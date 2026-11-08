"""Position-wise discrepancy analysis between groups of state sequences.

R-faithful port of TraMineR 2.2-14 ``seqdiff`` (``SOURCES/TraMineR/R/seqdiff.R``):
for each position window (driven by ``cmprange``) it computes a within-window
``seqdist`` and runs ``dissassoc`` on that distance matrix and the current group.
The complete statistics and discrepancy matrices are anchored to the exact R
oracle by ``tests/test_r_goldens_seqdiff.py``.

The public default is TraMineR's LCS distance with ``norm="auto"``. The current
distance layer implements its geometric-mean normalization. Explicit
``norm=None`` or ``norm="none"`` remains the maintained Python request for an
unnormalized distance.

Under ``with_missing=True``, the missing state ``nr`` is kept as a regular state
and structural ``void`` is deleted by ``seqdist`` through the window
``STSList`` metadata. With ``with_missing=False``, rows containing either marker
in the current window are excluded.
"""

from __future__ import annotations

import warnings
from collections.abc import Iterable, Mapping, Sequence
from numbers import Integral
from typing import Any

import numpy as np
import pandas as pd

from ._categorical_order import r_character_order, r_factor_levels
from ._deprecated import _caller_stacklevel
from ._rformat import r_as_character
from .sequences import STSList, _remap_issues


def _time_labels(x: STSList) -> list[str]:
    """R's ``colnames(seqdata)``: the explicit column names, else ``T<start..>``.

    Mirrors ``SOURCES/TraMineR/R/seqdef.R:249-257``.
    """
    if x.xtlab is not None and len(x.xtlab) >= x.seqlen:
        return [str(label) for label in x.xtlab[: x.seqlen]]
    return [f"T{x.start + offset}" for offset in range(x.seqlen)]


def _validate_cmprange(cmprange: object) -> tuple[int, int]:
    if cmprange is None or isinstance(cmprange, (str, bytes)):
        raise TypeError("cmprange must be a sequence of exactly two integer offsets")
    if not isinstance(cmprange, Iterable):
        raise TypeError("cmprange must be a sequence of exactly two integer offsets")
    values = list(cmprange)
    if len(values) != 2:
        raise ValueError("cmprange must contain exactly two integer offsets")
    for value in values:
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (Integral, np.integer)):
            raise TypeError("cmprange offsets must be integers")
    lower, upper = int(values[0]), int(values[1])
    if lower > upper:
        raise ValueError("cmprange lower offset must not exceed its upper offset")
    return lower, upper


def _resolve_seqdist_args(
    *,
    x: STSList,
    method: str,
    norm: str | None,
    seqdist_args: Mapping[str, Any] | None,
    seqdist_arg: Mapping[str, Any] | None,
) -> dict[str, Any]:
    if seqdist_args is not None and seqdist_arg is not None:
        raise ValueError("seqdist_args and deprecated seqdist_arg cannot be supplied together")

    supplied_map: Mapping[str, Any] | None
    if seqdist_arg is not None:
        warnings.warn(
            "seqdist_arg is deprecated; use seqdist_args instead",
            DeprecationWarning,
            stacklevel=_caller_stacklevel(),
        )
        supplied_map = seqdist_arg
    else:
        supplied_map = seqdist_args

    if supplied_map is not None and not isinstance(supplied_map, Mapping):
        raise TypeError("seqdist_args must be a mapping of seqdist options")

    if supplied_map is None:
        options: dict[str, Any] = {"method": method, "norm": norm}
    else:
        if method != "LCS" or norm != "auto":
            raise ValueError(
                "direct method or norm cannot be mixed with seqdist_args; "
                "put method and norm in the map"
            )
        options = dict(supplied_map)
        if "method" not in options:
            raise TypeError("seqdist_args must include the required seqdist method")

    if "sm" in options and "subm" in options:
        raise ValueError("seqdist_args cannot contain both sm and subm")
    if "sm" in options:
        options["subm"] = options.pop("sm")
    if isinstance(options.get("subm"), pd.DataFrame):
        options["subm"] = options["subm"].to_numpy(dtype=float, copy=True)

    for reserved_input in ("seqdata", "x"):
        if reserved_input in options:
            raise ValueError(
                f"seqdiff's current window owns the distance input; remove {reserved_input!r}"
            )

    if "void" in options:
        supplied_void = options.pop("void")
        if supplied_void != x.void_code:
            raise ValueError(
                "seqdiff reserves void for the current window's structural-void metadata"
            )

    if "refseq" in options or "ref_sets" in options:
        raise ValueError("seqdiff requires a square distance matrix; reference outputs are invalid")
    if bool(options.get("condensed", False)):
        raise ValueError("seqdiff requires a square distance matrix; condensed output is invalid")
    if options.get("full_matrix", True) is False:
        raise ValueError("seqdiff requires a square distance matrix; full_matrix=False is invalid")

    # R overwrites seqdist.args$with.missing from the direct seqdiff argument.
    options.pop("with_missing", None)
    options.pop("with.missing", None)
    return options


def _subset_window(
    x: STSList,
    row_indices: np.ndarray,
    column_indices: np.ndarray,
) -> STSList:
    rows = np.asarray(row_indices, dtype=int).reshape(-1)
    columns = np.asarray(column_indices, dtype=int).reshape(-1)
    return STSList(
        data=x.data[np.ix_(rows, columns)].copy(),
        alphabet=list(x.alphabet),
        labels=None if x.labels is None else list(x.labels),
        cpal=None if x.cpal is None else list(x.cpal),
        xtlab=None if x.xtlab is None else [x.xtlab[int(column)] for column in columns],
        xtstep=x.xtstep,
        tick_last=x.tick_last,
        nr=x.nr,
        void=x.void,
        nr_code=x.nr_code,
        void_code=x.void_code,
        weights=(None if x.weights is None else np.asarray(x.weights, dtype=float)[rows].copy()),
        start=x.start + (int(columns[0]) if columns.size else 0),
        ids=None if x.ids is None else [x.ids[int(row)] for row in rows],
        missing_color=x.missing_color,
        issues=_remap_issues(x.issues, rows),
    )


def _ordinary_group(
    group: object,
    n_sequences: int,
) -> tuple[np.ndarray, list[str]]:
    values = np.asarray(group)
    if values.ndim != 1 or values.shape[0] != n_sequences:
        raise ValueError("group must be one-dimensional with one value per sequence")

    if isinstance(group, pd.Series) and isinstance(group.dtype, pd.CategoricalDtype):
        levels = list(group.cat.remove_unused_categories().cat.categories)
    elif isinstance(group, pd.Categorical):
        # R's factor(group) drops unused levels (TX-AUD-20261010-026).
        levels = list(group.remove_unused_categories().categories)
    else:
        observed = values[~pd.isna(values)]
        observed_values = observed.tolist()
        if all(isinstance(value, (str, np.str_)) for value in observed_values):
            levels = r_character_order([str(value) for value in observed_values])
        else:
            # Numbers by value, booleans False first; text mixed with numbers is
            # refused, as by the other group functions (TX-AUD-20261010-096).
            levels = r_factor_levels(observed_values)

    def factor_label(level: object) -> str:
        # R's as.character: FALSE/TRUE, 1e+05, 0.3 (TX-AUD-20261010-026).
        text = r_as_character(level)
        return str(level) if text is None else text

    level_labels = [factor_label(level) for level in levels]
    labeled_values = np.full(n_sequences, None, dtype=object)
    for index, value in enumerate(values):
        if pd.isna(value):
            continue
        for level, label in zip(levels, level_labels, strict=True):
            if value == level:
                labeled_values[index] = label
                break
    return labeled_values, level_labels


def _sts_group_levels(group: STSList) -> tuple[dict[int, str], list[str]]:
    marker_codes = {int(code) for code in (group.nr_code, group.void_code) if code is not None}
    code_labels = {
        code: str(label) for code, label in enumerate(group.alphabet) if code not in marker_codes
    }
    return code_labels, list(code_labels.values())


def _sts_group_at_position(
    group: STSList, position: int, code_labels: dict[int, str]
) -> np.ndarray:
    codes = group.data[:, position - 1]
    values = np.full(group.nseq, None, dtype=object)
    for code, label in code_labels.items():
        values[codes == code] = label
    return values


def seqdiff(
    x: STSList,
    group: np.ndarray | pd.Series | STSList,
    *,
    cmprange: Sequence[int] = (0, 1),
    method: str = "LCS",
    norm: str | None = "auto",
    seqdist_args: Mapping[str, Any] | None = None,
    with_missing: bool = False,
    weighted: bool = True,
    squared: bool = False,
    seqdist_arg: Mapping[str, Any] | None = None,
) -> dict[str, pd.DataFrame]:
    """Run a sliding-window discrepancy analysis of sequences by group.

    ``seqdist_args`` is the complete R-style distance option map. It accepts
    Python ``subm`` and R ``sm`` spellings, while the direct ``with_missing``
    argument always wins. A supplied map must contain ``method``, as R
    ``seqdist`` has no method default. Direct ``method``/``norm`` remain
    conveniences when no map is supplied. ``cmprange`` must describe at least
    one complete, non-reversed window. ``seqdiff(weighted=...)`` controls only
    ``dissassoc``; an independent ``seqdist_args["weighted"]`` controls
    dataset-derived distance costs.
    """
    from .diss import dissassoc
    from .distances import seqdist
    from .distances.python_kernels import CostDiagnosticMessage

    if not isinstance(x, STSList):
        raise TypeError("seqdiff requires a state sequence object (seqdef)")

    c0, c1 = _validate_cmprange(cmprange)
    distance_options = _resolve_seqdist_args(
        x=x,
        method=method,
        norm=norm,
        seqdist_args=seqdist_args,
        seqdist_arg=seqdist_arg,
    )

    n_sequences, sequence_length = x.data.shape
    lo = max(1, 1 - c0)
    hi = min(sequence_length, sequence_length - c1)
    if lo > hi:
        raise ValueError("cmprange does not define any complete comparison window")
    positions = list(range(lo, hi + 1))

    group_sts = group if isinstance(group, STSList) else None
    if group_sts is not None:
        if group_sts.nseq != n_sequences:
            raise ValueError("group must have one value per sequence")
        if positions and group_sts.seqlen < max(positions):
            raise ValueError("an STSList group must have enough positions for every comparison")
        group_code_labels, level_labels = _sts_group_levels(group_sts)
        ordinary_values = None
    else:
        ordinary_values, level_labels = _ordinary_group(group, n_sequences)
        group_code_labels = {}

    discrepancy_weights = None
    if weighted and x.weights is not None:
        discrepancy_weights = np.asarray(x.weights, dtype=float)

    marker_codes = {int(code) for code in (x.nr_code, x.void_code) if code is not None}
    stat_columns = ["Pseudo F", "Pseudo Fbf", "Pseudo R2", "Bartlett", "Levene"]
    stat_rows: list[np.ndarray] = []
    discrepancy_rows: list[np.ndarray] = []

    for position in positions:
        columns = np.arange(position + c0 - 1, position + c1, dtype=int)
        if group_sts is not None:
            current_group = _sts_group_at_position(group_sts, position, group_code_labels)
        else:
            assert ordinary_values is not None
            current_group = ordinary_values

        keep = ~pd.isna(current_group)
        if not with_missing and marker_codes:
            raw_window = x.data[:, columns]
            keep &= ~np.isin(raw_window, list(marker_codes)).any(axis=1)
        row_indices = np.flatnonzero(keep)
        if row_indices.size == 0:
            # R hands seqdist a window without a sequence, and its matrix(nrow =
            # nbseq, ncol = maxsl) stops (seqdiff.R:51-62; TX-AUD-20261010-096).
            raise ValueError(
                "invalid 'ncol' value (too large or NA): a comparison window has no "
                "complete case"
            )
        window = _subset_window(x, row_indices, columns)
        current_group_subset = current_group[row_indices]

        call_options = dict(distance_options)
        distance_method = call_options.pop("method", "OM")
        with warnings.catch_warnings():
            # ``SOURCES/TraMineR/R/seqdiff.R:62`` wraps its internal seqdist call in
            # ``suppressMessages(...)``, which silences R's ``message()`` diagnostics
            # but not its ``warning()`` ones.  Silence the same channel here, so a
            # cost matrix that violates the triangle inequality is reported once by
            # a direct ``seqdist`` call rather than once per sliding window.
            warnings.simplefilter("ignore", CostDiagnosticMessage)
            # R's msg.warn notices are message()s too (" [!!] ...";
            # TraMineR-logging_helpers.R:20-25): "with.missing set as FALSE", "Empty
            # sequence(s)" (TX-AUD-20261019-059) and "'step' has been ignored because
            # 'breaks' is specified" (batch-9 review M-2).
            warnings.filterwarnings("ignore", message=r" \[!!\] ", category=UserWarning)
            distances = np.asarray(
                seqdist(
                    window,
                    method=distance_method,
                    with_missing=with_missing,
                    **call_options,
                ),
                dtype=float,
            )
        expected_shape = (len(row_indices), len(row_indices))
        if distances.ndim != 2 or distances.shape != expected_shape:
            raise ValueError(
                "seqdiff requires seqdist to return a square distance matrix "
                f"of shape {expected_shape}, got {distances.shape}"
            )

        window_weights = None if discrepancy_weights is None else discrepancy_weights[row_indices]
        result = dissassoc(
            distances,
            current_group_subset,
            weights=window_weights,
            R=0,
            weight_permutation="diss",
            squared=squared,
        )
        stat_rows.append(result.stat["t0"].reindex(stat_columns).to_numpy(dtype=float))

        discrepancy = result.groups["discrepancy"]
        row = np.zeros(len(level_labels) + 1, dtype=float)
        for index, label in enumerate(level_labels):
            if label in discrepancy.index:
                row[index] = float(discrepancy.loc[label])
        row[-1] = float(discrepancy.loc["Total"])
        discrepancy_rows.append(row)

    # R labels each comparison window with the time label of its position:
    # ``rownames(ret$stat) <- colnames(seqdata)[totrange]``
    # (SOURCES/TraMineR/R/seqdiff.R:29,32). A bare integer index dropped the
    # object's ``xtlab`` entirely, so calendar column names were lost on the way
    # out (cycle-6 §9). ``seqdef.R:254`` gives the default names
    # ``paste("T", start:(max(seql)+start-1))``, which is what an object without
    # explicit column names carries in R.
    row_labels = _time_labels(x)
    window_labels = [row_labels[position - 1] for position in positions]
    stat = pd.DataFrame(stat_rows, columns=stat_columns, index=window_labels)
    discrepancy = pd.DataFrame(
        discrepancy_rows,
        columns=level_labels + ["Total"],
        index=window_labels,
    )
    return {"stat": stat, "discrepancy": discrepancy}
