from __future__ import annotations

import contextlib
import io
import warnings
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from numbers import Real
from typing import Any, cast

import numpy as np
import pandas as pd

from ._categorical_order import r_character_order, stable_sorted_unique
from ._deprecated import _caller_stacklevel
from ._rformat import r_as_character
from .sequences import STSList, _require_stslist, alphabet


def _group_p_factor_label(value: object) -> str:
    """R's ``as.character`` of a level: ``TRUE``, ``1e+05``, ``0.3`` (TX-AUD-20261010-049)."""
    text = r_as_character(value)
    return str(value) if text is None else text


def _group_p_percent_label(value: float) -> str:
    if np.isnan(value):
        return "NaN"
    if np.isposinf(value):
        return "Inf"
    if np.isneginf(value):
        return "-Inf"
    rounded = float(np.round(value, 1))
    if rounded == 0:
        rounded = 0.0
    return f"{rounded:.1f}".removesuffix(".0")


def group_p(group: Sequence[object], weights: Sequence[float] | None = None) -> pd.Categorical:
    """Append weighted percentages to factor levels (extras ``group.p``).

    The input is converted to R-factor order: ordinary values are sorted as
    R's ``factor()`` sorts them (numbers by value, text in R's collation, which
    needs ICU; :func:`~traminex._categorical_order.stable_sorted_unique`),
    while an existing categorical preserves
    its declared order after unused levels are dropped. Missing group values do
    not contribute to the weighted table. The returned categorical is new and
    never mutates ``group`` or ``weights``.
    """
    group_dtype = getattr(group, "dtype", None)
    if isinstance(group, pd.Categorical) or isinstance(group_dtype, pd.CategoricalDtype):
        factor = pd.Categorical(group).remove_unused_categories()
    else:
        factor = pd.Categorical(group, categories=stable_sorted_unique(list(group)))
    codes = np.asarray(factor.codes, dtype=int).copy()
    if weights is None:
        resolved_weights: np.ndarray = np.ones(len(factor), dtype=float)
    else:
        resolved_weights = np.asarray(weights, dtype=float)
        if resolved_weights.ndim != 1 or resolved_weights.shape[0] != len(factor):
            raise ValueError("weights must have the same length as group")

    counts: np.ndarray = np.zeros(len(factor.categories), dtype=float)
    valid = codes >= 0
    np.add.at(counts, codes[valid], resolved_weights[valid])
    total = float(counts.sum())
    if np.isnan(resolved_weights[valid]).any():
        percentages = ["NA"] * len(counts)
    else:
        with np.errstate(divide="ignore", invalid="ignore"):
            raw_percentages = 100.0 * counts / total
        percentages = [_group_p_percent_label(value) for value in raw_percentages]

    levels = [
        f"{_group_p_factor_label(level)} ({percentage}%)"
        for level, percentage in zip(factor.categories, percentages, strict=True)
    ]
    return pd.Categorical.from_codes(codes, categories=levels, ordered=False)


def sorti(x: STSList, start: str = "end", sort_index: bool = True) -> np.ndarray:
    """Return the TraMineRextras lexicographic sequence ordering.

    ``start="beg"`` gives the first position highest priority; ``"end"``
    gives the final position highest priority.  Structural void is compared as
    missing, on a copy, matching R without mutating the caller.  Returned row
    indexes are zero-based; ``sort_index=False`` returns the inverse
    permutation used by sequence index plots.
    """
    if not isinstance(x, STSList):
        raise TypeError("x must be an STSList state sequence object")
    if start not in {"beg", "end"}:
        warnings.warn(
            'start should be one of "beg" or "end" ',
            UserWarning,
            stacklevel=2,
        )
        return np.arange(x.seqlen, dtype=np.int64)
    if x.seqlen == 0:
        raise ValueError("cannot sort zero-width sequence data")

    # R ranks cells over c(alphabet, nr, void) after replacing void by nr
    # (TraMineRextras sortv.R:5), whatever codes the markers are stored under
    # (TX-AUD-20261011-006). The ranks are a new int64 array: x is not touched.
    values = x._canonical_code_ranks(void_as_missing=True)[np.asarray(x.data)]

    priority = range(x.seqlen) if start == "beg" else range(x.seqlen - 1, -1, -1)
    keys = tuple(values[:, column] for column in reversed(tuple(priority)))
    indexes = np.asarray(np.lexsort(keys), dtype=np.int64)
    if sort_index:
        return indexes
    return np.asarray(np.argsort(indexes, kind="stable"), dtype=np.int64)


def sortv(x: STSList, start: str = "end") -> np.ndarray:
    """Return inverse lexicographic sequence ranks for plotting.

    The ranks are zero-based, where TraMineRextras' ``sortv`` returns one-based ranks
    (``order(sorti(...))``); both order the sequences alike, which is all a ``sortv=``
    argument uses. ``start`` is :func:`sorti`'s.
    """
    return sorti(x, start=start, sort_index=False)


def _rowmode_label(value: object) -> str:
    """R's ``as.character`` of a value: ``TRUE``, ``1e+05``, ``0.3`` (TX-AUD-20261010-049)."""
    text = r_as_character(value)
    return str(value) if text is None else text


def _rowmode_is_missing(value: object) -> bool:
    missing = pd.isna(value)
    return isinstance(missing, (bool, np.bool_)) and bool(missing)


def _rowmode_exclusions(except_: object) -> set[str]:
    if except_ is None:
        return set()
    values: list[Any]
    if isinstance(except_, (str, bytes)) or np.isscalar(except_):
        values = [except_]
    else:
        try:
            values = list(cast(Any, except_))
        except TypeError:
            values = [except_]
    return {_rowmode_label(value) for value in values if not _rowmode_is_missing(value)}


def _rowmode_scalar(v: Any, except_: object) -> str | None:
    categorical: pd.Categorical | None = None
    dtype = getattr(v, "dtype", None)
    if isinstance(v, pd.Categorical) or isinstance(dtype, pd.CategoricalDtype):
        categorical = pd.Categorical(v)
        values = list(categorical)
    else:
        values = np.asarray(v, dtype=object).tolist()

    observed = [value for value in values if not _rowmode_is_missing(value)]
    if not observed:
        return None

    counts: dict[str, int] = {}
    for value in observed:
        label = _rowmode_label(value)
        counts[label] = counts.get(label, 0) + 1

    for label in _rowmode_exclusions(except_):
        counts.pop(label, None)
    if not counts:
        return None

    if categorical is not None:
        ordered = [_rowmode_label(value) for value in categorical.categories]
    elif all(isinstance(value, (bool, np.bool_)) for value in observed):
        ordered = [_rowmode_label(value) for value in sorted(set(observed))]
    elif all(
        isinstance(value, Real) and not isinstance(value, (bool, np.bool_)) for value in observed
    ):
        ordered = [_rowmode_label(value) for value in sorted(set(observed))]
    else:
        ordered = r_character_order(list(counts))

    maximum = max(counts.values())
    return next(label for label in ordered if counts.get(label) == maximum)


def rowmode(v: Any, except_: Any = None) -> str | None | np.ndarray:
    """Return the TraMineRextras scalar mode as a character label.

    Missing values are omitted, exclusions are applied as whole values, and
    ties follow R table order (including pandas categorical level order and
    the supported R/ICU character collation). R's ``NULL`` and zero-length
    character no-result forms are both stably re-expressed as ``None``.

    A two-dimensional input retains the established Python extension: the
    corrected scalar operation is applied independently to each row and an
    object array is returned. Public R treats a matrix as one variable.
    """
    arr = np.asarray(v)
    if arr.ndim == 1:
        return _rowmode_scalar(v, except_)
    if arr.ndim == 2:
        return np.array([_rowmode_scalar(row, except_) for row in arr], dtype=object)
    raise ValueError("rowmode expects a 1D or 2D array-like")


_SEQTABSTOCC_DEFAULT_IDXS = tuple(range(1, 11))


class _OmittedSeqtabstoccIncludeMissing:
    def __repr__(self) -> str:
        return "False"


class _OmittedSeqtabstoccOption:
    def __repr__(self) -> str:
        return "None"


_SEQTABSTOCC_INCLUDE_MISSING_OMITTED = _OmittedSeqtabstoccIncludeMissing()
_SEQTABSTOCC_WITH_MISSING_OMITTED = _OmittedSeqtabstoccOption()
_SEQTABSTOCC_TLIM_OMITTED = _OmittedSeqtabstoccOption()
_SEQTABSTOCC_FORMAT_OMITTED = _OmittedSeqtabstoccOption()


def _seqtabstocc_missing_option(include_missing: object, with_missing: object) -> bool:
    include_supplied = include_missing is not _SEQTABSTOCC_INCLUDE_MISSING_OMITTED
    with_supplied = with_missing is not _SEQTABSTOCC_WITH_MISSING_OMITTED
    if include_supplied and with_supplied:
        raise ValueError("include_missing and with_missing cannot be specified together")
    value = with_missing if with_supplied else include_missing if include_supplied else False
    if not isinstance(value, (bool, np.bool_)):
        raise TypeError("include_missing/with_missing must be logical")
    return bool(value)


def _seqtabstocc_unique_row_names(patterns: Sequence[str]) -> list[str]:
    """Re-express R data-frame duplicate row-name suffixing."""
    used: set[str] = set()
    names: list[str] = []
    for pattern in patterns:
        candidate = pattern
        suffix = 1
        while candidate in used:
            candidate = f"{pattern}.{suffix}"
            suffix += 1
        used.add(candidate)
        names.append(candidate)
    return names


def _seqtabstocc_intermediate(x: STSList, *, include_missing: bool) -> STSList:
    """Build the unweighted, row-sorted DSS sequence used by the R wrapper."""
    excluded = {code for code in (x.void_code,) if code is not None}
    if not include_missing and x.nr_code is not None:
        excluded.add(int(x.nr_code))

    row_labels: list[list[str]] = []
    observed: list[str] = []
    for row in x.data:
        labels = r_character_order(
            [str(x._state_for_code(int(code))) for code in row if int(code) not in excluded]
        )
        row_labels.append(labels)
        observed.extend(labels)

    states = r_character_order(observed)
    void = x.void
    while void in states:
        void = f"{void}%"
    nr = x.nr
    while nr in states or nr == void:
        nr = f"{nr}*"
    alphabet = [*states, void]
    codebook = {state: code for code, state in enumerate(alphabet)}
    width = max((len(labels) for labels in row_labels), default=0)
    data = np.full((x.nseq, max(width, 1)), codebook[void], dtype=int)
    for row, labels in enumerate(row_labels):
        for column, label in enumerate(labels):
            data[row, column] = codebook[label]

    sorted_sequences = STSList(
        data=data,
        alphabet=alphabet,
        labels=list(alphabet),
        cpal=None,
        xtlab=[f"T{position}" for position in range(1, data.shape[1] + 1)],
        xtstep=1,
        tick_last=False,
        nr=nr,
        void=void,
        nr_code=None,
        void_code=codebook[void],
        weights=None,
        start=1,
        ids=None if x.ids is None else list(x.ids),
        missing_color=x.missing_color,
        issues=None,
    )

    from .stats_complexity import seqdss

    return seqdss(sorted_sequences)


def seqtabstocc(
    x: STSList,
    include_missing: bool = _SEQTABSTOCC_INCLUDE_MISSING_OMITTED,  # type: ignore[assignment]
    sep: str = "-",
    *,
    with_missing: bool | None = _SEQTABSTOCC_WITH_MISSING_OMITTED,  # type: ignore[assignment]
    idxs: int | Sequence[int] = _SEQTABSTOCC_DEFAULT_IDXS,
    weighted: bool = True,
    tlim: int | Sequence[int] | None = _SEQTABSTOCC_TLIM_OMITTED,  # type: ignore[assignment]
    format: str | None = _SEQTABSTOCC_FORMAT_OMITTED,  # type: ignore[assignment]
) -> pd.DataFrame:
    """Return frequencies of each sequence's co-occurring state set.

    This follows TraMineRextras' row sort, sequence redefinition, DSS collapse,
    and fixed ``seqtab(format="STS")`` pipeline. Structural void is always
    omitted; missing is included only when requested. The wrapper-created
    intermediate is deliberately unweighted, so source weights are ignored
    even when ``weighted=True``.

    ``include_missing`` is the established Python spelling. The R-equivalent
    ``with_missing`` compatibility spelling is accepted without warning; using
    both spellings is an error. ``idxs`` (R's one-based frequency ranks) and deprecated
    ``tlim`` are forwarded
    to the corrected :func:`seqtab`. Forwarding ``format`` conflicts with the
    wrapper's fixed STS format, as in R. ``sep`` is a retained Python extension.

    The stable pandas view retains ``pattern``, ``count``, ``prop``, and
    ``cumprop``. Exact R-style ``Freq``/``Percent`` evidence and seqtab metadata
    are available in :attr:`DataFrame.attrs`.
    """
    if not isinstance(x, STSList):
        raise TypeError("x must be an STSList state sequence object")
    if format is not _SEQTABSTOCC_FORMAT_OMITTED:
        raise TypeError("format conflicts with seqtabstocc's fixed format='STS'")
    resolved_missing = _seqtabstocc_missing_option(include_missing, with_missing)

    from .stats_basic import seqtab

    intermediate = _seqtabstocc_intermediate(x, include_missing=resolved_missing)
    if tlim is _SEQTABSTOCC_TLIM_OMITTED:
        table = seqtab(
            intermediate,
            idxs=idxs,
            weighted=weighted,
            format="STS",
            sep=sep,
            include_missing=True,
        )
    else:
        table = seqtab(
            intermediate,
            idxs=idxs,
            weighted=weighted,
            format="STS",
            tlim=tlim,
            sep=sep,
            include_missing=True,
        )
    patterns = _seqtabstocc_unique_row_names(table["sequence"].astype(str).tolist())

    result = table.rename(columns={"sequence": "pattern"}).copy()
    result["pattern"] = patterns
    result.attrs.update(table.attrs)
    freq = table.attrs["freq"].copy()
    freq.index = pd.Index(patterns, name=None)
    result.attrs.update(
        {
            "freq": freq,
            "format": "STS",
            "source_weights_ignored": x.weights is not None,
            "python_sep": sep,
        }
    )
    return result


_SEQINDIC_DEFAULT = ("visited", "trans", "entr", "cplx", "turb2n")
_SEQINDIC_COLUMN_ORDER = (
    "lgth",
    "nonm",
    "dlgth",
    "visited",
    "visitp",
    "recu",
    "trans",
    "transp",
    "entr",
    "meand",
    "dustd",
    "meand2",
    "dustd2",
    "nsubs",
    "volat",
    "cplx",
    "turb",
    "turbn",
    "turb2",
    "turb2n",
    "ppos",
    "nvolat",
    "vpos",
    "integr",
    "degrad",
    "bad",
    "prec",
    "insec",
)
_SEQINDIC_BINARY = ("ppos", "nvolat", "vpos", "integr")
_SEQINDIC_RANKED = ("degrad", "bad", "prec", "insec")
_SEQINDIC_IPOS_ALIASES = {
    "seqdata": "x",
    "pos.states": "pos_states",
    "neg.states": "neg_states",
    "with.missing": "include_missing",
    "with_missing": "include_missing",
}
_SEQINDIC_PREC_ALIASES = {
    "seqdata": "x",
    "state.order": "state_order",
    "state.equiv": "state_equiv",
    "with.missing": "with_missing",
    "include_missing": "with_missing",
    "weight.type": "weight_type",
    "spell.integr": "spell_integr",
    "border.effect": "border_effect",
}


def _translate_seqindic_args(
    arguments: Mapping[str, Any] | None,
    aliases: Mapping[str, str],
    *,
    name: str,
) -> dict[str, Any]:
    if arguments is None:
        return {}
    if not isinstance(arguments, Mapping):
        raise TypeError(f"{name} must be a mapping")
    translated: dict[str, Any] = {}
    origins: dict[str, str] = {}
    for key, value in arguments.items():
        if not isinstance(key, str):
            raise TypeError(f"{name} keys must be strings")
        normalized = aliases.get(key, key)
        if normalized in translated:
            raise ValueError(f"{name} contains synonymous keys {origins[normalized]!r} and {key!r}")
        translated[normalized] = value
        origins[normalized] = key
    return translated


@_require_stslist
def seqindic(
    x: STSList,
    *,
    indic: str | Sequence[str] = _SEQINDIC_DEFAULT,
    with_missing: bool = False,
    ipos_args: Mapping[str, Any] | None = None,
    prec_args: Mapping[str, Any] | None = None,
    w: float = 0.5,
) -> pd.DataFrame:
    """Table of per-sequence indicators (TraMineR ``seqindic``).

    Returns a DataFrame with one column per requested indicator (R column
    labels) and one row per sequence, indexed by the sequence ids as R names its
    rows (``1..n`` for an object without ids). The default ``indic`` matches R: ``Visited`` (# distinct visited
    states), ``Trans`` (transitions), ``Entr`` (longitudinal entropy, normalized),
    ``Cplx`` (complexity index) and ``Turb2n`` (normalized turbulence, type 2).

    Supports all TraMineR indicator groups: basic, diversity, complexity,
    binary, and ranked. ``ipos_args`` supplies arguments for binary indicators;
    ``prec_args`` supplies arguments for ranked indicators. R-style dotted keys
    and Python-style underscored keys are accepted without mutating either
    caller-owned mapping.

    This is a re-target: the previous implementation returned a one-hot state
    indicator matrix, which is not what TraMineR ``seqindic`` computes.
    """
    x = x._canonicalized()
    from .stats_basic import seqistatd, seqtransn
    from .stats_complexity import (
        _r_logical_pow,
        _seqindic_seqibad,
        _seqindic_seqidegrad,
        _seqindic_seqinsecurity,
        _seqindic_seqprecarity,
        seqdss,
        seqipos,
        seqivardur,
        seqivolatility,
        seqST,
    )
    from .stats_p2 import seqici, seqient, seqsubsn

    groups = {
        "basic": [
            "lgth",
            "nonm",
            "dlgth",
            "visited",
            "visitp",
            "recu",
            "trans",
            "transp",
            "meand",
            "meand2",
        ],
        "diversity": ["visitp", "entr", "dustd", "dustd2"],
        "complexity": [
            "transp",
            "nsubs",
            "volat",
            "cplx",
            "turb",
            "turbn",
            "turb2",
            "turb2n",
        ],
        "binary": list(_SEQINDIC_BINARY),
        "ranked": list(_SEQINDIC_RANKED),
    }
    indicators = [indic] if isinstance(indic, str) else list(indic)
    aliases = {"visit": "visited", "inpos": "integr", "prec2": "insec"}
    indicators = [aliases.get(name, name) for name in indicators]
    valid = set(_SEQINDIC_COLUMN_ORDER) | set(groups) | {"all"}
    invalid = [name for name in indicators if name not in valid]
    if invalid:
        raise ValueError(f"invalid values in indic: {', '.join(invalid)}")

    # R expands groups in this fixed sequence. Existing individual indicators
    # stay first and each group's members are appended with unique() semantics;
    # final columns are nevertheless emitted in the source's fixed order below.
    group_expansions = {
        "all": groups["basic"] + groups["diversity"] + groups["complexity"],
        **groups,
    }
    requested = indicators
    for group in ("all", "basic", "diversity", "complexity", "binary", "ranked"):
        if group not in requested:
            continue
        requested = [name for name in requested if name != group]
        for name in group_expansions[group]:
            if name not in requested:
                requested.append(name)

    resolved_ipos_args = _translate_seqindic_args(
        ipos_args, _SEQINDIC_IPOS_ALIASES, name="ipos_args"
    )
    if any(name in requested for name in _SEQINDIC_BINARY):
        if not resolved_ipos_args:
            raise ValueError(
                "at least one selected indicator requires a non empty ipos_args mapping"
            )
        if resolved_ipos_args.get("x") is not None:
            warnings.warn(
                "seqdata argument in ipos.args will be overwritten!",
                UserWarning,
                stacklevel=2,
            )
        if resolved_ipos_args.get("include_missing") is not None:
            warnings.warn(
                "with.missing argument in ipos.args will be overwritten!",
                UserWarning,
                stacklevel=2,
            )
        if resolved_ipos_args.get("index") is not None:
            warnings.warn(
                "index argument in ipos.args will be overwritten!",
                UserWarning,
                stacklevel=2,
            )
        resolved_ipos_args["x"] = x
        resolved_ipos_args["include_missing"] = with_missing

    resolved_prec_args = _translate_seqindic_args(
        prec_args, _SEQINDIC_PREC_ALIASES, name="prec_args"
    )
    if any(name in requested for name in _SEQINDIC_RANKED):
        if not resolved_prec_args:
            warnings.warn(
                "No prec.args list. Alphabet order used as state order.",
                UserWarning,
                stacklevel=2,
            )
        if resolved_prec_args.get("x") is not None:
            warnings.warn(
                "seqdata argument in prec.args will be overwritten!",
                UserWarning,
                stacklevel=2,
            )
        if resolved_prec_args.get("with_missing") is not None:
            warnings.warn(
                "with.missing argument in prec.args will be overwritten!",
                UserWarning,
                stacklevel=2,
            )
        # This intentionally mirrors a TraMineR 2.2-14 source typo: the warning
        # tests ipos.args$type, while prec.args$type itself is silently replaced.
        if resolved_ipos_args.get("type") is not None:
            warnings.warn(
                "type argument in prec.args will be overwritten!",
                UserWarning,
                stacklevel=2,
            )
        resolved_prec_args["x"] = x
        resolved_prec_args["with_missing"] = with_missing

    wm = with_missing
    K = len(x.alphabet)
    statl = [i for i in range(K) if i != x.void_code and (wm or i != x.nr_code)]

    istatd = np.asarray(seqistatd(x, include_missing=wm), dtype=float)
    nvisit = (istatd[:, statl] > 0).sum(axis=1).astype(float)

    nonvoid = (x.data != x.void_code) if x.void_code is not None else np.ones_like(x.data, bool)
    lgth = nonvoid.sum(axis=1).astype(float)
    nonm_mask = nonvoid & (x.data != x.nr_code) if x.nr_code is not None else nonvoid
    nonm = nonm_mask.sum(axis=1).astype(float)

    dss = seqdss(x, include_missing=wm)
    if dss.void_code is not None:
        dlgth = (dss.data != dss.void_code).sum(axis=1).astype(float)
    else:
        dlgth = np.full(x.nseq, dss.data.shape[1], dtype=float)

    # Reuse the approved alphabet boundary: it adds a phantom nr when needed,
    # without duplicating marker text that is deliberately an ordinary state.
    n_alpha = max(len(alphabet(x, with_missing=wm)), 1)
    vardur: dict[int, Any] = {}

    def _vardur(type_: int):
        if type_ not in vardur:
            vardur[type_] = seqivardur(x, type=type_, include_missing=wm)
        return vardur[type_]

    def _ipos(index: str, *, dss: bool | None = None, force_dss: bool = False) -> np.ndarray:
        call_args = dict(resolved_ipos_args)
        call_args["index"] = index
        if force_dss:
            call_args["dss"] = dss
        return np.asarray(seqipos(**call_args), dtype=float)

    def _ranked_call_args(allowed: set[str], *, missing_name: str) -> dict[str, Any]:
        call_args = {key: value for key, value in resolved_prec_args.items() if key in allowed}
        ranked_with_missing = bool(call_args.pop("with_missing", wm))
        call_args[missing_name] = ranked_with_missing
        return call_args

    def _check_private_spell_integr_collision() -> None:
        # seqindic retains this private-formal name, but seqprecarity.private
        # supplies spell.integr explicitly when it calls seqidegrad.private.
        # R therefore raises a duplicate-argument error unless an explicit
        # correction bypasses that nested call altogether.
        if "spell_integr" in resolved_prec_args and resolved_prec_args.get("correction") is None:
            raise TypeError('formal argument "spell.integr" matched by multiple actual arguments')

    def _column(name: str):
        if name == "lgth":
            return "Lgth", lgth
        if name == "nonm":
            return "NonM", nonm
        if name == "dlgth":
            return "Dlgth", dlgth
        if name == "visited":
            return "Visited", nvisit
        if name == "visitp":
            return "Visitp", nvisit / n_alpha
        if name == "recu":
            return "Recu", np.divide(
                dlgth,
                nvisit,
                out=np.full_like(dlgth, np.nan),
                where=nvisit > 0,
            )
        if name == "trans":
            return "Trans", seqtransn(x, with_missing=wm, norm=False)
        if name == "transp":
            return "Transp", seqtransn(x, with_missing=wm, norm=True)
        if name == "entr":
            return "Entr", seqient(x, norm=True, include_missing=wm)
        if name == "meand":
            return "MeanD", _vardur(1).meand
        if name == "dustd":
            return "Dustd", np.sqrt(np.asarray(_vardur(1), dtype=float))
        if name == "meand2":
            return "MeanD2", _vardur(2).meand
        if name == "dustd2":
            return "Dustd2", np.sqrt(np.asarray(_vardur(2), dtype=float))
        if name == "nsubs":
            return "Nsubs", np.asarray(seqsubsn(x, with_missing=wm), dtype=float).reshape(-1)
        if name == "volat":
            return "Volat", seqivolatility(x, w=w, include_missing=wm)
        if name == "cplx":
            return "Cplx", seqici(x, include_missing=wm)
        if name == "turb":
            return "Turb", seqST(x, norm=False, include_missing=wm, type=1)
        if name == "turbn":
            return "Turbn", seqST(x, norm=True, include_missing=wm, type=1)
        if name == "turb2":
            return "Turb2", seqST(x, norm=False, include_missing=wm, type=2)
        if name == "turb2n":
            return "Turb2n", seqST(x, norm=True, include_missing=wm, type=2)
        if name == "ppos":
            return "Ppos", _ipos("share", dss=False, force_dss=True)
        if name == "nvolat":
            return "Nvolat", _ipos("share", dss=True, force_dss=True)
        if name == "vpos":
            return "Vpos", _ipos("volatility")
        if name == "integr":
            return "Integr", _ipos("integr")
        if name == "degrad":
            args = _ranked_call_args(
                {
                    "x",
                    "state_order",
                    "state_equiv",
                    "stprec",
                    "with_missing",
                    "penalized",
                    "method",
                    "weight_type",
                    "pow",
                    "border_effect",
                },
                missing_name="with_missing",
            )
            # R's public seqidegrad treats logical pow as the spell.integr
            # switch and then coerces it to numeric before the private call.
            pow_value = args.get("pow")
            if isinstance(pow_value, (bool, np.bool_)):
                args["pow"], args["spell_integr"] = _r_logical_pow(pow_value, None)
            return "Degrad", _seqindic_seqidegrad(**args)
        if name == "bad":
            args = _ranked_call_args(
                {
                    "x",
                    "pow",
                    "state_order",
                    "state_equiv",
                    "stprec",
                    "with_missing",
                },
                missing_name="include_missing",
            )
            return "Bad", _seqindic_seqibad(**args)
        if name == "prec":
            _check_private_spell_integr_collision()
            args = _ranked_call_args(
                {
                    "x",
                    "correction",
                    "otto",
                    "a",
                    "b",
                    "method",
                    "state_order",
                    "state_equiv",
                    "stprec",
                    "with_missing",
                    "penalized",
                    "weight_type",
                    "pow",
                    "border_effect",
                },
                missing_name="include_missing",
            )
            return "Prec", _seqindic_seqprecarity(**args)
        if name == "insec":
            _check_private_spell_integr_collision()
            args = _ranked_call_args(
                {
                    "x",
                    "correction",
                    "state_order",
                    "state_equiv",
                    "stprec",
                    "with_missing",
                    "pow",
                    "spow",
                    "bound",
                    "method",
                    "penalized",
                    "weight_type",
                    "border_effect",
                },
                missing_name="with_missing",
            )
            # seqindic calls seqprecarity.private directly, whose default is
            # spow=0. The public seqinsecurity helper defaults spow to pow.
            args.setdefault("spow", 0.0)
            return "Insec", _seqindic_seqinsecurity(**args)
        raise AssertionError(f"unhandled validated seqindic indicator: {name}")

    out = {}
    for name in _SEQINDIC_COLUMN_ORDER:
        if name not in requested:
            continue
        # Preserve TraMineR 2.2-14's nesting bug: Visitp is computed only when
        # dlgth, visited, or recu also opens the surrounding source branch.
        if name == "visitp" and not any(gate in requested for gate in ("dlgth", "visited", "recu")):
            continue
        label, vals = _column(name)
        out[label] = np.asarray(vals, dtype=float)
    # R names the rows by rownames(seqdata): the sequence ids, 1..n for an object
    # without ids (TX-AUD-20261010-093).
    index = list(x.ids) if x.ids is not None else list(range(1, x.nseq + 1))
    return pd.DataFrame(out, index=index)


@_require_stslist
def seqindic_dyn(
    x: STSList,
    *,
    indic: str | Sequence[str] = "cplx",
    window_size: float = 0.2,
    sliding: bool = True,
    wstep: int = 1,
    with_missing: bool = False,
    endmiss_as_void: bool = False,
    silent_indic: bool = True,
    **kwargs,
) -> np.ndarray:
    """Dynamic (windowed) sequence indicator (TraMineRextras ``seqindic.dyn``).

    Computes :func:`seqindic` over a sliding (``[k-step+1 : k]``) or incremental
    (``[1 : k]``) window for each window-end position ``k`` in
    ``seq(step, maxl, by=wstep)``, returning the ``dynin`` matrix (n_sequences x
    n_windows). ``window_size < 1`` is a fraction of the longest sequence
    (``ceil(window_size * maxl)``); otherwise it is the window width. A cell is
    ``NaN`` when the sequence's window-end position is void (or missing, when
    ``endmiss_as_void`` and not ``with_missing``). The ``plot.dynin`` visualization
    is RP7. ``silent_indic=True`` suppresses messages from the per-window
    :func:`seqindic` calls. Provenance: TraMineRextras 0.6.9 ``seqindic-dyn.R``.

    ``indic`` names one indicator of :func:`seqindic`. As in R, an indicator
    group name (``"all"``, ``"basic"``, ``"diversity"``, ``"complexity"``,
    ``"binary"``, ``"ranked"``) is refused, and of a sequence of names only the
    first is used, with R's message as a ``UserWarning``; a name seqindic does
    not compute is refused by name (TX-AUD-20261010-051). Row ``i`` is sequence
    ``i``; column ``j`` is the window ending at position ``k_j``, which R names
    after that position's column label.
    """
    import math

    if not isinstance(indic, str):
        names = list(indic)
        if not names:
            raise ValueError("indic must name an indicator")
        if names[0] == "Develop":
            raise ValueError("indic 'Develop' is not supported by seqindic_dyn")
        if len(names) > 1:
            warnings.warn("Vector indic, only first value is used!", UserWarning, stacklevel=2)
        indic = str(names[0])
    if indic in {"all", "basic", "diversity", "complexity", "binary", "ranked"}:
        raise ValueError("Bad indic value, group name not supported!")
    if not isinstance(silent_indic, (bool, np.bool_)):
        raise TypeError("silent_indic must be a boolean")
    if wstep < 1:
        raise ValueError("wstep must be a strictly positive integer")
    if window_size <= 0:
        raise ValueError("bad window_size value")

    if x.nseq == 0:
        # R's max() of no lengths is -Inf, so its window count is NaN and matrix()
        # stops (seqindic-dyn.R:29, 72; TX-AUD-20261010-096).
        raise ValueError("invalid 'ncol' value (too large or NA): x has no sequences")
    if x.void_code is not None:
        lgth = (x.data != x.void_code).sum(axis=1)
    else:
        lgth = np.full(x.nseq, x.seqlen, dtype=int)
    maxl = int(lgth.max())
    step = math.ceil(window_size * maxl) if window_size < 1 else min(int(window_size), maxl)

    re_range = list(range(step, maxl + 1, wstep))
    if len(re_range) < 2:
        raise ValueError("There is only one window, wstep probably too large!")

    miss_codes: set[int] = set()
    if x.void_code is not None:
        miss_codes.add(x.void_code)
    if endmiss_as_void and not with_missing and x.nr_code is not None:
        miss_codes.add(x.nr_code)

    ind_dyn = np.full((x.nseq, len(re_range)), np.nan, dtype=float)
    for j, k in enumerate(re_range):
        lo = (k - step) if sliding else 0  # 0-based start of window
        sub = STSList(
            data=x.data[:, lo:k].copy(),
            alphabet=list(x.alphabet),
            labels=None if x.labels is None else list(x.labels),
            cpal=None if x.cpal is None else list(x.cpal),
            nr=x.nr,
            void=x.void,
            nr_code=x.nr_code,
            void_code=x.void_code,
            weights=None if x.weights is None else x.weights.copy(),
            start=x.start,
            xtlab=None if x.xtlab is None else list(x.xtlab[lo:k]),
            xtstep=x.xtstep,
            tick_last=x.tick_last,
            ids=None if x.ids is None else list(x.ids),
            missing_color=x.missing_color,
            issues=x.issues,
        )
        if silent_indic:
            with contextlib.redirect_stderr(io.StringIO()):
                frame = seqindic(sub, indic=[indic], with_missing=with_missing, **kwargs)
        else:
            frame = seqindic(sub, indic=[indic], with_missing=with_missing, **kwargs)
        if frame.shape[1] == 0:
            # R fails here with "replacement has length zero".
            raise ValueError(f"indic {indic!r} is not an indicator seqindic computes")
        ind_dyn[:, j] = np.asarray(frame.iloc[:, 0], dtype=float)
        # NA for sequences whose window-end position is void / end-missing
        end_codes = x.data[:, k - 1]
        wna = np.isin(end_codes, list(miss_codes)) if miss_codes else np.zeros(x.nseq, bool)
        ind_dyn[wna, j] = np.nan
    return ind_dyn


@_require_stslist
def seqasnum(x: STSList, *, with_missing: bool = False) -> np.ndarray:
    """Numeric coding of the sequences (TraMineR ``seqasnum``).

    Returns a float matrix of 0-based codes over the real-state alphabet -- code
    ``i-1`` for the ``i``-th real state, exactly as R. Cells that are not
    states are ``NaN``: the structural ``void`` symbol always, and the missing
    code (``nr``) unless ``with_missing`` (then ``nr`` is coded immediately
    after the real alphabet), independently of stored marker positions.
    """
    x = x._canonicalized()
    real_codes = [code for code in range(len(x.alphabet)) if code not in (x.nr_code, x.void_code)]
    lookup = np.full(len(x.alphabet), np.nan, dtype=float)
    lookup[real_codes] = np.arange(len(real_codes), dtype=float)
    if with_missing and x.nr_code is not None:
        lookup[x.nr_code] = len(real_codes)
    out = lookup[x.data]
    return out


@_require_stslist
def seqnum(x: STSList, *, with_missing: bool = False) -> STSList:
    """Recode the sequences over a numeric alphabet (TraMineR ``seqnum``).

    Returns a copy of ``x`` whose real-state labels are the numeric codes
    ``0..K-1`` (R's ``seqnum``: alphabet -> ``0:(nbstat-1)``); the underlying
    integer codes are unchanged; :func:`seqasnum` maps those codes to real-state
    ranks. The
    structural ``void`` symbol is preserved; the missing code (``nr``) is also
    relabelled numerically when ``with_missing`` (mirroring R), otherwise kept.
    """
    x = x._canonicalized()
    real_rank = 0
    nr = str(len(alphabet(x))) if with_missing else x.nr
    new_alph: list[str] = []
    for i, lab in enumerate(x.alphabet):
        if i == x.void_code:
            new_alph.append(lab)
        elif i == x.nr_code:
            new_alph.append(nr)
        else:
            new_alph.append(str(real_rank))
            real_rank += 1
    return replace(x, alphabet=new_alph, nr=nr)


@dataclass(frozen=True, eq=False)
class SeqTrateResult:
    """Labelled K by K by time result for ``seqtrate(time_varying=True)``."""

    values: np.ndarray
    from_states: tuple[str, ...]
    to_states: tuple[str, ...]
    time_labels: tuple[str, ...]

    @property
    def shape(self) -> tuple[int, int, int]:
        return self.values.shape

    @property
    def ndim(self) -> int:
        return self.values.ndim

    def __array__(self, dtype=None, copy=None) -> np.ndarray:
        out = np.asarray(self.values, dtype=dtype)
        return out.copy() if copy is True else out

    def __getitem__(self, key):
        return self.values[key]


def _seqtrate_states(
    x: STSList,
    sel_states: Sequence[str] | None,
    include_missing: bool,
) -> tuple[list[str], list[int | None]]:
    """Resolve R's ``sel.states`` labels to this STSList's integer codes."""
    if sel_states is None:
        marker_codes = {code for code in (x.nr_code, x.void_code) if code is not None}
        labels = [label for code, label in enumerate(x.alphabet) if code not in marker_codes]
    elif isinstance(sel_states, str):
        labels = [sel_states]
    else:
        labels = [str(label) for label in sel_states]

    # R appends nr even when it is absent from the data (and even when the user
    # already selected it), so do not de-duplicate this list.
    if include_missing:
        labels.append(str(x.nr))

    label_to_code: dict[str, int] = {}
    for code, label in enumerate(x.alphabet):
        label_to_code.setdefault(str(label), code)
    codes = [label_to_code.get(label) for label in labels]
    if include_missing:
        codes[-1] = x.nr_code
    return labels, codes


def _seqtrate_accumulate(
    x: STSList,
    state_codes: Sequence[int | None],
    origins: np.ndarray,
    destinations: np.ndarray,
    weights: np.ndarray,
    *,
    include_missing: bool,
    count: bool,
) -> np.ndarray:
    """Accumulate one stationary matrix using TraMineR's exact denominator."""
    nstates = len(state_codes)
    numerators: np.ndarray = np.zeros((nstates, nstates), dtype=float)
    denominators: np.ndarray = np.zeros(nstates, dtype=float)

    source = x.data[:, origins]
    target = x.data[:, destinations]
    eligible_target = np.ones(target.shape, dtype=bool)
    if x.void_code is not None:
        eligible_target &= target != x.void_code
    if not include_missing and x.nr_code is not None:
        eligible_target &= target != x.nr_code
    row_weights = np.asarray(weights, dtype=float)

    def _in_row_order(per_row_counts: np.ndarray) -> float:
        # Weighted per-row counts summed row by row from 0.0, exactly as the
        # former per-row loop did: np.add.accumulate is sequential (np.sum is
        # pairwise), so the result is bit-identical (audit S3-18).
        terms = np.concatenate(([0.0], row_weights * per_row_counts.astype(float)))
        return float(np.add.accumulate(terms)[-1])

    target_matches = {code: target == code for code in state_codes if code is not None}
    for source_index, source_code in enumerate(state_codes):
        if source_code is None:
            continue
        source_match = source == source_code
        denominators[source_index] = _in_row_order(
            np.count_nonzero(source_match & eligible_target, axis=1)
        )
        for target_index, target_code in enumerate(state_codes):
            if target_code is None:
                continue
            # R deliberately does not apply missingcond to PAB.  For the
            # default state set this is equivalent, but it matters when a
            # caller explicitly selects a missing or void marker.
            numerators[source_index, target_index] = _in_row_order(
                np.count_nonzero(source_match & target_matches[target_code], axis=1)
            )

    if count:
        # R's PA == 0 guard zeroes the complete row before its count branch.
        # This differs from returning the raw numerator when the caller selects
        # a missing/void target that is excluded from the denominator.
        numerators[denominators == 0, :] = 0
        return numerators
    return np.divide(
        numerators,
        denominators[:, None],
        out=np.zeros_like(numerators),
        where=denominators[:, None] != 0,
    )


def seqtrate(
    x: STSList,
    *,
    sel_states: Sequence[str] | None = None,
    time_varying: bool = False,
    weighted: bool = True,
    lag: int = 1,
    include_missing: bool = False,
    count: bool = False,
    statl: Sequence[str] | None = None,
    as_dataframe: bool = True,
) -> pd.DataFrame | np.ndarray | SeqTrateResult:
    """Transition probabilities or counts (TraMineR ``seqtrate``).

    ``sel_states`` preserves the caller's order.  Stationary output is K by K;
    time-varying output is K by K by the number of origin positions.  R's
    denominator includes every eligible destination (including states omitted
    from ``sel_states``), excludes structural void always, and excludes missing
    destinations unless ``include_missing``.  Consequently a selected subset's
    probability rows need not sum to one.  Empty source rows are exactly zero.

    ``statl`` is R's deprecated spelling for ``sel_states``.  Negative lag is
    supported for stationary output; TraMineR 2.2-14 refuses it for time-varying
    output, for which Python raises a clearer ``ValueError``.
    """
    if not isinstance(x, STSList):
        raise TypeError("x must be an STSList state sequence object")
    # Both marker representations are equivalent inputs (TX-AUD-20261010-094).
    x = x._canonicalized()
    if statl is not None:
        if sel_states is not None:
            raise ValueError("statl and sel_states cannot be supplied together")
        warnings.warn(
            "statl is deprecated; use sel_states instead",
            DeprecationWarning,
            stacklevel=_caller_stacklevel(),
        )
        sel_states = statl
    if not isinstance(lag, (int, np.integer)) or isinstance(lag, (bool, np.bool_)):
        raise TypeError("lag must be an integer")
    lag = int(lag)
    if time_varying and lag < 0:
        raise ValueError("negative lag is not supported when time_varying=True")
    if x.seqlen == 0 or abs(lag) >= x.seqlen:
        raise ValueError("absolute lag must be smaller than the sequence width")

    labels, state_codes = _seqtrate_states(x, sel_states, include_missing)
    if not labels:
        raise ValueError("sel_states must contain at least one state")
    if weighted and x.weights is not None:
        weights = np.asarray(x.weights, dtype=float)
        if weights.ndim != 1 or weights.shape[0] != x.nseq:
            raise ValueError(f"sequence weights must have length {x.nseq}")
    else:
        weights = np.ones(x.nseq, dtype=float)

    if lag < 0:
        origins: np.ndarray = np.arange(-lag, x.seqlen, dtype=int)
    else:
        origins = np.arange(0, x.seqlen - lag, dtype=int)
    destinations: np.ndarray = origins + lag

    if not time_varying:
        result = _seqtrate_accumulate(
            x,
            state_codes,
            origins,
            destinations,
            weights,
            include_missing=include_missing,
            count=count,
        )
        if not as_dataframe:
            return result
        return pd.DataFrame(result, index=labels, columns=labels)

    values: np.ndarray = np.empty((len(labels), len(labels), len(origins)), dtype=float)
    for position, (origin, destination) in enumerate(zip(origins, destinations)):
        values[:, :, position] = _seqtrate_accumulate(
            x,
            state_codes,
            np.array([origin], dtype=int),
            np.array([destination], dtype=int),
            weights,
            include_missing=include_missing,
            count=count,
        )
    if not as_dataframe:
        return values

    if x.xtlab is not None and len(x.xtlab) >= x.seqlen:
        all_time_labels = [str(label) for label in x.xtlab]
    else:
        all_time_labels = [f"T{x.start + offset}" for offset in range(x.seqlen)]
    time_labels = tuple(all_time_labels[int(origin)] for origin in origins)
    return SeqTrateResult(
        values=values,
        from_states=tuple(labels),
        to_states=tuple(labels),
        time_labels=time_labels,
    )


@_require_stslist
def seqhasmiss(x: STSList, per_sequence: bool = False):
    """Check if sequences contain missing (nr) or void codes."""
    missing_codes = []
    for code in (x.nr_code, x.void_code):
        if code is not None:
            missing_codes.append(int(code))
    if not missing_codes:
        return np.zeros(x.nseq, dtype=bool) if per_sequence else False
    arr = np.isin(x.data, np.array(missing_codes, dtype=int))
    if per_sequence:
        return arr.any(axis=1)
    return bool(arr.any())
