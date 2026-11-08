from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from ..sequences import STSList, _validate_sequence_weights
from ._validation import check_costs_requested, validate_smooth_window
from .python_kernels import _validate_dense_sequence_matrix, _warn_cost_quality

_SEQCOST_METHODS = ("CONSTANT", "TRATE", "FEATURES", "FUTURE", "INDELS", "INDELSLOG")
_SEQCOST_TRANSITIONS = ("previous", "next", "both")


@dataclass(frozen=True, eq=False)
class SeqCostMatrix:
    """NumPy-compatible labelled K by K by time substitution-cost array."""

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


def _seqcost_source_with_missing_state(seqdata: STSList, with_missing: bool) -> STSList:
    """Materialise the ``nr`` state in the codebook when ``with_missing`` asks for it.

    R's ``seqdata`` always carries an ``nr`` attribute, so
    ``alphabet(seqdata, with.missing = TRUE)`` can append it whether or not any
    cell is actually missing.  An ``STSList`` only allocates ``nr_code`` once a
    missing cell exists, so the same request needs the codebook widened first.
    The extra row/column is never reached by the source's own cells -- it exists
    so that a caller's ``refseq`` (or any later widening) can index it.

    Returns ``seqdata`` unchanged whenever nothing needs to be added.
    """
    if not with_missing or seqdata.nr_code is not None:
        return seqdata
    widened = seqdata.copy()
    widened.nr_code = len(widened.alphabet)
    widened.alphabet.append(widened.nr)
    if widened.labels is not None:
        widened.labels.append(widened.nr)
    if widened.cpal is not None:
        widened.cpal.append(widened.missing_color)
    return widened


def _seqcost_states(seqdata: STSList, with_missing: bool) -> tuple[list[str], list[int], bool]:
    """Return the cost-matrix states, honouring ``with_missing`` as 2.2-14 does.

    TraMineR 2.2-13 downgraded ``with.missing`` to ``FALSE`` here whenever the
    *source* carried no missing cell.  ``seqdist.R:112-113`` had already decided
    the flag across both the source and an external ``refseq``, so the downgrade
    silently undid a correct decision and returned an undersized ``sm`` -- an
    out-of-bounds read in R.  TraMineR 2.2-14 removed the block (the Peyrot
    issue); this follows it.  Pinned by
    ``tests/test_seqcost_refseq_missing_contract.py``.
    """
    marker_codes = {code for code in (seqdata.nr_code, seqdata.void_code) if code is not None}
    labels = [str(label) for code, label in enumerate(seqdata.alphabet) if code not in marker_codes]
    codes = [code for code in range(len(seqdata.alphabet)) if code not in marker_codes]

    nr_code = seqdata.nr_code
    if with_missing and nr_code is not None:
        labels.append(str(seqdata.nr))
        codes.append(int(nr_code))
    return labels, codes, bool(with_missing)


def _seqcost_time_labels(seqdata: STSList) -> list[str]:
    if seqdata.xtlab is not None and len(seqdata.xtlab) >= seqdata.seqlen:
        return [str(label) for label in seqdata.xtlab[: seqdata.seqlen]]
    return [f"T{seqdata.start + offset}" for offset in range(seqdata.seqlen)]


def _resolve_feature_columns(
    columns: Sequence[Any], selection: Any, *, component: str
) -> list[int]:
    selected = [selection] if isinstance(selection, (str, int, np.integer)) else list(selection)
    resolved: list[int] = []
    for value in selected:
        if isinstance(value, str):
            if value not in columns:
                raise ValueError(f"feature_type[{component!r}] has invalid column {value!r}")
            resolved.append(list(columns).index(value))
        elif isinstance(value, (int, np.integer)) and not isinstance(value, (bool, np.bool_)):
            # Match R's named daisy type list: numeric column references are 1-based.
            index = int(value) - 1
            if index < 0 or index >= len(columns):
                raise ValueError(
                    f"feature_type[{component!r}] column numbers must be in 1..{len(columns)}"
                )
            resolved.append(index)
        else:
            raise ValueError(
                f"feature_type[{component!r}] must contain column names or 1-based numbers"
            )
    return resolved


def _seqcost_gower(
    state_features: pd.DataFrame | np.ndarray,
    *,
    feature_weights: float | Sequence[float] | None,
    feature_type: Mapping[str, Any] | None,
) -> np.ndarray:
    """The cluster::daisy(metric='gower') path required by R ``seqcost``."""
    if isinstance(state_features, pd.DataFrame):
        frame = state_features.copy(deep=False)
    else:
        array = np.asarray(state_features)
        if array.ndim != 2:
            raise ValueError("state_features must be a 2D table")
        if not np.issubdtype(array.dtype, np.number):
            raise ValueError("non-numeric state_features must be supplied as a pandas DataFrame")
        frame = pd.DataFrame(array)

    nrows, nfeatures = frame.shape
    columns = list(frame.columns)
    raw_columns: list[np.ndarray] = []
    feature_types: list[str] = []
    for column in columns:
        series = frame[column]
        if isinstance(series.dtype, pd.CategoricalDtype):
            categorical = series.array
            values = categorical.codes.astype(float) + 1.0
            values[categorical.codes < 0] = np.nan
            raw_columns.append(values)
            feature_types.append("O" if series.dtype.ordered else "N")
        elif pd.api.types.is_bool_dtype(series.dtype):
            raw_columns.append(series.astype(float).to_numpy())
            feature_types.append("A")
        elif pd.api.types.is_numeric_dtype(series.dtype):
            raw_columns.append(pd.to_numeric(series, errors="coerce").to_numpy(dtype=float))
            feature_types.append("I")
        else:
            raise ValueError(
                f"invalid feature type for column {column!r}; use numeric or categorical data"
            )

    if feature_type is None:
        feature_type = {}
    if not isinstance(feature_type, Mapping):
        raise ValueError("feature_type must be a named mapping from daisy types to columns")
    allowed_components = {
        "asymm",
        "symm",
        "factor",
        "ordered",
        "logratio",
        "ordratio",
        "numeric",
        "integer",
    }
    unknown_components = sorted(set(feature_type) - allowed_components)
    if unknown_components:
        raise ValueError(f"unknown feature_type components: {unknown_components}")

    for component, selection in feature_type.items():
        indices = _resolve_feature_columns(columns, selection, component=component)
        for index in indices:
            values = raw_columns[index]
            if component in {"asymm", "symm"}:
                observed = values[~pd.isna(values)]
                unique = np.unique(observed)
                if unique.size > 2 or not np.all(np.isin(unique, [0.0, 1.0])):
                    raise ValueError(
                        "at least one binary feature has values outside {0, 1, missing}"
                    )
                feature_types[index] = "A" if component == "asymm" else "S"
            elif component == "ordratio":
                observed = values[~pd.isna(values)]
                levels = sorted(np.unique(observed).tolist())
                level_to_code = {level: code + 1.0 for code, level in enumerate(levels)}
                raw_columns[index] = np.array(
                    [np.nan if pd.isna(value) else level_to_code[value] for value in values],
                    dtype=float,
                )
                feature_types[index] = "T"
            elif component == "logratio":
                with np.errstate(divide="ignore", invalid="ignore"):
                    raw_columns[index] = np.log10(values.astype(float))
                feature_types[index] = "I"
            # daisy validates the remaining named components but relies on the
            # input column's own data class rather than coercing it here.

    if feature_weights is None:
        weights = np.ones(nfeatures, dtype=float)
    else:
        weights = np.asarray(feature_weights, dtype=float)
        if weights.ndim == 0 or weights.size == 1:
            weights = np.full(nfeatures, float(weights.reshape(-1)[0]), dtype=float)
        elif weights.ndim != 1 or weights.size != nfeatures:
            raise ValueError("feature weights must have length equal to the number of columns")

    numerator = np.zeros((nrows, nrows), dtype=float)
    denominator = np.zeros((nrows, nrows), dtype=float)
    for index, kind in enumerate(feature_types):
        values = np.asarray(raw_columns[index], dtype=float)
        missing = pd.isna(values)
        comparable = (~missing)[:, None] & (~missing)[None, :]

        if kind in {"N", "S", "A"}:
            contribution = (values[:, None] != values[None, :]).astype(float)
        else:
            observed = values[~missing]
            if observed.size:
                span = float(np.max(observed) - np.min(observed))
            else:
                span = 0.0
            # daisy replaces a zero range by one; the feature still contributes
            # its weight to the denominator with an all-zero distance.
            if span == 0:
                span = 1.0
            contribution = np.abs(values[:, None] - values[None, :]) / span

        if kind == "A":
            comparable &= ~((values[:, None] == 0) & (values[None, :] == 0))
        contribution = np.where(comparable, contribution, 0.0)
        delta = comparable.astype(float)
        numerator += weights[index] * delta * contribution
        denominator += weights[index] * delta

    result = np.divide(
        numerator,
        denominator,
        out=np.full_like(numerator, np.nan),
        where=denominator != 0,
    )
    np.fill_diagonal(result, 0.0)
    return result


def _seqcost_indels(
    seqdata: STSList,
    state_codes: Sequence[int],
    *,
    method: str,
    with_missing: bool,
    weighted: bool,
    time_varying: bool,
) -> tuple[np.ndarray, np.ndarray]:
    if time_varying:
        # TraMineR 2.2-14 deliberately omits weighted= here.  seqstatd therefore
        # uses attached sequence weights even when seqcost(weighted=False).
        from ..stats_basic import seqstatd

        distribution = seqstatd(seqdata, include_missing=with_missing)
        frequencies = np.asarray(distribution["Frequencies"], dtype=float)[:, state_codes].T
        frequencies[:, np.asarray(distribution["ValidStates"]) == 0] = np.nan
    else:
        if weighted and seqdata.weights is not None:
            sequence_weights = np.asarray(seqdata.weights, dtype=float)
        else:
            sequence_weights = np.ones(seqdata.nseq, dtype=float)
        if sequence_weights.ndim != 1 or sequence_weights.shape[0] != seqdata.nseq:
            raise ValueError(f"sequence weights must have length {seqdata.nseq}")
        all_counts = np.zeros(len(seqdata.alphabet), dtype=float)
        for row in range(seqdata.nseq):
            all_counts += sequence_weights[row] * np.bincount(
                seqdata.data[row], minlength=len(seqdata.alphabet)
            )
        total = float(all_counts.sum())
        if total == 0:
            frequencies = np.full(len(state_codes), np.nan, dtype=float)
        else:
            frequencies = all_counts[list(state_codes)] / total

    frequencies = np.where(np.isnan(frequencies), 1.0, frequencies)
    if method == "INDELSLOG":
        indels = np.log(2.0 / (1.0 + frequencies))
    else:
        with np.errstate(divide="ignore"):
            indels = 1.0 / frequencies
        indels = np.where(np.isinf(indels), np.finfo(float).max, indels)

    if time_varying:
        costs = np.zeros((len(state_codes), len(state_codes), seqdata.seqlen), dtype=float)
        with np.errstate(over="ignore"):
            for time in range(seqdata.seqlen):
                costs[:, :, time] = np.add.outer(indels[:, time], indels[:, time])
                np.fill_diagonal(costs[:, :, time], 0.0)
    else:
        with np.errstate(over="ignore"):
            costs = np.add.outer(indels, indels)
        np.fill_diagonal(costs, 0.0)
    return indels, costs


def _seqcost_trate_time(
    seqdata: STSList,
    *,
    cval: float,
    with_missing: bool,
    weighted: bool,
    transition: str,
    lag: int,
) -> np.ndarray:
    if lag < 0 or lag > 1:
        raise ValueError(
            "TraMineR 2.2-14 refuses time-varying TRATE seqcost when lag is negative or greater than 1"
        )
    from ..stats_utils import seqtrate

    rates = np.asarray(
        seqtrate(
            seqdata,
            time_varying=True,
            weighted=weighted,
            lag=lag,
            include_missing=with_missing,
            as_dataframe=False,
        ),
        dtype=float,
    )
    nstates = rates.shape[0]
    if nstates < 2:
        raise ValueError("time-varying TRATE seqcost requires at least two states")
    costs = np.zeros((nstates, nstates, seqdata.seqlen), dtype=float)
    symmetric_rates = rates + rates.swapaxes(0, 1)
    for time in range(seqdata.seqlen):
        if transition == "previous":
            slice_cost = cval if time == 0 else cval - symmetric_rates[:, :, time - 1]
        elif transition == "next":
            slice_cost = cval if time == seqdata.seqlen - 1 else cval - symmetric_rates[:, :, time]
        else:
            at_start = time == 0
            at_end = time == seqdata.seqlen - 1
            if at_start and at_end:
                slice_cost = np.full((nstates, nstates), cval, dtype=float)
            elif at_start:
                slice_cost = cval - 2.0 * symmetric_rates[:, :, time]
            elif at_end:
                slice_cost = cval - 2.0 * symmetric_rates[:, :, time - 1]
            else:
                slice_cost = cval - symmetric_rates[:, :, time - 1] - symmetric_rates[:, :, time]
        costs[:, :, time] = np.maximum(0.0, slice_cost)
        np.fill_diagonal(costs[:, :, time], 0.0)
    return costs


def _apply_fixed_missing_cost(
    costs: np.ndarray, miss_cost: float | Sequence[float], nstates: int
) -> None:
    values = np.asarray(miss_cost, dtype=float)
    if values.ndim == 0 or values.size == 1:
        values = np.full(nstates - 1, float(values.reshape(-1)[0]), dtype=float)
    elif values.ndim != 1 or values.size != nstates - 1:
        raise ValueError(f"miss_cost must be scalar or have length {nstates - 1}")
    if costs.ndim == 2:
        costs[-1, :-1] = values
        costs[:-1, -1] = values
    else:
        costs[-1, :-1, :] = values[:, None]
        costs[:-1, -1, :] = values[:, None]


def seqcost(
    seqdata: STSList,
    method: str,
    cval: float | None = None,
    *,
    with_missing: bool = False,
    miss_cost: float | Sequence[float] | None = None,
    time_varying: bool = False,
    weighted: bool = True,
    transition: str = "both",
    lag: int = 1,
    miss_cost_fixed: bool | None = None,
    state_features: pd.DataFrame | np.ndarray | None = None,
    feature_weights: float | Sequence[float] | None = None,
    feature_type: Mapping[str, Any] | None = None,
    proximities: bool = False,
    check_costs: bool = True,
) -> dict[str, object]:
    """Substitution/indel costs with the supported TraMineR ``seqcost`` controls.

    ``check_costs`` is R's ``check.costs``: ``True`` (R's default) is accepted,
    and traminex always checks the costs; ``False`` is refused, because R then
    returns negative distances for negative costs (TX-AUD-20261011-041).
    """
    check_costs_requested(check_costs)
    if not isinstance(seqdata, STSList):
        raise TypeError("seqdata must be an STSList state sequence object")
    if method not in _SEQCOST_METHODS:
        raise ValueError(f"method must be one of: {' '.join(_SEQCOST_METHODS)}")
    if transition not in _SEQCOST_TRANSITIONS:
        raise ValueError(f"transition must be one of: {' '.join(_SEQCOST_TRANSITIONS)}")
    if not isinstance(lag, (int, np.integer)) or isinstance(lag, (bool, np.bool_)):
        raise TypeError("lag must be an integer")
    lag = int(lag)

    # Both marker representations are the same sequences; the cost code below
    # indexes the codebook, which only the canonical form lays out as R does
    # (audit S1a-02).
    seqdata = seqdata._canonicalized()
    seqdata = _seqcost_source_with_missing_state(seqdata, with_missing)
    state_labels, state_codes, with_missing = _seqcost_states(seqdata, with_missing)
    nstates = len(state_labels)
    if nstates < 2 and (method == "TRATE" or (method == "CONSTANT" and time_varying)):
        raise ValueError(f"{method} costs require at least two states in this mode")
    default_cval = 4.0 if time_varying and method == "TRATE" and transition == "both" else 2.0
    resolved_cval = default_cval if cval is None else float(cval)
    resolved_miss_cost: float | Sequence[float] = resolved_cval if miss_cost is None else miss_cost
    if miss_cost_fixed is None:
        miss_cost_fixed = method not in {"INDELS", "INDELSLOG"}

    indel: float | np.ndarray = 1.0
    if method == "CONSTANT":
        if np.isnan(resolved_cval):
            raise ValueError("no value for the constant substitution-cost")
        if time_varying:
            costs = np.full((nstates, nstates, seqdata.seqlen), resolved_cval, dtype=float)
            for time in range(seqdata.seqlen):
                np.fill_diagonal(costs[:, :, time], 0.0)
        else:
            costs = np.full((nstates, nstates), resolved_cval, dtype=float)
            np.fill_diagonal(costs, 0.0)
    elif method == "FUTURE":
        if time_varying:
            raise ValueError("time_varying substitution cost is not implemented for method FUTURE")
        from ..stats_utils import seqtrate

        rates = np.asarray(
            seqtrate(
                seqdata,
                weighted=weighted,
                lag=lag,
                include_missing=with_missing,
                as_dataframe=False,
            ),
            dtype=float,
        )
        column_sums = rates.sum(axis=0)
        inverse = np.divide(
            1.0,
            column_sums,
            out=np.zeros_like(column_sums),
            where=column_sums != 0,
        )
        differences = rates[:, None, :] - rates[None, :, :]
        costs = np.sqrt(np.einsum("k,ijk->ij", inverse, differences**2))
        np.fill_diagonal(costs, 0.0)
        indel = 0.5 * float(np.max(costs))
    elif method == "FEATURES":
        if time_varying:
            raise ValueError(
                "time_varying substitution cost is not implemented for method FEATURES"
            )
        feature_shape = () if state_features is None else np.shape(state_features)
        if len(feature_shape) != 2 or feature_shape[0] != nstates:
            raise ValueError("state_features must contain exactly one row per state")
        if feature_shape[1] == 0:
            raise ValueError("state_features must contain at least one feature column")
        costs = _seqcost_gower(
            state_features,
            feature_weights=feature_weights,
            feature_type=feature_type,
        )
        indel = 0.5 * float(np.max(costs))
    elif method in {"INDELS", "INDELSLOG"}:
        indel, costs = _seqcost_indels(
            seqdata,
            state_codes,
            method=method,
            with_missing=with_missing,
            weighted=weighted,
            time_varying=time_varying,
        )
    elif time_varying:
        costs = _seqcost_trate_time(
            seqdata,
            cval=resolved_cval,
            with_missing=with_missing,
            weighted=weighted,
            transition=transition,
            lag=lag,
        )
    else:
        from ..stats_utils import seqtrate

        rates = np.asarray(
            seqtrate(
                seqdata,
                weighted=weighted,
                lag=lag,
                include_missing=with_missing,
                as_dataframe=False,
            ),
            dtype=float,
        )
        costs = resolved_cval - rates - rates.T
        np.fill_diagonal(costs, 0.0)
        indel = 0.5 * float(np.max(costs))

    if with_missing and miss_cost_fixed:
        _apply_fixed_missing_cost(costs, resolved_miss_cost, nstates)

    # R's seqcost runs checkcost on every matrix it generates, after the fixed
    # missing cost and before any conversion to proximities (seqcost.R:286;
    # TX-AUD-20261011-080). Its time-varying array is (K, K, T); the check takes
    # time first.
    _warn_cost_quality(
        np.moveaxis(costs, -1, 0) if time_varying else costs,
        indel=None if time_varying else indel,
        labels=state_labels,
        stacklevel=3,
    )

    component = "prox" if proximities else "sm"
    if proximities:
        with np.errstate(divide="ignore", invalid="ignore"):
            values = 1.0 - costs / np.max(costs)
    else:
        values = costs

    time_labels = _seqcost_time_labels(seqdata)
    if time_varying:
        labelled_values: object = SeqCostMatrix(
            values=np.asarray(values, dtype=float),
            from_states=tuple(state_labels),
            to_states=tuple(state_labels),
            time_labels=tuple(time_labels),
        )
    else:
        labelled_values = pd.DataFrame(values, index=state_labels, columns=state_labels)

    if isinstance(indel, np.ndarray):
        if indel.ndim == 2:
            labelled_indel: object = pd.DataFrame(indel, index=state_labels, columns=time_labels)
        else:
            labelled_indel = pd.Series(indel, index=state_labels, name="indel")
    else:
        labelled_indel = float(indel)
    return {"indel": labelled_indel, component: labelled_values}


def trate_substitution_matrix(
    seqmat: np.ndarray,
    *,
    smoothing: float = 0.0,
    mode: str = "symmetric",
    normalize: str = "none",
    weights: np.ndarray | None = None,
) -> np.ndarray:
    """Compute a transition-rate-based substitution matrix (TRATE family).

    With the defaults (``smoothing=0``, ``mode="symmetric"``, ``normalize="none"``)
    this follows the formula of R ``TraMineR::seqcost(method="TRATE")``:
    ``cost[i,j] = 2 - tr[i,j] - tr[j,i]`` where ``tr`` is the (optionally weighted)
    row-normalized transition-rate matrix of ``seqtrate``. As in R's ``seqtrate``,
    a *source state with no observed transitions* (``PA == 0``) yields a **zero**
    transition-rate row, not a uniform ``1/K`` row.

    This raw-array extension infers the state count from observed codes. It
    also applies that formula to singleton trajectories without transitions;
    the pinned R singleton TRATE call itself raises an error.

    - mode options:
      - "symmetric" (default): cost[i,j] = 2 - (T[i,j] + T[j,i]) for i != j; 0 on diag.
      - "row": cost[i,j] = 1 - T[i,j] for i != j; 0 on diag.
      - "col": cost[i,j] = 1 - T[j,i] for i != j; 0 on diag.
    - normalize:
      - "none" (default): return raw costs.
      - "mean_offdiag": divide by mean(off-diagonal) to scale typical off-diagonal to 1.
      - "clip_0_2": clip final values to [0,2]. Can be combined by chaining if needed.
    - smoothing: optional Laplace-style additive smoothing of transition counts.
      Defaults to 0 (R-faithful). A positive value regularizes only rows that have
      at least one observed transition; empty rows stay zero (R's ``PA == 0``).
    - weights: optional per-sequence weights for transition counting.
    """
    seqmat = _validate_dense_sequence_matrix(seqmat).astype(np.int64, copy=False)
    validated_weights = (
        None if weights is None else _validate_sequence_weights(weights, seqmat.shape[0])
    )
    K = int(seqmat.max()) + 1 if seqmat.size else 0
    a = seqmat[:, :-1]
    b = seqmat[:, 1:]
    if validated_weights is not None:
        # accumulate weighted transitions
        counts = np.zeros((K, K), dtype=float)
        for i in range(seqmat.shape[0]):
            idx = a[i] * K + b[i]
            np.add.at(counts, np.unravel_index(idx, (K, K)), validated_weights[i])
    else:
        idx = a.ravel() * K + b.ravel()
        counts = np.bincount(idx, minlength=K * K).astype(float).reshape(K, K)
    return _trate_costs_from_counts(counts, smoothing=smoothing, mode=mode, normalize=normalize)


def _trate_costs_from_counts(
    counts: np.ndarray,
    *,
    smoothing: float,
    mode: str,
    normalize: str,
) -> np.ndarray:
    """Convert adjacent-transition counts to the TRATE substitution matrix."""
    if not np.isfinite(smoothing) or smoothing < 0:
        raise ValueError("smoothing must be finite and non-negative")
    K = counts.shape[0]
    # Row-normalize. Match R's seqtrate: a source state with no observed
    # transitions (raw row sum 0, R's "PA == 0") yields a zero rate row, not a
    # uniform 1/K row. Smoothing (default 0 = R-faithful) regularizes only rows
    # that have at least one observed transition.
    raw_rowsum = counts.sum(axis=1, keepdims=True)
    denom = raw_rowsum + smoothing * K
    T = np.divide(
        counts + smoothing,
        denom,
        out=np.zeros_like(counts, dtype=float),
        where=raw_rowsum > 0,
    )
    # Vectorized cost computation
    if mode == "symmetric":
        C = 2.0 - (T + T.T)
        np.fill_diagonal(C, 0.0)
    elif mode == "row":
        C = 1.0 - T
        np.fill_diagonal(C, 0.0)
    elif mode == "col":
        C = 1.0 - T.T
        np.fill_diagonal(C, 0.0)
    else:
        raise ValueError(f"Unsupported TRATE mode: {mode}")
    if normalize == "mean_offdiag":
        off = C[~np.eye(K, dtype=bool)]
        m = off.mean() if off.size else 1.0
        if m > 0:
            C = C / m
    elif normalize == "none":
        pass
    elif normalize == "clip_0_2":
        C = np.clip(C, 0.0, 2.0)
    else:
        raise ValueError(f"Unsupported normalize: {normalize}")
    return C


def _trate_substitution_matrix_rows(
    rows: list[np.ndarray],
    *,
    n_states: int,
    smoothing: float = 0.0,
    mode: str = "symmetric",
    normalize: str = "none",
    weights: np.ndarray | None = None,
) -> np.ndarray:
    """TRATE costs from canonical variable-length rows (no structural void)."""
    K = int(n_states)
    if K < 0:
        raise ValueError("n_states must be non-negative")
    counts = np.zeros((K, K), dtype=float)
    if weights is None:
        row_weights = np.ones(len(rows), dtype=float)
    else:
        row_weights = _validate_sequence_weights(weights, len(rows))

    flat = counts.reshape(-1)
    for row, weight in zip(rows, row_weights):
        sequence = np.asarray(row, dtype=np.int64)
        if sequence.ndim != 1:
            raise ValueError("canonical sequence rows must be one-dimensional")
        if sequence.size == 0:
            continue
        if int(sequence.min()) < 0 or int(sequence.max()) >= K:
            raise ValueError("sequence state code is outside n_states")
        if sequence.size > 1:
            indices = sequence[:-1] * K + sequence[1:]
            np.add.at(flat, indices, float(weight))

    return _trate_costs_from_counts(counts, smoothing=smoothing, mode=mode, normalize=normalize)


def constant_substitution_matrix(seqmat: np.ndarray, cval: float = 2.0) -> np.ndarray:
    """Return a constant off-diagonal substitution matrix with 0 on the diagonal.

    The size K is inferred from seqmat. Off-diagonal entries are set to cval.
    """
    seqmat = _validate_dense_sequence_matrix(seqmat)
    K = int(seqmat.max()) + 1 if seqmat.size else 0
    if K <= 0:
        return np.zeros((0, 0), dtype=float)
    M = np.full((K, K), float(cval), dtype=float)
    np.fill_diagonal(M, 0.0)
    return M


def indelslog_costs(
    seqmat: np.ndarray,
    *,
    include_missing: bool = False,
    missing_codes: tuple[int, int] | None = None,
) -> tuple[float, np.ndarray]:
    """Compute INDELSLOG costs from global state frequencies (TraMineR seqcost INDELSLOG).

    Returns (indel, sm) where sm[i,j] = log(2/(1+p_i)) + log(2/(1+p_j)) for i!=j, and 0 on
    the diagonal. The scalar indel value is the mean of log(2/(1+p_i)) over included states,
    a scalar reduction of R's per-state INDELSLOG indel vector ``indels <- log(2/(1 + p))``
    (seqcost.R line 145).

    - include_missing: if True, include missing/void codes passed via missing_codes.
    - missing_codes: optional (nr_code, void_code) to explicitly include as valid states when include_missing=True.
    """
    seqmat = _validate_dense_sequence_matrix(seqmat).astype(np.int64, copy=False)
    if seqmat.size == 0:
        return 0.0, np.zeros((0, 0), dtype=float)
    K = int(seqmat.max()) + 1
    counts = np.bincount(seqmat.ravel(), minlength=K).astype(float)
    include_mask = np.ones(K, dtype=bool)
    if not include_missing and missing_codes is not None:
        for code in missing_codes:
            if code is not None and 0 <= int(code) < K:
                include_mask[int(code)] = False
    # Excluded marker states get a placeholder frequency; an absent real state
    # keeps p = 0, for which R's log(2/(1+p)) is exactly log 2 (audit S1a-10).
    eps = 1e-12
    p = np.divide(counts, counts.sum() if counts.sum() > 0 else 1.0)
    p = np.where(~include_mask, eps, p)
    # Per-state indel cost, matching R seqcost.R: indels <- log(2 / (1 + p))
    log_cost = np.log(2.0 / (1.0 + p))
    # build substitution matrix: sm[i,j] = log_cost[i] + log_cost[j]
    sm = np.add.outer(log_cost, log_cost)
    np.fill_diagonal(sm, 0.0)
    # scalar indel as the mean of the per-state log(2/(1+p)) costs
    indel = float(log_cost[include_mask].mean()) if include_mask.any() else float(log_cost.mean())
    return indel, sm


def indels_costs(
    seqmat: np.ndarray,
    *,
    include_missing: bool = False,
    missing_codes: tuple[int, int] | None = None,
) -> tuple[float, np.ndarray]:
    """Compute INDELS costs from global state frequencies (TraMineR seqcost INDELS).

    Returns (indel, sm) where sm[i,j] = (1/p_i) + (1/p_j) for i!=j, and 0 on the diagonal.
    The scalar indel value is the mean of 1/p_i over included states, a scalar reduction of
    R's per-state INDELS indel vector ``indels <- 1/p`` (seqcost.R line 147), with
    non-finite entries clamped to the largest finite float.

    - include_missing: if True, include missing/void codes passed via missing_codes.
    - missing_codes: optional (nr_code, void_code) to explicitly include as valid states when include_missing=True.
    """
    seqmat = _validate_dense_sequence_matrix(seqmat).astype(np.int64, copy=False)
    if seqmat.size == 0:
        return 0.0, np.zeros((0, 0), dtype=float)
    K = int(seqmat.max()) + 1
    counts = np.bincount(seqmat.ravel(), minlength=K).astype(float)
    include_mask = np.ones(K, dtype=bool)
    if not include_missing and missing_codes is not None:
        for code in missing_codes:
            if code is not None and 0 <= int(code) < K:
                include_mask[int(code)] = False
    # Excluded marker states get a placeholder frequency; an absent real state
    # keeps p = 0, so 1/p is Inf and is clamped to the largest double as R's
    # seqcost does -- the floor used to make it 1e12 (audit S1a-10).
    eps = 1e-12
    p = np.divide(counts, counts.sum() if counts.sum() > 0 else 1.0)
    p = np.where(~include_mask, eps, p)
    # Per-state indel cost, matching R seqcost.R: indels <- 1/p (inf clamped to max float)
    with np.errstate(divide="ignore"):
        inv = 1.0 / p
    inv = np.where(np.isinf(inv), np.finfo(float).max, inv)
    # build substitution matrix: sm[i,j] = inv[i] + inv[j]
    sm = np.add.outer(inv, inv)
    np.fill_diagonal(sm, 0.0)
    # scalar indel as the mean of the per-state 1/p costs
    indel = float(inv[include_mask].mean()) if include_mask.any() else float(inv.mean())
    return indel, sm


def gower_distance(
    features,
    *,
    weights=None,
    types=None,
    ordinal_levels=None,
) -> np.ndarray:
    """Gower mixed-type dissimilarity for the raw-array FEATURES extension.

    Computes the K x K Gower dissimilarity between rows of a per-object feature
    table. This extension drops constant numeric columns and returns zero for
    pairs without comparable positive-weight features. The separate STSList
    ``seqcost`` path follows TraMineR's ``cluster::daisy`` semantics, including
    constant columns in the denominator and returning NaN for incomparable pairs.

    d(i,j) = sum_f w_f * delta_ijf * d_ijf / sum_f w_f * delta_ijf

    where per feature f:
    - numeric/interval: d = |x_i - x_j| / range_f (zero-range columns are dropped);
      observed values must be finite.
    - ordinal: values mapped to level ranks scaled to [0,1], then treated as numeric.
    - nominal/factor: d = 0 if equal else 1.
    - asymmetric binary ("asymm"): d = 0/1, but joint absence (both 0) sets delta=0.
      Observed values must be 0 or 1.
    Missing values (NaN, None or pandas.NA) set delta=0 for affected pairs of
    every feature type and are excluded when inferring ordinal levels.

    Parameters
    ----------
    features : array-like
        2D table of shape ``(K, F)`` (object dtype allowed for mixed types).
    weights : sequence of float | None
        Finite, non-negative per-feature weights (default all 1), with at
        least one positive weight. Must be one-dimensional with length F.
    types : sequence of str | None
        Per-feature type in {"numeric"/"interval", "ordinal", "nominal"/"factor",
        "asymm"}. Defaults to all "numeric".
    ordinal_levels : sequence of sequence | None
        Ordered level lists for the ordinal columns (in ordinal-column order). If omitted,
        levels are inferred as the sorted unique values.

    Returns
    -------
    np.ndarray
        Symmetric ``(K, K)`` dissimilarity matrix in [0, 1] with zero diagonal.
    """
    arr = np.asarray(features, dtype=object)
    if arr.ndim != 2:
        raise ValueError("features must be 2D [K, F]")
    K, F = arr.shape
    if types is None:
        types = ["numeric"] * F
    if len(types) != F:
        raise ValueError("types length must equal number of feature columns")
    w = np.ones(F, dtype=float) if weights is None else np.asarray(weights, dtype=float)
    if w.ndim != 1 or w.size != F:
        raise ValueError("weights must be one-dimensional with length equal to feature columns")
    if not np.all(np.isfinite(w)) or np.any(w < 0) or (F > 0 and not np.any(w > 0)):
        raise ValueError("weights must be finite, non-negative, and have a positive total")

    if K == 0:
        return np.zeros((0, 0), dtype=float)

    # Double precision on every platform (TX-AUD-20261011-035). Overflow and lost
    # subnormal products are avoided by exact power-of-two scaling: each numeric
    # column is scaled so its largest magnitude lies in [0.5, 1) before its range
    # and differences are taken, and each pair's weights so that its largest
    # comparable weight does. A global rescaling could erase the only weight
    # comparable for a pair; what the per-pair scaling can push below the
    # subnormal range is negligible beside the pair's largest weight.
    terms = []  # (feature, comparable pairs, per-pair distance)
    ord_idx = 0
    with np.errstate(under="ignore"):
        for f in range(F):
            t = str(types[f]).lower()
            col = arr[:, f]
            miss = np.asarray(pd.isna(col), dtype=bool)
            comparable = (~miss)[:, None] & (~miss)[None, :]
            if t in ("numeric", "interval", "ordratio"):
                x = np.asarray(np.where(miss, np.nan, col), dtype=float)
                miss = np.isnan(x)
                comparable = (~miss)[:, None] & (~miss)[None, :]
                observed = x[~np.isnan(x)]
                if observed.size == 0:
                    continue
                if not np.all(np.isfinite(observed)):
                    raise ValueError("numeric Gower features must be finite or missing")
                # Subtract before normalization, preserving nearby large values;
                # the scaling keeps opposite-sign float64 extrema finite.
                exponent = int(np.frexp(np.max(np.abs(observed)))[1])
                x = np.ldexp(x, -exponent)
                observed = np.ldexp(observed, -exponent)
                rng = np.max(observed) - np.min(observed)
                if rng <= 0:
                    continue  # zero range contributes nothing and is excluded from denom
                dmat = np.abs(x[:, None] - x[None, :]) / rng
            elif t == "ordinal":
                if ordinal_levels is not None:
                    levels = list(ordinal_levels[ord_idx])
                else:
                    levels = sorted(set(col[~miss].tolist()))
                ord_idx += 1
                if len(levels) <= 1:
                    continue
                level_map = {lv: i for i, lv in enumerate(levels)}
                ranks = np.array(
                    [np.nan if absent else level_map[v] for v, absent in zip(col, miss)],
                    dtype=float,
                )
                # Rank differences are exact; one division rounds once.
                dmat = np.abs(ranks[:, None] - ranks[None, :]) / (len(levels) - 1)
            elif t in ("nominal", "factor", "categorical"):
                values = col.copy()
                values[miss] = None
                eq = values[:, None] == values[None, :]
                dmat = np.where(eq, 0.0, 1.0)
            elif t in ("asymm", "asymmetric", "binary"):
                b = np.asarray(np.where(miss, np.nan, col), dtype=float)
                miss = np.isnan(b)
                comparable = (~miss)[:, None] & (~miss)[None, :]
                if not np.all(np.isin(b[~miss], [0, 1])):
                    raise ValueError("asymmetric binary Gower features must be 0, 1 or missing")
                dmat = np.where(b[:, None] == b[None, :], 0.0, 1.0)
                both_zero = (b[:, None] == 0) & (b[None, :] == 0)
                comparable &= ~both_zero
            else:
                raise ValueError(f"Unknown Gower feature type: {types[f]!r}")
            terms.append((f, comparable, np.where(comparable, dmat, 0.0)))

        largest = np.zeros((K, K))
        for f, comparable, _ in terms:
            np.maximum(largest, np.where(comparable, w[f], 0.0), out=largest)
        pair_exponent = np.frexp(largest)[1]
        num = np.zeros((K, K))  # accumulated w_f * delta * d, scaled per pair
        den = np.zeros((K, K))  # accumulated w_f * delta, scaled per pair
        for f, comparable, dmat in terms:
            scaled = np.zeros((K, K))
            scaled[comparable] = np.ldexp(w[f], -pair_exponent[comparable])
            num += scaled * dmat
            den += scaled
        D = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
    np.fill_diagonal(D, 0.0)
    return D


def future_costs(
    seqmat: np.ndarray,
    *,
    weights: np.ndarray | None = None,
) -> tuple[float, np.ndarray]:
    """FUTURE substitution costs from common future (TraMineR seqcost.R FUTURE).

    Builds the row-normalized transition-rate matrix T (``seqtrate``), then applies
    the column-weighted Euclidean distance between rows (``chisqdista`` in seqcost.R):

        cost[i,j] = sqrt( sum_k Pdot_k * (T[i,k] - T[j,k])^2 ),  Pdot_k = 1 / colSum_k

    with ``Pdot_k = 0`` where a column sum is zero. Returns ``(indel, sm)`` with
    ``indel = 0.5 * max(sm)`` (seqcost.R line 105) and zero diagonal.
    """
    seqmat = _validate_dense_sequence_matrix(seqmat)
    validated_weights = (
        None if weights is None else _validate_sequence_weights(weights, seqmat.shape[0])
    )
    K = int(seqmat.max()) + 1 if seqmat.size else 0
    if K <= 0 or seqmat.shape[1] < 2:
        return 0.0, np.zeros((K, K), dtype=float)
    a = seqmat[:, :-1]
    b = seqmat[:, 1:]
    counts = np.zeros((K, K), dtype=float)
    if validated_weights is not None:
        for i in range(seqmat.shape[0]):
            np.add.at(counts, (a[i], b[i]), validated_weights[i])
    else:
        np.add.at(counts, (a.ravel(), b.ravel()), 1.0)
    row_sums = counts.sum(axis=1, keepdims=True)
    T = np.divide(counts, row_sums, out=np.zeros_like(counts), where=row_sums > 0)
    col_sums = T.sum(axis=0)
    # Pdot = 1 / colSums, with zero columns contributing 0 (R: cs[cs==0] <- Inf)
    pdot = np.where(col_sums == 0, 0.0, 1.0 / np.where(col_sums == 0, 1.0, col_sums))
    diff = T[:, None, :] - T[None, :, :]  # [K,K,K]
    sm = np.sqrt(np.einsum("k,ijk->ij", pdot, diff**2))
    np.fill_diagonal(sm, 0.0)
    indel = 0.5 * float(sm.max()) if sm.size else 0.0
    return indel, sm


def trate_substitution_matrix_timevarying(
    seqmat: np.ndarray,
    *,
    smoothing: float = 0.0,
    mode: str = "symmetric",
    normalize: str = "none",
    weights: np.ndarray | None = None,
    lag: int = 1,
) -> np.ndarray:
    """Time-varying TRATE substitution costs per position: returns array of shape [L, K, K].

    mode: 'symmetric'|'row'|'col' applies per time slice. weights are per-sequence.
    lag currently fixed to 1 for adjacent transitions. As in the time-constant
    helper and R's ``seqtrate``, a source state with no observed transitions at a
    position yields a zero rate row; ``smoothing`` defaults to 0 (R-faithful).
    """
    seqmat = _validate_dense_sequence_matrix(seqmat).astype(np.int64, copy=False)
    if not np.isfinite(smoothing) or smoothing < 0:
        raise ValueError("smoothing must be finite and non-negative")
    validated_weights = (
        None if weights is None else _validate_sequence_weights(weights, seqmat.shape[0])
    )
    n, L = seqmat.shape
    K = int(seqmat.max()) + 1 if seqmat.size else 0
    if lag != 1:
        raise ValueError(f"lag must be 1 for time-varying TRATE (got lag={lag})")
    if K == 0 or L == 0:
        # Empty data retain their time axis and still enforce the same option
        # vocabulary as a nonempty corpus.
        _trate_costs_from_counts(
            np.zeros((K, K)), smoothing=smoothing, mode=mode, normalize=normalize
        )
        return np.zeros((L, K, K), dtype=float)
    costs = np.zeros((L, K, K), dtype=float)
    for t in range(L):
        # build transitions impacting position t: previous and next
        counts = np.zeros((K, K), dtype=float)
        if t > 0:
            a = seqmat[:, t - 1]
            b = seqmat[:, t]
            if validated_weights is not None:
                idx = a * K + b
                np.add.at(counts, np.unravel_index(idx, (K, K)), validated_weights)
            else:
                idx = a * K + b
                np.add.at(counts, np.unravel_index(idx, (K, K)), 1.0)
        if t < L - 1:
            a = seqmat[:, t]
            b = seqmat[:, t + 1]
            if validated_weights is not None:
                idx = a * K + b
                np.add.at(counts, np.unravel_index(idx, (K, K)), validated_weights)
            else:
                idx = a * K + b
                np.add.at(counts, np.unravel_index(idx, (K, K)), 1.0)
        raw_rowsum = counts.sum(axis=1, keepdims=True)
        denom = raw_rowsum + smoothing * K
        T = np.divide(
            counts + smoothing,
            denom,
            out=np.zeros_like(counts, dtype=float),
            where=raw_rowsum > 0,
        )
        # Vectorized cost computation
        if mode == "symmetric":
            C = 2.0 - (T + T.T)
            np.fill_diagonal(C, 0.0)
        elif mode == "row":
            C = 1.0 - T
            np.fill_diagonal(C, 0.0)
        elif mode == "col":
            C = 1.0 - T.T
            np.fill_diagonal(C, 0.0)
        else:
            raise ValueError(f"Unsupported TRATE mode: {mode}")
        if normalize == "mean_offdiag":
            off = C[~np.eye(K, dtype=bool)]
            m = off.mean() if off.size else 1.0
            if m > 0:
                C = C / m
        elif normalize == "clip_0_2":
            C = np.clip(C, 0.0, 2.0)
        elif normalize == "none":
            pass
        else:
            raise ValueError(f"Unsupported normalize: {normalize}")
        costs[t] = C
    return costs


def indel_from_subm(subm: np.ndarray, *, strategy: str = "mean_offdiag_half") -> float:
    """Derive an indel cost from a substitution matrix.

    Strategies:
    - "mean_offdiag_half" (default): 0.5 * mean(off-diagonal costs)
    - "mean_offdiag": mean(off-diagonal costs)
    """
    M = np.asarray(subm, dtype=float)
    if M.ndim != 2 or M.shape[0] != M.shape[1]:
        raise ValueError("subm must be square [K,K]")
    K = M.shape[0]
    off = M[~np.eye(K, dtype=bool)]
    m = float(off.mean()) if off.size else 1.0
    if strategy == "mean_offdiag_half":
        return 0.5 * m
    if strategy == "mean_offdiag":
        return m
    raise ValueError(f"Unknown strategy: {strategy}")


def local_indel_by_pos(
    seqmat: np.ndarray,
    *,
    method: str = "change",
    base: float = 1.0,
    scale: float = 1.0,
    smooth_window: int | None = None,
) -> np.ndarray:
    """Compute per-position indel costs for OM (local variant).

    method options:
    - "change" (default):
      Estimate change rate at each position t as Pr[X_t != X_{t+1}] across sequences.
      indel[t] = clip(base * (1 - scale * change_rate[t]), min=1e-6)
    - "entropy":
      Compute normalized entropy H_t in [0,1] of the position-wise distribution; indel[t] = base * (1 - scale * H_t)
    If smooth_window is provided, apply moving average smoothing along positions.
    """
    seqmat = _validate_dense_sequence_matrix(seqmat)
    smooth_window = validate_smooth_window(smooth_window)
    n, L = seqmat.shape
    if L == 0:
        return np.array([], dtype=float)
    signal = np.zeros(L, dtype=float)
    if method == "change":
        if L > 1:
            diff = (seqmat[:, 1:] != seqmat[:, :-1]).astype(float)
            signal[:-1] = diff.mean(axis=0)
            signal[-1] = signal[-2] if L > 1 else 0.0
    elif method == "entropy":
        K = int(seqmat.max()) + 1 if seqmat.size else 0
        for t in range(L):
            vals, cnt = np.unique(seqmat[:, t], return_counts=True)
            s = cnt.sum()
            p = cnt / s if s > 0 else np.zeros_like(cnt, dtype=float)
            # Every share is positive (counts of observed codes): no constant in
            # the log (TX-AUD-20261010-095).
            H = float(-(p * np.log(p)).sum())
            Hn = H / (np.log(K) if K > 1 else 1.0)
            signal[t] = Hn
    else:
        raise ValueError("Unknown method for local_indel_by_pos: %r" % method)
    if smooth_window and smooth_window > 1:
        w = np.ones(smooth_window, dtype=float) / smooth_window
        left_pad = (smooth_window - 1) // 2
        right_pad = smooth_window // 2
        padded = np.pad(signal, (left_pad, right_pad), mode="edge")
        signal = np.convolve(padded, w, mode="valid")
    indel = base * (1.0 - scale * signal)
    indel = np.clip(indel, 1e-6, None)
    return indel
