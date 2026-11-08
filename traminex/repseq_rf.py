"""Relative-frequency groups derived from TraMineR 2.2-14 ``dissrf``/``seqrf``."""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Any, Dict, Literal, Optional, cast

import numpy as np
import pandas as pd
from scipy import linalg, stats

from ._categorical_order import r_order_keys
from ._rformat import r_number_text
from .diss import disscenter
from .sequences import STSList, _remap_issues


@dataclass(frozen=True)
class _BoundaryRecord:
    """One case's (possibly zero or fractional) allocation to an RF group."""

    case: int
    source_weight: float
    group: int
    assigned_weight: float


def _medoid_index(
    matrix: np.ndarray,
    *,
    weights: np.ndarray,
    squared: bool,
) -> int:
    """Return the first case at the minimum, as R's ``disscenter(medoids.index = "first")``.

    The contributions follow R's summation order (``_weighted_inertia_contrib``),
    so R's ties are exact ties and its first minimum is the exact first minimum
    (``disscenter.R:90-104``); no tolerance is needed (TX-AUD-20261011-023).
    """
    center_distances = np.asarray(
        disscenter(matrix, weights=weights, squared=squared),
        dtype=float,
    )
    minimum = float(np.min(center_distances))
    return int(np.flatnonzero(center_distances == minimum)[0])


def _symmetric_axis(gram: np.ndarray) -> tuple[np.ndarray, float]:
    """Return the leading classical-scaling eigenvector and eigenvalue."""
    values, vectors = linalg.eigh(
        gram,
        lower=True,
        driver="evr",
        check_finite=False,
    )
    leading = int(np.argmax(values))
    return vectors[:, leading], float(values[leading])


#: Relative tolerance under which two coordinate magnitudes count as tied.
_ORIENTATION_TIE = 1e-10


def _oriented(coordinates: np.ndarray) -> np.ndarray:
    """The axis with its coordinate of largest absolute value positive.

    A one-dimensional MDS axis has no sign of its own: R takes it from LAPACK,
    so its groups follow the numerical library (TX-REP-20261017-006).
    Coordinates within a relative ``_ORIENTATION_TIE`` of the largest magnitude
    are tied, and the first of them decides: on symmetric data two coordinates
    are equal in exact arithmetic and differ only in rounding, which would
    otherwise leave the choice to the build. Missing (zero-weight) coordinates
    are ignored, as R's ``which`` ignores NA.
    """
    magnitude = np.where(np.isfinite(coordinates), np.abs(coordinates), -np.inf)
    if not np.any(np.isfinite(magnitude)):
        return coordinates
    tied = magnitude >= np.max(magnitude) * (1.0 - _ORIENTATION_TIE)
    largest = int(np.flatnonzero(tied)[0])
    return -coordinates if coordinates[largest] < 0 else coordinates


def _cmdscale_1d(matrix: np.ndarray) -> np.ndarray:
    """One-dimensional ``stats::cmdscale`` coordinates, oriented by :func:`_oriented`."""
    squared = np.asarray(matrix, dtype=float) ** 2
    row_means = np.mean(squared, axis=1)
    centered = squared - row_means[:, None] - row_means[None, :] + float(np.mean(row_means))
    gram = -0.5 * centered
    spectrum = linalg.eigvalsh(gram, lower=True, driver="evr", check_finite=False)
    if spectrum.size > 1 and np.isclose(spectrum[-1], spectrum[-2], rtol=1e-12, atol=0.0):
        # A repeated leading root has no unique one-dimensional orientation.
        # Original order is the deterministic extension used by traminex. The
        # test is relative, as the axis is: R's cmdscale has none, and an
        # absolute tolerance made tiny distances (eigenvalues near 1e-20) count
        # as repeated roots (TX-AUD-20261011-052).
        return np.arange(matrix.shape[0], dtype=float)
    vector, value = _symmetric_axis(gram)
    return _oriented(vector * np.sqrt(max(value, 0.0)))


def _weighted_cmdscale_1d(diss: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """One-dimensional vegan ``wcmdscale`` (two weighted centering steps), oriented."""
    matrix = np.asarray(diss, dtype=float) ** 2
    resolved = np.asarray(weights, dtype=float)
    total = float(np.sum(resolved))
    roots = np.sqrt(resolved)

    first_means = np.sum(resolved[:, None] * matrix, axis=0) / total
    first = roots[:, None] * (matrix - first_means)
    transposed = first.T
    second_means = np.sum(resolved[:, None] * transposed, axis=0) / total
    weighted_centered = (roots[:, None] * (transposed - second_means)).T

    vector, value = _symmetric_axis(-0.5 * weighted_centered)
    with np.errstate(divide="ignore", invalid="ignore"):
        coordinates = vector * np.sqrt(abs(value)) / roots
    coordinates[~np.isfinite(coordinates)] = np.nan
    return _oriented(coordinates)


def _validate_dissimilarities(diss: np.ndarray) -> np.ndarray:
    matrix = np.asarray(diss, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("diss must be a square dissimilarity matrix")
    if matrix.shape[0] == 0:
        raise ValueError("diss must contain at least one case")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("diss must contain only finite values")
    return matrix


def _normalize_weights(ncase: int, weights: Optional[np.ndarray]) -> tuple[np.ndarray, bool]:
    if weights is None:
        return np.ones(ncase, dtype=float), False
    resolved = np.asarray(weights, dtype=float)
    if resolved.ndim != 1 or resolved.shape[0] != ncase:
        raise ValueError("weights must contain one value for each case")
    if not np.all(np.isfinite(resolved)) or np.any(resolved < 0):
        raise ValueError("weights must be finite and non-negative")
    total = float(np.sum(resolved))
    if total <= 0:
        raise ValueError("weights must have a positive sum")
    return ncase * resolved / total, True


def _resolve_axis(
    matrix: np.ndarray,
    sortv: object,
    *,
    weights: np.ndarray,
    weighted: bool,
    squared: bool,
) -> np.ndarray:
    ncase = matrix.shape[0]
    if sortv is None:
        return np.arange(ncase, dtype=float)
    if isinstance(sortv, str):
        if sortv != "mds":
            raise ValueError('the only scalar string allowed for sortv is "mds"')
        powered = matrix ** (2 if squared else 1)
        if weighted:
            return _weighted_cmdscale_1d(powered, weights)
        return _cmdscale_1d(powered)
    # R takes order(sortv) (dissrf.R:94): a factor by its level codes, text in R's
    # collation, missing values last (TX-REP-20261021-001).
    codes = getattr(sortv, "codes", None)
    if codes is None and isinstance(getattr(sortv, "dtype", None), pd.CategoricalDtype):
        codes = np.asarray(cast(pd.Series, sortv).cat.codes)  # a category-dtype Series
    values = np.asarray(codes if codes is not None else sortv)
    if values.ndim == 0 or values.size != ncase:
        raise ValueError("length of sortv must equal the number of cases")
    if codes is not None:
        keys = values.reshape(-1).astype(np.int64)
        keys[keys < 0] = np.iinfo(np.int64).max  # pandas codes a missing value -1
        return keys
    flat = values.reshape(-1).tolist()
    if any(isinstance(value, str) for value in flat):
        return cast(np.ndarray, r_order_keys(flat))
    try:
        numeric = values.reshape(-1).astype(float)
    except (TypeError, ValueError) as exc:
        raise ValueError("sortv must contain numeric values") from exc
    if not np.all(np.isfinite(numeric)):
        raise ValueError("sortv must contain only finite values")
    return numeric


def _proportional_records(
    sortorder: np.ndarray,
    weights: np.ndarray,
    k: int,
) -> list[_BoundaryRecord]:
    """Allocate sorted weighted intervals to equal groups, retaining boundary zeros."""
    sorted_weights = weights[sortorder]
    cumulative = np.cumsum(sorted_weights)
    group_size = float(np.sum(weights)) / k
    starts = np.zeros(k, dtype=int)
    ends = np.zeros(k, dtype=int)
    first_weights = np.zeros(k, dtype=float)
    last_weights = np.zeros(k, dtype=float)
    first_weights[0] = min(float(sorted_weights[0]), group_size)
    ends[-1] = len(sortorder) - 1
    last_weights[-1] = float(sorted_weights[-1])
    carried_weight = float(first_weights[0])

    for group in range(k - 1):
        boundary = (group + 1) * group_size
        boundary_position = min(
            int(np.searchsorted(cumulative, boundary, side="right")),
            len(sortorder) - 1,
        )
        ends[group] = boundary_position
        starts[group + 1] = boundary_position
        first_weights[group + 1] = min(
            float(cumulative[boundary_position]) - boundary,
            group_size,
        )
        if starts[group + 1] == starts[group]:
            carried_weight -= first_weights[group + 1]
        else:
            carried_weight = float(sorted_weights[boundary_position])
        if ends[group] > starts[group] + 1:
            assigned_before_last = float(np.sum(sorted_weights[starts[group] + 1 : ends[group]]))
        else:
            assigned_before_last = 0.0
        if ends[group] > starts[group]:
            assigned_before_last += first_weights[group]
            last_weights[group] = min(
                carried_weight - first_weights[group + 1],
                group_size - assigned_before_last,
            )

    records: list[_BoundaryRecord] = []
    for group in range(k):
        group_weights = sorted_weights[starts[group] : ends[group] + 1].copy()
        group_weights[-1] = last_weights[group]
        group_weights[0] = first_weights[group]
        for offset, assigned in enumerate(group_weights):
            position = starts[group] + offset
            records.append(
                _BoundaryRecord(
                    case=int(sortorder[position]),
                    source_weight=float(sorted_weights[position]),
                    group=group,
                    assigned_weight=float(assigned),
                )
            )
    return records


def _group_summaries(
    matrix: np.ndarray,
    records: list[_BoundaryRecord],
    k: int,
    *,
    squared: bool,
    exponent: float,
) -> tuple[np.ndarray, list[np.ndarray], list[np.ndarray], list[np.ndarray], float]:
    medoids = np.empty(k, dtype=int)
    index_list: list[np.ndarray] = []
    weights_list: list[np.ndarray] = []
    dist_list: list[np.ndarray] = []
    residual = 0.0

    for group in range(k):
        selected = [record for record in records if record.group == group]
        indexes = np.asarray([record.case for record in selected], dtype=int)
        group_weights = np.asarray([record.assigned_weight for record in selected], dtype=float)
        if indexes.size == 0 or float(np.sum(group_weights)) <= 0:
            raise ValueError("k creates an empty relative-frequency group")
        submatrix = matrix[np.ix_(indexes, indexes)]
        local_medoid = (
            0
            if indexes.size == 1
            else _medoid_index(submatrix, weights=group_weights, squared=squared)
        )
        distances = np.asarray(submatrix[:, local_medoid], dtype=float)
        medoids[group] = int(indexes[local_medoid])
        index_list.append(indexes)
        weights_list.append(group_weights)
        dist_list.append(distances)
        residual += float(np.sum(group_weights * distances**exponent))

    return medoids, index_list, weights_list, dist_list, residual


def _fit_statistics(
    matrix: np.ndarray,
    medoid: int,
    residual: float,
    *,
    exponent: float,
    k: int,
    weight_sum: float,
) -> tuple[float, float, float]:
    denominator = float(np.sum(matrix[:, medoid] ** exponent))
    with np.errstate(divide="ignore", invalid="ignore"):
        r_squared = float(1.0 - np.divide(residual, denominator))
        explained = np.divide(r_squared, k - 1)
        unexplained = np.divide(1.0 - r_squared, weight_sum - k)
        f_statistic = float(np.divide(explained, unexplained))
    if k > 1 and weight_sum > k:
        p_value = float(1.0 - stats.f.cdf(f_statistic, k - 1, weight_sum - k))
    else:
        # R (dissrf.R:262) warns here instead of silently reporting NA, with msg.warn's text:
        # its parts pasted with spaces after " [!!]", numbers as format() writes them
        # (TX-REP-20261023-008).
        warnings.warn(
            f" [!!] Negative  or zero df! k = {r_number_text(k, digits=7)}  wsum =  "
            f"{r_number_text(weight_sum, digits=7)}",
            stacklevel=2,
        )
        p_value = float("nan")
    return r_squared, f_statistic, p_value


def dissrf(
    diss: np.ndarray,
    *,
    k: Optional[int] = None,
    sortv="mds",
    weights: Optional[np.ndarray] = None,
    grp_meth: Literal["prop", "first", "random"] = "prop",
    squared: bool = False,
    pow: Optional[int] = None,
    random_state: int | np.random.Generator | None = None,
) -> Dict[str, Any]:
    """Partition sorted cases into TraMineR relative-frequency groups.

    ``medoids`` holds the zero-based row of each group's medoid and ``index_list`` the
    zero-based rows of each group (R's are one-based).
    The default ``sortv="mds"`` axis is oriented so that its coordinate of
    largest absolute value is positive (coordinates within a relative 1e-10 of
    it are tied, and the first decides): R takes the sign from LAPACK, so its
    groups follow the numerical library (TX-REP-20261017-006).
    """
    matrix = _validate_dissimilarities(diss)
    ncase = matrix.shape[0]
    resolved_weights, weighted = _normalize_weights(ncase, weights)
    if weighted and grp_meth != "prop":
        raise ValueError(
            f"Selected grp.meth '{grp_meth}' does not apply to weighted data! "
            "Use grp.meth='prop'."
        )
    if grp_meth not in {"prop", "first", "random"}:
        raise ValueError("grp_meth must be one of 'first', 'random', 'prop'")

    weight_sum = float(np.sum(resolved_weights))
    if k is None:
        k = min(int(np.floor(weight_sum / 10.0)), 100)
    if isinstance(k, bool) or not isinstance(k, (int, np.integer)):
        raise TypeError("k must be an integer or None")
    k = int(k)
    if not 1 <= k <= ncase:
        raise ValueError("k must be between 1 and the number of cases")
    if pow is None:
        exponent = 2.0 if squared else 1.0
    else:
        exponent = float(pow)
        if not np.isfinite(exponent) or exponent <= 0:
            raise ValueError("pow must be finite and positive")

    overall_medoid = _medoid_index(
        matrix,
        weights=resolved_weights,
        squared=bool(squared),
    )
    overall_distances = np.asarray(matrix[:, overall_medoid], dtype=float)
    if k == 1:
        indexes = np.arange(ncase, dtype=int)
        result: Dict[str, Any] = {
            "medoids": np.asarray([overall_medoid], dtype=int),
            "med_names": [overall_medoid],
            "dist_list": [overall_distances],
            "index_list": [indexes],
            "weights_list": [resolved_weights.copy()],
            "heights": np.ones(1, dtype=float),
            "at": np.asarray([0.5], dtype=float),
            "R2": 0.0,
            "Fstat": 0.0,
            "pvalue": float("nan"),
            "sizes": {"ncase": ncase, "wsum": weight_sum, "k": 1, "gsize": ncase},
            "grp_meth": "prop",
        }
        return result

    axis = _resolve_axis(
        matrix,
        sortv,
        weights=resolved_weights,
        weighted=weighted,
        squared=bool(squared),
    )
    sortorder = np.argsort(axis, kind="stable")

    if grp_meth == "prop":
        records = _proportional_records(sortorder, resolved_weights, k)
        medoids, index_list, weights_list, dist_list, residual = _group_summaries(
            matrix,
            records,
            k,
            squared=bool(squared),
            exponent=exponent,
        )
        r_squared, f_statistic, p_value = _fit_statistics(
            matrix,
            overall_medoid,
            residual,
            exponent=exponent,
            k=k,
            weight_sum=weight_sum,
        )
        return {
            "medoids": medoids,
            "med_names": medoids.tolist(),
            "dist_list": dist_list,
            "index_list": index_list,
            "weights_list": weights_list,
            "heights": np.full(k, 1.0 / k, dtype=float),
            "at": np.arange(k, dtype=float) + 0.5,
            "R2": r_squared,
            "Fstat": f_statistic,
            "pvalue": p_value,
            "sizes": {
                "ncase": ncase,
                "wsum": weight_sum,
                "k": k,
                "gsize": weight_sum / k,
            },
            "grp_meth": grp_meth,
        }

    base_size, remainder = divmod(ncase, k)
    group_sizes = np.full(k, base_size, dtype=int)
    if remainder:
        if grp_meth == "first":
            enlarged = np.arange(remainder, dtype=int)
        else:
            rng = (
                random_state
                if isinstance(random_state, np.random.Generator)
                else np.random.default_rng(random_state)
            )
            enlarged = rng.choice(k, size=remainder, replace=False)
        group_sizes[enlarged] += 1

    memberships = np.empty(ncase, dtype=int)
    cursor = 0
    records = []
    for group, group_size in enumerate(group_sizes):
        group_cases = sortorder[cursor : cursor + group_size]
        memberships[group_cases] = group
        cursor += int(group_size)
    for group in range(k):
        group_cases = np.flatnonzero(memberships == group)
        records.extend(_BoundaryRecord(int(case), 1.0, group, 1.0) for case in group_cases)

    medoids, index_list, weights_list, dist_list, residual = _group_summaries(
        matrix,
        records,
        k,
        squared=bool(squared),
        exponent=exponent,
    )
    case_medoids = medoids[memberships]
    case_distances = matrix[np.arange(ncase), case_medoids]
    heights = group_sizes.astype(float) / ncase
    at = (np.cumsum(heights) - heights / 2.0) * k
    r_squared, f_statistic, p_value = _fit_statistics(
        matrix,
        overall_medoid,
        residual,
        exponent=exponent,
        k=k,
        weight_sum=weight_sum,
    )
    return {
        "medoids": medoids,
        "med_names": medoids.tolist(),
        "dist_list": dist_list,
        "index_list": index_list,
        "weights_list": weights_list,
        "kmedoid_index": case_medoids.astype(int, copy=False),
        "kmedoid_dist": np.asarray(case_distances, dtype=float),
        "mdsk": memberships,
        "heights": heights,
        "at": at,
        "R2": r_squared,
        "Fstat": f_statistic,
        "pvalue": p_value,
        "sizes": {"ncase": ncase, "k": k, "ng": base_size, "r": remainder},
        "grp_meth": grp_meth,
    }


def _resolve_sortv(x: STSList, sortv):
    if not isinstance(sortv, str) or sortv not in {"from.start", "from.end"}:
        return sortv
    # seqrf.R:11-16 orders over the factor columns, c(alphabet, nr, void),
    # whatever codes the markers are stored under (TX-AUD-20261011-014).
    matrix = x._canonical_code_ranks(void_as_missing=False)[np.asarray(x.data)]
    reverse_time = sortv == "from.end"
    order = sorted(
        range(matrix.shape[0]),
        key=lambda row: tuple(matrix[row, ::-1] if reverse_time else matrix[row]),
    )
    ranks = np.empty(matrix.shape[0], dtype=int)
    ranks[np.asarray(order, dtype=int)] = np.arange(matrix.shape[0])
    return ranks


def _medoid_sequence_subset(
    source: STSList,
    medoids: np.ndarray,
    heights: np.ndarray,
) -> STSList:
    """Copy medoid rows and all compatible rendering metadata."""
    rows = np.asarray(medoids, dtype=int).reshape(-1)
    return STSList(
        data=np.asarray(source.data[rows], dtype=int).copy(),
        alphabet=list(source.alphabet),
        labels=None if source.labels is None else list(source.labels),
        cpal=None if source.cpal is None else list(source.cpal),
        xtlab=None if source.xtlab is None else list(source.xtlab),
        xtstep=source.xtstep,
        tick_last=source.tick_last,
        nr=source.nr,
        void=source.void,
        nr_code=source.nr_code,
        void_code=source.void_code,
        weights=np.asarray(heights, dtype=float).copy(),
        start=source.start,
        ids=None if source.ids is None else [source.ids[int(index)] for index in rows],
        missing_color=source.missing_color,
        issues=_remap_issues(source.issues, rows),
    )


def seqrf(
    x,
    diss: np.ndarray,
    *,
    k: Optional[int] = None,
    sortv="mds",
    weights: Optional[np.ndarray] = None,
    weighted: bool = True,
    grp_meth: Literal["prop", "first", "random"] = "prop",
    squared: bool = False,
    pow: Optional[int] = None,
    random_state: int | np.random.Generator | None = None,
) -> Dict[str, Any]:
    """Return RF medoids and representative sequences for an ``STSList``.

    ``medoids`` holds the zero-based row of each group's medoid (R's are one-based).
    The legacy ``seqtoplot`` array and RF mappings remain unchanged.
    ``medoid_sequences`` is the metadata-complete public ``STSList`` used by
    external renderers; its row weights are the RF group heights.
    The default ``sortv="mds"`` axis is oriented so that its coordinate of
    largest absolute value is positive (coordinates within a relative 1e-10 of
    it are tied, and the first decides): R takes the sign from LAPACK, so its
    groups follow the numerical library (TX-REP-20261017-006).
    """
    if not isinstance(x, STSList):
        raise TypeError("x must be a state sequence object")
    matrix = _validate_dissimilarities(diss)
    if matrix.shape[0] != x.nseq:
        raise ValueError("diss size must equal the number of sequences")
    resolved_sort = _resolve_sortv(x, sortv)
    resolved_weights = weights
    if not weighted:
        resolved_weights = None
    elif resolved_weights is None and x.weights is not None:
        resolved_weights = np.asarray(x.weights, dtype=float)

    rf = dissrf(
        matrix,
        k=k,
        sortv=resolved_sort,
        weights=resolved_weights,
        grp_meth=grp_meth,
        squared=squared,
        pow=pow,
        random_state=random_state,
    )
    medoids = np.asarray(rf["medoids"], dtype=int)
    heights = np.asarray(rf["heights"], dtype=float)
    return {
        "seqtoplot": np.asarray(x.data[medoids], dtype=int).copy(),
        "rf": rf,
        "medoids": medoids.copy(),
        "heights": heights.copy(),
        "medoid_sequences": _medoid_sequence_subset(x, medoids, heights),
    }
