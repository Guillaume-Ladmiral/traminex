"""Weighted multifactor analysis of a dissimilarity matrix."""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import List, Mapping, Optional, Sequence

import numpy as np
import pandas as pd

from ._rlinpack import dqrdc2_pivot
from .diss import _permutation_pvalues


@dataclass
class DissMultifactor:
    variables: List[str]
    pseudo_F: np.ndarray
    pseudo_R2: np.ndarray
    p_value: np.ndarray
    SS: np.ndarray
    mfac: pd.DataFrame


def _gower_matrix_weighted(
    diss: np.ndarray,
    squared: bool,
    weights: np.ndarray,
) -> np.ndarray:
    """Weighted Gower transform from TraMineR's matrix equation."""
    matrix = np.asarray(diss, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("diss must be a square dissimilarity matrix")
    case_weights = np.asarray(weights, dtype=float)
    if case_weights.shape != (matrix.shape[0],):
        raise ValueError("weights must have length equal to nrow(diss)")
    # Any total, as R's s <- weights/sum(weights) (dissmfacw.R:20): a zero total gives NaN
    # and a negative one finite values, without a warning (TX-AUD-20261011-005). dissmfacw
    # refuses such weights before it gets here.
    total = float(case_weights.sum())
    distances = matrix**2 if squared else matrix.copy()
    ones = np.ones(matrix.shape[0], dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        shares = case_weights / total
        left = np.eye(matrix.shape[0]) - np.outer(ones, shares)
        right = np.eye(matrix.shape[0]) - np.outer(shares, ones)
        return left @ (-0.5 * distances) @ right


# R's qr() tolerance: dqrdc2 sets aside a column whose remaining norm falls below this
# share of its own original norm (TX-REP-20261019-002).
_RANK_TOLERANCE = 1e-7


def _weighted_qr(design: np.ndarray, weights: np.ndarray) -> tuple[np.ndarray, np.ndarray, int]:
    """An orthonormal basis of the columns R's ``qr()`` keeps in the weighted design, and its rank.

    R decomposes ``W_sqrt * pred`` (dissmfacw.R:67); the rank and the kept columns are dqrdc2's.
    """
    matrix = np.asarray(design, dtype=float)
    case_weights = np.asarray(weights, dtype=float)
    if matrix.ndim != 2 or case_weights.shape != (matrix.shape[0],):
        raise ValueError("design and weights must have conformable rows")
    root_weights = np.sqrt(case_weights)
    weighted = root_weights[:, None] * matrix
    rank, pivot = dqrdc2_pivot(weighted, _RANK_TOLERANCE)
    q_matrix, _ = np.linalg.qr(weighted[:, pivot[:rank]])
    return q_matrix, root_weights, rank


def _design_rank(design: np.ndarray, weights: np.ndarray) -> int:
    return _weighted_qr(design, weights)[2]


def _weighted_hat_matrix_qr(design: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Weighted projection onto the column space of the design.

    R builds it from all ``ncol(X)`` columns of ``qr.Q(qr(X))`` (dissmfacw.R:67-69),
    which for a rank-deficient design adds one arbitrary direction; the
    projection here is onto the columns R's dqrdc2 keeps, so it is the same for
    any design spanning the same space (TX-AUD-20261010-040, TX-REP-20261019-002).
    """
    q_matrix, root_weights, rank = _weighted_qr(design, weights)
    basis = q_matrix[:, :rank]
    return np.outer(root_weights, root_weights) * (basis @ basis.T)


def _as_predictor_series(values: Sequence) -> pd.Series:
    if isinstance(values, pd.Series):
        return values.reset_index(drop=True).copy()
    array = np.asarray(values)
    if array.ndim != 1:
        raise ValueError("predictors must be one-dimensional")
    return pd.Series(array.copy())


def _is_numeric(series: pd.Series) -> bool:
    return bool(
        pd.api.types.is_numeric_dtype(series.dtype) and not pd.api.types.is_bool_dtype(series.dtype)
    )


def _build_design(
    factors: Mapping[str, Sequence],
) -> tuple[np.ndarray, np.ndarray, List[str]]:
    """Construct an R-model-matrix-style intercept and treatment contrasts."""
    names = list(factors)
    if not names:
        raise ValueError("factors must contain at least one predictor")
    series = [_as_predictor_series(factors[name]) for name in names]
    n_cases = len(series[0])
    if any(len(values) != n_cases for values in series):
        raise ValueError("all predictors must have the same length")

    columns = [np.ones(n_cases, dtype=float)]
    assignments = [0]
    for term, (name, values) in enumerate(zip(names, series), start=1):
        if _is_numeric(values):
            numeric = values.to_numpy(dtype=float)
            if not np.all(np.isfinite(numeric)):
                raise ValueError(f"numeric predictor {name!r} must be finite")
            columns.append(numeric)
            assignments.append(term)
            continue

        categories = sorted(pd.unique(values).tolist(), key=lambda value: str(value))
        if len(categories) < 2:
            raise ValueError(f"categorical predictor {name!r} must contain at least two levels")
        raw = values.to_numpy()
        for category in categories[1:]:
            columns.append((raw == category).astype(float))
            assignments.append(term)
    return (
        np.column_stack(columns),
        np.asarray(assignments, dtype=int),
        names,
    )


@np.errstate(divide="ignore", invalid="ignore")
def _deterministic_stats(
    X: np.ndarray,
    perm_X: np.ndarray,
    G: np.ndarray,
    var_list: np.ndarray,
    n_terms: int,
    p_list: np.ndarray,
    ss_backward: np.ndarray,
    sc_tot: float,
    total_weight: float,
    weights: np.ndarray,
    full_rank: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Evaluate term and complete-model statistics for one row permutation.

    The divisions follow NumPy (IEEE) semantics, as R's do: an all-zero
    dissimilarity matrix or a perfect fit yields NaN/Inf statistics, as in
    R's dissmfacw, rather than raising ``ZeroDivisionError``. Degrees of
    freedom are ranks (``full_rank`` for the model, ``p_list`` per term); a term
    with none has no test (NaN).
    """
    sc_tot = np.float64(sc_tot)
    total_weight = np.float64(total_weight)
    complete_index = n_terms
    pseudo_f = np.zeros(n_terms + 1, dtype=float)
    pseudo_r2 = np.zeros(n_terms + 1, dtype=float)
    explained_ss = np.zeros(n_terms + 1, dtype=float)
    column_count = full_rank

    complete_hat = _weighted_hat_matrix_qr(perm_X, weights)
    complete_explained = float(np.sum(complete_hat * G))
    complete_residual = sc_tot - complete_explained
    pseudo_r2[complete_index] = complete_explained / sc_tot
    # A model whose rank is the intercept's alone has no test.
    pseudo_f[complete_index] = (
        np.nan
        if column_count <= 1
        else (complete_explained / (column_count - 1))
        / (complete_residual / (total_weight - column_count))
    )
    explained_ss[complete_index] = complete_explained

    if n_terms == 1:
        pseudo_f[0] = pseudo_f[complete_index]
        pseudo_r2[0] = pseudo_r2[complete_index]
        explained_ss[0] = complete_explained
        return pseudo_f, pseudo_r2, explained_ss

    for term in range(1, n_terms + 1):
        mixed_design = X.copy()
        term_columns = var_list == term
        mixed_design[:, term_columns] = perm_X[:, term_columns]
        term_hat = _weighted_hat_matrix_qr(mixed_design, weights)
        model_explained = float(np.sum(term_hat * G))
        effect_ss = model_explained - float(ss_backward[term - 1])
        residual = sc_tot - model_explained
        if p_list[term - 1] == 0:
            pseudo_f[term - 1] = pseudo_r2[term - 1] = explained_ss[term - 1] = np.nan
            continue
        pseudo_f[term - 1] = (effect_ss / p_list[term - 1]) / (
            residual / (total_weight - column_count)
        )
        pseudo_r2[term - 1] = effect_ss / sc_tot
        explained_ss[term - 1] = effect_ss
    return pseudo_f, pseudo_r2, explained_ss


def _complete_predictors(
    factors: Mapping[str, Sequence],
    n_cases: int,
) -> tuple[dict[str, pd.Series], np.ndarray]:
    prepared: dict[str, pd.Series] = {}
    complete = np.ones(n_cases, dtype=bool)
    for name, values in factors.items():
        series = _as_predictor_series(values)
        if len(series) != n_cases:
            raise ValueError(f"predictor {name!r} has length {len(series)}, expected {n_cases}")
        prepared[name] = series
        complete &= ~series.isna().to_numpy()
    if not prepared:
        raise ValueError("factors must contain at least one predictor")
    return prepared, complete


def dissmfacw(
    diss: np.ndarray,
    factors: Mapping[str, Sequence],
    weights: Optional[Sequence[float]] = None,
    *,
    nperm: int = 1000,
    seed: Optional[int] = None,
    gower: bool = False,
    squared: bool = False,
) -> DissMultifactor:
    """Perform TraMineR-compatible Type-II multifactor discrepancy analysis.

    Permutation p-values follow R's ``TraMineR.permutation``: the share of
    the ``nperm`` rows, the observed one included, whose pseudo F is at least
    the observed one; NaN (R's NA) when ``nperm`` is 1 or when the statistic
    is NaN in the observed data or in any permutation.

    A rank-deficient design (collinear predictors, a constant) is projected onto
    its column space and its degrees of freedom are ranks, where R counts
    columns and projects onto one arbitrary extra direction: the model and the
    estimable terms have the values of the same design without its redundant
    columns, and a term adding no rank given the others has NaN statistics, with
    a ``UserWarning`` (TX-AUD-20261010-040). Terms are the named predictors of
    ``factors``; there is no formula interface or interaction term.
    """
    matrix = np.asarray(diss, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("diss must be a square dissimilarity matrix")
    n_cases = matrix.shape[0]
    prepared, complete = _complete_predictors(factors, n_cases)

    if weights is None:
        all_weights = np.ones(n_cases, dtype=float)
    else:
        all_weights = np.asarray(weights, dtype=float)
        if all_weights.shape != (n_cases,):
            raise ValueError(f"weights must have length {n_cases}")
    case_weights = all_weights[complete].copy()
    if not np.all(np.isfinite(case_weights)):
        raise ValueError("weights must be finite")
    if np.any(case_weights < 0):
        raise ValueError("weights must be non-negative")
    if float(case_weights.sum()) <= 0:
        raise ValueError("weights must have a positive sum")

    complete_factors = {
        name: values.loc[complete].reset_index(drop=True) for name, values in prepared.items()
    }
    design, assignments, names = _build_design(complete_factors)
    observed_matrix = matrix[np.ix_(complete, complete)]
    gram = (
        observed_matrix.copy()
        if gower
        else _gower_matrix_weighted(observed_matrix, bool(squared), case_weights)
    )
    total_weight = float(case_weights.sum())
    sc_total = float(np.sum(case_weights * np.diag(gram)))
    n_terms = len(names)
    # Degrees of freedom are ranks: R counts columns (dissmfacw.R:48-52, 79-90),
    # which a rank-deficient design overstates (TX-AUD-20261010-040).
    full_rank = _design_rank(design, case_weights)
    p_list = np.zeros(n_terms, dtype=float)
    backward_ss = np.zeros(n_terms, dtype=float)
    for term in range(1, n_terms + 1):
        reduced = design[:, assignments != term]
        p_list[term - 1] = full_rank - _design_rank(reduced, case_weights)
        reduced_hat = _weighted_hat_matrix_qr(reduced, case_weights)
        backward_ss[term - 1] = float(np.sum(reduced_hat * gram))
    if full_rank < design.shape[1]:
        untested = [name for name, df in zip(names, p_list) if df == 0]
        warnings.warn(
            f"dissmfacw: the design has rank {full_rank} for {design.shape[1]} columns; "
            f"terms adding no rank given the others have no test: {untested}",
            UserWarning,
            stacklevel=2,
        )

    permutation_count = int(nperm)
    if permutation_count < 1 or permutation_count != nperm:
        raise ValueError("nperm must be a positive integer")
    rng = np.random.default_rng(seed)
    observed_f, observed_r2, observed_ss = _deterministic_stats(
        design,
        design,
        gram,
        assignments,
        n_terms,
        p_list,
        backward_ss,
        sc_total,
        total_weight,
        case_weights,
        full_rank,
    )
    statistics = [observed_f]
    for _ in range(1, permutation_count):
        permutation = rng.permutation(design.shape[0])
        permuted_f, _, _ = _deterministic_stats(
            design,
            design[permutation],
            gram,
            assignments,
            n_terms,
            p_list,
            backward_ss,
            sc_total,
            total_weight,
            case_weights,
            full_rank,
        )
        statistics.append(permuted_f)
    # R's TraMineR.permutation (dissmfacw.R:124) reports NA p-values when
    # R <= 1 and NA for a statistic whose column holds a NaN, observed or
    # permuted (TraMineR-permutation.R:21, 25-26).
    if permutation_count > 1:
        p_values = _permutation_pvalues(np.vstack(statistics), observed_f)
    else:
        p_values = np.full(n_terms + 1, np.nan)
    variables = [*names, "Total"]
    table = pd.DataFrame(
        {
            "Variable": variables,
            "PseudoF": observed_f,
            "PseudoR2": observed_r2,
            "p_value": p_values,
        }
    )
    return DissMultifactor(
        variables=variables,
        pseudo_F=observed_f,
        pseudo_R2=observed_r2,
        p_value=p_values,
        SS=observed_ss,
        mfac=table,
    )
