"""Case-level marginality and gain indicators for dissimilarity data.

The implementation follows the weighted projection equations used by
TraMineRextras' ``dissindic``.  It is intentionally expressed directly in
terms of a Gower matrix and a categorical projection, rather than through the
multifactor-analysis implementation.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd


def _as_inputs(
    diss: np.ndarray,
    group: np.ndarray,
    weights: Optional[np.ndarray],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Validate inputs and return detached one-dimensional arrays."""
    matrix = np.asarray(diss, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("diss must be a square dissimilarity matrix")

    groups = np.asarray(group)
    if groups.ndim != 1 or groups.shape[0] != matrix.shape[0]:
        raise ValueError("group length must match nrow(diss)")

    if weights is None:
        case_weights = np.ones(matrix.shape[0], dtype=float)
    else:
        case_weights = np.asarray(weights, dtype=float)
        if case_weights.ndim != 1 or case_weights.shape[0] != matrix.shape[0]:
            raise ValueError("weights must have one value per observation")
        if not np.all(np.isfinite(case_weights)):
            raise ValueError("weights must be finite")
        # R computes with a zero weight (that case's residual is 0/0, NaN; TX-AUD-20261011-066).
        if np.any(case_weights < 0):
            raise ValueError("weights must be non-negative")

    return matrix, groups.copy(), case_weights.copy()


def _gower_from_distances(
    distances: np.ndarray,
    weights: np.ndarray,
    *,
    squared: bool,
) -> np.ndarray:
    """Weighted Gower transform from the equations in TraMineR."""
    transformed = distances**2 if squared else distances.copy()
    total = float(weights.sum())
    if not np.isfinite(total) or total <= 0:
        raise ValueError("weights must have a positive sum")
    shares = weights / total
    n_cases = distances.shape[0]
    left = np.eye(n_cases) - np.outer(np.ones(n_cases), shares)
    right = np.eye(n_cases) - np.outer(shares, np.ones(n_cases))
    return left @ (-0.5 * transformed) @ right


def _categorical_projection(groups: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Return the weighted projection for a complete categorical predictor."""
    # ``factorize(sort=True)`` gives the same column space as R's intercept plus
    # treatment contrasts while remaining well-defined for mixed scalar labels.
    codes, levels = pd.factorize(groups, sort=True)
    # A level whose weights are all zero has a zero weighted column: leaving it
    # out keeps the projection on the design's column space, so its cases take
    # no part (batch-12 review I-1). A QR of a zero column gives an arbitrary Q
    # column (R's qr.Q does too when that level is the reference one).
    weighted_levels = np.flatnonzero(np.bincount(codes, weights=weights > 0, minlength=len(levels)))
    if weighted_levels.size == 0:
        return np.zeros((groups.size, groups.size), dtype=float)
    design = np.equal.outer(codes, weighted_levels).astype(float)
    root_weights = np.sqrt(weights)
    q_matrix, _ = np.linalg.qr(root_weights[:, None] * design, mode="reduced")
    return np.outer(root_weights, root_weights) * (q_matrix @ q_matrix.T)


def dissindic(
    diss: np.ndarray,
    group: np.ndarray,
    *,
    gower: bool = False,
    squared: bool = False,
    weights: Optional[np.ndarray] = None,
) -> pd.DataFrame:
    """Compute TraMineRextras-compatible marginality and gain indicators."""
    matrix, groups, case_weights = _as_inputs(diss, group, weights)
    complete = ~pd.isna(groups)

    marginality = np.full(groups.shape[0], np.nan, dtype=float)
    gain = np.full(groups.shape[0], np.nan, dtype=float)
    if np.any(complete):
        observed_matrix = matrix[np.ix_(complete, complete)]
        observed_weights = case_weights[complete]
        gram = (
            observed_matrix.copy()
            if gower
            else _gower_from_distances(
                observed_matrix,
                observed_weights,
                squared=bool(squared),
            )
        )
        projection = _categorical_projection(groups[complete], observed_weights)
        weight_matrix = np.diag(observed_weights)
        with np.errstate(divide="ignore", invalid="ignore"):
            residual = np.diag((weight_matrix - projection) @ gram) / observed_weights
        null_residual = np.diag(gram)
        marginality[complete] = residual
        gain[complete] = null_residual - residual

    return pd.DataFrame(
        {
            "group": groups,
            "marginality": marginality,
            "gain": gain,
        }
    )
