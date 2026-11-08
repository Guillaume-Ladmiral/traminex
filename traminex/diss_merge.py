"""Greedy merging of groups using weighted average silhouette width."""

from __future__ import annotations

import warnings
from itertools import combinations
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd

from ._categorical_order import r_factor


def _renumber(labels: np.ndarray) -> np.ndarray:
    """R's ``as.integer(factor(group))``: level numbers from one, in R's level order.

    Text groups are numbered in R's collation (dissmergegroups.R:105,
    TX-AUD-20261010-039), numbers equal to 15 significant digits share a number
    (TX-AUD-20261019-046); callers refuse missing groups first.
    """
    return r_factor(np.asarray(labels))[1] + 1


def _wc_asw(D: np.ndarray, labels: np.ndarray, w: np.ndarray) -> float:
    """WeightedCluster's corrected weighted ASW statistic.

    The within-cluster denominator is the weighted cluster size minus one,
    exactly as in WeightedCluster's ``ClusterQualASWi`` implementation.
    """
    matrix = np.asarray(D, dtype=float)
    membership = _renumber(np.asarray(labels)) - 1
    weights = np.asarray(w, dtype=float)
    n_cases = matrix.shape[0]
    if matrix.ndim != 2 or matrix.shape != (n_cases, n_cases):
        raise ValueError("diss must be a square dissimilarity matrix")
    if membership.shape != (n_cases,) or weights.shape != (n_cases,):
        raise ValueError("diss, labels, and weights must have conformable lengths")

    n_groups = int(membership.max()) + 1 if n_cases else 0
    if n_groups < 2:
        raise ValueError("group must contain at least two different values")
    sizes = np.bincount(membership, weights=weights, minlength=n_groups)
    if np.any(sizes < 1.0):
        return float("nan")

    total = float(sizes.sum())
    weighted_score = 0.0
    for case in range(n_cases):
        if weights[case] <= 0:
            continue
        own = membership[case]
        sums = np.bincount(
            membership,
            weights=weights * matrix[case],
            minlength=n_groups,
        )
        within_sum = float(sums[own])
        within = 0.0 if sizes[own] <= 1.0 else within_sum / (sizes[own] - 1.0)
        between = min(
            float(sums[group] / sizes[group]) for group in range(n_groups) if group != own
        )
        scale = max(within, between)
        score = float("nan") if scale == 0 else (between - within) / scale
        weighted_score += weights[case] * score
    return float(weighted_score / total)


def _validated_inputs(
    diss: np.ndarray,
    group: np.ndarray,
    weights: Optional[np.ndarray],
    *,
    squared: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    matrix = np.asarray(diss, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("diss must be a square dissimilarity matrix")
    if squared:
        matrix = matrix**2
    else:
        matrix = matrix.copy()

    raw_groups = np.asarray(group)
    if raw_groups.ndim != 1 or raw_groups.shape[0] != matrix.shape[0]:
        raise ValueError("length of group non conformable with diss")
    if np.any(pd.isna(raw_groups)):
        raise ValueError("group must not contain missing values")

    if weights is None:
        case_weights = np.ones(matrix.shape[0], dtype=float)
    else:
        case_weights = np.asarray(weights, dtype=float)
        if case_weights.ndim != 1 or case_weights.shape[0] != matrix.shape[0]:
            raise ValueError("weights must have one value per observation")
        if np.any(~np.isfinite(case_weights)) or np.any(case_weights < 0):
            raise ValueError("weights must be finite and non-negative")
    if float(case_weights.sum()) <= 0:
        raise ValueError("weights must have a positive sum")
    return matrix, _renumber(raw_groups), case_weights.copy()


def dissmergegroups(
    diss: np.ndarray,
    group: np.ndarray,
    *,
    weights: Optional[np.ndarray] = None,
    measure: str = "ASW",
    crit: float = 0.2,
    ref: str = "max",
    min_group: int = 4,
    small: float = 0.05,
    squared: bool = False,
    silent: bool = True,
) -> Dict[str, Any]:
    """Greedily merge the pair of groups with the least ASW deterioration.

    ``silent`` is the one default here that differs from R's: ``dissmergegroups.R:77``
    defaults it to ``FALSE`` and prints a line per merge, while a library is
    better quiet, so this defaults to ``True``. It controls only the printed
    progress -- the returned ``group`` and ``quality`` are identical either way
    (cycle-7 W-07; pinned by
    ``tests/test_golden_default_coverage.py::test_dissmergegroups_is_quiet_by_default_but_otherwise_identical``).
    Every other default matches R's.

    The ASW is undefined when a group's weighted size is below 1. As in R,
    merging is then refused (``ValueError``; R stops at its first comparison),
    and when no merge is attempted (``min_group`` already reached) the grouping
    is returned with quality NaN and a ``RuntimeWarning``.

    When no merge leaves a positive quality (no pair of groups loses less than
    the reference quality), merging stops and the grouping is returned as it
    is. R's loop then ends only when the reference quality is above ``crit``
    times itself (positive, with ``crit`` below 1); otherwise, as for a
    negative reference with no merge of positive quality, it takes its merge
    branch with nothing to merge and does not end (dissmergegroups.R:127-168).
    This is a recorded correction (TX-AUD-20261011-016). A negative reference
    quality alone does not cause it: a merge of positive quality is made as in
    R.
    """
    if measure != "ASW":
        raise ValueError("only measure='ASW' is supported")
    if ref not in {"initial", "max", "previous"}:
        raise ValueError("Unknown ref value!")
    if not 0 <= float(crit) <= 1:
        raise ValueError("crit must be between 0 and 1")
    if int(min_group) < 1:
        raise ValueError("min_group must be at least one")
    if float(small) < 0:
        raise ValueError("small must be non-negative")

    matrix, membership, case_weights = _validated_inputs(
        diss,
        group,
        weights,
        squared=bool(squared),
    )
    n_groups = int(membership.max())
    if n_groups < 2:
        raise ValueError("group must contain at least two different values")

    population = float(case_weights.sum())
    minimum_size = float(small) * population if small < 1 else float(small)
    initial_sizes = np.bincount(membership - 1, weights=case_weights, minlength=n_groups)
    if np.any(initial_sizes < 1.0):
        # WeightedCluster's wcClusterQuality cannot compute the ASW then: it
        # warns and returns NA (dissmergegroups.R:124). R's merge loop compares
        # that NA and stops (lines 140, 155); with no merge to attempt, R
        # returns the grouping with quality NA. Merging only enlarges groups,
        # so only the initial grouping can be affected.
        message = "ASW can not be computed because at least one group has a weighted size below one"
        if n_groups > int(min_group):
            raise ValueError(f"{message}, so groups cannot be merged on it")
        warnings.warn(message, RuntimeWarning, stacklevel=2)
    initial_quality = _wc_asw(matrix, membership, case_weights)
    final_quality = initial_quality
    reference_quality = initial_quality

    while n_groups > int(min_group):
        sizes = np.bincount(
            membership - 1,
            weights=case_weights,
            minlength=n_groups,
        )
        if float(sizes.min()) > minimum_size:
            candidates = list(combinations(range(1, n_groups + 1), 2))
        else:
            smallest = int(np.argmin(sizes)) + 1
            candidates = [
                (min(smallest, other), max(smallest, other))
                for other in range(1, n_groups + 1)
                if other != smallest
            ]

        best_loss = reference_quality
        best_pair: tuple[int, int] | None = None
        for left, right in candidates:
            proposal = membership.copy()
            proposal[membership == right] = left
            quality = _wc_asw(matrix, proposal, case_weights)
            loss = reference_quality - quality
            if loss < best_loss:
                best_loss = loss
                best_pair = (left, right)

        if best_pair is None or best_loss > float(crit) * reference_quality:
            break

        left, right = best_pair
        if not silent:
            print(f"Merging groups {left} and {right}")
        membership[membership == right] = left
        membership = _renumber(membership)
        n_groups = int(membership.max())
        final_quality = reference_quality - best_loss
        if ref == "max":
            reference_quality = max(reference_quality, final_quality)
        elif ref == "previous":
            reference_quality = final_quality
        else:
            reference_quality = initial_quality

    return {"group": membership, "quality": float(final_quality)}
