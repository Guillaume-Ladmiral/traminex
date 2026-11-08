"""Multi-domain sequence analysis utilities.

This module provides tools for optimizing channel weights in multi-domain
(multi-channel) sequence data using clustering quality metrics.

These weight-search utilities operate over the **additive** per-channel distance
combination -- the non-R-faithful :func:`traminex.seqMD_distances` extension
(``D = sum_c w_c * seqdist(channel_c)``), *not* the faithful TraMineR CAT engine
(:func:`traminex.seqMD` / :func:`traminex.seqdistmc`). They take the lightweight
:class:`traminex.SeqMD` channel container.
"""

from __future__ import annotations

from itertools import product
from numbers import Real
from operator import index
from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd

from .sequences import SeqMD, _validate_analytical_weights, _validate_seqmd_channels


def _require_seqmd(obj: object) -> None:
    """Refuse an object that is not a SeqMD before its channels are read (TX-AUD-20261019-036)."""
    if not isinstance(obj, SeqMD):
        raise TypeError("obj must be a SeqMD multichannel sequence object")


def _validate_real(value, *, name: str, positive: bool) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real number")
    numeric = float(value)
    if not np.isfinite(numeric) or (numeric <= 0 if positive else numeric < 0):
        domain = "finite and strictly positive" if positive else "finite and non-negative"
        raise ValueError(f"{name} must be {domain}")
    return numeric


def _validate_positive_integer(value, *, name: str) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be a positive integer")
    try:
        numeric = index(value)
    except TypeError:
        raise TypeError(f"{name} must be a positive integer") from None
    if numeric <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return numeric


def _validate_metric(metric: str) -> None:
    if metric != "silhouette":
        raise ValueError(f"Unsupported metric: {metric}. Currently only 'silhouette' is supported.")


def _precompute_channel_distances(
    obj: SeqMD, method: str, per_channel_kwargs: Optional[List[dict]]
) -> List[np.ndarray]:
    """Compute one distance matrix per channel once, then reuse in optimizers."""
    from .distances import seqdist

    channels, _ = _validate_seqmd_channels(obj.channels)
    n_channels = len(channels)
    if per_channel_kwargs is None:
        kwargs_per_channel: List[dict] = [{} for _ in range(n_channels)]
    else:
        if len(per_channel_kwargs) != n_channels:
            raise ValueError(f"per_channel_kwargs must have {n_channels} elements")
        kwargs_per_channel = [dict(kws) for kws in per_channel_kwargs]

    out: List[np.ndarray] = []
    for idx, channel in enumerate(channels):
        out.append(
            np.asarray(seqdist(channel, method=method, **kwargs_per_channel[idx]), dtype=float)
        )
    return out


def _combine_channel_distances(
    channel_distances: List[np.ndarray], weights: List[float]
) -> np.ndarray:
    """Weighted sum of precomputed channel distance matrices."""
    if len(channel_distances) != len(weights):
        raise ValueError("weights length must match number of channel distance matrices")
    if not channel_distances:
        raise ValueError("channel_distances must be non-empty")

    expected_shape = np.asarray(channel_distances[0]).shape
    if len(expected_shape) != 2 or expected_shape[0] != expected_shape[1]:
        raise ValueError("channel distances must have the same square shape")
    D = np.zeros(expected_shape, dtype=float)
    for w, Dc in zip(weights, channel_distances):
        matrix = np.asarray(Dc, dtype=float)
        if matrix.shape != expected_shape:
            raise ValueError("channel distances must have the same square shape")
        D += float(w) * matrix
    return D


def seqMD_gridsearch_weights(
    obj: SeqMD,
    groups: Union[List[int], np.ndarray],
    *,
    method: str = "LCS",
    weight_grid: Optional[List[List[float]]] = None,
    weight_step: float = 0.1,
    metric: str = "silhouette",
    per_channel_kwargs: Optional[List[dict]] = None,
    return_all: bool = False,
) -> Union[List[float], Dict]:
    """Grid search for optimal channel weights in multi-domain sequences.

    Searches over a grid of channel weights to find the combination that
    best separates the given groups according to a clustering quality metric.

    Parameters
    ----------
    obj : SeqMD
        Multi-domain sequence object with multiple channels.
    groups : list or ndarray
        Group labels for each sequence (e.g., cluster assignments).
        Must have length equal to number of sequences.
    method : str, default='LCS'
        Distance method to use for each channel (e.g., 'OM', 'LCS', 'HAM').
    weight_grid : list of lists, optional
        Explicit weight grid to search. Each inner list should have length
        equal to number of channels. If None, generates grid automatically.
    weight_step : float, default=0.1
        Step size for automatic grid generation (ignored if weight_grid provided).
        Must be finite, strictly positive and at most 1. Weights are generated
        from 0 to 1 in steps of weight_step, then normalized.
    metric : str, default='silhouette'
        Quality metric to optimize. Currently only 'silhouette' is supported.
    per_channel_kwargs : list of dict, optional
        Per-channel keyword arguments for distance computation.
    return_all : bool, default=False
        If True, return dict with all weights and scores. If False, return only best weights.

    Returns
    -------
    list of float or dict
        If return_all=False: Best channel weights as list of floats.
        If return_all=True: Dict with keys:
            - 'best_weights': Best channel weights
            - 'best_score': Best metric score
            - 'all_weights': All tested weight combinations
            - 'all_scores': Corresponding metric scores

        Every weight vector is rescaled to sum to the number of channels
        before it is tested; ``best_weights`` and ``all_weights`` are reported
        on that scale, so ``best_weights`` is one of ``all_weights``.

    Raises
    ------
    ValueError
        If groups length doesn't match number of sequences, if groups holds a
        missing value, if weight_step is above 1, or if metric is unsupported.

    Notes
    -----
    The silhouette coefficient measures how well-separated clusters are, with values
    ranging from -1 (poor) to +1 (excellent). Higher values indicate better separation.

    For computational efficiency, consider:
    - Using coarser weight_step (e.g., 0.2 or 0.25)
    - Providing explicit weight_grid with fewer combinations
    - Using faster distance methods (e.g., HAM vs. OM)

    Examples
    --------
    >>> from traminex import seqdef, SeqMD, seqMD_gridsearch_weights
    >>> import numpy as np
    >>> # Create two-channel data
    >>> ch1 = seqdef(np.array([[0,1,2], [0,1,2], [2,1,0], [2,1,0]]))
    >>> ch2 = seqdef(np.array([[0,0,1], [0,0,1], [1,1,0], [1,1,0]]))
    >>> md = SeqMD([ch1, ch2], names=['work', 'family'])
    >>> groups = [0, 0, 1, 1]  # Two groups
    >>> # Search for optimal weights
    >>> best_w = seqMD_gridsearch_weights(md, groups, method='HAM', weight_step=0.25)
    >>> print(best_w)  # the weights sum to the number of channels
    [0.0, 2.0]
    """
    _require_seqmd(obj)
    channels, n_seqs = _validate_seqmd_channels(obj.channels)
    n_channels = len(channels)

    # Validate groups
    groups_arr: np.ndarray = np.asarray(groups)
    if bool(np.any(pd.isna(groups_arr))):
        # A missing label gives no silhouette (TX-AUD-20261019-041).
        raise ValueError("groups must not contain missing values")
    if len(groups_arr) != n_seqs:
        raise ValueError(
            f"groups length ({len(groups_arr)}) must equal number of sequences ({n_seqs})"
        )

    # Check that we have at least 2 groups
    unique_groups = np.unique(groups_arr)
    if len(unique_groups) < 2:
        raise ValueError("groups must contain at least 2 distinct values")

    _validate_metric(metric)

    # Generate or validate weight grid
    if weight_grid is None:
        step = _validate_real(weight_step, name="weight_step", positive=True)
        if step > 1:
            # The grid runs from 0 to 1 in steps of weight_step; a larger step
            # leaves no vector to test (TX-AUD-20261019-041).
            raise ValueError(
                "weight_step must be at most 1: the grid runs from 0 to 1 in steps of weight_step"
            )
        weight_grid = _generate_weight_grid(n_channels, step)
    else:
        try:
            raw_grid = list(weight_grid)
        except TypeError:
            raise TypeError("weight_grid must be a collection of weight vectors") from None
        if not raw_grid:
            raise ValueError("weight_grid must contain at least one weight vector")
        weight_grid = [
            _validate_analytical_weights(w, name="weight_grid entry", n_channels=n_channels)
            for w in raw_grid
        ]

    channel_distances = _precompute_channel_distances(obj, method, per_channel_kwargs)

    # Search over grid
    best_score = -np.inf
    best_weights = None
    all_scores = []
    tested_weights: List[List[float]] = []

    for weights in weight_grid:
        # Normalize weights to sum to n_channels (for interpretability)
        weights_normalized = [float(w) for w in weights]
        weight_sum = float(sum(weights_normalized))
        if weight_sum <= 0:
            raise ValueError("Each weight vector must have positive sum")
        if not np.isfinite(weight_sum):
            largest_weight = max(weights_normalized)
            weights_normalized = [w / largest_weight for w in weights_normalized]
            weight_sum = float(sum(weights_normalized))
        weights_normalized = [w / weight_sum * n_channels for w in weights_normalized]
        tested_weights.append(weights_normalized)

        # Compute weighted multi-domain distance matrix from cached channel distances
        D = _combine_channel_distances(channel_distances, weights_normalized)

        # Compute metric
        score = _compute_clustering_metric(D, groups_arr, metric)
        all_scores.append(score)

        # Track best
        if score > best_score:
            best_score = score
            best_weights = weights_normalized

    if best_weights is None:
        raise RuntimeError("No valid weight vector evaluated")

    if return_all:
        return {
            "best_weights": best_weights,
            "best_score": best_score,
            # The vectors as tested, on best_weights' scale (audit S2-16).
            "all_weights": tested_weights,
            "all_scores": all_scores,
        }
    else:
        return best_weights


def seqMD_optimize_weights(
    obj: SeqMD,
    groups: Union[List[int], np.ndarray],
    *,
    method: str = "LCS",
    metric: str = "silhouette",
    per_channel_kwargs: Optional[List[dict]] = None,
    init_weights: Optional[List[float]] = None,
    max_iter: int = 50,
    tol: float = 0.001,
) -> Dict:
    """Optimize channel weights using coordinate descent.

    Iteratively optimizes channel weights to maximize clustering quality.
    More efficient than grid search for many channels, but may find local optima.

    Parameters
    ----------
    obj : SeqMD
        Multi-domain sequence object with multiple channels.
    groups : list or ndarray
        Group labels for each sequence.
    method : str, default='LCS'
        Distance method to use for each channel.
    metric : str, default='silhouette'
        Quality metric to optimize.
    per_channel_kwargs : list of dict, optional
        Per-channel keyword arguments for distance computation.
    init_weights : list of float, optional
        Initial weights (default: equal weights). Values must be finite,
        non-negative, correctly shaped, and have positive sum. They are
        normalized to sum to the number of channels.
    max_iter : int, default=50
        Positive maximum number of optimization iterations. Booleans and
        non-integral values are rejected.
    tol : float, default=0.001
        Finite, non-negative convergence tolerance for score improvement.

    Returns
    -------
    dict
        Results dictionary with keys:
            - 'weights': Optimized channel weights (Python floats with a positive
              sum: an all-zero vector, whose combined distance is all zeros, is
              never a candidate; TX-AUD-20261011-018)
            - 'score': Final metric score
            - 'iterations': Number of iterations performed
            - 'converged': Whether optimization converged

    Notes
    -----
    Uses a simple coordinate descent approach: iteratively optimizes each
    channel weight while holding others fixed. Stops when improvement falls
    below tolerance or max_iter is reached.

    Examples
    --------
    >>> from traminex import seqdef, SeqMD, seqMD_optimize_weights
    >>> import numpy as np
    >>> ch1 = seqdef(np.array([[0,1,2], [0,1,2], [2,1,0], [2,1,0]]))
    >>> ch2 = seqdef(np.array([[0,0,1], [0,0,1], [1,1,0], [1,1,0]]))
    >>> md = SeqMD([ch1, ch2], names=['work', 'family'])
    >>> result = seqMD_optimize_weights(md, [0, 0, 1, 1], method='HAM')
    >>> print(result['weights'])
    [1.0, 1.0]
    >>> print(f"Silhouette score: {result['score']:.3f}")
    Silhouette score: 1.000
    """
    _require_seqmd(obj)
    channels, n_seqs = _validate_seqmd_channels(obj.channels)
    n_channels = len(channels)

    # Validate groups
    groups_arr: np.ndarray = np.asarray(groups)
    if bool(np.any(pd.isna(groups_arr))):
        # A missing label gives no silhouette (TX-AUD-20261019-041).
        raise ValueError("groups must not contain missing values")
    if len(groups_arr) != n_seqs:
        raise ValueError("groups length must equal number of sequences")
    if len(np.unique(groups_arr)) < 2:
        raise ValueError("groups must contain at least 2 distinct values")

    _validate_metric(metric)
    max_iter = _validate_positive_integer(max_iter, name="max_iter")
    tol = _validate_real(tol, name="tol", positive=False)

    # Initialize weights
    if init_weights is None:
        weights = [1.0] * n_channels
    else:
        weights = _validate_analytical_weights(
            init_weights,
            name="init_weights",
            n_channels=n_channels,
        )

    # Normalize to sum to n_channels
    weight_sum = float(sum(weights))
    if weight_sum <= 0:
        raise ValueError("init_weights must have positive sum")
    if not np.isfinite(weight_sum):
        largest_weight = max(weights)
        weights = [w / largest_weight for w in weights]
        weight_sum = float(sum(weights))
    weights = [w / weight_sum * n_channels for w in weights]

    channel_distances = _precompute_channel_distances(obj, method, per_channel_kwargs)

    # Initial score
    D = _combine_channel_distances(channel_distances, weights)
    score = _compute_clustering_metric(D, groups_arr, metric)
    prev_score = score

    converged = False
    for iteration in range(max_iter):
        improved = False

        # Try optimizing each weight
        for c in range(n_channels):
            best_w = weights[c]
            best_s = score

            # Search range for this weight (while keeping others fixed)
            for w_new in np.linspace(0.0, 2.0 * n_channels, 21):
                weights_try = weights.copy()
                weights_try[c] = float(w_new)
                if sum(weights_try) <= 0:
                    # All-zero weights give an all-zero distance, which
                    # seqMD_distances refuses (TX-AUD-20261011-018).
                    continue

                D = _combine_channel_distances(channel_distances, weights_try)
                s = _compute_clustering_metric(D, groups_arr, metric)

                if s > best_s:
                    best_s = s
                    best_w = w_new
                    improved = True

            weights[c] = float(best_w)
            score = best_s

        # Check convergence
        if not improved or (iteration > 0 and abs(score - prev_score) < tol):
            converged = True
            break

        prev_score = score

    return {
        "weights": [float(weight) for weight in weights],
        "score": float(score),
        "iterations": iteration + 1,
        "converged": converged,
    }


def _generate_weight_grid(n_channels: int, step: float) -> List[List[float]]:
    """Generate normalized weight grid for n channels.

    Generates all combinations of weights from 0 to 1 in steps of `step`,
    then keeps only those that sum to approximately 1 (normalized).

    Parameters
    ----------
    n_channels : int
        Number of channels.
    step : float
        Step size for weight values.

    Returns
    -------
    list of lists
        Weight grid, where each inner list has length n_channels.
    """
    # Generate candidate weight values
    values = np.arange(0.0, 1.0 + step / 2, step)

    # Generate all combinations
    grid = []
    for combo in product(values, repeat=n_channels):
        # Only keep if sum is non-zero
        if sum(combo) > 0:
            # Normalize to sum to 1
            normalized = [float(w / sum(combo)) for w in combo]
            grid.append(normalized)

    # Remove duplicates (due to normalization)
    unique_grid = []
    seen = set()
    for weights in grid:
        key = tuple(round(w, 6) for w in weights)
        if key not in seen:
            seen.add(key)
            unique_grid.append(weights)

    return unique_grid


def _compute_clustering_metric(D: np.ndarray, groups: np.ndarray, metric: str) -> float:
    """Compute clustering quality metric from distance matrix.

    Parameters
    ----------
    D : ndarray
        Symmetric distance matrix (n x n).
    groups : ndarray
        Group labels for each sequence.
    metric : str
        Metric name ('silhouette').

    Returns
    -------
    float
        Metric value (higher is better).
    """
    if metric == "silhouette":
        return _silhouette_from_distance_matrix(D, groups)
    else:
        raise ValueError(f"Unsupported metric: {metric}")


def _silhouette_from_distance_matrix(D: np.ndarray, groups: np.ndarray) -> float:
    """Compute average silhouette coefficient from distance matrix.

    Parameters
    ----------
    D : ndarray
        Symmetric distance matrix (n x n).
    groups : ndarray
        Group labels for each sequence.

    Returns
    -------
    float
        Average silhouette coefficient, range [-1, 1].

    Notes
    -----
    Silhouette coefficient for sample i:
        s(i) = (b(i) - a(i)) / max(a(i), b(i))
    where:
        a(i) = average distance from i to other samples in same cluster
        b(i) = minimum average distance from i to samples in other clusters
    """
    n = len(groups)
    unique_groups = np.unique(groups)
    n_groups = len(unique_groups)

    # Handle edge cases
    if n_groups < 2:
        return 0.0  # Silhouette undefined for single cluster

    silhouettes = []

    for i in range(n):
        group_i = groups[i]

        # a(i): mean distance to other points in same cluster
        same_cluster_mask = groups == group_i
        same_cluster_indices = np.where(same_cluster_mask)[0]

        if len(same_cluster_indices) == 1:
            # The standard silhouette convention defines singleton samples as 0.
            silhouettes.append(0.0)
            continue

        # Exclude self
        distances_same = D[i, same_cluster_mask]
        distances_same = distances_same[same_cluster_indices != i]
        a_i = distances_same.mean() if len(distances_same) > 0 else 0.0

        # b(i): minimum mean distance to other clusters
        b_i = np.inf
        for group_j in unique_groups:
            if group_j == group_i:
                continue
            other_cluster_mask = groups == group_j
            distances_other = D[i, other_cluster_mask]
            mean_dist = distances_other.mean()
            b_i = min(b_i, mean_dist)

        # Silhouette for this sample
        if max(a_i, b_i) > 0:
            s_i = (b_i - a_i) / max(a_i, b_i)
        else:
            s_i = 0.0

        silhouettes.append(s_i)

    return float(np.mean(silhouettes))
