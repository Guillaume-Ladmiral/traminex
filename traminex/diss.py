from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Literal, Optional, Sequence, Union

import numpy as np
import pandas as pd

from ._categorical_order import r_factor
from ._deprecated import _caller_stacklevel
from ._rsum import r_sum

ArrayLike = Union[np.ndarray, Sequence[float]]


def _as_square_matrix(diss: ArrayLike, *, squared: bool = False) -> np.ndarray:
    D = np.asarray(diss, dtype=float)
    if D.ndim != 2 or D.shape[0] != D.shape[1]:
        raise ValueError("diss must be a square dissimilarity matrix")
    if squared:
        D = D**2
    return D


def _resolve_weights(n: int, weights: Optional[ArrayLike]) -> np.ndarray:
    if weights is None:
        return np.ones(n, dtype=float)
    w = np.asarray(weights, dtype=float).reshape(-1)
    if w.shape[0] != n:
        raise ValueError("weights must have one value per observation")
    return w


def _declared_categories(group: Any) -> Optional[List[Any]]:
    """The declared categories of a pandas categorical group, else ``None``."""
    if isinstance(group, pd.Series) and isinstance(group.dtype, pd.CategoricalDtype):
        return list(group.cat.categories)
    if isinstance(group, pd.Categorical):
        return list(group.categories)
    return None


def _factor(
    values: np.ndarray, declared: Optional[List[Any]]
) -> tuple[List[str], np.ndarray, List[Any]]:
    """R's ``factor(values)``: level labels, 0-based codes (-1 missing), each level's first value.

    As :func:`traminex._categorical_order.r_factor`: R's order (text by R's
    collation, numbers by value, a pandas categorical's declared order with
    unused categories dropped), each level named as R's ``as.character`` writes
    it and numbers equal to 15 significant digits one level (TX-AUD-20261019-046).
    """
    return r_factor(values, declared)


def _c_accumulate(total: float, block: np.ndarray) -> float:
    """Continue a C ``double`` accumulator (``total += x``) over ``block`` in order."""
    if block.size == 0:
        return total
    return float(np.add.accumulate(np.concatenate(([total], block)))[-1])


_PRODUCT_BLOCK = 1 << 20  # products formed per vectorised block


def _compiled_pair_sum() -> Optional[Callable[..., tuple[float, float]]]:
    """The compiled :func:`_c_weighted_pair_sum` when Numba is dispatched, else ``None``.

    Read at each call, as ``TRAMINEX_NUMBA`` and the benchmark switch set it.
    """
    from .distances import python_kernels as kernels

    return kernels._c_weighted_pair_sum_nb if kernels._HAVE_NUMBA else None


def _c_weighted_pair_sum(D: np.ndarray, members: np.ndarray, w: np.ndarray) -> tuple[float, float]:
    """``tmrWeightedInertiaDist``'s running sum and weight total, in its C order.

    ``weightedinertia.cpp:17-31``: for each member ``i`` (outer) and each later
    member ``j`` (inner, by position in ``members``), ``result += w[j] * w[i] *
    D[j, i]``, and ``totweights += w[i]``. Each product keeps the kernel's
    association: factoring ``w[i]`` out would not be faithful (for ``w = (Inf, 1,
    1)`` against ``D[:, 0] = (0, 0, 2)`` R forms ``1 * Inf * 0 = NaN``), and a
    diagonal term never enters the sum, so an infinite weight forms no kept
    ``Inf * 0``. The products of a block of ``i`` columns are gathered with a
    boolean mask, which yields them in row-major order of the transposed block
    (``i`` outer, ``j`` inner), and ``np.add.accumulate`` adds them one by one,
    as the C loop does; a pairwise ``np.sum`` differs in the last bits. With Numba
    dispatched, the same loop runs compiled, with the same bytes; it was most of
    ``dissassoc``'s permutation test (TX-AUD-20261010-083).
    """
    compiled = _compiled_pair_sum()
    if compiled is not None:
        result, totweights = compiled(
            np.asarray(D, dtype=np.float64),
            np.asarray(members, dtype=np.int64),
            np.asarray(w, dtype=np.float64),
        )
        return float(result), float(totweights)
    size = members.size
    ws = w[members]
    result = 0.0
    step = max(1, _PRODUCT_BLOCK // max(size, 1))
    for start in range(0, size - 1, step):
        stop = min(start + step, size - 1)
        block = (ws[start:, None] * ws[None, start:stop]) * D[
            np.ix_(members[start:], members[start:stop])
        ]
        later = np.arange(size - start)[:, None] > np.arange(stop - start)[None, :]
        result = _c_accumulate(result, block.T[later.T])
    return result, _c_accumulate(0.0, ws)


def _c_weighted_inertia_dist(D: np.ndarray, members: np.ndarray, w: np.ndarray) -> float:
    """``tmrWeightedInertiaDist(var=FALSE)`` in its C accumulation order.

    :func:`_c_weighted_pair_sum` (``weightedinertia.cpp:17-34``); a positive
    total divides the result. ``members`` is not checked: the random-sampling
    test draws cases with replacement, which R's kernel takes as they come.
    """
    result, totweights = _c_weighted_pair_sum(D, members, w)
    return result / totweights if totweights > 0 else result


def _case_indices(indices: Optional[ArrayLike], n: int) -> np.ndarray:
    """Zero-based case positions for the public weighted-inertia helpers (all cases by default).

    R's kernels take integer case numbers. A non-integer array (such as weights
    given in the indices' place), a boolean mask, a negative or out-of-range
    position or a repeated one is refused (TX-AUD-20261011-036).
    """
    if indices is None:
        return np.arange(n)
    idx = np.asarray(indices).reshape(-1)
    if idx.size == 0:
        return np.zeros(0, dtype=int)
    if idx.dtype == np.bool_ or not np.issubdtype(idx.dtype, np.integer):
        raise TypeError("indices must be integer case positions (zero-based)")
    if int(idx.min()) < 0 or int(idx.max()) >= n:
        raise ValueError(f"indices must lie between 0 and {n - 1}")
    seen = np.zeros(n, dtype=bool)
    seen[idx] = True
    if int(np.count_nonzero(seen)) != idx.size:
        raise ValueError("indices must not repeat a case")
    return idx.astype(int)


def weighted_pairwise_sum(
    diss: ArrayLike,
    indices: Optional[ArrayLike] = None,
    weights: Optional[ArrayLike] = None,
) -> float:
    """Weighted sum of pairwise dissimilarities over a subset.

    Computes ``sum_{i<j in indices} w_i * w_j * D[i, j]`` (strictly off-diagonal,
    like TraMineR's C routines). This is the raw sum-of-squares building block of
    R TraMineR's weighted discrepancy analysis (``tmrWeightedInertiaDist`` before
    its normalization step).

    Parameters
    ----------
    diss : array-like
        Square dissimilarity matrix of shape (n, n).
    indices : array-like of int, optional
        Subset of observation indices (0-based, in range, not repeated).
        Defaults to all observations.
    weights : array-like, optional
        Per-observation weights of length n. Defaults to 1.
    """
    D = _as_square_matrix(diss)
    n = D.shape[0]
    w = _resolve_weights(n, weights)
    idx = _case_indices(indices, n)
    if idx.size <= 1:
        return 0.0
    # Non-finite case weights are R-faithfully propagated rather than refused
    # (see :func:`weighted_inertia`), so NaN/Inf arithmetic here is intended.
    # ``D.T[j, i]`` is ``D[i, j]``: the upper triangle, summed in R's order.
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        return _c_weighted_pair_sum(D.T, idx, w)[0]


def weighted_inertia(
    diss: ArrayLike,
    indices: Optional[ArrayLike] = None,
    weights: Optional[ArrayLike] = None,
) -> float:
    """Weighted inertia (sum of squares, Ward sense) of a subset of cases.

    R-faithful port of TraMineR's ``tmrWeightedInertiaDist(..., var=FALSE)``:
    ``weighted_pairwise_sum(D.T, indices, weights) / sum(weights[indices])``, with
    ``indices`` zero-based rows of ``diss`` (all rows by default; integers, in
    range and not repeated).
    The transpose follows R's lower-triangle access; for a symmetric matrix
    it leaves the value unchanged. Both sums run in the kernel's order (cases
    in the order of ``indices``), so the result equals R's to the last bit.

    The C kernel guards that division with ``if (totweights>0)`` and otherwise
    returns the **undivided** sum
    (``SOURCES/TraMineR/src/weightedinertia.cpp:34``), so a NaN, all-negative or
    zero-sum weight vector propagates rather than collapsing to a finite ``0.0``.
    Neither R nor traminex validates weights on this path. This is the SS
    quantity used by R ``dissassoc`` / ``disstree`` for group residuals.
    """
    D = _as_square_matrix(diss)
    n = D.shape[0]
    w = _resolve_weights(n, weights)
    # Summing non-finite weights is intended arithmetic on this path -- `Inf +
    # -Inf` is a legitimate NaN that R propagates -- so it must not raise a
    # numpy warning, matching the other weighted entry points of this module.
    idx = _case_indices(indices, n)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        return _c_weighted_inertia_dist(D, idx, w)


def dissvar(
    diss: ArrayLike,
    weights: Optional[ArrayLike] = None,
    squared: bool = False,
) -> float:
    """Discrepancy (pseudo-variance) of a dissimilarity matrix — R-faithful.

    Port of TraMineR 2.2-14 ``dissvar(diss, weights=, squared=)`` (Studer &
    Ritschard's discrepancy analysis):

    - unweighted: ``sum(diss) / (2 * n^2)``
    - weighted: ``tmrWeightedInertiaDist(..., var=TRUE)``, i.e.
      ``sum_{i<j} w_i w_j d_ij / W^2`` with ``W = sum(w)``

    Anchored to R by ``tests/test_r_goldens_dissvar.py`` (golden
    ``dissvar_phase2``). The previous traminex ``dissvar`` (variance of the
    lower-triangle distances, unweighted) was a simplified stand-in and is kept
    as :func:`dissvar_summary`.

    Parameters
    ----------
    diss : array-like
        Square dissimilarity matrix of shape (n, n).
    weights : array-like, optional
        Per-observation (e.g. survey) weights of length n.
    squared : bool, default False
        If True, compute on squared dissimilarities (R ``squared=TRUE``).
    """
    D = _as_square_matrix(diss, squared=squared)
    n = D.shape[0]
    if weights is None:
        # R's unweighted shortcut: sum(diss)/(2*nrow(diss)^2).
        return float(D.sum() / (2.0 * n * n))
    w = _resolve_weights(n, weights)
    # Non-finite case weights are R-faithfully propagated rather than refused
    # (see :func:`weighted_inertia`), so summing them is intended arithmetic --
    # `Inf + -Inf` is a legitimate NaN here -- and must not raise a numpy
    # warning, exactly as in the other weighted entry points of this module.
    with np.errstate(divide="ignore", invalid="ignore"):
        totw = _c_accumulate(0.0, w)  # the kernel's running total
        inertia = weighted_inertia(D, None, w)
        # Second guarded division of the C kernel's `var=TRUE` branch
        # (SOURCES/TraMineR/src/weightedinertia.cpp:36-38): when the total
        # weight is not strictly positive R returns the value undivided, so NaN
        # and the raw negative/zero-sum results survive rather than becoming a
        # finite 0.0.
        return inertia / totw if totw > 0 else inertia


def dissvar_grp(
    diss: ArrayLike,
    group: Optional[ArrayLike] = None,
    *,
    weights: Optional[ArrayLike] = None,
    squared: bool = False,
) -> Dict[str, float]:
    """Per-group discrepancy (pseudo-variance) — R-faithful.

    Port of TraMineRextras 0.6.9 ``dissvar.grp(diss, group, ...)``: the weighted
    discrepancy :func:`dissvar` computed within each group level over
    ``diss[ig, ig]``. (R's ``dissvar.grp`` forwards ``...`` to ``dissvar``
    without subsetting case weights; traminex subsets ``weights`` per group so
    the weighted path is usable.) Anchored to live R by
    ``tests/test_r_goldens_diss_mfa.py``.

    Parameters
    ----------
    diss : array-like
        Square dissimilarity matrix of shape (n, n).
    group : array-like, optional
        Group membership, one value per observation. ``None`` is one group of
        all cases, named ``"1"`` (R's ``group = NULL``, dissvar-grp.R:13-14).
    weights : array-like, optional
        Per-observation weights of length n (subset per group).
    squared : bool, default False
        Compute on squared dissimilarities.

    Returns
    -------
    dict
        ``{group_level: dissvar}`` for each group level, sorted, or in the
        declared category order for a pandas categorical (as R's ``factor``);
        each level is named as R names it (``"1"``, not ``"1.0"``).
    """
    D = _as_square_matrix(diss)
    n = D.shape[0]
    if group is None:
        group = np.ones(n)
    declared = _declared_categories(group)
    group = np.asarray(group)
    if group.shape[0] != n:
        raise ValueError("length(group) not compatible with size of diss")
    w = None if weights is None else np.asarray(weights, dtype=float)
    out: Dict[str, float] = {}
    labels, codes, _ = _factor(group, declared)
    for code, label in enumerate(labels):
        ig = np.flatnonzero(codes == code)
        sub = D[np.ix_(ig, ig)]
        wsub = None if w is None else w[ig]
        out[label] = float(dissvar(sub, weights=wsub, squared=squared))
    return out


def _weighted_inertia_contrib(D: np.ndarray, idx: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Port of TraMineR ``tmrWeightedInertiaContrib``, in its C order.

    For each case ``i`` in ``idx``: ``sum_{j in idx, j != i} w_j * D[., .] /
    sum_{j in idx} w_j``. ``weightedinertia.cpp:68-101`` adds case ``i``'s terms
    in index order to a running double from ``+0.0``: ``D[i, j] * w[j]`` for an
    earlier ``j`` and ``D[j, i] * w[j]`` for a later one. So identical cases tie
    exactly and the first is the medoid; a pairwise sum breaks such ties in the
    last bit (TX-AUD-20261011-023). Each row accumulates with
    ``np.add.accumulate`` from a ``+0.0`` column, the diagonal term set to
    ``+0.0``: a running total from ``+0.0`` is never ``-0.0``, so adding it is
    exact, as skipping it is in C.
    """
    idx = np.asarray(idx, dtype=int).reshape(-1)
    size = idx.size
    ws = w[idx]
    contrib = np.zeros(size, dtype=float)
    # Non-finite case weights are R-faithfully propagated here rather than
    # refused (see :func:`weighted_inertia`), so the resulting NaN/Inf arithmetic
    # is intended and must not raise a numpy warning. The total-weight sum is
    # inside the guard for the same reason: `Inf + -Inf` is a legitimate NaN.
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        step = max(1, _PRODUCT_BLOCK // max(size, 1))
        for start in range(0, size, step):
            rows = np.arange(start, min(start + step, size))
            # D[i, j] for j before i, D[j, i] for j after it (column-major dmat reads).
            earlier = np.arange(size)[None, :] < rows[:, None]
            terms = np.where(earlier, D[np.ix_(idx[rows], idx)], D[np.ix_(idx, idx[rows])].T)
            terms = terms * ws[None, :]
            terms[np.arange(rows.size), rows] = 0.0
            running = np.add.accumulate(
                np.concatenate((np.zeros((rows.size, 1)), terms), axis=1), axis=1
            )
            contrib[rows] = running[:, -1]
        totw = _c_accumulate(0.0, ws)
        if totw > 0:
            contrib = contrib / totw
    return contrib


def _r_weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    """R's ``weighted.mean(x, w)``: ``sum((x * w)[w != 0]) / sum(w)`` with R's long-double
    ``sum()`` (batch-12 review M-4), unguarded: a weight total of 0 gives NaN or Inf."""
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        return float(
            np.float64(r_sum((values * weights)[weights != 0])) / np.float64(r_sum(weights))
        )


def _disscenter_values(D: np.ndarray, grpint: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Distance of each case to its (weighted) group centre.

    Port of R ``disscenter(diss, group, weights)`` (trim=0): per group, the
    weighted inertia contribution minus half the weighted mean contribution.
    """
    n = D.shape[0]
    out = np.zeros(n, dtype=float)
    for g in np.unique(grpint):
        idx = np.flatnonzero(grpint == g)
        dc = _weighted_inertia_contrib(D, idx, w)
        wg = w[idx]
        # R's `dc - weighted.mean(dc, weights[cond])/2`
        # (SOURCES/TraMineR/R/disscenter.R:87) is unconditional, unlike the
        # guarded division inside the C contribution kernel.
        dc = dc - _r_weighted_mean(dc, wg) / 2.0
        out[idx] = dc
    return out


def _weighted_residuals_z(vals: np.ndarray, w: np.ndarray) -> float:
    """R ``WeightedResidualsZ``: sum_i w_i * (x_i - weighted.mean(x, w))^2.

    R's ``weighted.mean`` is ``sum((x * w)[w != 0]) / sum(w)``, without a guard: a weight total of 0 gives NaN or Inf,
    and a negative total is used as it is (dissassocweighted.R:4-6).
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        product = vals * w
        mean = np.sum(product[w != 0]) / np.sum(w)
        return float(np.sum(w * (vals - mean) ** 2))


def _unweighted_residuals_z(vals: np.ndarray) -> float:
    """R ``UnweightedResidualsZ``: sum_i (x_i - mean(x))^2."""
    if vals.size == 0:
        return 0.0
    return float(np.sum((vals - vals.mean()) ** 2))


_DISSASSOC_STAT_NAMES = ["Pseudo F", "Pseudo Fbf", "Pseudo R2", "Bartlett", "Levene"]


def _compute_weighted_test_values(
    group_terms: Sequence[tuple[float, float, float]],
    *,
    sc_tot: float,
    totweights: float,
    k: int,
    dissc_sc_tot: float,
) -> np.ndarray:
    """Port of R ``ComputeWeightedTestValues`` (dissassocweighted.R).

    ``group_terms`` holds per-group ``(ni, SCresi, zi)`` (weighted size,
    residual inertia, Levene residual sum).
    """
    # Use numpy float semantics throughout so degenerate cases yield inf/nan
    # exactly as R does, instead of raising ZeroDivisionError.
    sc_tot = np.float64(sc_tot)
    totweights = np.float64(totweights)
    k_float = np.float64(k)
    dissc_sc_tot = np.float64(dissc_sc_tot)
    sc_res = np.float64(0.0)
    lns = np.float64(0.0)
    nlnvi = np.float64(0.0)
    fbf_denom = np.float64(0.0)
    sumz = np.float64(0.0)
    s1ni = np.float64(0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        for ni, sc_resi, zi in group_terms:
            ni = np.float64(ni)
            sc_resi = np.float64(sc_resi)
            sumz += np.float64(zi)
            vari = sc_resi / ni
            lns += (ni - 1.0) * (vari / (totweights - k_float))
            nlnvi += (ni - 1.0) * np.log(vari)
            sc_res += sc_resi
            fbf_denom += (1.0 - ni / totweights) * vari
            s1ni += 1.0 / (ni - 1.0)
        c_calc = 1.0 + 1.0 / (3.0 * (k_float - 1.0)) * (s1ni - 1.0 / (totweights - k_float))
        t_calc = (totweights - k_float) * np.log(lns) - nlnvi
        bartlett = t_calc / c_calc
        sc_exp = sc_tot - sc_res
        pseudo_r2 = sc_exp / sc_tot
        pseudo_f = (sc_exp / (k_float - 1.0)) / (sc_res / (totweights - k_float))
        fpf = sc_exp / fbf_denom
        pseudo_w = ((dissc_sc_tot - sumz) / (k_float - 1.0)) / (sumz / (totweights - k_float))
    return np.array([pseudo_f, fpf, pseudo_r2, bartlett, pseudo_w], dtype=float)


def _permutation_pvalues(t: np.ndarray, t0: np.ndarray) -> np.ndarray:
    """R's ``sum(t[, i] >= t0[i]) / R`` for every column i.

    ``t`` holds all ``R`` rows, the observed statistics first. In R a NaN
    anywhere in a column, observed or permuted, makes the sum -- and so the
    p-value -- NA (``TraMineR-permutation.R:21``); NumPy's ``>=`` would
    silently count it as ``False`` instead.
    """
    exceed = (t >= t0[np.newaxis, :]).sum(axis=0) / float(t.shape[0])
    return np.where(np.isnan(t).any(axis=0), np.nan, exceed)


def _tramineR_permutation(
    n_items: int,
    R: int,
    statistic: Callable[[np.ndarray], np.ndarray],
    rng: np.random.Generator,
) -> Dict[str, Any]:
    """Port of R ``TraMineR.permutation``: t[0]=t0, pval = mean(t >= t0)."""
    identity = np.arange(n_items, dtype=int)
    t0 = np.asarray(statistic(identity), dtype=float)
    R = int(R)
    out: Dict[str, Any] = {"R": R, "t0": t0}
    if R > 1:
        t = np.empty((R, t0.size), dtype=float)
        t[0] = t0
        for r in range(1, R):
            t[r] = statistic(rng.permutation(n_items))
        out["t"] = t
        out["pval"] = _permutation_pvalues(t, t0)
    else:
        out["t"] = None
        out["pval"] = np.full(t0.size, np.nan)
    return out


@dataclass
class DissAssocResult:
    """Result of :func:`dissassoc` (mirrors R TraMineR's ``dissassoc`` object).

    Attributes
    ----------
    groups : pandas.DataFrame
        Per-group weighted size ``n`` and ``discrepancy``, plus a Total row.
    anova_table : pandas.DataFrame
        Pseudo ANOVA table (SS, df, MSE; rows Exp/Res/Total).
    stat : pandas.DataFrame
        Observed statistics ``t0`` and permutation ``p_value`` for
        Pseudo F, Pseudo Fbf, Pseudo R2, Bartlett and Levene.
    perms : dict
        Raw permutation output (``t0``, ``t``, ``pval``, ``R``).
    """

    groups: pd.DataFrame
    anova_table: pd.DataFrame
    stat: pd.DataFrame
    R: int
    weight_permutation: str
    perms: Dict[str, Any] = field(repr=False, default_factory=dict)

    def __str__(self) -> str:  # pragma: no cover - convenience output
        return (
            "Pseudo ANOVA table:\n"
            f"{self.anova_table}\n\n"
            f"Test values (p-values based on {self.R} permutations):\n"
            f"{self.stat}\n\n"
            "Discrepancy per level:\n"
            f"{self.groups}"
        )


def dissassoc(
    diss: ArrayLike,
    group: ArrayLike,
    weights: Optional[ArrayLike] = None,
    R: int = 1000,
    weight_permutation: Literal[
        "none",
        "diss",
        "replicate",
        "rounded-replicate",
        "group",
        "permute",
        "random-sampling",
    ] = "replicate",
    squared: bool = False,
    random_state: Optional[Union[int, np.random.Generator]] = None,
) -> DissAssocResult:
    """Discrepancy analysis of a dissimilarity matrix by group — R-faithful.

    Port of TraMineR 2.2-14 ``dissassoc(diss, group, weights=, R=,
    weight.permutation=, squared=)`` (Studer & Ritschard): the weighted pseudo
    ANOVA (pseudo F, pseudo R²), the Bartlett and Levene homogeneity tests,
    and the permutation test of the association.

    Observed statistics are anchored to live R by
    ``tests/test_r_goldens_dissassoc.py``; p-values are permutation-based
    (pval = proportion of permutations, including the observed one, with
    ``t >= t0`` — R's convention; NaN, R's NA, when the statistic is NaN in
    the observed data or in any permutation, or when ``R <= 1``) and
    therefore stochastic; pass ``random_state`` for reproducibility (Python
    addition, R relies on the global RNG).

    Supported ``weight_permutation`` variants are ``"none"``, ``"diss"``,
    ``"replicate"``, ``"rounded-replicate"``, ``"group"``, ``"permute"``,
    and ``"random-sampling"``. Without ``weights`` the unweighted test is used
    whatever the value (as in R); ``"none"`` with ``weights`` raises
    ``ValueError``, since R has no test for that combination and stops. R and
    NumPy have different random streams, so seeded sampled values are
    reproducible within Python while deterministic sufficient statistics are
    cross-language exact.

    Group levels (rows of ``groups``) are sorted, or taken in the declared
    category order of a pandas categorical, unused categories dropped, as
    R's ``factor``.

    The previous simplified traminex ``dissassoc`` (between/within effect
    size) is kept as :func:`dissassoc_effect`.
    """
    D = _as_square_matrix(diss)
    declared = _declared_categories(group)
    group = np.asarray(group)
    if group.shape[0] != D.shape[0]:
        raise ValueError("group must have one value per observation")

    # Remove missing group values (R: diss[!is.na(group), !is.na(group)]).
    valid = ~pd.isna(group)
    if not np.all(valid):
        keep = np.flatnonzero(valid)
        D = D[np.ix_(keep, keep)]
        group = group[keep]
        if weights is not None:
            weights = np.asarray(weights, dtype=float)[keep]
    n = D.shape[0]
    if squared:
        D = D**2

    use_replicate = weight_permutation in ("replicate", "rounded-replicate")
    unweighted = weights is None
    if unweighted:
        w = np.ones(n, dtype=float)
        use_replicate = False
        weight_permutation = "none"
    else:
        if weight_permutation == "none":
            # No branch of R's dissassocweighted handles weights with "none"
            # (dissassocweighted.R:342-359), so R stops; computing on would mix
            # weighted sums of squares with unweighted group sizes.
            raise ValueError(
                "weight_permutation='none' is the unweighted test and cannot be "
                "combined with weights; use 'replicate', 'rounded-replicate', "
                "'diss', 'group', 'permute' or 'random-sampling', or drop weights"
            )
        w = _resolve_weights(n, weights)
    if use_replicate:
        round_error = float(np.abs(np.round(w) - w).sum())
        if round_error > 0:
            if weight_permutation == "replicate":
                raise ValueError("to permute replicates, weights must be integers")
            w = np.round(w)
        w = w.astype(float)
    allowed_methods = {
        "none",
        "diss",
        "replicate",
        "rounded-replicate",
        "group",
        "permute",
        "random-sampling",
    }
    if weight_permutation not in allowed_methods:
        raise ValueError(
            "weight_permutation must be one of none, diss, replicate, "
            "rounded-replicate, group, permute, or random-sampling"
        )
    if weight_permutation == "random-sampling":
        if np.any(w < 0):
            raise ValueError("random-sampling weights must be nonnegative")
        weight_sum = float(w.sum())
        if weight_sum <= 0:
            raise ValueError("random-sampling weights must have a positive sum")
        w = w * n / weight_sum

    # R: grp <- factor(group[!is.na(group)]) (dissassocweighted.R:294) keeps a
    # factor's level order; plain values are sorted.
    labels, grpint, _ = _factor(group, declared)
    k = len(labels)
    totweights = float(w.sum())
    sc_tot = weighted_inertia(D, None, w)

    # Per-group discrepancy table.
    group_n = np.zeros(k + 1, dtype=float)
    group_disc = np.zeros(k + 1, dtype=float)
    sc_res = 0.0
    for g in range(k):
        idx = np.flatnonzero(grpint == g)
        ng = float(w[idx].sum())
        group_n[g] = ng
        # R computes every group whose weight total is not zero, negative totals included (dissassocweighted.R:313-330).
        if ng != 0:
            r = weighted_inertia(D, idx, w)
            group_disc[g] = r / ng
            sc_res += r
    group_n[k] = totweights
    # R divides in floating point: a zero total weight gives NaN (0/0) or an
    # infinity, where Python's float division raised (TX-AUD-20261010-096).
    with np.errstate(divide="ignore", invalid="ignore"):
        group_disc[k] = np.float64(sc_tot) / np.float64(totweights)
    groups_df = pd.DataFrame(
        {"n": group_n, "discrepancy": group_disc},
        index=labels + ["Total"],
    )
    ss = np.array([sc_tot - sc_res, sc_res, sc_tot], dtype=float)
    df_ = np.array([k - 1.0, totweights - k, totweights - 1.0], dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        mse = ss / df_  # R yields Inf/NaN on degenerate groupings (k == 1)
    anova = pd.DataFrame({"SS": ss, "df": df_, "MSE": mse}, index=["Exp", "Res", "Total"])

    rng = (
        random_state
        if isinstance(random_state, np.random.Generator)
        else np.random.default_rng(random_state)
    )
    indgrp: List[np.ndarray] = [np.flatnonzero(grpint == g) for g in range(k)]

    if weight_permutation == "none":
        dissc = _disscenter_values(D, grpint, w)
        dissc_sc_tot = _unweighted_residuals_z(dissc)

        def statistic(ind: np.ndarray) -> np.ndarray:
            terms = []
            for g in range(k):
                members = ind[indgrp[g]]
                ni = float(members.size)
                sc_resi = weighted_inertia(D, members, w)
                zi = _unweighted_residuals_z(dissc[members])
                terms.append((ni, sc_resi, zi))
            return _compute_weighted_test_values(
                terms, sc_tot=sc_tot, totweights=totweights, k=k, dissc_sc_tot=dissc_sc_tot
            )

        perms = _tramineR_permutation(n, R, statistic, rng)
    elif weight_permutation == "diss":
        dissc = _disscenter_values(D, grpint, w)
        dissc_sc_tot = _weighted_residuals_z(dissc, w)

        def statistic(ind: np.ndarray) -> np.ndarray:
            terms = []
            for g in range(k):
                members = ind[indgrp[g]]
                ni = float(w[indgrp[g]].sum())
                sc_resi = weighted_inertia(D, members, w)
                # The Levene residuals are taken at the *permuted* positions
                # (dissassocweighted.R:119-125: WeightedResidualsZ(dissc[wind],
                # ww[wind]) with wind the permuted members). Only the "group"
                # branch below uses the unpermuted positions (lines 96-98).
                zi = _weighted_residuals_z(dissc[members], w[members])
                terms.append((ni, sc_resi, zi))
            return _compute_weighted_test_values(
                terms, sc_tot=sc_tot, totweights=totweights, k=k, dissc_sc_tot=dissc_sc_tot
            )

        perms = _tramineR_permutation(n, R, statistic, rng)
    elif weight_permutation in ("group", "permute"):
        compute_t0 = True

        def statistic(ind: np.ndarray) -> np.ndarray:
            nonlocal compute_t0
            if weight_permutation == "permute":
                if compute_t0:
                    ww = w
                    compute_t0 = False
                else:
                    ww = w[rng.permutation(n)]
            else:
                ww = w[ind]
            sc_tot_p = weighted_inertia(D, None, ww)
            dissc_p = _disscenter_values(D, grpint[ind], ww)
            dissc_sc_tot_p = _weighted_residuals_z(dissc_p, ww)
            terms = []
            for g in range(k):
                members = ind[indgrp[g]]
                ni = float(ww[indgrp[g]].sum())
                sc_resi = weighted_inertia(D, members, ww)
                zi = _weighted_residuals_z(dissc_p[indgrp[g]], ww[indgrp[g]])
                terms.append((ni, sc_resi, zi))
            return _compute_weighted_test_values(
                terms, sc_tot=sc_tot_p, totweights=totweights, k=k, dissc_sc_tot=dissc_sc_tot_p
            )

        perms = _tramineR_permutation(n, R, statistic, rng)
    elif weight_permutation == "random-sampling":
        # R uses the ordinary weighted-dissimilarity statistic as t0, then two
        # independent probability-weighted bootstrap samples for every null
        # draw: one supplies distance objects and one supplies group profiles.
        dissc_t0 = _disscenter_values(D, grpint, w)
        dissc_sc_tot_t0 = _weighted_residuals_z(dissc_t0, w)
        t0_terms = []
        for g in range(k):
            members = indgrp[g]
            ni = float(w[members].sum())
            sc_resi = weighted_inertia(D, members, w)
            zi = _weighted_residuals_z(dissc_t0[members], w[members])
            t0_terms.append((ni, sc_resi, zi))
        t0 = _compute_weighted_test_values(
            t0_terms,
            sc_tot=sc_tot,
            totweights=totweights,
            k=k,
            dissc_sc_tot=dissc_sc_tot_t0,
        )

        unit_weights = np.ones(n, dtype=float)

        def sampled_statistic(diss_indices: np.ndarray, group_indices: np.ndarray) -> np.ndarray:
            sampled_groups = grpint[group_indices]
            # Drawn with replacement: the kernel takes the repeated cases as R's does.
            with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
                sc_tot_p = _c_weighted_inertia_dist(D, diss_indices, unit_weights)
            dissc_p = np.zeros(n, dtype=float)
            for g in range(k):
                condition = sampled_groups == g
                group_objects = diss_indices[condition]
                if group_objects.size:
                    dc = _weighted_inertia_contrib(D, group_objects, unit_weights)
                    dc -= float(dc.mean()) / 2.0
                    dissc_p[group_objects] = dc
            dissc_sc_tot_p = _unweighted_residuals_z(dissc_p)
            terms = []
            for g in range(k):
                condition = sampled_groups == g
                members = np.sort(diss_indices[condition])
                ni = float(members.size)
                with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
                    sc_resi = _c_weighted_inertia_dist(D, members, unit_weights)
                zi = _unweighted_residuals_z(dissc_p[condition])
                terms.append((ni, sc_resi, zi))
            return _compute_weighted_test_values(
                terms,
                sc_tot=sc_tot_p,
                totweights=totweights,
                k=k,
                dissc_sc_tot=dissc_sc_tot_p,
            )

        perms = {"R": int(R), "t0": t0}
        if R > 1:
            probability = w / float(w.sum())
            sampled = np.empty((int(R), t0.size), dtype=float)
            sampled[0] = t0
            for permutation in range(1, int(R)):
                diss_indices = np.asarray(
                    np.sort(rng.choice(n, size=n, replace=True, p=probability)), dtype=int
                )
                group_indices = np.asarray(
                    rng.choice(n, size=n, replace=True, p=probability), dtype=int
                )
                sampled[permutation] = sampled_statistic(diss_indices, group_indices)
            perms["t"] = sampled
            perms["pval"] = _permutation_pvalues(sampled, t0)
        else:
            perms["t"] = None
            perms["pval"] = np.full(t0.size, np.nan)
    else:  # replicate / rounded-replicate
        wint = w.astype(int)
        indrep = np.repeat(np.arange(n, dtype=int), wint)
        grouprep = np.repeat(grpint, wint)
        indgrp_rep: List[np.ndarray] = [np.flatnonzero(grouprep == g) for g in range(k)]
        dissc = _disscenter_values(D, grpint, w)
        dissc_sc_tot = _weighted_residuals_z(dissc, w)

        def statistic(ind: np.ndarray) -> np.ndarray:
            terms = []
            for g in range(k):
                nind = indrep[ind[indgrp_rep[g]]]
                wwt = np.bincount(nind, minlength=n).astype(float)
                wind = np.flatnonzero(wwt > 0)
                ni = float(nind.size)
                sc_resi = weighted_inertia(D, wind, wwt)
                zi = _weighted_residuals_z(dissc[wind], wwt[wind])
                terms.append((ni, sc_resi, zi))
            return _compute_weighted_test_values(
                terms, sc_tot=sc_tot, totweights=totweights, k=k, dissc_sc_tot=dissc_sc_tot
            )

        perms = _tramineR_permutation(indrep.size, R, statistic, rng)

    stat = pd.DataFrame(
        {"t0": perms["t0"], "p_value": perms["pval"]},
        index=_DISSASSOC_STAT_NAMES,
    )
    return DissAssocResult(
        groups=groups_df,
        anova_table=anova,
        stat=stat,
        R=int(R),
        weight_permutation=weight_permutation,
        perms=perms,
    )


def double_center_matrix(D: np.ndarray) -> np.ndarray:
    """Double-center a distance matrix (Gower centering approximation).

    .. note:: This helper was exported as ``disscenter`` before traminex
       0.0.6 — a **name clash** with R TraMineR's ``disscenter`` (distance of
       each case to its group centre), which is a different operation and now
       lives at :func:`disscenter`. See DUPLICATES.md A7.

    Parameters
    ----------
    D : np.ndarray
        Square distance matrix of shape (n, n)

    Returns
    -------
    np.ndarray
        Centered matrix C = D - row_mean - col_mean + grand_mean

    Notes
    -----
    This is a simplified centering operation used for visualization and
    exploratory analysis. For classical MDS, use `gower_matrix()` instead.
    """
    D = np.asarray(D, dtype=float)
    r = D.mean(axis=1, keepdims=True)
    c = D.mean(axis=0, keepdims=True)
    g = D.mean()
    return D - r - c + g


def disscenter_double_center(D: np.ndarray) -> np.ndarray:
    """Deprecated alias for :func:`double_center_matrix`.

    Introduced in 0.0.6 as a migration target for code that used the old
    ``disscenter`` (matrix double-centering); ``disscenter`` itself now has
    the R TraMineR meaning (distance of each case to its group centre). It
    stays exported while the public API is frozen at the 226-name 0.2.0
    surface (ECOSYSTEM.md pre-release sequencing); removal is an R6 decision.
    """
    import warnings

    warnings.warn(
        "disscenter_double_center is a deprecated alias; use double_center_matrix. "
        "Note that traminex.disscenter now means R TraMineR's disscenter "
        "(distance of each case to its group centre).",
        DeprecationWarning,
        stacklevel=_caller_stacklevel(),
    )
    return double_center_matrix(D)


def _weighted_inertia_contrib_ext(
    D: np.ndarray, idx: np.ndarray, others: np.ndarray, w: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Port of TraMineR ``tmrWeightedInertiaContribExt``.

    Returns the weighted mean distance to the group ``idx`` for the group
    members (excluding self) and for the ``others`` (non-members):
    ``sum_{j in idx} w_j * D[i, j] / sum_{j in idx} w_j``.

    The member half is guarded by ``if (totweights>0)`` in the C kernel; the
    outsider half is **not**, and the C kernel forms it as
    ``sum_j w_i w_j d_ij / (w_i * totweights)``
    (``SOURCES/TraMineR/src/weightedinertia.cpp:148-150``). The outsider's own
    weight ``w_i`` cancels algebraically but not numerically, so a zero- or
    non-finite-weight outsider yields NaN in R. That is reproduced here rather
    than cancelled away.
    """
    members = _weighted_inertia_contrib(D, idx, w)
    others = np.asarray(others, dtype=int).reshape(-1)
    idx = np.asarray(idx, dtype=int).reshape(-1)
    wg = w[idx]
    w_others = w[others]
    with np.errstate(divide="ignore", invalid="ignore"):
        # The group-weight sum belongs inside the guard for the same reason as
        # the arithmetic below: `Inf + -Inf` is a legitimate NaN that R keeps.
        totw = float(wg.sum())
        ext = (D[np.ix_(others, idx)] @ wg) * w_others / (w_others * totw)
    return members, ext


def disscenter(
    diss: ArrayLike,
    group: Optional[ArrayLike] = None,
    medoids_index: Optional[Union[str, bool]] = None,
    allcenter: bool = False,
    weights: Optional[ArrayLike] = None,
    squared: bool = False,
):
    """Distance of each case to its (weighted) group centre — R-faithful.

    Port of TraMineR 2.2-14 ``disscenter(diss, group=, medoids.index=,
    allcenter=, weights=, squared=)`` (Studer & Ritschard's discrepancy
    analysis). Anchored to live R by ``tests/test_r_goldens_disscenter.py``.

    .. warning:: **Name-clash resolution (DUPLICATES.md A7).** Before
       traminex 0.0.6, ``disscenter`` double-centred a matrix — a different
       operation. That helper is now :func:`double_center_matrix`.

    Parameters
    ----------
    diss : array-like
        Square dissimilarity matrix of shape (n, n).
    group : array-like, optional
        Group membership (one value per case). Default: a single group.
        Levels are sorted, or taken in the declared category order of a
        pandas categorical (unused categories dropped), as R's ``factor``.

        .. note:: **Extension: missing group values.** TraMineR stops when
           ``group`` holds a missing value, in every mode. traminex leaves
           such a case out of every group: its distance is 0 in the default
           mode, with ``allcenter=True`` it gets its distance to each group's
           centre (computed as for any other non-member), and it is never a
           medoid.
    medoids_index : {"first", "all", True, None}, optional
        If given, return the index (0-based) of each group's medoid instead
        of the distances: "first" returns the first index attaining the
        minimum, "all" returns every index. ``True`` means "first".

        .. note:: **Deliberate divergence.** R maps ``TRUE`` to ``"First"``
           and then rejects it (its own value check only accepts lowercase
           ``"first"``/``"all"``), so ``medoids.index=TRUE`` always errors in
           TraMineR 2.2-14. traminex honours the documented intent and
           treats ``True`` as ``"first"``.
    allcenter : bool, default False
        If True, return a DataFrame with the distance of *every* case to
        *each* group's centre (columns = group levels). Ignored when
        ``medoids_index`` is given (R behaviour).
    weights : array-like, optional
        Per-case weights of length n.
    squared : bool, default False
        Compute on squared dissimilarities.

    Returns
    -------
    numpy.ndarray, pandas.DataFrame, scalar or dict
        Distances per case (default), one column per group
        (``allcenter=True``), or the medoid index/indices — a scalar/array
        for a single group, else a dict keyed by group level.
    """
    if isinstance(medoids_index, bool):
        medoids_index = "first" if medoids_index else None
    if medoids_index is not None and medoids_index not in ("first", "all"):
        raise ValueError('medoids_index must be one of "first", "all" or None')
    ret_medoids = medoids_index is not None
    if ret_medoids:
        allcenter = False

    D = _as_square_matrix(diss, squared=squared)
    n = D.shape[0]
    if group is None:
        group = np.ones(n, dtype=int)
    declared = _declared_categories(group)
    group = np.asarray(group)
    if group.shape[0] != n:
        raise ValueError("group must have one value per observation")
    w = _resolve_weights(n, weights)

    # R: grp <- factor(group) (disscenter.R:47) keeps a factor's level order.
    labels, codes, firsts = _factor(group, declared)
    k = len(labels)

    if allcenter:
        out = np.zeros((n, k), dtype=float)
        for gi in range(k):
            members = np.flatnonzero(codes == gi)
            others = np.flatnonzero(codes != gi)
            d_mem, d_oth = _weighted_inertia_contrib_ext(D, members, others, w)
            wg = w[members]
            # R subtracts weighted.mean(dT[members], weights[members])/2
            # unconditionally (SOURCES/TraMineR/R/disscenter.R:80), so a zero or
            # non-finite group weight propagates NaN/Inf into every column entry.
            with np.errstate(divide="ignore", invalid="ignore"):
                mean_half = float(np.sum(d_mem * wg) / float(wg.sum())) / 2.0
                out[members, gi] = d_mem - mean_half
                out[others, gi] = d_oth - mean_half
        return pd.DataFrame(out, columns=labels)

    ret = np.zeros(n, dtype=float)
    medoids: Dict[Any, Any] = {}
    for gi, lv in enumerate(firsts):
        # The dict is keyed by group value; a level that merges values equal to
        # 15 significant digits is keyed by its first value.
        members = np.flatnonzero(codes == gi)
        dc = _weighted_inertia_contrib(D, members, w)
        wg = w[members]
        # Unconditional in R (SOURCES/TraMineR/R/disscenter.R:87); see
        # _disscenter_values for the same rule.
        dc = dc - _r_weighted_mean(dc, wg) / 2.0
        ret[members] = dc
        if ret_medoids:
            min_dc = float(dc.min())
            cands = members[dc == min_dc]
            medoids[lv] = int(cands.min()) if medoids_index == "first" else np.sort(cands)
    if ret_medoids:
        if k == 1:
            return medoids[firsts[0]]
        return medoids
    return ret


def gower_matrix(
    D: np.ndarray, *, squared: bool = True, weights: Optional[np.ndarray] = None
) -> np.ndarray:
    """Gower-centered Gram matrix — R-faithful.

    Port of TraMineR 2.2-14 ``gower_matrix(diss, squared=TRUE, weights=NULL)``
    (`SOURCES/TraMineR/R/dissmfacw.R:11`)::

        (I - 1 s^T) (D' / -2) (I - s 1^T),  s = w / sum(w),  D' = D^2 if squared

    With the defaults (``squared=True``, ``weights=None``) this is the classical
    double-centred Gram matrix ``-0.5 * J D^2 J`` used in classical MDS
    (``J = I - (1/n) 11^T``). Anchored to live R by
    ``tests/test_r_goldens_gower.py``.

    Parameters
    ----------
    D : np.ndarray
        Square symmetric distance matrix of shape (n, n).
    squared : bool, default True
        Square the dissimilarities before centring (R default ``squared=TRUE``).
    weights : np.ndarray, optional
        Per-case weights of length n. ``None`` means uniform weights, recovering
        the unweighted double-centring. Any total is taken, as in R: a zero total
        gives NaN and a negative one finite values, since the shares
        ``w / sum(w)`` still sum to 1 (TX-AUD-20261011-005).

    Returns
    -------
    np.ndarray
        Gram (inner product) matrix of shape (n, n). With the defaults its
        eigenvalues give the variance explained by each MDS dimension.
    """
    D = np.asarray(D, dtype=float)
    n = D.shape[0]
    if n == 0:
        return np.zeros((0, 0), dtype=float)
    if weights is None:
        w = np.ones(n, dtype=float)
    else:
        w = np.asarray(weights, dtype=float)
        if w.shape[0] != n:
            raise ValueError("weights must have length equal to nrow(D)")
    from .diss_mfa import _gower_matrix_weighted

    return _gower_matrix_weighted(D, squared, w)


def dissvar_summary(D: np.ndarray, labels: Optional[np.ndarray] = None) -> Dict[str, Any]:
    """Variance summaries for a distance matrix, optionally by group.

    .. note:: **Exploratory helper, not the TraMineR discrepancy.** This is the
       simplified, unweighted summary formerly exported as ``dissvar``. The
       R-faithful weighted discrepancy now lives at :func:`dissvar`.

    Parameters
    ----------
    D : np.ndarray
        Square distance matrix of shape (n, n)
    labels : np.ndarray, optional
        Group labels for each observation. If provided, computes
        within-group and between-group distance statistics.

    Returns
    -------
    dict
        Dictionary with keys:
        - 'overall_var': variance of all pairwise distances
        - 'within_mean': mean within-group distance (if labels provided)
        - 'between_mean': mean between-group distance (if labels provided)
    """
    D = np.asarray(D, dtype=float)
    n = D.shape[0]
    tril = D[np.tril_indices(n, k=-1)]
    out: Dict[str, Any] = {"overall_var": float(np.var(tril))}
    if labels is not None:
        labels = np.asarray(labels)
        mask_within = labels[:, None] == labels[None, :]
        within = D[np.tril_indices(n, k=-1)][(mask_within[np.tril_indices(n, k=-1)])]
        between = D[np.tril_indices(n, k=-1)][(~mask_within[np.tril_indices(n, k=-1)])]
        if within.size > 0:
            out["within_mean"] = float(within.mean())
        if between.size > 0:
            out["between_mean"] = float(between.mean())
    return out


def dissassoc_effect(D: np.ndarray, labels: np.ndarray) -> Dict[str, Any]:
    """Simple between/within effect size for a distance matrix.

    .. note:: **Exploratory helper, not the TraMineR test.** This is the
       simplified, unweighted effect size formerly exported as ``dissassoc``.
       The R-faithful weighted pseudo-F / pseudo-R² with permutation test now
       lives at :func:`dissassoc`.

    Quantifies how well group labels explain distance structure using a
    simple effect size: (mean_between - mean_within) / mean_between.

    Parameters
    ----------
    D : np.ndarray
        Square distance matrix of shape (n, n)
    labels : np.ndarray
        Group labels for each observation. Any number of unique group values is accepted.

    Returns
    -------
    dict
        Dictionary with keys:
        - 'mean_between': mean distance between groups
        - 'mean_within': mean distance within groups
        - 'effect': (mean_between - mean_within) / mean_between, where 1
          indicates zero within-group distance and negative values indicate
          larger within-group than between-group distances. Returns 0 when
          mean_between is not positive.
    """
    D = np.asarray(D, dtype=float)
    labels = np.asarray(labels)
    n = D.shape[0]
    tri_idx = np.tril_indices(n, k=-1)
    mask_within = labels[:, None] == labels[None, :]
    within = D[tri_idx][(mask_within[tri_idx])]
    between = D[tri_idx][(~mask_within[tri_idx])]
    mb = float(between.mean()) if between.size > 0 else 0.0
    mw = float(within.mean()) if within.size > 0 else 0.0
    effect = (mb - mw) / mb if mb > 0 else 0.0
    return {"mean_between": mb, "mean_within": mw, "effect": effect}


def dissassoc_test(
    D: np.ndarray, labels: np.ndarray, *, n_perm: int = 999, random_state: Optional[int] = None
) -> Dict[str, Any]:
    """Permutation test for the simple :func:`dissassoc_effect` effect size.

    .. note:: **Exploratory helper.** For the R-faithful TraMineR permutation
       test (pseudo F / pseudo R² / Bartlett / Levene), use :func:`dissassoc`.

    Returns a dict with observed effect, p_value (one-sided, effect_perm >= effect_obs), and n_perm.
    """
    rng = np.random.default_rng(random_state)
    obs = dissassoc_effect(D, labels)
    effect_obs = float(obs["effect"])
    labels = np.asarray(labels)
    count = 0
    for _ in range(int(n_perm)):
        perm = rng.permutation(labels)
        eff = dissassoc_effect(D, perm)["effect"]
        if eff >= effect_obs:
            count += 1
    p = (count + 1.0) / (n_perm + 1.0)
    return {"effect": effect_obs, "p_value": float(p), "n_perm": int(n_perm)}


def _zscore_cols(M: np.ndarray) -> np.ndarray:
    """Per-column z-score (R ``scale``: centre by mean, divide by sd, ddof=1)."""
    mu = M.mean(axis=0)
    sd = M.std(axis=0, ddof=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        return (M - mu) / sd


def _weighted_midrank(values: ArrayLike, weights: ArrayLike) -> np.ndarray:
    """TraMineR ``w.rank`` weighted midranks, including its final-tie fallback.

    A tie receives ``lower_weight + (tie_weight + 1) / 2``. TraMineR's loop
    leaves the first member of a tied maximum group unset and then replaces
    that value with its ordinary midrank; this helper preserves that observable
    behavior.
    """
    from scipy.stats import rankdata

    x = np.asarray(values, dtype=float).reshape(-1)
    w = np.asarray(weights, dtype=float).reshape(-1)
    if x.size != w.size:
        raise ValueError("weighted ranks require one weight per value")
    if not np.all(np.isfinite(x)):
        raise ValueError("weighted ranks require finite, non-constant distances")
    order = np.argsort(x, kind="stable")
    sorted_x = x[order]
    sorted_w = w[order]
    sorted_ranks = np.full(x.size, np.nan, dtype=float)
    lower_weight = 0.0
    start = 0
    while start < x.size:
        stop = start + 1
        while stop < x.size and sorted_x[stop] == sorted_x[start]:
            stop += 1
        tie_weight = float(sorted_w[start:stop].sum())
        rank = lower_weight + (tie_weight + 1.0) / 2.0
        if stop == x.size and stop - start > 1:
            sorted_ranks[start + 1 : stop] = rank
        else:
            sorted_ranks[start:stop] = rank
        lower_weight += tie_weight
        start = stop

    result = np.empty_like(sorted_ranks)
    result[order] = sorted_ranks
    missing = np.isnan(result)
    if np.any(missing):
        ordinary = rankdata(x, method="average")
        result[missing] = ordinary[missing]
    return result


def _wtd_corr_p(M: np.ndarray, ww: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Weighted Pearson correlation matrix + t-test p-values (R ``wtd.cor``).

    Mirrors TraMineR's vendored ``onecor.wtd``: each column is weight-standardised
    (``stdz``: subtract ``wtd.mean``, divide by ``sqrt(wtd.var)`` with ``ddof=1``);
    the correlation is the weighted slope ``sum(ww*Xa*Xb)/(sum(ww)-1)`` and the
    p-value is the slope's two-sided t-test from the weighted ``lm`` summary, with
    ``se = sqrt((1-r^2)/(m-2))`` and ``df = m-2`` (``m`` = number of pairs).
    """
    from scipy import stats as _spstats

    sw = float(ww.sum())
    mu = (ww[:, None] * M).sum(axis=0) / sw
    Mc = M - mu
    var = (ww[:, None] * Mc**2).sum(axis=0) / (sw - 1.0)
    Xs = Mc / np.sqrt(var)
    corr = (Xs.T @ (ww[:, None] * Xs)) / (sw - 1.0)
    corr = np.clip(corr, -1.0, 1.0)
    np.fill_diagonal(corr, 1.0)
    # R lm excludes zero-weight observations from its residual degrees of
    # freedom even though wtd.cor retains them in the input vectors.
    df = int(np.count_nonzero(ww > 0)) - 2
    with np.errstate(divide="ignore", invalid="ignore"):
        tval = corr * np.sqrt(df / (1.0 - corr**2))
    pval = 2.0 * _spstats.t.sf(np.abs(tval), df)
    pval = np.where(np.abs(corr) >= 1.0, 0.0, pval)
    return corr, pval


def dissdomassoc(
    domain_dists: "list[np.ndarray]",
    *,
    jointdiss: np.ndarray | None = None,
    what: "Sequence[str]" = ("pearson", "R2"),
    weights: np.ndarray | None = None,
    dnames: "Optional[Sequence[str]]" = None,
    w_rank: bool = False,
) -> Dict[str, Any]:
    """Measures of association between domains of sequence data -- R-faithful.

    Port of TraMineR 2.2-14 ``dissdomassoc`` (``dissdomassoc.R``): pairwise
    weighted Pearson / Spearman correlations between the (lower-triangle)
    dissimilarity vectors of each domain, their R2, and Cronbach's alpha.
    Anchored to live R by ``tests/test_r_goldens_dissdomassoc.py``.

    Parameters
    ----------
    domain_dists : list of np.ndarray
        One square (n, n) dissimilarity matrix per domain (>= 2 domains).
    jointdiss : np.ndarray, optional
        Dissimilarity of the combined (joint) domain; appended as an extra column.
    what : sequence of str, default ("pearson", "R2")
        Any of ``"pearson"``, ``"spearman"``, ``"R2"``, ``"cronbach"``,
        ``"cron.subsets"``, or ``"all"``.
    weights : np.ndarray, optional
        Per-case weights (length n). ``None`` -> uniform.
    dnames : sequence of str, optional
        Domain names (default ``dom.1 ... dom.k``).
    w_rank : bool, default False
        Use TraMineR's case-pair-weighted midranks for Spearman correlation.

    Returns
    -------
    dict
        Keys among ``Pearson``/``p.pearson``/``Pearson.Rsquare``/``Spearman``/
        ``p.spearman``/``Spearman.Rsquare``/``Cronbach``/``Cronbach.subsets``/
        ``dnames`` depending on ``what``.
    """
    from itertools import combinations

    if not isinstance(domain_dists, (list, tuple)) or len(domain_dists) < 2:
        raise ValueError("domain_dists should be a list of at least two distance matrices")

    whatlist = ["pearson", "spearman", "R2", "cronbach", "cron.subsets", "all"]
    if isinstance(what, str):
        what = (what,)
    what = list(what)
    if not all(value in whatlist for value in what):
        raise ValueError(
            "bad what values; allowed are pearson, spearman, R2, cronbach, cron.subsets, and all"
        )
    if "all" in what:
        what = ["pearson", "spearman", "R2", "cron.subsets"]
    if "R2" in what and not any(w in what for w in ("pearson", "spearman")):
        what = ["pearson"] + what

    mats = [np.asarray(D, dtype=float) for D in domain_dists]
    ndom = len(mats)
    ncases = mats[0].shape[0]
    # R's as.dist (dissdomassoc.R:55, 63) takes the lower triangle column by
    # column: for i < j in NumPy's row-major upper-triangle order that is
    # D[j, i], which also reads an asymmetric matrix as R does.
    i_idx, j_idx = np.triu_indices(ncases, k=1)
    cols = np.column_stack([D[j_idx, i_idx] for D in mats])

    w = np.ones(ncases) if weights is None else np.asarray(weights, dtype=float)
    if w.shape[0] != ncases:
        raise ValueError("length of weights not equal to number of cases!")
    ww = w[i_idx] * w[j_idx]

    if dnames is None:
        dnames = [f"dom.{i + 1}" for i in range(ndom)]
    dnames = list(dnames)

    if jointdiss is not None:
        jd = np.asarray(jointdiss, dtype=float)[j_idx, i_idx]
        cols = np.column_stack([cols, jd])
        ndomv: Any = np.array([ndom - 1] * ndom + [ndom], dtype=float)
    else:
        ndomv = float(ndom - 1)

    cols_s = _zscore_cols(cols)

    def _rsquare(corr: np.ndarray) -> np.ndarray:
        tmp = corr.copy()
        if jointdiss is not None:
            tmp[:ndom, ndom] = 0.0
        return (np.sum(tmp**2, axis=1) - 1.0) / ndomv

    res: Dict[str, Any] = {}

    if "pearson" in what:
        corr, pval = _wtd_corr_p(cols_s, ww)
        res["Pearson"] = corr
        res["p.pearson"] = pval
        if "R2" in what:
            res["Pearson.Rsquare"] = _rsquare(corr)

    if "spearman" in what:
        from scipy.stats import rankdata

        ordinary_ranks = np.column_stack([rankdata(cols_s[:, c]) for c in range(cols_s.shape[1])])
        if weights is not None and not np.all(w == 1) and w_rank:
            ranks = np.column_stack(
                [_weighted_midrank(cols_s[:, c], ww) for c in range(cols_s.shape[1])]
            )
        else:
            ranks = ordinary_ranks
        ranks_s = _zscore_cols(ranks)
        corr, pval = _wtd_corr_p(ranks_s, ww)
        res["Spearman"] = corr
        res["p.spearman"] = pval
        if "R2" in what:
            res["Spearman.Rsquare"] = _rsquare(corr)

    if any(c in what for c in ("cronbach", "cron.subsets")):
        sigmatot = float(np.var(cols_s[:, :ndom].sum(axis=1), ddof=1))
        res["Cronbach"] = (ndom / (ndom - 1)) * (1.0 - ndom / sigmatot)

    if "cron.subsets" in what and ndom > 2:
        subs: Dict[str, float] = {}
        for p in range(ndom - 1, 1, -1):
            for sub in combinations(range(ndom), p):
                sig = float(np.var(cols_s[:, list(sub)].sum(axis=1), ddof=1))
                label = "(" + ",".join(dnames[i] for i in sub) + ")"
                subs[label] = (p / (p - 1)) * (1.0 - p / sig)
        res["Cronbach.subsets"] = subs

    res["dnames"] = dnames
    return res
