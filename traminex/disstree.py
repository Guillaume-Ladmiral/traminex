"""R-faithful dissimilarity regression tree (TraMineR ``disstree``).

Port of TraMineR 2.2-14 ``disstree.R`` (Studer, Ritschard, Gabadinho &
Müller): a regression tree on a dissimilarity matrix that recursively picks,
among the candidate predictors, the binary partition of a predictor's levels
that minimises the residual weighted inertia (Ward sense), optionally gated
by a permutation test of the split's discrepancy association.

Anchored to live R by ``tests/test_r_goldens_disstree.py`` (goldens generated
by ``scripts/r_goldens/make_goldens_disstree.R``), including the weighted
path — see DUPLICATES.md A8 / migration plan Phase 2.4.

The scikit-learn helpers (``disstree_fit`` & co.: a classical MDS embedding
fed to a scikit-learn decision tree) are *exploratory*, not TraMineR's tree; they
live in :mod:`traminex.disstree_explore` and are re-exported here so the public
top-level API keeps resolving.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from itertools import combinations
from typing import Any, Dict, List, Literal, Optional, Union, cast

import numpy as np
import pandas as pd

from ._categorical_order import r_factor, r_string_to_double
from ._deprecated import MISSING, resolve_deprecated_alias
from ._rsum import r_sum as _r_sum
from .diss import (
    _PRODUCT_BLOCK,
    DissAssocResult,
    _c_accumulate,
    _c_weighted_inertia_dist,
    _resolve_weights,
    _tramineR_permutation,
    _weighted_inertia_contrib,
    dissassoc,
    dissvar,
    weighted_inertia,
)

# Re-export the exploratory sklearn-based helpers so the historical public
# names keep resolving from this module path as well (relabelled in
# disstree_explore.py; see DUPLICATES.md A8).
from .disstree_explore import (  # noqa: F401
    DissTreeResult,
    disstree2dot_explore,
    disstree_feature_importance,
    disstree_fit,
    disstree_leaf_summary,
    disstree_plot,
    disstree_predict,
)

# R-faithful leaf helpers (B-5) operate on the fitted DissTree below; re-exported
# here so the public names resolve from this module path too.
from .disstree_helpers import (  # noqa: F401
    disstree_assign,
    disstree_get_rules,
    disstreeleaf,
)

# R-faithful GraphViz DOT export / display (T20, B-12); the public ``disstree2dot``
# is now this faithful port (the sklearn one is ``disstree2dot_explore`` above).
from .disstree_viz import (  # noqa: F401
    disstree2dot,
    disstree2dotp,
    disstreedisplay,
)

ArrayLike = Union[np.ndarray, "pd.Series", list]


@dataclass
class DissTreeSplit:
    """A binary split of a node (R ``DissTreeSplit``)."""

    varindex: int
    prob: np.ndarray
    info: Dict[str, Any]
    index: Optional[np.ndarray] = None  # per-level side (1/2), categorical splits
    breaks: Optional[float] = None  # threshold, numeric-predictor splits only
    na_group: Optional[int] = None  # side (1/2) that received missing values
    labels: Optional[List[Any]] = None  # level values, categorical splits


@dataclass
class DissTreeNode:
    """A node of the dissimilarity tree (R ``DissTreeNode``)."""

    id: int
    info: Dict[str, Any]
    split: Optional[DissTreeSplit] = None
    kids: Optional[List["DissTreeNode"]] = None

    def __repr__(self) -> str:
        """A bounded summary: the id, n and the children's ids (the fields hold index arrays)."""
        kids = None if self.kids is None else [int(kid.id) for kid in self.kids]
        return f"DissTreeNode(id={self.id}, n={_summary_number(self.info.get('n'))}, kids={kids})"


def _summary_number(value: Any) -> str:
    """A number in a summary repr, to 5 significant digits (R's print digits), or its repr."""
    if isinstance(value, (int, float, np.integer, np.floating)) and not isinstance(value, bool):
        return f"{float(value):.5g}"
    return repr(value)


@dataclass
class DissTree:
    """Fitted dissimilarity tree (R ``disstree`` object).

    Attributes
    ----------
    root : DissTreeNode
        Tree root; nodes carry ``info`` (ind, n, vardis, depth, medoid,
        splitschedule) and ``split`` details (SCres, lpop/rpop, lvar/rvar,
        R2, pval).
    fitted : numpy.ndarray
        Leaf node id of each case.
    adjustment : DissAssocResult
        Global discrepancy analysis of the fitted leaf grouping.
    """

    root: DissTreeNode
    fitted: np.ndarray
    info: Dict[str, Any]
    data: pd.DataFrame
    weights: np.ndarray
    adjustment: DissAssocResult = field(repr=False, default=None)  # type: ignore[assignment]

    def __repr__(self) -> str:
        """A bounded summary, as R's print.disstree summarises a tree: node, leaf and depth counts, n and the global R2
        (``stat[3, 1]``, print.disstree.R). The fields are unchanged."""
        nodes = leaves = depth = 0
        stack = [self.root]
        while stack:
            node = stack.pop()
            nodes += 1
            depth = max(depth, int(node.info.get("depth", 0)))
            if node.kids is None:
                leaves += 1
            else:
                stack.extend(node.kids)
        r2 = None if self.adjustment is None else self.adjustment.stat.iloc[2, 0]
        return (
            f"DissTree(n={_summary_number(self.root.info.get('n'))}, nodes={nodes}, leaves={leaves}, "
            f"depth={depth}, R2={_summary_number(r2)})"
        )


@dataclass(frozen=True)
class _PredictorKind:
    """How R's ``DTNGroupFactorBinary`` sees one predictor column.

    R builds ``factor(pred, ordered = is.ordered(pred) || is.numeric(pred))``
    (``disstree.R:404``): a factor keeps its declared level order, an ordered
    factor or a numeric predictor is split over contiguous level prefixes, and
    only a numeric predictor reports a threshold (``breaks``, ``disstree.R:543``);
    every factor, ordered or not, reports level labels. A logical predictor is
    not numeric in R, so it is an unordered factor with levels FALSE/TRUE.
    """

    ordered: bool
    numeric: bool
    declared_levels: Optional[List[Any]]


def _predictor_kind(col: pd.Series) -> _PredictorKind:
    if isinstance(col.dtype, pd.CategoricalDtype):
        return _PredictorKind(bool(col.cat.ordered), False, list(col.cat.categories))
    if pd.api.types.is_bool_dtype(col.dtype):
        return _PredictorKind(False, False, None)
    numeric = bool(pd.api.types.is_numeric_dtype(col.dtype))
    return _PredictorKind(numeric, numeric, None)


def _node_factor(
    values: np.ndarray, declared_levels: Optional[List[Any]]
) -> tuple[List[str], np.ndarray]:
    """R's ``factor(pred)`` of the node-subset values (unused levels dropped): labels and codes.

    A categorical predictor keeps its declared category order; any other
    predictor is sorted as ``factor()`` sorts a plain vector: numbers by value,
    text by R's collation (TX-AUD-20261010-039). Levels are the values as R's
    ``as.character`` writes them, so numbers equal to 15 significant digits are
    one level (``r_factor``, TX-AUD-20261011-055). Codes are -1 for missing values.
    """
    labels, codes, _ = r_factor(values, declared_levels)
    return labels, codes


# The helpers below reproduce the rounding of R's DTNGroupFactorBinary
# (disstree.R:403-523). A level subset and its complement have mathematically
# equal residual inertia, so the last bit decides which becomes the left child
# (audit S4-04); NumPy's pairwise sums and matrix products round differently.


def _c_weighted_inter_inertia(
    D: np.ndarray, first: np.ndarray, second: np.ndarray, w: np.ndarray
) -> float:
    """``tmrWeightedInterInertia`` in its C accumulation order.

    ``weightedinertia.cpp:173-181``: ``result += w[i] * w[j] * D[i, j]`` for
    each ``i`` of the first group (outer) and each ``j`` of the second (inner).
    """
    if first.size == 0 or second.size == 0:
        return 0.0
    w_first, w_second = w[first], w[second]
    result = 0.0
    step = max(1, _PRODUCT_BLOCK // second.size)
    for start in range(0, first.size, step):
        rows = slice(start, start + step)
        block = (w_first[rows, None] * w_second[None, :]) * D[np.ix_(first[rows], second)]
        result = _c_accumulate(result, block.ravel())
    return result


def _group_factor_binary(
    D: np.ndarray,
    current_sc_res: float,
    values: np.ndarray,
    kind: _PredictorKind,
    min_size: float,
    varindex: int,
    ind: np.ndarray,
    w: np.ndarray,
) -> Optional[Dict[str, Any]]:
    """Best binary partition of one predictor (port of R ``DTNGroupFactorBinary``).

    Sizes, inertias and candidate residuals are computed with R's arithmetic
    (:func:`_r_sum`, :func:`_c_weighted_inertia_dist`,
    :func:`_c_weighted_inter_inertia`) so that ties between a level subset and
    its complement resolve as in R.
    """
    totpop = _r_sum(w[ind])
    levels, codes = _node_factor(values, kind.declared_levels)
    if len(levels) <= 1:
        return None
    is_ordered = kind.ordered
    nb_real = len(levels)
    na_mask = pd.isna(values)
    has_na = bool(na_mask.any())
    nb_grp = nb_real + 1 if has_na else nb_real

    grp_members: List[np.ndarray] = []
    level_masks: List[np.ndarray] = []
    grp_size = np.zeros(nb_grp, dtype=float)
    for gi in range(nb_real):
        cond = codes == gi
        level_masks.append(cond)
        grp_members.append(ind[cond])
        grp_size[gi] = _r_sum(w[ind[cond]])
    if has_na:
        grp_members.append(ind[na_mask])
        grp_size[nb_grp - 1] = _r_sum(w[ind[na_mask]])

    # Triangular inertia matrix: diagonal = within-group weighted pairwise
    # sums (R: the C inertia times the group size), below-diagonal [j, i] =
    # between-group weighted inter-inertia. Non-finite case weights propagate
    # NaN/Inf as in R, so that arithmetic must not raise numpy warnings.
    inertia_mat = np.zeros((nb_grp, nb_grp), dtype=float)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        for gi in range(nb_grp):
            inertia_mat[gi, gi] = _c_weighted_inertia_dist(D, grp_members[gi], w) * grp_size[gi]
            for gj in range(gi + 1, nb_grp):
                inertia_mat[gj, gi] = _c_weighted_inter_inertia(
                    D, grp_members[gi], grp_members[gj], w
                )
        sc_res_all = _r_sum(np.diag(inertia_mat) / grp_size)
    if np.isnan(sc_res_all) or np.isnan(current_sc_res):
        # R compares with `if (SCres > currentSCres)` (disstree.R:476), which
        # stops ("missing value where TRUE/FALSE needed") when either side is
        # NA/NaN. Skipping the predictor instead returned a tree with the root
        # alone as if nothing could be explained.
        raise ValueError(
            "disstree cannot compare candidate splits: a node's residual discrepancy "
            "is NaN, because the dissimilarity matrix has a missing value off its "
            "diagonal, a weight is missing, or all cases of a predictor level have "
            "weight 0 (R's disstree stops here too)"
        )
    if sc_res_all > current_sc_res:
        return None

    def inertia_of(co: List[int], pop: float) -> float:
        # R: sum(inertiaMat[co, co]) / pop -- column-major, long double.
        return _r_sum(inertia_mat[np.ix_(co, co)].ravel(order="F")) / pop

    best_sc_res = current_sc_res
    best: Optional[Dict[str, Any]] = None
    all_groups = list(range(nb_grp))
    max_grp = (nb_grp - 1) if is_ordered else int(np.ceil(nb_grp / 2.0))
    for p in range(1, max_grp + 1):
        if is_ordered:
            combos: List[tuple] = [tuple(range(p))]
            if has_na:
                combos.append(tuple(list(range(p)) + [nb_grp - 1]))
        else:
            combos = list(combinations(all_groups, p))
        for co in combos:
            popc = _r_sum(grp_size[list(co)])
            popothc = totpop - popc
            if popc > min_size and popothc > min_size:
                othc = [g for g in all_groups if g not in co]
                ico = inertia_of(list(co), popc)
                iothc = inertia_of(othc, popothc)
                sc_res = ico + iothc
                if sc_res < best_sc_res:
                    best_sc_res = sc_res
                    best = {
                        "co": list(co),
                        "othc": othc,
                        "ico": ico,
                        "iothc": iothc,
                        "popc": popc,
                        "popothc": popothc,
                    }
    if best is None:
        return None

    prob = np.array([best["popc"], best["popothc"]]) / totpop
    info = {
        "lpop": best["popc"],
        "rpop": best["popothc"],
        "lvar": best["ico"] / best["popc"],
        "rvar": best["iothc"] / best["popothc"],
        "SCres": best_sc_res,
    }
    co_real = [g for g in best["co"] if g < nb_real]
    variable = np.zeros(values.shape[0], dtype=bool)
    for g in co_real:
        variable |= level_masks[g]
    if kind.numeric:
        # R's breaks <- as.numeric(llgrp[...]) (disstree.R:544): the level text read back.
        breaks = r_string_to_double(levels[max(co_real)])
        spl = DissTreeSplit(
            varindex=varindex, prob=prob, info=info, index=np.array([1, 2]), breaks=breaks
        )
    else:
        index = np.array([1 if g in best["co"] else 2 for g in range(nb_real)], dtype=int)
        spl = DissTreeSplit(
            varindex=varindex, prob=prob, info=info, index=index, labels=list(levels)
        )
    if has_na:
        na_left = (nb_grp - 1) in best["co"]
        variable[na_mask] = na_left
        spl.na_group = 1 if na_left else 2
    return {"variable": variable, "spl": spl}


def _split_permutation_pval(
    D: np.ndarray,
    grp_mask: np.ndarray,
    indiv: np.ndarray,
    w: np.ndarray,
    R: int,
    weight_permutation: str,
    rng: np.random.Generator,
) -> float:
    """Permutation p-value of a binary split (port of R ``DTNdissassocweighted``)."""
    n_full = D.shape[0]

    if weight_permutation in ("diss", "group"):
        permut_group = weight_permutation == "group"

        def statistic(ind: np.ndarray) -> np.ndarray:
            ww = w
            if permut_group:
                ww = w.copy()
                ww[indiv] = w[indiv[ind]]
            g1 = indiv[ind[grp_mask]]
            g2 = indiv[ind[~grp_mask]]
            return np.array([-(weighted_inertia(D, g1, ww) + weighted_inertia(D, g2, ww))])

        perms = _tramineR_permutation(indiv.size, R, statistic, rng)
    elif weight_permutation == "none":

        def statistic(ind: np.ndarray) -> np.ndarray:
            g1 = indiv[ind[grp_mask]]
            g2 = indiv[ind[~grp_mask]]
            return np.array([-(weighted_inertia(D, g1) + weighted_inertia(D, g2))])

        perms = _tramineR_permutation(indiv.size, R, statistic, rng)
    else:  # replicate / rounded-replicate
        reps = w[indiv].astype(int)
        grp_rep = np.repeat(grp_mask, reps)
        indiv_rep = np.repeat(indiv, reps)

        def statistic(ind: np.ndarray) -> np.ndarray:
            g1n = indiv_rep[ind[grp_rep]]
            g2n = indiv_rep[ind[~grp_rep]]
            w1 = np.bincount(g1n, minlength=n_full).astype(float)
            w2 = np.bincount(g2n, minlength=n_full).astype(float)
            r1 = weighted_inertia(D, np.flatnonzero(w1 > 0), w1)
            r2 = weighted_inertia(D, np.flatnonzero(w2 > 0), w2)
            return np.array([-(r1 + r2)])

        perms = _tramineR_permutation(indiv_rep.size, R, statistic, rng)
    return float(perms["pval"][0])


def _add_split_schedule(tree: DissTree) -> None:
    """Port of R ``DTNaddCovariateSplitschedule`` (R2-based split ordering)."""
    sc_expl: Dict[int, float] = {}

    def collect(node: DissTreeNode, parent: Optional[DissTreeNode]) -> None:
        if parent is None:
            sc_expl[node.id] = node.info["vardis"] * node.info["n"]
        else:
            sc_tot = parent.info["vardis"] * parent.info["n"]
            assert parent.split is not None
            sc_expl[node.id] = parent.split.info["R2"] * sc_tot
        if node.kids is not None:
            collect(node.kids[0], node)
            collect(node.kids[1], node)

    collect(tree.root, None)
    ids = list(sc_expl.keys())
    vals = np.array([sc_expl[i] for i in ids], dtype=float)
    # R: as.integer(factor(1 - rank(SCexpl), ordered=TRUE)) — dense rank of
    # descending SCexpl (rank() uses average ranks for ties).
    order = vals.argsort(kind="mergesort")
    ranks = np.empty(vals.size, dtype=float)
    ranks[order] = np.arange(1, vals.size + 1, dtype=float)
    # rank(na.last = TRUE) gives each NaN its own trailing rank, which argsort
    # already did; only non-NaN ties are averaged (NaN == NaN selects nothing
    # and warned, audit S4-18).
    for v in np.unique(vals[~np.isnan(vals)]):
        m = vals == v
        ranks[m] = ranks[m].mean()
    neg = 1.0 - ranks
    codes = np.searchsorted(np.unique(neg), neg) + 1
    schedule = {nid: int(c) for nid, c in zip(ids, codes)}

    def assign(node: DissTreeNode) -> None:
        node.info["splitschedule"] = schedule[node.id]
        if node.kids is not None:
            assign(node.kids[0])
            assign(node.kids[1])

    assign(tree.root)


def disstree(
    diss: ArrayLike,
    predictors: pd.DataFrame,
    weights: Optional[ArrayLike] = None,
    min_size: float = 0.05,
    max_depth: int = 5,
    R: int = 1000,
    pval: float = 0.01,
    weight_permutation: str = "replicate",
    squared: bool = False,
    first: Optional[str] = None,
    random_state: Optional[Union[int, np.random.Generator]] = None,
    *,
    minSize=MISSING,
    maxdepth=MISSING,
) -> DissTree:
    """Weighted regression tree on a dissimilarity matrix — R-faithful.

    Port of TraMineR 2.2-14 ``disstree(formula, data, weights=, min.size=,
    max.depth=, R=, pval=, weight.permutation=, squared=, first=)``; the
    formula interface becomes an explicit ``predictors`` DataFrame (one
    column per covariate). As in R, a numeric column splits over thresholds
    (rules ``x <= t``); an ordered pandas categorical splits between two
    consecutive categories in their declared order, and any other column over
    level subsets (both give rules ``x %in% c(...)``). A categorical keeps its declared
    category order, other columns are sorted, and a boolean column is an
    unordered factor with levels ``FALSE``/``TRUE``, as a logical is in R.

    Splits minimise the residual weighted inertia; with ``R > 1`` each split
    must pass a permutation test of its discrepancy association at ``pval``
    (stochastic; pass ``random_state`` for reproducibility — Python
    addition). ``min_size`` is a fraction of the total weight when ``< 1``,
    else an absolute weighted size.

    Returns a :class:`DissTree`; a node's ``info["medoid"]`` and ``info["ind"]`` are
    zero-based rows of ``diss`` (R's are one-based). Node ids, per-node weighted sizes,
    discrepancies, split statistics (SCres, lpop/rpop, lvar/rvar, R2),
    medoids, the leaf assignment and the split schedule are anchored to live
    R by ``tests/test_r_goldens_disstree.py``. A level subset and its
    complement often have mathematically equal residual inertia; which one
    becomes the left child is then decided by rounding, so candidate splits
    are evaluated with R's summation order and precision (R's ``sum()``
    accumulates in long double) and resolve as in R where NumPy's
    ``longdouble`` matches R's, as on Linux x86-64
    (``tests/test_disstree_split_arithmetic.py``).

    As in R, a missing dissimilarity off the diagonal, a missing weight, or a
    predictor level whose cases all weigh 0 makes a node's residual
    discrepancy NaN, and then splits cannot be compared: ``ValueError`` is
    raised where R stops. If no comparison is made (``max_depth=1``, or no
    predictor with two levels) the root is returned alone. The diagonal of
    ``diss`` is ignored, as R's routines never read it.
    """
    min_size = resolve_deprecated_alias(
        min_size,
        minSize,
        canonical_name="min_size",
        alias_name="minSize",
        function_name="disstree",
        default=0.05,
    )
    max_depth = resolve_deprecated_alias(
        max_depth,
        maxdepth,
        canonical_name="max_depth",
        alias_name="maxdepth",
        function_name="disstree",
        default=5,
    )

    D = np.array(diss, dtype=float)  # a copy: the diagonal is reset below
    if D.ndim != 2 or D.shape[0] != D.shape[1]:
        raise ValueError("diss must be a square dissimilarity matrix")
    if squared:
        D = D**2
    # R's C kernels never read the diagonal (their loops run over j > i), so a
    # missing or non-zero self-dissimilarity does not affect R's tree; zeroing
    # it keeps traminex's matrix products from spreading such a value.
    np.fill_diagonal(D, 0.0)
    nobs = D.shape[0]
    predictors = pd.DataFrame(predictors).reset_index(drop=True)
    if predictors.shape[0] != nobs:
        raise ValueError("dissimilarity matrix and predictors must have the same size")

    if weights is None:
        w = np.ones(nobs, dtype=float)
        weight_permutation = "none"
    else:
        w = _resolve_weights(nobs, weights)
        if weight_permutation in ("replicate", "rounded-replicate"):
            round_error = float(np.abs(np.round(w) - w).sum())
            if round_error > 0:
                if weight_permutation == "replicate":
                    raise ValueError("to permute replicates, weights must be integers")
                w = np.round(w)
    if weight_permutation in ("permute", "random-sampling"):
        # Deliberate divergence: R's disstree runs these, but its split test
        # sends them to the replicate branch, which replicates cases by
        # as.integer(weights) -- real weights truncated, 0.5 counting as 0
        # (disstree.dissassoc.R:52-64). That is not a test of the weighted split.
        raise ValueError(
            f"weight_permutation={weight_permutation!r} is not supported by disstree: "
            "R's disstree accepts it, but its split test then replicates cases by "
            "their weights truncated to integers; use 'diss', 'group', 'replicate' "
            "or 'rounded-replicate'"
        )
    if weight_permutation not in ("none", "diss", "replicate", "rounded-replicate", "group"):
        # R's disstree stops on such a value too (at its final dissassoc).
        raise ValueError(
            f"unknown weight_permutation={weight_permutation!r}; use 'none', 'diss', "
            "'replicate', 'rounded-replicate' or 'group'"
        )

    pop = float(w.sum())
    if min_size < 1:
        min_size = pop * min_size
    if R <= 1:
        pval = 1.0
    elif pval < 1.0 / R:
        warnings.warn(
            f"minimum possible p-value with R={R} permutations is {1.0 / R}; "
            f"pval changed to {1.0 / R}",
            stacklevel=2,
        )
        pval = 1.0 / R

    first_index: Optional[int] = None
    if first is not None:
        if first not in predictors.columns:
            raise ValueError(f"unknown 'first' variable {first!r}; must be a predictor column")
        first_index = int(predictors.columns.get_loc(first))

    rng = (
        random_state
        if isinstance(random_state, np.random.Generator)
        else np.random.default_rng(random_state)
    )
    counter = {"id": 1}
    col_values = [predictors[c].to_numpy() for c in predictors.columns]
    col_kinds = [_predictor_kind(predictors[c]) for c in predictors.columns]

    def init_node(ind: np.ndarray, vardis: float, depth: int) -> DissTreeNode:
        dc = _weighted_inertia_contrib(D, ind, w)
        medoid = int(ind[int(np.argmin(dc))])
        node = DissTreeNode(
            id=counter["id"],
            info={
                "depth": depth,
                "vardis": float(vardis),
                "n": float(w[ind].sum()),
                "medoid": medoid,
                "ind": ind,
            },
        )
        counter["id"] += 1
        return node

    def build(ind: np.ndarray, vardis: float, depth: int, use_first: bool) -> DissTreeNode:
        node = init_node(ind, vardis, depth)
        sc_tot = vardis * node.info["n"]
        sc_res = sc_tot
        best: Optional[Dict[str, Any]] = None
        if depth >= max_depth:
            return node
        if use_first and first_index is not None:
            best = _group_factor_binary(
                D,
                sc_res,
                col_values[first_index][ind],
                col_kinds[first_index],
                min_size,
                first_index,
                ind,
                w,
            )
            if best is not None:
                sc_res = best["spl"].info["SCres"]
        else:
            for p in range(len(col_values)):
                spl = _group_factor_binary(
                    D, sc_res, col_values[p][ind], col_kinds[p], min_size, p, ind, w
                )
                if spl is not None and (
                    best is None or spl["spl"].info["SCres"] < best["spl"].info["SCres"]
                ):
                    best = spl
                    sc_res = spl["spl"].info["SCres"]
        if best is None:
            return node
        if R > 1:
            spval = _split_permutation_pval(D, best["variable"], ind, w, R, weight_permutation, rng)
            if spval > pval:
                return node
            best["spl"].info["pval"] = spval
        node.split = best["spl"]
        node.split.info["R2"] = 1.0 - (sc_res / sc_tot)
        mask = best["variable"]
        node.kids = [
            build(ind[mask], best["spl"].info["lvar"], depth + 1, False),
            build(ind[~mask], best["spl"].info["rvar"], depth + 1, False),
        ]
        return node

    vardis = dissvar(D, weights=w)
    root = build(np.arange(nobs, dtype=int), vardis, 1, use_first=True)
    fitted = disstreeleaf(root)
    adjustment_weight_permutation = cast(
        Literal[
            "none",
            "diss",
            "replicate",
            "rounded-replicate",
            "group",
            "permute",
            "random-sampling",
        ],
        "replicate" if weight_permutation == "none" else weight_permutation,
    )
    adjustment = dissassoc(
        D,
        fitted,
        weights=None if weight_permutation == "none" else w,
        R=R,
        weight_permutation=adjustment_weight_permutation,
        random_state=rng,
    )
    tree = DissTree(
        root=root,
        fitted=fitted,
        info={
            "method": "disstree",
            "n": pop,
            "parameters": {"min_size": min_size, "max_depth": max_depth, "R": R, "pval": pval},
            "object": None,
            "weight_permutation": weight_permutation,
        },
        data=predictors,
        weights=w,
        adjustment=adjustment,
    )
    _add_split_schedule(tree)
    return tree
