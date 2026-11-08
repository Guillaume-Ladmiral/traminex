"""R-faithful sequence regression tree (TraMineR ``seqtree``).

Port of TraMineR 2.2-14 ``seqtree.R``: a thin wrapper around
:func:`traminex.disstree.disstree` that takes a state-sequence object
(:class:`traminex.sequences.STSList`), reads the per-sequence (survey)
weights from it (R: ``attr(seqdata, "weights")``, set by
``seqdef(weights=)``), and computes the dissimilarity matrix with
``seqdist`` when one is not supplied.

The scikit-learn helpers (``seqtree_fit`` & co.: per-sequence features fed to
a scikit-learn decision tree) are *exploratory*, not TraMineR's tree; they live
in :mod:`traminex.seqtree_explore` and are re-exported here so the public
top-level API keeps resolving.
"""

from __future__ import annotations

from typing import Any, Dict, Optional, Union

import numpy as np

from ._deprecated import MISSING, resolve_deprecated_alias
from .disstree import DissTree, disstree
from .disstree_viz import seqtree2dot, seqtreedisplay  # noqa: F401  (R-faithful seq-tree viz, B-12)

# Re-export the exploratory sklearn-based helpers so the historical public
# names keep resolving from this module path as well (relabelled in
# seqtree_explore.py; see DUPLICATES.md A8).
from .seqtree_explore import (  # noqa: F401
    SeqTreeResult,
    seqtree_cv_fit,
    seqtree_fit,
    seqtree_leaf_summary,
    seqtree_plot,
    seqtree_predict,
)
from .sequences import STSList


def seqtree(
    seqdata: STSList,
    predictors,
    *,
    diss: Optional[np.ndarray] = None,
    weighted: bool = True,
    seqdist_args: Optional[Dict[str, Any]] = None,
    min_size: float = 0.05,
    max_depth: int = 5,
    R: int = 1000,
    pval: float = 0.01,
    weight_permutation: str = "replicate",
    squared: bool = False,
    first: Optional[str] = None,
    random_state: Optional[Union[int, np.random.Generator]] = None,
    minSize=MISSING,
    maxdepth=MISSING,
    seqdist_arg=MISSING,
) -> DissTree:
    """Weighted regression tree of state sequences — R-faithful.

    Port of TraMineR 2.2-14 ``seqtree(formula, data, weighted=, min.size=,
    max.depth=, R=, pval=, weight.permutation=, seqdist.args=, diss=,
    squared=, first=)``; the formula interface becomes ``seqdata`` plus an
    explicit ``predictors`` DataFrame.

    With ``weighted=True`` (default) the per-sequence weights stored on
    ``seqdata`` (``seqdef(..., weights=)``) drive the weighted-inertia
    splits, exactly as R reads ``attr(seqdata, "weights")``. When ``diss``
    is None, the dissimilarity matrix is computed with
    ``traminex.seqdist(seqdata, **seqdist_args)`` (default
    ``{"method": "LCS", "norm": "auto"}``, matching R's
    ``seqdist.args=list(method="LCS", norm="auto")``). An explicitly supplied
    dictionary replaces that default rather than being merged into it; an
    empty one calls seqdist with no arguments, as R's ``list()`` does, and
    seqdist refuses it for want of substitution costs (R stops for want of a
    method).

    A node's ``info["medoid"]`` and ``info["ind"]`` are zero-based rows of
    ``seqdata`` (R's are one-based), as in :func:`disstree`.
    """
    min_size = resolve_deprecated_alias(
        min_size,
        minSize,
        canonical_name="min_size",
        alias_name="minSize",
        function_name="seqtree",
        default=0.05,
    )
    max_depth = resolve_deprecated_alias(
        max_depth,
        maxdepth,
        canonical_name="max_depth",
        alias_name="maxdepth",
        function_name="seqtree",
        default=5,
    )
    seqdist_args = resolve_deprecated_alias(
        seqdist_args,
        seqdist_arg,
        canonical_name="seqdist_args",
        alias_name="seqdist_arg",
        function_name="seqtree",
        default=None,
    )

    if not isinstance(seqdata, STSList):
        raise TypeError("seqdata must be an STSList (see seqdef)")
    if diss is None:
        from .distances import seqdist

        # R passes seqdist.args as given (seqtree.R:15-16): only None takes the default,
        # and an empty dict calls seqdist without a method, which refuses it (TX-AUD-20261010-080).
        args = dict({"method": "LCS", "norm": "auto"} if seqdist_args is None else seqdist_args)
        diss = seqdist(seqdata, **args)
    weights = seqdata.weights if weighted else None
    tree = disstree(
        diss,
        predictors,
        weights=weights,
        min_size=min_size,
        max_depth=max_depth,
        R=R,
        pval=pval,
        weight_permutation=weight_permutation,
        squared=squared,
        first=first,
        random_state=random_state,
    )
    tree.info["method"] = "seqtree"
    tree.info["object"] = seqdata
    return tree
