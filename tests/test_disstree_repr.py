"""TX-AUD-20261019-065: the repr of a dissimilarity tree is a bounded summary.

R's print.disstree prints a summary (print.disstree.R: the global R2, then one line per node); the ledger maps it to
traminex's fields and repr. A tree's repr printed every node, its index arrays and the predictor frame (21,275
characters for the mvad tree with two predictors). It now names the node, leaf and depth counts, n and the global R2
(R's ``stat[3, 1]``); a node's repr names its id, n and children. The fields are unchanged.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd

import traminex as tx
from traminex.disstree import DissTree, DissTreeNode

STAT_ROWS = ["Pseudo F", "Pseudo Fbf", "Pseudo R2", "Bartlett", "Levene"]


def _node(node_id: int, depth: int, n: float, kids=None) -> DissTreeNode:
    info = {"depth": depth, "vardis": 1.0, "n": n, "medoid": 0, "ind": np.arange(int(n))}
    return DissTreeNode(id=node_id, info=info, kids=kids)


def _hand_tree() -> DissTree:
    # Five nodes, three leaves, depth 3: 1 -> (2 -> (4, 5), 3).
    inner = _node(2, 2, 5000, kids=[_node(4, 3, 3000), _node(5, 3, 2000)])
    root = _node(1, 1, 10000, kids=[inner, _node(3, 2, 5000)])
    stat = pd.DataFrame({"t0": [9.0, 9.5, 0.25, 1.0, 2.0], "p_value": np.nan}, index=STAT_ROWS)
    return DissTree(
        root=root,
        fitted=np.repeat([4, 5, 3], [3000, 2000, 5000]),
        info={"n": 10000.0},
        data=pd.DataFrame({"x": np.arange(10000)}),
        weights=np.ones(10000),
        adjustment=SimpleNamespace(stat=stat),  # type: ignore[arg-type]
    )


def test_a_tree_repr_names_its_counts_n_and_global_r2() -> None:
    tree = _hand_tree()
    assert repr(tree) == "DissTree(n=10000, nodes=5, leaves=3, depth=3, R2=0.25)"
    assert str(tree) == repr(tree)


def test_a_node_repr_names_its_id_n_and_children() -> None:
    tree = _hand_tree()
    assert repr(tree.root) == "DissTreeNode(id=1, n=10000, kids=[2, 3])"
    assert repr(tree.root.kids[1]) == "DissTreeNode(id=3, n=5000, kids=None)"


def test_a_fitted_tree_repr_is_bounded() -> None:
    rng = np.random.default_rng(7)
    groups = rng.integers(0, 4, size=400)
    y = groups * 5.0 + rng.normal(scale=0.5, size=400)
    diss = np.abs(y[:, None] - y[None, :])
    predictors = pd.DataFrame({"g": groups.astype(str), "z": rng.normal(size=400)})
    tree = tx.disstree(diss, predictors, R=1, min_size=0.05)
    text = repr(tree)
    assert len(text) < 300
    r2 = float(tree.adjustment.stat.iloc[2, 0])
    assert text.startswith("DissTree(n=400, nodes=") and text.endswith(f"R2={r2:.5g})")
