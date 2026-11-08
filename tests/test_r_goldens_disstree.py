"""R-golden tests for the weighted regression tree ``disstree``/``seqtree``
(Phase 2.4 / A8).

Anchors traminex's R-faithful tree to TraMineR 2.2-14 ``disstree`` /
``seqtree`` on the shared Phase 2 fixture
(``scripts/r_goldens/make_goldens_disstree.R``). R=1 disables the stochastic
permutation gate, so the grown tree is fully deterministic: node ids (preorder
creation), weighted sizes, discrepancies, medoids, split statistics
(SCres/lpop/rpop/lvar/rvar/R2/breaks), the split schedule, leaf memberships
and the adjustment t0 statistics are pinned. The permutation gate itself is
stochastic in R too and is checked structurally.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDEN_ATOL, load_golden

from traminex.disstree import disstree
from traminex.format import seqdef
from traminex.seqtree import seqtree


def _inputs():
    D = load_golden("phase2_input_D")
    wreal = load_golden("phase2_input_weights_real").ravel()
    wint = load_golden("phase2_input_weights_int").ravel()
    return D, wreal, wint


def _predictors(n=15):
    i = np.arange(1, n + 1)
    return pd.DataFrame(
        {
            "A": ["a", "a", "b", "b", "a", "b", "b", "b", "c", "c", "c", "c", "c", "a", "c"],
            "X": ((i * 7) % 15) + 0.1 * i,
        }
    )


def _seqtree_default_inputs():
    rows = [
        ["C", "A", "A"],
        ["B", "C", "B", "A", "C"],
        ["C"],
        ["B", "C", "C", "B", "B", "C"],
        ["B", "B", "A", "A", "B", "A"],
        ["A", "C", "C", "C", "C"],
        ["B", "A"],
        ["B", "B", "A", "B", "B", "C"],
        ["B", "B", "C", "C", "C"],
        ["B", "C", "B", "B", "C", "C"],
    ]
    predictors = pd.DataFrame(
        {
            "a": pd.Categorical(["x", "z", "y", "y", "z", "y", "y", "z", "x", "y"]),
            "b": [
                -0.23405312309912332,
                1.2180151083623492,
                -0.972122955655148,
                1.8720103381243531,
                1.0205646182785473,
                -2.6004555694046183,
                -1.1264420047778816,
                -0.3695439120196192,
                -0.3646856738463821,
                1.0806622434928328,
            ],
            "c": pd.Categorical(["v", "u", "u", "u", "v", "u", "v", "v", "v", "u"]),
        }
    )
    return seqdef(rows, alphabet=["A", "B", "C"]), predictors


def _fit_seqtree_internal(x, predictors, **kwargs):
    return seqtree(x, predictors, R=1, min_size=2, max_depth=3, **kwargs)


def _flatten(tree):
    """Mirror of the R flatten_tree helper (preorder; 1-based medoid/varindex)."""
    rows = []

    def rec(node):
        info = node.info
        sp = node.split
        if sp is not None:
            tail = [
                sp.varindex + 1,
                sp.info["SCres"],
                sp.info["lpop"],
                sp.info["rpop"],
                sp.info["lvar"],
                sp.info["rvar"],
                sp.info["R2"],
                np.nan if sp.breaks is None else float(sp.breaks),
            ]
            has_split = 1.0
        else:
            tail = [np.nan] * 8
            has_split = 0.0
        rows.append(
            [
                node.id,
                info["depth"],
                info["splitschedule"],
                info["n"],
                info["vardis"],
                info["medoid"] + 1,
                has_split,
            ]
            + tail
        )
        if node.kids is not None:
            rec(node.kids[0])
            rec(node.kids[1])

    rec(tree.root)
    return np.array(rows, dtype=float)


def test_disstree_unweighted_matches_r():
    D, _, _ = _inputs()
    tree = disstree(D, _predictors(), R=1, min_size=2, max_depth=4)
    np.testing.assert_allclose(
        _flatten(tree), load_golden("disstree_unw_nodes"), rtol=0, atol=R_GOLDEN_ATOL
    )
    np.testing.assert_array_equal(
        tree.fitted, load_golden("disstree_unw_fitted").ravel().astype(int)
    )


def test_disstree_weighted_matches_r():
    D, wreal, _ = _inputs()
    tree = disstree(
        D, _predictors(), weights=wreal, R=1, min_size=2, max_depth=4, weight_permutation="diss"
    )
    np.testing.assert_allclose(
        _flatten(tree), load_golden("disstree_wreal_nodes"), rtol=0, atol=R_GOLDEN_ATOL
    )
    np.testing.assert_array_equal(
        tree.fitted, load_golden("disstree_wreal_fitted").ravel().astype(int)
    )
    np.testing.assert_allclose(
        tree.adjustment.stat["t0"].to_numpy(),
        load_golden("disstree_wreal_adjust_t0").ravel(),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


def test_seqtree_weighted_from_seqdef_matches_r():
    # Weights flow from seqdef(weights=) through the sequence object, exactly
    # like R's attr(seqdata, "weights") — the weighted path of A8.
    D, _, wint = _inputs()
    n = 15
    states = np.empty((n, 8), dtype=object)
    for r in range(1, n + 1):
        for c in range(1, 9):
            states[r - 1, c - 1] = ["A", "B", "C"][(r + c) % 3]
    x = seqdef(states, alphabet=["A", "B", "C"], weights=wint)
    tree = seqtree(
        x,
        _predictors(),
        diss=D,
        R=1,
        min_size=2,
        max_depth=4,
        weight_permutation="replicate",
    )
    np.testing.assert_allclose(
        _flatten(tree), load_golden("seqtree_wint_nodes"), rtol=0, atol=R_GOLDEN_ATOL
    )
    np.testing.assert_array_equal(
        tree.fitted, load_golden("seqtree_wint_fitted").ravel().astype(int)
    )
    assert tree.info["method"] == "seqtree"


def test_seqtree_omitted_distance_args_matches_r_default_normalization():
    x, predictors = _seqtree_default_inputs()

    actual = _fit_seqtree_internal(x, predictors)
    explicit_auto = _fit_seqtree_internal(
        x,
        predictors,
        seqdist_args={"method": "LCS", "norm": "auto"},
    )

    assert actual.root.info["vardis"] == pytest.approx(
        0.24051168840622211, rel=0, abs=R_GOLDEN_ATOL
    )
    np.testing.assert_array_equal(actual.fitted, [3, 2, 2, 2, 3, 2, 3, 3, 3, 2])
    assert actual.root.info["vardis"] == pytest.approx(
        explicit_auto.root.info["vardis"], rel=0, abs=R_GOLDEN_ATOL
    )
    np.testing.assert_array_equal(actual.fitted, explicit_auto.fitted)


def test_seqtree_explicit_distance_args_replace_omitted_defaults():
    x, predictors = _seqtree_default_inputs()

    explicit_raw = _fit_seqtree_internal(x, predictors, seqdist_args={"method": "LCS"})
    with pytest.deprecated_call(match="seqdist_arg"):
        alias_raw = _fit_seqtree_internal(x, predictors, seqdist_arg={"method": "LCS"})

    assert explicit_raw.root.info["vardis"] == pytest.approx(2.15, rel=0, abs=R_GOLDEN_ATOL)
    np.testing.assert_array_equal(explicit_raw.fitted, [5, 2, 4, 5, 2, 4, 4, 2, 5, 5])
    assert alias_raw.root.info["vardis"] == pytest.approx(
        explicit_raw.root.info["vardis"], rel=0, abs=R_GOLDEN_ATOL
    )
    np.testing.assert_array_equal(alias_raw.fitted, explicit_raw.fitted)


def test_disstree_permutation_gate_structure():
    # With R>1 the permutation gate is stochastic (in R too); check it is
    # seeded-reproducible, that strong structure passes the gate, and that a
    # pure-noise predictor is rejected (tree stays a stump).
    D, _, _ = _inputs()
    pred = _predictors()
    t1 = disstree(D, pred, R=100, pval=0.05, min_size=2, max_depth=3, random_state=11)
    t2 = disstree(D, pred, R=100, pval=0.05, min_size=2, max_depth=3, random_state=11)
    np.testing.assert_array_equal(t1.fitted, t2.fitted)
    assert t1.root.split is not None  # strong signal splits
    assert t1.root.split.info["pval"] <= 0.05

    rng = np.random.default_rng(5)
    noise = pd.DataFrame({"N": rng.permutation(["u", "v", "w"] * 5)})
    t3 = disstree(D, noise, R=200, pval=0.01, min_size=2, max_depth=3, random_state=12)
    assert t3.root.split is None  # noise rejected by the gate
