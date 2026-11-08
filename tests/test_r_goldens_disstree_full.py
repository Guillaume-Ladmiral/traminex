"""R-golden tests for the disstree leaf helpers (RP2 / B-5).

Anchors traminex's R-faithful ``disstreeleaf(label=True)`` /
``disstree_get_rules`` / ``disstree_assign`` to TraMineR 2.2-14 on the shared
Phase-2 fixture and the same deterministic tree as
``tests/test_r_goldens_disstree.py`` (R=1 -> no permutation gate). Goldens by
``scripts/r_goldens/make_goldens_disstree_full.R``:

- ``disstree_full_leaflabels.txt`` — one root->leaf path label per case
- ``disstree_full_rules.txt``      — the R-syntax classification rules
- ``disstree_full_assign.csv``     — disstree.assign(rules, pred): 1-based index
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from r_goldens_utils import R_GOLDEN_ATOL, load_golden, load_golden_lines

from traminex.disstree import (
    disstree,
    disstree_assign,
    disstree_get_rules,
    disstreeleaf,
)


def _inputs():
    D = load_golden("phase2_input_D")
    return D


def _predictors(n=15):
    i = np.arange(1, n + 1)
    return pd.DataFrame(
        {
            "A": ["a", "a", "b", "b", "a", "b", "b", "b", "c", "c", "c", "c", "c", "a", "c"],
            "X": ((i * 7) % 15) + 0.1 * i,
        }
    )


def _tree():
    D = _inputs()
    return disstree(D, _predictors(), R=1, min_size=2, max_depth=4)


def test_disstreeleaf_labels_match_r():
    tree = _tree()
    labels = disstreeleaf(tree, label=True, collapse="; ")
    golden = load_golden_lines("disstree_full_leaflabels")
    assert [str(s) for s in labels] == golden


def test_disstreeleaf_ids_unchanged():
    # label=False still equals the existing fitted-leaf golden (regression guard).
    tree = _tree()
    np.testing.assert_array_equal(
        disstreeleaf(tree), load_golden("disstree_unw_fitted").ravel().astype(int)
    )
    np.testing.assert_array_equal(disstreeleaf(tree), tree.fitted)


def test_disstree_get_rules_match_r():
    tree = _tree()
    rules = disstree_get_rules(tree, collapse="; ")
    golden = load_golden_lines("disstree_full_rules")
    assert rules == golden


def test_disstree_assign_matches_r():
    tree = _tree()
    rules = disstree_get_rules(tree, collapse="; ")
    assigned = disstree_assign(rules, _predictors())
    np.testing.assert_allclose(
        assigned, load_golden("disstree_full_assign").ravel(), rtol=0, atol=R_GOLDEN_ATOL
    )


def test_disstree_assign_consistent_with_leaf_partition():
    # Each case's assigned rule must reproduce its leaf membership: the rule the
    # case matches is exactly the rule built from that case's leaf label.
    tree = _tree()
    rules = disstree_get_rules(tree, collapse="; ")
    assigned = disstree_assign(rules, _predictors())
    labels = [str(s) for s in disstreeleaf(tree, label=True, collapse="; ")]
    # map each leaf label to its rule index (sorted order, as get_rules)
    sorted_labels = sorted(set(labels))
    for case_label, idx in zip(labels, assigned):
        assert int(idx) - 1 == sorted_labels.index(case_label)
