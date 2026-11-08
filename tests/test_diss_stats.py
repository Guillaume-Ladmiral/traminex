import numpy as np
import pytest

from traminex.diss import (
    dissassoc_effect,
    dissassoc_test,
    dissvar,
    weighted_inertia,
    weighted_pairwise_sum,
)


def test_dissassoc_effect_and_perm():
    # dissassoc_effect is the relabelled exploratory helper formerly exported
    # as dissassoc (the R-faithful dissassoc is golden-tested separately).
    D = np.array(
        [
            [0.0, 1.0, 2.0, 2.0],
            [1.0, 0.0, 2.0, 2.0],
            [2.0, 2.0, 0.0, 1.0],
            [2.0, 2.0, 1.0, 0.0],
        ]
    )
    labels = np.array([0, 0, 1, 1])
    res = dissassoc_effect(D, labels)
    assert res["effect"] >= 0.0
    test = dissassoc_test(D, labels, n_perm=199, random_state=0)
    assert 0.0 <= test["p_value"] <= 1.0


@pytest.mark.parametrize("distance", [2.0, 1e308, np.nextafter(0.0, 1.0)])
def test_weighted_pairwise_sum_counts_each_pair_once_without_overflow(distance):
    matrix = np.array([[np.nan, distance], [distance, np.inf]])
    matrix.setflags(write=False)
    weights = np.ones(2)
    weights.setflags(write=False)
    # One off-diagonal pair, weight 1. No averaging of two copies of it.
    assert weighted_pairwise_sum(matrix, weights=weights) == distance
    assert weighted_inertia(matrix, weights=weights) == distance / 2
    assert dissvar(matrix, weights=weights) == distance / 4
    np.testing.assert_array_equal(matrix, [[np.nan, distance], [distance, np.inf]])
    np.testing.assert_array_equal(weights, [1, 1])


def test_weighted_pairwise_sum_and_r_inertia_keep_their_triangle_order():
    matrix = np.array([[np.nan, 2.0], [100.0, np.inf]])
    weights = np.array([2.0, 3.0])
    assert weighted_pairwise_sum(matrix, weights=weights) == 12
    # R's column-major C kernel reads the lower triangle for inertia.
    assert weighted_inertia(matrix, weights=weights) == 600 / 5
    assert dissvar(matrix, weights=weights) == (600 / 5) / 5
    # A reordered subset changes which input cell is above the diagonal.
    assert weighted_pairwise_sum(matrix, indices=[1, 0], weights=weights) == 600
