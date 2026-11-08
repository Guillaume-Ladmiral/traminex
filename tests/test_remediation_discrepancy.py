"""TX-AUD-20261005-017: self distances never enter weighted discrepancy.

The TraMineR inertia kernels accumulate only pairs with j > i. Expectations
below are hand sums for distances 1, 2, 3 and case weights 1, 2, 3, also
scaled by 1/4 to exercise fractional weights. The pinned R probe independently
confirms variance 26/36 and centres (11, 17, 11)/18.
"""

from __future__ import annotations

import numpy as np
import pytest

from traminex.diss import disscenter, dissvar, weighted_inertia, weighted_pairwise_sum

DIAGONALS = [0.0, 7.0, np.nan, np.inf, 1e200]


def _distances(diagonal):
    distances = np.array([[0.0, 1.0, 2.0], [1.0, 0.0, 3.0], [2.0, 3.0, 0.0]])
    np.fill_diagonal(distances, diagonal)
    return distances


@pytest.mark.parametrize("diagonal", DIAGONALS)
@pytest.mark.parametrize(
    "weights,pair_sum,inertia",
    [([1.0, 2.0, 3.0], 26.0, 26 / 6), ([0.25, 0.5, 0.75], 26 / 16, 26 / 24)],
)
def test_weighted_sum_and_inertia_ignore_self_distances(diagonal, weights, pair_sum, inertia):
    distances = _distances(diagonal)
    original = distances.copy()
    assert weighted_pairwise_sum(distances, weights=weights) == pytest.approx(pair_sum)
    assert weighted_inertia(distances, weights=weights) == pytest.approx(inertia)
    np.testing.assert_array_equal(distances, original)


@pytest.mark.parametrize("diagonal", DIAGONALS)
@pytest.mark.parametrize("weights", [[1.0, 2.0, 3.0], [0.25, 0.5, 0.75]])
@pytest.mark.parametrize("squared,expected", [(False, 26 / 36), (True, 68 / 36)])
def test_weighted_variance_ignores_self_distances(diagonal, weights, squared, expected):
    distances = _distances(diagonal)
    original = distances.copy()
    # R silently overflows when squaring a huge diagonal before its C kernel
    # ignores the resulting Inf; this test concerns the analytical result.
    with np.errstate(over="ignore"):
        assert dissvar(distances, weights=weights, squared=squared) == pytest.approx(expected)
    np.testing.assert_array_equal(distances, original)


@pytest.mark.parametrize("diagonal", DIAGONALS)
@pytest.mark.parametrize(
    "indices,pair_sum,inertia",
    [([0, 2], 6.0, 1.5), ([2, 0], 6.0, 1.5), ([1], 0.0, 0.0), ([], 0.0, 0.0)],
)
def test_subset_inertia_ignores_self_distances(diagonal, indices, pair_sum, inertia):
    distances = _distances(diagonal)
    distances.setflags(write=False)
    assert weighted_pairwise_sum(distances, indices, [1.0, 2.0, 3.0]) == pytest.approx(pair_sum)
    assert weighted_inertia(distances, indices, [1.0, 2.0, 3.0]) == pytest.approx(inertia)


@pytest.mark.parametrize("diagonal", DIAGONALS)
@pytest.mark.parametrize("weights", [[1.0, 2.0, 3.0], [0.25, 0.5, 0.75]])
@pytest.mark.parametrize(
    "squared,expected",
    [(False, [11 / 18, 17 / 18, 11 / 18]), (True, [4 / 9, 25 / 9, 16 / 9])],
)
def test_group_centres_ignore_self_distances(diagonal, weights, squared, expected):
    distances = _distances(diagonal)
    original = distances.copy()
    with np.errstate(over="ignore"):
        result = disscenter(distances, weights=weights, squared=squared)
    np.testing.assert_allclose(result, expected, rtol=0, atol=1e-14)
    np.testing.assert_array_equal(distances, original)


@pytest.mark.parametrize("diagonal", DIAGONALS)
@pytest.mark.parametrize("weights", [[1.0, 2.0, 3.0], [0.25, 0.5, 0.75]])
def test_grouped_centres_and_singletons_ignore_self_distances(diagonal, weights):
    distances = _distances(diagonal)
    distances.setflags(write=False)
    result = disscenter(distances, group=["a", "a", "b"], weights=weights)
    np.testing.assert_allclose(result, [4 / 9, 1 / 9, 0.0], rtol=0, atol=1e-14)


@pytest.mark.parametrize("diagonal", DIAGONALS)
@pytest.mark.parametrize("weights", [[1.0, 2.0, 3.0], [0.25, 0.5, 0.75]])
def test_all_centres_ignore_self_distances(diagonal, weights):
    distances = _distances(diagonal)
    distances.setflags(write=False)
    result = disscenter(distances, group=["a", "a", "b"], weights=weights, allcenter=True)
    assert list(result.columns) == ["a", "b"]
    np.testing.assert_allclose(
        result.to_numpy(), [[4 / 9, 2.0], [1 / 9, 3.0], [22 / 9, 0.0]], rtol=0, atol=1e-14
    )


@pytest.mark.parametrize("diagonal", DIAGONALS)
def test_medoid_selection_ignores_self_distances(diagonal):
    distances = _distances(diagonal)
    distances.setflags(write=False)
    assert disscenter(distances, weights=[1.0, 2.0, 3.0], medoids_index="first") == 0
    np.testing.assert_array_equal(
        disscenter(distances, weights=[1.0, 2.0, 3.0], medoids_index="all"), [0, 2]
    )
    assert disscenter(
        distances, group=["a", "a", "b"], weights=[1.0, 2.0, 3.0], medoids_index="first"
    ) == {"a": 1, "b": 2}


@pytest.mark.parametrize("off_diagonal", [np.nan, np.inf])
@pytest.mark.parametrize("function", [weighted_pairwise_sum, weighted_inertia, dissvar])
def test_nonfinite_off_diagonal_is_preserved(off_diagonal, function):
    distances = _distances(0.0)
    distances[0, 1] = distances[1, 0] = off_diagonal
    result = function(distances, weights=[1.0, 2.0, 3.0])
    if np.isnan(off_diagonal):
        assert np.isnan(result)
    else:
        assert np.isposinf(result)


@pytest.mark.parametrize("function", [weighted_pairwise_sum, weighted_inertia, dissvar])
def test_zero_weight_does_not_hide_infinite_off_diagonal(function):
    distances = _distances(0.0)
    distances[0, 1] = distances[1, 0] = np.inf
    assert np.isnan(function(distances, weights=[0.0, 2.0, 3.0]))


def test_unweighted_variance_keeps_r_whole_matrix_shortcut():
    # TraMineR's weights=NULL branch includes the diagonal in sum(diss).
    assert dissvar(_distances(2.0)) == pytest.approx(1.0)
