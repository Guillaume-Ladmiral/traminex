"""TX-AUD-20261011-036: the public weighted-inertia helpers check their indices.

``weighted_pairwise_sum(diss, indices, weights)`` and ``weighted_inertia(diss, indices, weights)`` take zero-based
integer case positions, in range and not repeated. They accepted float positions (truncated), negative ones (counted
from the end), repeated ones and a float weight vector given in the indices' place (``weighted_inertia(D, w)``), and
returned a number for each. The signatures stay: the weightedcluster sibling calls them positionally.
"""

from __future__ import annotations

import numpy as np
import pytest

import traminex as tx

D = np.abs(np.subtract.outer(np.arange(6.0), np.arange(6.0))) ** 1.5
W = np.array([1.0, 2.0, 0.5, 1.5, 1.0, 3.0])
HELPERS = [tx.weighted_pairwise_sum, tx.weighted_inertia]


@pytest.mark.parametrize("function", HELPERS)
@pytest.mark.parametrize(
    "indices, error",
    [
        ([1.7, 0.2], TypeError),
        (np.array([0.0, 1.0]), TypeError),
        (np.array([True, False, True, False, True, False]), TypeError),
        ([-1, 0], ValueError),
        ([0, 6], ValueError),
        ([0, 0, 1], ValueError),
    ],
)
def test_invalid_indices_are_refused(function, indices, error) -> None:
    with pytest.raises(error, match="indices"):
        function(D, indices, W)


def test_weights_given_in_the_indices_place_are_refused() -> None:
    with pytest.raises(TypeError, match="indices"):
        tx.weighted_inertia(D, W)


@pytest.mark.parametrize("function", HELPERS)
def test_valid_indices_are_accepted(function) -> None:
    every = function(D, None, W)
    assert function(D, list(range(6)), W) == every
    assert function(D, np.arange(6, dtype=np.int32), W) == every
    assert function(D, [], W) == 0.0
    assert function(D, [4], W) == 0.0
    # The order of the positions is R's summation order, so it may change the last bits only.
    np.testing.assert_allclose(function(D, [5, 0, 2], W), function(D, [0, 2, 5], W), rtol=1e-15)
