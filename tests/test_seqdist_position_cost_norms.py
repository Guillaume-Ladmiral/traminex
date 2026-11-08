"""R-vocabulary norms with position-dependent OM costs use the position maxima.

``subm_by_pos``/``indel_by_pos`` are a traminex extension (``docs/OMLOC.md``),
whose documented normalization bound uses the largest position-specific
substitution and indel costs. ``seqdist``'s R-vocabulary normalization ignored
them and used the plain scalar ``indel``/``subm``, so "normalized" distances
exceeded 1 -- 2.5 for ``maxdist``/``gmean`` against the documented 0.625
(2026-09-22 audit S1b-02). R's formulas now take the position maxima wherever
they take the indel and the largest substitution cost.
"""

from __future__ import annotations

import numpy as np
import pytest

import traminex as tx
from traminex.distances import pairwise

X = np.array([[0, 0, 1, 1], [2, 2, 0, 0], [1, 2, 2, 1]])
L = X.shape[1]
SUBM_BY_POS = np.stack(
    [np.full((3, 3), cost) - np.diag(np.full(3, cost)) for cost in (1.0, 3.0, 2.0, 6.0)]
)
INDEL_BY_POS = np.array([0.5, 1.0, 4.0, 2.0])
COSTS = {"subm": 2.0, "indel": 1.0, "subm_by_pos": SUBM_BY_POS, "indel_by_pos": INDEL_BY_POS}


def test_maxdist_equals_the_documented_low_level_max_bound():
    documented = pairwise(X, method="OM", norm="max", **COSTS)
    for spelling in ("maxdist", "max"):
        np.testing.assert_allclose(
            tx.seqdist(X, method="OM", norm=spelling, **COSTS), documented, rtol=0, atol=1e-12
        )


def test_maxlength_divides_by_the_largest_position_indel():
    raw = pairwise(X, method="OM", **COSTS)
    expected = raw / (L * INDEL_BY_POS.max())
    np.testing.assert_allclose(
        tx.seqdist(X, method="OM", norm="maxlength", **COSTS), expected, rtol=0, atol=1e-12
    )


@pytest.mark.parametrize("norm", ["maxlength", "gmean", "maxdist", "YujianBo"])
def test_normalized_distances_stay_within_zero_and_one(norm):
    result = tx.seqdist(X, method="OM", norm=norm, **COSTS)
    assert np.all(result >= 0) and np.all(result <= 1 + 1e-12), result
