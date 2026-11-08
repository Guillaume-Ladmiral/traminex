"""R-golden tests for gower_matrix (RP4 T15).

Anchors the (squared, weights)-parametrised double-centring to live TraMineR
TraMineR 2.2-14 ``gower_matrix`` (see scripts/r_goldens/make_goldens_gower.R).
"""

from __future__ import annotations

import numpy as np
import pytest
from r_goldens_utils import R_GOLDEN_ATOL, load_golden

from traminex import gower_matrix


def _inputs():
    D = load_golden("phase2_input_D")
    w = load_golden("phase2_input_weights_real").ravel()
    return D, w


# tag -> (squared, weighted)
_CASES = {
    "sq_uw": (True, False),
    "raw_uw": (False, False),
    "sq_w": (True, True),
    "raw_w": (False, True),
}


@pytest.mark.parametrize("tag", list(_CASES))
def test_gower_matrix_matches_live_r(tag):
    squared, weighted = _CASES[tag]
    D, w = _inputs()
    weights = w if weighted else None
    out = gower_matrix(D, squared=squared, weights=weights)
    golden = load_golden(f"gower_{tag}")
    np.testing.assert_allclose(out, golden, rtol=0, atol=R_GOLDEN_ATOL)


def test_gower_matrix_default_is_squared_unweighted_double_centering():
    # Back-compat: gower_matrix(D) keeps the classical -0.5 J D^2 J MDS Gram.
    D, _ = _inputs()
    n = D.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    expected = -0.5 * (J @ (D**2) @ J)
    np.testing.assert_allclose(gower_matrix(D), expected, rtol=0, atol=R_GOLDEN_ATOL)
