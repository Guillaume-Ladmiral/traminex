"""R-golden parity for DHD/HAM with a caller-supplied substitution matrix.

Phase 1 / A1 (see MIGRATION_PLAN.md) adds an explicit-``subm`` path to
``traminex.distances.pairwise(method="DHD")`` so traminex can faithfully serve
weightedcluster's position-cost-weighted Hamming distances. TraMineR honours a
user ``sm`` in ``seqdist`` (a static ``(K,K)`` matrix behaves as
``method="HAM"``; a time-varying ``(K,K,L)`` array as ``method="DHD"``); traminex
previously ignored it and always rebuilt transition-rate costs.

Goldens are produced by ``scripts/r_goldens/make_goldens_dhd_subm.R`` against
TraMineR 2.2-14 (norm="none", i.e. the raw position-wise substitution sum).
"""

from __future__ import annotations

import numpy as np
from r_goldens_utils import R_GOLDEN_ATOL, assert_matches_golden

from traminex.distances import pairwise

# The 6x8 ABC sequences shared with the other distance goldens, 0-based.
_SEQS = np.array(
    [
        [0, 0, 0, 1, 1, 2, 2, 2],
        [0, 1, 0, 1, 2, 1, 2, 0],
        [2, 2, 2, 2, 0, 0, 0, 0],
        [0, 0, 1, 1, 1, 1, 2, 2],
        [1, 1, 1, 1, 1, 1, 1, 1],
        [2, 0, 2, 0, 1, 2, 1, 0],
    ],
    dtype=np.int32,
)
_K = 3
_L = _SEQS.shape[1]


def test_dhd_static_subm_matches_r_ham_with_sm():
    """A static (K,K) subm reproduces R seqdist(method='HAM', sm=...)."""
    sm = np.full((_K, _K), 2.0)
    np.fill_diagonal(sm, 0.0)
    D = pairwise(_SEQS, method="DHD", subm=sm, norm="none")
    assert_matches_golden(D, "DHD_subm_static_sm2_none")


def test_dhd_timevarying_subm_matches_r_dhd_with_sm():
    """A time-varying (L,K,K) subm reproduces R seqdist(method='DHD', sm=...)."""
    sm_t = np.full((_L, _K, _K), 2.0)
    for t in range(_L):
        np.fill_diagonal(sm_t[t], 0.0)
    # Position index 4 (R's position 5) has expensive off-diagonal costs.
    sm_t[4] = np.array([[0.0, 5.0, 5.0], [5.0, 0.0, 5.0], [5.0, 5.0, 0.0]])
    D = pairwise(_SEQS, method="DHD", subm=sm_t, norm="none")
    assert_matches_golden(D, "DHD_subm_timevarying_none")


def test_dhd_without_subm_still_trate_faithful():
    """Bare DHD (no subm) is unchanged: TraMineR transition-rate time-varying costs."""
    D_default = pairwise(_SEQS, method="DHD")
    # Explicit-subm path must not have altered the no-subm behaviour.
    assert D_default.shape == (_SEQS.shape[0], _SEQS.shape[0])
    assert np.allclose(D_default, D_default.T, rtol=0, atol=R_GOLDEN_ATOL)
    assert np.allclose(np.diag(D_default), 0.0, rtol=0, atol=R_GOLDEN_ATOL)
