"""Parity tests for TWED (Time-Warp Edit Distance) vs TraMineR TWEDdistance.cpp.

Marteau (2009) DP with the substitution matrix as local cost, integer timestamps,
stiffness ``nu`` and penalty ``h`` (lambda in the C++ source: params[["lambda"]] <- h).
"""

import numpy as np
import pytest

import traminex as tx
from traminex.distances.python_kernels import twed_distance


def test_twed_hand_computed_small_case():
    # a=[0,1], b=[0,2], nu=1, h=1, scalar subm=2. Hand DP gives fmat[2,2] = 2.0.
    a = np.array([0, 1], dtype=np.int32)
    b = np.array([0, 2], dtype=np.int32)
    d = twed_distance(a, b, nu=1.0, lam=1.0, subm=2.0)
    assert d == pytest.approx(2.0)


def test_twed_identical_sequences_zero():
    a = np.array([0, 1, 2, 1, 0], dtype=np.int32)
    assert twed_distance(a, a, nu=1.0, lam=1.0, subm=2.0) == pytest.approx(0.0)


def test_twed_symmetric():
    rng = np.random.default_rng(0)
    a = rng.integers(0, 4, size=8).astype(np.int32)
    b = rng.integers(0, 4, size=8).astype(np.int32)
    dab = twed_distance(a, b, nu=0.5, lam=1.0, subm=2.0)
    dba = twed_distance(b, a, nu=0.5, lam=1.0, subm=2.0)
    assert dab == pytest.approx(dba)


def test_twed_seqdist_pairwise():
    seq = np.array([[0, 0, 1, 1, 2], [0, 1, 1, 2, 2], [2, 2, 1, 0, 0]], dtype=np.int32)
    D = tx.seqdist(seq, method="TWED", nu=0.5, h=1.0, subm=2.0)
    assert D.shape == (3, 3)
    assert np.allclose(np.diag(D), 0.0)
    assert np.allclose(D, D.T)
    assert np.all(D >= 0)


def test_twed_requires_nu():
    seq = np.array([[0, 1, 2], [2, 1, 0]], dtype=np.int32)
    with pytest.raises((ValueError, TypeError)):
        tx.seqdist(seq, method="TWED", h=1.0)


def test_twed_matrix_subm_matches_scalar_when_constant():
    a = np.array([0, 1, 2], dtype=np.int32)
    b = np.array([2, 1, 0], dtype=np.int32)
    K = 3
    M = np.full((K, K), 2.0)
    np.fill_diagonal(M, 0.0)
    d_mat = twed_distance(a, b, nu=0.7, lam=0.5, subm=M)
    d_sca = twed_distance(a, b, nu=0.7, lam=0.5, subm=2.0)
    assert d_mat == pytest.approx(d_sca)
