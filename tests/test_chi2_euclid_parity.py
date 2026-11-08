"""Parity tests for CHI2 (pdotj-weighted) and EUCLID distances vs TraMineR seqdist-CHI2.R.

The exact references use a single interval (step = L), where each sequence's vector is
its state time-proportions P[i] and pdotj is the overall (mean) state proportion.

    CHI2 (norm=False):  d(i,j) = sqrt( sum_k (P[i,k]-P[j,k])^2 / pdotj[k] )   over pdotj>0
    EUCLID (norm=False): d(i,j) = sqrt( sum_k (P[i,k]-P[j,k])^2 )            over present states
"""

import numpy as np
import pytest

import traminex as tx


def _proportions(seqmat, K):
    n, L = seqmat.shape
    P = np.zeros((n, K), dtype=float)
    for i in range(n):
        c = np.bincount(seqmat[i], minlength=K).astype(float)
        P[i] = c / c.sum()
    return P


def _chi2_single_ref(seqmat):
    K = int(seqmat.max()) + 1
    P = _proportions(seqmat, K)
    pdotj = P.mean(axis=0)
    keep = pdotj > 0
    n = seqmat.shape[0]
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            diff = P[i, keep] - P[j, keep]
            D[i, j] = D[j, i] = np.sqrt(np.sum(diff * diff / pdotj[keep]))
    return D


def _euclid_single_ref(seqmat):
    K = int(seqmat.max()) + 1
    P = _proportions(seqmat, K)
    present = P.sum(axis=0) > 0
    n = seqmat.shape[0]
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            diff = P[i, present] - P[j, present]
            D[i, j] = D[j, i] = np.sqrt(np.sum(diff * diff))
    return D


SEQ = np.array([[0, 0, 1, 1, 2], [0, 1, 1, 2, 2], [2, 2, 1, 0, 0], [0, 0, 0, 1, 2]], dtype=np.int32)


def test_chi2_single_interval_matches_pdotj_weighted_reference():
    D = tx.seqdist(SEQ, method="CHI2", step=SEQ.shape[1], norm="none")
    assert np.allclose(D, _chi2_single_ref(SEQ))


def test_euclid_single_interval_matches_euclidean_reference():
    D = tx.seqdist(SEQ, method="EUCLID", step=SEQ.shape[1], norm="none")
    assert np.allclose(D, _euclid_single_ref(SEQ))


def test_euclid_norm_divides_single_interval_by_sqrt2():
    # EUCLID norm sets the divisor maxd=2 (single interval -> /sqrt(2))
    D_raw = tx.seqdist(SEQ, method="EUCLID", step=SEQ.shape[1], norm="none")
    D_norm = tx.seqdist(SEQ, method="EUCLID", step=SEQ.shape[1], norm="auto")
    assert np.allclose(D_norm, D_raw / np.sqrt(2.0))


@pytest.mark.parametrize("method", ["CHI2", "EUCLID"])
def test_chi2_euclid_properties_multi_interval(method):
    D = tx.seqdist(SEQ, method=method, step=2, overlap=False)
    assert D.shape == (4, 4)
    assert np.allclose(np.diag(D), 0.0)
    assert np.allclose(D, D.T)
    assert np.all(D >= 0)


@pytest.mark.parametrize("method", ["CHI2", "EUCLID"])
def test_chi2_euclid_identical_rows_zero(method):
    X = np.array([[0, 1, 2, 1], [0, 1, 2, 1]], dtype=np.int32)
    D = tx.seqdist(X, method=method, step=2, overlap=True)
    assert D[0, 1] == pytest.approx(0.0, abs=1e-12)


def test_chi2_default_step_is_per_position():
    # R parity: default step=1 (one interval per position)
    D_default = tx.seqdist(SEQ, method="CHI2")
    D_step1 = tx.seqdist(SEQ, method="CHI2", step=1)
    assert np.allclose(D_default, D_step1)
