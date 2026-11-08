"""Parity tests for NMS (Number of Matching Subsequences) vs TraMineR NMSdistance.cpp.

kvect[k] = number of common length-(k+1) subsequences (Elzinga count kernel). With
distMethod=2 and kweights all 1, dist = sqrt(Ival + Jval - 2*Aval) where I/J are the
self subsequence counts and A the common counts. The final sqrt is applied by
TraMineR's R wrapper (seqdist.R), on top of what NMSdistance.cpp computes.

Hand cases:
- a=[0,1], b=[0,1] identical -> 0.
- a=[0,1], b=[1,0]: self=[2,1] each (k0=2, k1=1), common=[2,0] -> sqrt(3+3-2*2) = sqrt(2).
"""

import numpy as np
import pytest

import traminex as tx
from traminex.distances.python_kernels import subseq_counts


def test_subseq_counts_basic():
    a = np.array([0, 1])
    # self-counts of [0,1]: one length-1 of each symbol (2), one length-2 "01" (1)
    kv = subseq_counts(a, a, 2)
    assert list(kv) == [2.0, 1.0]


def test_subseq_counts_common_disjoint_order():
    a = np.array([0, 1])
    b = np.array([1, 0])
    kv = subseq_counts(a, b, 2)
    # length-1 common: "0" (1x1) + "1" (1x1) = 2; length-2 common: none
    assert list(kv) == [2.0, 0.0]


def test_nms_identical_zero():
    seq = np.array([[0, 1, 2, 1], [0, 1, 2, 1]], dtype=np.int32)
    D = tx.seqdist(seq, method="NMS")
    assert D[0, 1] == pytest.approx(0.0)


def test_nms_hand_case_two_sequences():
    seq = np.array([[0, 1], [1, 0]], dtype=np.int32)
    D = tx.seqdist(seq, method="NMS")
    assert D[0, 1] == pytest.approx(np.sqrt(2.0))


def test_nms_properties():
    seq = np.array(
        [[0, 0, 1, 1, 2], [0, 1, 1, 2, 2], [2, 2, 1, 0, 0], [0, 1, 2, 1, 0]], dtype=np.int32
    )
    D = tx.seqdist(seq, method="NMS")
    assert D.shape == (4, 4)
    assert np.allclose(np.diag(D), 0.0)
    assert np.allclose(D, D.T)
    assert np.all(D >= -1e-9)


def test_nms_soft_match_identity_equals_hard():
    seq = np.array([[0, 1, 2], [2, 1, 0], [0, 0, 1]], dtype=np.int32)
    K = 3
    ident = np.eye(K, dtype=float)
    D_hard = tx.seqdist(seq, method="NMS")
    D_soft = tx.seqdist(seq, method="NMS", prox=ident)
    assert np.allclose(D_hard, D_soft)


def test_nms_soft_match_uses_traminer_soft_recurrence():
    """A non-identity proximity must use NMSMSTSoftdistanceII, not hard NMS."""
    seq = np.array([[0, 0], [1, 0]], dtype=np.int32)
    prox = np.array([[1.0, 0.5], [0.5, 1.0]])

    distance = tx.seqdist(seq, method="NMS", prox=prox)

    assert distance[0, 1] == pytest.approx(np.sqrt(2.0))


@pytest.mark.parametrize("method", ["NMS", "SVRspell"])
@pytest.mark.parametrize("size", [2, 4])
def test_soft_subsequence_proximity_requires_exact_state_dimensions(method, size):
    seq = np.array([[0, 1], [1, 2], [2, 0]], dtype=np.int32)

    with pytest.raises(ValueError, match=r"prox.*3.*3"):
        tx.seqdist(seq, method=method, prox=np.eye(size))


@pytest.mark.parametrize("method", ["NMS", "SVRspell"])
def test_soft_subsequence_proximity_rejects_nonfinite_values(method):
    seq = np.array([[0, 1], [1, 2], [2, 0]], dtype=np.int32)

    with pytest.raises(ValueError, match=r"prox.*finite"):
        tx.seqdist(seq, method=method, prox=np.full((3, 3), np.nan))


@pytest.mark.parametrize("method", ["NMS", "SVRspell"])
def test_soft_subsequence_proximity_reports_indefinite_matrix(method):
    seq = np.array([[0, 1], [1, 2], [2, 0]], dtype=np.int32)
    prox = np.array([[1.0, 2.0, 0.0], [2.0, 1.0, 0.0], [0.0, 0.0, 1.0]])

    with pytest.warns(UserWarning, match="not positive semi-definite"):
        distance = tx.seqdist(seq, method=method, prox=prox)

    assert np.all(np.isfinite(distance))
