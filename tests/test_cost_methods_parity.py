"""Parity tests for substitution-cost methods added for TraMineR seqcost.R parity:
FEATURES (Gower-based) and FUTURE (chi-squared on common future / transition rates).

Reference formulas are replicated independently from seqcost.R.
"""

import numpy as np
import pytest

import traminex as tx
from traminex.distances import seqsubm
from traminex.distances.costs import future_costs


def _future_reference(seqmat: np.ndarray) -> np.ndarray:
    """Replicate seqcost.R chisqdista(seqtrate(...)) for method='FUTURE'."""
    K = int(seqmat.max()) + 1
    a = seqmat[:, :-1].ravel()
    b = seqmat[:, 1:].ravel()
    C = np.zeros((K, K), dtype=float)
    np.add.at(C, (a, b), 1.0)
    rs = C.sum(axis=1, keepdims=True)
    T = np.divide(C, rs, out=np.zeros_like(C), where=rs > 0)
    cs = T.sum(axis=0)
    Pdot = np.where(cs == 0, 0.0, 1.0 / np.where(cs == 0, 1.0, cs))
    D = np.zeros((K, K), dtype=float)
    for i in range(K):
        for j in range(i + 1, K):
            d = float(np.sqrt(np.sum(Pdot * (T[i] - T[j]) ** 2)))
            D[i, j] = D[j, i] = d
    return D


def test_future_costs_matches_chisqdista_reference():
    seqmat = np.array(
        [[0, 0, 1, 1, 2], [0, 1, 1, 2, 2], [2, 2, 1, 0, 0], [0, 0, 0, 1, 2]],
        dtype=np.int32,
    )
    indel, sm = future_costs(seqmat)
    ref = _future_reference(seqmat)
    assert np.allclose(sm, ref)
    assert np.allclose(np.diag(sm), 0.0)
    assert np.allclose(sm, sm.T)
    # seqcost.R: ret$indel <- 0.5 * max(costs)
    assert indel == pytest.approx(0.5 * float(ref.max()))


def test_seqsubm_future_returns_cost_matrix():
    seqmat = np.array([[0, 0, 1, 1, 2], [0, 1, 1, 2, 2], [2, 2, 1, 0, 0]], dtype=np.int32)
    M = seqsubm(seqmat, method="FUTURE")
    _, ref = future_costs(seqmat)
    assert np.allclose(M, ref)


def test_seqsubm_features_gower_numeric():
    # 3 states, one numeric feature -> gower = |delta|/range
    seqmat = np.array([[0, 1, 2, 0], [2, 1, 0, 1]], dtype=np.int32)
    state_features = np.array([[0.0], [5.0], [10.0]])
    M = seqsubm(seqmat, method="FEATURES", state_features=state_features, feature_types=["numeric"])
    assert M.shape == (3, 3)
    assert np.allclose(np.diag(M), 0.0)
    assert M[0, 2] == pytest.approx(1.0)  # |0-10|/10
    assert M[0, 1] == pytest.approx(0.5)
    assert M[1, 2] == pytest.approx(0.5)


def test_seqsubm_features_requires_one_row_per_state():
    seqmat = np.array([[0, 1, 2], [2, 1, 0]], dtype=np.int32)  # 3 states
    bad_features = np.array([[0.0], [1.0]])  # only 2 rows
    with pytest.raises(ValueError, match="one row per state"):
        seqsubm(seqmat, method="FEATURES", state_features=bad_features, feature_types=["numeric"])


def test_seqsubm_features_pipeline_into_om():
    seqmat = np.array([[0, 0, 1, 1, 2], [0, 1, 1, 2, 2], [2, 2, 1, 0, 0]], dtype=np.int32)
    sf = np.array([[0.0, "x"], [1.0, "y"], [2.0, "x"]], dtype=object)
    M = seqsubm(seqmat, method="FEATURES", state_features=sf, feature_types=["numeric", "nominal"])
    D = tx.seqdist(seqmat, method="OM", indel=0.5 * float(M.max()), subm=M)
    assert D.shape == (3, 3)
    assert np.allclose(D, D.T)
    assert np.all(D >= 0)
