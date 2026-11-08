"""Parity tests for OMslen (OMvdistance.cpp) and OMspell (OMPerdistance.cpp).

Key invariants used as ground truth (R semantics; see also
tests/test_r_goldens_distances.py for live TraMineR 2.2-14 goldens):
- OMslen weights each position by spell_length^(-h) (seqdist.R: dur.mat^(-1*h),
  h default 0.5). With all spells of length 1 the weights are 1, both
  link="gmean" (sm*sqrt(w_i*w_j)) and link="mean"
  (sm/2*(w_i+w_j), per seqdist.R and OMvdistance.h) reduce to plain OM.
- OMspell runs on the DSS with durations dur^tpow - 1 (seqdist.R subtracts 1);
  with timecost=0 and no consecutive repeats (DSS == sequence) OMspell == OM.
Plus a hand-computed OMspell DP case.
"""

import numpy as np
import pytest

import traminex as tx
from traminex.distances._spells import spell_encode, spell_length_per_position
from traminex.distances.python_kernels import omslen_distance, omspell_distance


def test_spell_encode_basic():
    s, d = spell_encode(np.array([0, 0, 1, 2, 2]))
    assert list(s) == [0, 1, 2]
    assert list(d) == [2.0, 1.0, 2.0]
    assert list(spell_length_per_position(np.array([0, 0, 1, 2, 2]))) == [2.0, 2.0, 1.0, 2.0, 2.0]


def _no_repeat_rows(n, L, K, seed):
    rng = np.random.default_rng(seed)
    rows = []
    for _ in range(n):
        row = [rng.integers(0, K)]
        while len(row) < L:
            nxt = rng.integers(0, K)
            if nxt != row[-1]:
                row.append(nxt)
        rows.append(row)
    return np.array(rows, dtype=np.int32)


def test_omslen_equals_om_when_no_consecutive_repeats():
    seqmat = _no_repeat_rows(5, 7, 3, seed=1)
    indel, subm = 1.0, 2.0
    D_om = tx.seqdist(seqmat, method="OM", indel=indel, subm=subm)
    # unit spells -> weights 1: gmean substitution cost = sm (equals OM) ...
    D_gmean = tx.seqdist(seqmat, method="OMslen", indel=indel, subm=subm, link="gmean")
    assert np.allclose(D_gmean, D_om), "OMslen link=gmean should equal OM with unit spells"
    # ... and mean passes sm/2 before applying the sum w_i + w_j = 2.
    D_mean = tx.seqdist(seqmat, method="OMslen", indel=indel, subm=subm, link="mean")
    assert np.allclose(D_mean, D_om), "OMslen link=mean should equal OM for unit spells"


def test_omspell_equals_om_when_timecost_zero_and_no_repeats():
    seqmat = _no_repeat_rows(5, 7, 3, seed=2)
    indel, subm = 1.0, 2.0
    D_om = tx.seqdist(seqmat, method="OM", indel=indel, subm=subm)
    D_spell = tx.seqdist(seqmat, method="OMspell", indel=indel, subm=subm, expcost=0.0)
    assert np.allclose(D_spell, D_om)


def test_omspell_hand_computed_case():
    # a=[0,0,1] DSS([0,1],[2,1]); b=[0,1,1] DSS([0,1],[1,2]); indel=1, subm=2, timecost=0.5
    a = np.array([0, 0, 1], dtype=np.int32)
    b = np.array([0, 1, 1], dtype=np.int32)
    d = omspell_distance(a, b, indel=1.0, subm=2.0, timecost=0.5)
    assert d == pytest.approx(1.0)


@pytest.mark.parametrize("method,kw", [("OMslen", {"link": "mean"}), ("OMspell", {"expcost": 0.5})])
def test_omspell_omslen_properties(method, kw):
    seqmat = np.array(
        [[0, 0, 1, 1, 2], [0, 1, 1, 2, 2], [2, 2, 1, 0, 0], [0, 0, 0, 1, 2]], dtype=np.int32
    )
    D = tx.seqdist(seqmat, method=method, indel=1.0, subm=2.0, **kw)
    assert D.shape == (4, 4)
    assert np.allclose(np.diag(D), 0.0)
    assert np.allclose(D, D.T)
    assert np.all(D >= 0)


def test_omslen_differs_from_om_with_repeats():
    seqmat = np.array([[0, 0, 0, 1, 1], [0, 1, 1, 1, 1]], dtype=np.int32)
    D_om = tx.seqdist(seqmat, method="OM", indel=1.0, subm=2.0)
    D_slen = tx.seqdist(seqmat, method="OMslen", indel=1.0, subm=2.0, link="mean")
    assert not np.allclose(D_slen, D_om)


def test_omslen_symmetric_kernel():
    a = np.array([0, 0, 1, 2], dtype=np.int32)
    b = np.array([0, 1, 1, 2], dtype=np.int32)
    da = omslen_distance(
        a,
        b,
        indel=1.0,
        subm=2.0,
        link="gmean",
        dur_a=spell_length_per_position(a),
        dur_b=spell_length_per_position(b),
    )
    db = omslen_distance(
        b,
        a,
        indel=1.0,
        subm=2.0,
        link="gmean",
        dur_a=spell_length_per_position(b),
        dur_b=spell_length_per_position(a),
    )
    assert da == pytest.approx(db)


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ Empty\ sequence\(s\)\ \ \d+\ \ may\ cause\ inconsistent\ results!:UserWarning"
)
@pytest.mark.parametrize(
    "method,kwargs", [("OMslen", {"link": "mean"}), ("OMspell", {"expcost": 0.5})]
)
def test_infinite_vector_indel_makes_every_nonidentity_pair_nan(method, kwargs):
    sequences = tx.seqdef(
        np.asarray([["B", "%"], ["B", "B"], ["%", "%"]], dtype=object),
        alphabet=["A", "B"],
    )

    actual = tx.seqdist(
        sequences,
        method=method,
        subm=2.0,
        indel=np.asarray([np.inf, 1.0]),
        norm="none",
        **kwargs,
    )

    np.testing.assert_array_equal(np.diag(actual), np.zeros(3))
    assert np.isnan(actual[np.triu_indices(3, 1)]).all()
    assert np.isnan(actual[np.tril_indices(3, -1)]).all()
