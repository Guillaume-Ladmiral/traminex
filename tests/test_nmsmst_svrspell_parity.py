"""Parity tests for NMSMST (NMSMSTdistance.cpp) and SVRspell (NMSDURSoftdistance.cpp).

These kernels count, over every pair of matching length-(k+1) subsequence embeddings
(p in a, q in b with a[p]==b[q] elementwise):
- NMSMST: sum of minimum shared times  sum_l min(dur_a[p_l], dur_b[q_l])
- SVRspell: product of summed times     (sum_l dur_a[p_l]) * (sum_l dur_b[q_l])     (identity softmatch)

We validate the C++ ports against independent brute-force embedding enumeration.
"""

from itertools import combinations

import numpy as np
import pytest

import traminex as tx
from traminex.distances._spells import spell_encode
from traminex.distances.python_kernels import nmsmst_counts, svrspell_counts


def _nmsmst_ref(a, da, b, db, ksize):
    kv = np.zeros(ksize)
    for length in range(1, ksize + 1):
        tot = 0.0
        for p in combinations(range(len(a)), length):
            for q in combinations(range(len(b)), length):
                if all(a[p[l]] == b[q[l]] for l in range(length)):
                    tot += sum(min(da[p[l]], db[q[l]]) for l in range(length))
        kv[length - 1] = tot
    return kv


def _svr_ref(a, da, b, db, ksize):
    kv = np.zeros(ksize)
    for length in range(1, ksize + 1):
        tot = 0.0
        for p in combinations(range(len(a)), length):
            for q in combinations(range(len(b)), length):
                if all(a[p[l]] == b[q[l]] for l in range(length)):
                    wa = sum(da[p[l]] for l in range(length))
                    wb = sum(db[q[l]] for l in range(length))
                    tot += wa * wb
        kv[length - 1] = tot
    return kv


@pytest.mark.parametrize(
    "a,da,b,db",
    [
        ([0, 1, 2], [1.0, 1.0, 1.0], [0, 2, 1], [1.0, 1.0, 1.0]),
        ([0, 1], [2.0, 3.0], [0, 1], [1.0, 4.0]),
        ([0, 1, 0], [2.0, 1.0, 3.0], [0, 1, 2], [1.0, 2.0, 1.0]),
    ],
)
def test_nmsmst_counts_match_bruteforce(a, da, b, db):
    a = np.array(a)
    b = np.array(b)
    da = np.array(da)
    db = np.array(db)
    ks = min(len(a), len(b))
    assert np.allclose(nmsmst_counts(a, b, da, db, ks), _nmsmst_ref(a, da, b, db, ks))


@pytest.mark.parametrize(
    "a,da,b,db",
    [
        ([0, 1, 2], [1.0, 1.0, 1.0], [0, 2, 1], [1.0, 1.0, 1.0]),
        ([0, 1], [2.0, 3.0], [0, 1], [1.0, 4.0]),
        ([0, 1, 0], [2.0, 1.0, 3.0], [0, 1, 2], [1.0, 2.0, 1.0]),
    ],
)
def test_svrspell_counts_match_bruteforce_identity(a, da, b, db):
    a = np.array(a)
    b = np.array(b)
    da = np.array(da)
    db = np.array(db)
    ks = min(len(a), len(b))
    K = int(max(a.max(), b.max())) + 1
    got = svrspell_counts(a, b, da, db, ks, np.eye(K))
    assert np.allclose(got, _svr_ref(a, da, b, db, ks))


def test_nmsmst_seqdist_matches_bruteforce_distance():
    seqmat = np.array([[0, 0, 1, 2], [0, 1, 1, 2], [2, 1, 0, 0]], dtype=np.int32)
    D = tx.seqdist(seqmat, method="NMSMST")
    # brute-force distance via spell encoding, kweights all 1
    spells = [spell_encode(seqmat[i]) for i in range(3)]
    seqs = [s for s, _ in spells]
    durs = [d for _, d in spells]

    def selfval(i):
        m = len(seqs[i])
        return _nmsmst_ref(seqs[i], durs[i], seqs[i], durs[i], m).sum()

    for i in range(3):
        for j in range(i + 1, 3):
            ks = min(len(seqs[i]), len(seqs[j]))
            common = _nmsmst_ref(seqs[i], durs[i], seqs[j], durs[j], ks).sum()
            # final sqrt applied by TraMineR's R wrapper (seqdist.R)
            expected = np.sqrt(selfval(i) + selfval(j) - 2 * common)
            assert D[i, j] == pytest.approx(expected)


@pytest.mark.parametrize("method", ["NMSMST", "SVRspell"])
def test_nmsmst_svrspell_properties_with_spells(method):
    seqmat = np.array(
        [[0, 0, 0, 1, 2], [0, 1, 1, 2, 2], [2, 2, 1, 1, 0], [0, 0, 1, 1, 1]], dtype=np.int32
    )
    D = tx.seqdist(seqmat, method=method)
    assert D.shape == (4, 4)
    assert np.allclose(np.diag(D), 0.0)
    assert np.allclose(D, D.T)
    assert np.all(D >= -1e-9)


@pytest.mark.parametrize("method", ["NMSMST", "SVRspell"])
def test_nmsmst_svrspell_identical_zero(method):
    seq = np.array([[0, 0, 1, 2, 2], [0, 0, 1, 2, 2]], dtype=np.int32)
    D = tx.seqdist(seq, method=method)
    assert D[0, 1] == pytest.approx(0.0)


def test_compiled_spell_counts_equal_the_python_kernels_bit_for_bit():
    # 2026-10-25 study, finding F-06: with Numba dispatched, nmsmst_counts and svrspell_counts run
    # compiled copies of the Python loops; their counts must be the same bytes, overflow included.
    from traminex.distances import python_kernels as pk
    from traminex.distances._spells import spell_length_power

    if not pk._NUMBA_IMPORTED:
        pytest.skip("numba is not installed")
    rng = np.random.default_rng(20261107)
    for trial in range(400):
        n_states = int(rng.integers(1, 6))
        (a, da), (b, db) = (
            spell_encode(rng.integers(0, n_states, size=int(rng.integers(0, 30)))) for _ in range(2)
        )
        tpow = float(rng.choice([1.0, 0.5, 2.0, 0.0, -1.0, 400.0]))
        da, db = spell_length_power(da, tpow), spell_length_power(db, tpow)
        ks = min(len(a), len(b))
        soft = np.eye(n_states)
        if trial % 2:
            soft = np.where(soft > 0, 1.0, rng.uniform(0.0, 1.0, size=soft.shape))
        arrays = pk._spell_kernel_arrays(a, b, da, db)
        with np.errstate(over="ignore", invalid="ignore"):
            pairs = (
                (pk._nmsmst_counts(a, b, da, db, ks, None), pk._nmsmst_counts_nb(*arrays, ks)),
                (
                    pk._svrspell_counts(a, b, da, db, ks, soft, None),
                    pk._svrspell_counts_nb(*arrays, ks, soft),
                ),
            )
        for reference, compiled in pairs:
            assert reference.tobytes() == compiled.tobytes(), trial


@pytest.mark.parametrize("method", ["NMSMST", "SVRspell"])
def test_seqdist_spell_methods_are_the_same_bytes_with_and_without_numba(method, monkeypatch):
    from traminex.distances import python_kernels as pk

    if not pk._NUMBA_IMPORTED:
        pytest.skip("numba is not installed")
    rng = np.random.default_rng(20261108)
    seqmat = rng.integers(0, 4, size=(12, 25)).astype(np.int32)
    monkeypatch.setattr(pk, "_HAVE_NUMBA", True)
    compiled = tx.seqdist(seqmat, method=method, tpow=0.5)
    monkeypatch.setattr(pk, "_HAVE_NUMBA", False)
    python = tx.seqdist(seqmat, method=method, tpow=0.5)
    assert compiled.tobytes() == python.tobytes()
