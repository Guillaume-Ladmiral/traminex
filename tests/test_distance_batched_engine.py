"""The row engine evaluates its dynamic programs for many pairs at once (TX-AUD-20261010-022).

OM with per-state indels (INDELS, INDELSLOG, indel vectors, OMstran), OMslen,
OMspell, OMloc, TWED, LCS on ragged rows and NMS ran one pure-Python dynamic
program per pair, 80-200 times slower than TraMineR, and a ``refseq`` computed
the whole matrix to return one column. The batched kernels in
``traminex.distances._batched`` evaluate the same recurrences cell by cell for a
batch of pairs, so their results must equal the per-pair reference bit for bit;
``python_kernels._USE_BATCHED = False`` runs that reference.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

import traminex as tx
from traminex.distances import python_kernels as pk

ALPHABET = ["A", "B", "C", "D"]
SM4 = np.array([[0, 1, 2, 3], [1, 0, 1, 2], [2, 1, 0, 1], [3, 2, 1, 0]], dtype=float)
ASYM = np.array([[0, 1, 2, 3], [2, 0, 1, 2], [1, 3, 0, 1], [2, 1, 2, 0]], dtype=float)
PROX = np.array(
    [[1, 0.3, 0.2, 0.1], [0.3, 1, 0.3, 0.2], [0.2, 0.3, 1, 0.3], [0.1, 0.2, 0.3, 1]], dtype=float
)


def _random_rows(
    seed: int, n: int, width: int, *, ragged: bool, runs: bool
) -> list[list[str | None]]:
    rng = np.random.default_rng(seed)
    rows = []
    for _ in range(n):
        if runs:  # long spells, so spell-based methods see durations > 1
            row: list[str | None] = []
            while len(row) < width:
                row += [ALPHABET[int(rng.integers(4))]] * int(rng.integers(1, 5))
            row = row[:width]
        else:
            row = [ALPHABET[int(rng.integers(4))] for _ in range(width)]
        if ragged:
            keep = int(rng.integers(1, width + 1))
            row = row[:keep] + [None] * (width - keep)
        rows.append(row)
    rows.append(list(rows[0]))  # a duplicate row, which the asymmetric route deduplicates
    return rows


def _seqs(seed: int, *, ragged: bool, runs: bool = False) -> tx.STSList:
    return tx.seqdef(
        _random_rows(seed, 14, 9, ragged=ragged, runs=runs), alphabet=ALPHABET, right="DEL"
    )


CASES = {
    "om_indelslog": dict(method="OM", subm="INDELSLOG"),
    "om_indel_vector": dict(method="OM", subm=SM4, indel=np.array([1.0, 2.0, 0.5, 1.5])),
    "om_asymmetric": dict(method="OM", subm=ASYM, indel=1.5),
    "om_scalar_maxdist": dict(method="OM", subm=SM4, indel=1.0, norm="maxdist"),
    "lcs_gmean": dict(method="LCS", norm="gmean"),
    "omslen_mean": dict(method="OMslen", subm=SM4, indel=1.0, h=0.5),
    "omslen_gmean_yujianbo": dict(
        method="OMslen", subm=SM4, indel=1.0, link="gmean", norm="YujianBo"
    ),
    "omslen_asymmetric": dict(method="OMslen", subm=ASYM, indel=1.0),
    "omspell": dict(method="OMspell", subm=SM4, indel=1.0, expcost=0.5, tpow=1.5),
    "omspell_asymmetric_maxlength": dict(method="OMspell", subm=ASYM, indel=1.0, norm="maxlength"),
    "omloc": dict(method="OMloc", subm=SM4, expcost=0.4, context=0.3),
    "omloc_gmean": dict(method="OMloc", subm="CONSTANT", expcost=0.25, norm="gmean"),
    "twed": dict(method="TWED", subm=SM4, nu=0.3, h=0.2),
    "twed_asymmetric_yujianbo": dict(method="TWED", subm=ASYM, nu=0.5, h=0.5, norm="YujianBo"),
    "omstran": dict(method="OMstran", subm=SM4, indel=1.0, otto=0.2),
    "omstran_trate": dict(method="OMstran", subm="TRATE", otto=0.5, transindel="prob"),
    "nms": dict(method="NMS"),
    "nms_kweights": dict(method="NMS", kweights=[1.0, 0.5, 0.25, 2.0]),
    "nms_soft": dict(method="NMS", prox=PROX),
}
RAGGED_ONLY = {"lcs_gmean"}  # with Numba, dense LCS runs its kernel, not the row engine


def _both(seqs: tx.STSList, options: dict, monkeypatch) -> tuple[np.ndarray, np.ndarray]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        monkeypatch.setattr(pk, "_USE_BATCHED", False)
        reference = tx.seqdist(seqs, **options)
        monkeypatch.setattr(pk, "_USE_BATCHED", True)
        batched = tx.seqdist(seqs, **options)
    return reference, batched


@pytest.mark.parametrize("ragged", [False, True])
@pytest.mark.parametrize("case", sorted(CASES))
def test_batched_equals_per_pair_reference(case: str, ragged: bool, monkeypatch) -> None:
    if case in RAGGED_ONLY and not ragged and pk._HAVE_NUMBA:
        pytest.skip("dense LCS runs the compiled or Numba kernel, not the row engine")
    seqs = _seqs(20261013 + len(case), ragged=ragged, runs=case.startswith(("omspell", "omslen")))
    reference, batched = _both(seqs, CASES[case], monkeypatch)
    assert np.array_equal(reference, batched, equal_nan=True)


def test_batched_equals_reference_on_empty_and_single_state_rows(monkeypatch) -> None:
    seqs = tx.seqdef(
        [
            ["A", "A", "A"],
            [None, None, None],
            ["B", None, None],
            ["A", "A", "A"],
            ["A", None, None],
        ],
        alphabet=["A", "B"],
        right="DEL",
    )
    for options in (
        dict(method="OM", subm="INDELS"),
        dict(method="TWED", subm=2.0, nu=0.5, h=0.5),
        dict(method="LCS", norm="maxdist"),
        dict(method="OMspell", subm=2.0, indel=1.0),
        dict(method="NMS"),
        dict(method="NMS", prox=np.array([[1.0, 0.5], [0.5, 1.0]])),
    ):
        reference, batched = _both(seqs, options, monkeypatch)
        assert np.array_equal(reference, batched, equal_nan=True), options


@pytest.mark.parametrize(
    "options",
    [dict(method="OM", subm=SM4, indel=1.0), dict(method="OM", subm="CONSTANT", norm="maxlength")],
)
def test_dense_om_without_numba_equals_reference(options: dict, monkeypatch) -> None:
    monkeypatch.setattr(pk, "_HAVE_NUMBA", False)
    reference, batched = _both(_seqs(5, ragged=False), options, monkeypatch)
    assert np.array_equal(reference, batched, equal_nan=True)


@pytest.mark.parametrize("norm", ["none", "gmean", "maxdist"])
def test_dense_lcs_without_numba_equals_reference(norm: str, monkeypatch) -> None:
    monkeypatch.setattr(pk, "_HAVE_NUMBA", False)
    reference, batched = _both(_seqs(5, ragged=False), dict(method="LCS", norm=norm), monkeypatch)
    assert np.array_equal(reference, batched, equal_nan=True)


@pytest.mark.parametrize("soft", [False, True])
def test_batched_nms_counts_keep_the_reference_addition_order(soft: bool) -> None:
    """Counts above 2**53 round, so only the reference's order gives its bits."""
    rng = np.random.default_rng(3)
    rows = []
    for _ in range(8):
        row: list[int] = []
        while len(row) < int(rng.integers(40, 49)):  # unequal lengths: padded batches
            row += [int(rng.integers(2))] * int(rng.integers(3, 12))
        rows.append(np.asarray(row[:48]))
    left, right = np.triu_indices(len(rows))
    lengths = np.asarray([len(row) for row in rows])
    sizes = np.minimum(lengths[left], lengths[right])
    prox = PROX[:2, :2] if soft else None
    if soft:
        batched = pk._batched.soft_counts_raw(rows, left, right, sizes, prox)
    else:
        batched = pk._batched.subseq_counts_raw(rows, left, right, sizes)
    assert max(float(vector.max()) for vector in batched) > 2.0**53
    for i, j, size, vector in zip(left, right, sizes, batched):
        expected = pk.subseq_counts(rows[i], rows[j], int(size), prox=prox)
        assert np.array_equal(vector, expected), (i, j)


REFERENCES = (
    "om_distance",
    "lcs_distance",
    "omslen_distance",
    "omspell_distance",
    "_omloc_distance",
    "twed_distance",
    "subseq_counts",
)


@pytest.mark.parametrize("numba", [False, True])
@pytest.mark.parametrize("ragged", [False, True])
def test_row_engine_does_not_call_the_per_pair_reference(
    ragged: bool, numba: bool, monkeypatch
) -> None:
    """No route reaches a per-pair reference kernel, with Numba dispatch on or off.

    The dispatcher reads ``_HAVE_NUMBA`` at call time, so the ordinary suite also runs
    the core install's path: without Numba, dense LCS called the per-pair reference
    and this test failed only under ``TRAMINEX_NUMBA=off`` (TX-AUD-20261019-022).
    """
    if numba and not pk._NUMBA_IMPORTED:
        pytest.skip("Numba is not installed")
    monkeypatch.setattr(pk, "_HAVE_NUMBA", numba)
    seqs = _seqs(7, ragged=ragged, runs=True)
    expected = {}
    for case, options in CASES.items():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            expected[case] = tx.seqdist(seqs, **options)

    def refuse(*args, **kwargs):
        raise AssertionError("the per-pair reference was called")

    for name in REFERENCES:
        monkeypatch.setattr(pk, name, refuse)
    for case, options in CASES.items():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            np.testing.assert_array_equal(tx.seqdist(seqs, **options), expected[case])


@pytest.mark.parametrize(
    ("refseq", "pairs"),
    [(0, 14), (3, 14), ([[1, 2], [4, 5, 6]], 6)],
)
@pytest.mark.parametrize("case", sorted(CASES))
def test_refseq_evaluates_only_the_pairs_it_returns(
    case: str, refseq, pairs: int, monkeypatch
) -> None:
    """Every row-engine method evaluates at most the pairs a ``refseq`` returns (batch-2 review M-6).

    The bound is ``<=``, not ``==``: a pair of identical rows is answered without a dynamic
    program, and an index set below a method's crossover runs the per-pair reference instead
    (:func:`python_kernels._worth_batching`). With asymmetric costs R returns d(sequence,
    reference), which the mirrored full matrix holds only above its diagonal, so those
    results are compared with the per-pair reference; the set block (set 1 before set 2)
    is the full matrix's in both cases.
    """
    seqs = _seqs(11, ragged=True, runs=True)
    rows = _random_rows(11, 14, 9, ragged=True, runs=True)
    # refseq=0 is the most frequent sequence: here row 1, the only repeated one.
    assert [rows.count(row) for row in rows] == [2] + [1] * 13 + [2]
    options = CASES[case]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        full = tx.seqdist(seqs, **options)
        monkeypatch.setattr(pk, "_USE_BATCHED", False)
        reference = tx.seqdist(seqs, refseq=refseq, **options)
        monkeypatch.setattr(pk, "_USE_BATCHED", True)
    evaluated: list[int] = []
    original = pk._evaluate_raw_batches

    def counting(indices_i, indices_j, raw_batch):
        evaluated.append(len(indices_i))
        return original(indices_i, indices_j, raw_batch)

    monkeypatch.setattr(pk, "_evaluate_raw_batches", counting)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = tx.seqdist(seqs, refseq=refseq, **options)
    assert sum(evaluated) <= pairs
    np.testing.assert_array_equal(result, reference)
    if isinstance(refseq, list):
        np.testing.assert_array_equal(
            result, full[np.ix_(np.asarray(refseq[0]) - 1, np.asarray(refseq[1]) - 1)]
        )
    else:
        assert result.shape == (seqs.nseq,)
        if options.get("subm") is not ASYM:
            np.testing.assert_array_equal(result, full[:, max(refseq - 1, 0)])
