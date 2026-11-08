"""A handful of pairs uses the per-pair reference (batch-2 review M-4, TX-AUD-20261010-022).

The batched engine sweeps every cell of the dynamic-programming table with a
few NumPy calls, whatever the number of pairs, so for one pair of length 300 it
took 0.8 to 1.9 seconds where the per-pair reference takes 0.13 to 0.23 (OM,
OMslen, OMloc, OMstran, TWED, LCS, OMspell). Below each method's measured
crossover, the row engine now evaluates the per-pair reference, which gives the
same bits. NMS stays batched: its batched counts are faster even for one pair.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

import traminex as tx
from traminex.distances import python_kernels as pk

SM4 = np.array([[0, 1, 2, 3], [1, 0, 1, 2], [2, 1, 0, 1], [3, 2, 1, 0]], float)
ASYM = np.array([[0, 1, 2, 3], [2, 0, 1, 2], [2, 1, 0, 1], [3, 2, 2, 0]], float)

# method options and whether rows are of unequal length (LCS batches only then)
CASES = {
    "OM_INDELSLOG": (dict(method="OM", subm="INDELSLOG"), False),
    "OM_asymmetric": (dict(method="OM", subm=ASYM, indel=1.0), False),
    "LCS_ragged": (dict(method="LCS"), True),
    "OMslen": (dict(method="OMslen", subm=SM4, indel=1.0), False),
    "OMspell": (dict(method="OMspell", subm=SM4, indel=1.0), False),
    "OMloc": (dict(method="OMloc", subm="CONSTANT", expcost=0.5), False),
    "OMstran": (dict(method="OMstran", subm="CONSTANT", indel=1, otto=0.1), False),
    "TWED": (dict(method="TWED", subm=SM4, nu=0.5, h=0.5), False),
}


def _seqs(n: int, ragged: bool) -> tx.STSList:
    rng = np.random.default_rng(3)
    rows = []
    for k in range(n):
        target = 12 - (k % 4 if ragged else 0)
        row: list[object] = []
        while len(row) < target:
            row += [str("ABCD"[int(rng.integers(4))])] * int(rng.integers(1, 4))
        rows.append(row[:target] + [None] * (12 - target))
    return tx.seqdef(rows, alphabet=list("ABCD"), right="DEL")


def _distances(seqs, options, monkeypatch, **extra) -> tuple[np.ndarray, int]:
    calls: list[int] = []
    original = pk._evaluate_raw_batches

    def counting(left, right, raw_batch):
        calls.append(len(left))
        return original(left, right, raw_batch)

    monkeypatch.setattr(pk, "_evaluate_raw_batches", counting)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = np.asarray(tx.seqdist(seqs, **options, **extra))
    monkeypatch.setattr(pk, "_evaluate_raw_batches", original)
    return result, len(calls)


def _reference(seqs, options, monkeypatch, **extra) -> np.ndarray:
    monkeypatch.setattr(pk, "_USE_BATCHED", False)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = np.asarray(tx.seqdist(seqs, **options, **extra))
    monkeypatch.setattr(pk, "_USE_BATCHED", True)
    return result


@pytest.mark.parametrize("case", sorted(CASES))
def test_one_pair_uses_the_per_pair_reference(case: str, monkeypatch) -> None:
    options, ragged = CASES[case]
    seqs = _seqs(2, ragged)
    result, batched_calls = _distances(seqs, options, monkeypatch)
    assert batched_calls == 0
    np.testing.assert_array_equal(result, _reference(seqs, options, monkeypatch))


@pytest.mark.parametrize("case", sorted(CASES))
def test_many_pairs_stay_batched(case: str, monkeypatch) -> None:
    options, ragged = CASES[case]
    seqs = _seqs(8, ragged)  # 28 pairs, above every crossover
    result, batched_calls = _distances(seqs, options, monkeypatch)
    assert batched_calls >= 1
    np.testing.assert_array_equal(result, _reference(seqs, options, monkeypatch))


def test_a_refseq_pair_uses_the_per_pair_reference(monkeypatch) -> None:
    options, _ = CASES["TWED"]
    seqs = _seqs(8, False)
    result, batched_calls = _distances(seqs, options, monkeypatch, refseq=[[1], [2]])
    assert batched_calls == 0
    full = _reference(seqs, options, monkeypatch)
    np.testing.assert_array_equal(result, full[np.ix_([0], [1])])


@pytest.mark.parametrize("prox", [False, True], ids=["hard", "soft"])
def test_nms_stays_batched_for_one_pair(prox: bool, monkeypatch) -> None:
    options = dict(method="NMS")
    if prox:
        matrix = np.full((4, 4), 0.3)
        np.fill_diagonal(matrix, 1.0)
        options["prox"] = matrix
    seqs = _seqs(2, False)
    result, batched_calls = _distances(seqs, options, monkeypatch)
    assert batched_calls == 1
    np.testing.assert_array_equal(result, _reference(seqs, options, monkeypatch))
