"""``seqalign(include_missing=True)`` sizes ``sm`` with the missing state.

``seqalign.R:9-14`` names the rows and columns of ``sm`` with
``alphabet(seqdata)`` plus ``nr`` whenever ``with.missing=TRUE``, so it takes
the ``(K + 1) x (K + 1)`` matrix that ``seqsubm(with.missing=TRUE)`` builds and
refuses a ``K x K`` one, even on data without a missing cell.  traminex added
the missing state only when the object carried an nr code, so on missing-free
data it rejected its own ``seqsubm(..., with_missing=True)`` output ("sm must
be 3x3") and accepted the matrix R refuses (2026-09-22 audit S3-05).

Expected values: ``scripts/r_goldens/make_goldens_audit20260922_s3a.R``
(section S3-05), live TraMineR 2.2-14.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import traminex as tx
from traminex.sequences import STSList

GOLDENS = Path(__file__).parent / "r_goldens" / "audit20260922_s3a.txt"


def _r() -> dict[str, tuple[str, str]]:
    out = {}
    for line in GOLDENS.read_text().splitlines():
        if line.startswith("#"):
            continue
        case, kind, payload, *_ = line.split("\t")
        out[case] = (kind, payload)
    return out


R = _r()


def _r_payload(case: str) -> list[str]:
    kind, payload = R[case]
    assert kind == "value", R[case]
    return payload.split(" ")


def _sa() -> STSList:
    rows = [["A", "B", "C"], ["A", "C", "C"], ["C", "A", "B"]]
    return tx.seqdef(np.array(rows, dtype=object), alphabet=["A", "B", "C"])


SA_VAR = np.array([0, 1, 3, 5, 1, 0, 2, 5, 3, 2, 0, 5, 5, 5, 5, 0], dtype=float).reshape(4, 4).T


def _check(result, case: str, rows: int, cols: int) -> None:
    assert list(result.operation) == _r_payload(f"{case}_operation")
    assert list(result.seq1) == _r_payload(f"{case}_seq1")
    assert list(result.seq2) == _r_payload(f"{case}_seq2")
    np.testing.assert_allclose(
        result.cost, np.array(_r_payload(f"{case}_cost"), dtype=float), rtol=0, atol=0
    )
    expected = np.array(_r_payload(f"{case}_costmatrix"), dtype=float).reshape(cols, rows).T
    np.testing.assert_allclose(result.cost_matrix, expected, rtol=0, atol=0)


def test_seqsubm_with_missing_builds_rs_four_by_four_matrix():
    sm = np.asarray(tx.seqsubm(_sa(), "CONSTANT", 2.0, with_missing=True), dtype=float)
    expected = np.array(_r_payload("s3_05_align_const_sm"), dtype=float).reshape(4, 4).T
    np.testing.assert_array_equal(sm, expected)


def test_seqalign_accepts_the_with_missing_cost_matrix_as_r():
    sm = np.asarray(tx.seqsubm(_sa(), "CONSTANT", 2.0, with_missing=True), dtype=float)
    result = tx.seqalign(_sa(), (0, 1), indel=1.0, sm=sm, include_missing=True)
    _check(result, "s3_05_align_const_12", 4, 4)


def test_seqalign_indexes_a_variable_cost_matrix_as_r():
    result = tx.seqalign(_sa(), (0, 2), indel=4.0, sm=SA_VAR, include_missing=True)
    _check(result, "s3_05_align_var_13", 4, 4)


def test_seqalign_refuses_a_k_by_k_matrix_with_missing_as_r():
    assert R["s3_05_align_k_by_k_sm"][0] == "error"
    sm = np.asarray(tx.seqsubm(_sa(), "CONSTANT", 2.0), dtype=float)
    with pytest.raises(ValueError, match="sm must be 4×4 matrix"):
        tx.seqalign(_sa(), (0, 1), indel=1.0, sm=sm, include_missing=True)


def test_default_costs_follow_the_same_size():
    result = tx.seqalign(_sa(), (0, 1), indel=1.0, include_missing=True)
    _check(result, "s3_05_align_const_12", 4, 4)


def test_without_missing_the_matrix_stays_k_by_k():
    sm = np.asarray(tx.seqsubm(_sa(), "CONSTANT", 2.0), dtype=float)
    result = tx.seqalign(_sa(), (0, 1), indel=1.0, sm=sm)
    assert list(result.operation) == _r_payload("s3_05_align_const_12_operation")
    with pytest.raises(ValueError, match="sm must be 3×3 matrix"):
        tx.seqalign(_sa(), (0, 1), indel=1.0, sm=np.zeros((4, 4)))
