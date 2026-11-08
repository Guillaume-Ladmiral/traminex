"""``seqcomp`` and ``seqmpos`` compare states by label across two objects.

TraMineR's ``seqcomp`` (``x[1:lx] == y[1:ly]``) and ``seqmpos``
(``seq1 == seq2``) compare the two objects' factor columns, and R's
``Ops.factor`` compares factor *labels*.  traminex compared the integer codes,
so two objects whose alphabets differ in content or order were compared code
by code: ``A-B-B`` against ``B-A-A`` was "identical" with 3 matching positions,
while ``A-B-B`` coded under ``["B", "A"]`` did not match itself (2026-09-22
audit S3-01).  The per-object void and missing exclusions of ``seqmpos`` apply
to each object's own markers.

Expected values: ``scripts/r_goldens/make_goldens_audit20260922_s3a.R``
(section S3-01), live TraMineR 2.2-14.
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


def _r_value(case: str):
    kind, payload = R[case]
    assert kind == "value"
    if payload in ("TRUE", "FALSE"):
        return payload == "TRUE"
    return int(payload)


def _sts(rows, **kwargs) -> STSList:
    return tx.seqdef(np.array(rows, dtype=object), **kwargs)


def _row(x: STSList, index: int) -> STSList:
    out = x.copy()
    out.data = x.data[index : index + 1]
    return out


OBJECTS = {
    "x": lambda: _sts([["A", "B", "B"]]),
    "y": lambda: _sts([["B", "A", "A"]]),
    "z": lambda: _sts([["A", "B", "B"]], alphabet=["B", "A"]),
    "w": lambda: _sts([["A", "C", "C"]]),
    "xs": lambda: _sts([["A", "B", "B"]], alphabet=["A", "B", "C"]),
    "m": lambda: _sts([["A", None, "B"]]),
    "n": lambda: _sts([["A", "B", "B"]]),
    "a": lambda: _sts([["A", "B", "C"], ["C", "C", "A"]], alphabet=["A", "B", "C"]),
    "b": lambda: _sts([["A", "C", "C"], ["C", "B", "A"]], alphabet=["C", "B", "A"]),
    "v1": lambda: _sts([["A", "B", "B"], ["A", "B", None]]),
    "v2": lambda: _sts([["B", "B", "B"], ["A", "B", None]], alphabet=["B", "A"]),
    # A real state spelled like m's nr text, under a different nr text.
    "st": lambda: _sts([["A", "*", "B"]], alphabet=["A", "*", "B"], nr="?"),
    # A real state spelled like vr's void text, under a different void text.
    "sv": lambda: _sts([["A", "B", "%"]], alphabet=["A", "B", "%"], void="#"),
    "vr": lambda: _sts([["A", "B", None]]),
}


@pytest.mark.parametrize(
    "case,left,right",
    [
        ("s3_01_seqcomp_x_y", "x", "y"),
        ("s3_01_seqcomp_x_z", "x", "z"),
        ("s3_01_seqcomp_x_w", "x", "w"),
        ("s3_01_seqcomp_x_xs", "x", "xs"),
        ("s3_01_seqcomp_m_n", "m", "n"),
        ("s3_01_seqcomp_m_m", "m", "m"),
    ],
)
def test_seqcomp_compares_state_labels_as_r(case, left, right):
    assert tx.seqcomp(OBJECTS[left](), OBJECTS[right]()) is _r_value(case)


def test_seqcomp_ragged_rows_with_different_codebooks_match_r():
    left = _row(OBJECTS["v1"](), 1)
    right = _row(OBJECTS["v2"](), 1)
    assert tx.seqcomp(left, right) is _r_value("s3_01_seqcomp_v1r2_v2r2")


@pytest.mark.parametrize(
    "case,left,right,kwargs",
    [
        ("s3_01_seqmpos_x_y", "x", "y", {}),
        ("s3_01_seqmpos_x_z", "x", "z", {}),
        ("s3_01_seqmpos_x_w", "x", "w", {}),
        ("s3_01_seqmpos_x_xs", "x", "xs", {}),
        ("s3_01_seqmpos_m_n", "m", "n", {}),
        ("s3_01_seqmpos_m_n_wm", "m", "n", {"with_missing": True}),
        ("s3_01_seqmpos_m_m_wm", "m", "m", {"with_missing": True}),
        ("s3_01_seqmpos_a_b", "a", "b", {}),
        ("s3_01_seqmpos_v1_v2", "v1", "v2", {}),
        ("s3_01_seqmpos_st_m", "st", "m", {}),
        ("s3_01_seqmpos_m_st", "m", "st", {}),
        ("s3_01_seqmpos_st_m_wm", "st", "m", {"with_missing": True}),
        ("s3_01_seqmpos_sv_vr", "sv", "vr", {}),
        ("s3_01_seqmpos_vr_sv", "vr", "sv", {}),
    ],
)
def test_seqmpos_counts_matching_state_labels_as_r(case, left, right, kwargs):
    result = tx.seqmpos(OBJECTS[left](), OBJECTS[right](), **kwargs)
    assert isinstance(result, int)
    assert result == _r_value(case)


def test_both_marker_representations_compare_identically():
    """``A-*-B-%`` held with external marker codes equals its in-alphabet twin."""
    external = STSList(
        np.array([[0, 3, 1, 4]]), alphabet=["A", "B"], nr="*", void="%", nr_code=3, void_code=4
    )
    twin = STSList(
        np.array([[0, 2, 1, 3]]),
        alphabet=["A", "B", "*", "%"],
        nr="*",
        void="%",
        nr_code=2,
        void_code=3,
    )
    # Same alphabet as ``twin`` but reordered, so every code differs.
    permuted = STSList(
        np.array([[3, 1, 2, 0]]),
        alphabet=["%", "*", "B", "A"],
        nr="*",
        void="%",
        nr_code=1,
        void_code=0,
    )
    objects = (external, twin, permuted)
    for left in objects:
        for right in objects:
            assert tx.seqcomp(left, right) is True
            assert tx.seqmpos(left, right) == 2
            assert tx.seqmpos(left, right, with_missing=True) == 3
