"""Event-sequence gaps are printed as R's ``format()`` prints them.

``Sequence::sprint`` renders each gap, and the trailing gap to the end of
observation, by calling R's ``format()`` (``TMRNumberFormat``,
``eventseq.cpp:90-100``), which keeps seven significant digits and switches
to scientific notation whenever that is narrower: 100000 prints as ``1e+05``
and 0.0001 as ``1e-04``, while 123456 stays fixed and 1/3 prints as
``0.3333333``. traminex used Python's ``%.7g``, printing ``100000`` and
``0.0001``, so ``as_character`` and the row labels of ``seqeapplysub`` differed
from R (2026-09-22 audit S5-08).

Expected values: ``scripts/r_goldens/make_goldens_audit20260922_s5.R``
(TraMineR 2.2-14).
"""

from __future__ import annotations

from pathlib import Path

import pytest

import traminex as tx

GOLDEN = Path(__file__).parent / "r_goldens" / "audit20260922_s5_event_objects.txt"
R = {
    fields[1]: fields[2].split(" ")
    for fields in (line.split("\t") for line in GOLDEN.read_text().splitlines())
    if fields[0] == "gap_format"
}


@pytest.fixture
def eseq():
    return tx.seqecreate(
        id=[1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5],
        timestamp=[0, 100000, 0, 0.0001, 3, 0, 123456, 0, 1 / 3, 0, 1234567.8],
        event=["a", "b", "a", "c", "b", "a", "b", "a", "b", "a", "b"],
    )


def test_gaps_print_like_r(eseq):
    assert eseq.as_character() == R["created"]


def test_the_trailing_gap_prints_like_r(eseq):
    eseq.lengths = [200000, 5, 123460, 1, 2e6]
    assert eseq.as_character() == R["with_length"]


def test_seqeapplysub_rows_carry_the_same_labels(eseq):
    eseq.lengths = [200000, 5, 123460, 1, 2e6]
    counts = tx.seqeapplysub(tx.seqefsub(eseq, str_subseq=["(a)-(b)"]), method="presence")
    assert list(counts.index) == R["applysub_rows"]
