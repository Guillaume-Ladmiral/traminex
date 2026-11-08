"""``seqecontain`` also tests the subsequences of a ``seqefsub`` result.

``seqecontain.R:6`` starts with ``if (is.subseqelist(eseq)) eseq <- eseq$subseq``,
so a subsequence list is accepted and each subsequence is tested for the
given events, with the same any/all (``unknown.exclude``) rule. traminex
raised ``TypeError`` for it (2026-09-22 audit S5-06).

Expected values: ``scripts/r_goldens/make_goldens_audit20260922_s5.R``
(TraMineR 2.2-14).
"""

from __future__ import annotations

from pathlib import Path

import pytest

import traminex as tx

GOLDEN = Path(__file__).parent / "r_goldens" / "audit20260922_s5_event_objects.txt"
LINES = [line.split("\t") for line in GOLDEN.read_text().splitlines() if not line.startswith("#")]
(SUBSEQ,) = [fields[1] for fields in LINES if fields[0] == "contain_subseq"]
CASES = [fields[1:] for fields in LINES if fields[0] == "contain"]


def _subsequences() -> tx.SubSeqElist:
    eseq = tx.seqecreate(
        id=[1, 1, 2, 2, 2], timestamp=[0, 5, 0, 1, 3], event=["a", "b", "a", "c", "b"]
    )
    return tx.seqefsub(eseq, min_support=1)


def test_the_subsequences_are_r_s():
    assert _subsequences().subseq == SUBSEQ.split(" ")


@pytest.mark.parametrize("events,exclude,expected", CASES)
def test_each_subsequence_is_tested_for_the_events(events, exclude, expected):
    result = tx.seqecontain(_subsequences(), events.split(" "), unknown_exclude=(exclude == "TRUE"))
    assert result.tolist() == [value == "TRUE" for value in expected.split(" ")]
