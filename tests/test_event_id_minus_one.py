"""An event sequence whose id is -1 is an ordinary sequence (TX-AUD-20261019-006).

The event engine used -1, a valid id, as two sentinels ported from TraMineR's C++: "no sequence
counted yet" (TreeEventNode::lastID) and "generic pattern" (Sequence::isGeneric). A sequence
with id -1 was therefore never counted and printed without its gaps; R crashes on the same
input ("attempt access index 0/0 in VECTOR_ELT"). Expected values are worked by hand: both
sequences hold A at time 0 and B at time 1, so (A), (B) and (A)-(B) each occur in both.
"""

from __future__ import annotations

import pytest

import traminex as tx

TIMES = [0, 1, 0, 1]
EVENTS = ["A", "B", "A", "B"]


def _sequences(ids):
    return tx.seqecreate(id=ids, timestamp=TIMES, event=EVENTS)


@pytest.mark.parametrize("ids", [[-1, -1, 2, 2], [2, 2, -1, -1], [1, 1, 2, 2]])
def test_subsequences_are_counted_for_an_id_of_minus_one(ids) -> None:
    result = tx.seqefsub(_sequences(ids), min_support=2)
    assert list(result.subseq) == ["(A)", "(A)-(B)", "(B)"]
    assert result.data["Support"].tolist() == [1.0, 1.0, 1.0]
    assert result.data["Count"].tolist() == [2.0, 2.0, 2.0]


def test_count_method_two_is_unchanged() -> None:
    constraint = tx.seqeconstraint(count_method=2)
    minus_one = tx.seqefsub(_sequences([-1, -1, 2, 2]), min_support=2, constraint=constraint)
    ordinary = tx.seqefsub(_sequences([1, 1, 2, 2]), min_support=2, constraint=constraint)
    assert list(minus_one.subseq) == list(ordinary.subseq)
    assert minus_one.data.to_dict("list") == ordinary.data.to_dict("list")


def test_a_sequence_with_id_minus_one_prints_its_gaps() -> None:
    assert [str(sequence) for sequence in _sequences([-1, -1, 2, 2])] == ["(A)-1-(B)", "(A)-1-(B)"]
