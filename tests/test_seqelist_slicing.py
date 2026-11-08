"""Slicing a ``SeqElist`` keeps exactly the events of the selected sequences.

R's ``[.seqelist`` returns the selected event sequences, and with them exactly
their events (``seqe2TSE`` of the subset). traminex rebuilt the flat
``id``/``time``/``label`` views of a slice by cutting the flat rows into runs
of identical id text and assuming one run per sequence. A sequence with no
event (an id whose only row is the end event) has no run, so ``el[1:3]``
raised ``IndexError`` and ``el[0:2]`` returned the next sequence's events;
fractional ids that R's ``as.integer`` merges into one sequence (1.2 and 1.7)
make two runs, so ``el[1:2]`` returned the wrong event (2026-09-22 audit
S5-03). The flat rows now remember which sequence they belong to.

Expected values: ``scripts/r_goldens/make_goldens_audit20260922_s5.R``
(TraMineR 2.2-14).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import traminex as tx
from traminex.events import SeqElist

GOLDEN = Path(__file__).parent / "r_goldens" / "audit20260922_s5_event_objects.txt"
SLICES = [
    line.split("\t")[1:] for line in GOLDEN.read_text().splitlines() if line.startswith("slice\t")
]


def _source(name: str) -> SeqElist:
    if name == "with_end":
        return tx.seqecreate(
            id=[1, 1, 2, 3, 3],
            timestamp=[0, 1, 5, 0, 2],
            event=["a", "b", "END", "a", "c"],
            end_event="END",
        )
    return tx.seqecreate(id=[1.2, 1.7, 2.0], timestamp=[0, 1, 0], event=["a", "b", "c"])


@pytest.mark.parametrize("name,r_index,ids,lengths,events", SLICES)
def test_a_slice_keeps_exactly_the_events_of_its_sequences(name, r_index, ids, lengths, events):
    positions = [int(value) - 1 for value in r_index.split(" ")]
    assert positions == list(range(positions[0], positions[-1] + 1))
    sliced = _source(name)[positions[0] : positions[-1] + 1]

    assert tx.seqeid(sliced).tolist() == [int(value) for value in ids.split(" ")]
    np.testing.assert_array_equal(tx.seqelength(sliced), np.array(lengths.split(" "), float))
    rows = [] if events == "none" else [row.split(",") for row in events.split(" ")]
    # The flat id keeps the caller's id text; R reports as.integer of it.
    assert [int(float(value)) for value in sliced.id] == [int(row[0]) for row in rows]
    np.testing.assert_array_equal(sliced.time, np.array([row[1] for row in rows], dtype=float))
    assert sliced.label == [row[2] for row in rows]


def test_a_slice_of_a_slice_still_maps_rows_to_sequences():
    eseq = _source("with_end")
    twice = eseq[1:3][1:2]
    assert tx.seqeid(twice).tolist() == [3]
    assert twice.id == ["3", "3"]
    assert twice.label == ["a", "c"]


def test_flat_views_that_cannot_be_matched_to_sequences_are_refused():
    """A hand-built list whose flat rows do not map to its sequences is refused."""
    eseq = _source("with_end")
    rebuilt = SeqElist(
        eseq.sequences,
        eseq.dictionary,
        flat_id=list(eseq.id),
        flat_time=np.array(eseq.time),
        flat_label=list(eseq.label),
    )
    with pytest.raises(ValueError, match="flat"):
        rebuilt[0:2]


def test_a_group_subset_slices_exactly_too():
    """The seqedplot group views are rebuilt per sequence and slice exactly."""
    from traminex.events._duration import _subset_seqelist

    subset = _subset_seqelist(_source("with_end"), [1, 2])
    sliced = subset[1:2]
    assert tx.seqeid(sliced).tolist() == [3]
    assert sliced.label == ["a", "c"]
    assert subset[0:1].label == []
