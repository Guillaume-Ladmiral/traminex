"""Stored event vectors cannot be changed in place behind their object's back.

``EventSequence`` re-validates the whole (id, time, label) triple on every
assignment, and cycle-10 N-19 froze its ``time`` array because an in-place
write went straight past that validation. The ``id`` and ``label`` lists had
the same hole: ``ev.label.append("C")`` left four labels against three times
(and ``frequent_events`` counted the extra label), and ``ev.id[0] = 123``
stored a non-string the constructor rejects. ``SeqElist``'s flat
``id``/``time``/``label`` views -- properties with no setter, which must
describe the same events as the per-sequence objects (the cycle-8 F-16
contract) -- handed out their stored list and array, so a write changed what
``events_transition_matrix`` saw while ``seqefsub`` kept reading the
sequences (2026-09-22 audit S5-05). The stored vectors are now read-only;
reading them is unchanged and assignment stays the way to change an
``EventSequence``.
"""

from __future__ import annotations

import copy
import pickle

import numpy as np
import pytest

import traminex as tx
from traminex.events import SeqElist


def _events() -> tx.EventSequence:
    return tx.EventSequence(
        id=["p", "p", "q"], time=np.array([0.0, 1.0, 0.0]), label=["A", "B", "A"]
    )


@pytest.mark.parametrize(
    "write",
    [
        lambda ev: ev.label.append("C"),
        lambda ev: ev.label.extend(["C"]),
        lambda ev: ev.label.insert(0, "C"),
        lambda ev: ev.label.pop(),
        lambda ev: ev.label.remove("A"),
        lambda ev: ev.label.clear(),
        lambda ev: ev.label.sort(),
        lambda ev: ev.label.reverse(),
        lambda ev: ev.label.__setitem__(0, "Z"),
        lambda ev: ev.label.__delitem__(0),
        lambda ev: ev.id.__setitem__(0, 123),
        lambda ev: ev.id.append("r"),
    ],
)
def test_event_sequence_id_and_label_refuse_in_place_writes(write):
    events = _events()
    with pytest.raises(TypeError, match="read-only"):
        write(events)
    assert events.id == ["p", "p", "q"]
    assert events.label == ["A", "B", "A"]
    assert tx.frequent_events(events) == [("A", 2, 2 / 3), ("B", 1, 1 / 3)]


def test_event_sequence_vectors_still_read_and_reassign_as_lists():
    events = _events()
    assert isinstance(events.label, list) and events.label == ["A", "B", "A"]
    assert events.label + ["C"] == ["A", "B", "A", "C"]
    events.label = ["X", "Y", "Z"]
    assert events.label == ["X", "Y", "Z"]
    with pytest.raises(ValueError, match="same length"):
        events.label = ["X"]
    with pytest.raises(ValueError, match="same length"):
        events.label += ["W"]
    # copy, deepcopy and pickle rebuild the vectors without tripping the guard
    for clone in (copy.copy(events), copy.deepcopy(events), pickle.loads(pickle.dumps(events))):
        assert clone.id == ["p", "p", "q"] and clone.label == ["X", "Y", "Z"]


def test_seqelist_flat_views_refuse_in_place_writes():
    eseq = tx.seqecreate(id=[1, 1, 2], timestamp=[0, 1, 0], event=["a", "b", "a"])
    before = tx.events_transition_matrix(eseq)
    with pytest.raises(ValueError, match="read-only"):
        eseq.time[0] = 99.0
    with pytest.raises(TypeError, match="read-only"):
        eseq.label[0] = "b"
    with pytest.raises(TypeError, match="read-only"):
        eseq.id.append("3")
    after = tx.events_transition_matrix(eseq)
    assert before[0] == after[0]
    np.testing.assert_array_equal(before[1], after[1])
    assert eseq.label == ["a", "b", "a"] and eseq.id == ["1", "1", "2"]
    np.testing.assert_array_equal(eseq.time, [0.0, 1.0, 0.0])
    # a slice's views are read-only too
    with pytest.raises(ValueError, match="read-only"):
        eseq[0:1].time[0] = 5.0


def test_building_a_seqelist_does_not_freeze_the_caller_s_vectors():
    eseq = tx.seqecreate(id=[1, 1], timestamp=[0, 1], event=["a", "b"])
    times = np.array([0.0, 1.0])
    ids = ["1", "1"]
    rebuilt = SeqElist(
        eseq.sequences, eseq.dictionary, flat_id=ids, flat_time=times, flat_label=["a", "b"]
    )
    times[0] = 7.0
    ids.append("2")
    np.testing.assert_array_equal(rebuilt.time, [0.0, 1.0])
    assert rebuilt.id == ["1", "1"]
