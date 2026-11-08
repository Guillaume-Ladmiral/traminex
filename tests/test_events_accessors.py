"""Behavioural tests for the RE-EXPRESSED event/state type guards.

R's ``is.eseq`` / ``is.seqelist`` / ``is.stslist`` S3 predicates are re-expressed
Pythonically as ``isinstance`` checks; ``levels`` is the object's ``dictionary``
attribute (pinned numerically by test_r_goldens_events_obj.py). No R golden --
these are type idioms, not numeric ports.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from traminex import (
    is_eseq,
    is_seqelist,
    is_stslist,
    seqdef,
    seqecontain,
    seqecreate,
    seqedist,
    seqefsub,
    seqeid,
    seqelength,
    seqentrans,
    seqeweight,
)
from traminex.events._duration import _subset_seqelist


def _eseq():
    return seqecreate(
        pd.DataFrame({"id": ["a", "a", "b"], "timestamp": [0, 1, 0], "event": ["Y", "X", "X"]})
    )


def test_is_seqelist_and_is_eseq():
    e = _eseq()
    assert is_seqelist(e)
    assert not is_eseq(e)
    assert not is_stslist(e)
    # An individual element is an eseq, not a seqelist.
    assert is_eseq(e[0])
    assert not is_seqelist(e[0])


def test_is_stslist():
    sd = seqdef(np.array([["A", "B"], ["B", "A"]], dtype=object))
    assert is_stslist(sd)
    assert not is_seqelist(sd)
    assert not is_eseq(sd)


def test_guards_reject_foreign_objects():
    assert not is_eseq(42)
    assert not is_seqelist("not a seqelist")
    assert not is_stslist(None)


def test_levels_is_sorted_dictionary():
    # R levels() == the sorted-unique event dictionary.
    assert list(_eseq().dictionary) == ["X", "Y"]


def test_seqelist_slice_preserves_container_metadata_and_copy_semantics():
    events = _eseq()
    events.weights = [2.0, 3.0]
    events.lengths = [4.0, 5.0]

    selected = events[::-1]

    assert is_seqelist(selected)
    assert selected.dictionary == events.dictionary
    assert selected.id == ["b", "a", "a"]
    np.testing.assert_array_equal(selected.time, [0.0, 0.0, 1.0])
    assert selected.label == ["X", "Y", "X"]
    np.testing.assert_array_equal(seqeid(selected), [2, 1])
    np.testing.assert_array_equal(seqeweight(selected), [3.0, 2.0])
    np.testing.assert_array_equal(seqelength(selected), [5.0, 4.0])

    selected.weights = [7.0, 8.0]
    selected.lengths = [9.0, 10.0]
    np.testing.assert_array_equal(seqeweight(events), [2.0, 3.0])
    np.testing.assert_array_equal(seqelength(events), [4.0, 5.0])


def test_seqelist_slice_survives_event_analysis_consumer_chain():
    selected = _eseq()[::-1]

    distances = seqedist(selected, idcost=[1.0, 1.0], vparam=0.5)
    np.testing.assert_array_equal(distances.shape, [2, 2])
    np.testing.assert_array_equal(seqecontain(selected, ["X"]), [True, True])

    subsequences = seqefsub(selected, min_support=1, max_k=2)
    transitions = seqentrans(subsequences)
    assert len(subsequences.data) > 0
    assert len(transitions.data) == len(subsequences.data)


def test_seqelist_slice_preserves_rebuilt_group_flat_views():
    # Group subsets rebuild their flat views from the selected sequences
    # (cycle-8 audit F-16): a counting helper reading them must see the
    # group's events rather than an empty vector. Only the numeric id code
    # survives at the object level, so the rebuilt flat id is its string.
    grouped = _subset_seqelist(_eseq(), [1])

    selected = grouped[:]

    assert is_seqelist(selected)
    assert len(selected) == 1
    assert selected.id == ["2"]
    np.testing.assert_array_equal(selected.time, [0.0])
    assert selected.label == ["X"]
    np.testing.assert_array_equal(seqeid(selected), [2])
