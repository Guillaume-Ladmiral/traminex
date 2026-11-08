"""Regressions for the 2026-09-02 audit's event-layer findings."""

from __future__ import annotations

import pytest

from traminex.events import (
    events_to_state_sequences,
    events_transition_matrix,
    frequent_event_pairs,
    frequent_events,
    seqecreate,
)
from traminex.events._helpers import EventSequence


def test_fractional_event_times_reach_the_state_grid():
    events = EventSequence(id=["p1", "p1"], time=[0.0, 2.5], label=["home", "work"])
    sts = events_to_state_sequences(events)
    assert any("work" in row for row in sts.to_strings())


def test_integer_event_times_keep_the_historical_grid():
    events = EventSequence(id=["p1", "p1"], time=[0.0, 2.0], label=["home", "work"])
    sts = events_to_state_sequences(events)
    assert sts.seqlen == 3
    assert sts.to_strings() == ["home-home-work"]


def test_seqelist_flat_views_match_the_per_sequence_events():
    # The counting helpers deliberately duck-type SeqElist (pinned by the
    # phase-1 tie-order tests), so the flat id/time/label views must describe
    # the same events as the per-sequence objects: no end_event rows, and
    # group subsets must expose their group's events, not empty vectors.
    elist = seqecreate(
        id=["p1", "p1", "p1", "p2", "p2"],
        timestamp=[1.0, 2.0, 3.0, 1.0, 2.0],
        event=["hire", "promo", "END", "hire", "END"],
        end_event="END",
    )
    assert "END" not in elist.label
    assert len(elist.label) == sum(len(s.codes) for s in elist.sequences)
    counted = dict((label, count) for label, count, _ in frequent_events(elist))
    assert counted == {"hire": 2, "promo": 1}

    from traminex.events._duration import _subset_seqelist

    subset = _subset_seqelist(elist, [0])
    assert list(subset.label) == ["hire", "promo"]
    assert frequent_events(subset) == [("hire", 1, 0.5), ("promo", 1, 0.5)]


def test_event_sequence_revalidates_on_mutation():
    events = EventSequence(id=["p1", "p1"], time=[0.0, 1.0], label=["a", "b"])
    with pytest.raises(ValueError, match="same length"):
        events.label = ["a"]
    with pytest.raises((TypeError, ValueError)):
        events.time = "not-a-vector"
    events.label = ["a", "c"]  # consistent updates stay allowed
    assert events.label == ["a", "c"]


def test_state_sequence_rows_and_helper_orderings_are_documented():
    # Sorted row order is the pinned contract for events_to_state_sequences
    # (test_events_to_state_sequences_uses_event_ids_as_fresh_row_ids); the
    # counting helpers group by first appearance. Both sides now say so.
    events = EventSequence(
        id=["zeta", "zeta", "alpha", "alpha"],
        time=[0.0, 1.0, 0.0, 1.0],
        label=["x", "y", "u", "v"],
    )
    sts = events_to_state_sequences(events)
    assert sts.ids == ["alpha", "zeta"]
    assert "sorted" in (events_to_state_sequences.__doc__ or "")
    assert "typo" in (events_transition_matrix.__doc__ or "").lower()


def test_support_fraction_denominators_are_documented():
    from traminex.events import _helpers

    assert "denominator" in (frequent_events.__doc__ or "").lower()
    assert "denominator" in (frequent_event_pairs.__doc__ or "").lower()
    assert "denominator" in (_helpers.frequent_event_pairs_timegap.__doc__ or "").lower()
