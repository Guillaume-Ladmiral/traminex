import numpy as np

from traminex.events import EventSequence, events_transition_degrees, frequent_event_pairs_timegap


def test_events_degrees_and_timegap():
    events = EventSequence(
        id=["p1", "p1", "p2", "p2", "p2"],
        time=np.array([1, 2, 1, 3, 4]),
        label=["A", "B", "A", "B", "C"],
    )
    labs, indeg, outdeg = events_transition_degrees(events, normalize=False)
    assert len(labs) == indeg.shape[0] == outdeg.shape[0]
    pairs = frequent_event_pairs_timegap(events, max_gap=1, min_support=1)
    assert isinstance(pairs, list)
