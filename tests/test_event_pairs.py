import numpy as np

from traminex.events import EventSequence, frequent_event_pairs


def test_frequent_event_pairs():
    events = EventSequence(
        id=["p1", "p1", "p2", "p2", "p2"],
        time=np.array([1, 2, 1, 2, 3]),
        label=["A", "B", "A", "B", "C"],
    )
    pairs = frequent_event_pairs(events, min_support=1)
    assert isinstance(pairs, list)
    assert len(pairs) > 0
