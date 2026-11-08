# Example: event analytics (weighted degrees and time-gap pairs)
import numpy as np

from traminex.events import (
    EventSequence,
    events_transition_degrees,
    events_transition_matrix,
    frequent_event_pairs_timegap,
)

# Build small example events
ids = ["p1", "p1", "p1", "p2", "p2", "p3"]
times = np.array([1, 2, 5, 1, 2, 1])
labels = ["A", "B", "C", "A", "C", "B"]
E = EventSequence(id=ids, time=times, label=labels)

labs, M = events_transition_matrix(E, normalize=False)
print("Labels:", labs)
print("Transition matrix:\n", M)

labs, indeg, outdeg = events_transition_degrees(E, normalize=False, id_weights={"p2": 2.0})
print("Weighted outdegrees:", outdeg)

pairs = frequent_event_pairs_timegap(E, max_gap=2, min_support=1)
print("Time-gap pairs:", pairs)
