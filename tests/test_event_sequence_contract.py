"""Boundary and cross-helper contracts for the flat event convenience layer."""

from __future__ import annotations

import numpy as np
import pytest

from traminex.events import (
    EventSequence,
    events_to_state_sequences,
    events_transition_degrees,
    events_transition_matrix,
    frequent_event_pairs,
    frequent_event_pairs_timegap,
    frequent_events,
)


@pytest.mark.parametrize(
    ("ids", "times", "labels"),
    [
        (["p1", "p2"], [0.0], ["A"]),
        (["p1"], [0.0, 1.0], ["A"]),
        (["p1"], [0.0], ["A", "B"]),
    ],
)
def test_event_sequence_refuses_every_unequal_length_permutation(ids, times, labels):
    with pytest.raises(ValueError, match="id, time, and label must have the same length"):
        EventSequence(id=ids, time=np.asarray(times), label=labels)


@pytest.mark.parametrize(
    ("times", "error", "message"),
    [
        (np.array(0.0), ValueError, "time must be one-dimensional"),
        (np.array([[0.0]]), ValueError, "time must be one-dimensional"),
        (np.array([True]), TypeError, "time must contain numeric values"),
        (np.array(["0"]), TypeError, "time must contain numeric values"),
        (np.array([0, "1"], dtype=object), TypeError, "time must contain numeric values"),
        (np.array([np.nan]), ValueError, "time must contain only finite values"),
        (np.array([np.inf]), ValueError, "time must contain only finite values"),
        (np.array([-np.inf]), ValueError, "time must contain only finite values"),
    ],
)
def test_event_sequence_refuses_malformed_time_vectors(times, error, message):
    with pytest.raises(error, match=message):
        EventSequence(id=["p1"] * times.size, time=times, label=["A"] * times.size)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("id", [1]),
        ("id", [None]),
        ("id", [True]),
        ("id", "p1"),
        ("label", [1]),
        ("label", [None]),
        ("label", [True]),
        ("label", "A"),
    ],
)
def test_event_sequence_refuses_malformed_id_and_label_vectors(field, value):
    kwargs = {"id": ["p1"], "time": np.array([0.0]), "label": ["A"]}
    kwargs[field] = value

    with pytest.raises(TypeError, match=rf"{field} must be a collection of strings"):
        EventSequence(**kwargs)


def test_event_sequence_normalizes_justified_vector_inputs_with_defensive_copies():
    ids = ["p1", "p1"]
    times = [0, 0.5]
    labels = ["A", "B"]

    events = EventSequence(id=ids, time=times, label=labels)  # type: ignore[arg-type]
    ids[0] = "changed"
    times[0] = 99
    labels[0] = "changed"

    assert events.id == ["p1", "p1"]
    np.testing.assert_array_equal(events.time, np.array([0.0, 0.5]))
    assert events.label == ["A", "B"]


def test_empty_event_population_is_consistent_across_every_convenience_helper():
    events = EventSequence(id=[], time=np.array([], dtype=float), label=[])

    assert frequent_events(events) == []
    assert frequent_event_pairs(events) == []
    assert frequent_event_pairs_timegap(events, max_gap=1) == []
    labels, matrix = events_transition_matrix(events, normalize=True)
    degree_labels, indegree, outdegree = events_transition_degrees(events, normalize=True)
    assert labels == degree_labels == []
    assert matrix.shape == (0, 0)
    assert indegree.shape == outdegree.shape == (0,)

    automatic = events_to_state_sequences(events)
    supplied = events_to_state_sequences(events, time_grid=np.array([0.0, 1.0]))
    assert automatic.data.shape == (0, 0)
    assert supplied.data.shape == (0, 2)
    assert automatic.alphabet == supplied.alphabet == []
    assert automatic.ids == supplied.ids == []


def test_unsorted_fractional_times_preserve_stable_input_order_on_ties():
    events = EventSequence(
        id=["p1", "p1", "p1"],
        time=np.array([1.5, 0.5, 0.5]),
        label=["C", "B", "A"],
    )

    pairs = {pair: count for pair, count, _ in frequent_event_pairs(events)}
    assert pairs == {("A", "C"): 1, ("B", "A"): 1}

    labels, matrix = events_transition_matrix(events)
    index = {label: position for position, label in enumerate(labels)}
    assert matrix[index["B"], index["A"]] == 1.0
    assert matrix[index["A"], index["C"]] == 1.0

    states = events_to_state_sequences(events, time_grid=np.array([0.5, 1.0, 1.5]))
    assert states.to_strings() == ["A-A-C"]


def test_event_helpers_share_one_population_and_adjacent_pair_denominator():
    events = EventSequence(
        id=["p1", "p2", "p1", "p1", "p2"],
        time=np.array([2.0, 1.0, 0.0, 1.0, 2.0]),
        label=["C", "A", "A", "B", "B"],
    )

    frequencies = frequent_events(events)
    pairs = frequent_event_pairs(events)
    labels, matrix = events_transition_matrix(events)
    degree_labels, indegree, outdegree = events_transition_degrees(events)
    states = events_to_state_sequences(events, time_grid=np.array([0.0, 1.0, 2.0]))

    assert sum(count for _, count, _ in frequencies) == len(events.label) == 5
    assert sum(count for _, count, _ in pairs) == 3
    assert matrix.sum() == indegree.sum() == outdegree.sum() == 3.0
    assert labels == degree_labels == sorted(set(events.label))
    assert states.ids == ["p1", "p2"]


def _weighted_events() -> EventSequence:
    return EventSequence(
        id=["p1", "p1", "p2", "p2", "p3"],
        time=np.array([0.0, 1.0, 0.0, 1.0, 0.0]),
        label=["A", "B", "B", "A", "C"],
    )


def test_transition_weights_preserve_defaults_zero_fractional_and_extra_ids():
    events = _weighted_events()

    labels, matrix = events_transition_matrix(
        events,
        id_weights={"p1": 0.5, "p2": 0.0, "unused": 7.0},
    )
    index = {label: position for position, label in enumerate(labels)}
    assert matrix[index["A"], index["B"]] == 0.5
    assert matrix[index["B"], index["A"]] == 0.0
    assert matrix.sum() == 0.5

    _, defaulted = events_transition_matrix(events, id_weights={"p1": 0.5})
    assert defaulted.sum() == 1.5

    _, normalized = events_transition_matrix(
        events,
        normalize=True,
        id_weights={"p1": 0.5, "p2": 0.0},
    )
    assert normalized[index["A"], index["B"]] == 1.0
    assert np.all(normalized[index["B"]] == 0.0)

    _, indegree, outdegree = events_transition_degrees(
        events,
        id_weights={"p1": 0.5, "p2": 0.0},
    )
    assert indegree.sum() == outdegree.sum() == 0.5


def test_all_zero_effective_transition_weights_define_finite_zero_rows():
    events = _weighted_events()

    _, raw = events_transition_matrix(events, id_weights={"p1": 0.0, "p2": 0.0})
    _, normalized = events_transition_matrix(
        events,
        normalize=True,
        id_weights={"p1": 0.0, "p2": 0.0},
    )

    np.testing.assert_array_equal(raw, np.zeros_like(raw))
    np.testing.assert_array_equal(normalized, np.zeros_like(normalized))
    assert np.all(np.isfinite(normalized))


@pytest.mark.parametrize(
    ("weight", "error", "message"),
    [
        (-1.0, ValueError, "non-negative"),
        (np.nan, ValueError, "finite"),
        (np.inf, ValueError, "finite"),
        (-np.inf, ValueError, "finite"),
        ("1", TypeError, "real number"),
        (True, TypeError, "real number"),
        ([1.0], TypeError, "real number"),
        (np.array(1.0), TypeError, "real number"),
    ],
)
@pytest.mark.parametrize("consumer", [events_transition_matrix, events_transition_degrees])
def test_transition_consumers_reject_every_invalid_scalar_weight(
    weight,
    error,
    message,
    consumer,
):
    with pytest.raises(error, match=message):
        consumer(_weighted_events(), id_weights={"p1": weight})


def test_transition_matrix_validates_unused_weights_and_mapping_shape_eagerly():
    events = EventSequence(id=["p1"], time=np.array([0.0]), label=["A"])

    with pytest.raises(ValueError, match="finite"):
        events_transition_matrix(events, id_weights={"unused": np.nan})
    with pytest.raises(TypeError, match="id_weights must be a mapping"):
        events_transition_matrix(events, id_weights=[("p1", 1.0)])  # type: ignore[arg-type]


def test_the_validated_time_vector_cannot_be_mutated_in_place():
    """In-place writes must not bypass the validation reassignment enforces.

    ``EventSequence.__setattr__`` re-validates the whole (id, time, label)
    triple on assignment, so ``e.time = [0.0, nan]`` is refused with "time must
    contain only finite values". The same NaN written *through* the array --
    ``e.time[1] = nan`` -- was accepted, leaving the object holding a value its
    own validator rejects (cycle-10 N-19). The normalized vector is already a
    fresh copy, so freezing it cannot affect the caller's array.
    """
    import numpy as np

    import traminex as tx

    sequence = tx.EventSequence(id=["a", "a"], time=np.array([0.0, 1.0]), label=["x", "y"])

    with pytest.raises(ValueError, match="read-only|assignment destination"):
        sequence.time[1] = float("nan")

    assert np.all(np.isfinite(sequence.time))
    # Reassignment remains the supported way to change it, and still validates.
    sequence.time = np.array([0.0, 2.0])
    assert sequence.time.tolist() == [0.0, 2.0]
    with pytest.raises(ValueError, match="finite"):
        sequence.time = np.array([0.0, float("nan")])


def test_freezing_time_does_not_freeze_the_caller_s_array():
    """The caller keeps a writable array; only the stored copy is frozen."""
    import numpy as np

    import traminex as tx

    supplied = np.array([0.0, 1.0])
    tx.EventSequence(id=["a", "a"], time=supplied, label=["x", "y"])
    supplied[0] = 5.0  # must not raise
    assert supplied.tolist() == [5.0, 1.0]
