"""Event sequences compare by value (TX-AUD-20261011-059).

R defines no equality for event sequences (``Ops.seqelist`` is commented out, seqe.R:41). The generated dataclass
``__eq__`` of ``Eseq`` compared its NumPy fields inside a tuple, so two distinct sequences of more than one event raised
``ValueError`` (the truth value of an array is ambiguous). Equality is now by value over the compared fields: the id,
times, event codes, dictionary, weight and stored length, and the engine's gap chain, which the input order of
fractional times decides (batch-11 review M-6). An ``Eseq`` is mutable, so it stays unhashable.
"""

from __future__ import annotations

import itertools

import pytest

import traminex as tx
from traminex.events._object import _copy_eseq


def _first(timestamps, events, ids=(1, 1, 1)):
    return tx.seqecreate(id=list(ids), timestamp=list(timestamps), event=list(events)).sequences[0]


def test_equal_event_sequences_compare_equal() -> None:
    first, second = _first([0, 1, 3], "abc"), _first([0, 1, 3], "abc")
    assert first is not second
    assert first == second
    assert not first != second
    assert _copy_eseq(first) == first


@pytest.mark.parametrize(
    "timestamps,events,ids",
    [([0, 1, 4], "abc", (1, 1, 1)), ([0, 1, 3], "abd", (1, 1, 1)), ([0, 1, 3], "abc", (2, 2, 2))],
)
def test_a_different_time_event_or_id_makes_them_unequal(timestamps, events, ids) -> None:
    assert _first([0, 1, 3], "abc") != _first(timestamps, events, ids)


def test_a_different_weight_or_length_makes_them_unequal() -> None:
    first, second = _first([0, 1, 3], "abc"), _first([0, 1, 3], "abc")
    second.weight = 2.0
    assert first != second
    second.weight, second.obs_time = first.weight, 10.0
    assert first != second


def test_an_event_sequence_is_not_equal_to_another_type_and_is_unhashable() -> None:
    first = _first([0, 1, 3], "abc")
    assert first != first.as_character()
    with pytest.raises(TypeError):
        hash(first)


def _by_input_order(timestamps, events):
    """One Eseq per input order of the events, keyed by its printed label."""
    built = {}
    for order in itertools.permutations(range(len(events))):
        sequence = _first(
            [timestamps[i] for i in order], [events[i] for i in order], ids=[1] * len(events)
        )
        built.setdefault(sequence.as_character(), sequence)
    return built


def test_input_orders_that_print_differently_are_unequal() -> None:
    # With fractional times the input order decides how the gaps round, so the printed gaps and the counts; equality
    # compares the gap chain the engine reads (batch-11 review M-6).
    built = _by_input_order([1.3, 0.8, 0.1, 1.3], "abcd")
    assert len(built) > 1
    first, second = list(built.values())[:2]
    assert first != second


def test_input_orders_of_whole_number_times_stay_equal() -> None:
    built = {}
    for order in itertools.permutations(range(3)):
        times, events = [[0, 2, 2][i] for i in order], ["a", "b", "c"]
        built[order] = _first(times, [events[i] for i in order])
    assert len({s.as_character() for s in built.values()}) == 1
    assert all(s == built[(0, 1, 2)] for s in built.values())
