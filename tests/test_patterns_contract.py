"""The documented contract of the traminex-only pattern miners (W-03).

``traminex.patterns.frequent_subsequences`` and
``traminex.events._helpers.frequent_ngrams`` have **no TraMineR counterpart**,
so they can never be closed by comparison with R. Authority tier 3 applies
instead: their own documented contract, which lives in their docstrings.

Every expectation below is worked out **by hand** from that contract on inputs
small enough to check by inspection -- never copied from what the code
currently prints. Each case shows the derivation in a comment so a reader can
redo the arithmetic without running anything.

The pre-existing ``test_patterns.py``, ``test_patterns_support.py`` and
``test_patterns_noncontig.py`` assert only that the output is non-empty or that
one list is no shorter than another. They are kept, but they pin no value; this
module is what holds the contract.
"""

from __future__ import annotations

from math import comb

import numpy as np
import pandas as pd
import pytest

from traminex import seqdef
from traminex.events._helpers import frequent_ngrams
from traminex.patterns import frequent_subsequences


# --------------------------------------------------------------------------- #
# The worked corpus: two sequences of length three over the alphabet {a, b}.
#
#     sequence 1:  a b a
#     sequence 2:  a b b
#
# nseq = 2, seqlen = 3.
# --------------------------------------------------------------------------- #
def _corpus():
    return seqdef(pd.DataFrame([["a", "b", "a"], ["a", "b", "b"]]))


def test_the_worked_corpus_is_what_the_derivations_assume():
    x = _corpus()
    assert list(x.alphabet) == ["a", "b"]
    assert (x.nseq, x.seqlen) == (2, 3)


# --------------------------------------------------------------------------- #
# 1. contiguous=True -- one occurrence per sliding window.
# --------------------------------------------------------------------------- #
def test_contiguous_occurrence_counts_and_denominator():
    """By hand.

    Windows of sequence 1 (a b a): k=1 -> a, b, a; k=2 -> (a,b), (b,a).
    Windows of sequence 2 (a b b): k=1 -> a, b, b; k=2 -> (a,b), (b,b).

    Occurrence counts: a=3, b=3, (a,b)=2, (b,a)=1, (b,b)=1.
    Denominator = nseq * sum(seqlen - k + 1 for k in 1..2)
                = 2 * ((3-1+1) + (3-2+1)) = 2 * (3 + 2) = 10.
    min_support=2 keeps a, b and (a,b).
    Order: counts 3, 3, 2; the two 3s tie and break on the tuple, ("a",)<("b",).
    """
    got = frequent_subsequences(_corpus(), max_len=2, contiguous=True, min_support=2)

    assert got == [
        (("a",), 3, 3 / 10),
        (("b",), 3, 3 / 10),
        (("a", "b"), 2, 2 / 10),
    ]


def test_contiguous_min_support_is_an_inclusive_floor():
    """(b,a) and (b,b) occur exactly once, so min_support=1 must keep them and
    min_support=2 must drop them."""
    kept = frequent_subsequences(_corpus(), max_len=2, contiguous=True, min_support=1)
    dropped = frequent_subsequences(_corpus(), max_len=2, contiguous=True, min_support=2)

    assert (("b", "a"), 1, 1 / 10) in kept
    assert (("b", "b"), 1, 1 / 10) in kept
    assert [item for item in dropped if item[1] < 2] == []
    # exactly the two singleton patterns were removed
    assert len(kept) - len(dropped) == 2


def test_contiguous_support_fractions_sum_to_one_when_nothing_is_filtered():
    """Every window is counted exactly once, so with min_support=1 the
    fractions partition the 10 slots."""
    got = frequent_subsequences(_corpus(), max_len=2, contiguous=True, min_support=1)

    assert sum(count for _, count, _ in got) == 10
    assert sum(fraction for _, _, fraction in got) == pytest.approx(1.0)


# --------------------------------------------------------------------------- #
# 2. support_type="sequences" -- the denominator is the number of sequences.
# --------------------------------------------------------------------------- #
def test_sequences_support_counts_sequences_not_occurrences():
    """By hand.

    Sequence 1 contains a, b, (a,b), (b,a). Sequence 2 contains a, b, (a,b),
    (b,b). So the per-sequence counts are a=2, b=2, (a,b)=2, (b,a)=1, (b,b)=1.
    Denominator = nseq = 2, so the kept patterns all have fraction 1.0.
    Order: all three tie at 2, so the tuple decides:
    ("a",) < ("a","b") < ("b",)  -- a prefix sorts before its extension.
    """
    got = frequent_subsequences(
        _corpus(), max_len=2, contiguous=True, min_support=2, support_type="sequences"
    )

    assert got == [
        (("a",), 2, 1.0),
        (("a", "b"), 2, 1.0),
        (("b",), 2, 1.0),
    ]


def test_sequences_support_fraction_never_exceeds_one():
    got = frequent_subsequences(
        _corpus(), max_len=3, contiguous=False, min_support=1, support_type="sequences"
    )

    assert all(0.0 < fraction <= 1.0 for _, _, fraction in got)
    assert all(count <= 2 for _, count, _ in got)


def test_the_two_support_types_disagree_where_a_pattern_repeats_in_one_sequence():
    """Guards the two cases above against passing for the same reason.

    "a" occurs twice inside sequence 1, so its occurrence count (3) exceeds its
    sequence count (2).
    """
    occurrences = dict(
        (pattern, count)
        for pattern, count, _ in frequent_subsequences(
            _corpus(), max_len=1, contiguous=True, min_support=1
        )
    )
    sequences = dict(
        (pattern, count)
        for pattern, count, _ in frequent_subsequences(
            _corpus(), max_len=1, contiguous=True, min_support=1, support_type="sequences"
        )
    )

    assert occurrences[("a",)] == 3
    assert sequences[("a",)] == 2


# --------------------------------------------------------------------------- #
# 3. contiguous=False -- one occurrence per combination of positions.
# --------------------------------------------------------------------------- #
def test_non_contiguous_counts_position_sets_and_uses_a_binomial_denominator():
    """By hand.

    Sequence 1 (a b a), positions 0,1,2:
      k=1 -> a, b, a
      k=2 -> (0,1)=(a,b), (0,2)=(a,a), (1,2)=(b,a)
    Sequence 2 (a b b):
      k=1 -> a, b, b
      k=2 -> (0,1)=(a,b), (0,2)=(a,b), (1,2)=(b,b)

    Counts: a=3, b=3, (a,b)=3, (a,a)=1, (b,a)=1, (b,b)=1.
    Note (a,b)=3 while the corpus holds only 2 sequences: position sets are
    counted, so sequence 2 contributes it twice.
    Denominator = nseq * sum(C(3,k) for k in 1..2) = 2 * (3 + 3) = 12.
    Order: all three kept patterns tie at 3, so ("a",) < ("a","b") < ("b",).
    """
    got = frequent_subsequences(_corpus(), max_len=2, contiguous=False, min_support=2)

    assert got == [
        (("a",), 3, 3 / 12),
        (("a", "b"), 3, 3 / 12),
        (("b",), 3, 3 / 12),
    ]


def test_non_contiguous_can_count_a_pattern_more_often_than_contiguous_does():
    """(a,b) appears in 2 windows but in 3 position sets."""
    contiguous = dict(
        (pattern, count)
        for pattern, count, _ in frequent_subsequences(
            _corpus(), max_len=2, contiguous=True, min_support=1
        )
    )
    combinations = dict(
        (pattern, count)
        for pattern, count, _ in frequent_subsequences(
            _corpus(), max_len=2, contiguous=False, min_support=1
        )
    )

    assert contiguous[("a", "b")] == 2
    assert combinations[("a", "b")] == 3
    # (a,a) is not contiguous anywhere but is a valid position set
    assert ("a", "a") not in contiguous
    assert combinations[("a", "a")] == 1


def test_non_contiguous_denominator_matches_the_documented_formula():
    x = _corpus()
    for max_len in (1, 2, 3, 5):
        got = frequent_subsequences(x, max_len=max_len, contiguous=False, min_support=1)
        expected_denominator = x.nseq * sum(
            comb(x.seqlen, k) for k in range(1, min(max_len, x.seqlen) + 1)
        )
        total = sum(count for _, count, _ in got)
        assert total == expected_denominator, max_len
        assert sum(fraction for _, _, fraction in got) == pytest.approx(1.0)


# --------------------------------------------------------------------------- #
# 4. max_len beyond the sequence length contributes nothing, either side.
# --------------------------------------------------------------------------- #
def test_lengths_longer_than_the_sequences_change_neither_numerator_nor_denominator():
    """By hand.

    With max_len=5 and seqlen=3 the contiguous denominator is
    2 * (3 + 2 + 1 + 0 + 0) = 12, and no pattern longer than 3 can exist.
    """
    x = _corpus()
    got = frequent_subsequences(x, max_len=5, contiguous=True, min_support=1)

    assert max(len(pattern) for pattern, _, _ in got) == 3
    assert sum(count for _, count, _ in got) == 12
    assert (("a",), 3, 3 / 12) in got
    # max_len=3 and max_len=5 must agree exactly
    assert got == frequent_subsequences(x, max_len=3, contiguous=True, min_support=1)


# --------------------------------------------------------------------------- #
# 5. The ragged case -- padding is an ordinary state and inflates the
#    denominator. Documented, not silently accepted.
# --------------------------------------------------------------------------- #
def test_ragged_sequences_count_the_void_marker_as_a_state():
    """By hand.

        sequence 1:  a b a b
        sequence 2:  a b % %      ("%" is the void padding)

    nseq = 2, seqlen = 4 -- the full matrix width, for both sequences.

    Contiguous windows, k=1: sequence 1 gives a,b,a,b; sequence 2 gives
    a,b,%,%. So a=3, b=3, %=2.
    k=2: sequence 1 gives (a,b),(b,a),(a,b); sequence 2 gives (a,b),(b,%),(%,%).
    So (a,b)=3, (b,a)=1, (b,%)=1, (%,%)=1.

    Denominator = 2 * ((4-1+1) + (4-2+1)) = 2 * (4 + 3) = 14 -- it assumes both
    sequences fill all four positions, which sequence 2 does not.
    """
    x = seqdef(pd.DataFrame([["a", "b", "a", "b"], ["a", "b", None, None]]))
    assert list(x.alphabet) == ["a", "b", "%"]
    assert (x.nseq, x.seqlen) == (2, 4)

    got = frequent_subsequences(x, max_len=2, contiguous=True, min_support=1)

    assert got == [
        (("a",), 3, 3 / 14),
        (("a", "b"), 3, 3 / 14),
        (("b",), 3, 3 / 14),
        (("%",), 2, 2 / 14),
        (("%", "%"), 1, 1 / 14),
        (("b", "%"), 1, 1 / 14),
        (("b", "a"), 1, 1 / 14),
    ]


def test_the_ragged_denominator_exceeds_the_observed_positions():
    """The documented over-count, asserted rather than assumed.

    Sequence 2 holds 2 real states, not 4, so an occurrence-based denominator
    built from the matrix width is larger than the data supports.
    """
    x = seqdef(pd.DataFrame([["a", "b", "a", "b"], ["a", "b", None, None]]))
    got = frequent_subsequences(x, max_len=1, contiguous=True, min_support=1)

    padded_denominator = x.nseq * x.seqlen
    real_positions = 4 + 2

    assert sum(count for _, count, _ in got) == padded_denominator == 8
    assert real_positions < padded_denominator
    # and the padding itself is reported as though it were a state
    assert (("%",), 2, 2 / 8) in got


# --------------------------------------------------------------------------- #
# 6. Degenerate inputs.
# --------------------------------------------------------------------------- #
def test_an_unknown_support_type_is_refused():
    with pytest.raises(ValueError, match="support_type"):
        frequent_subsequences(_corpus(), support_type="proportion")


def test_a_min_support_above_every_count_returns_an_empty_list():
    assert frequent_subsequences(_corpus(), max_len=2, min_support=99) == []


def test_a_single_sequence_of_one_state_is_handled():
    """By hand: one window, denominator 1 * max(1-1+1, 0) = 1."""
    x = seqdef(pd.DataFrame([["a"]]))

    assert frequent_subsequences(x, max_len=1, min_support=1) == [(("a",), 1, 1.0)]


# --------------------------------------------------------------------------- #
# 7. frequent_ngrams.
# --------------------------------------------------------------------------- #
def test_frequent_ngrams_counts_windows_over_the_position_denominator():
    """By hand.

    n=2 over (a b a) and (a b b): sequence 1 gives (a,b),(b,a); sequence 2
    gives (a,b),(b,b). Counts (a,b)=2, (b,a)=1, (b,b)=1.
    Denominator = nseq * (seqlen - n + 1) = 2 * (3 - 2 + 1) = 4.
    """
    got = frequent_ngrams(_corpus(), n=2, min_support=1)

    assert got == [
        (("a", "b"), 2, 2 / 4),
        (("b", "a"), 1, 1 / 4),
        (("b", "b"), 1, 1 / 4),
    ]


def test_frequent_ngrams_unigrams_match_the_contiguous_length_one_counts():
    """n=1: denominator 2 * 3 = 6, counts a=3, b=3."""
    assert frequent_ngrams(_corpus(), n=1, min_support=1) == [
        (("a",), 3, 3 / 6),
        (("b",), 3, 3 / 6),
    ]


@pytest.mark.parametrize("n", [0, -1, 4, 10])
def test_frequent_ngrams_returns_empty_outside_the_usable_range(n):
    """n <= 0 or n > seqlen yields [] rather than raising."""
    assert frequent_ngrams(_corpus(), n=n, min_support=1) == []


def test_frequent_ngrams_min_support_is_an_inclusive_floor():
    assert frequent_ngrams(_corpus(), n=2, min_support=2) == [(("a", "b"), 2, 2 / 4)]
    assert frequent_ngrams(_corpus(), n=2, min_support=3) == []


def test_frequent_ngrams_agrees_with_frequent_subsequences_at_the_same_length():
    """Both count contiguous windows, so the counts must coincide; only the
    denominators differ, because frequent_subsequences sums over 1..max_len."""
    ngrams = {
        pattern: count for pattern, count, _ in frequent_ngrams(_corpus(), n=2, min_support=1)
    }
    subsequences = {
        pattern: count
        for pattern, count, _ in frequent_subsequences(
            _corpus(), max_len=2, contiguous=True, min_support=1
        )
        if len(pattern) == 2
    }

    assert ngrams == subsequences


def test_frequent_ngrams_on_a_uniform_corpus_has_one_gram_at_fraction_one():
    """By hand: (a,a) fills every one of 2 * (3-2+1) = 4 window positions."""
    x = seqdef(pd.DataFrame([["a", "a", "a"], ["a", "a", "a"]]))

    assert frequent_ngrams(x, n=2, min_support=1) == [(("a", "a"), 4, 1.0)]


def test_numpy_int_matrices_are_accepted_and_labelled_by_the_alphabet():
    """seqdef over integer codes names the states "0"/"1"; the contract is the
    same, so this only pins the labelling."""
    x = seqdef(np.array([[0, 1, 0], [0, 1, 1]]))

    assert frequent_ngrams(x, n=2, min_support=2) == [(("0", "1"), 2, 0.5)]
