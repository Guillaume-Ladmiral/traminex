"""traminex-only helpers do not treat the void padding as a state (audit S3-19).

Void is the padding that gives unequal-length sequences a common width; it is
never a state anywhere in the package (``seqlength``, ``seqdss``, ``seqtrate``
and the rest drop it). Three traminex-only helpers, which have no R
counterpart and are judged against their own documented maths, counted it:

* ``entropy`` (global Shannon entropy of the state distribution) counted void
  cells as a state and ``log(K)`` included it, so data with a single observed
  state had an entropy of 0.64 instead of 0;
* ``seqtransn_count`` ("number of state changes per sequence") and so
  ``turbulence_proxy`` counted the step from the last state into the padding
  as a change, which ``seqtree_fit`` then used as a feature;
* ``_transition_count_matrix`` (behind ``transition_heatmap``, documented as
  equivalent to ``seqtrate`` when row-normalised) had transitions into and
  out of void.

Expected values follow from the definitions: the same cells without padding
give the same answer.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import traminex as tx
from traminex.stats_basic import _transition_count_matrix

N = None


def _sts(rows, alphabet):
    return tx.seqdef(np.array(rows, dtype=object), alphabet=alphabet)


def test_entropy_of_a_single_observed_state_is_zero_despite_padding():
    x = _sts([["A", "A", "A"], ["A", N, N]], ["A"])
    assert x.void_code is not None
    assert tx.entropy(x, normalized=False) == pytest.approx(0.0, abs=1e-9)
    assert tx.entropy(x) == 0.0


def test_padding_does_not_change_the_entropy():
    ragged = _sts([["A", "B", "B"], ["A", N, N]], ["A", "B"])
    # The same four state cells (A, B, B, A) without padding.
    rectangular = _sts([["A", "B"], ["B", "A"]], ["A", "B"])
    assert tx.entropy(ragged, normalized=False) == pytest.approx(math.log(2), abs=1e-9)
    assert tx.entropy(ragged, normalized=False) == tx.entropy(rectangular, normalized=False)
    assert tx.entropy(ragged) == tx.entropy(rectangular)


def test_entering_the_padding_is_not_a_state_change():
    y = _sts([["A", "B", "B"], ["A", N, N]], ["A", "B"])
    np.testing.assert_array_equal(tx.seqtransn_count(y, normalized=False), [1.0, 0.0])
    # Normalised by the object width minus one, as documented.
    np.testing.assert_allclose(tx.seqtransn_count(y), [0.5, 0.0])
    assert tx.turbulence_proxy(y) == pytest.approx(0.25)
    # The per-sequence counts agree with R's seqtransn on the same object.
    np.testing.assert_array_equal(tx.seqtransn_count(y, normalized=False), tx.seqtransn(y))


def test_transition_matrix_has_no_transition_into_or_out_of_the_padding():
    y = _sts([["A", "B", "B"], ["A", N, N]], ["A", "B"])
    void = y.void_code
    counts = _transition_count_matrix(y)
    assert counts.shape == (len(y.alphabet), len(y.alphabet))
    assert counts[void].sum() == 0 and counts[:, void].sum() == 0
    assert counts[0, 1] == 1 and counts[1, 1] == 1 and counts.sum() == 2
    rates = _transition_count_matrix(y, normalize=True)
    np.testing.assert_allclose(rates[:2, :2], [[0.0, 1.0], [0.0, 1.0]])


def test_exploratory_tree_features_leave_missing_and_void_out():
    """seqtree_fit's features are over the alphabet's states (TX-AUD-20261010-081).

    The features had a column for the missing and void markers and counted
    their cells: in ``B B % %`` the step into the padding was a transition
    ``tr_B_%`` and the padding a state of the entropy. A marker cell is not a
    state, so it has no column, joins no transition and is outside the entropy.
    """
    from traminex.seqtree_explore import _features_from_sequences

    x = _sts([["A", "A", "B", "B"], ["A", np.nan, "B", N], ["B", "B", N, N]], ["A", "B"])
    features, names = _features_from_sequences(x)
    assert names == [
        "dur_A",
        "dur_B",
        "tr_A_A",
        "tr_A_B",
        "tr_B_A",
        "tr_B_B",
        "entropy_norm",
        "transition_count_norm",
    ]
    expected = np.array(
        [
            # durations; transitions between two states; entropy of the states, over log(2)
            [2, 2, 1, 1, 0, 1, 1.0],
            [1, 1, 0, 0, 0, 0, 1.0],
            [0, 2, 0, 0, 0, 1, 0.0],
        ]
    )
    np.testing.assert_allclose(features[:, :7], expected, rtol=1e-9, atol=1e-9)
    # Changes between two adjacent state cells over seqlen - 1 (batch-12 review M-2): A A B B one,
    # A * B % and B B % % none.
    np.testing.assert_allclose(features[:, 7], [1 / 3, 0, 0], rtol=0, atol=1e-15)


def test_exploratory_tree_docs_say_what_they_are():
    """No 'analog' wording: an MDS embedding or sequence features fed to a scikit-learn tree (TX-AUD-20261010-081)."""
    import inspect

    from traminex import disstree_explore, seqtree_explore

    for text in (
        disstree_explore.__doc__,
        seqtree_explore.__doc__,
        inspect.getdoc(disstree_explore.disstree_fit),
        inspect.getdoc(seqtree_explore.seqtree_fit),
    ):
        words = " ".join(text.split()).lower()
        assert "analog" not in words and "-like" not in words and "approximat" not in words, words
        assert "scikit-learn" in words, words
