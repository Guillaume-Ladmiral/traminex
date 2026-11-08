"""``seqRLCP`` agrees with ``seqdist(method="RLCP")`` on ragged sequences.

``seqRLCP`` is traminex-only and documents itself as equivalent to
``seqdist(method="RLCP")``, but it passed the raw codes -- void included -- to
the kernel, so a ragged object compared its padding as a state (2026-09-22
audit S1a-08). It now takes ``seqdist``'s raw distances (with the missing state
kept as a state, as before) and keeps its own documented ``norm="length"``:
the low-level ``l1 + l2`` denominator, ``2 * L`` for equal lengths.
"""

from __future__ import annotations

import numpy as np

import traminex as tx

RAGGED = tx.seqdef(
    np.array([list("ABCB"), ["A", "C", "B", None], list("BBCB")], dtype=object),
    alphabet=["A", "B", "C"],
)


def test_ragged_objects_delete_void_as_seqdist_does():
    np.testing.assert_array_equal(tx.seqRLCP(RAGGED), tx.seqdist(RAGGED, method="RLCP"))


def test_length_norm_uses_the_sum_of_the_two_lengths():
    raw = np.asarray(tx.seqdist(RAGGED, method="RLCP"), dtype=float)
    lengths = np.array([4.0, 3.0, 4.0])
    np.testing.assert_allclose(
        tx.seqRLCP(RAGGED, norm="length"), raw / np.add.outer(lengths, lengths), atol=1e-12
    )


def test_a_missing_state_is_still_compared_as_a_state():
    with_missing = tx.seqdef(
        np.array([list("ABCB"), ["A", None, "C", "B"]], dtype=object), alphabet=["A", "B", "C"]
    )
    np.testing.assert_array_equal(
        tx.seqRLCP(with_missing), tx.seqdist(with_missing, method="RLCP", with_missing=True)
    )


def test_plain_matrices_are_unchanged():
    X = np.array([[0, 1, 2, 1], [1, 1, 2, 1], [2, 2, 0, 1]])
    from traminex.distances import pairwise

    np.testing.assert_array_equal(tx.seqRLCP(X), pairwise(X, method="RLCP"))
