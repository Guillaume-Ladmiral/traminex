"""JSD_TRANS between sequences without a transition is 0, and the documents say so (TX-AUD-20261011-046).

A sequence with fewer than two positions has no transition, so it has no transition distribution; the kernels
return 0 for any pair involving one. No document stated it.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

import traminex as tx
from traminex.distances import python_kernels

REPO = Path(__file__).resolve().parents[1]


def test_four_length_1_sequences_give_0() -> None:
    seqmat = np.array([[0], [1], [2], [1]], dtype=np.int32)
    np.testing.assert_array_equal(
        np.asarray(tx.seqdist(seqmat, method="JSD_TRANS"), dtype=float), np.zeros((4, 4))
    )


def test_a_pair_involving_a_sequence_without_a_transition_gives_0() -> None:
    longer = np.array([0, 1, 1, 0])
    assert python_kernels.jsd_transition_distance(np.array([0]), longer, 2) == 0.0
    assert python_kernels.jsd_transition_distance(longer, np.array([1]), 2) == 0.0


def test_the_kernel_docstring_and_distances_md_state_it() -> None:
    doc = " ".join((python_kernels.jsd_transition_distance.__doc__ or "").split())
    assert "fewer than two positions has no transition; the distance is then 0" in doc
    text = " ".join((REPO / "docs" / "DISTANCES.md").read_text(encoding="utf-8").split())
    assert "A sequence with fewer than two positions has no transition" in text
