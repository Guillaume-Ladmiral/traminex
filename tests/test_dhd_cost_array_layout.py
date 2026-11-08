"""DHD's explicit cost array is time-first, and a TraMineR-layout array says so.

``seqdist(method="DHD", subm=array)`` takes the time-varying costs time-first,
``(T, K, K)``, while ``seqcost(..., time_varying=True)["sm"]`` returns R's
``(K, K, T)``. So R's own workflow -- passing ``seqcost``'s array to
``seqdist`` -- failed with a shape error that did not name the layout, and when
``T == K`` with a misleading "diagonals must be zero" (2026-09-22 audit
S1a-03). The time-first convention is kept (it is pinned by the DHD goldens and
``seqMD``); the refusal now names the layout and the one-line conversion.
"""

from __future__ import annotations

import numpy as np
import pytest

import traminex as tx

ROWS = [list("AABCA"), list("ABBCC"), list("CCABA"), list("BACCB")]


def _sequences(width: int = 5):
    return tx.seqdef(
        np.array([row[:width] for row in ROWS], dtype=object), alphabet=["A", "B", "C"]
    )


def test_seqcosts_r_layout_array_is_refused_with_the_conversion():
    sequences = _sequences()
    r_layout = np.asarray(tx.seqcost(sequences, "TRATE", time_varying=True)["sm"], dtype=float)
    assert r_layout.shape == (3, 3, 5)

    with pytest.raises(ValueError, match=r"np\.moveaxis\(sm, -1, 0\)"):
        tx.seqdist(sequences, method="DHD", subm=r_layout)

    converted = tx.seqdist(sequences, method="DHD", subm=np.moveaxis(r_layout, -1, 0))
    assert np.asarray(converted).shape == (4, 4)


def test_a_square_cube_in_r_layout_gets_the_layout_hint_not_a_diagonal_complaint():
    """With T == K both readings have the right shape; only the diagonals tell."""
    sequences = _sequences(width=3)
    r_layout = np.asarray(tx.seqcost(sequences, "TRATE", time_varying=True)["sm"], dtype=float)
    assert r_layout.shape == (3, 3, 3)
    assert np.any(np.diagonal(r_layout, axis1=1, axis2=2) != 0), "time-first reading is invalid"

    with pytest.raises(ValueError, match=r"\(K, K, T\) layout.*np\.moveaxis"):
        tx.seqdist(sequences, method="DHD", subm=r_layout)
