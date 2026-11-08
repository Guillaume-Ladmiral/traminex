"""``seqeapplysub`` names the missing window size for the CWIN count method.

The CWIN method counts the windows of ``window.size`` that contain a
subsequence, so it needs that size. ``seqefsub`` refuses CWIN without one
(``seqefsub.R:23-25``); ``seqeapplysub.R`` has no such check, and R aborts
the whole session (``std::bad_array_new_length``, exit status 134). traminex
raised ``ValueError: cannot convert float NaN to integer`` from inside the
counting engine (2026-09-22 audit S5-11). It now raises the same named
refusal as ``seqefsub``; the parity ledger records this where R crashes.
"""

from __future__ import annotations

import numpy as np
import pytest

import traminex as tx


def _subsequences():
    eseq = tx.seqecreate(id=[1, 1, 2, 2], timestamp=[0, 2, 0, 3], event=["a", "b", "a", "b"])
    return tx.seqefsub(eseq, min_support=1)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"method": "CWIN"},
        {"method": 3},
        {"constraint": tx.seqeconstraint(count_method="CWIN")},
        {"constraint": tx.seqeconstraint(count_method=3), "rules": True},
    ],
)
def test_cwin_without_a_window_size_is_refused_by_name(kwargs):
    with pytest.raises(ValueError, match=r"count method 3 or 'CWIN' requires a window\.size"):
        tx.seqeapplysub(_subsequences(), **kwargs)


def test_cwin_with_a_window_size_still_counts():
    counts = tx.seqeapplysub(
        _subsequences(), constraint=tx.seqeconstraint(window_size=2, count_method="CWIN")
    )
    assert counts.shape == (2, 3)
    assert np.all(np.isfinite(counts.to_numpy()))
