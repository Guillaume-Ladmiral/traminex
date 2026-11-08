"""R-golden tests for the faithful event-sequence object (eseq / seqelist).

Pins ``as.character(seqecreate(tse))`` -- the TraMineR event-label strings --
and ``levels(eseq)`` (the event dictionary) against live TraMineR 2.2-14, for a
TSE exercising a leading gap, simultaneous events (comma group, sorted by event
code), a repeated event and a singleton.

Generator: scripts/r_goldens/make_goldens_events_obj.R.
"""

from __future__ import annotations

import pandas as pd
from r_goldens_utils import load_golden_lines

from traminex import seqecreate


def _tse() -> pd.DataFrame:
    # Mirrors scripts/r_goldens/make_goldens_events_obj.R.
    return pd.DataFrame(
        {
            "id": ["i1", "i1", "i1", "i2", "i2", "i3", "i3", "i4", "i4", "i4", "i5"],
            "timestamp": [0, 2, 2, 1, 3, 0, 0, 0, 1, 2, 4],
            "event": ["A", "B", "C", "A", "B", "B", "A", "C", "C", "A", "B"],
        }
    )


def test_eseq_as_character_matches_r():
    e = seqecreate(_tse())
    expected = load_golden_lines("eseq_aschar")
    assert e.as_character() == expected
    # Per-sequence __str__ matches each label too.
    assert [str(s) for s in e] == expected
    assert str(e) == "\n".join(expected)


def test_eseq_levels_matches_r():
    e = seqecreate(_tse())
    assert list(e.dictionary) == load_golden_lines("eseq_levels")


def test_eseq_simultaneous_events_sorted_by_code():
    # i3 input order is B,A at the same timestamp; R sorts ties by event code -> (A,B).
    e = seqecreate(_tse())
    assert e.sequences[2].as_character() == "(A,B)"


def test_seqecreate_keeps_flat_views_for_approx_helpers():
    # The non-R helpers still see flat per-event input-order views.
    e = seqecreate(_tse())
    assert e.id[:3] == ["i1", "i1", "i1"]
    assert e.label[5:7] == ["B", "A"]  # i3 input order preserved (not code-sorted)
    assert len(e) == 5
