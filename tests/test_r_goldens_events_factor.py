"""R-golden tests: seqecreate honors a factor / pandas Categorical event column's
LEVEL ORDER as the event dictionary (R keeps factor levels; it does not sort them).

This is the RP5b precursor: the event-mining goldens are generated on actcal.tse,
whose event column is a factor with non-alphabetical levels, and the factor codes
drive the comma-group ordering of simultaneous events. A plain (string) event
column still falls back to sorted-unique (R factor() of a character vector).

Generator: scripts/r_goldens/make_goldens_events_factor.R.
"""

from __future__ import annotations

import pandas as pd
from r_goldens_utils import load_golden_lines

from traminex import seqecreate


def _factor_tse() -> pd.DataFrame:
    # Mirrors scripts/r_goldens/make_goldens_events_factor.R: levels Gamma<Beta<Alpha.
    return pd.DataFrame(
        {
            "id": ["p1", "p1", "p1", "p2", "p2", "p3", "p3"],
            "timestamp": [0, 2, 2, 1, 3, 0, 0],
            "event": pd.Categorical(
                ["Beta", "Gamma", "Alpha", "Beta", "Alpha", "Gamma", "Alpha"],
                categories=["Gamma", "Beta", "Alpha"],
            ),
        }
    )


def test_seqecreate_factor_levels_drive_dictionary():
    e = seqecreate(_factor_tse())
    # Dictionary = the factor levels in declared order (NOT sorted alphabetically).
    assert list(e.dictionary) == load_golden_lines("eseq_factor_levels")
    assert list(e.dictionary) == ["Gamma", "Beta", "Alpha"]


def test_seqecreate_factor_simultaneous_sorted_by_factor_code():
    e = seqecreate(_factor_tse())
    # Same-timestamp events comma-group by FACTOR CODE (Gamma=1, Alpha=3) -> (Gamma,Alpha),
    # not alphabetically (which would give (Alpha,Gamma)).
    assert e.as_character() == load_golden_lines("eseq_factor_aschar")
    assert e.sequences[0].as_character() == "(Beta)-2-(Gamma,Alpha)"
    assert e.sequences[2].as_character() == "(Gamma,Alpha)"


def test_seqecreate_string_event_still_sorted_unique():
    # A plain (non-categorical) event column keeps the R factor()-of-character behaviour:
    # dictionary = sorted unique.
    tse = pd.DataFrame(
        {
            "id": ["i1", "i1", "i2"],
            "timestamp": [0, 1, 0],
            "event": ["C", "A", "B"],
        }
    )
    e = seqecreate(tse)
    assert list(e.dictionary) == ["A", "B", "C"]
