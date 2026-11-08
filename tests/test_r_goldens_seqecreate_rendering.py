"""seqecreate renders event values and matches end_event as R does.

Repairs audit findings TX-AUD-20261010-002 (a numeric ``end_event`` was
ignored) and TX-AUD-20261010-016 (numeric and logical event labels were
rendered with Python ``str``). R builds ``factor(event)``, so labels are R's
``as.character`` of the values, and finds ``end.event`` with ``==`` against
those labels (seqecreate.R:63-73). Fixtures:
``scripts/r_goldens/make_goldens_seqecreate_rendering.R``.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
import pytest

from traminex.events import seqecreate, seqelength

GOLDEN = Path(__file__).resolve().parent / "r_goldens" / "seqecreate_rendering.csv"


def _load() -> dict[str, dict[str, list[str]]]:
    cases: dict[str, dict[str, list[tuple[int, str]]]] = defaultdict(lambda: defaultdict(list))
    with GOLDEN.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            cases[row["case"]][row["kind"]].append((int(row["position"]), row["value"]))
    return {
        case: {kind: [value for _, value in sorted(items)] for kind, items in kinds.items()}
        for case, kinds in cases.items()
    }


R = _load()

# The arguments of each TraMineR::seqecreate call in the generator. R numeric
# vectors are doubles (Python float); the *_integer cases are R integers.
CALLS: dict[str, dict[str, object]] = {
    "numeric_end": dict(
        id=[1, 1, 1, 2, 2],
        timestamp=[0, 2, 5, 1, 3],
        event=[1.0, 2.0, 99.0, 2.0, 99.0],
        end_event=99.0,
    ),
    "labels_1e5": dict(
        id=[1, 1, 1, 2, 2],
        timestamp=[0, 2, 5, 1, 3],
        event=[1.0, 100000.0, 7.0, 100000.0, 7.0],
        end_event=7.0,
    ),
    "logical": dict(id=[1, 1, 2], timestamp=[0, 2, 1], event=[True, False, True]),
    "end_1e5": dict(
        id=[1, 1, 2, 2], timestamp=[0, 2, 1, 4], event=[1.0, 100000.0, 2.0, 100000.0], end_event=1e5
    ),
    "absent_end": dict(id=[1, 1, 2], timestamp=[0, 2, 1], event=[1.0, 2.0, 2.0], end_event=42.0),
    "logical_end": dict(
        id=[1, 1, 2, 2], timestamp=[0, 3, 1, 2], event=[True, False, True, False], end_event=False
    ),
    "numeric_end_integer": dict(
        id=[1, 1, 1, 2, 2], timestamp=[0, 2, 5, 1, 3], event=[1, 2, 99, 2, 99], end_event=99
    ),
    "labels_1e5_integer": dict(
        id=[1, 1, 1, 2, 2], timestamp=[0, 2, 5, 1, 3], event=[1, 100000, 7, 100000, 7], end_event=7
    ),
    "double_end": dict(
        id=[1, 1, 2, 2], timestamp=[0, 3, 1, 2], event=[0.5, 1.0, 0.5, 1.0], end_event=1.0
    ),
}


def test_every_generated_case_is_exercised():
    assert set(CALLS) == set(R)


@pytest.mark.parametrize("case", sorted(CALLS))
def test_event_dictionary_matches_r_levels(case: str):
    assert list(seqecreate(**CALLS[case]).dictionary) == R[case]["level"]


@pytest.mark.parametrize("case", sorted(CALLS))
def test_event_sequences_match_r(case: str):
    assert seqecreate(**CALLS[case]).as_character() == R[case]["sequence"]


@pytest.mark.parametrize("case", sorted(CALLS))
def test_observation_lengths_match_r(case: str):
    lengths = np.asarray(seqelength(seqecreate(**CALLS[case])), dtype=float)
    np.testing.assert_array_equal(lengths, [float(value) for value in R[case]["length"]])


def test_string_events_keep_their_dictionary():
    e = seqecreate(id=[1, 1, 2], timestamp=[0, 1, 0], event=["b", "a", "b"], end_event="a")
    assert list(e.dictionary) == ["a", "b"]
    assert e.as_character() == ["(b)-1", "(b)"]
