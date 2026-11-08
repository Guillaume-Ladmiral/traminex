"""dissmergegroups stops when no merge is found, where R loops (TX-AUD-20261011-016).

R's dissmergegroups records a pair only when its loss is below the reference
quality, that is, when the merged quality is positive. When no pair qualifies and
the reference quality is 0 or below, the loop takes the merge branch with nothing
to merge, printing "Merging groups and" without end (dissmergegroups.R:127-168).
traminex stops and returns the grouping unchanged, a recorded correction. On a
grouping that does merge, including one with a negative reference quality and a
merge of positive quality (batch-4 review M-1), traminex gives R's result.
Fixture: ``scripts/r_goldens/make_goldens_dissmergegroups_negative.R`` (R run under
setTimeLimit).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

import traminex as tx

HERE = Path(__file__).resolve().parent / "r_goldens"
GOLDEN = pd.read_csv(
    HERE / "dissmergegroups_negative.csv", dtype=str, keep_default_na=False
).set_index("case")
DISS = np.loadtxt(HERE / "dissmergegroups_negative_diss.csv", delimiter=",")


def _ints(text: str) -> np.ndarray:
    return np.array([int(value) for value in text.split("|")])


def test_r_does_not_finish_on_a_negative_quality_grouping() -> None:
    row = GOLDEN.loc["negative"]
    assert float(row["reference_asw"]) < 0
    assert row["status"] == "error: reached elapsed time limit"
    assert row["empty_merge_seen"] == "TRUE"


def test_traminex_stops_and_keeps_the_grouping() -> None:
    group = _ints(GOLDEN.loc["negative", "input"])
    result = tx.dissmergegroups(DISS, group)
    np.testing.assert_array_equal(np.asarray(result["group"]), group)


def test_a_merging_grouping_gives_r_result() -> None:
    row = GOLDEN.loc["control"]
    assert row["status"] == "returned" and row["empty_merge_seen"] == "FALSE"
    result = tx.dissmergegroups(DISS, _ints(row["input"]))
    np.testing.assert_array_equal(np.asarray(result["group"]), _ints(row["groups"]))


def test_a_negative_reference_with_a_positive_merge_merges_as_in_r() -> None:
    row = GOLDEN.loc["negative_merging"]
    assert float(row["reference_asw"]) < 0
    assert row["status"] == "returned" and row["empty_merge_seen"] == "FALSE"
    result = tx.dissmergegroups(DISS, _ints(row["input"]))
    np.testing.assert_array_equal(np.asarray(result["group"]), _ints(row["groups"]))
