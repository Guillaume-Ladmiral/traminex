"""Group labels read as R's as.character writes them (TX-AUD-20261010-026, -049).

seqdiff names its discrepancy columns after factor(group) levels, and
TraMineRextras' group.p and rowmode label with as.character. traminex used
Python's str(): False/True, 100000, 1000000000000000, 0.30000000000000004, and
kept a categorical's unused level as a column. A NaN group value is missing
here, where R's factor() keeps NaN as a level of its own; the expected values
are R's for the same groups with NA in place of NaN. Fixture:
``scripts/r_goldens/make_goldens_group_label_rendering.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "group_label_rendering.csv",
    dtype=str,
    keep_default_na=False,
)
STRINGS = ["AABB", "ABBB", "AAAB", "BBBA", "ABAB", "BAAA", "AABA", "BBAB"]
SEQS = tx.seqdef(np.array([list(s) for s in STRINGS], dtype=object), alphabet=["A", "B"])
GROUPS = {
    "logical": np.array([True, False] * 4),
    "doubles": np.array([100000, 2, 0.1, 1e15] * 2),
    "float_sum": np.array([0.1 + 0.2, 2.0] * 4),
    "unused_level": pd.Categorical(["b", "a"] * 4, categories=["a", "b", "zzz_unused"]),
    "edge_numbers": np.array([-0.0, 1e-20, 123456, 2] * 2),
}
MODES = {
    "big": [1e5, 1e5, 2.0],
    "huge": [1e15, 1e15, 3.0],
    "tiny": [1e-20, 1e-20, 2.0],
    "plain": [123456.0, 123456.0, 1.0],
    "truth": [True, True, False],
    "nan": [np.nan, np.nan, 1.0],
}
NAN_GROUP = np.array([np.nan, 1, 2, np.nan, 1, 2, np.nan, 1])


def _r(function: str, case: str) -> list[str]:
    rows = GOLDEN[(GOLDEN["function_name"] == function) & (GOLDEN["case"] == case)]
    return rows.sort_values("position", key=lambda c: c.astype(int))["label"].tolist()


@pytest.mark.parametrize("case", sorted(GROUPS))
def test_seqdiff_columns_are_r_factor_labels(case: str) -> None:
    columns = list(tx.seqdiff(SEQS, GROUPS[case])["discrepancy"].columns)
    assert columns == _r("seqdiff", case)


@pytest.mark.parametrize("case", ["doubles", "edge_numbers", "logical"])
def test_group_p_levels_are_r_labels(case: str) -> None:
    assert list(tx.group_p(GROUPS[case]).categories) == _r("group.p", case)


@pytest.mark.parametrize("case", sorted(MODES))
def test_rowmode_is_the_r_label(case: str) -> None:
    assert tx.rowmode(np.array(MODES[case], dtype=object)) == _r("rowmode", case)[0]


def test_a_nan_group_is_missing_where_r_keeps_a_nan_level() -> None:
    assert "NaN" in _r("seqdiff", "nan_level")
    assert any(label.startswith("NaN (") for label in _r("group.p", "nan_level"))
    columns = list(tx.seqdiff(SEQS, NAN_GROUP)["discrepancy"].columns)
    assert columns == _r("seqdiff", "nan_as_na")
    assert list(tx.group_p(NAN_GROUP).categories) == _r("group.p", "nan_as_na")
