"""R-golden coverage for TraMineRextras 0.6.9 ``sorti`` and ``sortv``."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDENS_DIR

import traminex
from traminex import stats_utils
from traminex.sequences import STSList


def _cases() -> dict[str, STSList]:
    return {
        "declared_order": STSList(
            data=np.array(
                [
                    [1, 2, 0],
                    [2, 1, 0],
                    [1, 2, 2],
                    [2, 1, 2],
                    [1, 2, 0],
                ],
                dtype=int,
            ),
            alphabet=["C", "B", "A", "*", "%"],
            labels=["C", "B", "A", "*", "%"],
            nr_code=3,
            void_code=4,
        ),
        "ragged_void": STSList(
            data=np.array([[0, 1, 3], [0, 2, 3], [1, 0, 1], [3, 3, 3]], dtype=int),
            alphabet=["A", "B", "*", "%"],
            labels=["A", "B", "*", "%"],
            nr_code=2,
            void_code=3,
        ),
        "all_void": STSList(
            data=np.full((3, 2), 3, dtype=int),
            alphabet=["A", "B", "*", "%"],
            labels=["A", "B", "*", "%"],
            nr_code=2,
            void_code=3,
        ),
    }


def _golden(case_name: str, start: str, result_type: str) -> np.ndarray:
    records = pd.read_csv(R_GOLDENS_DIR / "sortv_r3_cases.csv")
    selected = records.loc[
        (records["case"] == case_name)
        & (records["start"] == start)
        & (records["result_type"] == result_type)
    ]
    # R row indexes are one-based; this is the sole representation conversion.
    return selected.sort_values("position")["value"].to_numpy(dtype=int) - 1


@pytest.mark.parametrize("case_name", list(_cases()))
@pytest.mark.parametrize("start", ["beg", "end"])
def test_sorti_indexes_and_inverse_match_r_without_mutation(case_name: str, start: str) -> None:
    seqdata = _cases()[case_name]
    before = seqdata.data.copy()

    indexes = stats_utils.sorti(seqdata, start=start)
    inverse = stats_utils.sorti(seqdata, start=start, sort_index=False)

    np.testing.assert_array_equal(indexes, _golden(case_name, start, "sorti"))
    np.testing.assert_array_equal(inverse, _golden(case_name, start, "inverse"))
    np.testing.assert_array_equal(inverse, np.argsort(indexes, kind="stable"))
    np.testing.assert_array_equal(seqdata.data, before)
    assert indexes.dtype == np.int64
    assert inverse.dtype == np.int64


@pytest.mark.parametrize("case_name", list(_cases()))
@pytest.mark.parametrize("start", ["beg", "end"])
def test_sortv_matches_r_inverse_sort_values(case_name: str, start: str) -> None:
    result = stats_utils.sortv(_cases()[case_name], start=start)
    np.testing.assert_array_equal(result, _golden(case_name, start, "sortv"))


def test_sorti_invalid_start_matches_r_warning_and_zero_based_result() -> None:
    conditions = (R_GOLDENS_DIR / "sortv_r3_conditions.txt").read_text(encoding="utf-8")
    assert 'invalid_start|warning|start should be one of "beg" or "end"' in conditions
    assert "invalid_start|result|1,2,3" in conditions

    with pytest.warns(UserWarning, match='start should be one of "beg" or "end"'):
        result = stats_utils.sorti(_cases()["declared_order"], start="middle")
    np.testing.assert_array_equal(result, np.arange(3, dtype=np.int64))


def test_sorting_zero_width_matches_r_error_boundary() -> None:
    conditions = (R_GOLDENS_DIR / "sortv_r3_conditions.txt").read_text(encoding="utf-8")
    assert "zero_width_sorti|simpleError|second argument must be a list" in conditions
    assert "zero_width_sortv|simpleError|second argument must be a list" in conditions
    empty = STSList(data=np.empty((3, 0), dtype=int), alphabet=["A", "B"])

    with pytest.raises(ValueError, match="zero-width"):
        stats_utils.sorti(empty)
    with pytest.raises(ValueError, match="zero-width"):
        stats_utils.sortv(empty)


def test_sorting_helpers_are_public_r_derived_names() -> None:
    assert traminex.sorti is stats_utils.sorti
    assert traminex.sortv is stats_utils.sortv
