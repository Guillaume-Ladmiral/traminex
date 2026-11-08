from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDEN_ATOL

from traminex import STSList, toPersonPeriod

GOLDEN_DIR = Path(__file__).parent / "r_goldens"
CASES = pd.read_csv(
    GOLDEN_DIR / "to_person_period_cases.csv",
    keep_default_na=False,
)
CASE_NAMES = CASES["case"].drop_duplicates().tolist()


def _fixture(name: str) -> tuple[STSList, str]:
    if name == "custom_ids_ragged_void":
        return (
            STSList(
                data=np.array(
                    [
                        [1, 0, 4, 4],
                        [2, 1, 0, 4],
                        [0, 2, 1, 2],
                    ],
                    dtype=int,
                ),
                alphabet=["A", "B", "C", "*", "%"],
                nr="*",
                void="%",
                nr_code=3,
                void_code=4,
                start=7,
                ids=["zeta", "alpha", "mu"],
            ),
            "timestamp",
        )
    if name == "interior_missing":
        return (
            STSList(
                data=np.array(
                    [
                        [0, 3, 1, 2],
                        [1, 2, 3, 0],
                        [2, 0, 1, 2],
                    ],
                    dtype=int,
                ),
                alphabet=["A", "B", "C", "*", "%"],
                nr="*",
                void="%",
                nr_code=3,
                void_code=4,
                start=4,
                ids=["zeta", "alpha", "mu"],
            ),
            "timestamp",
        )
    if name == "no_missing":
        return (
            STSList(
                data=np.array(
                    [
                        [2, 1, 0],
                        [0, 2, 1],
                        [1, 0, 2],
                    ],
                    dtype=int,
                ),
                alphabet=["A", "B", "C"],
                start=12,
                ids=["zeta", "alpha", "mu"],
            ),
            "timestamp",
        )
    if name == "custom_time_col":
        return (
            STSList(
                data=np.array([[1, 0], [0, 2]], dtype=int),
                alphabet=["A", "B", "C"],
                start=20,
                ids=["person-b", "person-a"],
            ),
            "period",
        )
    raise AssertionError(f"unknown golden case {name!r}")


def _expected(name: str) -> tuple[pd.DataFrame, str]:
    rows = CASES.loc[CASES["case"] == name].copy()
    rows["row"] = rows["row"].astype(int)
    rows["timestamp"] = rows["timestamp"].astype(int)
    rows = rows.sort_values("row", kind="stable")
    time_col = str(rows["time_col"].iloc[0])
    expected = rows[["id", "state", "timestamp"]].rename(columns={"timestamp": time_col})
    expected["state"] = expected["state"].mask(expected["state"] == "<NA>", pd.NA)
    return expected.reset_index(drop=True), time_col


def _state_values(series: pd.Series) -> list[object]:
    return ["<NA>" if pd.isna(value) else value for value in series]


def _assert_case(name: str) -> pd.DataFrame:
    sequence, requested_time_col = _fixture(name)
    expected, golden_time_col = _expected(name)
    assert requested_time_col == golden_time_col

    actual = toPersonPeriod(sequence, time_col=requested_time_col)

    assert list(actual.columns) == ["id", "state", requested_time_col]
    assert actual["id"].tolist() == expected["id"].tolist()
    assert actual[requested_time_col].tolist() == expected[requested_time_col].tolist()
    assert _state_values(actual["state"]) == _state_values(expected["state"])
    return actual


def test_golden_artifact_covers_the_full_person_period_contract():
    assert set(CASE_NAMES) == {
        "custom_ids_ragged_void",
        "interior_missing",
        "no_missing",
        "custom_time_col",
    }
    assert set(CASES["contract"]) == {"R_EXACT"}
    assert set(CASES["columns"]) == {"id;state;timestamp", "id;state;period"}


@pytest.mark.parametrize("name", CASE_NAMES)
def test_person_period_cases_match_r(name: str):
    _assert_case(name)


def test_custom_ids_and_rows_are_sorted_by_id_timestamp_state():
    actual = _assert_case("custom_ids_ragged_void")
    sort_probe = actual.assign(
        _state=actual["state"].map(lambda value: "<NA>" if pd.isna(value) else value)
    )
    sorted_probe = sort_probe.sort_values(["id", "timestamp", "_state"], kind="stable").reset_index(
        drop=True
    )
    pd.testing.assert_frame_equal(
        sort_probe.reset_index(drop=True),
        sorted_probe,
        check_exact=False,
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )
    assert actual["id"].drop_duplicates().tolist() == ["alpha", "mu", "zeta"]


def test_ragged_void_rows_are_removed_and_each_timestamp_starts_at_zero():
    actual = _assert_case("custom_ids_ragged_void")
    assert "%" not in actual["state"].tolist()
    assert len(actual) == 9
    assert actual.groupby("id", sort=False)["timestamp"].apply(list).to_dict() == {
        "alpha": [0, 1, 2],
        "mu": [0, 1, 2, 3],
        "zeta": [0, 1],
    }


def test_nonresponse_marker_becomes_missing_state():
    actual = _assert_case("interior_missing")
    missing_rows = actual.loc[actual["state"].isna(), ["id", "timestamp"]]
    assert missing_rows.to_dict("records") == [
        {"id": "alpha", "timestamp": 2},
        {"id": "zeta", "timestamp": 1},
    ]
    assert "*" not in _state_values(actual["state"])


def test_no_missing_case_has_dense_zero_based_timestamps():
    actual = _assert_case("no_missing")
    assert actual["state"].notna().all()
    assert actual.groupby("id", sort=False)["timestamp"].apply(list).to_dict() == {
        "alpha": [0, 1, 2],
        "mu": [0, 1, 2],
        "zeta": [0, 1, 2],
    }


def test_time_column_can_be_renamed_without_changing_r_values():
    actual = _assert_case("custom_time_col")
    assert list(actual.columns) == ["id", "state", "period"]
    assert actual["period"].tolist() == [0, 1, 0, 1]


def test_marker_text_states_are_preserved_by_code_aware_filtering():
    sequence = STSList(
        data=np.array([[0, 1, 2, 3]], dtype=int),
        alphabet=["*", "%", "<NR>", "<VOID>"],
        nr="<NR>",
        void="<VOID>",
        nr_code=2,
        void_code=3,
        ids=["p"],
    )

    actual = toPersonPeriod(sequence)

    assert _state_values(actual["state"]) == ["*", "%", "<NA>"]
    assert actual["timestamp"].tolist() == [0, 1, 2]


def test_equal_length_tuple_ids_remain_scalar_identities():
    sequence = STSList(
        data=np.array([[0], [1]], dtype=int),
        alphabet=["A", "B"],
        ids=[(2, "b"), (1, "a")],
    )

    actual = toPersonPeriod(sequence)

    assert actual["id"].tolist() == [(1, "a"), (2, "b")]
    assert actual["state"].tolist() == ["B", "A"]


def test_missing_ids_fall_back_to_one_based_integers():
    sequence = STSList(
        data=np.array([[0], [1]], dtype=int),
        alphabet=["A", "B"],
    )
    assert toPersonPeriod(sequence)["id"].tolist() == [1, 2]


def test_default_ids_follow_r_rowname_order_beyond_nine():
    sequence = STSList(
        data=np.zeros((12, 1), dtype=int),
        alphabet=["A"],
    )
    assert toPersonPeriod(sequence)["id"].tolist() == [1, 10, 11, 12, 2, 3, 4, 5, 6, 7, 8, 9]


def test_duplicate_ids_use_r_factor_level_state_order():
    sequence = STSList(
        data=np.array([[0], [1]], dtype=int),
        alphabet=["Z", "A"],
        ids=["p", "p"],
    )
    assert toPersonPeriod(sequence)["state"].tolist() == ["Z", "A"]


@pytest.mark.parametrize("time_col", ["id", "state"])
def test_time_column_cannot_overwrite_required_r_columns(time_col: str):
    sequence, _ = _fixture("no_missing")
    with pytest.raises(ValueError, match="time_col must differ"):
        toPersonPeriod(sequence, time_col=time_col)


@pytest.mark.parametrize("time_col", ["__id_order", "__state_order"])
def test_time_column_can_use_former_internal_sort_names(time_col: str):
    sequence, _ = _fixture("custom_time_col")
    actual = toPersonPeriod(sequence, time_col=time_col)
    assert list(actual.columns) == ["id", "state", time_col]
    assert actual[time_col].tolist() == [0, 1, 0, 1]
