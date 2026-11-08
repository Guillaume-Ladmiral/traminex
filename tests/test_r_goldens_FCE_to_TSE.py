"""TraMineRextras 0.6.9 goldens for ``FCE_to_TSE``."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from traminex.format import FCE_to_TSE

GOLDEN_DIR = Path(__file__).resolve().parent / "r_goldens"
CASES = pd.read_csv(
    GOLDEN_DIR / "fce_to_tse_cases.csv",
    dtype=str,
    keep_default_na=False,
)


def _records(kind: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in (GOLDEN_DIR / "fce_to_tse_errors.txt").read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        record_kind, name, message = line.split("|", 2)
        if record_kind == kind:
            assert name not in out
            out[name] = message.replace(r"\n", "\n")
    return out


R_ERRORS = _records("ERROR")
SELECTED_ERRORS = _records("SELECTED_ERROR")
DIVERGENCES = _records("DIVERGENCE")
CASE_NAMES = tuple(CASES["case"].drop_duplicates())


def _fce() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "pid": ["z", "a", "m"],
            "eB": [2.5, np.nan, 1.25],
            "eA": [2.5, 3.75, np.nan],
            "eC": [np.nan, 3.75, 1.25],
        }
    )


def _empty_observed() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "pid": ["z", "a"],
            "e1": [np.nan, np.nan],
            "e2": [np.nan, np.nan],
        }
    )


def _eventless_mixed() -> pd.DataFrame:
    return pd.DataFrame({"pid": ["z", "a", "m"], "e1": [2.0, np.nan, 1.0]})


def _call_case(name: str) -> pd.DataFrame:
    data = _fce()
    if name == "null_id_character_cols_implicit":
        return FCE_to_TSE(data, cols=["eB", "eA", "eC"])
    if name == "single_row_null_id_uses_index":
        return FCE_to_TSE(data.iloc[[0]], cols=["eB", "eA"])
    if name == "id_column_numeric_cols_implicit":
        return FCE_to_TSE(data, id="pid", cols=[2, 3, 4])
    if name == "id_numeric_selector":
        return FCE_to_TSE(data, id=1, cols=["eB", "eA", "eC"])
    if name == "explicit_character_ids":
        return FCE_to_TSE(data, id=["q", "r", "s"], cols=["eB", "eA"])
    if name == "explicit_integer_ids":
        return FCE_to_TSE(data, id=[30, 10, 20], cols=["eB", "eA"])
    if name == "explicit_events_first_event":
        return FCE_to_TSE(
            data,
            id="pid",
            cols=["eB", "eA", "eC"],
            eventlist=["beta", "alpha", "gamma"],
            firstEvent="birth",
        )
    if name == "all_missing_without_first_event":
        return FCE_to_TSE(_empty_observed(), id="pid", cols=["e1", "e2"])
    if name == "eventless_person_does_not_corrupt_ids":
        return FCE_to_TSE(_eventless_mixed(), id="pid", cols=["e1"])
    if name == "all_missing_with_first_event":
        return FCE_to_TSE(
            _empty_observed(),
            id="pid",
            cols=["e1", "e2"],
            firstEvent="origin",
        )
    raise AssertionError(f"unhandled FCE golden case: {name}")


def _case(name: str) -> pd.DataFrame:
    rows = CASES.loc[CASES["case"] == name].copy()
    assert not rows.empty
    rows["row"] = rows["row"].astype(int)
    return rows.sort_values("row", kind="stable").reset_index(drop=True)


def _render_cell(value: object) -> str:
    try:
        missing = pd.isna(value)
    except (TypeError, ValueError):
        missing = False
    if isinstance(missing, (bool, np.bool_)) and bool(missing):
        return "<NA>"
    if isinstance(value, (int, np.integer)) and not isinstance(value, (bool, np.bool_)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        numeric = float(value)
        if numeric.is_integer():
            return str(int(numeric))
        return format(numeric, ".15g")
    return str(value)


def _assert_case(actual: pd.DataFrame, name: str, *, time_column: str = "time") -> None:
    expected = _case(name)
    assert expected["columns"].unique().tolist() == ["id;time;event"]
    assert len(expected["contract"].unique()) == 1
    assert actual.columns.tolist() == ["id", time_column, "event"]
    assert len(actual) == int(expected["nrow"].iloc[0])
    if actual.empty:
        assert expected["row"].tolist() == [0]
        assert expected[["id", "time", "event"]].values.tolist() == [
            ["<EMPTY>", "<EMPTY>", "<EMPTY>"]
        ]
        return

    actual_rows = [
        [_render_cell(value) for value in row]
        for row in actual[["id", time_column, "event"]].itertuples(index=False, name=None)
    ]
    assert actual_rows == expected[["id", "time", "event"]].values.tolist()

    id_type = expected["id_type"].unique().tolist()
    assert len(id_type) == 1
    if id_type[0] == "numeric":
        assert pd.api.types.is_numeric_dtype(actual["id"].dtype)
    else:
        nonmissing_ids = actual.loc[actual["id"].notna(), "id"]
        assert all(isinstance(value, str) for value in nonmissing_ids)


def test_golden_artifact_covers_the_complete_selector_and_empty_matrix():
    assert set(CASE_NAMES) == {
        "null_id_character_cols_implicit",
        "single_row_null_id_uses_index",
        "id_column_numeric_cols_implicit",
        "id_numeric_selector",
        "explicit_character_ids",
        "explicit_integer_ids",
        "explicit_events_first_event",
        "all_missing_without_first_event",
        "eventless_person_does_not_corrupt_ids",
        "all_missing_with_first_event",
    }
    assert set(CASES["columns"]) == {"id;time;event"}
    assert set(CASES["time_type"]) == {"numeric"}
    assert set(CASES["event_type"]) == {"character"}
    assert set(CASES["contract"]) == {"R_EXACT", "SELECTED_DIVERGENCE"}
    assert set(CASES.loc[CASES["contract"] == "SELECTED_DIVERGENCE", "case"]) == {
        "all_missing_without_first_event",
        "eventless_person_does_not_corrupt_ids",
        "single_row_null_id_uses_index",
    }


def test_divergence_artifact_documents_each_selected_contract():
    assert set(DIVERGENCES) == {
        "all_missing_without_first_event",
        "eventless_person_does_not_corrupt_ids",
        "single_row_null_id_uses_index",
        "eventlist_too_short",
        "eventlist_too_long",
    }
    assert "1:0" in DIVERGENCES["all_missing_without_first_event"]
    assert "ids=z,m" in DIVERGENCES["eventless_person_does_not_corrupt_ids"]
    assert all("selected Python" in message for message in DIVERGENCES.values())


@pytest.mark.parametrize("name", CASE_NAMES)
def test_fce_to_tse_cases_match_r(name: str):
    _assert_case(_call_case(name), name)


def test_all_missing_input_is_a_true_empty_result_not_the_r_dummy_row():
    actual = _call_case("all_missing_without_first_event")
    assert actual.empty
    _assert_case(actual, "all_missing_without_first_event")


def test_eventless_person_contributes_no_rows_without_corrupting_other_ids():
    actual = _call_case("eventless_person_does_not_corrupt_ids")
    assert actual["id"].tolist() == ["z", "m"]
    _assert_case(actual, "eventless_person_does_not_corrupt_ids")


def test_first_event_is_at_zero_and_ties_sort_by_time_then_event_per_person():
    actual = _call_case("explicit_events_first_event")
    assert actual.loc[actual["event"] == "birth", "time"].tolist() == [0, 0, 0]
    assert actual["id"].drop_duplicates().tolist() == ["z", "a", "m"]
    assert actual.loc[actual["id"] == "z", "event"].tolist() == ["birth", "alpha", "beta"]
    _assert_case(actual, "explicit_events_first_event")


def test_timestamp_compatibility_renames_only_the_time_column():
    actual = FCE_to_TSE(
        _fce(),
        id="pid",
        cols=["eB", "eA", "eC"],
        eventlist=["beta", "alpha", "gamma"],
        firstEvent="birth",
        time_col="timestamp",
    )
    _assert_case(actual, "explicit_events_first_event", time_column="timestamp")


def test_integer_id_column_is_not_promoted_by_float_event_times():
    data = pd.DataFrame({"pid": [20, 10], "event_time": [2.5, 1.25]})
    actual = FCE_to_TSE(data, id="pid", cols=["event_time"])
    assert actual["id"].tolist() == [20, 10]
    assert pd.api.types.is_integer_dtype(actual["id"].dtype)
    assert actual["time"].tolist() == [2.5, 1.25]


def _error_calls() -> dict[str, object]:
    data = _fce()
    return {
        "unknown_character_cols": lambda: FCE_to_TSE(
            data,
            id="pid",
            cols=["not_a_column"],
        ),
        "numeric_col_out_of_bounds": lambda: FCE_to_TSE(data, id="pid", cols=[99]),
        "unknown_id_selector": lambda: FCE_to_TSE(data, id="not_an_id", cols=["eB"]),
        "empty_cols": lambda: FCE_to_TSE(data, id="pid", cols=[]),
        "zero_row_input": lambda: FCE_to_TSE(
            data.iloc[0:0],
            id="pid",
            cols=["eB", "eA"],
        ),
        "eventlist_too_short": lambda: FCE_to_TSE(
            data,
            id="pid",
            cols=["eB", "eA"],
            eventlist=["only"],
        ),
        "eventlist_too_long": lambda: FCE_to_TSE(
            data,
            id="pid",
            cols=["eB"],
            eventlist=["kept", "ignored"],
        ),
    }


@pytest.mark.parametrize("name", sorted(R_ERRORS))
def test_fce_to_tse_oracle_refusals(name: str):
    calls = _error_calls()
    with pytest.raises(
        (IndexError, KeyError, TypeError, ValueError), match=re.escape(R_ERRORS[name])
    ):
        calls[name]()


@pytest.mark.parametrize("name", sorted(SELECTED_ERRORS))
def test_fce_to_tse_selected_refusals(name: str):
    calls = _error_calls()
    assert set(SELECTED_ERRORS) == {"eventlist_too_short", "eventlist_too_long"}
    with pytest.raises(ValueError, match=re.escape(SELECTED_ERRORS[name])):
        calls[name]()
