"""TraMineRextras 0.6.9 goldens for event-driven ``TSE_to_STS``."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from traminex import STSList, TSE_to_STS, alphabet

GOLDEN_DIR = Path(__file__).resolve().parent / "r_goldens"
CASES = pd.read_csv(
    GOLDEN_DIR / "tse_to_sts_cases.csv",
    dtype=str,
    keep_default_na=False,
)


def _records(kind: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in (GOLDEN_DIR / "tse_to_sts_errors.txt").read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        record_kind, name, message = line.split("|", 2)
        if record_kind == kind:
            assert name not in out
            out[name] = message.replace(r"\n", "\n")
    return out


R_ERRORS = _records("ERROR")
DIVERGENCES = _records("DIVERGENCE")
CASE_NAMES = tuple(CASES["case"].drop_duplicates())
# R's alphabet(seqdef(result)) per case: the default alphabet since owner decision D2.
_ALPHABET_ROWS = pd.read_csv(
    GOLDEN_DIR / "tse_to_sts_alphabets.csv", dtype={"state": str}, keep_default_na=False
)
ALPHABETS = {
    name: group.sort_values("position")["state"].tolist()
    for name, group in _ALPHABET_ROWS.groupby("case", sort=False)
}


def _canonical_tse() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "id": ["z", "z", "z", "a", "a", "a"],
            "time": [2.75, 0, 2.75, 1, 3, 1],
            "event": ["B", "A", "C", "A", "C", "B"],
        }
    )


def _explicit_stm() -> pd.DataFrame:
    # TraMineRextras::seqe2stm(c("A", "B"), dropList=list(B="A"),
    #                          firstState="origin")
    return pd.DataFrame(
        {
            "A": ["A", "A", "A.B", "A.B"],
            "B": ["B", "B", "B", "B"],
        },
        index=["origin", "A", "A.B", "B"],
    )


def _call_case(name: str) -> STSList:
    if name == "default_stm_default_selectors":
        return TSE_to_STS(_canonical_tse(), tmax=5)
    if name == "default_stm_named_selectors":
        named = _canonical_tse().rename(
            columns={"id": "person", "time": "moment", "event": "action"}
        )
        return TSE_to_STS(
            named,
            id_col="person",
            time_col="moment",
            event_col="action",
            tmax=5,
        )
    if name == "explicit_stm_ties_unknown_first_state":
        data = pd.DataFrame(
            {
                "id": ["p", "p", "p", "q", "q"],
                "time": [0, 0, 3, 1, 4],
                "event": ["B", "A", "Z", "B", "A"],
            }
        )
        return TSE_to_STS(
            data,
            event_col="event",
            stm=_explicit_stm(),
            tmin=1,
            tmax=5,
            first_state="origin",
        )
    if name == "numeric_selectors_tmin_slice":
        data = pd.DataFrame(
            {
                "noise": ["x", "x", "y", "y"],
                "action": ["A", "B", "B", "A"],
                "pid": [20, 20, 10, 10],
                "moment": [0, 2, 0, 3],
            }
        )
        return TSE_to_STS(
            data,
            id_col=3,
            time_col=4,
            event_col=2,
            tmin=2,
            tmax=4,
        )
    if name == "late_event_clamps_to_tmax":
        data = pd.DataFrame(
            {
                "id": ["p", "p"],
                "time": [0, 99],
                "event": ["A", "B"],
            }
        )
        return TSE_to_STS(data, tmax=4)
    if name == "infinite_event_clamps_to_tmax":
        data = pd.DataFrame({"id": ["p"], "time": [np.inf], "event": ["A"]})
        stm = pd.DataFrame({"A": ["A", "A"]}, index=["origin", "A"])
        return TSE_to_STS(data, stm=stm, tmax=4, first_state="origin")
    if name == "negative_boundary_event_applies_from_start":
        data = pd.DataFrame(
            {
                "id": ["p", "q"],
                "time": [-0.25, -1],
                "event": ["A", "A"],
            }
        )
        stm = pd.DataFrame({"A": ["A", "A"]}, index=["origin", "A"])
        return TSE_to_STS(data, stm=stm, tmax=3, first_state="origin")
    if name == "marker_text_states_remain_real":
        data = pd.DataFrame(
            {
                "id": ["p", "q"],
                "time": [0, 0],
                "event": ["*", "%"],
            }
        )
        return TSE_to_STS(data, tmax=2)
    if name == "auto_stm_custom_first_state_safe":
        data = pd.DataFrame(
            {
                "id": ["p", "p"],
                "time": [0, 2],
                "event": ["A", "B"],
            }
        )
        return TSE_to_STS(data, tmax=4, first_state="origin")
    if name == "one_event_auto_stm_transitions":
        data = pd.DataFrame(
            {
                "id": ["p", "p"],
                "time": [0, 2],
                "event": ["A", "A"],
            }
        )
        return TSE_to_STS(data, tmax=4)
    if name == "noncontiguous_ids_group_stably":
        data = pd.DataFrame(
            {
                "id": ["p", "q", "p", "q"],
                "time": [0, 0, 2, 2],
                "event": ["A", "B", "B", "A"],
            }
        )
        return TSE_to_STS(data, tmax=4)
    raise AssertionError(f"unhandled TSE golden case: {name}")


def _case(name: str) -> pd.DataFrame:
    rows = CASES.loc[CASES["case"] == name].copy()
    assert not rows.empty
    rows["row"] = rows["row"].astype(int)
    rows["column_index"] = rows["column_index"].astype(int)
    return rows.sort_values(["row", "column_index"], kind="stable").reset_index(drop=True)


def _decoded(x: STSList) -> np.ndarray:
    return np.asarray(x.alphabet, dtype=object)[x.data]


def test_state_at_time_keeps_ids_aligned_in_first_appearance_order():
    # R's route keeps the persons in first-appearance order (TX-AUD-20261019-039).
    data = pd.DataFrame(
        {
            "id": ["case-z", "case-a", "case-z", "case-a"],
            "time": [1, 1, 2, 2],
            "state": ["A", "B", "C", "A"],
        }
    )

    actual = TSE_to_STS(data)

    assert actual.ids == ["case-z", "case-a"]
    assert dict(zip(actual.ids, _decoded(actual).tolist(), strict=True)) == {
        "case-a": ["B", "A"],
        "case-z": ["A", "C"],
    }


def _assert_case(actual: STSList, name: str) -> None:
    expected = _case(name)
    expected_nrow = int(expected["nrow"].iloc[0])
    expected_ncol = int(expected["ncol"].iloc[0])
    expected_ids = expected.groupby("row", sort=False)["id"].first().tolist()
    expected_columns = (
        expected.loc[expected["row"] == 1]
        .sort_values("column_index", kind="stable")["column"]
        .tolist()
    )
    expected_states = expected["state"].to_numpy(dtype=object).reshape(expected_nrow, expected_ncol)

    assert isinstance(actual, STSList)
    assert actual.data.shape == (expected_nrow, expected_ncol)
    np.testing.assert_array_equal(_decoded(actual), expected_states)
    assert actual.ids == expected_ids
    assert actual.xtlab == expected_columns
    assert actual.start == int(expected_columns[0][1:])
    # Owner decision D2 (2026-10-19): R's seqdef order of the states, recorded from R;
    # until then this compared with first appearance (TX-AUD-20261019-071).
    assert alphabet(actual) == ALPHABETS[name]


def test_golden_artifact_covers_the_event_contract_matrix():
    assert set(CASE_NAMES) == {
        "default_stm_default_selectors",
        "default_stm_named_selectors",
        "explicit_stm_ties_unknown_first_state",
        "numeric_selectors_tmin_slice",
        "late_event_clamps_to_tmax",
        "infinite_event_clamps_to_tmax",
        "negative_boundary_event_applies_from_start",
        "marker_text_states_remain_real",
        "auto_stm_custom_first_state_safe",
        "one_event_auto_stm_transitions",
        "noncontiguous_ids_group_stably",
    }
    assert set(CASES["contract"]) == {"R_EXACT", "SELECTED_DIVERGENCE"}
    assert set(CASES.loc[CASES["contract"] == "SELECTED_DIVERGENCE", "case"]) == set(DIVERGENCES)
    assert set(R_ERRORS) == {
        "tmax_required",
        "tmin_below_one",
        "unknown_named_id_selector",
        "unknown_numeric_event_selector",
        "empty_auto_stm",
    }


def test_divergence_artifact_documents_each_selected_safety_contract():
    assert set(DIVERGENCES) == {
        "auto_stm_custom_first_state_safe",
        "one_event_auto_stm_transitions",
        "noncontiguous_ids_group_stably",
    }
    assert "subscript out of bounds" in DIVERGENCES["auto_stm_custom_first_state_safe"]
    assert "states=None,None,None,None" in DIVERGENCES["one_event_auto_stm_transitions"]
    assert "subscript out of bounds" in DIVERGENCES["noncontiguous_ids_group_stably"]
    assert all("selected Python" in message for message in DIVERGENCES.values())


@pytest.mark.parametrize("name", CASE_NAMES)
def test_event_driven_cases_match_r(name: str):
    _assert_case(_call_case(name), name)


def test_raw_event_time_becomes_state_at_floor_time_plus_one():
    actual = _call_case("default_stm_default_selectors")
    states = _decoded(actual)
    # z has A at raw time 0, then B and C at raw time 2.75.
    assert states[0].tolist() == ["A", "A", "B.A.C", "B.A.C", "B.A.C"]


def test_tied_events_sort_by_event_name_and_unknown_events_are_ignored():
    actual = _call_case("explicit_stm_ties_unknown_first_state")
    states = _decoded(actual)
    # p's input tie is B then A, but lexical processing applies A then B. The
    # explicit B transition drops A, and the later unknown Z has no effect.
    assert states[0].tolist() == ["B", "B", "B", "B", "B"]


def test_tmin_slices_r_columns_and_numeric_ids_become_row_labels():
    actual = _call_case("numeric_selectors_tmin_slice")
    assert actual.ids == ["20", "10"]
    assert actual.xtlab == ["a2", "a3", "a4"]
    assert actual.start == 2
    _assert_case(actual, "numeric_selectors_tmin_slice")


def test_event_after_horizon_is_clamped_into_final_cell():
    actual = _call_case("late_event_clamps_to_tmax")
    assert _decoded(actual).tolist() == [["A", "A", "A", "A.B"]]
    _assert_case(actual, "late_event_clamps_to_tmax")


def test_positive_infinite_event_is_clamped_into_final_cell():
    actual = _call_case("infinite_event_clamps_to_tmax")
    assert _decoded(actual).tolist() == [["origin", "origin", "origin", "A"]]


def test_negative_boundary_event_applies_from_first_cell():
    actual = _call_case("negative_boundary_event_applies_from_start")
    assert _decoded(actual).tolist() == [["A", "A", "A"], ["A", "A", "A"]]


def test_marker_text_states_are_not_reclassified_as_missing_or_void():
    actual = _call_case("marker_text_states_remain_real")
    assert _decoded(actual).tolist() == [["*", "*"], ["%", "%"]]
    assert actual.nr_code is None
    assert actual.void_code is None
    assert alphabet(actual) == ["*", "%"]
    assert actual.nr not in actual.alphabet
    assert actual.void not in actual.alphabet


def test_auto_stm_honors_custom_first_state_without_r_index_failure():
    actual = _call_case("auto_stm_custom_first_state_safe")
    assert _decoded(actual).tolist() == [["A", "A", "A.B", "A.B"]]
    _assert_case(actual, "auto_stm_custom_first_state_safe")


def test_one_event_auto_stm_keeps_dimensions_and_applies_the_event():
    actual = _call_case("one_event_auto_stm_transitions")
    assert _decoded(actual).tolist() == [["A", "A", "A", "A"]]
    _assert_case(actual, "one_event_auto_stm_transitions")


def test_noncontiguous_ids_are_grouped_once_in_stable_first_appearance_order():
    actual = _call_case("noncontiguous_ids_group_stably")
    assert actual.ids == ["p", "q"]
    assert _decoded(actual).tolist() == [
        ["A", "A", "A.B", "A.B"],
        ["B", "B", "A.B", "A.B"],
    ]
    _assert_case(actual, "noncontiguous_ids_group_stably")


def test_explicit_state_column_keeps_state_at_time_compatibility_path():
    data = pd.DataFrame(
        {
            "id": ["p", "p", "q", "q"],
            "time": [1, 2, 1, 2],
            "state": ["left", "right", "down", "down"],
            # The event column is deliberately incompatible with the states.
            "event": ["A", "B", "A", "B"],
        }
    )

    actual = TSE_to_STS(data, id_col="id", time_col="time", state_col="state")

    assert _decoded(actual).tolist() == [["left", "right"], ["down", "down"]]
    assert "A.B" not in actual.alphabet


def test_event_mode_rejects_overrides_of_derived_metadata():
    with pytest.raises(TypeError, match="derives ids, xtlab, and start"):
        TSE_to_STS(_canonical_tse(), tmax=5, ids=["x", "y"])


def test_event_mode_preserves_explicitly_disabled_markers():
    actual = TSE_to_STS(_canonical_tse(), tmax=5, nr=None, void=None)
    assert actual.nr is None
    assert actual.void is None
    assert actual.nr_code is None
    assert actual.void_code is None


def _error_call(name: str) -> STSList:
    data = _canonical_tse()
    if name == "tmax_required":
        return TSE_to_STS(data)
    if name == "tmin_below_one":
        return TSE_to_STS(data, tmin=0, tmax=5)
    if name == "unknown_named_id_selector":
        return TSE_to_STS(data, id_col="missing", event_col="event", tmax=5)
    if name == "unknown_numeric_event_selector":
        return TSE_to_STS(data, id_col=1, time_col=2, event_col=99, tmax=5)
    if name == "empty_auto_stm":
        return TSE_to_STS(data.iloc[0:0], tmax=3)
    raise AssertionError(f"unhandled TSE error case: {name}")


@pytest.mark.parametrize("name", tuple(R_ERRORS))
def test_event_driven_errors_match_r(name: str):
    with pytest.raises(ValueError) as exc_info:
        _error_call(name)
    assert str(exc_info.value) == R_ERRORS[name]
