"""R-golden coverage for TraMineRextras 0.6.9 ``seqsamm``."""

from __future__ import annotations

from dataclasses import is_dataclass

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDEN_ATOL, R_GOLDENS_DIR

import traminex
from traminex.samm import SAMMResult
from traminex.sequences import STSList


def _ragged_fixture() -> tuple[STSList, pd.DataFrame]:
    x = STSList(
        data=np.array(
            [
                [0, 0, 1, 1, 2, 2, 2],
                [0, 1, 1, 2, 4, 4, 4],
                [2, 2, 2, 2, 2, 2, 2],
                [0, 0, 0, 3, 1, 1, 4],
            ],
            dtype=int,
        ),
        alphabet=["Alpha", "Beta", "Gamma"],
        labels=["Alpha", "Beta", "Gamma"],
        cpal=["#111111", "#222222", "#333333"],
        xtlab=["t1", "t2", "t3", "t4", "t5", "t6", "t7"],
        xtstep=2,
        tick_last=True,
        nr="*",
        void="%",
        nr_code=3,
        void_code=4,
        weights=np.array([1.0, 2.0, 3.0, 4.0]),
        start=5,
        ids=["one", "two", "three", "four"],
        missing_color="#999999",
        issues={"example": [2, 4]},
    )
    covar = pd.DataFrame(
        {
            "age": np.array([20, 30, 40, 50], dtype=np.int64),
            "sector": pd.Categorical(["x", "y", "x", "z"], categories=["z", "x", "y"]),
        },
        index=[10, 20, 30, 40],
    )
    return x, covar


def _constant_fixture() -> STSList:
    return STSList(
        data=np.array([[0, 0, 0, 0, 0], [1, 1, 1, 1, 1]], dtype=int),
        alphabet=["A", "B"],
        labels=["A", "B"],
        cpal=["#7FC97F", "#BEAED4"],
        nr="*",
        void="%",
        nr_code=2,
        void_code=3,
        start=3,
    )


def _repeated_fixture() -> tuple[STSList, pd.DataFrame]:
    x = STSList(
        data=np.array([[0, 1, 0, 0, 2, 2, 0], [0, 0, 1, 1, 0, 2, 2]], dtype=int),
        alphabet=["State A", "State B", "State C"],
        labels=["State A", "State B", "State C"],
        cpal=["#AA0000", "#00AA00", "#0000AA"],
        nr="*",
        void="%",
        nr_code=3,
        void_code=4,
        start=4,
    )
    covar = pd.DataFrame(
        {
            "score": np.array([10, 20], dtype=np.int64),
            "cohort": pd.Categorical(["old", "new"], categories=["new", "old"]),
        }
    )
    return x, covar


def _frame_records(case_name: str) -> pd.DataFrame:
    records = pd.read_csv(R_GOLDENS_DIR / "samm_r3_frame.csv", dtype=str, keep_default_na=False)
    return records.loc[records["case"] == case_name].reset_index(drop=True)


def _golden_frame(case_name: str) -> pd.DataFrame:
    records = _frame_records(case_name)
    columns = records.loc[records["record_type"] == "column"].sort_values(
        "value", key=lambda values: values.astype(int)
    )
    column_order = columns["column"].tolist()
    row_order = list(dict.fromkeys(records.loc[records["record_type"] == "cell", "row"]))
    result = pd.DataFrame(index=row_order)
    cells = records.loc[records["record_type"] == "cell"]
    for record in columns.itertuples(index=False):
        values = [
            cells.loc[(cells["row"] == row) & (cells["column"] == record.column), "value"].iloc[0]
            for row in row_order
        ]
        level_record = records.loc[
            (records["record_type"] == "levels") & (records["column"] == record.column),
            "value",
        ]
        if not level_record.empty:
            categories = level_record.iloc[0].split("\x1e")
            result[record.column] = pd.Categorical(values, categories=categories)
        elif record.detail == "integer":
            result[record.column] = pd.Series(values, index=row_order, dtype="int64")
        elif record.detail == "double":
            result[record.column] = pd.Series(values, index=row_order, dtype="float64")
        elif record.detail == "logical":
            result[record.column] = pd.Series(
                [value == "TRUE" for value in values], index=row_order, dtype="bool"
            )
        else:
            raise AssertionError(f"unsupported R column type: {record.detail}")
    return result.loc[:, column_order].reset_index(drop=True)


def _metadata(case_name: str, section: str) -> dict[str, str]:
    records = pd.read_csv(R_GOLDENS_DIR / "samm_r3_metadata.csv", dtype=str, keep_default_na=False)
    selected = records.loc[(records["case"] == case_name) & (records["section"] == section)]
    return dict(zip(selected["key"], selected["value"], strict=True))


def _golden_sequence_codes(case_name: str) -> np.ndarray:
    shape = _metadata(case_name, "sequences")
    codes = _metadata(case_name, "sequence_code")
    result = np.empty((int(shape["nrow"]), int(shape["ncol"])), dtype=int)
    for key, value in codes.items():
        row, column = (int(part) - 1 for part in key.split(","))
        result[row, column] = int(value)
    return result


def test_seqsamm_frame_matches_every_r_row_and_preserves_callers():
    x, covar = _ragged_fixture()
    data_snapshot = x.data.copy()
    weights_snapshot = x.weights.copy()
    covar_snapshot = covar.copy(deep=True)

    result = traminex.seqsamm(x, 3, covar=covar)

    assert isinstance(result, SAMMResult)
    pd.testing.assert_frame_equal(
        result.frame,
        _golden_frame("ragged_covariates"),
        check_exact=False,
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )
    np.testing.assert_array_equal(x.data, data_snapshot)
    np.testing.assert_array_equal(x.weights, weights_snapshot)
    pd.testing.assert_frame_equal(
        covar, covar_snapshot, check_exact=False, rtol=0, atol=R_GOLDEN_ATOL
    )


def test_seqsamm_result_and_embedded_sequence_metadata_match_r():
    x, covar = _ragged_fixture()
    result = traminex.seqsamm(x, 3, covar=covar)
    expected = _metadata("ragged_covariates", "sequences_attribute")

    assert result.spell == tuple(_metadata("ragged_covariates", "result")["spell"].split("\x1e"))
    assert result.state_columns == tuple(
        _metadata("ragged_covariates", "result")["sname"].split("\x1e")
    )
    assert result.sublength == 3
    np.testing.assert_array_equal(
        result.typology,
        _metadata("ragged_covariates", "result")["typology"].split("\x1e"),
    )
    np.testing.assert_array_equal(
        result.sequences.data, _golden_sequence_codes("ragged_covariates")
    )
    assert result.sequences.alphabet == expected["alphabet"].split("\x1e")
    assert result.sequences.labels == expected["labels"].split("\x1e")
    assert result.sequences.cpal == expected["cpal"].split("\x1e")
    assert result.sequences.start == int(expected["start"])
    assert result.sequences.xtstep == int(expected["xtstep"])
    assert result.sequences.tick_last is True
    assert result.sequences.nr == expected["nr"]
    assert result.sequences.void == expected["void"]
    assert result.sequences.nr_code == 3
    assert result.sequences.void_code == 4
    assert result.sequences.missing_color == expected["missing.color"]
    assert result.sequences.xtlab is None
    assert result.sequences.weights is None
    assert result.sequences.ids is None
    assert result.sequences.issues is None


def test_seqsamm_zero_transition_case_matches_r_and_is_not_empty():
    result = traminex.seqsamm(_constant_fixture(), 2)

    pd.testing.assert_frame_equal(
        result.frame, _golden_frame("no_transitions"), check_exact=False, rtol=0, atol=R_GOLDEN_ATOL
    )
    assert len(result.frame) == 6
    assert not result.frame["transition"].any()
    np.testing.assert_array_equal(result.sequences.data, _golden_sequence_codes("no_transitions"))
    np.testing.assert_array_equal(result.typology, np.repeat("None", 6))


@pytest.mark.parametrize("sublength", [True, 1, 3.0, 7, 8])
def test_seqsamm_refuses_invalid_sublengths(sublength):
    x, _ = _ragged_fixture()
    transcript = (R_GOLDENS_DIR / "samm_r3_errors.txt").read_text(encoding="utf-8")
    assert "sublength_one|simpleError" in transcript
    assert "sublength_width|simpleError" in transcript
    assert "sublength_fractional|success" in transcript

    with pytest.raises((TypeError, ValueError), match="sublength"):
        traminex.seqsamm(x, sublength)


def test_seqsamm_validates_inputs_and_covariate_shape():
    x, covar = _ragged_fixture()
    transcript = (R_GOLDENS_DIR / "samm_r3_errors.txt").read_text(encoding="utf-8")
    assert "covar_short|success" in transcript
    with pytest.raises(TypeError, match="STSList"):
        traminex.seqsamm(x.data, 3)
    with pytest.raises(TypeError, match="DataFrame"):
        traminex.seqsamm(x, 3, covar=np.asarray(covar))
    with pytest.raises(ValueError, match="one row per sequence"):
        traminex.seqsamm(x, 3, covar=covar.iloc[:-1])
    with pytest.raises(ValueError, match="reserved"):
        traminex.seqsamm(x, 3, covar=pd.DataFrame({"id": range(x.nseq)}))


def test_seqsamm_is_public_but_result_type_remains_module_level():
    assert traminex.seqsamm.__module__ == "traminex.samm"
    assert is_dataclass(SAMMResult)
    assert "seqsamm" in traminex.__all__
    assert "SAMMResult" not in traminex.__all__


def _consumer_sequence_records(case_name: str) -> pd.DataFrame:
    records = pd.read_csv(
        R_GOLDENS_DIR / "samm_consumers_r3_sequences.csv",
        dtype=str,
        keep_default_na=False,
    )
    return records.loc[records["case"] == case_name].reset_index(drop=True)


def _assert_consumer_sequence(result: STSList, case_name: str) -> None:
    records = _consumer_sequence_records(case_name)
    metadata = dict(
        zip(
            records.loc[records["record_type"] == "metadata", "column"],
            records.loc[records["record_type"] == "metadata", "value"],
            strict=True,
        )
    )
    expected = np.empty((int(metadata["nrow"]), int(metadata["ncol"])), dtype=int)
    for record in records.loc[records["record_type"] == "code"].itertuples(index=False):
        expected[int(record.row) - 1, int(record.column) - 1] = int(record.value)
    np.testing.assert_array_equal(result.data, expected)

    attributes = dict(
        zip(
            records.loc[records["record_type"] == "attribute", "column"],
            records.loc[records["record_type"] == "attribute", "value"],
            strict=True,
        )
    )
    assert result.start == int(attributes["start"])
    assert result.nr == attributes["nr"]
    assert result.void == attributes["void"]
    assert result.alphabet == attributes["alphabet"].split("\x1e")
    assert result.labels == attributes["labels"].split("\x1e")
    assert result.cpal == attributes["cpal"].split("\x1e")
    assert result.missing_color == attributes["missing.color"]
    assert result.xtstep == int(attributes["xtstep"])
    assert result.tick_last is (attributes["tick.last"] == "TRUE")


def _golden_eha_frame(case_name: str) -> pd.DataFrame:
    records = pd.read_csv(
        R_GOLDENS_DIR / "samm_consumers_r3_eha.csv",
        dtype=str,
        keep_default_na=False,
    )
    records = records.loc[records["case"] == case_name]
    columns = records.loc[records["record_type"] == "column"].sort_values(
        "value", key=lambda values: values.astype(int)
    )
    row_order = list(dict.fromkeys(records.loc[records["record_type"] == "cell", "row"]))
    cells = records.loc[records["record_type"] == "cell"]
    result = pd.DataFrame(index=row_order)
    for record in columns.itertuples(index=False):
        values = [
            cells.loc[(cells["row"] == row) & (cells["column"] == record.column), "value"].iloc[0]
            for row in row_order
        ]
        level_record = records.loc[
            (records["record_type"] == "levels") & (records["column"] == record.column),
            "value",
        ]
        if not level_record.empty:
            result[record.column] = pd.Categorical(
                values, categories=level_record.iloc[0].split("\x1e")
            )
        elif record.detail == "integer":
            result[record.column] = pd.Series(values, index=row_order, dtype="int64")
        elif record.detail == "double":
            result[record.column] = pd.Series(values, index=row_order, dtype="float64")
        elif record.detail == "logical":
            result[record.column] = pd.Series(
                [value == "TRUE" for value in values], index=row_order, dtype="bool"
            )
        else:
            raise AssertionError(f"unsupported R column type: {record.detail}")
    return result.loc[:, columns["column"].tolist()].reset_index(drop=True)


def test_seqsammseq_transition_and_empty_routes_match_r():
    ragged_x, covar = _ragged_fixture()
    ragged = traminex.seqsamm(ragged_x, 3, covar=covar)
    repeated_x, repeat_covar = _repeated_fixture()
    repeated = traminex.seqsamm(repeated_x, 2, covar=repeat_covar)

    cases = {
        "alpha": traminex.seqsammseq(ragged, "Alpha"),
        "gamma_empty": traminex.seqsammseq(ragged, "Gamma"),
        "unknown_empty": traminex.seqsammseq(ragged, "Unknown"),
        "repeated_a": traminex.seqsammseq(repeated, "State A"),
    }

    for case_name, result in cases.items():
        _assert_consumer_sequence(result, case_name)
    assert cases["alpha"].nseq == int(
        ((ragged.frame["s.1"] == "Alpha") & ragged.frame["transition"]).sum()
    )
    assert cases["gamma_empty"].data.shape == (0, 3)
    assert cases["unknown_empty"].data.shape == (0, 3)


def test_seqsammeha_all_r_routes_match_exact_frames():
    ragged_x, covar = _ragged_fixture()
    ragged = traminex.seqsamm(ragged_x, 3, covar=covar)
    repeated_x, repeat_covar = _repeated_fixture()
    repeated = traminex.seqsamm(repeated_x, 2, covar=repeat_covar)
    alpha_types = ["late", "early", "late"]
    factor_types = pd.Categorical(alpha_types, categories=["late", "early", "unused"], ordered=True)
    mixed_types = ["a", "A", "a"]
    repeat_types = ["route2", "route1", "route2", "route3"]

    cases = {
        "alpha_character": traminex.seqsammeha(ragged, "Alpha", alpha_types),
        "alpha_factor": traminex.seqsammeha(ragged, "Alpha", factor_types),
        "alpha_spell": traminex.seqsammeha(ragged, "Alpha", alpha_types, persper=False),
        "alpha_mixed_case": traminex.seqsammeha(ragged, "Alpha", mixed_types),
        "gamma_empty": traminex.seqsammeha(
            ragged,
            "Gamma",
            pd.Categorical([], categories=["unused", "future"], ordered=True),
        ),
        "repeated_person": traminex.seqsammeha(repeated, "State A", repeat_types),
        "repeated_spell": traminex.seqsammeha(repeated, "State A", repeat_types, persper=False),
    }

    for case_name, result in cases.items():
        pd.testing.assert_frame_equal(
            result, _golden_eha_frame(case_name), check_exact=False, rtol=0, atol=R_GOLDEN_ATOL
        )
    assert list(cases["alpha_factor"].columns[-3:]) == [
        "SAMMlate",
        "SAMMearly",
        "SAMMunused",
    ]
    assert list(cases["alpha_mixed_case"]["SAMMtypology"].cat.categories) == [
        "None",
        "a",
        "A",
    ]
    assert cases["repeated_spell"].groupby(["id", "begin"], observed=False).size().eq(1).all()


def test_seqsammeha_locale_order_matches_r():
    """Typology levels take R's collation (TX-AUD-20261011-015); they used a portable key."""
    x, covar = _ragged_fixture()
    samm = traminex.seqsamm(x, 3, covar=covar)

    result = traminex.seqsammeha(samm, "Alpha", ["a", "á", "b"])
    expected = _golden_eha_frame("alpha_locale")

    assert list(expected["SAMMtypology"].cat.categories) == ["None", "a", "á", "b"]
    assert list(result["SAMMtypology"].cat.categories) == ["None", "a", "á", "b"]


def test_seqsammeha_preserves_samm_metadata_and_frames():
    x, covar = _ragged_fixture()
    samm = traminex.seqsamm(x, 3, covar=covar)
    frame_snapshot = samm.frame.copy(deep=True)
    data_snapshot = samm.sequences.data.copy()
    typology_snapshot = samm.typology.copy()
    state_columns_snapshot = samm.state_columns

    traminex.seqsammeha(samm, "Alpha", ["late", "early", "late"])
    traminex.seqsammseq(samm, "Alpha")

    pd.testing.assert_frame_equal(
        samm.frame, frame_snapshot, check_exact=False, rtol=0, atol=R_GOLDEN_ATOL
    )
    np.testing.assert_array_equal(samm.sequences.data, data_snapshot)
    np.testing.assert_array_equal(samm.typology, typology_snapshot)
    assert samm.state_columns == state_columns_snapshot


def test_samm_consumers_validate_lengths_and_types():
    transcript = (R_GOLDENS_DIR / "samm_r3_errors.txt").read_text(encoding="utf-8")
    assert "seqsammeha_typology_short|warning|" in transcript
    assert "seqsammeha_typology_long|warning|" in transcript
    assert "seqsammseq_not_samm|success" in transcript
    x, covar = _ragged_fixture()
    samm = traminex.seqsamm(x, 3, covar=covar)

    with pytest.raises(ValueError, match="typology"):
        traminex.seqsammeha(samm, "Alpha", ["one", "two"])
    with pytest.raises(ValueError, match="typology"):
        traminex.seqsammeha(samm, "Alpha", ["one", "two", "three", "four"])
    with pytest.raises(ValueError, match="spell"):
        traminex.seqsammeha(samm, "Unknown", [])
    with pytest.raises(TypeError, match="SAMMResult"):
        traminex.seqsammseq(samm.frame, "Alpha")
    with pytest.raises(TypeError, match="SAMMResult"):
        traminex.seqsammeha(samm.frame, "Alpha", [])
    with pytest.raises(TypeError, match="spell"):
        traminex.seqsammseq(samm, 1)
    with pytest.raises(TypeError, match="persper"):
        traminex.seqsammeha(samm, "Alpha", ["one", "two", "three"], persper=1)


def test_samm_consumers_are_public_with_frozen_signatures():
    import inspect

    assert traminex.seqsammseq.__module__ == "traminex.samm"
    assert traminex.seqsammeha.__module__ == "traminex.samm"
    assert str(inspect.signature(traminex.seqsammseq)) == (
        "(samm: 'SAMMResult', spell: 'str') -> 'STSList'"
    )
    assert str(inspect.signature(traminex.seqsammeha)) == (
        "(samm: 'SAMMResult', spell: 'str', typology: 'Sequence[object]', "
        "persper: 'bool' = True) -> 'pd.DataFrame'"
    )
