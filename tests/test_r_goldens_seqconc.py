"""Complete public TraMineR 2.2-14 evidence for ``seqconc``."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from corestats_fixtures import fixture_a, fixture_b, fixture_c
from r_goldens_utils import R_GOLDEN_ATOL

from traminex import seqconc, seqfind, seqpm, seqtab

GOLDEN_DIR = Path(__file__).resolve().parent / "r_goldens"


def _golden(case: str) -> pd.DataFrame:
    rows = pd.read_csv(
        GOLDEN_DIR / "seqconc_public_cases.csv",
        dtype=str,
        keep_default_na=False,
    )
    selected = rows.loc[rows["case"] == case].reset_index(drop=True)
    assert not selected.empty
    return selected


def _assert_result(case: str, actual: np.ndarray) -> None:
    expected = _golden(case)
    assert isinstance(actual, np.ndarray)
    assert actual.ndim == 1
    assert actual.tolist() == expected["value"].tolist()
    assert tuple(actual.row_names) == tuple(expected["row_name"])
    assert actual.vname == expected["column_name"].iloc[0]


@pytest.mark.parametrize(
    ("case", "builder"),
    [
        ("sts_A", fixture_a),
        ("sts_B", fixture_b),
        ("sts_C", fixture_c),
    ],
)
def test_seqconc_stslist_uses_object_void_and_preserves_array_behavior(case: str, builder) -> None:
    x = builder()
    before = x.copy()
    actual = seqconc(x, void="ignored")
    _assert_result(case, actual)
    assert list(actual) == actual.tolist()
    np.testing.assert_array_equal(x.data, before.data)
    assert x.ids == before.ids
    assert x.alphabet == before.alphabet


def test_seqconc_raw_vector_default_void_and_fallback_rows_match_r() -> None:
    _assert_result("raw_vector", seqconc(["A", None, "B"]))


def test_seqconc_raw_matrix_literal_separator_name_and_void_match_r() -> None:
    data = np.array(
        [["A", None, "VOID"], ["B", "C", "VOID"]],
        dtype=object,
    )
    before = data.copy()
    actual = seqconc(data, sep=":", vname="Joined", void="VOID")
    _assert_result("raw_matrix_void", actual)
    np.testing.assert_array_equal(data, before)


@pytest.mark.parametrize(
    ("case", "selector"),
    [
        ("dataframe_named_var", ["third", "first"]),
        ("dataframe_integer_var", [3, 1]),
    ],
)
def test_seqconc_dataframe_selectors_preserve_requested_order_and_index(
    case: str, selector
) -> None:
    data = pd.DataFrame(
        {
            "first": ["A", "B"],
            "second": ["C", "D"],
            "third": ["E", "F"],
        },
        index=["id-1", "id-2"],
    )
    before = data.copy(deep=True)
    kwargs = {"sep": "|", "var": selector}
    if case == "dataframe_named_var":
        kwargs["vname"] = "Selected"
    _assert_result(case, seqconc(data, **kwargs))
    pd.testing.assert_frame_equal(data, before, check_exact=False, rtol=0, atol=R_GOLDEN_ATOL)


def test_seqconc_missing_only_one_row_zero_columns_and_zero_selector_match_r() -> None:
    _assert_result("missing_only", seqconc(np.array([[None]], dtype=object)))
    _assert_result(
        "one_row",
        seqconc(
            pd.DataFrame(
                [["A", None, "B"]],
                index=["only"],
                columns=["x", "y", "z"],
            )
        ),
    )
    _assert_result("zero_columns", seqconc(np.empty((2, 0), dtype=object)))
    frame = pd.DataFrame({"first": ["A", "B"], "second": ["C", "D"]}, index=["id-1", "id-2"])
    _assert_result("zero_selector", seqconc(frame, var=0))


def test_seqconc_public_errors_match_r_boundaries() -> None:
    errors = pd.read_csv(
        GOLDEN_DIR / "seqconc_public_errors.txt",
        sep="|",
        dtype=str,
        keep_default_na=False,
    )
    assert set(errors["case"]) == {"explicit_null_void", "unknown_column", "vname_length"}
    data = pd.DataFrame({"first": ["A"], "second": ["B"]})
    with pytest.raises((TypeError, ValueError), match="void|length zero"):
        seqconc(data, void=None)
    with pytest.raises(ValueError, match="undefined columns"):
        seqconc(data, var="missing")
    with pytest.raises((TypeError, ValueError), match="vname|length"):
        seqconc(data, vname=["A", "B"])


def test_seqconc_corrected_result_remains_compatible_with_direct_consumers() -> None:
    x = fixture_a()
    assert seqfind(x, x).tolist() == [0, 1, 2, 3]
    assert not seqtab(x, idxs=0, format="STS").empty
    assert seqpm(x, "AA", sep="")["MTab"]["nbocc"].iloc[0] == 1
