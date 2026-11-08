"""Public TraMineRextras 0.6.9 evidence for ``rowmode``."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from traminex import rowmode

GOLDEN = Path(__file__).resolve().parent / "r_goldens" / "rowmode_public_cases.csv"


def _expected(case: str) -> pd.Series:
    rows = pd.read_csv(GOLDEN, dtype=str, keep_default_na=False)
    selected = rows.loc[rows["case"] == case]
    assert len(selected) == 1
    return selected.iloc[0]


@pytest.mark.parametrize(
    ("case", "values", "except_"),
    [
        ("numeric_tie", [2, 1, 2, 1], None),
        ("numeric_decimal", [2.5, 1.2, 2.5], None),
        ("logical_tie", [True, False, True, False], None),
        ("missing_character", [None, "B", "B", "A"], None),
        ("scalar_except", ["AA", "AA", "A"], "AA"),
        ("multiple_except", ["A", "B", "C", "C"], np.array(["C", "B"])),
        ("character_locale_tie", ["a", "B", "a", "B"], None),
        ("character_accent_tie", ["á", "ä", "á", "ä"], None),
    ],
)
def test_rowmode_scalar_values_match_public_r(case: str, values: object, except_: object) -> None:
    expected = _expected(case)
    actual = rowmode(values, except_=except_)
    assert isinstance(actual, str)
    assert actual == expected["value"]
    assert expected["r_type"] == "character"
    assert expected["r_length"] == "1"


def test_rowmode_pandas_categorical_preserves_factor_level_tie_order() -> None:
    values = pd.Categorical(["A", "B"], categories=["B", "A"])
    expected = _expected("factor_tie")
    assert rowmode(values) == expected["value"] == "B"


@pytest.mark.parametrize(
    ("case", "values", "except_"),
    [
        ("all_excluded", ["A", "A"], "A"),
        ("empty_character", [], None),
        ("all_na_character", [None, np.nan], None),
    ],
)
def test_rowmode_no_result_is_stably_reexpressed_as_none(
    case: str, values: object, except_: object
) -> None:
    expected = _expected(case)
    assert expected["r_length"] == "0"
    assert rowmode(values, except_=except_) is None


def test_rowmode_retains_explicit_row_wise_python_extension() -> None:
    values = np.array([["B", "A", "A"], ["B", "B", "A"]], dtype=object)
    actual = rowmode(values)
    assert isinstance(actual, np.ndarray)
    assert actual.tolist() == ["A", "B"]
    assert _expected("matrix_flattened")["value"] == "A"


def test_rowmode_rejects_more_than_two_dimensions() -> None:
    with pytest.raises(ValueError, match="1D or 2D"):
        rowmode(np.zeros((1, 1, 1)))
