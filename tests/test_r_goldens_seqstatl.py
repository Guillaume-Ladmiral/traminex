"""Complete public TraMineR 2.2-14 evidence for ``seqstatl``."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from corestats_fixtures import fixture_a, fixture_b, fixture_c
from r_goldens_utils import R_GOLDEN_ATOL

from traminex import seqstatl

GOLDEN_DIR = Path(__file__).resolve().parent / "r_goldens"


def _expected(case: str) -> list[str] | list[float]:
    rows = pd.read_csv(
        GOLDEN_DIR / "seqstatl_public_cases.csv",
        dtype=str,
        keep_default_na=False,
    )
    selected = rows.loc[rows["case"] == case].sort_values("position")
    assert not selected.empty
    assert selected["r_type"].nunique() == 1
    r_type = selected["r_type"].iloc[0]
    if r_type == "double":
        return [float(value) for value in selected["value"]]
    assert r_type == "character"
    return list(selected["value"])


@pytest.mark.parametrize(
    ("case", "factory"),
    [
        ("sts_A", fixture_a),
        ("sts_B", fixture_b),
        ("sts_C", fixture_c),
    ],
)
def test_seqstatl_stslist_markers_match_r(case: str, factory) -> None:
    data = factory()
    before_data = data.data.copy()
    before_alphabet = list(data.alphabet)

    assert seqstatl(data) == _expected(case)

    np.testing.assert_array_equal(data.data, before_data)
    assert data.alphabet == before_alphabet


def _raw_extended() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "t1": ["é", "A", None, "-"],
            "t2": ["Z", "ä", "10", "!"],
            "t3": ["2", "A", "2", "é"],
        }
    )


@pytest.mark.parametrize(
    ("case", "selector"),
    [
        ("raw_extended", None),
        ("named_var", ["t3", "t1"]),
        ("integer_var", [3, 1]),
    ],
)
def test_seqstatl_raw_dataframe_and_var_match_r(
    case: str,
    selector: list[str] | list[int] | None,
) -> None:
    data = _raw_extended()
    before = data.copy(deep=True)

    assert seqstatl(data, var=selector) == _expected(case)

    pd.testing.assert_frame_equal(data, before, check_exact=False, rtol=0, atol=R_GOLDEN_ATOL)


@pytest.mark.parametrize(
    ("case", "data"),
    [
        ("compressed_dash", ["A-B-C", "B-D"]),
        ("compressed_colon", ["A:é:!", "Z:ä:A"]),
    ],
)
def test_seqstatl_autodetects_compressed_inputs(case: str, data: list[str]) -> None:
    before = list(data)
    assert seqstatl(data) == _expected(case)
    assert data == before


def test_seqstatl_sps_uses_public_format_route() -> None:
    data = ["(A,2)-(B,1)", "(B,1)-(C,2)"]
    before = list(data)
    assert seqstatl(data, format="SPS") == _expected("sps")
    assert data == before


@pytest.mark.parametrize("format", ["STS", "SPELL", "DSS"])
def test_seqstatl_accepts_only_source_format_spellings(format: str) -> None:
    data = np.array([["B", "A"], ["C", "B"]], dtype=object)
    before = data.copy()
    assert seqstatl(data, format=format) == _expected(f"format_{format}")
    np.testing.assert_array_equal(data, before)


def test_seqstatl_rejects_lowercase_format() -> None:
    errors = pd.read_csv(
        GOLDEN_DIR / "seqstatl_public_errors.txt",
        sep="|",
        dtype=str,
        keep_default_na=False,
    )
    expected = errors.loc[errors["case"] == "invalid_lowercase_format"].iloc[0]
    assert "Input format must be one of" in expected["message"]
    with pytest.raises(ValueError, match="Input format must be one of"):
        seqstatl(np.array([["A"]], dtype=object), format="sts")


@pytest.mark.parametrize(
    ("case", "data"),
    [
        (
            "numeric_integer_text",
            np.array([["3", "1", "10"], ["2", "10", "1"]], dtype=object),
        ),
        (
            "numeric_actual",
            np.array([[3, 1, 10], [2, 10, 1]], dtype=int),
        ),
        (
            "numeric_decimal_text",
            np.array([["2.5", "1.2"], ["10", "2.5"]], dtype=object),
        ),
        (
            "mixed_numeric_text",
            np.array([["10", "A"], ["2", "!"]], dtype=object),
        ),
        (
            "raw_missing",
            np.array([[None, "A"], ["B", None]], dtype=object),
        ),
        (
            "raw_void_literal",
            np.array([["%", None], ["A", "%"]], dtype=object),
        ),
    ],
)
def test_seqstatl_preserves_r_output_type_and_character_order(
    case: str,
    data: np.ndarray,
) -> None:
    before = data.copy()
    actual = seqstatl(data)
    expected = _expected(case)

    assert actual == expected
    assert all(type(value) is type(want) for value, want in zip(actual, expected))
    np.testing.assert_array_equal(data, before)
