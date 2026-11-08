"""Complete public TraMineR 2.2-14 evidence for ``seqdecomp``."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDEN_ATOL

from traminex import seqconc, seqdecomp, seqdim

GOLDEN_DIR = Path(__file__).resolve().parent / "r_goldens"


def _golden(case: str) -> pd.DataFrame:
    rows = pd.read_csv(
        GOLDEN_DIR / "seqdecomp_public_cases.csv",
        dtype=str,
        keep_default_na=False,
    )
    selected = rows.loc[rows["case"] == case].reset_index(drop=True)
    assert not selected.empty
    return selected


def _expected_matrix(case: str) -> tuple[np.ndarray, tuple[str, ...], tuple[str, ...]]:
    expected = _golden(case)
    nrow = int(expected["nrow"].iloc[0])
    ncol = int(expected["ncol"].iloc[0])
    matrix = np.empty((nrow, ncol), dtype=object)
    for row in expected.itertuples(index=False):
        value = None if row.is_na == "TRUE" else row.value
        matrix[int(row.row_position) - 1, int(row.column_position) - 1] = value
    row_names = tuple(expected.drop_duplicates("row_position")["row_name"])
    column_names = tuple(expected.drop_duplicates("column_position")["column_name"])
    return matrix, row_names, column_names


def _assert_result(case: str, actual: np.ndarray) -> None:
    expected, row_names, column_names = _expected_matrix(case)
    assert isinstance(actual, np.ndarray)
    np.testing.assert_array_equal(actual, expected)
    assert tuple(actual.row_names) == row_names
    assert tuple(actual.column_names) == column_names


@pytest.mark.parametrize(
    ("case", "values", "kwargs"),
    [
        ("standard", ["A-B-C", "D-E"], {"sep": "-"}),
        ("colon", ["A:B", "C:D:E"], {"sep": ":"}),
        ("empty_separator", ["ABC", "DE"], {"sep": ""}),
        ("regex_pipe", ["|A||B|", "C|D"], {"sep": r"\|"}),
        ("regex_dot", ["ABC", "D"], {"sep": "."}),
        ("multiple_missing", ["A-X-NA", "B-C"], {"miss": ["X", "NA"]}),
        ("leading_separator", ["-A-B", "C-D"], {}),
        ("trailing_separator", ["A-B-", "C-D"], {}),
        ("repeated_separator", ["A--B", "C-D"], {}),
        ("empty_string_mixed", ["", "A-B"], {}),
    ],
)
def test_seqdecomp_tokenization_padding_and_missing_values_match_r(
    case: str, values: list[str], kwargs: dict[str, object]
) -> None:
    before = list(values)
    _assert_result(case, seqdecomp(values, **kwargs))
    assert values == before


def test_seqdecomp_dataframe_selector_vnames_and_index_match_r() -> None:
    data = pd.DataFrame(
        {"first": ["A-B", "C-D"], "second": ["E-F", "G-H"]},
        index=["r1", "r2"],
    )
    before = data.copy(deep=True)
    actual = seqdecomp(
        data,
        var="second",
        vnames=["left", "right"],
    )
    _assert_result("dataframe_selected", actual)
    pd.testing.assert_frame_equal(data, before, check_exact=False, rtol=0, atol=R_GOLDEN_ATOL)


@pytest.mark.parametrize(
    ("case", "selector"),
    [
        ("matrix_column_major", None),
        ("matrix_var_order", [2, 1]),
    ],
)
def test_seqdecomp_matrix_uses_r_column_major_extraction(
    case: str, selector: list[int] | None
) -> None:
    data = np.array([["A-B", "C-D"], ["E-F", "G-H"]], dtype=object)
    before = data.copy()
    _assert_result(case, seqdecomp(data, var=selector))
    np.testing.assert_array_equal(data, before)


def test_seqdecomp_public_errors_match_r_boundaries() -> None:
    errors = pd.read_csv(
        GOLDEN_DIR / "seqdecomp_public_errors.txt",
        sep="|",
        dtype=str,
        keep_default_na=False,
    )
    assert set(errors["case"]) == {
        "empty_string_only",
        "dataframe_multicolumn_rows",
        "vnames_length",
        "unknown_column",
    }
    data = pd.DataFrame(
        {"first": ["A-B", "C-D"], "second": ["E-F", "G-H"]},
        index=["r1", "r2"],
    )
    with pytest.raises(ValueError, match="dimnames|empty"):
        seqdecomp("")
    with pytest.raises(ValueError, match="dimnames|row"):
        seqdecomp(data)
    with pytest.raises(ValueError, match="vnames|length"):
        seqdecomp(["A-B", "C-D"], vnames=["only"])
    with pytest.raises(ValueError, match="undefined columns"):
        seqdecomp(data, var="missing")


def test_seqdecomp_preserves_seqconc_metadata_and_seqdim_integration() -> None:
    data = pd.DataFrame(
        [["A", "B"], ["C", "D"]],
        index=["person-1", "person-2"],
        columns=["t1", "t2"],
    )
    compressed = seqconc(data, vname="Trajectory")
    actual = seqdecomp(compressed, var="Trajectory")
    assert tuple(actual.row_names) == ("person-1", "person-2")
    np.testing.assert_array_equal(seqdim(compressed), [2, 2])
