"""R-golden parity for R2 seqfpos and seqfposend single-state branches."""

from __future__ import annotations

import csv
import re

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDEN_ATOL, R_GOLDENS_DIR, load_golden_lines

from traminex import seqdef, seqfpos, seqfposend

RAW = np.asarray(
    [
        ["A", "A", "B", "C", "C"],
        ["B", "A", "A", "C", "B"],
        ["B", "B", "C", "C", "C"],
        ["*", "A", "*", "B", "B"],
        ["A", "B", "*", "*", "*"],
        ["*", "*", "*", "*", "*"],
    ],
    dtype=object,
)


def _fixture():
    return seqdef(RAW, alphabet=["A", "B", "C", "D"], missing="*")


def _rows():
    with (R_GOLDENS_DIR / "seqfpos_r2_cases.csv").open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _values(rows):
    return np.asarray(
        [float("nan") if row["value"] == "NA" else float(row["value"]) for row in rows]
    )


def test_seqfpos_single_state_matrix_matches_r():
    rows = [row for row in _rows() if row["function"] == "seqfpos"]
    cases = sorted({row["case"] for row in rows})
    assert cases == ["raw_A", "raw_B", "raw_C", "raw_D", "raw_Z", "raw_missing"]

    x = _fixture()
    for case in cases:
        expected_rows = [row for row in rows if row["case"] == case]
        state = expected_rows[0]["state"]
        actual = seqfpos(x, state)
        np.testing.assert_allclose(
            actual,
            _values(expected_rows),
            rtol=0,
            atol=R_GOLDEN_ATOL,
            equal_nan=True,
        )


def test_seqfpos_singular_missing_search_ignores_compatibility_flag():
    x = _fixture()
    np.testing.assert_allclose(
        seqfpos(x, "*", include_missing=False),
        seqfpos(x, "*", include_missing=True),
        equal_nan=True,
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


def test_seqfposend_single_state_matrix_matches_r():
    rows = [row for row in _rows() if row["function"] == "seqfposend"]
    cases = sorted({row["case"] for row in rows})
    assert len(cases) == 6 * 2 * 4 * 2
    pattern = re.compile(
        r"^spells_.+_wm(?P<wm>TRUE|FALSE)_lead(?P<lead>-?\d+)_origin(?P<origin>TRUE|FALSE)$"
    )

    x = _fixture()
    for case in cases:
        expected_rows = [row for row in rows if row["case"] == case]
        route = pattern.match(case)
        assert route is not None
        actual = seqfposend(
            x,
            expected_rows[0]["state"],
            include_missing=route["wm"] == "TRUE",
            lead=int(route["lead"]),
            from_seq_start=route["origin"] == "TRUE",
        )
        np.testing.assert_allclose(
            actual,
            _values(expected_rows),
            rtol=0,
            atol=R_GOLDEN_ATOL,
            equal_nan=True,
        )


def test_seqfpos_invalid_state_and_type_dispositions_are_frozen():
    assert load_golden_lines("seqfpos_r2_errors") == [
        "seqfpos_non_stslist=ERROR:data is not a sequence object, use 'seqdef' function to create one",
        "seqfposend_non_stslist=ERROR:data is not a sequence object, use 'seqdef' function to create one",
        "seqfpos_invalid_state=NA,NA,NA,NA,NA,NA",
        "seqfposend_invalid_state=0,0,0,0,0,0",
    ]
    assert np.isnan(seqfpos(_fixture(), "Z")).all()
    np.testing.assert_array_equal(seqfposend(_fixture(), "Z"), np.zeros(6))
    with pytest.raises(TypeError, match="state sequence object"):
        seqfpos(RAW, "A")
    with pytest.raises(TypeError, match="state sequence object"):
        seqfposend(RAW, "A")


def test_vectorized_dataframe_compatibility_and_argument_exclusion():
    x = _fixture()
    assert isinstance(seqfpos(x), pd.DataFrame)
    assert isinstance(seqfposend(x), pd.DataFrame)
    with pytest.raises(ValueError, match="state and states"):
        seqfpos(x, "A", states=["B"])
    with pytest.raises(ValueError, match="state and states"):
        seqfposend(x, "A", states=["B"])
