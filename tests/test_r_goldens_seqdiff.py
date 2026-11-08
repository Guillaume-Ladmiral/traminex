"""R-golden tests for seqdiff (RP4 T16).

Anchors traminex's sliding-window discrepancy analysis to live TraMineR 2.2-14
``seqdiff`` (see scripts/r_goldens/make_goldens_seqdiff.R): per position window,
the Pseudo F / Pseudo Fbf / Pseudo R2 / Bartlett / Levene statistics and the
per-group discrepancy of an LCS distance matrix (norm="none").
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDEN_ATOL, load_golden

from traminex import seqdef, seqdiff

_ALPH = ["1", "2", "3", "4"]
_GOLDEN_DIR = Path(__file__).parent / "r_goldens"
_PUBLIC_CASES = pd.read_csv(_GOLDEN_DIR / "seqdiff_public_cases.csv")


def _fixture():
    codes = load_golden("seqdiff_input_codes").astype(int)
    w = load_golden("seqdiff_input_weights").ravel()
    chars = np.array([[_ALPH[c] for c in row] for row in codes], dtype=object)
    return seqdef(chars, alphabet=_ALPH, weights=w)


@pytest.mark.parametrize("tag,weighted", [("uw", False), ("w", True)])
def test_seqdiff_matches_live_r(tag, weighted):
    x = _fixture()
    g = load_golden("seqdiff_input_group").ravel().astype(int)
    out = seqdiff(x, g, method="LCS", norm=None, weighted=weighted)
    np.testing.assert_allclose(
        out["stat"].to_numpy(dtype=float),
        load_golden(f"seqdiff_{tag}_stat"),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )
    np.testing.assert_allclose(
        out["discrepancy"].to_numpy(dtype=float),
        load_golden(f"seqdiff_{tag}_disc"),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


def _assert_public_case(out, case):
    expected_case = _PUBLIC_CASES[_PUBLIC_CASES["case"] == case]
    assert not expected_case.empty
    for kind in ("stat", "discrepancy"):
        expected = expected_case[expected_case["kind"] == kind]
        row_positions = expected["row_position"].drop_duplicates().astype(int).tolist()
        column_labels = expected.sort_values("column_position")["column"].drop_duplicates().tolist()
        values = np.array(
            [
                [
                    expected.loc[
                        (expected["row_position"] == row) & (expected["column_position"] == column),
                        "value",
                    ].item()
                    for column in range(1, len(column_labels) + 1)
                ]
                for row in row_positions
            ],
            dtype=float,
        )
        # The golden records R's own row labels; compare against those rather
        # than against the position numbers, which was a Python-authored
        # expectation that let the two diverge unnoticed (cycle-6 §9).
        r_row_labels = expected["r_row_label"].drop_duplicates().tolist()
        assert r_row_labels == [f"T{position}" for position in row_positions]
        assert out[kind].index.tolist() == r_row_labels
        assert list(out[kind].columns) == column_labels
        np.testing.assert_allclose(
            out[kind].to_numpy(dtype=float),
            values,
            rtol=0,
            atol=R_GOLDEN_ATOL,
        )


def test_seqdiff_public_default_and_modern_map_match_live_r():
    x = _fixture()
    group = load_golden("seqdiff_input_group").ravel().astype(int)

    _assert_public_case(seqdiff(x, group), "default_lcs_auto")
    _assert_public_case(
        seqdiff(
            x,
            group,
            seqdist_args={"method": "LCS", "norm": "none"},
            weighted=False,
        ),
        "lcs_none_map",
    )


@pytest.mark.parametrize("substitution_keyword", ["sm", "subm"])
def test_seqdiff_public_custom_om_map_matches_live_r(substitution_keyword):
    x = _fixture()
    group = load_golden("seqdiff_input_group").ravel().astype(int)
    subm = pd.DataFrame(
        np.array(
            [
                [0.0, 0.5, 1.5, 2.0],
                [0.5, 0.0, 2.0, 1.0],
                [1.5, 2.0, 0.0, 0.75],
                [2.0, 1.0, 0.75, 0.0],
            ]
        ),
        index=_ALPH,
        columns=_ALPH,
    )
    original_subm = subm.copy(deep=True)

    out = seqdiff(
        x,
        group,
        seqdist_args={
            "method": "OM",
            "norm": "none",
            substitution_keyword: subm,
            "indel": 0.75,
        },
        weighted=False,
    )
    _assert_public_case(out, "om_custom")
    pd.testing.assert_frame_equal(
        subm, original_subm, check_exact=False, rtol=0, atol=R_GOLDEN_ATOL
    )


@pytest.mark.parametrize(
    "distance_weighted,case",
    [
        (True, "chi2_weighted_distance"),
        (False, "chi2_unweighted_distance"),
    ],
)
def test_seqdiff_public_chi2_distance_weight_control_matches_live_r(
    distance_weighted,
    case,
):
    x = _fixture()
    group = load_golden("seqdiff_input_group").ravel().astype(int)

    out = seqdiff(
        x,
        group,
        seqdist_args={
            "method": "CHI2",
            "step": 2,
            "overlap": True,
            "weighted": distance_weighted,
        },
        weighted=False,
    )
    _assert_public_case(out, case)


def test_seqdiff_public_ordinary_missing_group_matches_live_r():
    x = _fixture()
    group = load_golden("seqdiff_input_group").ravel().astype(float)
    group[[2, 7]] = np.nan

    _assert_public_case(seqdiff(x, group), "ordinary_group_missing")


def test_seqdiff_public_ordinary_character_group_order_matches_live_r():
    x = _fixture()
    group = np.array(["a", "B"] * 5, dtype=object)

    _assert_public_case(seqdiff(x, group), "ordinary_group_locale")


def test_seqdiff_public_sts_group_matches_live_r():
    x = _fixture()
    group = seqdef(
        np.array(
            [
                ["z", "*", "m", "z", "a", "m"],
                ["a", "a", "m", "%", "%", "%"],
                ["m", "z", "z", "a", "a", "m"],
                ["z", "m", "a", "a", "m", "z"],
                ["a", "z", "m", "z", "a", "m"],
                ["m", "a", "z", "m", "z", "a"],
                ["z", "m", "a", "z", "m", "a"],
                ["a", "z", "m", "a", "z", "m"],
                ["m", "a", "z", "m", "a", "z"],
                ["z", "a", "m", "z", "m", "a"],
            ],
            dtype=object,
        ),
        alphabet=["z", "a", "m"],
        nr="*",
        void="%",
    )

    _assert_public_case(
        seqdiff(
            x,
            group,
            seqdist_args={"method": "LCS", "norm": "none"},
        ),
        "sts_group",
    )


def test_seqdiff_public_missing_override_matches_live_r():
    matrix = np.array(
        [
            ["1", "2", "3", "4", "1", "2"],
            ["2", "1", "*", "3", "4", "1"],
            ["3", "3", "1", "2", "2", "3"],
            ["4", "1", "2", "*", "2", "1"],
            ["1", "1", "2", "3", "4", "2"],
            ["2", "3", "4", "1", "1", "3"],
            ["3", "2", "1", "4", "3", "%"],
            ["4", "4", "1", "1", "1", "%"],
            ["1", "3", "2", "4", "3", "%"],
            ["2", "2", "3", "3", "4", "%"],
            ["3", "1", "2", "2", "1", "%"],
            ["4", "2", "1", "3", "2", "%"],
        ],
        dtype=object,
    )
    weights = 0.5 + 0.1 * np.arange(1, len(matrix) + 1)
    x = seqdef(
        matrix,
        alphabet=_ALPH,
        nr="*",
        void="%",
        weights=weights,
    )
    group = np.array([1, 2] * 6)

    out = seqdiff(
        x,
        group,
        seqdist_args={
            "method": "LCS",
            "norm": "none",
            "with.missing": False,
        },
        with_missing=True,
    )
    _assert_public_case(out, "missing_override")


def test_seqdiff_public_deprecated_alias_and_conflict_evidence():
    x = _fixture()
    group = load_golden("seqdiff_input_group").ravel().astype(int)

    with pytest.deprecated_call(match="seqdist_arg"):
        out = seqdiff(
            x,
            group,
            seqdist_arg={"method": "LCS", "norm": "none"},
            weighted=False,
        )
    _assert_public_case(out, "deprecated_alias")

    warnings_text = (_GOLDEN_DIR / "seqdiff_public_warnings.txt").read_text()
    errors = pd.read_csv(
        _GOLDEN_DIR / "seqdiff_public_errors.txt",
        sep="|",
        dtype=str,
        keep_default_na=False,
    )
    assert "seqdist_arg is deprecated" in warnings_text
    assert set(errors["case"]) == {
        "both_maps",
        "empty_map",
        "cmprange_reversed",
        "cmprange_outside",
    }
    assert (
        "seqdist.args and seqdist_arg cannot be specified together"
        in errors.loc[errors["case"] == "both_maps", "message"].item()
    )

    with pytest.raises(TypeError, match="method"):
        seqdiff(x, group, seqdist_args={})


# --- cycle-6 §9 (S9): R labels the result rows with the sequence's time labels ---
def test_result_rows_carry_the_time_labels_r_uses():
    """``SOURCES/TraMineR/R/seqdiff.R:29,32`` sets

        rownames(ret$stat)        <- colnames(seqdata)[totrange]
        rownames(ret$discrepancy) <- colnames(seqdata)[totrange]

    so R labels each comparison window with the *time label* of its position.
    traminex used a bare 1..n integer index and dropped `xtlab` entirely, so a
    sequence object built with calendar column names lost them on the way out
    while the parity ledger recorded the row as EXACT.
    """
    years = ["1990", "1991", "1992", "1993", "1994", "1995"]
    sequences = seqdef(
        np.array(
            [
                ["A", "A", "B", "B", "C", "C"],
                ["A", "B", "B", "C", "C", "C"],
                ["B", "B", "B", "C", "C", "A"],
                ["A", "A", "A", "B", "B", "C"],
                ["C", "C", "B", "B", "A", "A"],
                ["B", "A", "B", "A", "B", "A"],
            ],
            dtype=object,
        ),
        alphabet=["A", "B", "C"],
        cnames=years,
    )
    group = np.array(["g1", "g1", "g1", "g2", "g2", "g2"])

    result = seqdiff(sequences, group)
    # cmprange=(0, 1) compares positions 1..5, so R labels them 1990..1994.
    assert list(result["stat"].index) == years[:5]
    assert list(result["discrepancy"].index) == years[:5]

    wider = seqdiff(sequences, group, cmprange=(0, 2))
    assert list(wider["stat"].index) == years[:4]
    assert list(wider["discrepancy"].index) == years[:4]


def test_result_rows_fall_back_to_r_default_time_labels():
    """`seqdef` without `cnames` gives R's default column names T1..TL."""
    sequences = seqdef(
        np.array(
            [
                ["A", "A", "B", "B"],
                ["A", "B", "B", "A"],
                ["B", "B", "A", "A"],
                ["B", "A", "A", "B"],
            ],
            dtype=object,
        ),
        alphabet=["A", "B"],
    )
    result = seqdiff(sequences, np.array(["g1", "g1", "g2", "g2"]))
    assert list(result["stat"].index) == ["T1", "T2", "T3"]
    assert list(result["discrepancy"].index) == ["T1", "T2", "T3"]
