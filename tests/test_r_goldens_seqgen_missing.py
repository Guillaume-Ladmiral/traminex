"""R-derived coverage for TraMineRextras 0.6.9 ``seqgen.missing``."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDENS_DIR

import traminex
from traminex import sequences
from traminex.sequences import STSList


def _rich_fixture() -> STSList:
    return STSList(
        data=np.array(
            [
                [0, 1, 0, 1, 0, 1],
                [1, 0, 1, 0, 3, 3],
                [0, 1, 0, 3, 3, 3],
                [1, 0, 1, 0, 1, 3],
                [0, 1, 0, 1, 0, 1],
                [1, 0, 1, 0, 1, 0],
            ],
            dtype=int,
        ),
        alphabet=["A", "B", "*", "%"],
        labels=["Alpha", "Beta", "Missing", "Void"],
        cpal=["#111111", "#222222"],
        xtlab=["2000", "2001", "2002", "2003", "2004", "2005"],
        xtstep=2,
        tick_last=True,
        nr_code=2,
        void_code=3,
        weights=np.array([1, 2, 3, 4, 5, 6], dtype=float),
        start=2000,
        ids=["a", "b", "c", "d", "e", "f"],
        missing_color="#999999",
        issues={"empty": [3]},
    )


def _reference(
    x: STSList,
    *,
    p_cases: float,
    p_left: float,
    p_gaps: float,
    p_right: float,
    mt_left: str,
    mt_gaps: str,
    mt_right: str,
    seed: int,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    result = x.data.copy()
    length = int(np.max(np.sum(result != x.void_code, axis=1), initial=0))
    selected_count = int(np.rint(p_cases * x.nseq))
    selected = np.sort(rng.choice(x.nseq, size=selected_count, replace=False))
    left_draws = rng.uniform(0, p_left, size=x.nseq)
    gap_draws = rng.uniform(0, p_gaps, size=x.nseq)
    right_draws = rng.uniform(0, p_right, size=x.nseq)
    code: Mapping[str, int] = {"nr": int(x.nr_code), "void": int(x.void_code)}

    for row in selected:
        gap_count = int(np.rint(gap_draws[row] * length))
        gaps = rng.choice(length, size=gap_count, replace=False)
        result[row, gaps] = code[mt_gaps]
        left_count = int(np.rint(left_draws[row] * length))
        if left_count:
            result[row, :left_count] = code[mt_left]
        right_count = int(np.rint(right_draws[row] * length))
        if right_count:
            result[row, length - right_count : length] = code[mt_right]
    return result


def _metadata(x: STSList) -> tuple[object, ...]:
    return (
        tuple(x.alphabet),
        None if x.labels is None else tuple(x.labels),
        None if x.cpal is None else tuple(x.cpal),
        None if x.xtlab is None else tuple(x.xtlab),
        x.xtstep,
        x.tick_last,
        x.nr,
        x.void,
        x.nr_code,
        x.void_code,
        None if x.weights is None else tuple(x.weights),
        x.start,
        None if x.ids is None else tuple(x.ids),
        x.missing_color,
        x.issues,
    )


def test_r_manifest_freezes_rounding_selection_and_overlap_execution() -> None:
    manifest = pd.read_csv(
        R_GOLDENS_DIR / "seqgen_missing_r3_manifest.csv",
        dtype=str,
        keep_default_na=False,
    )
    rounded = manifest.loc[manifest["case"] == "rounded_case_counts"]
    assert dict(zip(rounded["index"], rounded["value"], strict=True)) == {
        "0": "0",
        "0.1": "0",
        "0.3": "2",
        "0.5": "2",
        "1": "5",
    }
    summary = manifest.loc[
        (manifest["case"] == "seeded_overlap") & (manifest["record_type"] == "summary")
    ]
    facts = dict(zip(summary["index"], summary["value"], strict=True))
    assert facts == {
        "n": "6",
        "length": "6",
        "rounded_cases": "4",
        "selected_ids": "1,2,3,6",
        "assignment_order": "gaps,left,right",
        "overwrite_precedence": "right>left>gaps",
    }
    output = manifest.loc[
        (manifest["case"] == "seeded_overlap") & (manifest["record_type"] == "output")
    ].set_index("index")["value"]
    assert output["1,1"] == "%"  # left overwrites an earlier gap
    assert output["1,4"] == "*"  # right overwrites an earlier gap
    assert output["2,6"] == "*"  # an unoverlapped gap remains nonresponse


def test_seeded_algorithm_matches_r_derived_deterministic_mechanics() -> None:
    x = _rich_fixture()
    before = x.data.copy()
    kwargs = {
        "p_cases": 2 / 3,
        "p_left": 0.9,
        "p_gaps": 0.9,
        "p_right": 0.9,
        "mt_left": "void",
        "mt_gaps": "nr",
        "mt_right": "nr",
    }

    result = sequences.seqgen_missing(x, **kwargs, random_state=1407)

    np.testing.assert_array_equal(result.data, _reference(x, **kwargs, seed=1407))
    np.testing.assert_array_equal(x.data, before)
    assert _metadata(result) == _metadata(x)
    assert result is not x
    assert result.data is not x.data


def test_seed_and_generator_routes_are_reproducible_and_isolated() -> None:
    x = _rich_fixture()
    kwargs = dict(p_cases=1, p_left=0.7, p_gaps=0.5, p_right=0.6)
    first = sequences.seqgen_missing(x, **kwargs, random_state=88)
    second = sequences.seqgen_missing(x, **kwargs, random_state=88)
    np.testing.assert_array_equal(first.data, second.data)

    generator_a = np.random.default_rng(91)
    generator_b = np.random.default_rng(91)
    a1 = sequences.seqgen_missing(x, **kwargs, random_state=generator_a)
    a2 = sequences.seqgen_missing(x, **kwargs, random_state=generator_a)
    b1 = sequences.seqgen_missing(x, **kwargs, random_state=generator_b)
    b2 = sequences.seqgen_missing(x, **kwargs, random_state=generator_b)
    np.testing.assert_array_equal(a1.data, b1.data)
    np.testing.assert_array_equal(a2.data, b2.data)

    np.random.seed(9321)
    expected = np.random.random(5)
    np.random.seed(9321)
    sequences.seqgen_missing(x, **kwargs, random_state=94)
    np.testing.assert_array_equal(np.random.random(5), expected)


def test_probability_boundaries_zero_width_and_ragged_length() -> None:
    x = _rich_fixture()
    unchanged = sequences.seqgen_missing(
        x,
        p_cases=1,
        p_left=0,
        p_gaps=0,
        p_right=0,
        random_state=1,
    )
    np.testing.assert_array_equal(unchanged.data, x.data)
    no_cases = sequences.seqgen_missing(x, p_cases=0, random_state=1)
    np.testing.assert_array_equal(no_cases.data, x.data)

    # max(seqlength) is 6 despite ragged rows, as frozen by the R manifest.
    result = sequences.seqgen_missing(
        x,
        p_cases=1,
        p_left=1,
        p_gaps=1,
        p_right=1,
        random_state=12,
    )
    assert result.data.shape == x.data.shape

    empty = STSList(data=np.empty((3, 0), dtype=int), alphabet=["A", "B"])
    empty_result = sequences.seqgen_missing(
        empty,
        p_cases=1,
        p_left=1,
        p_gaps=1,
        p_right=1,
        random_state=4,
    )
    assert empty_result.data.shape == (3, 0)


@pytest.mark.parametrize("name", ["p_cases", "p_left", "p_gaps", "p_right"])
@pytest.mark.parametrize("value", [-0.01, 1.01, np.nan, np.inf])
def test_probability_validation_is_bounded_and_finite(name: str, value: float) -> None:
    transcript = (R_GOLDENS_DIR / "seqgen_missing_r3_conditions.txt").read_text(encoding="utf-8")
    assert "p_cases_above_one|simpleError|cannot take a sample larger" in transcript
    assert "p_left_below_zero|warning|NAs produced" in transcript
    with pytest.raises(ValueError, match=name):
        sequences.seqgen_missing(_rich_fixture(), **{name: value}, random_state=3)


@pytest.mark.parametrize("name", ["mt_left", "mt_gaps", "mt_right"])
def test_missing_type_validation_matches_r_error_boundary(name: str) -> None:
    transcript = (R_GOLDENS_DIR / "seqgen_missing_r3_conditions.txt").read_text(encoding="utf-8")
    assert "invalid_missing_type|simpleError|replacement has length zero" in transcript
    with pytest.raises(ValueError, match=name):
        sequences.seqgen_missing(_rich_fixture(), **{name: "bad"}, random_state=3)


def test_missing_code_is_added_without_reclassifying_real_states() -> None:
    x = STSList(data=np.zeros((4, 8), dtype=int), alphabet=["A"])
    result = sequences.seqgen_missing(
        x,
        p_cases=1,
        p_left=1,
        p_right=1,
        random_state=7,
    )
    assert x.alphabet == ["A"]
    assert x.nr_code is None
    assert result.nr_code == 1
    assert result.alphabet == ["A", "*"]
    assert np.any(result.data == result.nr_code)


def test_seqgen_missing_is_the_public_r_derived_name() -> None:
    assert traminex.seqgen_missing is sequences.seqgen_missing
