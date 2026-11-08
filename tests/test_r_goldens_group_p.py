"""R-golden coverage for TraMineRextras 0.6.9 ``group.p``."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDENS_DIR

import traminex
from traminex import stats_utils


def _inputs() -> dict[str, tuple[Sequence[object], Sequence[float] | None]]:
    return {
        "character_unit": (["b", "a", "b", "c"], None),
        "numeric_integer": ([10, 2, 10, 3], [2, 1, 7, 2]),
        "factor_unused_na": (
            pd.Categorical(["b", "a", "b", None], categories=["z", "b", "a", "c"]),
            None,
        ),
        "fractional": (
            pd.Categorical(["b", "a", "b", "c"], categories=["c", "b", "a", "z"]),
            [0.1, 0.2, 0.3, 0.4],
        ),
        "third_rounding": (["x", "y", "z"], [1, 1, 1]),
        "zero_total": (
            pd.Categorical(["b", "a"], categories=["b", "a", "z"]),
            [0, 0],
        ),
        "na_group": (["b", None, "a"], [1, 2, 3]),
        "na_weight": (["b", "a", "b"], [1, np.nan, 3]),
        "mixed_case": (["a", "A"], None),
        "signed_zero": (["a", "b"], [1, -1]),
        "infinite_weight": (["a", "b"], [np.inf, 1]),
    }


def _golden() -> pd.DataFrame:
    return pd.read_csv(
        R_GOLDENS_DIR / "group_p_r3_cases.csv",
        dtype=str,
        keep_default_na=False,
    )


def _snapshot(values: Sequence[object]) -> tuple[object, ...]:
    if isinstance(values, pd.Categorical):
        return (tuple(values.categories), tuple(values.codes), values.ordered)
    return tuple(values)


@pytest.mark.parametrize("case_name", list(_inputs()))
def test_group_p_values_levels_codes_and_nonmutation_match_r(case_name: str) -> None:
    group, weights = _inputs()[case_name]
    group_before = _snapshot(group)
    weights_before = None if weights is None else tuple(weights)

    result = stats_utils.group_p(group, weights)

    case = _golden().loc[lambda frame: frame["case"] == case_name]
    expected_values = case.loc[lambda frame: frame["record_type"] == "value"]
    expected_levels = case.loc[lambda frame: frame["record_type"] == "level"]
    assert isinstance(result, pd.Categorical)
    assert list(result.categories) == expected_levels["level"].tolist()
    assert result.codes.tolist() == [
        -1 if code == "NA" else int(code) - 1 for code in expected_values["code"]
    ]
    assert ["<NA>" if pd.isna(value) else str(value) for value in result] == (
        expected_values["value"].tolist()
    )
    assert _snapshot(group) == group_before
    if weights is not None:
        np.testing.assert_equal(np.asarray(weights), np.asarray(weights_before))


def test_group_p_rejects_weight_length_mismatch_recorded_by_r() -> None:
    transcript = (R_GOLDENS_DIR / "group_p_r3_errors.txt").read_text(encoding="utf-8")
    assert "weight_length|simpleError|variable lengths differ" in transcript
    with pytest.raises(ValueError, match="weights must have the same length as group"):
        stats_utils.group_p(["a", "b"], [1])


def test_group_p_locale_order_matches_r() -> None:
    """Text levels take R's collation (TX-AUD-20261011-015); they used a portable key."""
    labels = ["_", "-", "!", "?", ".", "*", "%", "a", "A", "á", "ä", "b"]
    case = _golden().loc[lambda frame: frame["case"] == "locale_labels"]
    r_levels = case.loc[lambda frame: frame["record_type"] == "level", "level"].tolist()
    assert r_levels == [f"{label} (8.3%)" for label in labels]

    result = stats_utils.group_p(labels[::-1])

    assert list(result.categories) == r_levels


def test_group_p_is_the_public_r_derived_name() -> None:
    assert traminex.group_p is stats_utils.group_p
