"""Public TraMineR 2.2-14 evidence for the ``seqsubm`` wrapper."""

from __future__ import annotations

import copy
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from corestats_fixtures import fixture_a, fixture_b
from r_goldens_utils import R_GOLDEN_ATOL
from test_r_goldens_seqcost import STATIONARY_CASES, TIMEVARYING_CASES

from traminex.distances import seqsubm
from traminex.distances.costs import SeqCostMatrix, seqcost

GOLDEN_DIR = Path(__file__).resolve().parent / "r_goldens"
ATOL = 1e-10


def _features() -> pd.DataFrame:
    return pd.DataFrame({"constant": [5, 5, 5, 5], "varying": [0, 10, 0, 10]})


PUBLIC_CASES: dict[str, tuple[Any, dict[str, Any]]] = {
    "trate_unweighted": (fixture_a, {"method": "TRATE", "weighted": False}),
    "trate_weighted": (
        lambda: fixture_a(weighted=True),
        {"method": "TRATE", "weighted": True},
    ),
    "trate_time_missing": (
        lambda: fixture_b(weighted=True),
        {
            "method": "TRATE",
            "time_varying": True,
            "with_missing": True,
            "weighted": True,
            "transition": "both",
            "lag": 1,
            "miss_cost": [7, 8, 9],
            "miss_cost_fixed": True,
        },
    ),
    "features": (
        fixture_a,
        {
            "method": "FEATURES",
            "state_features": _features(),
            "feature_weights": [2, 1],
            "feature_type": {"numeric": ["constant", "varying"]},
        },
    ),
}


def _expected(case: str) -> tuple[list[str], list[str], np.ndarray]:
    table = pd.read_csv(
        GOLDEN_DIR / "seqsubm_public_cases.csv",
        dtype={"case": str, "options": str, "state_i": str, "state_j": str, "time_label": str},
        keep_default_na=False,
    )
    rows = table.loc[table["case"] == case]
    assert not rows.empty
    states = rows["state_i"].drop_duplicates().tolist()
    assert rows["state_j"].drop_duplicates().tolist() == states
    times = rows.loc[rows["time_label"] != "", "time_label"].drop_duplicates().tolist()
    state_index = {state: index for index, state in enumerate(states)}
    if times:
        time_index = {label: index for index, label in enumerate(times)}
        expected = np.full((len(states), len(states), len(times)), np.nan)
        for row in rows.itertuples(index=False):
            expected[
                state_index[row.state_i],
                state_index[row.state_j],
                time_index[row.time_label],
            ] = float(row.value)
    else:
        expected = np.full((len(states), len(states)), np.nan)
        for row in rows.itertuples(index=False):
            expected[state_index[row.state_i], state_index[row.state_j]] = float(row.value)
    return states, times, expected


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ \ \ \ \ \ for\ example,\ cost\ is\ zero\ for\ the\ pair\ \(\ A\ ,\ C\ \):traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ \[!\]\ \ \d+\ substitution\ cost\ value\(s\)\ are\ \d+!:traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.parametrize("case", list(PUBLIC_CASES))
def test_public_seqsubm_matches_direct_r_golden(case: str):
    builder, kwargs = PUBLIC_CASES[case]
    actual = seqsubm(builder(), **kwargs)
    states, times, expected = _expected(case)

    if times:
        assert isinstance(actual, SeqCostMatrix)
        assert list(actual.from_states) == states
        assert list(actual.to_states) == states
        assert list(actual.time_labels) == times
    else:
        assert isinstance(actual, pd.DataFrame)
        assert actual.index.tolist() == states
        assert actual.columns.tolist() == states
    np.testing.assert_allclose(np.asarray(actual), expected, rtol=0, atol=ATOL)


@pytest.mark.parametrize(
    ("case", "spec"),
    [*STATIONARY_CASES.items(), *TIMEVARYING_CASES.items()],
)
def test_public_seqsubm_is_seqcost_sm_for_every_cost_case(case: str, spec):
    builder, kwargs = spec
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        expected = seqcost(builder(), **kwargs)["sm"]
        actual = seqsubm(builder(), **kwargs)

    assert type(actual) is type(expected)
    np.testing.assert_allclose(
        np.asarray(actual),
        np.asarray(expected),
        rtol=0,
        atol=ATOL,
        equal_nan=True,
        err_msg=case,
    )
    if isinstance(expected, pd.DataFrame):
        assert actual.index.equals(expected.index)
        assert actual.columns.equals(expected.columns)
    else:
        assert actual.from_states == expected.from_states
        assert actual.to_states == expected.to_states
        assert actual.time_labels == expected.time_labels


def test_public_seqsubm_proximities_returns_r_null_analogue():
    outcomes = pd.read_csv(
        GOLDEN_DIR / "seqsubm_public_outcomes.txt",
        sep="|",
        dtype=str,
        keep_default_na=False,
    )
    assert outcomes.to_dict("records") == [
        {
            "case": "proximities_null",
            "call": "TraMineR::seqsubm(fx$A,method=TRATE,proximities=TRUE)",
            "result": "NULL",
        }
    ]
    assert seqsubm(fixture_a(), method="TRATE", proximities=True) is None


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ At\ least\ one\ substitution\ cost\ greater\ than\ twice\ the\ indel\ cost\.:traminex.distances.python_kernels.CostDiagnosticMessage"
)
def test_public_seqsubm_requires_method_and_accepts_positional_cval():
    with pytest.raises(TypeError, match="method"):
        seqsubm(fixture_a())
    actual = seqsubm(fixture_a(), "CONSTANT", 3.5)
    expected = seqcost(fixture_a(), "CONSTANT", 3.5)["sm"]
    pd.testing.assert_frame_equal(actual, expected, check_exact=False, rtol=0, atol=R_GOLDEN_ATOL)


def test_public_seqsubm_propagates_seqcost_validation():
    x = fixture_a()
    with pytest.raises(ValueError) as seqcost_error:
        seqcost(x, "constant")
    with pytest.raises(
        type(seqcost_error.value), match=str(seqcost_error.value).replace("[", r"\[")
    ):
        seqsubm(x, "constant")


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ At\ least\ one\ substitution\ cost\ greater\ than\ twice\ the\ indel\ cost\.:traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ \ at\ least,\ one\ substitution\ cost\ doesn't\ respect\ the\ triangle\ inequality\.:traminex.distances.python_kernels.CostDiagnosticMessage"
)
def test_public_and_raw_seqsubm_inputs_are_not_mutated():
    x = fixture_b(weighted=True)
    data_before = x.data.copy()
    weights_before = x.weights.copy()
    features = pd.DataFrame({"score": [0.0, 0.5, 1.0, 0.25]})
    features_before = features.copy(deep=True)
    options = {
        "method": "FEATURES",
        "with_missing": True,
        "state_features": features,
        "feature_weights": [2.0],
        "feature_type": {"numeric": ["score"]},
        "miss_cost": [7.0, 8.0, 9.0],
        "miss_cost_fixed": True,
    }
    options_before = copy.deepcopy(options)

    seqsubm(x, **options)

    np.testing.assert_array_equal(x.data, data_before)
    np.testing.assert_array_equal(x.weights, weights_before)
    pd.testing.assert_frame_equal(
        features, features_before, check_exact=False, rtol=0, atol=R_GOLDEN_ATOL
    )
    assert options.keys() == options_before.keys()
    for key in options.keys() - {"state_features"}:
        assert options[key] == options_before[key]
    pd.testing.assert_frame_equal(
        options["state_features"],
        options_before["state_features"],
        check_exact=False,
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )

    raw = np.array([[0, 1, 2], [2, 1, 0]], dtype=np.int32)
    raw_before = raw.copy()
    seqsubm(raw)
    np.testing.assert_array_equal(raw, raw_before)
