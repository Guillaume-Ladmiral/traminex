"""Full TraMineR 2.2-14 result-contract parity for ``seqcost``.

The oracle artifacts are produced by ``make_goldens_seqcost.R`` with explicit
row-major loops.  They cover all six methods, scalar/vector/matrix indels,
time-varying costs, fixed and derived missing costs, and proximity output.
"""

from __future__ import annotations

import warnings
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from corestats_fixtures import fixture_a, fixture_b, fixture_c

import traminex as tx
from traminex.api_manifest import PUBLIC_API
from traminex.distances import seqcost
from traminex.format import seqdef
from traminex.sequences import STSList

GOLDEN_DIR = Path(__file__).resolve().parent / "r_goldens"
ATOL = 1e-10

Builder = Callable[[], STSList]
Case = tuple[Builder, dict[str, Any]]


def _absent_states() -> STSList:
    return seqdef(np.full((3, 3), "A", dtype=object), alphabet=["A", "B", "C"])


def _one_state() -> STSList:
    return seqdef(np.full((2, 3), "A", dtype=object), alphabet=["A"])


def _short() -> STSList:
    return seqdef(
        np.array([["A", "B"], ["B", "A"], ["A", "A"]], dtype=object),
        alphabet=["A", "B"],
    )


def _custom_time() -> STSList:
    x = fixture_a()
    x.xtlab = [f"age{i}" for i in range(10, 16)]
    return x


def _features_mixed() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "interval": [0, 5, 10, 10],
            "category": pd.Categorical(["x", "y", "x", "z"]),
            "flag": [0, 1, 0, 1],
            "ordered_score": [10, 20, 40, 30],
            "ratio": [1, 10, 100, 1000],
        }
    )


def _features_constant() -> pd.DataFrame:
    return pd.DataFrame({"constant": [5, 5, 5, 5], "varying": [0, 10, 0, 10]})


def _features_missing() -> pd.DataFrame:
    return pd.DataFrame({"score": [0, 0.5, 1, 0.25]})


FEATURE_MIXED_KWARGS = {
    "state_features": _features_mixed(),
    "feature_weights": [2, 1, 3, 1, 0.5],
    "feature_type": {
        "asymm": "flag",
        "ordratio": "ordered_score",
        "logratio": "ratio",
    },
}


STATIONARY_CASES: dict[str, Case] = {
    "constant_default": (fixture_a, {"method": "CONSTANT"}),
    "constant_custom": (fixture_a, {"method": "CONSTANT", "cval": 3.5}),
    "constant_no_observed_missing": (
        fixture_a,
        {"method": "CONSTANT", "with_missing": True},
    ),
    "trate_unweighted_lag1": (
        fixture_a,
        {"method": "TRATE", "weighted": False, "lag": 1},
    ),
    "trate_weighted_lag1": (
        lambda: fixture_a(weighted=True),
        {"method": "TRATE", "weighted": True, "lag": 1},
    ),
    "trate_unweighted_lag2": (
        fixture_a,
        {"method": "TRATE", "weighted": False, "lag": 2},
    ),
    "trate_unweighted_lag_neg1": (
        fixture_a,
        {"method": "TRATE", "weighted": False, "lag": -1},
    ),
    "trate_missing_fixed": (
        lambda: fixture_b(weighted=True),
        {
            "method": "TRATE",
            "weighted": True,
            "with_missing": True,
            "cval": 3,
            "miss_cost": [7, 8, 9],
            "miss_cost_fixed": True,
        },
    ),
    "trate_missing_derived": (
        fixture_b,
        {
            "method": "TRATE",
            "weighted": False,
            "with_missing": True,
            "cval": 3,
            "miss_cost_fixed": False,
        },
    ),
    "future_unweighted_lag2": (
        fixture_a,
        {"method": "FUTURE", "weighted": False, "lag": 2},
    ),
    "future_weighted_lag1": (
        lambda: fixture_a(weighted=True),
        {"method": "FUTURE", "weighted": True, "lag": 1},
    ),
    "future_missing_fixed": (
        fixture_b,
        {
            "method": "FUTURE",
            "weighted": False,
            "with_missing": True,
            "miss_cost": 7,
            "miss_cost_fixed": True,
        },
    ),
    "future_missing_derived": (
        fixture_b,
        {
            "method": "FUTURE",
            "weighted": False,
            "with_missing": True,
            "miss_cost_fixed": False,
        },
    ),
    "features_mixed": (fixture_a, {"method": "FEATURES", **FEATURE_MIXED_KWARGS}),
    "features_constant_numeric": (
        fixture_a,
        {"method": "FEATURES", "state_features": _features_constant()},
    ),
    "features_missing_fixed": (
        fixture_b,
        {
            "method": "FEATURES",
            "with_missing": True,
            "state_features": _features_missing(),
            "miss_cost": 2,
            "miss_cost_fixed": True,
        },
    ),
    "features_missing_derived": (
        fixture_b,
        {
            "method": "FEATURES",
            "with_missing": True,
            "state_features": _features_missing(),
            "miss_cost_fixed": False,
        },
    ),
    "indels_unweighted": (fixture_a, {"method": "INDELS", "weighted": False}),
    "indels_weighted": (
        lambda: fixture_a(weighted=True),
        {"method": "INDELS", "weighted": True},
    ),
    "indels_missing_unweighted": (
        fixture_b,
        {"method": "INDELS", "weighted": False, "with_missing": True},
    ),
    "indels_missing_weighted": (
        lambda: fixture_b(weighted=True),
        {"method": "INDELS", "weighted": True, "with_missing": True},
    ),
    "indels_ragged_missing_weighted": (
        lambda: fixture_c(weighted=True),
        {"method": "INDELS", "weighted": True, "with_missing": True},
    ),
    "indels_missing_fixed": (
        fixture_b,
        {
            "method": "INDELS",
            "weighted": False,
            "with_missing": True,
            "miss_cost": [7, 8, 9],
            "miss_cost_fixed": True,
        },
    ),
    "indels_absent_states": (_absent_states, {"method": "INDELS"}),
    "indelslog_unweighted": (fixture_a, {"method": "INDELSLOG", "weighted": False}),
    "indelslog_weighted": (
        lambda: fixture_a(weighted=True),
        {"method": "INDELSLOG", "weighted": True},
    ),
    "indelslog_missing_weighted": (
        lambda: fixture_b(weighted=True),
        {"method": "INDELSLOG", "weighted": True, "with_missing": True},
    ),
    "indelslog_ragged_missing_weighted": (
        lambda: fixture_c(weighted=True),
        {"method": "INDELSLOG", "weighted": True, "with_missing": True},
    ),
    "indelslog_missing_fixed": (
        fixture_b,
        {
            "method": "INDELSLOG",
            "weighted": False,
            "with_missing": True,
            "miss_cost": [7, 8, 9],
            "miss_cost_fixed": True,
        },
    ),
    "indelslog_absent_states": (_absent_states, {"method": "INDELSLOG"}),
}


TIMEVARYING_CASES: dict[str, Case] = {
    "constant_time_default": (
        fixture_a,
        {"method": "CONSTANT", "time_varying": True},
    ),
    "constant_time_missing_vector": (
        fixture_b,
        {
            "method": "CONSTANT",
            "time_varying": True,
            "with_missing": True,
            "cval": 3,
            "miss_cost": [7, 8, 9],
        },
    ),
    "trate_time_previous": (
        fixture_a,
        {
            "method": "TRATE",
            "time_varying": True,
            "transition": "previous",
            "weighted": False,
        },
    ),
    "trate_time_next": (
        fixture_a,
        {
            "method": "TRATE",
            "time_varying": True,
            "transition": "next",
            "weighted": False,
        },
    ),
    "trate_time_both": (
        fixture_a,
        {
            "method": "TRATE",
            "time_varying": True,
            "transition": "both",
            "weighted": False,
        },
    ),
    "trate_time_both_weighted": (
        lambda: fixture_a(weighted=True),
        {
            "method": "TRATE",
            "time_varying": True,
            "transition": "both",
            "weighted": True,
        },
    ),
    "trate_time_missing_fixed": (
        fixture_b,
        {
            "method": "TRATE",
            "time_varying": True,
            "transition": "both",
            "weighted": False,
            "with_missing": True,
            "miss_cost": [7, 8, 9],
            "miss_cost_fixed": True,
        },
    ),
    "trate_time_missing_derived": (
        fixture_b,
        {
            "method": "TRATE",
            "time_varying": True,
            "transition": "both",
            "weighted": False,
            "with_missing": True,
            "miss_cost_fixed": False,
        },
    ),
    "trate_time_short_boundary": (
        _short,
        {
            "method": "TRATE",
            "time_varying": True,
            "transition": "both",
            "weighted": False,
        },
    ),
    "trate_time_custom_labels": (
        _custom_time,
        {
            "method": "TRATE",
            "time_varying": True,
            "transition": "both",
            "weighted": False,
        },
    ),
}

for method in ("INDELS", "INDELSLOG"):
    prefix = method.lower()
    TIMEVARYING_CASES.update(
        {
            f"{prefix}_time_unweighted_object": (
                fixture_a,
                {"method": method, "time_varying": True, "weighted": False},
            ),
            f"{prefix}_time_weighted_object_flag_false": (
                lambda: fixture_a(weighted=True),
                {"method": method, "time_varying": True, "weighted": False},
            ),
            f"{prefix}_time_weighted_object_flag_true": (
                lambda: fixture_a(weighted=True),
                {"method": method, "time_varying": True, "weighted": True},
            ),
            f"{prefix}_time_missing_weighted": (
                lambda: fixture_b(weighted=True),
                {
                    "method": method,
                    "time_varying": True,
                    "weighted": True,
                    "with_missing": True,
                },
            ),
            f"{prefix}_time_ragged_missing_weighted": (
                lambda: fixture_c(weighted=True),
                {
                    "method": method,
                    "time_varying": True,
                    "weighted": True,
                    "with_missing": True,
                },
            ),
            f"{prefix}_time_absent_states": (
                _absent_states,
                {"method": method, "time_varying": True},
            ),
        }
    )


PROXIMITY_CASES: dict[str, Case] = {
    "prox_constant_stationary": (
        fixture_a,
        {"method": "CONSTANT", "proximities": True},
    ),
    "prox_trate_weighted_stationary": (
        lambda: fixture_a(weighted=True),
        {"method": "TRATE", "weighted": True, "proximities": True},
    ),
    "prox_trate_missing_fixed_stationary": (
        fixture_b,
        {
            "method": "TRATE",
            "weighted": False,
            "with_missing": True,
            "miss_cost": 7,
            "proximities": True,
        },
    ),
    "prox_features_stationary": (
        fixture_a,
        {
            "method": "FEATURES",
            "state_features": _features_constant(),
            "proximities": True,
        },
    ),
    "prox_indelslog_stationary": (
        fixture_a,
        {"method": "INDELSLOG", "weighted": False, "proximities": True},
    ),
    "prox_trate_time_both": (
        fixture_a,
        {
            "method": "TRATE",
            "time_varying": True,
            "transition": "both",
            "weighted": False,
            "proximities": True,
        },
    ),
    "prox_one_state_nan": (
        _one_state,
        {"method": "CONSTANT", "proximities": True},
    ),
    "prox_indels_time_absent_nan": (
        _absent_states,
        {"method": "INDELS", "time_varying": True, "proximities": True},
    ),
}


def _load(name: str) -> pd.DataFrame:
    return pd.read_csv(GOLDEN_DIR / name, dtype=str, keep_default_na=False)


def _number(value: str) -> float:
    if value == "DBL_MAX":
        return np.finfo(float).max
    if value in {"NA", "NaN"}:
        return np.nan
    if value == "Inf":
        return np.inf
    if value == "-Inf":
        return -np.inf
    return float(value)


def _ordered(values: pd.Series) -> list[str]:
    return values.loc[values != ""].drop_duplicates().tolist()


def _expected_matrix(rows: pd.DataFrame, component: str) -> tuple[list[str], list[str], np.ndarray]:
    values = rows.loc[rows["component"] == component]
    states = _ordered(values["state_i"])
    time_labels = _ordered(values["time_label"])
    assert _ordered(values["state_j"]) == states
    expected_rows = len(states) * len(states) * (len(time_labels) or 1)
    assert len(values) == expected_rows
    assert not values[["state_i", "state_j", "time_label"]].duplicated().any()
    state_index = {state: i for i, state in enumerate(states)}
    if time_labels:
        time_index = {label: i for i, label in enumerate(time_labels)}
        expected = np.full((len(states), len(states), len(time_labels)), np.nan, dtype=float)
        for row in values.itertuples(index=False):
            expected[
                state_index[row.state_i],
                state_index[row.state_j],
                time_index[row.time_label],
            ] = _number(row.value)
    else:
        expected = np.full((len(states), len(states)), np.nan, dtype=float)
        for row in values.itertuples(index=False):
            expected[state_index[row.state_i], state_index[row.state_j]] = _number(row.value)
    return states, time_labels, expected


def _expected_indel(
    rows: pd.DataFrame, states: list[str], time_labels: list[str]
) -> float | np.ndarray:
    values = rows.loc[rows["component"] == "indel"]
    if values.iloc[0]["state_i"] == "":
        assert len(values) == 1
        return _number(values.iloc[0]["value"])
    state_index = {state: i for i, state in enumerate(states)}
    if time_labels:
        time_index = {label: i for i, label in enumerate(time_labels)}
        assert len(values) == len(states) * len(time_labels)
        assert not values[["state_i", "time_label"]].duplicated().any()
        expected = np.full((len(states), len(time_labels)), np.nan, dtype=float)
        for row in values.itertuples(index=False):
            expected[state_index[row.state_i], time_index[row.time_label]] = _number(row.value)
        return expected
    assert len(values) == len(states)
    assert not values[["state_i"]].duplicated().any()
    expected = np.full(len(states), np.nan, dtype=float)
    for row in values.itertuples(index=False):
        expected[state_index[row.state_i]] = _number(row.value)
    return expected


def _assert_case(table: pd.DataFrame, case: str, spec: Case) -> None:
    builder, kwargs = spec
    rows = table.loc[table["case"] == case]
    assert not rows.empty, f"missing R golden case {case}"
    component = "prox" if kwargs.get("proximities") else "sm"
    states, time_labels, expected_matrix = _expected_matrix(rows, component)
    expected_indel = _expected_indel(rows, states, time_labels)

    result = seqcost(builder(), **kwargs)
    assert set(result) == {"indel", component}

    matrix = result[component]
    if time_labels:
        assert matrix.shape == expected_matrix.shape
        assert list(matrix.from_states) == states
        assert list(matrix.to_states) == states
        assert list(matrix.time_labels) == time_labels
    else:
        assert isinstance(matrix, pd.DataFrame)
        assert matrix.index.tolist() == states
        assert matrix.columns.tolist() == states
    np.testing.assert_allclose(
        np.asarray(matrix), expected_matrix, rtol=0, atol=ATOL, equal_nan=True
    )

    actual_indel = result["indel"]
    if isinstance(expected_indel, np.ndarray) and expected_indel.ndim == 1:
        assert isinstance(actual_indel, pd.Series)
        assert actual_indel.index.tolist() == states
    elif isinstance(expected_indel, np.ndarray):
        assert isinstance(actual_indel, pd.DataFrame)
        assert actual_indel.index.tolist() == states
        assert actual_indel.columns.tolist() == time_labels
    if isinstance(expected_indel, np.ndarray):
        np.testing.assert_allclose(
            np.asarray(actual_indel), expected_indel, rtol=0, atol=ATOL, equal_nan=True
        )
    else:
        assert float(actual_indel) == pytest.approx(expected_indel, rel=0, abs=ATOL, nan_ok=True)


def test_case_maps_cover_every_oracle_result():
    assert set(_load("seqcost_stationary_cases.csv")["case"]) == set(STATIONARY_CASES)
    assert set(_load("seqcost_timevarying_cases.csv")["case"]) == set(TIMEVARYING_CASES)
    assert set(_load("seqcost_proximity_cases.csv")["case"]) == set(PROXIMITY_CASES)


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ At\ least\ one\ substitution\ cost\ greater\ than\ twice\ the\ indel\ cost\.:traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ \ \ \ \ \ for\ example,\ cost\ is\ zero\ for\ the\ pair\ \(\ A\ ,\ C\ \):traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ \ at\ least,\ one\ substitution\ cost\ doesn't\ respect\ the\ triangle\ inequality\.:traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ \[!\]\ \ \d+\ substitution\ cost\ value\(s\)\ are\ \d+!:traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ at\ least\ one\ indel\ cost\ does\ not\ respect\ the\ triangle\ inequality\.:traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.parametrize("case", list(STATIONARY_CASES))
def test_stationary_seqcost_matches_r(case: str):
    _assert_case(_load("seqcost_stationary_cases.csv"), case, STATIONARY_CASES[case])


@pytest.mark.parametrize("case", list(TIMEVARYING_CASES))
def test_timevarying_seqcost_matches_r(case: str):
    _assert_case(_load("seqcost_timevarying_cases.csv"), case, TIMEVARYING_CASES[case])


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ At\ least\ one\ substitution\ cost\ greater\ than\ twice\ the\ indel\ cost\.:traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ \ \ \ \ \ for\ example,\ cost\ is\ zero\ for\ the\ pair\ \(\ A\ ,\ C\ \):traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ \[!\]\ \ \d+\ substitution\ cost\ value\(s\)\ are\ \d+!:traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ max\ substitution\ cost\ less\ or\ equal\ to\ \d+!:traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.parametrize("case", list(PROXIMITY_CASES))
def test_seqcost_proximities_match_r(case: str):
    _assert_case(_load("seqcost_proximity_cases.csv"), case, PROXIMITY_CASES[case])


def test_no_observed_missing_still_adds_the_requested_missing_state():
    """``with_missing=True`` widens the matrix even with no observed missing cell.

    TraMineR 2.2-14 downgraded the flag here and warned; 2.2-14 removed that
    block because ``seqdist`` may legitimately have decided ``with.missing``
    from an external ``refseq`` instead (the Peyrot issue).  The paired golden
    case is ``constant_no_observed_missing`` in
    ``tests/r_goldens/seqcost_stationary_cases.csv``, regenerated under 2.2-14:
    a 5x5 matrix over ``A,B,C,D,*``.  See
    ``tests/test_seqcost_refseq_missing_contract.py`` for the end-to-end route.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        result = seqcost(fixture_a(), "CONSTANT", with_missing=True)
    assert result["sm"].shape == (5, 5)


def test_oracle_error_artifact_covers_every_refusal():
    errors = pd.read_csv(GOLDEN_DIR / "seqcost_stationary_errors.txt", sep="|", dtype=str)
    assert set(errors["case"]) == {
        "non_stslist",
        "missing_method",
        "invalid_method",
        "trate_one_state_stationary",
        "trate_one_state_timevarying",
        "constant_one_state_timevarying",
        "invalid_transition",
        "constant_na",
        "features_missing",
        "features_wrong_rows",
        "features_zero_columns",
        "features_bad_weights",
        "features_unnamed_type",
        "features_bad_type_column",
        "features_bad_binary",
        "future_timevarying",
        "features_timevarying",
        "trate_time_lag2_previous",
        "trate_time_lag2_next",
        "trate_time_lag2_both",
        "trate_time_negative_lag",
    }


def test_python_refusals_preserve_r_supported_boundary():
    x = fixture_a()
    with pytest.raises(TypeError, match="STSList|sequence object"):
        seqcost(np.array([["A"]]), "CONSTANT")
    with pytest.raises(TypeError):
        seqcost(x)  # type: ignore[call-arg]
    with pytest.raises(ValueError, match="method must be one of"):
        seqcost(x, "constant")
    for time_varying in (False, True):
        with pytest.raises(ValueError, match="at least two states"):
            seqcost(_one_state(), "TRATE", time_varying=time_varying)
    with pytest.raises(ValueError, match="at least two states"):
        seqcost(_one_state(), "CONSTANT", time_varying=True)
    with pytest.raises(ValueError, match="transition must be one of"):
        seqcost(x, "CONSTANT", transition="bad")
    with pytest.raises(ValueError, match="constant substitution-cost"):
        seqcost(x, "CONSTANT", cval=np.nan)
    with pytest.raises(ValueError, match="one row per state"):
        seqcost(x, "FEATURES")
    with pytest.raises(ValueError, match="one row per state"):
        seqcost(x, "FEATURES", state_features=pd.DataFrame({"x": [1, 2, 3]}))
    with pytest.raises(ValueError, match="at least one feature"):
        seqcost(x, "FEATURES", state_features=pd.DataFrame(index=range(4)))
    with pytest.raises(ValueError, match="2D table|one row per state"):
        seqcost(x, "FEATURES", state_features=np.array(1.0))
    with pytest.raises(ValueError, match="weights.*length|length.*weights"):
        seqcost(
            x,
            "FEATURES",
            state_features=_features_constant(),
            feature_weights=[1, 2, 3],
        )
    with pytest.raises(ValueError, match="named.*mapping"):
        seqcost(
            x,
            "FEATURES",
            state_features=_features_mixed(),
            feature_type=["numeric"],
        )
    with pytest.raises(ValueError, match="invalid column|not_a_column"):
        seqcost(
            x,
            "FEATURES",
            state_features=_features_mixed(),
            feature_type={"asymm": "not_a_column"},
        )
    with pytest.raises(ValueError, match="binary"):
        seqcost(
            x,
            "FEATURES",
            state_features=_features_mixed(),
            feature_type={"asymm": "category"},
        )
    with pytest.raises(ValueError, match="time.varying|time_varying"):
        seqcost(x, "FUTURE", time_varying=True)
    with pytest.raises(ValueError, match="time.varying|time_varying"):
        seqcost(x, "FEATURES", time_varying=True, state_features=_features_constant())
    for transition in ("previous", "next", "both"):
        with pytest.raises(ValueError, match="lag"):
            seqcost(x, "TRATE", time_varying=True, transition=transition, lag=2)
    with pytest.raises(ValueError, match="lag"):
        seqcost(x, "TRATE", time_varying=True, lag=-1)


def test_seqcost_is_exported_exactly_once_at_both_public_surfaces():
    assert tx.seqcost is seqcost
    assert PUBLIC_API.count("seqcost") == 1
    assert tx.__all__.count("seqcost") == 1
