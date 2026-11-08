"""R-golden tests for the weighted ``dissassoc`` (Phase 2.2 / A6).

Anchors traminex's ``dissassoc`` to TraMineR 2.2-14 ``dissassoc(diss, group,
weights=, R=, weight.permutation=, squared=)`` on the shared Phase 2 fixture
(``scripts/r_goldens/make_goldens_dissassoc.R``).

The observed statistics (t0: Pseudo F, Pseudo Fbf, Pseudo R2, Bartlett,
Levene), the pseudo-ANOVA table and the per-group discrepancies are pinned to
R. P-values are permutation-based and stochastic in R as well, so they are
checked structurally (range, reproducibility under a fixed seed, and obvious
signal/no-signal behaviour), per MIGRATION_PLAN.md Phase 2.2.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest
from r_goldens_utils import R_GOLDEN_ATOL, load_golden, load_golden_lines

from traminex.diss import dissassoc

GROUPS = np.repeat([1, 2, 3], 5)
R2_CASES = Path(__file__).with_name("r_goldens") / "dissassoc_r2_cases.csv"


def _inputs():
    D = load_golden("phase2_input_D")
    wreal = load_golden("phase2_input_weights_real").ravel()
    wint = load_golden("phase2_input_weights_int").ravel()
    return D, wreal, wint


def _assert_case(res, tag):
    np.testing.assert_allclose(
        res.anova_table.to_numpy(),
        load_golden(f"dissassoc_{tag}_anova"),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )
    np.testing.assert_allclose(
        res.groups.to_numpy(), load_golden(f"dissassoc_{tag}_groups"), rtol=0, atol=R_GOLDEN_ATOL
    )
    np.testing.assert_allclose(
        res.stat["t0"].to_numpy(),
        load_golden(f"dissassoc_{tag}_t0").ravel(),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


def test_dissassoc_unweighted_matches_r():
    D, _, _ = _inputs()
    res = dissassoc(D, GROUPS, R=1)
    assert res.weight_permutation == "none"
    _assert_case(res, "unw")
    assert np.all(np.isnan(res.stat["p_value"].to_numpy()))  # R: R=1 -> NA pvals


def test_dissassoc_weighted_diss_matches_r():
    D, wreal, _ = _inputs()
    res = dissassoc(D, GROUPS, weights=wreal, R=1, weight_permutation="diss")
    _assert_case(res, "wreal")


def test_dissassoc_weighted_replicate_matches_r():
    D, _, wint = _inputs()
    res = dissassoc(D, GROUPS, weights=wint, R=1, weight_permutation="replicate")
    _assert_case(res, "wint_rep")


def test_dissassoc_weighted_squared_matches_r():
    D, wreal, _ = _inputs()
    res = dissassoc(D, GROUPS, weights=wreal, R=1, weight_permutation="diss", squared=True)
    np.testing.assert_allclose(
        res.stat["t0"].to_numpy(),
        load_golden("dissassoc_wreal_sq_t0").ravel(),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


def test_dissassoc_permutation_structure():
    # Strong group structure in the fixture -> minimal possible p-value for
    # Pseudo F; deterministic under a fixed seed; pval includes the observed
    # permutation (R convention), so it is >= 1/R.
    D, wreal, wint = _inputs()
    res1 = dissassoc(D, GROUPS, weights=wreal, R=100, weight_permutation="diss", random_state=42)
    res2 = dissassoc(D, GROUPS, weights=wreal, R=100, weight_permutation="diss", random_state=42)
    p1 = res1.stat["p_value"].to_numpy()
    np.testing.assert_array_equal(p1, res2.stat["p_value"].to_numpy())
    assert np.all(p1[~np.isnan(p1)] >= 1.0 / 100.0)
    assert np.all(p1[~np.isnan(p1)] <= 1.0)
    assert p1[0] == 1.0 / 100.0  # Pseudo F: observed grouping is the best
    # Under weight_permutation="diss" R takes the Levene residuals at the
    # permuted positions (dissassocweighted.R:119-125), so the permuted Levene
    # statistic varies. This line used to pin a constant column (p-value 1),
    # which was audit finding S4-01; R's permuted statistics are compared cell
    # by cell in tests/test_dissassoc_permutation_mechanics.py.
    assert np.unique(res1.perms["t"][:, 4]).size > 1

    # Replicate and unweighted variants run end-to-end with permutations.
    res3 = dissassoc(D, GROUPS, weights=wint, R=50, weight_permutation="replicate", random_state=7)
    assert res3.stat["p_value"]["Pseudo F"] == 1.0 / 50.0
    res4 = dissassoc(D, GROUPS, R=50, random_state=7)
    assert res4.stat["p_value"]["Pseudo F"] == 1.0 / 50.0
    res5 = dissassoc(D, GROUPS, weights=wreal, R=50, weight_permutation="group", random_state=7)
    assert res5.stat["p_value"]["Pseudo F"] == 1.0 / 50.0

    # No-signal labels -> large Pseudo F p-value.
    rng = np.random.default_rng(0)
    noise_groups = rng.permutation(GROUPS)
    res6 = dissassoc(D, noise_groups, R=200, random_state=3)
    # t0 of noise grouping should not be the extreme of its own permutation
    # distribution with overwhelming probability.
    assert res6.stat["p_value"]["Pseudo F"] > 0.05


def test_dissassoc_drops_missing_groups():
    D, wreal, _ = _inputs()
    groups = GROUPS.astype(float)
    groups_with_nan = groups.copy()
    # Make case 0 missing; compare against explicit removal.
    groups_with_nan[0] = np.nan
    res_nan = dissassoc(D, groups_with_nan, weights=wreal, R=1, weight_permutation="diss")
    keep = np.arange(1, 15)
    res_sub = dissassoc(
        D[np.ix_(keep, keep)], groups[keep], weights=wreal[keep], R=1, weight_permutation="diss"
    )
    np.testing.assert_allclose(
        res_nan.stat["t0"].to_numpy(), res_sub.stat["t0"].to_numpy(), rtol=0, atol=1e-12
    )


def _r2_rows():
    with R2_CASES.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _case_rows(case, section):
    return [row for row in _r2_rows() if row["case"] == case and row["section"] == section]


def _matrix(case, section):
    rows = _case_rows(case, section)
    result = np.full(
        (max(int(row["row"]) for row in rows), max(int(row["col"]) for row in rows)),
        np.nan,
    )
    for row in rows:
        result[int(row["row"]) - 1, int(row["col"]) - 1] = float(row["value"])
    return result


def _r2_config(case):
    row = next(row for row in _r2_rows() if row["case"] == case)
    return {
        "method": row["method"],
        "weights_kind": row["weights_kind"],
        "R": int(row["R"]),
        "squared": row["squared"] == "TRUE",
        "group_kind": row["group_kind"],
    }


_R2_CASE_NAMES = [
    "uw_none_r31",
    "uw_requested_permute_r0",
    "uw_requested_random_sampling_r1",
    "wint_diss_r0",
    "wreal_diss_r31_sq",
    "wint_replicate_r1",
    "wint_replicate_r31",
    "wreal_rounded_replicate_r1",
    "wreal_rounded_replicate_r31_two",
    "wreal_group_r1",
    "wreal_group_r31",
    "wint_permute_r1",
    "wint_permute_r31",
    "wreal_permute_r31_sq_two",
    "wint_random_sampling_r1",
    "wint_random_sampling_r31",
    "wreal_random_sampling_r31_sq_two",
]


def _run_r2(case, random_state=4401):
    D, wreal, wint = _inputs()
    config = _r2_config(case)
    weights = {"none": None, "integer": wint, "fractional": wreal}[config["weights_kind"]]
    groups = np.resize(np.asarray([1, 2]), D.shape[0]) if config["group_kind"] == "two" else GROUPS
    return dissassoc(
        D,
        groups,
        weights=weights,
        R=config["R"],
        weight_permutation=config["method"],
        squared=config["squared"],
        random_state=random_state,
    )


@pytest.mark.parametrize("case", _R2_CASE_NAMES)
def test_dissassoc_r2_deterministic_components_match_r(case):
    assert {row["case"] for row in _r2_rows()} == set(_R2_CASE_NAMES)
    out = _run_r2(case)
    np.testing.assert_allclose(
        out.groups.to_numpy(), _matrix(case, "groups"), rtol=0, atol=R_GOLDEN_ATOL
    )
    np.testing.assert_allclose(
        out.anova_table.to_numpy(), _matrix(case, "anova"), rtol=0, atol=R_GOLDEN_ATOL
    )
    np.testing.assert_allclose(
        out.stat["t0"].to_numpy(),
        _matrix(case, "t0").ravel(),
        rtol=0,
        atol=R_GOLDEN_ATOL,
        equal_nan=True,
    )
    attribute = _case_rows(case, "attribute")[0]["text"]
    assert out.weight_permutation == attribute.split("=", 1)[1]
    assert list(out.stat.index) == ["Pseudo F", "Pseudo Fbf", "Pseudo R2", "Bartlett", "Levene"]
    assert list(out.stat.columns) == ["t0", "p_value"]
    assert out.perms["R"] == _r2_config(case)["R"]


@pytest.mark.parametrize("case", [case for case in _R2_CASE_NAMES if _r2_config(case)["R"] > 1])
def test_dissassoc_r2_sampled_results_have_valid_schema_and_probabilities(case):
    out = _run_r2(case)
    assert out.perms["t"].shape == (out.R, 5)
    np.testing.assert_allclose(
        out.perms["t"][0], out.perms["t0"], equal_nan=True, rtol=0, atol=R_GOLDEN_ATOL
    )
    pvalues = out.stat["p_value"].to_numpy()
    finite = np.isfinite(pvalues)
    assert np.all((pvalues[finite] >= 1 / out.R) & (pvalues[finite] <= 1))


@pytest.mark.parametrize("method", ["permute", "random-sampling"])
def test_dissassoc_new_modes_are_reproducible_and_do_not_mutate_inputs(method):
    D, wreal, _ = _inputs()
    D_before = D.copy()
    weights_before = wreal.copy()
    first = dissassoc(D, GROUPS, weights=wreal, R=31, weight_permutation=method, random_state=991)
    second = dissassoc(D, GROUPS, weights=wreal, R=31, weight_permutation=method, random_state=991)
    different = dissassoc(
        D, GROUPS, weights=wreal, R=31, weight_permutation=method, random_state=992
    )
    np.testing.assert_array_equal(first.perms["t"], second.perms["t"])
    assert not np.array_equal(first.perms["t"][1:], different.perms["t"][1:])
    np.testing.assert_array_equal(D, D_before)
    np.testing.assert_array_equal(wreal, weights_before)


def test_dissassoc_r2_oracle_conditions_and_refusals_are_frozen():
    conditions = {row["case"]: row["text"] for row in _r2_rows() if row["section"] == "condition"}
    assert conditions == {
        "wreal_rounded_replicate_r1": "MESSAGE:Weigths loss : 3.5 (0.179487179487179)",
        "wreal_rounded_replicate_r31_two": "MESSAGE:Weigths loss : 3.5 (0.179487179487179)",
        "wint_random_sampling_r1": "MESSAGE: [>] weights corrected to match sample size. Mean correction 1",
        "wint_random_sampling_r31": "MESSAGE: [>] weights corrected to match sample size. Mean correction 1",
        "wreal_random_sampling_r31_sq_two": "MESSAGE: [>] weights corrected to match sample size. Mean correction 0.3",
    }
    assert load_golden_lines("dissassoc_r2_errors") == [
        "invalid_method=ERROR:attempt to set 'rownames' on an object with no dimensions",
        "fractional_replicate=ERROR: [!] to permute replicate, you should specify integer weights",
        "weights_short=ERROR:missing value where TRUE/FALSE needed",
        "group_short=VALUE",
        "negative_R=VALUE",
        "negative_sampling_weight=ERROR:negative probability",
        "zero_sampling_weights=ERROR:missing value where TRUE/FALSE needed",
    ]
    D, wreal, _ = _inputs()
    with pytest.raises(ValueError, match="weight_permutation"):
        dissassoc(D, GROUPS, weights=wreal, R=1, weight_permutation="bogus")
    with pytest.raises(ValueError, match="integers"):
        dissassoc(D, GROUPS, weights=wreal, R=1, weight_permutation="replicate")
    with pytest.raises(ValueError, match="weights"):
        dissassoc(D, GROUPS, weights=wreal[:-1], R=1, weight_permutation="diss")
