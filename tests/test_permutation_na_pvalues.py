"""Permutation p-values are NaN wherever TraMineR's are NA.

Audit 2026-09-22, findings S4-03 and S4-07. traminex counted ``t >= t0`` with
NumPy, where a comparison with NaN is simply false, so a statistic that was NaN
in the observed data got a p-value of 0 (reading as maximally significant) and
a NaN among the permuted statistics lowered the count; ``dissmfacw`` with
``nperm=1`` reported p-values of 1.

TraMineR 2.2-14 computes ``sum(t[, i] >= t0[i]) / R``
(``SOURCES/TraMineR/R/TraMineR-permutation.R:21``), which is NA as soon as the
column holds an NA, observed or permuted, and returns NA p-values when
``R <= 1`` (lines 25-26). ``dissassoc`` (all weighting modes) and ``dissmfacw``
(``dissmfacw.R:124``) use that function. Expected values come from live R
through ``scripts/r_goldens/make_goldens_audit20260922_s4.R``; the permutations
are replayed from the fixture, as in ``test_dissassoc_permutation_mechanics``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from r_goldens_utils import R_GOLDEN_ATOL
from test_dissassoc_permutation_mechanics import ScriptedGenerator

from traminex import dissmfacw
from traminex.diss import dissassoc

FIXTURE = Path(__file__).with_name("r_goldens") / "permutation_na_pvalues.txt"


def _rows() -> dict[str, dict[str, list[list[str]]]]:
    cases: dict[str, dict[str, list[list[str]]]] = {}
    for line in FIXTURE.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        case, key, values = line.split("\t")
        cases.setdefault(case, {}).setdefault(key, []).append(values.split(" "))
    return cases


CASES = _rows()


def _numbers(tokens: list[str]) -> np.ndarray:
    return np.array([np.nan if token == "NA" else float(token) for token in tokens])


def _matrix(tokens: list[str]) -> np.ndarray:
    values = _numbers(tokens)
    return values.reshape(int(round(np.sqrt(values.size))), -1)


def _replay(case: str):
    rows = CASES[case]
    script = [np.array([int(token) - 1 for token in row]) for row in rows["perm"]]
    generator = ScriptedGenerator(script)
    result = dissassoc(
        _matrix(rows["D"][0]),
        np.array(rows["group"][0], dtype=object),
        R=len(script) + 1,
        random_state=generator,
    )
    assert generator.calls == len(script)
    return result, np.vstack([_numbers(row) for row in rows["t"]]), _numbers(rows["pval"][0])


def test_observed_nan_statistic_has_na_pvalue():
    result, expected_t, expected_pval = _replay("singleton")
    np.testing.assert_allclose(result.perms["t"], expected_t, rtol=1e-12, atol=R_GOLDEN_ATOL)
    assert np.isnan(expected_pval[3])  # R: NA for the Bartlett column
    np.testing.assert_array_equal(result.perms["pval"], expected_pval)
    np.testing.assert_array_equal(result.stat["p_value"].to_numpy(), expected_pval)


def test_permuted_nan_statistic_has_na_pvalue():
    result, expected_t, expected_pval = _replay("permuted_nan")
    np.testing.assert_allclose(result.perms["t"], expected_t, rtol=1e-12, atol=R_GOLDEN_ATOL)
    assert not np.isnan(expected_t[0, 3]) and np.isnan(expected_t[2, 3])
    np.testing.assert_array_equal(result.perms["pval"], expected_pval)


def test_random_sampling_observed_nan_statistic_has_na_pvalue():
    rows = CASES["random_sampling"]
    singleton = CASES["singleton"]
    D = _matrix(singleton["D"][0])
    group = np.array(singleton["group"][0], dtype=object)
    result = dissassoc(
        D,
        group,
        weights=np.ones(group.size),
        R=20,
        weight_permutation="random-sampling",
        random_state=20260922,
    )
    np.testing.assert_allclose(
        result.perms["t0"], _numbers(rows["t0"][0]), rtol=1e-12, atol=R_GOLDEN_ATOL
    )
    assert np.isnan(_numbers(rows["bartlett_pval"][0])[0])
    assert np.isnan(result.perms["pval"][3])


def test_dissmfacw_single_permutation_has_na_pvalues():
    rows = CASES["dissmfacw_R1"]
    singleton = CASES["singleton"]
    D = _matrix(singleton["D"][0])
    factors = {"g": singleton["group"][0], "h": [1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1]}
    result = dissmfacw(D, factors, nperm=1)
    np.testing.assert_allclose(
        result.pseudo_F, _numbers(rows["PseudoF"][0]), rtol=1e-12, atol=R_GOLDEN_ATOL
    )
    np.testing.assert_allclose(
        result.pseudo_R2, _numbers(rows["PseudoR2"][0]), rtol=1e-12, atol=R_GOLDEN_ATOL
    )
    assert np.isnan(_numbers(rows["p_value"][0])).all()
    assert np.isnan(result.p_value).all()
    assert result.mfac["p_value"].isna().all()


def test_dissmfacw_nan_statistic_has_na_pvalue_with_permutations():
    rows = CASES["dissmfacw_zero"]
    groups = CASES["singleton"]["group"][0]
    factors = {"g": groups, "h": [1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1]}
    result = dissmfacw(np.zeros((13, 13)), factors, nperm=5, seed=1)
    assert np.isnan(_numbers(rows["PseudoF"][0])).all()
    assert np.isnan(result.pseudo_F).all()
    assert np.isnan(_numbers(rows["p_value"][0])).all()
    assert np.isnan(result.p_value).all()
