"""Ranks follow R's LINPACK dqrdc2 rule (TX-REP-20261019-002, batch-11 ruling R11-m).

R's ``qr()`` keeps a column unless its remaining norm falls below ``tol`` times its own original norm, which it then
moves to the end (``src/appl/dqrdc2.f``). traminex's ``dissmfacw`` counted pivots of a pivoted QR above ``1e-7`` of
the largest, so a column scaled by 1e9 made the others vanish (rank 1 or 2 where R has 3); the binomial fit behind
``seqerulesdisc`` chose its columns once with ``np.linalg.matrix_rank`` and refused a design that its least-squares
step found rank deficient. ``dissmfacw`` now uses dqrdc2 with R's ``qr()`` tolerance 1e-7, and the fit applies it as
``glm.fit`` does: with ``tol = 1e-11`` to the weighted design at every iteration, aliasing the columns past the rank.
Fixture: ``scripts/r_goldens/make_goldens_rank_and_convergence.R``.
"""

from __future__ import annotations

import csv
import warnings
from pathlib import Path

import numpy as np
import pytest

import traminex as tx
from traminex._rlinpack import dqrdc2_pivot
from traminex.discrete_rules import _binomial_irls

GOLDEN_PATH = Path(__file__).resolve().parent / "r_goldens" / "rank_and_convergence.csv"
with open(GOLDEN_PATH, encoding="utf-8", newline="") as handle:
    GOLDEN = {(row["case"], row["field"]): row["value"] for row in csv.DictReader(handle)}
QR_CASES = sorted(case for case, field in GOLDEN if field == "qr_rank")


def _matrix(text: str) -> np.ndarray:
    return np.array([[float(cell) for cell in row.split(",")] for row in text.split("|")])


def _vector(text: str) -> np.ndarray:
    return np.array([np.nan if cell == "NA" else float(cell) for cell in text.split()])


def test_the_golden_covers_the_designs() -> None:
    assert len(QR_CASES) == 20


@pytest.mark.parametrize("case", QR_CASES)
def test_dqrdc2_pivot_gives_r_s_qr_rank_and_pivot(case: str) -> None:
    rank, pivot = dqrdc2_pivot(_matrix(GOLDEN[(case, "design")]), 1e-7)
    assert rank == int(GOLDEN[(case, "qr_rank")])
    assert [p + 1 for p in pivot] == [int(p) for p in GOLDEN[(case, "qr_pivot")].split()]


@pytest.mark.parametrize("case", ["mfac_scaled", "mfac_scaled_weighted"])
def test_dissmfacw_on_a_scaled_design_equals_r(case: str) -> None:
    factors = {
        "a": _vector(GOLDEN[("mfac_scaled", "a")]),
        "b": _vector(GOLDEN[("mfac_scaled", "b")]),
    }
    weights = _vector(GOLDEN[("mfac_scaled", "weights")]) if case.endswith("weighted") else None
    result = tx.dissmfacw(_matrix(GOLDEN[("mfac_scaled", "diss")]), factors, weights, nperm=1)
    assert list(result.variables) == GOLDEN[(case, "variables")].split()
    np.testing.assert_allclose(
        result.pseudo_F, _vector(GOLDEN[(case, "pseudo_f")]), rtol=0, atol=1e-10
    )
    np.testing.assert_allclose(
        result.pseudo_R2, _vector(GOLDEN[(case, "pseudo_r2")]), rtol=0, atol=1e-10
    )


@pytest.mark.parametrize("link", ["logit", "cloglog"])
@pytest.mark.parametrize("k", [10, 12])
def test_the_binomial_fit_aliases_columns_as_glm_does(k: int, link: str) -> None:
    y = _vector(GOLDEN[("glm_response", "y")])
    design = _matrix(GOLDEN[(f"glm_near_collinear_1e-{k}", "design")])
    fit = _binomial_irls(y, design, link=link)
    coefficients = _vector(GOLDEN[(f"glm_near_collinear_1e-{k}_{link}", "coefficients")])
    pvalues = _vector(GOLDEN[(f"glm_near_collinear_1e-{k}_{link}", "pvalues")])
    np.testing.assert_array_equal(np.isnan(fit.coefficients), np.isnan(coefficients))
    np.testing.assert_array_equal(np.isnan(fit.pvalues), np.isnan(pvalues))
    if k == 12:
        # At 1e-10 R keeps the third column; its coefficients (about 6e9) are too ill-conditioned to compare.
        kept = ~np.isnan(coefficients)
        np.testing.assert_allclose(fit.coefficients[kept], coefficients[kept], rtol=2e-7, atol=2e-8)
        np.testing.assert_allclose(fit.pvalues[kept], pvalues[kept], rtol=2e-6, atol=2e-8)


# TX-AUD-20261019-085: on separated data R's glm.fit stops after maxit = 25 iterations, warns "glm.fit: algorithm did
# not converge" and returns the last iterate; the fit raised RuntimeError. The golden's data are the reproduction's.
@pytest.mark.parametrize("link", ["logit", "cloglog"])
def test_a_fit_that_does_not_converge_warns_and_keeps_the_last_iterate(link: str) -> None:
    case = f"glm_separated_{link}"
    assert (GOLDEN[(case, "converged")], GOLDEN[(case, "iterations")]) == ("FALSE", "25")
    r_warning = GOLDEN[(case, "warnings")].split(" | ")[0]
    assert r_warning == "glm.fit: algorithm did not converge"
    x = _vector(GOLDEN[("glm_separated", "x")])
    y = _vector(GOLDEN[("glm_separated", "y")])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        fit = _binomial_irls(y, np.column_stack([np.ones(len(x)), x]), link=link)
    # R's second warning there (fitted probabilities 0 or 1) is checked with the others below.
    assert str(caught[0].message) == r_warning and caught[0].category is UserWarning
    np.testing.assert_allclose(
        fit.coefficients, _vector(GOLDEN[(case, "coefficients")]), rtol=2e-7, atol=2e-8
    )


# TX-REP-20261022-004: outside [-30, 30] R's logit linkinv and mu.eta (C code) give mu = eps/(1+eps) or 1/(1+eps)
# and mu.eta = eps; traminex used expit and mu(1-mu) floored at eps, so the working weights of separated rows were up
# to about 155 times R's and the logit fits drifted (2.8e-6 on gap's slope). After the fix the fits agree to about
# 1e-11; the bar rtol 1e-8, atol 1e-9 leaves room for last-bit arithmetic on separated data.
@pytest.mark.parametrize("link", ["logit", "cloglog"])
@pytest.mark.parametrize("data", ["gap", "rand30", "separated"])
def test_separated_fits_follow_r_s_binomial_family(data: str, link: str) -> None:
    x = _vector(GOLDEN[(f"glm_{data}", "x")])
    y = _vector(GOLDEN[(f"glm_{data}", "y")])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit = _binomial_irls(y, np.column_stack([np.ones(len(x)), x]), link=link)
    expected = _vector(GOLDEN[(f"glm_{data}_{link}", "coefficients")])
    np.testing.assert_allclose(fit.coefficients, expected, rtol=1e-8, atol=1e-9)


# TX-REP-20261022-005: after the fit glm.fit warns "glm.fit: fitted probabilities numerically 0 or 1 occurred" when a
# fitted mean is above 1 - 10 eps or below 10 eps (after its non-convergence warning); traminex did not. Every golden
# fit whose warnings R recorded: the same messages in the same order (none on the near-collinear designs).
WARNING_CASES = sorted(case for case, field in GOLDEN if field == "warnings")


@pytest.mark.parametrize("case", WARNING_CASES)
def test_the_binomial_fit_warns_as_glm_does(case: str) -> None:
    data, link = case.rsplit("_", 1)
    if data.startswith("glm_near_collinear"):
        design, y = _matrix(GOLDEN[(data, "design")]), _vector(GOLDEN[("glm_response", "y")])
    else:
        x, y = _vector(GOLDEN[(data, "x")]), _vector(GOLDEN[(data, "y")])
        design = np.column_stack([np.ones(len(x)), x])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        _binomial_irls(y, design, link=link)
    expected = [message for message in GOLDEN[(case, "warnings")].split(" | ") if message]
    assert [str(w.message) for w in caught if issubclass(w.category, UserWarning)] == expected


def test_the_golden_records_r_s_warnings_for_ten_fits() -> None:
    assert len(WARNING_CASES) == 10
