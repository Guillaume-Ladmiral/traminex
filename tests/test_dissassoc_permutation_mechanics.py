"""Permutation mechanics of ``dissassoc`` replayed against TraMineR 2.2-14.

Audit 2026-09-22, finding S4-01: under ``weight_permutation="diss"`` traminex
computed the Levene residuals at the *unpermuted* group positions, so the
permuted Levene statistic never changed and its p-value was always 1. R's
``dissassocweighted.permdiss`` takes them at the permuted positions
(``SOURCES/TraMineR/R/dissassocweighted.R:119-125``: ``wind <-
sort.int(ind[indgrp[[i]]])`` then ``WeightedResidualsZ(dissc[wind], ww[wind])``).
The unpermuted positions belong to the *group* branch only (lines 96-98).

R and NumPy draw different random permutations, so the R generator
(``scripts/r_goldens/make_goldens_audit20260922_s4.R``) replaces ``sample`` with a
fixed, deterministic permutation script that it writes into the fixture; this
test replays the same script through a NumPy ``Generator`` and compares every
permuted statistic with R's, cell by cell.

Finding S4-06: with ``weights`` and ``weight_permutation="none"`` traminex
returned a mixture of weighted sums of squares and unweighted group counts
(a negative Bartlett statistic, for example). No branch of R's
``dissassocweighted`` handles that combination
(``dissassocweighted.R:342-359``), so R stops; traminex now refuses it.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from r_goldens_utils import R_GOLDEN_ATOL

from traminex.diss import dissassoc

FIXTURES = Path(__file__).with_name("r_goldens")


class ScriptedGenerator(np.random.Generator):
    """A generator whose ``permutation`` replays a fixed 0-based script."""

    def __init__(self, permutations: list[np.ndarray]):
        super().__init__(np.random.PCG64(0))
        self._permutations = permutations
        self.calls = 0

    def permutation(self, x):  # type: ignore[override]
        script = self._permutations[self.calls]
        self.calls += 1
        assert int(x) == script.size
        return script.copy()


def _value(token: str) -> float:
    return float("nan") if token == "NA" else float(token)


def _read_fixture(name: str) -> dict[str, list[list[str]]]:
    rows: dict[str, list[list[str]]] = {}
    for line in (FIXTURES / name).read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        key, values = line.split("\t", 1)
        rows.setdefault(key, []).append(values.split(" "))
    return rows


def _numbers(tokens: list[str]) -> np.ndarray:
    return np.array([_value(token) for token in tokens], dtype=float)


def test_diss_permutation_matches_r_for_every_statistic():
    rows = _read_fixture("dissassoc_diss_scripted_permutations.txt")
    D = _numbers(rows["D"][0])
    D = D.reshape(int(round(np.sqrt(D.size))), -1)
    group = np.array([None if token == "NA" else token for token in rows["group"][0]], object)
    weights = _numbers(rows["weights"][0])
    R = int(rows["R"][0][0])
    script = [np.array([int(token) - 1 for token in row]) for row in rows["perm"]]
    expected_t = np.vstack([_numbers(row) for row in rows["t"]])
    expected_pval = _numbers(rows["pval"][0])

    generator = ScriptedGenerator(script)
    result = dissassoc(
        D, group, weights=weights, R=R, weight_permutation="diss", random_state=generator
    )

    assert generator.calls == R - 1
    np.testing.assert_allclose(result.perms["t"], expected_t, rtol=1e-12, atol=R_GOLDEN_ATOL)
    np.testing.assert_array_equal(result.perms["pval"], expected_pval)
    # The Levene column is genuinely permuted (it was constant before the fix).
    assert np.unique(result.perms["t"][:, 4]).size > 1


def test_weights_with_the_unweighted_permutation_are_refused():
    rows = _read_fixture("dissassoc_weights_none.txt")
    assert rows["none"][0][0] == "error"  # R stops on this combination
    permutation_rows = _read_fixture("dissassoc_diss_scripted_permutations.txt")
    D = _numbers(permutation_rows["D"][0])
    D = D.reshape(int(round(np.sqrt(D.size))), -1)
    group = np.array(
        [None if token == "NA" else token for token in permutation_rows["group"][0]], object
    )
    weights = _numbers(permutation_rows["weights"][0])

    with pytest.raises(ValueError, match="weight_permutation='none'"):
        dissassoc(D, group, weights=weights, R=1, weight_permutation="none")

    # The weighted analysis R does run on these data is unaffected.
    result = dissassoc(D, group, weights=weights, R=1, weight_permutation="diss")
    np.testing.assert_allclose(
        result.stat["t0"].to_numpy(), _numbers(rows["diss_t0"][0]), rtol=1e-12, atol=R_GOLDEN_ATOL
    )
    # Without weights, "none" remains the unweighted analysis.
    unweighted = dissassoc(D, group, R=1, weight_permutation="none")
    assert unweighted.weight_permutation == "none"
