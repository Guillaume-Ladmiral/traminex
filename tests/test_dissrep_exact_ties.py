"""Exactly tied density scores are ordered by R's tie rule on every build (TX-AUD-20261019-017).

dissrep scores a candidate by the summed weights of its neighbours and orders the candidates
stably (R's order(); dissrep.R:40-43, 84). When two sums are equal in exact arithmetic, the
floating-point summation order of the build's BLAS decided which came first: on mvad, R's
OpenBLAS and NumPy disagreed and the fourth representative differed. Scores that the summation
error could have split are summed again with correct rounding, so an exact tie is a tie and the
lower row comes first.
"""

from __future__ import annotations

import warnings

import numpy as np

import traminex as tx


def _tie():
    # Rows 0 and 1 each have three neighbours (distance 0.05 < 0.1 * dmax): rows 0, 2, 3 with
    # weights 0.3, 0.2, 0.1 and rows 1, 4, 5 with weights 0.1, 0.2, 0.3. The exact sums are equal
    # (the same three numbers); left-to-right addition gives 0.6 and 0.6000000000000001.
    # Every other row has fewer or lighter neighbours (0.5, 0.4, 0.3, 0.4).
    weights = np.array([0.3, 0.1, 0.2, 0.1, 0.2, 0.3])
    diss = np.ones((6, 6))
    np.fill_diagonal(diss, 0.0)
    for i, j in ((0, 2), (0, 3), (1, 4), (1, 5)):
        diss[i, j] = diss[j, i] = 0.05
    return diss, weights


def test_an_exact_density_tie_keeps_the_lower_row_first() -> None:
    diss, weights = _tie()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = tx.dissrep(diss, weights=weights, nrep=1)
    assert np.asarray(result["representatives"]).tolist() == [0]
    assert result["scores"][0] == result["scores"][1]


def test_the_frequency_and_distance_criteria_use_the_same_rule() -> None:
    diss, weights = _tie()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for criterion in ("freq", "dist"):
            result = tx.dissrep(diss, weights=weights, nrep=1, criterion=criterion)
            scores = np.asarray(result["scores"], dtype=float)
            order = np.argsort(-scores if criterion == "freq" else scores, kind="stable")
            assert np.asarray(result["representatives"]).tolist() == [int(order[0])]
