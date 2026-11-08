"""Numeric edges of ``seqdist`` against live TraMineR 2.2-14.

Expected values: ``scripts/r_goldens/make_goldens_seqdist_numeric_edges.R``
(2026-09-22 audit S1a-05, S1b-08, S1b-10, S1b-12).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

import traminex as tx

GOLDENS = Path(__file__).parent / "r_goldens"


def _r(case: str) -> np.ndarray:
    for line in (GOLDENS / "seqdist_numeric_edges.txt").read_text().splitlines():
        name, _, values = line.partition("\t")
        if name == case:
            flat = np.array(values.split(" "), dtype=float)
            size = int(round(np.sqrt(flat.size)))
            return flat.reshape(size, size).T if size * size == flat.size else flat
    raise KeyError(case)


FIVE = tx.seqdef(
    np.array([list("ABBCA"), list("BBAAC"), list("CAABB")], dtype=object), alphabet=["A", "B", "C"]
)


def test_chi2_with_step_beyond_the_length_and_overlap_builds_rs_three_intervals():
    """seqdist-CHI2.R shortcuts to one interval only when step == ncol; with a
    larger step and overlap it builds [1, 1+step/2], [1, ncol] and the shifted
    interval. traminex collapsed them to one, giving all zeros here."""
    np.testing.assert_allclose(
        tx.seqdist(FIVE, method="CHI2", step=6, overlap=True),
        _r("chi2_step6_overlap"),
        rtol=1e-12,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        tx.seqdist(FIVE, method="CHI2", step=6), _r("chi2_step6_no_overlap"), atol=1e-12
    )


def test_om_with_an_infinite_indel_and_maxdist_is_nan_off_the_diagonal():
    """R's bound is abs(n - m) * indel + maxscost * min(m, n): 0 * Inf is NaN."""
    expected = _r("om_inf_indel_maxdist")
    result = np.asarray(
        tx.seqdist(FIVE, method="OM", subm="CONSTANT", indel=np.inf, norm="maxdist"), dtype=float
    )
    np.testing.assert_array_equal(np.isnan(result), np.isnan(expected))
    np.testing.assert_array_equal(np.diag(result), np.diag(expected))


def test_omslen_halves_the_cost_before_weighting_it():
    """seqdist.R:737 halves a substitution cost before weighting by the spell
    lengths, so INDELS's largest-double cost does not overflow to Inf."""
    source = tx.seqdef(
        np.array([list("ABBC"), list("AACC"), list("BBAA")], dtype=object), alphabet=["A", "B", "C"]
    )
    reference = tx.seqdef(np.array([["A", "B", None, "C"]], dtype=object), alphabet=["A", "B", "C"])
    result = tx.seqdist(source, method="OMSLEN", subm="INDELS", refseq=reference, with_missing=True)
    np.testing.assert_allclose(result, _r("omslen_indels_missing_ref"), rtol=1e-15, atol=0)


def test_twed_keeps_the_infinite_border_where_rs_stand_in_caps_the_distance():
    """Deliberate divergence. TWEDdistance.cpp fills the impossible border cells
    with ``1000 * maxscost`` (an "arbitrary high value"). With costs of 0.02 and
    stiffness 5 that is 20 -- lower than the real warping cost -- so R returns 20
    for the length-1 sequence against every longer one. traminex keeps the
    TWED definition's infinite border (docs/API_PARITY.md)."""
    ragged = tx.seqdef(
        np.array(
            [
                list("ABCABCA") + [None] * 3,
                ["C"] + [None] * 9,
                ["A", "B"] + [None] * 8,
                ["C"] * 10,
            ],
            dtype=object,
        ),
        alphabet=["A", "B", "C"],
    )
    small = np.array([[0, 0.01, 0.02], [0.01, 0, 0.015], [0.02, 0.015, 0]])

    def twed(nu: float) -> np.ndarray:
        return np.asarray(
            tx.seqdist(ragged, method="TWED", subm=small, nu=nu, h=0.5, norm="none"), dtype=float
        )

    # nu = 1: the stand-in (20) is far above every real warping cost, so it
    # never binds and the two implementations agree exactly.
    np.testing.assert_allclose(twed(1.0), _r("twed_stand_in_nu1"), rtol=1e-12, atol=1e-12)
    # nu = 5: a finite border can only lower a path's cost, so R <= traminex
    # everywhere, strictly where the stand-in binds (every pair of the length-1
    # sequence, where R returns exactly 1000 * maxscost = 20).
    r_values, result = _r("twed_stand_in_nu5"), twed(5.0)
    assert np.all(r_values <= result + 1e-12)
    np.testing.assert_array_equal(r_values[1, [0, 3]], [20.0, 20.0])
    assert np.all(result[1, [0, 3]] > 20.0)
