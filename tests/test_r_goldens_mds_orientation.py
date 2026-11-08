"""The one-dimensional MDS axis has a documented orientation (TX-REP-20261017-006).

seqrf, dissrf and seqplot's rf groups sort the sequences by the first MDS
coordinate, whose sign is arbitrary: R takes it from LAPACK, and on these 14
sequences R (OpenBLAS 0.3.26) gives the negation of traminex (SciPy's LAPACK
3.9.0), so the groups, medoids and R2 differed (k = 3: R2 0.2 against R's
0.2667). traminex now orients the axis by a rule -- the coordinate of largest
absolute value is positive; coordinates within a relative 1e-10 of it are tied
and the first decides; missing (zero-weight) coordinates are ignored -- and the
expected values are R's seqrf given the axis oriented by that rule. The
tolerance matters on symmetric data, where two coordinates are equal in exact
arithmetic and differ only in rounding (line5, line7, mirror6; batch-5 review
I-3). On this build SciPy's sign already equals the rule's for the 14
sequences, so the seqrf comparisons also run with the eigenvectors negated,
which a missing orientation fails (batch-5 review M-6).
Fixture: ``scripts/r_goldens/make_goldens_mds_orientation.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx
import traminex.repseq_rf as rf_module

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "mds_orientation.csv",
    dtype=str,
    keep_default_na=False,
)
STRINGS = ["AAABBC", "AABBCC", "ABBBCC", "AAAABB", "BBBCCC", "ABCABC", "AACCBB",
           "CCBBAA", "ABBBBB", "AAACCC", "BCCCCC", "AAAAAC", "BBAACC", "CABBBA"]  # fmt: skip
WEIGHTS = np.array([1, 2, 1, 1, 0.5, 1, 1, 2, 1, 1, 1, 1, 0, 1], dtype=float)
ROWS = np.array([list(s) for s in STRINGS], dtype=object)
PLAIN = tx.seqdef(ROWS, alphabet=["A", "B", "C"])
WEIGHTED = tx.seqdef(ROWS, alphabet=["A", "B", "C"], weights=WEIGHTS)
DISS = np.asarray(tx.seqdist(PLAIN, method="LCS"), dtype=float)
SYMMETRIC = {
    "line5": ["AAAA", "AAAB", "AABB", "ABBB", "BBBB"],
    "line7": ["AAAAAA", "AAAAAB", "AAAABB", "AAABBB", "AABBBB", "ABBBBB", "BBBBBB"],
    "mirror6": ["AAAB", "BAAA", "AABB", "BBAA", "ABBB", "BBBA"],
}
SYMMETRIC_SEQS = {
    name: tx.seqdef(np.array([list(s) for s in strings], dtype=object), alphabet=["A", "B"])
    for name, strings in SYMMETRIC.items()
}
SYMMETRIC_DISS = {
    name: np.asarray(tx.seqdist(seqs, method="HAM"), dtype=float)
    for name, seqs in SYMMETRIC_SEQS.items()
}
RUNS = {
    "unweighted_first_k3": (PLAIN, DISS, 3, "first"),
    "unweighted_first_k4": (PLAIN, DISS, 4, "first"),
    "unweighted_prop_k3": (PLAIN, DISS, 3, "prop"),
    "unweighted_prop_k4": (PLAIN, DISS, 4, "prop"),
    "weighted_prop_k3": (WEIGHTED, DISS, 3, "prop"),
    **{
        f"{name}_first_k2": (SYMMETRIC_SEQS[name], SYMMETRIC_DISS[name], 2, "first")
        for name in SYMMETRIC
    },
}


def _r(case: str, quantity: str) -> np.ndarray:
    rows = GOLDEN[(GOLDEN["case"] == case) & (GOLDEN["quantity"] == quantity)]
    values = rows.sort_values("position", key=lambda c: c.astype(int))["value"].str.strip()
    return np.array([np.nan if value == "NA" else float(value) for value in values])


@pytest.mark.parametrize("library_sign", [1, -1])
def test_the_unweighted_axis_follows_the_rule(library_sign: int) -> None:
    """The oriented axis is R's raw axis up to its sign, whichever sign R's LAPACK gives.

    The raw sign follows the numerical library (generator header), so asserting it
    failed on an oracle whose LAPACK returns the other sign; ``library_sign`` -1 stands
    for one (TX-REP-20261019-006, batch-5 review M-8).
    """
    r_axis = library_sign * _r("unweighted", "r_axis")
    oriented = _r("unweighted", "oriented_axis")
    assert np.allclose(oriented, r_axis, rtol=0, atol=1e-12) or np.allclose(
        oriented, -r_axis, rtol=0, atol=1e-12
    )
    np.testing.assert_allclose(
        rf_module._cmdscale_1d(DISS), _r("unweighted", "oriented_axis"), rtol=1e-9, atol=1e-12
    )


def test_the_weighted_axis_follows_the_rule() -> None:
    weights, _ = rf_module._normalize_weights(len(WEIGHTS), WEIGHTS)
    np.testing.assert_allclose(
        rf_module._weighted_cmdscale_1d(DISS, weights),
        _r("weighted", "oriented_axis"),
        rtol=1e-9,
        atol=1e-12,
    )


@pytest.mark.parametrize("name", sorted(SYMMETRIC))
def test_symmetric_axes_follow_the_rule(name: str) -> None:
    np.testing.assert_allclose(
        rf_module._cmdscale_1d(SYMMETRIC_DISS[name]),
        _r(name, "oriented_axis"),
        rtol=1e-9,
        atol=1e-12,
    )


def _assert_seqrf_equals_r(case: str) -> None:
    seqs, diss, k, grp_meth = RUNS[case]
    result = tx.seqrf(seqs, diss, k=k, grp_meth=grp_meth)["rf"]
    assert [int(m) + 1 for m in result["medoids"]] == [int(m) for m in _r(case, "medoids")]
    np.testing.assert_allclose(result["heights"], _r(case, "heights"), rtol=1e-12)
    assert result["R2"] == pytest.approx(_r(case, "R2")[0], rel=1e-9)
    assert result["Fstat"] == pytest.approx(_r(case, "Fstat")[0], rel=1e-9)


def _negate_eigenvectors(monkeypatch) -> None:
    real = rf_module.linalg.eigh

    def negated(*args, **kwargs):
        values, vectors = real(*args, **kwargs)
        return values, -vectors

    monkeypatch.setattr(rf_module.linalg, "eigh", negated)


@pytest.mark.parametrize("case", sorted(RUNS))
def test_default_seqrf_equals_r_on_the_oriented_axis(case: str) -> None:
    _assert_seqrf_equals_r(case)


@pytest.mark.parametrize("case", sorted(RUNS))
def test_default_seqrf_equals_r_with_the_eigenvectors_negated(case: str, monkeypatch) -> None:
    _negate_eigenvectors(monkeypatch)
    _assert_seqrf_equals_r(case)


def test_the_orientation_does_not_depend_on_the_eigenvector_sign(monkeypatch) -> None:
    expected = rf_module._cmdscale_1d(DISS)
    _negate_eigenvectors(monkeypatch)
    np.testing.assert_allclose(rf_module._cmdscale_1d(DISS), expected, rtol=0, atol=0)


def test_ties_and_missing_coordinates() -> None:
    np.testing.assert_array_equal(
        rf_module._oriented(np.array([-2.0, 1.0, 2.0])), [2.0, -1.0, -2.0]
    )
    np.testing.assert_array_equal(
        rf_module._oriented(np.array([np.nan, -3.0, 1.0])), [np.nan, 3.0, -1.0]
    )


def test_coordinates_within_the_tolerance_are_tied() -> None:
    # Equal in exact arithmetic, apart in the last bits: the first decides.
    np.testing.assert_array_equal(
        rf_module._oriented(np.array([2.0 - 4e-16, 0.0, -2.0])), [2.0 - 4e-16, 0.0, -2.0]
    )
    np.testing.assert_array_equal(
        rf_module._oriented(np.array([-(2.0 - 4e-16), 0.0, 2.0])), [2.0 - 4e-16, 0.0, -2.0]
    )
    # Outside the tolerance the larger magnitude still decides.
    np.testing.assert_array_equal(
        rf_module._oriented(np.array([2.0 - 1e-6, 0.0, -2.0])), [-(2.0 - 1e-6), 0.0, 2.0]
    )
