"""An alphabet state no compared sequence contains changes no distance (TX-AUD-20261019-002).

R golden: ``scripts/r_goldens/make_goldens_indels_unobserved.R`` (``indels_unobserved.csv``).
INDELS gives a state of frequency 0 the stand-in indel ``.Machine$double.xmax`` (seqcost.R).
In R that cost enters the normalization bound and OMloc's local costs although no sequence
contains the state, so every normalized distance is 0 or NaN (TWED's raw distances are Inf).
traminex corrects this (DIVERGENCE row of seqdist): the expected values are R's own results
with that state left out of the alphabet (``observed`` rows), where no stand-in exists.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDEN_ATOL, R_GOLDENS_DIR

import traminex as tx

GOLDEN = pd.read_csv(R_GOLDENS_DIR / "indels_unobserved.csv", dtype={"value": str})
ROWS = [list("AABB"), list("ABAB"), list("BBBA"), list("AAAA")]
EXTRA = {"TWED": {"h": 0.5, "nu": 0.5}, "OMstran": {"otto": 0.5}}
CALLS = sorted(
    {
        (method, norm)
        for method, norm in zip(GOLDEN["method"], GOLDEN["norm"])
        if not str(method).startswith("reference")
    }
)


def _quiet(function, *args, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return function(*args, **kwargs)


def _r(alphabet: str, method: str, norm: str) -> np.ndarray:
    rows = GOLDEN[
        (GOLDEN["alphabet"] == alphabet) & (GOLDEN["method"] == method) & (GOLDEN["norm"] == norm)
    ]
    assert not rows.empty and not rows["value"].str.startswith("ERROR").any(), (
        alphabet,
        method,
        norm,
    )
    return rows["value"].astype(float).to_numpy()


def _upper(matrix) -> np.ndarray:
    values = np.asarray(matrix, dtype=float)
    return values[np.triu_indices(values.shape[0], k=1)]


def _seqdist(alphabet: list[str], method: str, norm: str) -> np.ndarray:
    x = _quiet(tx.seqdef, ROWS, alphabet=alphabet)
    return _upper(
        _quiet(tx.seqdist, x, method=method, subm="INDELS", norm=norm, **EXTRA.get(method, {}))
    )


@pytest.mark.parametrize("method, norm", CALLS)
def test_unobserved_state_changes_no_distance(method: str, norm: str) -> None:
    expected = _r("observed", method, norm)
    np.testing.assert_allclose(
        _seqdist(["A", "B", "C"], method, norm), expected, rtol=1e-12, atol=R_GOLDEN_ATOL
    )


@pytest.mark.parametrize("method, norm", CALLS)
def test_indels_unchanged_when_every_state_is_observed(method: str, norm: str) -> None:
    expected = _r("observed", method, norm)
    np.testing.assert_allclose(
        _seqdist(["A", "B"], method, norm), expected, rtol=1e-12, atol=R_GOLDEN_ATOL
    )


def test_r_defect_is_recorded() -> None:
    # R: every normalized OM, HAM and DHD distance is 0 or NaN with the unobserved state.
    for method in ("OM", "HAM", "DHD"):
        for norm in ("auto", "maxlength", "maxdist", "YujianBo", "gmean"):
            values = _r("with_unobserved", method, norm)
            assert np.all((values == 0) | np.isnan(values)), (method, norm)
            assert np.all(_r("observed", method, norm) > 0), (method, norm)
    assert np.all(np.isinf(_r("with_unobserved", "TWED", "none")))
    # OMstran's raw values change too in R: it rescales its costs by the largest indel
    # (seqdist-OMstran.R:73-77).
    assert not np.allclose(
        _r("with_unobserved", "OMstran", "none"), _r("observed", "OMstran", "none")
    )


def test_costs_returned_by_seqcost_stay_r_s() -> None:
    x = _quiet(tx.seqdef, ROWS, alphabet=["A", "B", "C"])
    costs = _quiet(tx.seqcost, x, "INDELS")
    assert np.asarray(costs["indel"], dtype=float)[2] == np.finfo(float).max
    assert np.all(np.asarray(costs["sm"], dtype=float)[2, :2] >= np.finfo(float).max)


def test_a_state_the_reference_contains_keeps_its_cost() -> None:
    x = _quiet(tx.seqdef, ROWS, alphabet=["A", "B", "C"])
    reference = _quiet(tx.seqdef, [list("CCCC")], alphabet=["A", "B", "C"])
    actual = np.asarray(_quiet(tx.seqdist, x, method="OM", subm="INDELS", refseq=reference), float)
    expected = _r("reference_contains_it", "OM", "none")
    assert np.all(np.isinf(expected))
    np.testing.assert_array_equal(actual.ravel(), expected)


REFERENCE_CALLS = sorted(
    {
        (method, norm)
        for alphabet, method, norm in zip(GOLDEN["alphabet"], GOLDEN["method"], GOLDEN["norm"])
        if alphabet == "external_reference_observed"
    }
)


@pytest.mark.parametrize("alphabet", [["A", "B", "C"], ["A", "B"]])
@pytest.mark.parametrize("method, norm", REFERENCE_CALLS)
def test_an_external_reference_route_is_corrected(
    method: str, norm: str, alphabet: list[str]
) -> None:
    # Batch-7 review I-2: seqdist reaches these distances through its combined-rows route.
    x = _quiet(tx.seqdef, ROWS, alphabet=alphabet)
    reference = _quiet(tx.seqdef, [list("ABBA")], alphabet=alphabet)
    actual = _quiet(tx.seqdist, x, method=method, subm="INDELS", norm=norm, refseq=reference)
    expected = _r("external_reference_observed", method, norm)
    np.testing.assert_allclose(
        np.asarray(actual, float).ravel(), expected, rtol=1e-12, atol=R_GOLDEN_ATOL
    )


@pytest.mark.parametrize("alphabet", [["A", "B", "C"], ["A", "B"]])
@pytest.mark.parametrize("norm", ["none", "maxdist"])
def test_an_index_set_route_is_corrected(norm: str, alphabet: list[str]) -> None:
    x = _quiet(tx.seqdef, ROWS, alphabet=alphabet)
    actual = _quiet(tx.seqdist, x, method="DHD", subm="INDELS", norm=norm, refseq=([1, 2], [3, 4]))
    # R's rows: (1, 3), (2, 3), (1, 4), (2, 4), the 2 x 2 matrix by columns.
    expected = _r("index_sets_observed", "DHD", norm)
    np.testing.assert_allclose(
        np.asarray(actual, float).ravel(order="F"), expected, rtol=1e-12, atol=R_GOLDEN_ATOL
    )


def test_r_defect_is_recorded_on_the_reference_routes() -> None:
    for method in ("OM", "HAM", "DHD"):
        assert np.all(_r("external_reference_with_unobserved", method, "maxdist") == 0), method
    assert np.all(_r("index_sets_with_unobserved", "DHD", "maxdist") == 0)


@pytest.mark.parametrize("norm", ["maxlength", "maxdist"])
def test_dhd_bound_uses_the_states_present_at_each_position(norm: str) -> None:
    """Worked by hand. Sequences (A,A), (B,A), (A,C): position 1 holds A, A, B (frequencies
    2/3, 1/3), position 2 holds A, A, C. INDELS indels are 1/frequency: 1.5 and 3, so
    cost(A,B) at 1 and cost(A,C) at 2 are 4.5; C at 1 and B at 2 have the stand-in. Raw DHD:
    d12 = 4.5, d13 = 4.5, d23 = 9. The bound sums the largest cost among the states present
    at each position: 4.5 + 4.5 = 9, so the normalized distances are 0.5, 0.5 and 1."""
    x = _quiet(tx.seqdef, [["A", "A"], ["B", "A"], ["A", "C"]], alphabet=["A", "B", "C"])
    raw = _upper(_quiet(tx.seqdist, x, method="DHD", subm="INDELS"))
    np.testing.assert_allclose(raw, [4.5, 4.5, 9.0], rtol=1e-12)
    normalized = _upper(_quiet(tx.seqdist, x, method="DHD", subm="INDELS", norm=norm))
    np.testing.assert_allclose(normalized, [0.5, 0.5, 1.0], rtol=1e-12)
