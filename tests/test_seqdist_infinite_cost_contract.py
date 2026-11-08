"""An infinite substitution cost must not silently produce an all-inf matrix (F-019).

``seqdist`` derives an automatic indel from an explicit cost matrix as R does --
``SOURCES/TraMineR/R/seqdist.R:359-365``: ``indel <- max(sm)/2`` over the **whole**
matrix, not the observed sub-matrix.  That is load-bearing and is pinned below:
changing only the cost between two *unobserved* states changes R's normalized
distances (``norm="maxdist"``: 1/3 for an unobserved cost of 2, and 0.04 /
0.0267 / 0.0133 for 50), so restricting the maximum to observed states would
break parity.

The consequence used to be that an ``Inf`` in an unused cell made the automatic
indel infinite and every distance ``inf``, with only numpy overflow warnings.
R refuses the same three inputs (``Error: missing value where TRUE/FALSE
needed``, raised inside ``checkcost`` once ``tol <- 1e-07 * max(sma)`` is ``Inf``).
Every R result here is TraMineR 2.2-14's, recorded by
``scripts/r_goldens/make_goldens_audit20260922_s10.R`` in
``tests/r_goldens/inline_r_anchors.txt``.
traminex need not reproduce R's internal error text, but it must not silently
return a meaningless all-infinite matrix, so it refuses with a message that names
the cause.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest
from r_goldens_utils import r_anchor, r_anchor_values

import traminex as tx
from traminex.distances import pairwise

_RAW = np.array([["A", "B", "A"], ["B", "A", "B"], ["A", "A", "B"]], dtype=object)
_ALPHABET = ["A", "B", "C", "D", "E"]  # D and E are never observed


def _sequences():
    return tx.seqdef(_RAW.copy(), alphabet=_ALPHABET)


def _costs(row: int, column: int, value: float) -> np.ndarray:
    costs = np.full((5, 5), 2.0)
    np.fill_diagonal(costs, 0.0)
    costs[row, column] = costs[column, row] = value
    return costs


@pytest.mark.parametrize(
    ("tag", "row", "column"),
    [
        ("both observed", 0, 1),
        ("both unobserved", 3, 4),
        ("observed to unobserved", 2, 3),
    ],
)
def test_an_infinite_substitution_cost_is_refused(tag, row, column):
    r_error = r_anchor(f"infinite_cost_{row + 1}_{column + 1}", "error")
    assert r_error == ["missing value where TRUE/FALSE needed"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with pytest.raises(ValueError, match="finite"):
            tx.seqdist(_sequences(), method="OM", subm=_costs(row, column, np.inf))


def test_an_infinite_cost_no_observed_state_reaches_is_computed_with_a_given_indel():
    """R stops here too; traminex computes, since no distance can reach the cell.

    With the indel given, nothing is derived from the infinite cost, and the
    distances equal those with that cell finite (docs/API_PARITY.md; audit
    S1a-13, which found the documentation claiming this input is refused).
    """
    r_error = r_anchor("infinite_unobserved_cost_given_indel", "error")
    assert r_error == ["missing value where TRUE/FALSE needed"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        infinite = tx.seqdist(_sequences(), method="OM", subm=_costs(3, 4, np.inf), indel=1.0)
        finite = tx.seqdist(_sequences(), method="OM", subm=_costs(3, 4, 2.0), indel=1.0)
    np.testing.assert_array_equal(infinite, finite)
    assert np.all(np.isfinite(infinite))


@pytest.mark.parametrize("unobserved_cost", [2.0, 50.0])
def test_the_automatic_indel_still_uses_the_whole_matrix_as_r_does(unobserved_cost):
    """TraMineR 2.2-14's d(1,2), d(1,3), d(2,3)."""
    expected = r_anchor_values(f"unobserved_cost_{unobserved_cost:g}")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        distances = np.asarray(
            tx.seqdist(
                _sequences(),
                method="OM",
                subm=_costs(3, 4, unobserved_cost),
                norm="maxdist",
            )
        )
    actual = [distances[0, 1], distances[0, 2], distances[1, 2]]
    np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-12)


def test_the_indels_route_with_an_external_reference_keeps_its_infinite_distances():
    """This legitimate infinity must survive.

    TraMineR 2.2-13 returned ``Inf Inf Inf`` here. 2.2-14 aborts inside
    ``checkcost`` instead and returns ``Inf Inf Inf`` only with
    ``check.costs = FALSE``; traminex keeps computing, the documented INDELS
    divergence (docs/API_PARITY.md).
    """
    assert r_anchor("indels_refseq_unobserved", "error") == [
        "missing value where TRUE/FALSE needed"
    ]
    r_without_check = r_anchor_values("indels_refseq_unobserved_no_check")
    external = tx.seqdef(np.array([["D", "D", "D"]], dtype=object), alphabet=_ALPHABET)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        distances = np.asarray(
            tx.seqdist(_sequences(), method="OM", subm="INDELS", refseq=external)
        )
    assert distances.shape == r_without_check.shape == (3,)
    assert np.all(np.isinf(r_without_check))
    np.testing.assert_array_equal(distances, r_without_check)


def test_the_triangle_inequality_scan_does_not_leak_a_numpy_warning():
    """A non-finite cost matrix must not print ``invalid value encountered``.

    ``_first_triangle_inequality_violation`` reproduces R's ``checktriangleineq``
    C loop, where ``Inf - Inf`` is a silent NaN. In NumPy the same subtraction
    raises ``RuntimeWarning: invalid value encountered in scalar subtract``,
    which reached the user before the (correct) R-ordered diagnostic
    (2026-09-20 T-07, cycle-10 N-18). The analytical result was never affected.
    """
    import warnings

    import numpy as np

    from traminex.distances.python_kernels import _first_triangle_inequality_violation

    costs = np.array([[0.0, np.inf, 2.0], [np.inf, 0.0, 2.0], [2.0, 2.0, 0.0]])
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        _first_triangle_inequality_violation(costs, tol=1e-7)


_VOID = 9
_RAGGED = np.array([[0, 1, 2, _VOID], [1, 0, _VOID, _VOID], [2, 2, 1, 0]], dtype=np.int32)
_DENSE = np.array([[0, 1, 2], [1, 0, 0], [2, 2, 1]], dtype=np.int32)
_ROUTES = {
    "void": lambda subm: pairwise(_RAGGED, "OM", subm=subm, indel=1.0, void=_VOID),
    # One indel per active state: the rows hold the three codes 0, 1, 2.
    "indel_vector": lambda subm: pairwise(_DENSE, "OM", subm=subm, indel=np.array([1.0, 1.0, 1.0])),
    "omslen": lambda subm: pairwise(_DENSE, "OMslen", subm=subm, indel=1.0),
    "twed": lambda subm: pairwise(_DENSE, "TWED", subm=subm, nu=0.5, h=0.5),
}


def _with_cost(size: int, row: int, column: int, value: float) -> np.ndarray:
    costs = np.ones((size, size))
    np.fill_diagonal(costs, 0.0)
    costs[row, column] = costs[column, row] = value
    return costs


@pytest.mark.parametrize("value", [np.nan, np.inf], ids=["nan", "inf"])
@pytest.mark.parametrize("route", sorted(_ROUTES))
def test_every_low_level_route_refuses_a_reachable_non_finite_cost(route: str, value: float):
    # TX-AUD-20261019-010: the dense OM route refused it, these routes computed silently; R stops
    # on every such matrix with only "missing value where TRUE/FALSE needed".
    with pytest.raises(ValueError, match="substitution costs must be finite"):
        _ROUTES[route](_with_cost(3, 0, 2, value))


@pytest.mark.parametrize("route", sorted(_ROUTES))
def test_a_non_finite_cost_no_row_reaches_is_still_accepted(route: str):
    # State 3 occurs in no row, so no distance can reach its cost (the documented contract).
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = np.asarray(_ROUTES[route](_with_cost(4, 0, 3, np.nan)), dtype=float)
    assert np.all(np.isfinite(result))
