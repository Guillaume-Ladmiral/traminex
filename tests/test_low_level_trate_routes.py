"""OM's low-level subm="TRATE" works on every pairwise route (TX-AUD-20261019-024).

``docs/DISTANCES.md`` documents ``pairwise(seqmat, method="OM", subm="TRATE", ...)``. The dense route resolved the
string, but the void route and a per-state indel vector branched off first and refused it with the message meant
for the other methods ("cost-method strings are resolved by seqdist"). Every route now resolves OM's TRATE costs,
from the rows without their void cells, and returns what seqdist returns for the same data.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

import traminex as tx
from traminex.distances import pairwise

DENSE = tx.seqdef(
    np.array([list("AABBCC"), list("ABBCCD"), list("CCBAAD"), list("DDCBAA")], dtype=object),
    alphabet=["A", "B", "C", "D"],
)
#: Ragged rows: dropping the void cells changes the transitions counted (Review Focus 2).
RAGGED = tx.seqdef(
    np.array(
        [
            ["A", "A", "B", "C", None, None],
            ["A", "B", "B", "C", "C", "D"],
            ["C", "C", "B", "A", None, None],
            ["D", "C", "B", "A", "A", None],
        ],
        dtype=object,
    ),
    alphabet=["A", "B", "C", "D"],
    right="DEL",
)


def _quiet(function, *args, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return np.asarray(function(*args, **kwargs), dtype=float)


def test_dense_route_equals_seqdist() -> None:
    np.testing.assert_allclose(
        _quiet(pairwise, np.asarray(DENSE.data), method="OM", subm="TRATE", indel=1.0),
        _quiet(tx.seqdist, DENSE, method="OM", subm="TRATE", indel=1.0),
        rtol=0,
        atol=1e-12,
    )


def test_void_route_equals_seqdist() -> None:
    assert RAGGED.void_code is not None and np.any(RAGGED.data == RAGGED.void_code)
    np.testing.assert_allclose(
        _quiet(
            pairwise,
            np.asarray(RAGGED.data),
            method="OM",
            subm="TRATE",
            indel=1.0,
            void=RAGGED.void_code,
        ),
        _quiet(tx.seqdist, RAGGED, method="OM", subm="TRATE", indel=1.0),
        rtol=0,
        atol=1e-12,
    )


def test_indel_vector_route_equals_seqdist() -> None:
    indel = np.array([1.0, 2.0, 1.5, 1.0])
    np.testing.assert_allclose(
        _quiet(pairwise, np.asarray(DENSE.data), method="OM", subm="TRATE", indel=indel),
        _quiet(tx.seqdist, DENSE, method="OM", subm="TRATE", indel=indel),
        rtol=0,
        atol=1e-12,
    )


@pytest.mark.parametrize("method", ["OMSLEN", "OMSPELL", "TWED"])
def test_other_methods_still_leave_cost_strings_to_seqdist(method: str) -> None:
    extra = {"nu": 0.5, "h": 0.5} if method == "TWED" else {}
    with pytest.raises(ValueError, match="resolved by seqdist"):
        pairwise(np.asarray(DENSE.data), method=method, subm="TRATE", **extra)


#: One position: no transition, so no transition rates (batch-9 review M-3).
WIDTH_1 = tx.seqdef(np.array([["A"], ["B"], ["A"]], dtype=object), alphabet=["A", "B"])
#: Stored width 2, one real cell per row: seqdist and live R compute cost-2 distances
#: (evidence/batch9/fixpass/m3_void_width_{py,r}.out), so the void route does too.
VOID_WIDTH_2 = tx.seqdef(
    np.array([["A", None], ["B", None], ["A", None]], dtype=object),
    alphabet=["A", "B"],
    right="DEL",
)


@pytest.mark.parametrize("indel", [1.0, np.array([1.0, 1.0])], ids=["dense", "indel_vector"])
def test_every_raw_route_applies_the_formula_without_a_transition(indel) -> None:
    """``docs/DISTANCES.md`` documents the raw-array extension: without a transition the TRATE formula
    still applies (distinct states cost 2), where seqdist refuses as R does. Every low-level route keeps
    that contract (batch-9 review M-3, ruled not a defect)."""
    with pytest.raises(ValueError, match="absolute lag must be smaller than the sequence width"):
        _quiet(tx.seqdist, WIDTH_1, method="OM", subm="TRATE", indel=indel)
    np.testing.assert_array_equal(
        _quiet(pairwise, np.asarray(WIDTH_1.data), method="OM", subm="TRATE", indel=indel),
        [[0.0, 2.0, 0.0], [2.0, 0.0, 2.0], [0.0, 2.0, 0.0]],
    )


def test_void_route_on_a_wider_stored_width_computes_as_seqdist_does() -> None:
    np.testing.assert_allclose(
        _quiet(
            pairwise,
            np.asarray(VOID_WIDTH_2.data),
            method="OM",
            subm="TRATE",
            indel=1.0,
            void=VOID_WIDTH_2.void_code,
        ),
        _quiet(tx.seqdist, VOID_WIDTH_2, method="OM", subm="TRATE", indel=1.0),
        rtol=0,
        atol=1e-12,
    )
