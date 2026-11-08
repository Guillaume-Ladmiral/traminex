"""traminex must emit TraMineR 2.2-14's non-fatal cost diagnostics, word for word.

``SOURCES/TraMineR/R/TraMineR-checkcost.R`` (called from ``seqdist.R:380-390``)
checks a caller-supplied substitution-cost matrix, or a DHD time-varying array,
and reports five non-fatal conditions: a largest cost of 0 or below (after
which it checks nothing else), zero costs between distinct states, a
substitution triangle-inequality violation, an indel triangle-inequality
violation or a cost above twice the indel, and asymmetry. The first four go
through ``message()``; only asymmetry is a ``warning()``.

traminex had ported the TraMineR 2.2-13 version: the twice-the-indel and
indel-triangle texts were 2.2-13's and went out as ordinary warnings, the
triangle message lacked 2.2-14's layout, asymmetry read "The substitution cost
matrix is not symmetric", the zero-cost and ``max(sm) <= 0`` diagnostics were
missing (an all-zero matrix instead produced a spurious triangle message),
states were named by integer code, and time-varying arrays got no diagnostic
but the asymmetry one. This module's earlier version pinned those 2.2-13 texts
while its docstring said they had been confirmed against 2.2-14 (2026-09-22
audit S10-01, S1a-11).

Expected diagnostics: live TraMineR 2.2-14 via
``scripts/r_goldens/make_goldens_audit20260922_s10.R``
(``tests/r_goldens/checkcost_diagnostics.txt``). R's ``message()`` diagnostics
arrive as :class:`~traminex.distances.python_kernels.CostDiagnosticMessage`
warnings, its ``warning()`` as ``UserWarning``; a Python warning carries no
equivalent of the newline ``message()`` appends, so that one newline is removed
before comparing. Everything else is compared exactly, including order.
"""

from __future__ import annotations

import csv
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pytest

import traminex as tx
from traminex.distances.python_kernels import CostDiagnosticMessage

FIXTURE = Path(__file__).parent / "r_goldens" / "checkcost_diagnostics.txt"

_RAW = [list("ABCABC"), list("BACCAB"), list("CCABBA")]
_RAW_MISSING = [["A", "B", None, "A"], ["B", "A", "C", "C"], ["C", None, "A", "B"]]
_RAW_SEQDIFF = _RAW + [list("AABBCC"), list("BBBAAC"), list("CABCAB")]


def _matrix(*values: float) -> np.ndarray:
    return np.array(values, dtype=float).reshape(3, 3)


SYM = _matrix(0, 2, 2, 2, 0, 2, 2, 2, 0)
GT = _matrix(0, 5, 5, 5, 0, 5, 5, 5, 0)
TRI = _matrix(0, 10, 2, 10, 0, 2, 2, 2, 0)
TRIF = _matrix(0, 10.5, 1 / 3, 10.5, 0, 2.25, 1 / 3, 2.25, 0)
ASY = _matrix(0, 3, 2, 1, 0, 2, 2, 2, 0)
ZAB = _matrix(0, 0, 2, 0, 0, 2, 2, 2, 0)
ZBC = _matrix(0, 2, 2, 2, 0, 0, 2, 0, 0)
ZTWO = _matrix(0, 0, 0, 0, 0, 2, 0, 2, 0)
ZALL = np.zeros((3, 3))


def _with_missing_costs(row: int, col: int, value: float) -> np.ndarray:
    costs = np.full((4, 4), 2.0)
    np.fill_diagonal(costs, 0.0)
    costs[row, col] = costs[col, row] = value
    return costs


def _dhd_array(*edits: tuple[int, int, int, float], symmetric: bool = True) -> np.ndarray:
    """A (time, state, state) array: R's ``array(2, c(3, 3, 6))`` time-first.

    Each edit is R's one-based ``(row, col, position)`` and a value.
    """
    costs = np.full((6, 3, 3), 2.0)
    for time in range(6):
        np.fill_diagonal(costs[time], 0.0)
    for row, col, position, value in edits:
        costs[position - 1, row - 1, col - 1] = value
        if symmetric:
            costs[position - 1, col - 1, row - 1] = value
    return costs


def _sequences():
    return tx.seqdef(np.array(_RAW, dtype=object), alphabet=["A", "B", "C"])


def _seqdist(method: str, subm, **options):
    return lambda: tx.seqdist(_sequences(), method=method, subm=subm, **options)


def _seqdist_missing(subm, indel: float):
    def call():
        data = tx.seqdef(np.array(_RAW_MISSING, dtype=object), alphabet=["A", "B", "C"])
        return tx.seqdist(data, method="OM", subm=subm, indel=indel, with_missing=True)

    return call


def _seqdiff(subm, indel):
    def call():
        data = tx.seqdef(np.array(_RAW_SEQDIFF, dtype=object), alphabet=["A", "B", "C"])
        return tx.seqdiff(
            data,
            np.array([1, 1, 1, 2, 2, 2]),
            cmprange=(0, 1),
            seqdist_args={"method": "OM", "sm": subm, "indel": indel},
        )

    return call


TWED_OPTIONS = {"indel": 6.0, "nu": 0.5, "h": 0.5}
CASES = {
    "om_twice_indel": _seqdist("OM", GT, indel=1.0),
    "om_twice_indel_auto": _seqdist("OM", GT),
    "om_triangle": _seqdist("OM", TRI, indel=6.0),
    "om_triangle_fractional": _seqdist("OM", TRIF, indel=6.0),
    "om_indel_vector_triangle": _seqdist("OM", SYM, indel=np.array([0.1, 5.0, 5.0])),
    "om_indel_vector_no_twice": _seqdist("OM", GT, indel=np.array([1.0, 1.0, 1.0])),
    "om_asymmetric": _seqdist("OM", ASY, indel=2.0),
    "om_well_formed": _seqdist("OM", SYM, indel=1.0),
    "om_zero_ab": _seqdist("OM", ZAB, indel=1.0),
    "om_zero_bc": _seqdist("OM", ZBC, indel=1.0),
    "om_zero_two_pairs": _seqdist("OM", ZTWO, indel=1.0),
    "om_all_zero": _seqdist("OM", ZALL, indel=1.0),
    "omloc_triangle": _seqdist("OMloc", TRI),
    "omloc_zero_ab": _seqdist("OMloc", ZAB),
    "omslen_triangle": _seqdist("OMslen", TRI, indel=6.0),
    "omslen_zero_ab": _seqdist("OMslen", ZAB, indel=6.0),
    "omspell_triangle": _seqdist("OMspell", TRI, indel=6.0),
    "omspell_zero_ab": _seqdist("OMspell", ZAB, indel=6.0),
    "omstran_triangle": _seqdist("OMstran", TRI, indel=6.0, otto=0.5),
    "omstran_zero_ab": _seqdist("OMstran", ZAB, indel=6.0, otto=0.5),
    "twed_triangle": _seqdist("TWED", TRI, **TWED_OPTIONS),
    "twed_zero_ab": _seqdist("TWED", ZAB, **TWED_OPTIONS),
    "ham_triangle": _seqdist("HAM", TRI),
    "ham_zero_ab": _seqdist("HAM", ZAB),
    "om_missing_zero": _seqdist_missing(_with_missing_costs(2, 3, 0.0), indel=1.0),
    "om_missing_triangle": _seqdist_missing(_with_missing_costs(0, 3, 9.0), indel=5.0),
    "dhd_zero": _seqdist("DHD", _dhd_array((1, 2, 3, 0.0), (1, 3, 5, 0.0))),
    "dhd_triangle": _seqdist("DHD", _dhd_array((1, 2, 4, 10.0), (1, 2, 6, 9.0))),
    "dhd_asymmetric": _seqdist("DHD", _dhd_array((1, 2, 2, 3.0), (1, 2, 5, 3.0), symmetric=False)),
    "dhd_all_zero": _seqdist("DHD", np.zeros((6, 3, 3))),
    "seqdiff_twice_indel": _seqdiff(GT, 1.0),
    "seqdiff_indel_vector_triangle": _seqdiff(SYM, [0.1, 5.0, 5.0]),
    "seqdiff_asymmetric": _seqdiff(ASY, 2.0),
}

#: R never runs checkcost for HAM (a dangling ``else`` in seqdist.R:380-386);
#: traminex deliberately gives HAM the triangle-inequality diagnostic
#: (docs/API_PARITY.md), so these two are asserted separately below.
HAM_CASES = {"ham_triangle", "ham_zero_ab"}

_NEGATIVE = np.full((3, 3), -1.0)
np.fill_diagonal(_NEGATIVE, 0.0)
#: ``max(sm) <= 0`` makes R's checkcost return before its refusals, so R
#: computes with these; traminex refuses them (docs/API_PARITY.md).
REFUSED_CASES = {
    "om_negative_max_zero": _seqdist("OM", _NEGATIVE, indel=1.0),
    "om_wrong_size_all_zero": _seqdist("OM", np.zeros((2, 2)), indel=1.0),
}


def _r_diagnostics() -> dict[str, list[tuple[str, str]]]:
    recorded: dict[str, list[tuple[str, str]]] = defaultdict(list)
    with FIXTURE.open(newline="") as handle:
        rows = csv.reader((line for line in handle if not line.startswith("#")), delimiter="\t")
        for case, _order, channel, text in rows:
            if channel == "none":
                recorded[case] = []
                continue
            decoded = text.replace("\\n", "\n")
            if channel == "message":
                assert decoded.endswith("\n"), (case, decoded)
                decoded = decoded[:-1]
            recorded[case].append((channel, decoded))
    return dict(recorded)


R = _r_diagnostics()


def _traminex_diagnostics(call) -> list[tuple[str, str]]:
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        call()
    out = []
    for item in caught:
        if issubclass(item.category, CostDiagnosticMessage):
            channel = "message"
        elif item.category is UserWarning:
            channel = "warning"
        else:
            channel = item.category.__name__
        out.append((channel, str(item.message)))
    return out


def test_every_recorded_case_is_exercised():
    assert set(CASES) | set(REFUSED_CASES) == set(R)


@pytest.mark.parametrize("case", sorted(set(CASES) - HAM_CASES))
def test_diagnostics_match_traminer_2_2_14(case):
    assert _traminex_diagnostics(CASES[case]) == R[case]


def test_ham_gets_only_the_triangle_diagnostic_r_never_runs():
    """R emits nothing for HAM; traminex's documented choice emits the triangle one."""
    assert R["ham_triangle"] == [] and R["ham_zero_ab"] == []
    assert _traminex_diagnostics(CASES["ham_triangle"]) == R["om_triangle"]
    assert _traminex_diagnostics(CASES["ham_zero_ab"]) == []


@pytest.mark.parametrize("case", sorted(REFUSED_CASES))
def test_a_non_positive_maximum_does_not_bypass_the_refusals(case):
    """R reports the maximum and computes anyway; traminex refuses (documented)."""
    assert R[case] == [("message", " [!!] max substitution cost less or equal to 0!")]
    with pytest.raises(ValueError, match="negative|must be 3x3"):
        REFUSED_CASES[case]()


@pytest.mark.parametrize(
    "case",
    ["om_twice_indel", "om_triangle", "om_indel_vector_triangle", "om_zero_ab", "om_all_zero"],
)
def test_diagnostics_are_not_fatal(case):
    """Every one of them still returns the distance matrix."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = CASES[case]()
    assert result.shape == (3, 3)
    assert np.all(np.isfinite(result))
