"""seqMD CAT costs when a combined state has a void component (audit S2-19).

A combined state can carry a void component: the padding of a shorter
channel, or a trailing void in a channel with no missing values. R builds
each CAT cost by matching every component in its channel's cost alphabet
(``seqMD.R:289-335``); a void component has no entry, so ``match()`` gives NA
and every cost pairing that state with another is NA, while the diagonal
keeps its initial 0 and a scalar indel is unaffected. ``what="cost"``
returns that matrix; ``what="diss"`` stops with "NA indel and/or substitution
CAT costs, cannot compute MD distances!" (``seqMD.R:357-359``).

traminex raised "combined state component '%' is not in channel 1 costs"
for both, so ``what="cost"`` refused a matrix R returns. It now returns R's
matrix with NaN for NA, and ``what="diss"`` raises R's message.

Source of truth: live TraMineR 2.2-14; values in
``tests/r_goldens/audit20260922_s2b_cases.txt``
(``scripts/r_goldens/make_goldens_audit20260922_s2b.R``).
"""

from __future__ import annotations

import ctypes.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from traminex import STSList, seqcost, seqdef, seqMD

pytestmark = [
    pytest.mark.skipif(
        ctypes.util.find_library("icui18n") is None,
        reason="the combined alphabet order needs the system ICU library",
    ),
    # Every fixture has cases of unequal effective length, which seqMD (like
    # R) reports with this warning.
    pytest.mark.filterwarnings("ignore:Cases with sequences of different length"),
]

FIXTURE = Path(__file__).resolve().parent / "r_goldens" / "audit20260922_s2b_cases.txt"


def _r(case: str) -> dict[str, str]:
    fields: dict[str, str] = {}
    for line in FIXTURE.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        name, field, value = line.split("\t")
        if name == case:
            fields[field] = value
    assert fields, case
    return fields


def _numbers(text: str) -> np.ndarray:
    return np.array([np.nan if value == "NA" else float(value) for value in text.split("|")])


def _assert_cost_matches_r(costs: dict, case: str) -> None:
    expected = _r(case)
    assert costs["alphabet"] == expected["alphabet"].split("|")
    shape = tuple(int(size) for size in expected["dim"].split("|"))
    np.testing.assert_allclose(
        np.asarray(costs["sm"], dtype=float),
        _numbers(expected["sm"]).reshape(shape, order="F"),
        rtol=0,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        np.atleast_1d(np.asarray(costs["indel"], dtype=float)),
        _numbers(expected["indel"]),
        rtol=0,
        atol=1e-12,
    )


def _unequal_channels() -> list[STSList]:
    shorter = seqdef(
        pd.DataFrame({"a": ["A", "B"], "b": [None, "A"], "c": ["B", "B"]}),
        alphabet=["A", "B"],
        right=None,
    )
    longer = seqdef(
        pd.DataFrame(
            {"a": ["x", "y"], "b": ["y", "y"], "c": ["x", "y"], "d": ["y", "x"], "e": ["x", None]}
        ),
        alphabet=["x", "y"],
    )
    return [shorter, longer]


def _void_channels() -> list[STSList]:
    first = seqdef(
        pd.DataFrame({"a": ["A", "B"], "b": ["A", "B"], "c": ["B", None]}),
        alphabet=["A", "B"],
    )
    second = seqdef(
        pd.DataFrame({"a": ["x", "y"], "b": ["y", "y"], "c": ["x", "y"]}),
        alphabet=["x", "y"],
    )
    return [first, second]


def test_unequal_widths_cost_matches_r():
    _assert_cost_matches_r(
        seqMD(_unequal_channels(), sm="CONSTANT", what="cost"), "unequal_constant"
    )


def test_void_component_cost_matches_r():
    _assert_cost_matches_r(seqMD(_void_channels(), sm="CONSTANT", what="cost"), "void_constant")


def test_void_component_cost_with_mean_link_and_weights_matches_r():
    costs = seqMD(_void_channels(), sm="CONSTANT", what="cost", link="mean", cweight=[2, 0.5])
    _assert_cost_matches_r(costs, "void_mean_cweight")


def test_void_component_time_varying_cost_matches_r():
    channels = _void_channels()
    matrices = [np.asarray(seqcost(c, "TRATE", time_varying=True)["sm"]) for c in channels]
    _assert_cost_matches_r(seqMD(channels, sm=matrices, what="cost"), "void_trate_time_varying")


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ \ \ \ \ \ for\ example,\ cost\ is\ zero\ for\ the\ pair\ \(\ A\ ,\ B\ \):traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ \ \ \ \ \ for\ example,\ cost\ is\ zero\ for\ the\ pair\ \(\ x\ ,\ y\ \):traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ \[!\]\ At\ pos\ \d+,\ \d+\ substitution\ cost\ value\(s\)\ is\ \d+!:traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.parametrize(
    "case, channels, kwargs",
    [
        ("unequal_diss_om", _unequal_channels, {"method": "OM", "sm": "CONSTANT"}),
        ("void_diss_ham", _void_channels, {"method": "HAM"}),
        ("void_diss_dhd", _void_channels, {"method": "DHD"}),
    ],
)
def test_void_component_distances_are_refused_with_r_message(case, channels, kwargs):
    with pytest.raises(ValueError) as excinfo:
        seqMD(channels(), what="diss", **kwargs)
    assert str(excinfo.value) == _r(case)["error"]


def test_a_per_state_indel_vector_is_undefined_for_a_void_component():
    """R's list indel is unusable in 2.2-14 ("subscript out of bounds"); the
    traminex per-state vector follows the same NA rule as the costs."""
    costs = seqMD(_void_channels(), sm="CONSTANT", indel=[[1.0, 2.0], [3.0, 4.0]], what="cost")
    alphabet = costs["alphabet"]
    indel = np.asarray(costs["indel"], dtype=float)
    assert np.isnan(indel[alphabet.index("%+y")])
    assert indel[alphabet.index("A+x")] == pytest.approx(4.0)
    assert indel[alphabet.index("B+y")] == pytest.approx(6.0)
