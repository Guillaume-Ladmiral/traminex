"""seqMD's combined (CAT) costs are symmetric even for an asymmetric channel matrix.

TraMineR computes the combined cost of states ``i < j`` once, from each
channel's ``sm[ipos, jpos]``, and writes it to both ``[i, j]`` and ``[j, i]``
(``seqMD.R:307-316``, and the same loop for time-varying costs). traminex
summed full channel blocks, so below the diagonal it used the transposed
channel costs: for an asymmetric channel matrix its CAT costs, and the
``what="diss"`` distances computed from them, differed from R's (2026-09-22
audit S10-13).

Expected values: live TraMineR 2.2-14 via
``scripts/r_goldens/make_goldens_audit20260922_s10b.R``
(``tests/r_goldens/seqMD_cat_costs.txt``).
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest
from r_goldens_utils import load_golden_lines

from traminex import seqdef, seqMD

ASY = np.array([[0.0, 3.0, 2.0], [1.0, 0.0, 2.0], [2.0, 2.0, 0.0]])
SYM = np.array([[0.0, 2.0, 2.0], [2.0, 0.0, 2.0], [2.0, 2.0, 0.0]])


def _recorded() -> dict[tuple[str, str], list[str]]:
    rows: dict[tuple[str, str], list[str]] = {}
    for line in load_golden_lines("seqMD_cat_costs"):
        if not line.startswith("#"):
            case, kind, *fields = line.split("\t")
            rows[(case, kind)] = fields
    return rows


R = _recorded()


def _r_array(case: str) -> np.ndarray:
    shape = tuple(int(value) for value in R[(case, "dim")])
    return np.array([float(value) for value in R[(case, "values")]]).reshape(shape, order="F")


def _channels():
    a = seqdef([list("ABAB"), list("BABA"), list("AABB")], alphabet=["A", "B", "U"])
    x = seqdef([list("XYXY"), list("YXYX"), list("XXYY")], alphabet=["X", "Y", "V"])
    return [a, x]


def _arrays():
    base = np.full((3, 3, 4), 2.0)
    for time in range(4):
        np.fill_diagonal(base[:, :, time], 0.0)
    asymmetric = base.copy()
    asymmetric[0, 1, 0] = 3.0  # R's A_asym[1, 2, 1] and [1, 2, 3]
    asymmetric[0, 1, 2] = 3.0
    return asymmetric, base


def _quiet(**kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return seqMD(_channels(), **kwargs)


@pytest.mark.parametrize("case", ["asym_stationary", "asym_time_varying"])
def test_cat_costs_of_an_asymmetric_channel_are_rs(case):
    if case == "asym_stationary":
        actual = _quiet(what="cost", sm=[ASY, SYM], indel=[1, 1])
    else:
        actual = _quiet(method="DHD", what="cost", sm=list(_arrays()), indel=[1, 1])
    np.testing.assert_allclose(actual["sm"], _r_array(case), rtol=0, atol=1e-12)
    assert actual["alphabet"] == R[(case, "alphabet")]


def test_distances_from_an_asymmetric_channel_are_rs():
    actual = _quiet(method="OM", what="diss", sm=[ASY, SYM], indel=[1, 1])
    np.testing.assert_allclose(actual, _r_array("asym_diss_om"), rtol=0, atol=1e-12)
