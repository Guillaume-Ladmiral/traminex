"""``seqMD(what="cost")`` ignores ``method``, as TraMineR does.

In R, whether seqMD's costs vary over time depends only on the dimensions of
``sm`` (``seqMD.R:37``); the switch that makes DHD time-varying sits inside
``if (what == "diss")`` (``seqMD.R:81-93``). traminex made ``method="DHD"``
time-varying for ``what="cost"`` as well, so ``seqMD(what="cost",
method="DHD", sm="TRATE")`` returned a states-by-states-by-time array where R
returns the stationary matrix, and an explicit 2-D ``sm`` was refused
(2026-09-22 audit S10-12).

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


def _recorded() -> dict[tuple[str, str], list[str]]:
    rows: dict[tuple[str, str], list[str]] = {}
    for line in load_golden_lines("seqMD_cat_costs"):
        if line.startswith("#"):
            continue
        case, kind, *fields = line.split("\t")
        rows[(case, kind)] = fields
    return rows


R = _recorded()


def _channels():
    missing = seqdef([["A", None, "B", "A"], list("BABA"), list("ABAB")], alphabet=["A", "B", "U"])
    other = seqdef([list("XYXY"), list("YXYX"), list("XXYY")], alphabet=["X", "Y", "V"])
    return [missing, other]


@pytest.mark.parametrize("method", [None, "DHD", "OM", "HAM"])
def test_cost_is_rs_stationary_matrix_whatever_the_method(method):
    case = f"trate_method_{method or 'none'}"
    shape = tuple(int(value) for value in R[(case, "dim")])
    expected = np.array([float(value) for value in R[(case, "values")]]).reshape(shape, order="F")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        actual = seqMD(
            _channels(), method=method, what="cost", sm="TRATE", with_missing=[True, False]
        )

    assert np.asarray(actual["sm"]).shape == shape == (5, 5)
    np.testing.assert_allclose(actual["sm"], expected, rtol=0, atol=1e-12)
    np.testing.assert_allclose(
        np.atleast_1d(actual["indel"]), [float(v) for v in R[(case, "indel")]], rtol=0, atol=1e-12
    )
    assert actual["alphabet"] == R[(case, "alphabet")]


def test_an_explicit_matrix_with_method_dhd_is_not_forced_to_be_time_varying():
    """R takes a 2-D ``sm`` with ``what="cost"`` whatever ``method`` says."""
    sm = [np.full((4, 4), 2.0) - 2.0 * np.eye(4), np.full((3, 3), 2.0) - 2.0 * np.eye(3)]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with_dhd = seqMD(
            _channels(), method="DHD", what="cost", sm=sm, indel=1, with_missing=[True, False]
        )
        without = seqMD(_channels(), what="cost", sm=sm, indel=1, with_missing=[True, False])
    np.testing.assert_array_equal(with_dhd["sm"], without["sm"])
