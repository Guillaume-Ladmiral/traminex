"""check_costs is read as R's if() reads check.costs (TX-REP-20261019-009, batch-6 review M-8).

R evaluates ``if (check.costs)``: TRUE and FALSE, a number (0 is FALSE), and the strings as.logical accepts
("TRUE", "T", "true", "FALSE", "F", "false"); any other string stops with "argument is not interpretable as
logical", NA with "missing value where TRUE/FALSE needed" and NULL with "argument is of length zero"
(``scripts/r_goldens/make_goldens_seqdist_check_order.R``: with a non-zero diagonal R refuses when it reads TRUE
and computes when it reads FALSE). traminex read Python truthiness, so "FALSE" and "no" ran the checks and None was
reported as "check_costs=False". traminex keeps its refusal of R's FALSE (TX-AUD-20261011-041).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx
from traminex.distances._validation import CHECK_COSTS_FALSE_MESSAGE

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "seqdist_check_order.csv",
    dtype=str,
    keep_default_na=False,
).set_index("case")
X = tx.seqdef(
    np.array([list("AABC"), list("ABBC"), list("CCBA"), list("BBCA")], dtype=object),
    alphabet=["A", "B", "C"],
)
#: The golden's check.costs values in Python: R's NA is a float NaN, R's NULL is None.
VALUES = {
    "string_TRUE": "TRUE",
    "string_T": "T",
    "string_true": "true",
    "string_FALSE": "FALSE",
    "string_F": "F",
    "string_false": "false",
    "string_no": "no",
    "number_0": 0,
    "number_1": 1,
    "number_2_5": 2.5,
    "NA": float("nan"),
    "NULL": None,
    "TRUE": True,
    "FALSE": False,
}


def _r_reading(name: str) -> str:
    """ "TRUE", "FALSE" or R's error for one value, from the golden."""
    if f"check_costs_bad_{name}" in GOLDEN.index:
        return "FALSE" if GOLDEN.loc[f"check_costs_bad_{name}", "ok"] == "TRUE" else "TRUE"
    row = GOLDEN.loc[f"check_costs_{name}"]
    return "TRUE" if row["ok"] == "TRUE" else row["error"]


CALLS = {
    "seqdist": lambda value: tx.seqdist(X, method="OM", subm="CONSTANT", check_costs=value),
    "seqcost": lambda value: tx.seqcost(X, "CONSTANT", check_costs=value),
    "seqsubm": lambda value: tx.seqsubm(X, method="CONSTANT", check_costs=value),
    "seqMD": lambda value: tx.seqMD(
        [X, X], method="OM", sm=["CONSTANT", "CONSTANT"], what="diss", check_costs=value
    ),
}


@pytest.mark.parametrize("entry", sorted(CALLS))
@pytest.mark.parametrize("name", sorted(VALUES))
def test_check_costs_is_read_as_r_reads_it(name: str, entry: str) -> None:
    reading = _r_reading(name)
    call = CALLS[entry]
    if reading == "TRUE":
        np.testing.assert_array_equal(
            np.asarray(call(VALUES[name]) if entry != "seqcost" else call(VALUES[name])["sm"]),
            np.asarray(call(True) if entry != "seqcost" else call(True)["sm"]),
        )
    elif reading == "FALSE":
        with pytest.raises(ValueError, match="check_costs=False is not supported"):
            call(VALUES[name])
    else:
        with pytest.raises(ValueError, match=reading):
            call(VALUES[name])


def test_r_reads_the_values_as_expected() -> None:
    readings = {name: _r_reading(name) for name in VALUES}
    assert (
        readings["string_FALSE"]
        == readings["string_F"]
        == readings["string_false"]
        == readings["number_0"]
        == "FALSE"
    )
    assert readings["string_T"] == readings["string_true"] == readings["number_2_5"] == "TRUE"
    assert readings["string_no"] == "argument is not interpretable as logical"
    assert readings["NA"] == "missing value where TRUE/FALSE needed"
    assert readings["NULL"] == "argument is of length zero"


def test_numpy_booleans_are_read_as_booleans() -> None:
    with pytest.raises(ValueError, match="check_costs=False is not supported"):
        tx.seqdist(X, method="OM", subm="CONSTANT", check_costs=np.bool_(False))
    np.testing.assert_array_equal(
        tx.seqdist(X, method="OM", subm="CONSTANT", check_costs=np.bool_(True)),
        tx.seqdist(X, method="OM", subm="CONSTANT"),
    )
    assert "check_costs=False" in CHECK_COSTS_FALSE_MESSAGE


def test_zero_dimensional_boolean_arrays_are_read_as_booleans() -> None:
    """A 0-d array is a scalar: np.array(True) is R's TRUE (batch-9 review I-2)."""
    np.testing.assert_array_equal(
        tx.seqdist(X, method="OM", subm="CONSTANT", check_costs=np.array(True)),
        tx.seqdist(X, method="OM", subm="CONSTANT"),
    )
    with pytest.raises(ValueError, match="check_costs=False is not supported"):
        tx.seqdist(X, method="OM", subm="CONSTANT", check_costs=np.array(False))


GOOD = np.array([[0, 2, 2], [2, 0, 2], [2, 2, 0]], dtype=float)
#: Routes on which R ignores check.costs (seqdist.R:370-460: it reads the value only for a matrix sm of the OM
#: family and TWED, an array for DHD and costs seqcost generates; seqMD.R:201 returns what = "MDseq" before
#: seqMD.R:266 reads it): the golden rows, their call in traminex and the value passed.
UNREAD = {
    "check_costs_unread_lcs_no": (lambda v: tx.seqdist(X, method="LCS", check_costs=v), "no"),
    "check_costs_unread_lcs_NA": (
        lambda v: tx.seqdist(X, method="LCS", check_costs=v),
        float("nan"),
    ),
    "check_costs_unread_lcs_NULL": (lambda v: tx.seqdist(X, method="LCS", check_costs=v), None),
    "check_costs_unread_ham_matrix_no": (
        lambda v: tx.seqdist(X, method="HAM", subm=GOOD, check_costs=v),
        "no",
    ),
    "check_costs_unread_ham_default_no": (
        lambda v: tx.seqdist(X, method="HAM", check_costs=v),
        "no",
    ),
    "check_costs_unread_dhd_default_no": (
        lambda v: tx.seqdist(X, method="DHD", check_costs=v),
        "no",
    ),
    "check_costs_unread_chi2_no": (lambda v: tx.seqdist(X, method="CHI2", check_costs=v), "no"),
    "check_costs_unread_seqmd_mdseq_no": (
        lambda v: tx.seqMD([X, X], what="MDseq", check_costs=v),
        "no",
    ),
}
#: R's if() reading of each value, from the golden's OM rows, where R reads it.
READING = {"no": "string_no", "NA": "NA", "NULL": "NULL"}


@pytest.mark.parametrize("case", sorted(UNREAD))
def test_traminex_reads_check_costs_where_r_ignores_it(case: str) -> None:
    """A traminex rule: R computes, traminex reads the value as R would and refuses (batch-9 review I-2).

    traminex checks the costs on every route and refuses FALSE everywhere (TX-AUD-20261011-041), so it reads the
    value everywhere too, with the message R gives where it does read it.
    """
    assert GOLDEN.loc[case, "ok"] == "TRUE", GOLDEN.loc[case, "error"]
    call, value = UNREAD[case]
    reading = _r_reading(READING[case.rsplit("_", 1)[1]])
    with pytest.raises(ValueError, match=reading):
        call(value)
    np.asarray(call(True))
