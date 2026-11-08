"""seqdist checks its arguments in R's order (TX-REP-20261019-010, batch-6 review M-9).

R checks norm, indel, sm, prox, link and step in that order for every method (seqdist.R:127-190), so with two
invalid arguments R names the earlier one. traminex checked link, step, prox, sm and indel, then norm, and named
link where R names norm, and prox where R names sm. The rows of
``scripts/r_goldens/make_goldens_seqdist_check_order.R`` pair every two adjacent checks and the four pairs of the
batch-6 review; each refusal must name the argument R names, the first one R's message quotes.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx
from traminex.distances import pairwise

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "seqdist_check_order.csv",
    dtype=str,
    keep_default_na=False,
).set_index("case")
X = tx.seqdef(
    np.array([list("AABC"), list("ABBC"), list("CCBA"), list("BBCA")], dtype=object),
    alphabet=["A", "B", "C"],
)
#: The generator's calls in traminex's spelling (subm for R's sm).
CALLS = {
    "order_norm_indel": dict(subm="CONSTANT", norm="bogus", indel="bogus"),
    "order_indel_sm": dict(indel="bogus", subm="bogus"),
    "order_sm_prox": dict(subm="bogus", prox="bogus"),
    "order_prox_link": dict(subm="CONSTANT", prox="bogus", link="bogus"),
    "order_link_step": dict(subm="CONSTANT", link="bogus", step=0),
    "order_norm_link": dict(subm="CONSTANT", norm="bogus", link="bogus"),
    "order_indel_step": dict(subm="CONSTANT", indel="bogus", step=0),
    "order_norm_step": dict(subm="CONSTANT", norm="bogus", step=0),
    # A logical or complex sm is refused at R's sm step, after norm and indel (batch-9 review M-1, M-5).
    "order_norm_logical_sm": dict(norm="bogus", subm=True),
    "order_indel_logical_sm": dict(indel="bogus", subm=True),
    "order_norm_complex_sm": dict(norm="bogus", subm=np.complex128(2)),
    # R does not trim a norm spelling (batch-9 review M-4).
    "order_norm_spaces_indel": dict(subm="CONSTANT", norm=" maxlength ", indel="bogus"),
}


def test_the_golden_rows_are_r_refusals() -> None:
    assert all(GOLDEN.loc[case, "ok"] == "FALSE" for case in CALLS)


@pytest.mark.parametrize("case", sorted(CALLS))
def test_the_refusal_names_the_argument_r_names(case: str) -> None:
    named = re.search(r"'([A-Za-z_.]+)'", GOLDEN.loc[case, "error"]).group(1)
    with pytest.raises(ValueError) as excinfo:
        tx.seqdist(X, method="OM", **CALLS[case])
    message = str(excinfo.value)
    others = {"norm", "indel", "sm", "prox", "link", "step"} - {named}
    assert re.search(rf"\b{named}\b", message), (named, message)
    # The message is about R's argument alone: no later argument is named instead.
    assert not any(re.search(rf"^'?{other}\b", message) for other in others), message


@pytest.mark.parametrize("subm", [np.complex128(2), complex(2, 0), np.array(2 + 0j)])
def test_a_complex_sm_is_refused_as_r_refuses_it(subm) -> None:
    """R's sm is a matrix, an array or a method name (seqdist.R:153-168); traminex adds real numbers, and a
    complex number is not one: it computed with its real part (batch-9 review M-5)."""
    assert GOLDEN.loc["sm_complex", "ok"] == "FALSE"
    with pytest.raises(ValueError, match=re.escape(GOLDEN.loc["sm_complex", "error"])):
        tx.seqdist(X, method="OM", subm=subm)
    with pytest.raises(ValueError, match=re.escape(GOLDEN.loc["sm_complex", "error"])):
        pairwise(np.array([[0, 1, 2], [1, 1, 2]]), method="OM", subm=subm, indel=1.0)
