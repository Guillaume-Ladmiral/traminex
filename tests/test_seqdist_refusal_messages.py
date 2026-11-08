"""seqdist's refusals for HAM's norm, DHD's CONSTANT costs and OMstran's generated costs say what is wrong
(TX-AUD-20261019-079).

For an invalid norm, seqdist(method="HAM") said "'gmean'/'maxdist' are not implemented", yet both compute and
match R (R: "'norm' must be one of: auto, none, maxlength, gmean, maxdist, YujianBo"); since seqdist checks the
norm spelling first, as R does (TX-REP-20261019-010), the message is the norm vocabulary's. For DHD with
sm = "CONSTANT", R says "'sm = \"CONSTANT\"' is not relevant for DHD, consider HAM instead" (seqdist.R:253;
``scripts/r_goldens/make_goldens_seqdist_check_order.R`` row ``dhd_constant``); traminex now says the same.
"""

from __future__ import annotations

import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "seqdist_check_order.csv",
    dtype=str,
    keep_default_na=False,
).set_index("case")
X = tx.seqdef(
    np.array([list("AABC"), list("ABBC"), list("CCBA"), list("BBCA")], dtype=object),
    alphabet=["A", "B", "C"],
)


def test_ham_invalid_norm_message_does_not_deny_gmean_and_maxdist() -> None:
    with pytest.raises(ValueError) as excinfo:
        tx.seqdist(X, method="HAM", norm="bogus")
    assert "not implemented" not in str(excinfo.value)
    assert "gmean" in str(excinfo.value) and "maxdist" in str(excinfo.value)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for norm in ("gmean", "maxdist"):
            assert np.all(np.isfinite(tx.seqdist(X, method="HAM", norm=norm)))


def test_dhd_constant_refusal_is_r_message() -> None:
    assert GOLDEN.loc["dhd_constant", "ok"] == "FALSE"
    with pytest.raises(ValueError, match=re.escape(GOLDEN.loc["dhd_constant", "error"])):
        tx.seqdist(X, method="DHD", subm="CONSTANT")


def test_omstran_names_the_non_finite_generated_substitution_costs() -> None:
    """With two states the data never show, seqcost's INDELS indels are the largest double, and the
    substitution cost between those two states overflows to infinity (R: seqcost.R:152-153; R itself stops
    in checkcost). OMstran's refusal blamed the indel (batch-9 review I-1)."""
    x = tx.seqdef(
        np.array([list("AABA"), list("ABBA"), list("BBAA")], dtype=object),
        alphabet=["A", "B", "C", "D"],
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert not np.all(np.isfinite(np.asarray(tx.seqcost(x, "INDELS")["sm"], dtype=float)))
        with pytest.raises(ValueError) as excinfo:
            tx.seqdist(x, method="OMstran", subm="INDELS", otto=0.5)
    message = str(excinfo.value)
    assert "indel cost should be positive" not in message
    assert re.search(r"substitution costs\b.*not finite", message), message
