"""check_costs: R's default is accepted; False is refused (TX-AUD-20261011-041).

traminex refused ``check_costs`` at every value in seqdist, seqcost, seqsubm,
seqMD and seqdistmc, while the ledger said only False was refused. R's
``check.costs = TRUE`` runs checkcost; with FALSE, R computes with a negative
indel, a non-zero diagonal or a wrongly sized matrix, and for negative costs it
returns negative "distances", which no dissimilarity can be; it does so even
with TRUE for an all-negative matrix, because checkcost returns after its
"max <= 0" message before the negative-cost check (TraMineR-checkcost.R:5-8,
52-54). traminex therefore keeps its cost checks mandatory: ``check_costs=True``
is accepted, ``check_costs=False`` refused. Fixture
``scripts/r_goldens/make_goldens_check_costs.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "check_costs.csv",
    dtype=str,
    keep_default_na=False,
)
S3 = tx.seqdef(
    np.array([["A", "B", "A"], ["B", "B", "A"], ["A", "A", "B"]], dtype=object), alphabet=["A", "B"]
)
RAW = np.array([[0, 1, 0], [1, 1, 0], [0, 0, 1]])
NEG = np.array([[0.0, -1.0], [-1.0, 0.0]])
CALLS = {
    "seqdist": lambda **kw: tx.seqdist(S3, method="OM", subm="CONSTANT", **kw),
    "seqcost": lambda **kw: tx.seqcost(S3, "CONSTANT", **kw)["sm"],
    "seqsubm": lambda **kw: tx.seqsubm(S3, method="CONSTANT", **kw),
    "seqsubm_raw": lambda **kw: tx.seqsubm(RAW, method="TRATE", **kw),
    "seqMD": lambda **kw: tx.seqMD(
        [S3, S3], method="OM", sm=["CONSTANT", "CONSTANT"], what="diss", **kw
    ),
    "seqdistmc": lambda **kw: tx.seqdistmc(
        [S3, S3], method="OM", sm=["CONSTANT", "CONSTANT"], **kw
    ),
}


def _r(case: str, check_costs: str) -> pd.Series:
    rows = GOLDEN[(GOLDEN["case"] == case) & (GOLDEN["check_costs"] == check_costs)]
    assert len(rows) == 1
    return rows.iloc[0]


def _values(row: pd.Series) -> np.ndarray:
    return np.array([float(v) for v in row["value"].split()])


@pytest.mark.parametrize("name", sorted(CALLS))
def test_check_costs_true_is_accepted(name: str) -> None:
    np.testing.assert_array_equal(
        np.asarray(CALLS[name](check_costs=True)), np.asarray(CALLS[name]())
    )


@pytest.mark.parametrize("name", sorted(CALLS))
def test_check_costs_false_is_refused(name: str) -> None:
    with pytest.raises(ValueError, match="check_costs=False is not supported"):
        CALLS[name](check_costs=False)


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ max\ substitution\ cost\ less\ or\ equal\ to\ \d+!:traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.parametrize("case", ["negative_indel", "negative_sm", "seqmd_negative_sm"])
def test_r_returns_negative_distances_where_traminex_refuses(case: str) -> None:
    assert _r(case, "FALSE")["ok"] == "TRUE" and _values(_r(case, "FALSE")).min() < 0
    with pytest.raises(ValueError):
        if case == "negative_indel":
            tx.seqdist(S3, method="OM", subm=np.array([[0.0, 2.0], [2.0, 0.0]]), indel=-1)
        elif case == "negative_sm":
            tx.seqdist(S3, method="OM", subm=NEG, indel=1)
        else:
            tx.seqMD([S3, S3], method="OM", sm=[NEG, NEG], indel=1, what="diss")


#: The generator's three cases R refuses with check.costs = TRUE, and the part of
#: traminex's refusal that names the fault.
TRUE_REFUSALS = {
    "nonzero_diagonal": (dict(subm=np.array([[1.0, 2.0], [2.0, 1.0]]), indel=1), "diagonal"),
    "negative_indel": (dict(subm=np.array([[0.0, 2.0], [2.0, 0.0]]), indel=-1), "indel"),
    "wrong_size": (
        dict(subm=np.array([[0, 2, 2], [2, 0, 2], [2, 2, 0]], dtype=float), indel=1),
        "2x2",
    ),
}


@pytest.mark.parametrize("case", sorted(TRUE_REFUSALS))
def test_r_refuses_these_with_check_costs_true(case: str) -> None:
    """R refuses each; so must traminex (only R's side was asserted: batch-6 review M-7)."""
    assert _r(case, "TRUE")["ok"] == "FALSE"
    kwargs, fault = TRUE_REFUSALS[case]
    with pytest.raises(ValueError, match=fault):
        tx.seqdist(S3, method="OM", check_costs=True, **kwargs)


def test_valid_costs_give_r_values() -> None:
    expected = _values(_r("triangle_violation", "TRUE")).reshape(3, 3, order="F")
    with pytest.warns(UserWarning):
        actual = tx.seqdist(S3, method="OM", subm=np.array([[0.0, 5.0], [5.0, 0.0]]), indel=1)
    np.testing.assert_allclose(actual, expected, atol=1e-12, rtol=0)
    trate = _values(_r("seqcost_trate", "TRUE")).reshape(2, 2, order="F")
    np.testing.assert_allclose(
        tx.seqcost(S3, "TRATE", check_costs=True)["sm"], trate, atol=1e-12, rtol=0
    )
