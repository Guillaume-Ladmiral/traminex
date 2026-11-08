"""seqerulesdisc leaves out a rule whose predictor is aliased (TX-AUD-20261011-010).

When a rule's predictor event is constant in its model data, R's glm aliases the
coefficient and seqerulesdisc reports the next row, AGE's, under the
predictor's name (seqerulesdisc.R:81-87: ``A(HR=1.4,p=0) => B``). traminex raised
"rank deficient" and returned no rule at all. The coefficient has no estimate,
so its p-value is missing and, as R's own filter does for a missing p-value
(:96-98), that rule is left out; the others are returned with R's values.
Fixture: ``scripts/r_goldens/make_goldens_erulesdisc_aliased.R`` (R's rules and
which predictors R's glm aliases).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "erulesdisc_aliased.csv",
    dtype=str,
    keep_default_na=False,
)


def _rules(constant_a: bool) -> list[str]:
    ids = np.arange(1, 49)
    data = pd.DataFrame(
        {
            "IDPERS": ids,
            "A": np.full(48, 18) if constant_a else 19 + ids % 4,
            "AST": np.ones(48, int) if constant_a else (ids % 4 != 0).astype(int),
            "B": 21 + ids % 4,
            "BST": (ids % 5 != 0).astype(int),
            "C": 23 + ids % 3,
            "CST": (ids % 3 != 0).astype(int),
        }
    )
    ages = pd.DataFrame({"IDPERS": ids, "AGE": 18}), pd.DataFrame({"IDPERS": ids, "AGE": 26})
    discrete = tx.createdatadiscrete(list(ids), data, ["A", "B", "C"], *ages)
    events = [
        (row.IDPERS, getattr(row, event), event)
        for row in data.itertuples()
        for event in "ABC"
        if getattr(row, event + "ST") == 1
    ]
    eseq = tx.seqecreate(
        id=[e[0] for e in events], timestamp=[e[1] for e in events], event=[e[2] for e in events]
    )
    subseq = tx.seqefsub(eseq, str_subseq=["(A)-(B)", "(B)-(C)", "(A)-(C)"])
    result = tx.seqerulesdisc(subseq, discrete, data, pvalue=1, adjust=False)
    return result["Hazards"].tolist() if len(result) else []


def _r(case: str, kind: str) -> pd.DataFrame:
    return GOLDEN[(GOLDEN["case"] == case) & (GOLDEN["kind"] == kind)]


def test_control_rules_are_r_rules() -> None:
    assert _rules(constant_a=False) == _r("control", "r_rule")["key"].tolist()


def test_a_rule_with_an_aliased_predictor_is_left_out() -> None:
    aliased = _r("constant_a", "aliased")
    kept = set(aliased.loc[aliased["value"] == "FALSE", "key"])
    expected = [
        rule
        for rule in _r("constant_a", "r_rule")["key"]
        if f"{rule.split('(')[0]}=>{rule.rsplit(' ', 1)[-1]}" in kept
    ]
    assert expected == ["B(HR=1.73,p=0.2591) => C"]
    with pytest.warns(UserWarning, match="aliased"):
        assert _rules(constant_a=True) == expected
