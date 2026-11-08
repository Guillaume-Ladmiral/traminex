"""disstree_assign parses R-style covariate names (TX-AUD-20261010-006).

R evaluates every rule with ``eval(parse(text = rule))`` (disstree.R:143-176):
dotted names (``N.Eastern``) and backquoted names (```my var```) work. For a
non-syntactic covariate R's own ``disstree.get.rules`` writes the bare name
(``my var <= 5``), which R's ``disstree.assign`` then cannot read (both errors
are recorded); traminex assigns such rules, and the expected value is the rule
of each case's terminal node. Fixtures:
``scripts/r_goldens/make_goldens_disstree_assign_names.R``.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx
from traminex.disstree_helpers import disstree_assign

GOLDENS = Path(__file__).resolve().parent / "r_goldens"


def _load() -> dict[str, dict[str, list[str]]]:
    cases: dict[str, dict[str, list[tuple[int, str]]]] = defaultdict(lambda: defaultdict(list))
    with (GOLDENS / "disstree_assign_names.csv").open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            cases[row["case"]][row["kind"]].append((int(row["position"]), row["value"]))
    return {
        case: {kind: [value for _, value in sorted(items)] for kind, items in kinds.items()}
        for case, kinds in cases.items()
    }


R = _load()
R_ERRORS = dict(
    line.split("|", 1)
    for line in (GOLDENS / "disstree_assign_names_errors.txt")
    .read_text(encoding="utf-8")
    .splitlines()
    if line and not line.startswith("#")
)

# The profiles and covariates of the generator's disstree.assign calls.
MANUAL = {
    "manual_dotted": (pd.DataFrame({"N.Eastern": ["a", "b", "a"]}), ["N.Eastern"]),
    "manual_backquoted": (pd.DataFrame({"my var": ["x", "y", "x"]}), ["my var"]),
    "manual_comparison": (pd.DataFrame({"N.x": [1, 2, 3]}), ["N.x"]),
}
FITTED_PROFILE = pd.DataFrame(
    {
        "N.Eastern": ["a", "a", "b", "b", "a", "b", "a", "b"],
        "my var": [1, 1, 5, 5, 2, 6, 2, 6],
    }
)


@pytest.mark.parametrize("case", sorted(MANUAL))
def test_r_style_names_assign_as_in_r(case: str):
    profile, covar = MANUAL[case]
    actual = disstree_assign(R[case]["rule"], profile, covar=covar)
    np.testing.assert_array_equal(actual, [float(value) for value in R[case]["assign"]])


def test_rules_r_writes_for_a_non_syntactic_name_are_assigned_by_leaf():
    # R cannot read its own rules here (recorded errors); the definition-derived
    # value is the rule index of each case's terminal node.
    assert R_ERRORS["covar_attribute"].startswith("covariates missing in profile")
    assert "unexpected symbol" in R_ERRORS["covar_bare"]
    case = R["fitted_nonsyntactic"]
    actual = disstree_assign(case["rule"], FITTED_PROFILE, covar=["N.Eastern", "my var"])
    np.testing.assert_array_equal(actual, [float(value) for value in case["leaf"]])


def test_traminex_rules_for_a_non_syntactic_name_round_trip():
    rows = ["AABBC", "ABBCC", "CCBAA", "BBCCA", "AAABB", "CCCBB", "AABBB", "CCCBA"]
    seqs = tx.seqdef(np.array([list(row) for row in rows], dtype=object), alphabet=["A", "B", "C"])
    diss = tx.seqdist(seqs, method="LCS")
    tree = tx.disstree(diss, FITTED_PROFILE, R=1, pval=1.0, min_size=0.1, random_state=1)
    assigned = disstree_assign(tx.disstree_get_rules(tree), FITTED_PROFILE)
    leaves = tx.disstreeleaf(tree)
    # Cases share a rule exactly when they share a terminal node.
    pairs = {(int(leaf), int(rule)) for leaf, rule in zip(leaves, assigned)}
    assert len(pairs) == len({leaf for leaf, _ in pairs}) == len({rule for _, rule in pairs})
    assert len(pairs) > 1
