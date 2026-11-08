"""Numeric ``disstree`` rules print their threshold as TraMineR does.

Audit 2026-09-22, finding S4-05. ``disstree_get_rules`` printed a numeric
threshold with C's ``%.5g`` (``1.2346e+05`` for 123456), where TraMineR prints
``format(x, digits = getOption("digits") - 2)``
(``SOURCES/TraMineR/R/disstree.R:63-70``), which keeps fixed notation unless it
is wider than scientific notation (``123456``). Because ``disstree_assign``
evaluates the printed rule, a case at 123457 was assigned to the ``<=`` leaf.
Expected rules and assignments come from live R through
``scripts/r_goldens/make_goldens_audit20260922_s4.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from traminex import disstree, disstree_assign, disstree_get_rules

FIXTURE = Path(__file__).with_name("r_goldens") / "disstree_threshold_rules.txt"


def _cases() -> dict[str, dict[str, str]]:
    cases: dict[str, dict[str, str]] = {}
    for line in FIXTURE.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        case, key, value = line.split("\t")
        cases.setdefault(case, {})[key] = value
    return cases


CASES = _cases()
POS = np.array([float(token) for token in CASES.pop("pos")["values"].split(" ")])
D = np.abs(POS[:, None] - POS[None, :])


@pytest.mark.parametrize("base", sorted(CASES, key=float))
def test_threshold_rule_and_assignment_match_r(base):
    expected = CASES[base]
    value = float(base)
    predictor = pd.DataFrame({"inc": [value] * 6 + [value + abs(value)] * 6})

    tree = disstree(D, predictor, R=1, min_size=0, max_depth=2)
    rules = disstree_get_rules(tree)

    assert tree.root.split.breaks == float(expected["breaks"])
    assert " || ".join(rules) == expected["rules"]
    profile = pd.DataFrame({"inc": [float(token) for token in expected["profile"].split(" ")]})
    expected_assign = [float(token) for token in expected["assign"].split(" ")]
    assert disstree_assign(rules, profile).tolist() == expected_assign
