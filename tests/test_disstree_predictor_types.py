"""``disstree`` splits pandas categorical and boolean predictors as TraMineR does.

Audit 2026-09-22, finding S4-02. ``disstree`` took the levels of every
predictor from ``np.unique`` (sorted values), ignoring the declared order of a
pandas categorical, and emitted a numeric threshold (``breaks``) for every
ordered column, so an ordered categorical with string categories raised
``ValueError: could not convert string to float`` and one with integer
categories split over numeric order instead of the declared order. A boolean
column was treated as numeric (``b <= 0``).

TraMineR 2.2-14 (``SOURCES/TraMineR/R/disstree.R:404``) builds ``factor(pred,
ordered = is.ordered(pred) || is.numeric(pred))``, which keeps a factor's level
order and drops unused levels; it emits ``breaks`` only for a numeric predictor
(``disstree.R:543``) and splits every factor, ordered or not, over level labels.
A logical predictor is not numeric in R, so it is an unordered factor with
levels ``FALSE``/``TRUE``. Expected values come from live R through
``scripts/r_goldens/make_goldens_audit20260922_s4.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from traminex import disstree, disstree_assign, disstree_get_rules
from traminex._rformat import r_as_character

FIXTURE = Path(__file__).with_name("r_goldens") / "disstree_predictor_types.txt"


def _fixture() -> tuple[np.ndarray, dict[str, dict[str, list[list[str]]]]]:
    cases: dict[str, dict[str, list[list[str]]]] = {}
    pos = None
    for line in FIXTURE.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if fields[0] == "pos":
            pos = np.array([float(token) for token in fields[1].split(" ")])
            continue
        cases.setdefault(fields[0], {}).setdefault(fields[1], []).append(fields[2:])
    assert pos is not None
    return pos, cases


POS, CASES = _fixture()
D = np.abs(POS[:, None] - POS[None, :])
N = POS.size

PREDICTORS = {
    "ordered_string": pd.Categorical(
        ["low"] * 3 + ["high"] * 3 + ["mid"] * 6, categories=["low", "mid", "high"], ordered=True
    ),
    "ordered_integer": pd.Categorical(
        [3] * 3 + [2] * 3 + [1] * 6, categories=[3, 1, 2], ordered=True
    ),
    "logical": np.array([True] * 6 + [False] * 6),
    "factor_level_order": pd.Categorical(list("aaabbbcccddd"), categories=["d", "c", "b", "a"]),
    "ordered_missing_unused": pd.Categorical(
        ["x"] * 3 + [None] + ["y"] * 2 + ["z"] * 6, categories=["w", "z", "y", "x"], ordered=True
    ),
}


def _splits(node):
    if node.split is None:
        return []
    return [node] + _splits(node.kids[0]) + _splits(node.kids[1])


@pytest.mark.parametrize("case", sorted(PREDICTORS))
def test_tree_structure_matches_r(case):
    values = PREDICTORS[case]
    expected = CASES[case]
    assert [r_as_character(v) or "NA" for v in values] == expected["values"][0][0].split(" ")

    tree = disstree(D, pd.DataFrame({"v": values}), R=1, min_size=0, max_depth=3)

    observed = _splits(tree.root)
    assert len(observed) == len(expected["split"])
    for node, row in zip(observed, expected["split"]):
        fields = dict(item.split("=", 1) for item in row[1:])
        split = node.split
        assert node.id == int(row[0])
        assert split.varindex + 1 == int(fields["varindex"])
        assert split.breaks is None and fields["breaks"] == "NA"
        assert "|".join(r_as_character(label) for label in split.labels) == fields["labels"]
        assert ",".join(str(side) for side in split.index) == fields["index"]
        na_group = "NA" if split.na_group is None else str(split.na_group)
        assert na_group == fields["nagroup"]
        assert split.info["lpop"] == float(fields["lpop"])
        assert split.info["rpop"] == float(fields["rpop"])
        assert split.info["R2"] == pytest.approx(float(fields["R2"]), rel=1e-12)
    assert tree.fitted.tolist() == [int(token) for token in expected["fitted"][0][0].split(" ")]
    assert " || ".join(disstree_get_rules(tree)) == expected["rules"][0][0]


@pytest.mark.parametrize(
    "case", ["ordered_string", "ordered_integer", "logical", "factor_level_order"]
)
def test_assign_on_the_r_rules_matches_r(case):
    rules = CASES[case]["rules"][0][0].split(" || ")
    assigned = disstree_assign(rules, pd.DataFrame({"v": PREDICTORS[case]}))
    expected = [float(token) for token in CASES[case]["assign"][0][0].split(" ")]
    assert assigned.tolist() == expected
