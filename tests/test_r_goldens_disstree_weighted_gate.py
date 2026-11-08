"""``disstree`` and ``seqtree`` against TraMineR 2.2-14, on what no test here ran before.

The supervision repository's measurement of 2026-11-02 (its
``docs/TRAMINEX_COMPONENT_GAPS.md``) found that no suite of the family ran the weighted
branches of the split permutation test (``traminex/disstree.py``,
``_split_permutation_pval``; R: ``SOURCES/TraMineR/R/disstree.dissassoc.R:52-64``), though
``disstree(D, predictors, weights=w)`` enters one with its default arguments. Only
pyweightedcluster's tests ran ``squared``, ``first`` and the rounding of weights.

R and NumPy draw different random permutations, so the generator
(``scripts/r_goldens/make_goldens_disstree_weighted_gate.R``) replaces ``sample`` with a
deterministic script and writes each answer into the fixture, in the order R asked. These
tests replay those rows through a NumPy ``Generator``. The whole tree must then equal R's,
every p-value exactly, since a p-value is a count over ``R``; and every row must be used,
so the order and the size of R's requests are pinned too.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDEN_ATOL, R_GOLDENS_DIR

from traminex import disstree, seqdef, seqtree

GATE = {"R": 30, "pval": 1.0}
DISSTREE_CASES = {
    "gate_diss": ("real", {**GATE, "weight_permutation": "diss"}),
    "gate_group": ("real", {**GATE, "weight_permutation": "group"}),
    "gate_replicate": ("integer", {**GATE, "weight_permutation": "replicate"}),
    "gate_rounded": ("real", {**GATE, "weight_permutation": "rounded-replicate"}),
    "gate_none_weighted": ("real", {**GATE, "weight_permutation": "none"}),
    "gate_rejects": ("integer", {"R": 30, "pval": 0.05, "weight_permutation": "replicate"}),
    "weak_replicate": (
        "integer",
        {"R": 100, "pval": 1.0, "weight_permutation": "replicate", "max_depth": 2},
    ),
    "weak_rounded": (
        "real",
        {"R": 100, "pval": 1.0, "weight_permutation": "rounded-replicate", "max_depth": 2},
    ),
    "squared_unweighted": ("none", {"R": 1, "squared": True}),
    "squared_weighted": ("real", {"R": 1, "squared": True, "weight_permutation": "diss"}),
    "first_phase": ("none", {"R": 1, "first": "phase"}),
    "first_cohort": ("real", {"R": 1, "first": "cohort", "weight_permutation": "diss"}),
    "first_no_split": ("none", {"R": 1, "first": "cohort", "min_size": 8}),
    "first_gate": ("real", {**GATE, "first": "cohort", "weight_permutation": "group"}),
    "rounded_deterministic": ("real", {"R": 1, "weight_permutation": "rounded-replicate"}),
    "rounded_to_zero": ("rounds_to_zero", {"R": 1, "weight_permutation": "rounded-replicate"}),
}
COLUMNS = {
    "rounded_to_zero": ["phase"],
    "weak_replicate": ["noise"],
    "weak_rounded": ["noise"],
}
SEQTREE_CASES = {
    "seqtree_first_diss": ("real", {**GATE, "first": "phase", "weight_permutation": "diss"}),
    "seqtree_squared": ("integer", {"R": 1, "squared": True, "weight_permutation": "replicate"}),
}


def _rows() -> dict[str, dict[str, list[str]]]:
    rows: dict[str, dict[str, list[str]]] = {}
    text = (R_GOLDENS_DIR / "disstree_weighted_gate.txt").read_text(encoding="utf-8")
    for line in text.splitlines():
        if not line.startswith("#"):
            case, key, value = line.split("\t")
            rows.setdefault(case, {}).setdefault(key, []).append(value)
    return rows


ROWS = _rows()
INPUT = {key: value[0] for key, value in ROWS["input"].items()}


def _numbers(text: str) -> np.ndarray:
    return np.array([float("nan") if token == "NA" else float(token) for token in text.split(" ")])


D = _numbers(INPUT["D"]).reshape(14, 14)
WEIGHTS = {name: _numbers(INPUT[name]) for name in ("real", "integer", "rounds_to_zero")}
ALL_PREDICTORS = pd.DataFrame(
    {
        "phase": INPUT["phase"].split(" "),
        "cohort": _numbers(INPUT["cohort"]),
        "noise": INPUT["noise"].split(" "),
    }
)
PREDICTORS = ALL_PREDICTORS[["phase", "cohort"]]


class ScriptedGenerator(np.random.Generator):
    """A generator whose ``permutation`` replays R's 1-based rows, refusing another size."""

    def __init__(self, rows: list[str]):
        super().__init__(np.random.PCG64(0))
        self.script = [np.array(row.split(" "), dtype=int) - 1 for row in rows]
        self.calls = 0

    def permutation(self, x):  # type: ignore[override]
        row = self.script[self.calls]
        self.calls += 1
        assert int(x) == row.size
        return row.copy()


def _flatten(tree) -> np.ndarray:
    """The generator's ``flatten_tree``: preorder, 1-based medoid and varindex."""
    rows = []

    def visit(node) -> None:
        split = node.split
        tail = [np.nan] * 9
        if split is not None:
            tail = [
                split.varindex + 1,
                split.info["SCres"],
                split.info["lpop"],
                split.info["rpop"],
                split.info["lvar"],
                split.info["rvar"],
                split.info["R2"],
                split.info.get("pval", np.nan),
                np.nan if split.breaks is None else float(split.breaks),
            ]
        info = node.info
        rows.append(
            [node.id, info["depth"], info["splitschedule"], info["n"], info["vardis"]]
            + [info["medoid"] + 1, float(split is not None)]
            + tail
        )
        if node.kids is not None:
            visit(node.kids[0])
            visit(node.kids[1])

    visit(tree.root)
    return np.array(rows, dtype=float)


def _assert_tree_equals_r(tree, case: str, generator: ScriptedGenerator) -> None:
    expected = ROWS[case]
    assert generator.calls == len(generator.script), "R asked for a different number"
    nodes = np.array([_numbers(row) for row in expected["node"]])
    observed = _flatten(tree)
    assert observed.shape == nodes.shape
    # The p-values, column 14, are counts over R: exactly R's, not merely close.
    np.testing.assert_array_equal(observed[:, 14], nodes[:, 14])
    np.testing.assert_allclose(observed, nodes, rtol=0, atol=R_GOLDEN_ATOL)
    np.testing.assert_array_equal(tree.fitted, _numbers(expected["fitted"][0]).astype(int))
    stat = tree.adjustment.stat
    np.testing.assert_allclose(
        stat["t0"].to_numpy(), _numbers(expected["adjust_t0"][0]), rtol=0, atol=R_GOLDEN_ATOL
    )
    np.testing.assert_array_equal(stat["p_value"].to_numpy(), _numbers(expected["adjust_p"][0]))


def _weights(name: str):
    return None if name == "none" else WEIGHTS[name]


@pytest.mark.parametrize("case", sorted(DISSTREE_CASES))
def test_disstree_matches_r(case):
    name, options = DISSTREE_CASES[case]
    assert ROWS[case]["weights"] == [name]
    columns = COLUMNS.get(case, ["phase", "cohort"])
    generator = ScriptedGenerator(ROWS[case].get("perm", []))
    tree = disstree(
        D,
        ALL_PREDICTORS[columns],
        weights=_weights(name),
        random_state=generator,
        **{"max_depth": 3, **options},
    )
    _assert_tree_equals_r(tree, case, generator)


@pytest.mark.parametrize("case", sorted(SEQTREE_CASES))
def test_seqtree_matches_r(case):
    name, options = SEQTREE_CASES[case]
    assert ROWS[case]["weights"] == [name]
    states = [list(row) for row in INPUT["states"].split(" ")]
    sequences = seqdef(states, alphabet=["A", "B", "C"], weights=_weights(name))
    generator = ScriptedGenerator(ROWS[case].get("perm", []))
    tree = seqtree(sequences, PREDICTORS, diss=D, max_depth=3, random_state=generator, **options)
    _assert_tree_equals_r(tree, case, generator)


def test_every_recorded_tree_has_a_test():
    recorded = {case for case, rows in ROWS.items() if "node" in rows}
    assert recorded == set(DISSTREE_CASES) | set(SEQTREE_CASES)


@pytest.mark.parametrize("case", ["weak_replicate", "weak_rounded"])
def test_the_weak_cases_stay_where_a_count_is_sensitive(case):
    # Otherwise these would lose what they are for: a wrong statistic moves a mid-range count.
    root = _numbers(ROWS[case]["node"][0])
    assert 0.1 < root[14] < 0.9


def test_the_threshold_case_rejects_a_split_the_open_case_keeps():
    # Otherwise gate_rejects would not pin that the p-value decides the tree's shape.
    assert len(ROWS["gate_rejects"]["node"]) < len(ROWS["gate_replicate"]["node"])


def test_real_weights_are_refused_for_replicates_as_in_r():
    assert ROWS["refuse_replicate_real"]["refusal"] == [
        " [!] To permute replicate, you should specify integer weights"
    ]
    with pytest.raises(ValueError, match="to permute replicates, weights must be integers"):
        disstree(D, PREDICTORS, weights=WEIGHTS["real"], R=1, weight_permutation="replicate")


def test_an_unknown_first_variable_is_refused():
    # R stops too, though not on purpose: which(...) is empty, and its check of the index
    # fails on the missing value ("missing value where TRUE/FALSE needed").
    assert ROWS["refuse_first_absent"]["refusal"] == ["missing value where TRUE/FALSE needed"]
    with pytest.raises(ValueError, match="unknown 'first' variable 'absent'"):
        disstree(D, PREDICTORS, R=1, first="absent")


def test_weights_rounded_to_zero_stop_where_r_stops():
    # A candidate group left with no weight gives a NaN residual, and R stops on comparing it.
    assert ROWS["refuse_rounded_to_zero"]["refusal"] == ["missing value where TRUE/FALSE needed"]
    with pytest.raises(ValueError, match="cannot compare candidate splits"):
        disstree(
            D,
            PREDICTORS,
            weights=WEIGHTS["rounds_to_zero"],
            R=1,
            weight_permutation="rounded-replicate",
        )


def test_the_defaults_enter_the_replicate_branch():
    # disstree(D, predictors, weights=w) is the call that made this gap matter: its defaults
    # are R's (SOURCES/TraMineR/R/disstree.R:238-239), R=1000 and "replicate". R draws other
    # permutations, so this checks only that the call runs and gives counts over 1000.
    tree = disstree(D, PREDICTORS, weights=WEIGHTS["integer"], max_depth=2, random_state=7)
    assert tree.info["parameters"]["R"] == 1000
    assert tree.info["weight_permutation"] == "replicate"
    if tree.root.split is not None:
        count = tree.root.split.info["pval"] * 1000
        assert count == pytest.approx(round(count), abs=1e-9)
    assert tree.adjustment.stat.loc["Pseudo F", "p_value"] * 1000 == pytest.approx(
        round(tree.adjustment.stat.loc["Pseudo F", "p_value"] * 1000), abs=1e-9
    )
