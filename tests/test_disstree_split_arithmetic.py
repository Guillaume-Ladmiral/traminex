"""``disstree`` evaluates candidate splits with TraMineR's arithmetic.

Audit 2026-09-22, finding S4-04. When a level subset and its complement have
mathematically equal residual inertia -- which happens for every candidate of an
unordered predictor, once as ``co`` and once as its complement -- the first one
enumerated is kept unless the other is smaller by rounding, so the last bit
decides which child is on the left, and so the node ids and the fitted leaves.
traminex summed with NumPy's pairwise ``sum`` and matrix products, TraMineR
(``SOURCES/TraMineR/R/disstree.R:403-523``) with R's ``sum()``, which
accumulates in long double in column-major order, and with the C kernels
``tmrWeightedInertiaDist`` / ``tmrWeightedInterInertia``
(``SOURCES/TraMineR/src/weightedinertia.cpp``), which accumulate doubles in
member order. On the audit's example R put ``{a, b}`` on the left and
traminex ``{c}``.

Expected values come from live R through
``scripts/r_goldens/make_goldens_audit20260922_s4.R`` and are compared exactly.
R's ``sum()`` uses the platform's long double; the fixture was produced on
Linux x86-64 (80-bit extended precision), so the exact comparison runs only
where NumPy's ``longdouble`` has that precision too.
"""

from __future__ import annotations

import importlib
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from traminex import disstree

# ``traminex.disstree`` is the re-exported function; the module holds the block size.
disstree_module = importlib.import_module("traminex.disstree")
diss_module = importlib.import_module("traminex.diss")
FIXTURE = Path(__file__).with_name("r_goldens") / "disstree_split_arithmetic.txt"

pytestmark = pytest.mark.skipif(
    np.finfo(np.longdouble).nmant != 63,
    reason="the R goldens were produced with an x87 80-bit long-double sum()",
)


def _cases() -> dict[str, dict[str, list[list[str]]]]:
    cases: dict[str, dict[str, list[list[str]]]] = {}
    for line in FIXTURE.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        cases.setdefault(fields[0], {}).setdefault(fields[1], []).append(fields[2:])
    return cases


CASES = _cases()


def _splits(node):
    if node.split is None:
        return []
    return [node] + _splits(node.kids[0]) + _splits(node.kids[1])


@pytest.mark.parametrize("block", [None, 3])
@pytest.mark.parametrize("case", sorted(CASES))
def test_splits_match_r_bit_for_bit(case, block, monkeypatch):
    if block is not None:
        # Carry the sequential accumulator across many small product blocks.
        monkeypatch.setattr(disstree_module, "_PRODUCT_BLOCK", block)
        monkeypatch.setattr(
            diss_module, "_PRODUCT_BLOCK", block
        )  # the within-group sum lives there
        # The blocks are those of the NumPy path; the compiled one has none.
        monkeypatch.setattr(diss_module, "_compiled_pair_sum", lambda: None)
    rows = CASES[case]
    D = np.array([float(token) for token in rows["D"][0][0].split(" ")])
    D = D.reshape(int(round(np.sqrt(D.size))), -1)
    f = np.array(rows["f"][0][0].split(" "), dtype=object)
    if rows["ordered"][0][0] == "TRUE":
        f = pd.Categorical(f, categories=sorted(set(f)), ordered=True)
    weights_text = rows["weights"][0][0]
    weights = None if weights_text == "none" else [float(v) for v in weights_text.split(" ")]

    tree = disstree(
        D,
        pd.DataFrame({"f": f}),
        weights=weights,
        R=1,
        min_size=0.1,
        max_depth=3,
        weight_permutation="diss",
    )

    observed = _splits(tree.root)
    assert len(observed) == len(rows["split"])
    for node, row in zip(observed, rows["split"]):
        fields = dict(item.split("=", 1) for item in row[1:])
        split = node.split
        assert node.id == int(row[0])
        assert "|".join(split.labels) == fields["labels"]
        assert ",".join(str(side) for side in split.index) == fields["index"]
        for key in ("lpop", "rpop", "lvar", "rvar", "SCres"):
            assert split.info[key] == float(fields[key]), key
    expected_fitted = [int(token) for token in rows["fitted"][0][0].split(" ")]
    assert tree.fitted.tolist() == expected_fitted
