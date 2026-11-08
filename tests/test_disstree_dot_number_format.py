"""Tree DOT labels format numbers with R's ``format(x, digits = d)`` rule.

What was wrong (2026-09-22 audit S6-03): ``disstree2dot`` and its siblings
formatted node sizes, discrepancies and R2 with Python's ``"%.{d}g"``, which
switches to scientific notation as soon as a value has more integer digits than
``d``: a node of 1200 cases was labelled ``n: 1.2e+03`` where TraMineR writes
``n: 1200``. The global-quality lines also used the caller's ``digits`` (3)
where R uses ``getOption("digits") - 2`` (5).

Source of truth: TraMineR 2.2-14 ``disstree2dot.R`` -- ``formd`` is
``format(x, digits = digits)`` for node labels and
``format(x, digits = getOption("digits") - 2)`` for the quality lines. R's
``format`` prefers fixed notation unless it is wider than scientific notation.
Expected strings are live R output in ``tests/r_goldens/r_format_digits.csv``
(``scripts/r_goldens/make_goldens_audit20260922_s6.R``).
"""

from __future__ import annotations

import csv
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import load_golden

from traminex.disstree import disstree, disstree2dot, disstree2dotp
from traminex.disstree_viz import _fmt, _quality_lines

FIXTURE = Path(__file__).with_name("r_goldens") / "r_format_digits.csv"


def _r_format_rows():
    with FIXTURE.open(newline="") as handle:
        header = handle.readline()
        rows = list(csv.DictReader(handle))
    assert header == '# getOption("digits") = 7\n'
    return [(float(row["value"]), int(row["digits"]), row["formatted"]) for row in rows]


R_FORMAT = _r_format_rows()
R_FORMAT_BY_KEY = {(value, digits): text for value, digits, text in R_FORMAT}


@pytest.mark.parametrize(
    ("value", "digits", "expected"),
    R_FORMAT,
    ids=[f"{value!r}-d{digits}" for value, digits, _ in R_FORMAT],
)
def test_fmt_matches_r_format(value, digits, expected):
    assert _fmt(value, digits) == expected


def _tree():
    D = load_golden("phase2_input_D")
    i = np.arange(1, 16)
    predictors = pd.DataFrame(
        {
            "A": ["a", "a", "b", "b", "a", "b", "b", "b", "c", "c", "c", "c", "c", "a", "c"],
            "X": ((i * 7) % 15) + 0.1 * i,
        }
    )
    return disstree(D, predictors, R=1, min_size=2, max_depth=4)


def test_node_label_of_a_large_node_keeps_fixed_notation(tmp_path, monkeypatch):
    # The audit's reproduction: a node of 1200 cases.
    monkeypatch.chdir(tmp_path)
    tree = _tree()
    tree.root.info["n"] = 1200
    disstree2dot(tree, "tree", use_title=False)

    root = next(
        line
        for line in (tmp_path / "tree.dot").read_text().splitlines()
        if line.startswith('"tree_none"[')
    )
    assert f'label="n: {R_FORMAT_BY_KEY[(1200.0, 3)]} s2: ' in root
    assert "e+03" not in root


def test_quality_lines_use_r_five_significant_digits():
    stat = pd.DataFrame(
        {"t0": [12.34, 0.567, 1234.5], "p_value": [0.0005, 0.5, 0.04]},
        index=["Pseudo F", "Pseudo R2", "Levene"],
    )
    tree = SimpleNamespace(adjustment=SimpleNamespace(stat=stat))

    assert _quality_lines(tree) == [
        f"Pseudo F: {R_FORMAT_BY_KEY[(12.34, 5)]} ***",
        f"Pseudo R2: {R_FORMAT_BY_KEY[(0.567, 5)]} ",
        f"Levene: {R_FORMAT_BY_KEY[(1234.5, 5)]} *",
    ]


def test_disstree2dotp_quality_digits_do_not_follow_the_node_digits(tmp_path, monkeypatch):
    # R's disstree2dotp passes `digits` to the node labels only.
    monkeypatch.chdir(tmp_path)
    tree = _tree()
    disstree2dotp(tree, "tree", digits=2, with_quality=True, use_title=False)

    quality = next(
        line
        for line in (tmp_path / "tree.dot").read_text().splitlines()
        if "node_qualityimage" in line
    )
    assert quality == (
        '"node_qualityimage"[shape=box, label="Global quality\\n'
        + "\\n".join(_quality_lines(tree))
        + '"];'
    )
