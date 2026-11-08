"""Tree displays show R's node text when no node image is drawn (TX-AUD-20261010-014).

R's disstreedisplay always draws a node image and writes the node text into its
title: the branch label and "n: .. s2: .." in ``main``, "Split: .. R2:.." in
``sub`` (disstree2dot.R:365-383); the DOT label is then " ". traminex draws no
node image unless one is asked for, so it left every node box empty. A node
without an image now carries R's title text as its DOT label. Fixture:
``scripts/r_goldens/make_goldens_disstree_node_labels.R`` (R's title() calls).
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
from r_goldens_utils import load_golden

import traminex as tx
from traminex.disstree import disstree, disstree2dot
from traminex.disstree_viz import seqtree2dot

TITLES = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "disstree_node_labels.csv",
    dtype=str,
    keep_default_na=False,
)
NODE = re.compile(r'^"(?P<name>[^"]+)"\[shape=box,\s+label="(?P<label>.*)"\];$')


def _tree():
    D = load_golden("phase2_input_D")
    i = np.arange(1, 16)
    pred = pd.DataFrame(
        {
            "A": ["a", "a", "b", "b", "a", "b", "b", "b", "c", "c", "c", "c", "c", "a", "c"],
            "X": ((i * 7) % 15) + 0.1 * i,
        }
    )
    return disstree(D, pred, R=1, min_size=2, max_depth=4)


def _r_label(row: pd.Series) -> str:
    text = row["main"] if not row["sub"] else f"{row['main']}\n{row['sub']}"
    return text.replace("\n", "\\n")


def _node_labels(lines: list[str]) -> dict[str, str]:
    labels = {}
    for line in lines:
        match = NODE.match(line)
        if match and not match["name"].startswith("node_"):
            labels[match["name"]] = match["label"]
    return labels


def test_node_labels_are_r_titles(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    disstree2dot(_tree(), "tree")
    labels = _node_labels((tmp_path / "tree.dot").read_text(encoding="utf-8").splitlines())
    assert labels == {row["node_name"]: _r_label(row) for _, row in TITLES.iterrows()}


def test_seqtree_without_images_has_the_node_text(tmp_path, monkeypatch) -> None:
    # The same tree drawn as a sequence tree (15 sequences, one per case).
    monkeypatch.chdir(tmp_path)
    rows = [list("ABAB"), list("AABB"), list("BBAA")] * 5
    seqdata = tx.seqdef(np.array(rows, dtype=object), alphabet=["A", "B"])
    seqtree2dot(_tree(), "tree", seqdata=seqdata)
    labels = _node_labels((tmp_path / "tree.dot").read_text(encoding="utf-8").splitlines())
    assert labels == {row["node_name"]: _r_label(row) for _, row in TITLES.iterrows()}
