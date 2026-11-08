"""Structure golden for the R-faithful ``disstree2dot`` (RP2 / T20).

Per decision D4 the tree-viz layer is pinned by a **structure test**, not a
numeric golden: the emitted GraphViz ``.dot`` (topology + node/edge labels) must
match TraMineR 2.2-14 ``disstree2dot`` on the text/topology path
(``image.fun=NULL``, ``use.title=FALSE``), plus the R2 portable image callback
protocol (structural DOT and callback manifest, never raster pixels). Goldens by
``scripts/r_goldens/make_goldens_disstree_dot.R`` (filename root ``"tree"``).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from r_goldens_utils import load_golden, load_golden_lines

from traminex.disstree import disstree, disstree2dot


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


def _emit(tmp_path, monkeypatch, **kw) -> list[str]:
    monkeypatch.chdir(tmp_path)
    disstree2dot(_tree(), "tree", use_title=False, **kw)
    return (tmp_path / "tree.dot").read_text(encoding="utf-8").splitlines()


def test_disstree2dot_plain_matches_r(tmp_path, monkeypatch):
    assert _emit(tmp_path, monkeypatch) == load_golden_lines("disstree_dot_plain", ".dot")


def test_disstree2dot_show_depth_matches_r(tmp_path, monkeypatch):
    assert _emit(tmp_path, monkeypatch, show_depth=True) == load_golden_lines(
        "disstree_dot_depth", ".dot"
    )


def _manifest(case: str) -> pd.DataFrame:
    rows = pd.read_csv(Path(__file__).with_name("r_goldens") / "disstree_dot_r2_image_calls.csv")
    return rows.loc[rows["case"] == case].reset_index(drop=True)


def test_disstree2dot_explicit_image_calls_and_dot_match_r(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    image_data = pd.DataFrame(
        {
            "row_id": [f"row-{i:02d}" for i in range(1, 16)],
            "value": np.arange(1, 16) * 10,
        }
    )
    calls: list[tuple[pd.DataFrame, object]] = []

    def emit(subset, output_path):
        calls.append((subset.copy(), output_path))
        output_path.write_bytes(b"portable image")

    disstree2dot(
        _tree(),
        "tree",
        image_fun=emit,
        image_data=image_data,
        image_format="png",
        use_title=False,
        show_depth=True,
        legend_text="<B>R2 legend</B>",
        image_legend="legend.png",
        image_quality="quality.png",
    )

    expected = _manifest("explicit_all")
    assert len(calls) == len(expected)
    assert [type(path).__name__ for _, path in calls] == ["PosixPath"] * len(expected)
    assert [path.name for _, path in calls] == expected["image_filename"].tolist()
    assert [";".join(subset["row_id"]) for subset, _ in calls] == expected["subset_values"].tolist()
    assert (tmp_path / "tree.dot").read_text(encoding="utf-8").splitlines() == (
        load_golden_lines("disstree_dot_r2_images", ".dot")
    )


def test_disstree2dot_default_data_matches_one_based_r_manifest(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    tree = _tree()

    def make_one_based(node):
        node.info["ind"] = np.asarray(node.info["ind"]) + 1
        if node.split is not None:
            make_one_based(node.kids[0])
            make_one_based(node.kids[1])

    make_one_based(tree.root)
    calls = []

    def emit(subset, output_path):
        calls.append((subset.copy(), output_path.name))
        output_path.touch()

    disstree2dot(
        tree,
        "tree",
        image_fun=emit,
        image_format="jpg",
        use_title=True,
        show_depth=True,
    )

    expected = _manifest("default_all")
    assert [name for _, name in calls] == expected["image_filename"].tolist()
    assert [";".join(map(str, subset.iloc[:, 0])) for subset, _ in calls] == expected[
        "subset_values"
    ].tolist()


def test_disstree2dot_only_leaf_matches_r_call_order(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    image_data = pd.DataFrame(
        {
            "row_id": [f"row-{i:02d}" for i in range(1, 16)],
            "value": np.arange(15),
        }
    )
    calls = []

    def emit(subset, output_path):
        calls.append((subset.copy(), output_path.name))
        output_path.touch()

    disstree2dot(
        _tree(),
        "tree",
        image_fun=emit,
        image_data=image_data,
        only_leaf=True,
        image_format="png",
        use_title=True,
    )

    expected = _manifest("explicit_leaf")
    assert [name for _, name in calls] == expected["image_filename"].tolist()
    assert [";".join(subset["row_id"]) for subset, _ in calls] == expected["subset_values"].tolist()
