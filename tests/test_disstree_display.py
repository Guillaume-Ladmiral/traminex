"""Structure + graceful-fallback tests for the tree viz layer (RP2 / B-12, D4).

Per decision D4 the viz layer is gated by **structure tests, not numeric
goldens**: ``disstree2dotp`` adds a global-fit quality node, ``seqtree2dot`` adds
an alphabet legend node, and ``disstreedisplay`` / ``seqtreedisplay`` render via
GraphViz ``dot`` as an **optional** extra — degrading gracefully (clear error,
``.dot`` still written) when ``dot`` is absent. The suite must pass with GraphViz
not installed, so the render path is exercised via a monkeypatched locator.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import load_golden

import traminex.disstree_viz as viz
from traminex.disstree import disstree, disstree2dot, disstree2dotp, disstreedisplay
from traminex.disstree_viz import _fmt, _stars
from traminex.format import seqdef
from traminex.seqtree import seqtree, seqtree2dot, seqtreedisplay


def _predictors(n=15):
    i = np.arange(1, n + 1)
    return pd.DataFrame(
        {
            "A": ["a", "a", "b", "b", "a", "b", "b", "b", "c", "c", "c", "c", "c", "a", "c"],
            "X": ((i * 7) % 15) + 0.1 * i,
        }
    )


def _tree():
    D = load_golden("phase2_input_D")
    return disstree(D, _predictors(), R=1, min_size=2, max_depth=4)


def _seqtree():
    D = load_golden("phase2_input_D")
    wint = load_golden("phase2_input_weights_int").ravel()
    n = 15
    states = np.empty((n, 8), dtype=object)
    for r in range(1, n + 1):
        for c in range(1, 9):
            states[r - 1, c - 1] = ["A", "B", "C"][(r + c) % 3]
    x = seqdef(states, alphabet=["A", "B", "C"], weights=wint)
    return seqtree(x, _predictors(), diss=D, R=1, min_size=2, max_depth=4)


def _emit(tmp_path, monkeypatch, fn, *args, **kw) -> list[str]:
    monkeypatch.chdir(tmp_path)
    fn(*args, **kw)
    return (tmp_path / "tree.dot").read_text(encoding="utf-8").splitlines()


def test_stars_thresholds():
    assert _stars(0.0005) == "***"
    assert _stars(0.005) == "**"
    assert _stars(0.02) == "*"
    assert _stars(0.5) == ""
    assert _stars(0.0) == "NA"


def test_disstree2dotp_no_quality_matches_plain_topology(tmp_path, monkeypatch):
    # without the quality node, disstree2dotp == disstree2dot text path.
    monkeypatch.chdir(tmp_path)
    tree = _tree()
    disstree2dot(tree, "a", use_title=False)
    disstree2dotp(tree, "b", with_quality=False, use_title=False)
    a = (tmp_path / "a.dot").read_text().splitlines()
    b = (tmp_path / "b.dot").read_text().splitlines()
    # node names differ only by the filename prefix (a_* vs b_*); compare stripped.
    assert [ln.replace('"a', '"b') for ln in a] == b


def test_disstree2dotp_quality_node(tmp_path, monkeypatch):
    tree = _tree()
    lines = _emit(tmp_path, monkeypatch, disstree2dotp, tree, "tree", with_quality=True)
    qual = [ln for ln in lines if "node_qualityimage" in ln]
    assert len(qual) == 1
    stat = tree.adjustment.stat
    for name in ("Pseudo F", "Pseudo R2", "Levene"):
        # R formats the quality values with format(x, digits = getOption("digits") - 2),
        # i.e. 5 significant digits (disstree2dot.R); "%.3g" was pinned here until
        # the 2026-09-22 audit S6-03. _fmt is checked against live R in
        # tests/test_disstree_dot_number_format.py.
        val = _fmt(float(stat.loc[name, "t0"]), 5)
        assert f"{name}: {val} " in qual[0]


def test_seqtree2dot_legend_node(tmp_path, monkeypatch):
    tree = _seqtree()
    lines = _emit(tmp_path, monkeypatch, seqtree2dot, tree, "tree", with_legend=True)
    legend = [ln for ln in lines if "node_image.legend" in ln]
    assert len(legend) == 1
    assert "Legend: A, B, C" in legend[0]
    # topology (the tree nodes/edges) is still present
    assert any('"tree_none"->"tree_none_left"' in ln for ln in lines)


def test_seqtree2dot_requires_seqdata(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError):
        seqtree2dot(_tree(), "tree")  # a disstree has no tree.info['object']


@pytest.mark.parametrize(
    "image_data",
    [
        pd.DataFrame({"value": range(14)}),
        np.arange(14),
        list(range(14)),
    ],
)
def test_disstree2dot_rejects_image_data_with_wrong_row_count(tmp_path, monkeypatch, image_data):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError, match="15.*root"):
        disstree2dot(
            _tree(),
            "tree",
            image_fun=lambda subset, output_path: output_path.touch(),
            image_data=image_data,
        )
    assert not (tmp_path / "tree.dot").exists()


def test_disstree2dot_subsets_numpy_and_sequence_positionally(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    expected_leaf_positions = [[10, 12, 14], [0, 1, 4, 13], [2, 3, 5, 6, 7], [8, 9, 11]]

    numpy_calls = []

    def emit_numpy(subset, output_path):
        numpy_calls.append(subset.copy())
        output_path.touch()

    disstree2dot(
        _tree(),
        "numpy-tree",
        image_fun=emit_numpy,
        image_data=np.arange(30).reshape(15, 2),
        only_leaf=True,
    )
    assert [subset[:, 0].tolist() for subset in numpy_calls] == [
        [2 * i for i in positions] for positions in expected_leaf_positions
    ]

    sequence_calls = []

    def emit_sequence(subset, output_path):
        sequence_calls.append(subset)
        output_path.touch()

    data = [f"item-{i}" for i in range(15)]
    disstree2dot(
        _tree(),
        "sequence-tree",
        image_fun=emit_sequence,
        image_data=data,
        only_leaf=True,
    )
    assert sequence_calls == [[data[i] for i in positions] for positions in expected_leaf_positions]


def test_disstree2dot_requires_callback_to_create_regular_file(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(RuntimeError, match="image callback.*regular file"):
        disstree2dot(_tree(), "missing", image_fun=lambda subset, output_path: None)
    assert not (tmp_path / "missing.dot").exists()

    def directory_output(subset, output_path):
        output_path.mkdir()

    with pytest.raises(RuntimeError, match="image callback.*regular file"):
        disstree2dot(_tree(), "directory", image_fun=directory_output)
    assert not (tmp_path / "directory.dot").exists()


def test_disstree2dot_accepts_empty_regular_callback_files(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    disstree2dot(
        _tree(),
        "tree",
        image_fun=lambda subset, output_path: output_path.touch(),
        only_leaf=True,
    )
    assert (tmp_path / "tree.dot").is_file()


def test_disstree2dot_image_title_and_relative_reference_contract(tmp_path):
    output_dir = tmp_path / "nested"
    output_dir.mkdir()
    legend = output_dir / "legend panel.png"
    quality = output_dir / "quality.png"
    legend.touch()
    quality.touch()
    callback_paths: list[Path] = []

    def emit(subset, output_path):
        callback_paths.append(output_path)
        output_path.touch()

    disstree2dot(
        _tree(),
        str(output_dir / "tree"),
        image_fun=emit,
        only_leaf=True,
        use_title=True,
        image_legend=legend,
        image_quality=quality,
        legend_text="portable legend",
    )
    dot = (output_dir / "tree.dot").read_text(encoding="utf-8")
    assert callback_paths
    assert all(path.parent == output_dir for path in callback_paths)
    assert 'image="tree_none_left.png"' in dot
    assert 'xlabel="n: 3 s2: 0.706"' in dot
    assert any(line.endswith('tree_none"[shape=box,  label=" "];') for line in dot.splitlines())
    assert 'image="legend panel.png"' in dot
    assert 'image="quality.png"' in dot
    assert '"node_legend.text"[shape=box, label=<portable legend>]' in dot


def test_disstree2dot_auxiliary_nodes_do_not_require_callback(tmp_path, monkeypatch):
    lines = _emit(
        tmp_path,
        monkeypatch,
        disstree2dot,
        _tree(),
        "tree",
        use_title=False,
        legend_text="text",
        image_legend="legend.png",
        image_quality="quality.png",
    )
    assert any("node_image.legend" in line and 'image="legend.png"' in line for line in lines)
    assert any("node_qualityimage" in line and 'image="quality.png"' in line for line in lines)
    assert any("node_legend.text" in line and "label=<text>" in line for line in lines)


def test_disstreedisplay_graceful_without_graphviz(tmp_path, monkeypatch):
    monkeypatch.setattr(viz, "_locate_dot", lambda: None)
    with pytest.raises(RuntimeError, match="GraphViz"):
        disstreedisplay(_tree(), str(tmp_path / "tree"))
    # the .dot is still written even though rendering failed (graceful)
    assert (tmp_path / "tree.dot").exists()


def test_seqtreedisplay_graceful_without_graphviz(tmp_path, monkeypatch):
    monkeypatch.setattr(viz, "_locate_dot", lambda: None)
    with pytest.raises(RuntimeError, match="GraphViz"):
        seqtreedisplay(_seqtree(), str(tmp_path / "tree"))
    assert (tmp_path / "tree.dot").exists()


def test_default_display_paths_share_one_session_owned_temp_directory(tmp_path, monkeypatch):
    monkeypatch.setattr(viz.tempfile, "tempdir", str(tmp_path))

    def fake_render(dot_path, out_path, image_format):
        assert Path(dot_path).is_file()
        Path(out_path).touch()

    monkeypatch.setattr(viz, "_render_dot", fake_render)
    try:
        outputs = [
            disstreedisplay(_tree()),
            disstreedisplay(_tree()),
            seqtreedisplay(_seqtree()),
        ]
        directories = {Path(output).parent for output in outputs}

        assert len(directories) == 1
        assert directories == set(tmp_path.iterdir())
        assert all(Path(output).is_file() for output in outputs)
    finally:
        owner = getattr(viz, "_DISPLAY_TEMP_DIRECTORY", None)
        if owner is not None:
            owner.cleanup()
            viz._DISPLAY_TEMP_DIRECTORY = None


# TX-AUD-20261010-082: GVPATH may name GraphViz's installation directory, its bin directory or the dot executable, and
# without dot the error names GVPATH, the paths tried and the DOT file. A stand-in executable plays dot.
def _stand_in_dot(prefix: Path) -> Path:
    executable = prefix / "bin" / "dot"
    executable.parent.mkdir(parents=True)
    executable.write_text("#!/bin/sh\nexit 0\n")
    executable.chmod(0o755)
    return executable


@pytest.mark.parametrize("form", ["prefix", "bin", "executable"])
def test_gvpath_names_the_installation_its_bin_directory_or_dot(tmp_path, monkeypatch, form):
    executable = _stand_in_dot(tmp_path / "graphviz")
    setting = {"prefix": tmp_path / "graphviz", "bin": executable.parent, "executable": executable}[
        form
    ]
    monkeypatch.setenv("GVPATH", str(setting))
    monkeypatch.setenv("PATH", str(tmp_path / "empty"))
    assert viz._locate_dot() == str(executable)


def test_the_missing_dot_error_names_gvpath_the_paths_tried_and_the_dot_file(tmp_path, monkeypatch):
    prefix = tmp_path / "graphviz"
    prefix.mkdir()
    monkeypatch.setenv("GVPATH", str(prefix))
    monkeypatch.setenv("PATH", str(tmp_path / "empty"))
    with pytest.raises(RuntimeError) as excinfo:
        disstreedisplay(_tree(), str(tmp_path / "tree"))
    message = str(excinfo.value)
    assert f"GVPATH={str(prefix)!r}" in message
    assert str(prefix / "bin" / "dot") in message and str(prefix / "dot") in message
    assert str(tmp_path / "tree.dot") in message
    assert (tmp_path / "tree.dot").exists()
