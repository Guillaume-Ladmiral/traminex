"""Graphviz calls and DOT text survive unusual file and column names.

What was wrong (2026-09-22 audit S6-17):

(a) ``disstreedisplay``/``seqtreedisplay`` passed the output and ``.dot`` paths to
    ``dot`` as given, so a relative name starting with ``-`` was read as an
    option: ``disstreedisplay(tree, "-Tsvg")`` failed with "Missing argument for
    -o flag" and left a stray ``-Tsvg.dot`` behind.
(b) Column names, alphabet labels and image paths went into DOT quoted strings
    unescaped, so a ``"`` in a predictor name produced a DOT syntax error.

Source of truth: Graphviz treats any argument starting with ``-`` as an option,
and in a DOT quoted string ``\\"`` is the only escape (a bare ``"`` ends the
string). TraMineR's ``disstree2dot.R`` pastes names raw too, so for such names
R writes invalid DOT; traminex escapes them. Names without ``"`` produce
byte-identical DOT (the R-golden DOT tests are unchanged).
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import load_golden

import traminex.disstree_viz as viz
from traminex.disstree import disstree, disstree2dot, disstreedisplay


def _tree(column: str = "X"):
    D = load_golden("phase2_input_D")
    i = np.arange(1, 16)
    predictors = pd.DataFrame(
        {
            "A": ["a", "a", "b", "b", "a", "b", "b", "b", "c", "c", "c", "c", "c", "a", "c"],
            column: ((i * 7) % 15) + 0.1 * i,
        }
    )
    return disstree(D, predictors, R=1, min_size=2, max_depth=4)


def _real_dot(monkeypatch):
    """Use a real ``dot``: found by traminex, or beside the running Python."""
    if viz._locate_dot() is None:
        environment = Path(sys.executable).resolve().parent.parent
        if (environment / "bin" / "dot").is_file():
            monkeypatch.setenv("GVPATH", os.fspath(environment))
    if viz._locate_dot() is None:
        pytest.skip("GraphViz dot is not available")


def test_render_command_never_passes_an_option_like_path(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    calls = []

    def fake_run(arguments, **kwargs):
        calls.append(list(arguments))
        Path(arguments[3]).touch()
        return type("Completed", (), {"returncode": 0, "stdout": "", "stderr": ""})()

    monkeypatch.setattr(viz, "_locate_dot", lambda: "/fake/dot")
    monkeypatch.setattr(viz.subprocess, "run", fake_run)

    disstreedisplay(_tree(), "-Tsvg")

    (arguments,) = calls
    assert arguments[:3] == ["/fake/dot", "-Tpng", "-o"]
    assert all(not argument.startswith("-") for argument in arguments[3:])
    assert Path(arguments[3]) == tmp_path / "-Tsvg.png"
    assert Path(arguments[4]) == tmp_path / "-Tsvg.dot"


def test_a_leading_dash_filename_renders_with_real_dot(tmp_path, monkeypatch):
    _real_dot(monkeypatch)
    monkeypatch.chdir(tmp_path)

    output = disstreedisplay(_tree(), "-Tsvg")

    assert (tmp_path / "-Tsvg.png").is_file()
    assert Path(output).resolve() == (tmp_path / "-Tsvg.png").resolve()


def test_double_quotes_in_names_are_escaped(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    tree = _tree('q"uote')

    disstree2dot(tree, "tree", use_title=False, image_legend='le"gend.png')

    dot = (tmp_path / "tree.dot").read_text()
    assert 'Split: q\\"uote R2:' in dot
    assert 'image="le\\"gend.png"' in dot
    assert 'q"uote' not in dot.replace('q\\"uote', "")


def test_double_quotes_in_names_give_valid_dot(tmp_path, monkeypatch):
    # The audit's reproduction: the text-label DOT of a tree split on a column
    # whose name contains a double quote made dot fail with a syntax error.
    _real_dot(monkeypatch)
    monkeypatch.chdir(tmp_path)
    disstree2dot(_tree('q"uote'), "quoted", use_title=False)

    viz._render_dot("quoted.dot", "quoted.svg", "svg")

    svg = (tmp_path / "quoted.svg").read_text()
    assert "q&quot;uote" in svg or 'q"uote' in svg


@pytest.mark.parametrize("basename", ['tree ü "root"', "tree & <root>", "plain"])
def test_nested_callback_images_are_embedded_by_real_dot(tmp_path, monkeypatch, basename):
    """Graphviz resolves image references relative to its working directory."""
    import xml.etree.ElementTree as ET

    plt = pytest.importorskip("matplotlib.pyplot")
    _real_dot(monkeypatch)
    monkeypatch.chdir(tmp_path)
    directory = tmp_path / "nested panels"
    stem = directory / basename
    images = []

    def emit(subset, path):
        plt.imsave(path, np.zeros((2, 2, 3)))
        images.append(path)

    disstree2dot(_tree(), str(stem), image_fun=emit, only_leaf=True)
    output = directory / "rendered.svg"
    viz._render_dot(str(stem) + ".dot", str(output), "svg")

    embedded = ET.parse(output).findall(".//{http://www.w3.org/2000/svg}image")
    assert images
    assert len(embedded) == len(images)
    assert all(
        (output.parent / element.attrib["{http://www.w3.org/1999/xlink}href"]).is_file()
        for element in embedded
    )
    assert all(path.is_file() for path in images)
    assert Path.cwd() == tmp_path
