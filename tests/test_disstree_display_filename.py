"""disstreedisplay/seqtreedisplay accept R's image-path ``filename``.

What was wrong (2026-09-22 audit S6-18): the display functions always treated
``filename`` as a stem and appended the format, so the R-style call
``disstreedisplay(tree, "tree.png")`` wrote ``tree.png.png`` and
``tree.png.dot``.

Source of truth: TraMineR 2.2-14 ``disstree2dot.R`` -- ``disstreedisplay`` and
``seqtreedisplay`` copy the rendered image to ``filename`` itself
(``file.copy(..., filename)``). traminex now does the same when ``filename``
already ends with the image format's extension, and keeps its documented stem
behaviour otherwise (a name without that extension gains it, and the ``.dot``
file is written beside the image).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import load_golden

import traminex.disstree_viz as viz
from traminex.disstree import disstree, disstreedisplay
from traminex.format import seqdef
from traminex.seqtree import seqtree, seqtreedisplay


def _predictors():
    i = np.arange(1, 16)
    return pd.DataFrame(
        {
            "A": ["a", "a", "b", "b", "a", "b", "b", "b", "c", "c", "c", "c", "c", "a", "c"],
            "X": ((i * 7) % 15) + 0.1 * i,
        }
    )


def _tree():
    return disstree(load_golden("phase2_input_D"), _predictors(), R=1, min_size=2, max_depth=4)


def _seqtree():
    states = np.array([[["A", "B", "C"][(r + c) % 3] for c in range(8)] for r in range(15)])
    x = seqdef(states, alphabet=["A", "B", "C"])
    return seqtree(x, _predictors(), diss=load_golden("phase2_input_D"), R=1, min_size=2)


@pytest.fixture
def rendered(monkeypatch):
    calls = []

    def fake_render(dot_path, out_path, image_format):
        calls.append((Path(dot_path), Path(out_path), image_format))
        Path(out_path).write_bytes(b"image")

    monkeypatch.setattr(viz, "_render_dot", fake_render)
    return calls


@pytest.mark.parametrize("display", [disstreedisplay, seqtreedisplay])
@pytest.mark.parametrize(
    ("filename", "image_format", "image", "dot"),
    [
        ("tree.png", "png", "tree.png", "tree.dot"),
        ("tree.PNG", "png", "tree.PNG", "tree.dot"),
        ("tree.svg", "svg", "tree.svg", "tree.dot"),
        ("tree", "png", "tree.png", "tree.dot"),
        ("tree.v2", "png", "tree.v2.png", "tree.v2.dot"),
    ],
)
def test_image_path_follows_r_when_filename_has_the_extension(
    tmp_path, rendered, display, filename, image_format, image, dot
):
    tree = _tree() if display is disstreedisplay else _seqtree()

    output = display(tree, str(tmp_path / filename), image_format=image_format)

    assert Path(output) == tmp_path / image
    assert rendered == [(tmp_path / dot, tmp_path / image, image_format)]
    assert sorted(path.name for path in tmp_path.iterdir()) == sorted({image, dot})
