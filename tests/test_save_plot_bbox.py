"""save_plot honours the ``bbox_inches="standard"`` option it documents.

What was wrong (2026-09-22 audit S6-19): the docstring offered
``bbox_inches="standard"`` for default margins, but the value went straight to
``Figure.savefig``, which raised ``AttributeError`` (matplotlib 3.11).

Source of truth: the function's own docstring -- "standard" keeps matplotlib's
default margins, i.e. the full figure canvas that ``bbox_inches=None`` saves.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.image import imread

from traminex.plotting import save_plot


def _saved_shape(tmp_path, name, **options):
    fig, ax = plt.subplots(figsize=(4, 3))
    ax.plot([0, 1], [0, 1])
    path = tmp_path / name
    save_plot(fig, str(path), dpi=20, **options)
    plt.close(fig)
    return imread(path).shape[:2]


def test_standard_keeps_the_full_canvas(tmp_path):
    assert _saved_shape(tmp_path, "standard.png", bbox_inches="standard") == (60, 80)
    assert _saved_shape(tmp_path, "none.png", bbox_inches=None) == (60, 80)


def test_tight_remains_the_default(tmp_path):
    assert _saved_shape(tmp_path, "tight.png") != (60, 80)
