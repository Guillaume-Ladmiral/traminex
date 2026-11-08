"""Legend colours match the drawn line colours when no palette is given.

What was wrong (2026-09-22 audit S6-10): with an ``STSList`` built without
``cpal`` (or ``SeqImplic(cpal=None)``), ``seqsplot`` and ``plot_seqimplic`` let
matplotlib's colour cycle colour the lines (C0, C1, C2, ...) but built legend
handles with ``color=None``, which matplotlib draws in its default line colour,
so every legend entry showed the first colour. In ``seqsplot``'s group mode the
cycle also restarted in each panel and skipped unobserved states, so one state
could change colour from panel to panel.

Source of truth: a legend is only correct if each entry shows the colour of the
curve it names. Without a palette, state curves now take the shared default
state palette every other renderer uses (``ensure_cmap``), and ``seqsplot``'s
per-state group curves beyond the eight-colour default take explicit cycle
colours.
"""

from __future__ import annotations

import dataclasses

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgba

from traminex.plotting import plot_seqimplic, seqsplot
from traminex.sequences import STSList
from traminex.stats_p2 import seqimplic
from traminex.visualization.palette import _state_color, ensure_cmap

DATA = np.array(
    [
        [0, 0, 1, 1, 2, 2],
        [0, 1, 1, 2, 2, 2],
        [1, 1, 1, 2, 3, 3],
        [0, 0, 0, 0, 1, 1],
        [2, 2, 2, 3, 3, 3],
        [0, 0, 1, 1, 2, 2],
        [3, 3, 3, 3, 3, 3],
        [0, 0, 1, 1, 2, 2],
    ]
)


def _no_palette(data=DATA):
    return STSList(np.array(data), alphabet=["A", "B", "C", "D"])


def _legend_colours(legend):
    handles = legend.legend_handles if hasattr(legend, "legend_handles") else legend.legendHandles
    return {
        text.get_text(): to_rgba(handle.get_color())
        for text, handle in zip(legend.get_texts(), handles)
    }


def test_seqsplot_group_mode_legend_matches_the_state_curves():
    x = _no_palette()
    # Group g2 never visits B, so its panel skips a state.
    group = np.array(["g1", "g1", "g2", "g1", "g2", "g1", "g2", "g1"])
    state_cmap = ensure_cmap(x, None)

    fig, axes = seqsplot(x, group=group)

    legend = _legend_colours(fig.legends[0])
    expected = {state: to_rgba(_state_color(state_cmap, code)) for code, state in enumerate("ABCD")}
    assert legend == expected
    for axis in axes:
        for line in axis.lines:
            assert to_rgba(line.get_color()) == expected[line.get_label()]
    plt.close(fig)


def test_seqsplot_per_state_legend_matches_group_curves_beyond_eight_groups():
    x = _no_palette(np.vstack([DATA] * 2)[:9])
    group = np.arange(9)

    fig, axes = seqsplot(x, group=group, per_state=True)

    legend = _legend_colours(fig.legends[0])
    assert len(set(legend.values())) == 9
    for axis in np.ravel(axes):
        for line in axis.lines:
            assert to_rgba(line.get_color()) == legend[line.get_label()]
    plt.close(fig)


def test_plot_seqimplic_legend_matches_the_lines_without_a_palette():
    x = _no_palette()
    result = dataclasses.replace(
        seqimplic(x, np.array(["g1", "g2"] * 4)),
        cpal=None,
    )

    fig, axes = plot_seqimplic(result)

    legend = _legend_colours(fig.legends[0])
    assert len(set(legend.values())) == 4
    for axis in np.ravel(axes):
        for line in axis.lines[:4]:
            assert to_rgba(line.get_color()) == legend[line.get_label()]
    plt.close(fig)
