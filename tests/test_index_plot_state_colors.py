"""Every index-plot cell must be drawn in its own state's palette colour.

``imshow`` with a listed colormap but no value range stretches the palette over
whichever codes happen to be displayed.  When the rows on screen do not span
the whole codebook -- a sampled subset, an unobserved last state, missing
values present elsewhere in the object -- every state was drawn in a
neighbour's colour, e.g. state ``C`` in the missing-value grey (2026-09-22 audit
S6-01).  R's ``plot.stslist`` draws each state with ``cpal[state]``; the
legend and every bar renderer here assume the same mapping.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.colors import to_hex

from traminex import seqdef
from traminex.plotting import seqplot, sequence_index_plot, sequence_index_plot_grouped


def _drawn_colours(ax) -> dict[int, str]:
    image = ax.images[0]
    codes = np.asarray(image.get_array())
    rgba = image.to_rgba(codes)
    return {int(code): to_hex(rgba[codes == code][0]) for code in np.unique(codes)}


def _expected(x, codes) -> dict[int, str]:
    return {int(code): to_hex(x.cpal[int(code)]) for code in codes}


def _with_missing():
    rows = [
        ["A", "A", "B", "C"],
        ["A", "B", "B", "C"],
        ["A", None, "D", "D"],
        ["B", "B", "C", "C"],
    ]
    return seqdef(np.array(rows, dtype=object), alphabet=["A", "B", "C", "D"])


@pytest.mark.parametrize(
    "rows",
    [[0, 1, 3], [0, 1], [3]],
    ids=["no-missing-rows", "two-rows", "one-row-three-states"],
)
def test_a_displayed_subset_keeps_each_states_colour(rows):
    """Rows without missing values from an object whose codebook ends in ``*``."""
    x = _with_missing()
    assert x.nr_code == 4 and x.cpal is not None
    subset = x.copy()
    subset.data = x.data[rows]
    subset.weights = None
    fig, ax = sequence_index_plot(subset)
    try:
        drawn = _drawn_colours(ax)
        assert drawn == _expected(subset, drawn)
    finally:
        plt.close(fig)


def test_an_unobserved_last_state_does_not_shift_the_palette():
    x = seqdef(np.array([["A", "B", "B"], ["B", "A", "A"]], dtype=object), alphabet=["A", "B", "C"])
    x.cpal = ["#1b9e77", "#d95f02", "#7570b3"]
    fig, ax = sequence_index_plot(x)
    try:
        assert _drawn_colours(ax) == {0: "#1b9e77", 1: "#d95f02"}
    finally:
        plt.close(fig)


def test_a_single_displayed_state_is_drawn_in_its_own_colour():
    x = seqdef(np.array([["B", "B", "B"]], dtype=object), alphabet=["A", "B", "C"])
    x.cpal = ["#1b9e77", "#d95f02", "#7570b3"]
    fig, ax = sequence_index_plot(x)
    try:
        assert _drawn_colours(ax) == {1: "#d95f02"}
    finally:
        plt.close(fig)


def test_sampling_through_seqplot_keeps_each_states_colour():
    """The public route the audit reproduced: missing-bearing rows weighted out."""
    x = _with_missing()
    weights = np.array([0.0 if (x.data[r] == x.nr_code).any() else 1.0 for r in range(x.nseq)])
    x = seqdef(
        np.array(
            [[x.alphabet[c] if c != x.nr_code else None for c in row] for row in x.data],
            dtype=object,
        ),
        alphabet=["A", "B", "C", "D"],
        weights=weights,
    )
    fig, ax = seqplot(x, type="i", sampling=3, sample_meth="random", random_state=1)
    try:
        drawn = _drawn_colours(ax)
        assert x.nr_code not in drawn, "the weighted-out missing rows must not be sampled"
        assert drawn == _expected(x, drawn)
    finally:
        plt.close(fig)


def test_the_grouped_index_plot_keeps_each_states_colour():
    x = seqdef(np.array([["A", "B", "B"], ["B", "A", "A"]], dtype=object), alphabet=["A", "B", "C"])
    x.cpal = ["#1b9e77", "#d95f02", "#7570b3"]
    fig, ax, _order = sequence_index_plot_grouped(x, np.array([0, 1]))
    try:
        assert _drawn_colours(ax) == {0: "#1b9e77", 1: "#d95f02"}
    finally:
        plt.close(fig)
