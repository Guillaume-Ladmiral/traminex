"""RP7a: faithful seqfplot (plot.stslist.freq) structure test.

seqfplot draws the top distinct sequences as horizontal coloured state-segment
bands, thickness proportional to (weighted) frequency -- R
barplot(seqgbar(...), width=table$Percent, horiz=TRUE). The frequencies come from
seqtab (R-golden'd in seqtab_S); this asserts the band heights map to them and the
cell colours map to the state palette.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import numpy as np
from r_goldens_utils import R_GOLDEN_ATOL, load_golden

from traminex.format import seqdef
from traminex.visualization.palette import ensure_cmap
from traminex.visualization.plots.frequency import _state_color, seqfplot, seqfrequency_plot

# Same fixture as tests/test_r_goldens_seqtab.py (-> seqtab_S golden).
_S = np.array(
    [
        ["A", "A", "B", "C"],
        ["A", "A", "B", "C"],
        ["A", "B", "B", "C"],
        ["A", "A", "B", "C"],
        ["A", "B", "B", "C"],
    ],
    dtype=object,
)


def _fixture_s():
    return seqdef(_S, alphabet=["A", "B", "C"])


def _bands(ax):
    """Group rectangle patches into bands keyed by their y0; return [(y0, height, [colors])]."""
    from collections import defaultdict

    groups = defaultdict(list)
    for p in ax.patches:
        y0 = round(float(p.get_y()), 9)
        groups[y0].append(p)
    bands = []
    for y0 in sorted(groups):
        rects = sorted(groups[y0], key=lambda r: r.get_x())
        heights = {round(float(r.get_height()), 9) for r in rects}
        assert len(heights) == 1  # a band has uniform thickness
        colors = [tuple(np.round(r.get_facecolor()[:3], 6)) for r in rects]
        bands.append((y0, heights.pop(), colors))
    return bands


def test_seqfplot_band_heights_match_seqtab_props():
    x = _fixture_s()
    fig, ax = seqfplot(x)
    bands = _bands(ax)
    # two distinct sequences -> two bands; one cell rectangle per time position (L=4)
    assert len(bands) == 2
    assert all(len(b[2]) == 4 for b in bands)
    # band thickness (bottom-up == descending frequency) == seqtab proportions
    heights_bottom_up = [b[1] for b in bands]
    golden = load_golden("seqtab_S")  # columns [count, prop*100], desc
    expected_props = golden[:, 1] / 100.0
    np.testing.assert_allclose(heights_bottom_up, expected_props, rtol=0, atol=R_GOLDEN_ATOL)


def test_seqfplot_cell_colours_match_state_palette():
    x = _fixture_s()
    fig, ax = seqfplot(x)
    cm = ensure_cmap(x, None)
    # most frequent sequence is "A-A-B-C" (codes 0,0,1,2), drawn as the bottom band
    bottom = _bands(ax)[0]
    expected = [tuple(np.round(_state_color(cm, c)[:3], 6)) for c in (0, 0, 1, 2)]
    assert bottom[2] == expected


def test_seqfplot_top_limits_number_of_bands():
    x = _fixture_s()
    fig, ax = seqfplot(x, top=1)
    assert len(_bands(ax)) == 1


def test_seqfplot_top_none_includes_every_distinct_sequence():
    mat = np.array(
        [list(f"{value:04b}".replace("0", "A").replace("1", "B")) for value in range(12)],
        dtype=object,
    )
    x = seqdef(mat, alphabet=["A", "B"])

    fig_default, ax_default = seqfplot(x)
    fig_all, ax_all = seqfplot(x, top=None)

    assert len(_bands(ax_default)) == 10
    assert len(_bands(ax_all)) == 12

    fig_bars, ax_bars = seqfrequency_plot(x, top=None, seaborn=False)
    assert len(ax_bars.patches) == 12


def test_seqfplot_representatives_align_duplicates_weights_missing_and_void():
    x = seqdef(
        np.array(
            [
                ["A", "*", "B", "%"],
                ["A", "*", "B", "%"],
                ["B", "B", "%", "%"],
            ],
            dtype=object,
        ),
        alphabet=["A", "B"],
        nr="*",
        void="%",
        weights=[2.0, 3.0, 1.0],
    )

    fig, ax = seqfplot(x)
    bands = _bands(ax)
    np.testing.assert_allclose([band[1] for band in bands], [5 / 6, 1 / 6])

    cm = ensure_cmap(x, None)
    expected = [tuple(np.round(_state_color(cm, int(code))[:3], 6)) for code in x.data[0]]
    assert bands[0][2] == expected


def test_visualization_subpackage_all_is_not_the_public_plotting_contract():
    """cycle-6 §8.3: the numbers in traminex/visualization/__init__.py's comment.

    The comment there claimed the list was the "Public visualization API". It is
    not: `traminex.plotting` is, and the two differ in both directions. The
    membership is deliberately left alone -- nothing user-visible depends on it --
    so this test pins the figures the corrected comment quotes rather than the
    list itself.
    """
    import traminex.plotting as plotting
    import traminex.visualization as visualization
    from traminex.api_manifest import PUBLIC_API_SET

    public_plot_names = set(plotting.__all__) & PUBLIC_API_SET
    assert len(plotting.__all__) == 38
    assert len(public_plot_names) == 37
    assert len(visualization.__all__) == 23
    assert len(public_plot_names - set(visualization.__all__)) == 19
    assert set(visualization.__all__) - PUBLIC_API_SET == {
        "compute_even_ticks",
        "determine_layout",
        "figure_to_buffer",
        "legend_figure",
        "save_plot",
    }
