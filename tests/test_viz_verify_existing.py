"""RP7a VERIFY: existing traminex plots draw the R-golden'd series.

Structure tests (NOT pixel baselines) keyed to the already-pinned series. Each
test loads the *same* golden series the plot consumes and asserts the matplotlib
artists/backing data map to it (atol=1e-10). This de-risks the wrong-target trap
(the RP3-RP6 pattern) before any new plot is ported on top of these.

VERIFY verdicts (vs live TraMineR 2.2-14 plot structure, read from
SOURCES/TraMineR/R/{plot.stslist,plot.stslist.statd,plot.stslist.freq,seqlegend}.R):

- sequence_index_plot  <-> seqIplot/plot.stslist        : FAITHFUL (raw state grid;
  R draws coloured bars, traminex imshow -- same data->artist). sortv added later
  via the seqplot dispatcher helpers.
- state_distribution_plot <-> seqdplot/plot.stslist.statd: FAITHFUL (stacked bars =
  R barplot(seqstatd$Frequencies)).
- chronogram             <-> (heatmap of seqstatd)       : correct SERIES, a
  traminex-original heatmap rendering; the faithful seqdplot artist is
  state_distribution_plot.
- state_legend           <-> seqlegend                   : FAITHFUL (one swatch per
  alphabet state, colour=cpal, label=state label -- the alphabet here, as the
  fixture has no labels; "missing" and void handling: test_state_legend_labels).
- sequence_index_plot_grouped / cluster_chronograms      : FAITHFUL per-group series.
- seqfrequency_plot      <-> seqfplot/plot.stslist.freq   : draws the correct series
  (seqtab freq, desc) but as a simplified bar chart. The faithful coloured
  state-segment seqfplot is ADDED as a separate symbol (maintainer decision,
  2026-06-21); this test pins that seqfrequency_plot consumes seqtab.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import numpy as np
from corestats_fixtures import ALPH_A, SEQS_A, fixture_a
from r_goldens_utils import R_GOLDEN_ATOL, load_golden

from traminex.format import seqdef
from traminex.visualization.legend import state_legend
from traminex.visualization.palette import ensure_cmap
from traminex.visualization.plots.chronogram import chronogram, cluster_chronograms
from traminex.visualization.plots.distribution import state_distribution_plot
from traminex.visualization.plots.frequency import seqfrequency_plot
from traminex.visualization.plots.index import (
    sequence_index_plot,
    sequence_index_plot_grouped,
)


def _img(ax):
    """The backing array of the first imshow on an axes."""
    return np.asarray(ax.images[0].get_array())


# --- sequence_index_plot <-> seqIplot/plot.stslist -------------------------


def test_index_plot_draws_raw_state_matrix():
    x = fixture_a()
    fig, ax = sequence_index_plot(x)
    arr = _img(ax)
    # rows == n sequences, columns == sequence length L
    assert arr.shape == x.data.shape == (4, 6)
    # the heatmap IS the raw state-code matrix
    assert np.array_equal(arr, x.data)
    # distinct colours == alphabet size (one state code per alphabet symbol)
    assert len(np.unique(arr)) == len(x.alphabet) == 4


# --- state_distribution_plot <-> seqdplot/plot.stslist.statd ----------------


def test_state_distribution_stacked_heights_match_seqstatd_golden():
    x = fixture_a()
    fig, ax = state_distribution_plot(x, mode="d")
    # one stacked BarContainer per alphabet state (band count == alphabet)
    containers = ax.containers
    assert len(containers) == len(x.alphabet) == 4
    # heights[t, k] == seqstatd Frequencies (L x K), the pinned series
    heights = np.column_stack([[p.get_height() for p in c] for c in containers])
    assert heights.shape == (6, 4)
    np.testing.assert_allclose(heights, load_golden("seqstatd_A_freq"), rtol=0, atol=R_GOLDEN_ATOL)


# --- chronogram (heatmap of the seqstatd distribution) ---------------------


def test_chronogram_heatmap_is_seqstatd_distribution():
    x = fixture_a()
    fig, ax = chronogram(x)  # matplotlib (non-seaborn) path -> imshow
    arr = _img(ax)
    # K states x L positions == seqstatd Frequencies transposed
    np.testing.assert_allclose(arr, load_golden("seqstatd_A_freq").T, rtol=0, atol=R_GOLDEN_ATOL)


# --- state_legend <-> seqlegend --------------------------------------------


def test_state_legend_one_swatch_per_state():
    import matplotlib.pyplot as plt

    x = fixture_a()
    fig, ax = plt.subplots()
    patches, labels = state_legend(x, ax=ax)
    # one swatch per alphabet state, labels == alphabet
    assert len(patches) == len(x.alphabet) == 4
    assert list(labels) == ALPH_A
    # swatch colours == the resolved palette (distinct, one per state)
    cm = ensure_cmap(x, None)
    expected = [cm.colors[i % len(cm.colors)] for i in range(len(x.alphabet))]
    got = [tuple(p.get_facecolor()[:3]) for p in patches]
    for g, e in zip(got, expected):
        np.testing.assert_allclose(g, e[:3], atol=1e-12)


# --- seqfrequency_plot <-> seqtab (series VERIFY; artist is simplified) ------


def test_seqfrequency_plot_bar_widths_match_seqtab_props_descending():
    seqs = np.array(
        [
            ["A", "A", "B", "C"],
            ["A", "A", "B", "C"],
            ["A", "B", "B", "C"],
            ["A", "A", "B", "C"],
            ["A", "B", "B", "C"],
        ],
        dtype=object,
    )
    x = seqdef(seqs, alphabet=["A", "B", "C"])
    fig, ax = seqfrequency_plot(x, seaborn=False)
    widths = np.array([p.get_width() for p in ax.containers[0]])
    # bind to the R golden (seqtab_S: [count, prop*100], desc) rather than an
    # in-test seqtab recompute -- the bars must equal the R value, not "whatever
    # seqtab returns".
    expected = load_golden("seqtab_S")[:, 1] / 100.0
    np.testing.assert_allclose(widths, expected, rtol=0, atol=R_GOLDEN_ATOL)
    assert np.all(np.diff(expected) <= 1e-12)  # descending order


# --- grouped variants: per-group series -------------------------------------


def test_grouped_index_plot_sorts_rows_by_cluster():
    x = fixture_a()
    labels = np.array([0, 1, 0, 1])
    fig, ax, order = sequence_index_plot_grouped(x, labels)
    # stable sort groups sequences by cluster label (the sortv-by-cluster order)
    assert list(order) == [0, 2, 1, 3]
    arr = _img(ax)
    assert np.array_equal(arr, x.data[order])


def test_cluster_chronograms_one_panel_per_cluster_with_group_series():
    x = fixture_a()
    labels = np.array([0, 1, 0, 1])
    fig, axes = cluster_chronograms(x, labels)
    panels = [a for a in np.asarray(axes).ravel() if a.images]
    assert len(panels) == 2  # one chronogram per unique cluster label
    # cluster 0 == sequences 0 and 2; its panel == seqstatd of that subset
    from traminex.stats import seqstatd

    sub0 = seqdef(np.array([SEQS_A[0], SEQS_A[2]], dtype=object), alphabet=ALPH_A)
    expected0 = np.asarray(seqstatd(sub0, include_missing=False)["Frequencies"]).T
    np.testing.assert_allclose(_img(panels[0]), expected0, rtol=0, atol=R_GOLDEN_ATOL)
