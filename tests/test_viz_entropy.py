"""RP7a: seqHtplot / seqdHplot (plot.stslist.statd type='Ht'/'dH') structure test.

seqHtplot draws the transversal-entropy line == seqstatd$Entropy (golden
seqstatd_A_ve row 1). seqdHplot = the stacked seqstatd distribution + that entropy
overlay (delegates to state_distribution_plot mode='dH').
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import numpy as np
from corestats_fixtures import fixture_a
from r_goldens_utils import R_GOLDEN_ATOL, load_golden

from traminex.visualization.plots.summary import seqdHplot, seqHtplot


def test_seqHtplot_line_is_seqstatd_entropy():
    x = fixture_a()
    fig, ax = seqHtplot(x)
    ent = load_golden("seqstatd_A_ve")[1]  # row 1 == Entropy
    (line,) = ax.lines
    np.testing.assert_allclose(line.get_ydata(), ent, rtol=0, atol=R_GOLDEN_ATOL)
    # plotted at time positions 1..L
    np.testing.assert_array_equal(line.get_xdata(), np.arange(1, len(ent) + 1))


def test_seqHtplot_grouped_draws_per_group_entropy_lines():
    # folded seqplot.tentrop: one entropy line per group == seqstatd of the subset
    from corestats_fixtures import ALPH_A, SEQS_A

    from traminex.format import seqdef
    from traminex.stats import seqstatd

    x = fixture_a()
    group = np.array([0, 1, 0, 1])
    fig, ax = seqHtplot(x, group=group)
    assert len(ax.lines) == 2
    for line, lev in zip(ax.lines, [0, 1]):
        rows = np.where(group == lev)[0]
        sub = seqdef(np.array([SEQS_A[r] for r in rows], dtype=object), alphabet=ALPH_A)
        exp = np.asarray(seqstatd(sub, include_missing=False)["Entropy"], dtype=float)
        np.testing.assert_allclose(line.get_ydata(), exp, rtol=0, atol=R_GOLDEN_ATOL)


def test_seqdHplot_has_distribution_bars_and_entropy_overlay():
    x = fixture_a()
    fig, ax = seqdHplot(x)
    # stacked bars == seqstatd Frequencies (4 state bands)
    containers = ax.containers
    assert len(containers) == 4
    heights = np.column_stack([[p.get_height() for p in c] for c in containers])
    np.testing.assert_allclose(heights, load_golden("seqstatd_A_freq"), rtol=0, atol=R_GOLDEN_ATOL)
    # the entropy overlay line == seqstatd Entropy
    ent = load_golden("seqstatd_A_ve")[1]
    (line,) = ax.lines
    np.testing.assert_allclose(line.get_ydata(), ent, rtol=0, atol=R_GOLDEN_ATOL)
