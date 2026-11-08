"""RP7a: seqmtplot (plot.stslist.meant) structure test.

Bars = seqmeant Mean (mean time per state); 95% CI whiskers =
Mean +/- qt(.975, n-1) * SE (R plot.stslist.meant). Asserted against the
seqmeant_A_ms golden (cols [Mean, Var, Stdev, SE]).
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import numpy as np
from corestats_fixtures import fixture_a
from r_goldens_utils import R_GOLDEN_ATOL, load_golden
from scipy.stats import t as student_t

from traminex.visualization.plots.summary import seqmtplot


def test_seqmtplot_bar_heights_match_seqmeant_mean():
    x = fixture_a()
    fig, ax = seqmtplot(x)
    golden = load_golden("seqmeant_A_ms")  # [Mean, Var, Stdev, SE]
    # bars are Rectangle patches (ordered left-to-right == state order)
    heights = np.array([p.get_height() for p in ax.patches])
    assert len(heights) == 4  # one bar per state A,B,C,D
    np.testing.assert_allclose(heights, golden[:, 0], rtol=0, atol=R_GOLDEN_ATOL)


def test_seqmtplot_ci_whiskers_match_t_interval():
    from matplotlib.container import ErrorbarContainer

    x = fixture_a()
    # R draws the whiskers only when asked (serr = TRUE; TX-AUD-20261010-043).
    fig, ax = seqmtplot(x, errbar=True)
    golden = load_golden("seqmeant_A_ms")
    means, se = golden[:, 0], golden[:, 3]
    q = float(student_t.ppf(0.975, x.nseq - 1))
    # the errorbar container's vertical bar-segments span Mean +/- q*SE
    err = [c for c in ax.containers if isinstance(c, ErrorbarContainer)][0]
    segs = err.lines[2][0].get_segments()
    lows = np.array([s[0, 1] for s in segs])
    highs = np.array([s[1, 1] for s in segs])
    np.testing.assert_allclose(lows, means - q * se, rtol=0, atol=R_GOLDEN_ATOL)
    np.testing.assert_allclose(highs, means + q * se, rtol=0, atol=R_GOLDEN_ATOL)


def test_seqmtplot_ylim_is_zero_to_length():
    x = fixture_a()
    fig, ax = seqmtplot(x)
    assert ax.get_ylim() == (0.0, float(x.data.shape[1]))
