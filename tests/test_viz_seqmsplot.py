"""RP7a: seqmsplot (plot.stslist.modst) structure test.

One bar per time position: height = the modal-state frequency, colour = the modal
state (R prof.freq[modal, t] = Frequencies[t]; barplot space=0, ylim 0..1.2).
Asserted against the seqmodst_A_cf golden (row 0 = modal codes, row 1 = freqs).
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import numpy as np
from corestats_fixtures import fixture_a
from r_goldens_utils import R_GOLDEN_ATOL, load_golden

from traminex.visualization.palette import ensure_cmap
from traminex.visualization.plots.summary import _state_color, seqmsplot


def test_seqmsplot_bar_heights_match_modal_frequencies():
    x = fixture_a()
    fig, ax = seqmsplot(x)
    golden = load_golden("seqmodst_A_cf")  # [modal codes; modal freqs]
    heights = np.array([p.get_height() for p in ax.patches])
    assert len(heights) == 6  # one bar per time position
    np.testing.assert_allclose(heights, golden[1], rtol=0, atol=R_GOLDEN_ATOL)


def test_seqmsplot_bar_colours_match_modal_states():
    x = fixture_a()
    fig, ax = seqmsplot(x)
    golden = load_golden("seqmodst_A_cf")
    modal_codes = golden[0].astype(int)
    cm = ensure_cmap(x, None)
    expected = [tuple(np.round(_state_color(cm, c)[:3], 6)) for c in modal_codes]
    got = [tuple(np.round(p.get_facecolor()[:3], 6)) for p in ax.patches]
    assert got == expected


def test_seqmsplot_ylim_zero_to_1p2():
    x = fixture_a()
    fig, ax = seqmsplot(x)
    assert ax.get_ylim() == (0.0, 1.2)
