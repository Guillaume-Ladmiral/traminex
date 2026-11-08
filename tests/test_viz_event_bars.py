"""RP7a: plot_subseqelist / plot_subseqelistchisq structure tests.

plot_subseqelist bars == seqefsub Support (golden seqefsub_syn_k3). The
discriminating-subsequence plot draws one panel per group with bars ==
seqecmpgroup Freq.<group> (golden seqecmpgroup_g2 columns Freq.0 / Freq.1).
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
from r_goldens_utils import R_GOLDEN_ATOL, R_GOLDENS_DIR, load_golden, load_golden_lines

from traminex import seqecmpgroup, seqecreate, seqefsub
from traminex.visualization.plots.events_viz import plot_subseqelist, plot_subseqelistchisq


def _syn_eseq():
    tse = pd.DataFrame(
        {
            "id": ["s1", "s1", "s1", "s1", "s2", "s2", "s2", "s3", "s3", "s3"],
            "timestamp": [0, 1, 1, 3, 0, 2, 4, 0, 1, 2],
            "event": ["A", "B", "C", "A", "A", "B", "C", "A", "B", "C"],
        }
    )
    return seqecreate(tse)


def test_plot_subseqelist_bars_match_seqefsub_support():
    fsub = seqefsub(_syn_eseq(), min_support=1, max_k=3)
    fig, ax = plot_subseqelist(fsub)
    golden = load_golden("seqefsub_syn_k3")  # [Support, Count]
    widths = np.array([p.get_width() for p in ax.patches])
    assert len(widths) == len(fsub.subseq)
    np.testing.assert_allclose(widths, golden[:, 0], rtol=0, atol=R_GOLDEN_ATOL)


def _actcal_fsub():
    df = pd.read_csv(R_GOLDENS_DIR / "actcal_tse.csv")
    levels = load_golden_lines("eseq_actcal_levels")
    df["event"] = pd.Categorical(df["event"].astype(str), categories=levels)
    ea = seqecreate(df, time_col="time")
    return seqefsub(ea, pmin_support=0.05, max_k=2)


def test_plot_subseqelistchisq_panel_bars_match_group_freqs():
    fsub = _actcal_fsub()
    group = np.arange(1, len(fsub.eseq) + 1) % 2
    cmp = seqecmpgroup(fsub, group)
    fig, axes = plot_subseqelistchisq(cmp)
    cols = load_golden_lines("seqecmpgroup_g2_cols")
    golden = load_golden("seqecmpgroup_g2")
    assert len(axes) == len(cmp.labels) == 2
    for ax, g in zip(axes, cmp.labels):
        widths = np.array([p.get_width() for p in ax.patches])
        expected = golden[:, cols.index(f"Freq.{g}")]
        np.testing.assert_allclose(widths, expected, rtol=0, atol=R_GOLDEN_ATOL)
