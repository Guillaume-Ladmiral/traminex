"""seqeordplot disposition (D-ord): deprecation shim forwarding to seqpcplot.

seqeordplot was suppressed upstream in TraMineRextras 0.4.3 (2018-08-22) and
superseded by seqpcplot; traminex provides no analytical port, only a thin
forwarder that warns and delegates to seqpcplot. These tests fix that the
DeprecationWarning fires and that a figure (the seqpcplot result) is returned.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import warnings

import matplotlib.pyplot as plt
import pandas as pd
import pytest
from matplotlib.colors import to_rgba

import traminex
from traminex.events import seqecreate


def _eseq():
    tse = pd.DataFrame(
        {
            "id": ["s1", "s1", "s2", "s2", "s3"],
            "timestamp": [0, 2, 0, 1, 3],
            "event": ["A", "B", "B", "A", "C"],
        }
    )
    return seqecreate(tse)


def test_seqeordplot_is_public_and_top_level():
    assert traminex.seqeordplot is traminex.plotting.seqeordplot


def test_seqeordplot_warns_and_returns_figure():
    from traminex.plotting import seqeordplot

    with pytest.warns(DeprecationWarning, match="seqpcplot"):
        fig, ax = seqeordplot(_eseq())
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_seqeordplot_forwards_to_seqpcplot():
    # The shim's output equals seqpcplot's: same number of polylines drawn.
    from traminex.plotting import seqeordplot, seqpcplot

    e = _eseq()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        fig_ord, ax_ord = seqeordplot(e)
    fig_pc, ax_pc = seqpcplot(e)
    assert len(ax_ord.lines) == len(ax_pc.lines) == len(e)
    plt.close(fig_ord)
    plt.close(fig_pc)


def test_seqeordplot_warns_and_forwards_palette():
    from traminex.plotting import seqeordplot

    with pytest.warns(DeprecationWarning, match="seqpcplot"):
        fig, ax = seqeordplot(_eseq(), cpal=("red", "blue"))

    assert [line.get_color() for line in ax.lines] == [
        to_rgba("red"),
        to_rgba("blue"),
        to_rgba("red"),
    ]
    plt.close(fig)
