"""seqplot passes R's ``weighted`` option to every plot type.

What was wrong (2026-09-22 audit S6-20): ``seqplot(type="d"/"Ht"/"dH"/"f",
weighted=False)`` raised ``TypeError`` (the routed renderer has no ``weighted``
parameter), while TraMineR honours it.

Source of truth: TraMineR 2.2-14 ``seqplot.R`` hands every option that the
series function accepts to it (``fargs <- olist[names(olist) %in% flist]``);
``seqstatd`` (types d, Ht, dH) and ``seqtab`` (type f) both take ``weighted``,
so ``weighted = FALSE`` draws the unweighted series. The expected series are
traminex's R-golden'd ``seqstatd``/``seqtab`` with ``weighted=False``.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest
from corestats_fixtures import fixture_a

from traminex import seqstatd, seqtab
from traminex.visualization.plots.dispatch import seqplot


def _bar_heights(ax):
    return np.column_stack([[patch.get_height() for patch in c] for c in ax.containers])


@pytest.mark.parametrize("weighted", [False, True])
@pytest.mark.parametrize("plot_type", ["d", "dH"])
def test_distribution_types_follow_weighted(plot_type, weighted):
    x = fixture_a(weighted=True)

    fig, ax = seqplot(x, type=plot_type, weighted=weighted)

    expected = seqstatd(x, weighted=weighted)["Frequencies"]
    np.testing.assert_allclose(_bar_heights(ax), expected, rtol=0, atol=1e-12)
    plt.close(fig)


@pytest.mark.parametrize("weighted", [False, True])
def test_entropy_type_follows_weighted(weighted):
    x = fixture_a(weighted=True)

    fig, ax = seqplot(x, type="Ht", weighted=weighted)

    expected = np.asarray(seqstatd(x, weighted=weighted)["Entropy"], dtype=float)
    np.testing.assert_allclose(ax.lines[0].get_ydata(), expected, rtol=0, atol=1e-12)
    plt.close(fig)


@pytest.mark.parametrize("weighted", [False, True])
def test_frequency_type_follows_weighted(weighted):
    x = fixture_a(weighted=True)

    fig, ax = seqplot(x, type="f", weighted=weighted)

    expected = seqtab(x, idxs=0, normalize=True, weighted=weighted)["prop"].to_numpy()[:10]
    heights = sorted({round(patch.get_height(), 12) for patch in ax.patches}, reverse=True)
    assert heights == sorted(set(np.round(expected, 12)), reverse=True)
    plt.close(fig)


def test_weights_attached_to_the_caller_are_untouched():
    x = fixture_a(weighted=True)
    weights = x.weights.copy()

    fig, _ = seqplot(x, type="d", weighted=False)

    np.testing.assert_array_equal(x.weights, weights)
    plt.close(fig)
