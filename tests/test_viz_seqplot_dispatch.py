"""RP7a: seqplot(type=...) dispatcher structure test.

Mirrors seqplot.R: each type routes to the matching renderer. Asserted by
checking the routed artists equal those of the per-type plot it dispatches to.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pytest
from corestats_fixtures import fixture_a
from matplotlib.collections import QuadMesh
from r_goldens_utils import R_GOLDEN_ATOL, load_golden

from traminex.format import seqdef
from traminex.visualization.plots.dispatch import seqplot

# Same 5-sequence fixture as tests/test_viz_seqfplot.py (-> seqtab_S golden).
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


def test_seqplot_type_i_is_index_plot():
    x = fixture_a()
    fig, ax = seqplot(x, type="i")
    assert np.array_equal(np.asarray(ax.images[0].get_array()), x.data)


@pytest.mark.parametrize("plot_type", ["i", "I"])
def test_seqplot_index_routes_forward_sampling_keywords(plot_type):
    x = fixture_a(weighted=True)
    fig, ax = seqplot(
        x,
        type=plot_type,
        sampling=2,
        sample_meth="random",
        random_state=24,
    )
    # Weighted rows are drawn as a mesh with R's heights (TX-AUD-20261019-048).
    drawn = (
        ax.images[0] if ax.images else next(c for c in ax.collections if isinstance(c, QuadMesh))
    )
    assert np.asarray(drawn.get_array()).size == 2 * x.seqlen
    plt.close(fig)


def test_seqplot_type_d_is_state_distribution():
    x = fixture_a()
    fig, ax = seqplot(x, type="d")
    heights = np.column_stack([[p.get_height() for p in c] for c in ax.containers])
    np.testing.assert_allclose(heights, load_golden("seqstatd_A_freq"), rtol=0, atol=R_GOLDEN_ATOL)


def test_seqplot_type_mt_is_mean_time():
    x = fixture_a()
    fig, ax = seqplot(x, type="mt")
    heights = np.array([p.get_height() for p in ax.patches])
    np.testing.assert_allclose(
        heights, load_golden("seqmeant_A_ms")[:, 0], rtol=0, atol=R_GOLDEN_ATOL
    )


def test_seqplot_type_f_is_seqfplot():
    # type="f" must route to the faithful seqfplot (Rectangle state-bands, NOT the
    # bar-container seqfrequency_plot) AND the band thicknesses must equal the
    # seqtab_S proportions -- bound to the golden so a heights bug is caught
    # (a bare patches>0/containers==0 check is vacuous to a heights mutation).
    x = seqdef(_S, alphabet=["A", "B", "C"])
    fig, ax = seqplot(x, type="f")
    assert len(ax.containers) == 0  # bands, not a seaborn/bar renderer
    groups = defaultdict(list)
    for p in ax.patches:
        groups[round(float(p.get_y()), 9)].append(p)
    heights = [float(groups[y0][0].get_height()) for y0 in sorted(groups)]
    expected = load_golden("seqtab_S")[:, 1] / 100.0  # desc props
    np.testing.assert_allclose(heights, expected, rtol=0, atol=R_GOLDEN_ATOL)


def test_seqplot_unknown_type_raises():
    with pytest.raises(ValueError):
        seqplot(fixture_a(), type="zzz")
