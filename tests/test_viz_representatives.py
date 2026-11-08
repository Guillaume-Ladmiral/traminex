"""RP7a: seqrplot (plot.stslist.rep) + seqrfplot (plot.seqrf) structure tests.

seqrplot draws the representative sequences as coloured bands, thickness = the
coverage share na% from seqrep/dissrep (golden dissrep_density_uw_stats).
seqrfplot draws (left) the RF-group medoid sequences as bands, thickness = group
heights, and (right) one box per group of the within-group residual distances --
the medoids/heights are golden'd in seqrf_uw_repcodes / seqrf_uw_heights, the
residuals come from dissrf.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.colors import to_rgba
from r_goldens_utils import R_GOLDEN_ATOL, load_golden

from traminex.format import seqdef
from traminex.visualization.palette import ensure_cmap
from traminex.visualization.plots import representatives
from traminex.visualization.plots.representatives import (
    _state_color,
    seqrfplot,
    seqrplot,
)


def _band_heights(ax):
    """Rectangle bands grouped bottom-up by y0; returns [(height, [colors])].

    Groups by a coarse y0 key but returns the full-precision band height so the
    golden comparison keeps atol=1e-10.
    """
    groups = defaultdict(list)
    for p in ax.patches:
        groups[round(float(p.get_y()), 6)].append(p)
    bands = []
    for y0 in sorted(groups):
        rects = sorted(groups[y0], key=lambda r: r.get_x())
        hs = [float(r.get_height()) for r in rects]
        assert max(hs) - min(hs) < 1e-9  # a band has uniform thickness
        colors = [tuple(np.round(r.get_facecolor()[:3], 6)) for r in rects]
        bands.append((hs[0], colors))
    return bands


# --- seqrplot ---------------------------------------------------------------


def test_seqrplot_band_heights_match_dissrep_coverage():
    D = load_golden("phase2_input_D")  # 15x15 phase-2 distance matrix
    n = D.shape[0]
    # codes do not affect rep selection (only D does); build a valid 15-seq STSList
    chars = np.array(
        [[["A", "B", "C", "D"][(i + j) % 4] for j in range(6)] for i in range(n)],
        dtype=object,
    )
    x = seqdef(chars, alphabet=["A", "B", "C", "D"])
    fig, ax = seqrplot(x, diss=D, criterion="density", coverage=0.5, weighted=False)
    stats = load_golden("dissrep_density_uw_stats")  # rows r1..rk,Total; col 1 == na(%)
    idx = load_golden("dissrep_density_uw_idx").ravel().astype(int)
    nrep = len(idx)
    expected = stats[:nrep, 1] / 100.0
    bands = _band_heights(ax)
    assert len(bands) == nrep
    np.testing.assert_allclose([b[0] for b in bands], expected, rtol=0, atol=R_GOLDEN_ATOL)


# --- seqrfplot --------------------------------------------------------------

_RF_ALPH = ["1", "2", "3", "4"]


def _seqrf_fixture():
    codes = load_golden("seqrf_input_codes").astype(int)
    w = load_golden("seqrf_input_weights").ravel()
    chars = np.array([[_RF_ALPH[c] for c in row] for row in codes], dtype=object)
    return seqdef(chars, alphabet=_RF_ALPH, weights=w)


def test_seqrfplot_left_medoid_bands_match_seqrf_golden():
    x = _seqrf_fixture()
    D = load_golden("seqrf_input_diss")
    sortv = load_golden("seqrf_input_sortv").ravel()
    fig, (axL, axR) = seqrfplot(x, D, k=3, sortv=sortv, weighted=False)
    heights = load_golden("seqrf_uw_heights").ravel()
    repcodes = load_golden("seqrf_uw_repcodes").astype(int)
    bands = _band_heights(axL)
    assert len(bands) == 3
    np.testing.assert_allclose([b[0] for b in bands], heights, rtol=0, atol=R_GOLDEN_ATOL)
    # band i cell colours == the i-th medoid sequence's state colours
    cm = ensure_cmap(x, None)
    for i, band in enumerate(bands):
        expected = [tuple(np.round(_state_color(cm, int(c))[:3], 6)) for c in repcodes[i]]
        assert band[1] == expected


def test_seqrfplot_right_boxes_match_weighted_box_golden():
    x = _seqrf_fixture()
    D = load_golden("seqrf_input_diss")
    sortv = load_golden("seqrf_input_sortv").ravel()
    fig, (axL, axR) = seqrfplot(x, D, k=3, sortv=sortv, weighted=False)
    # R plot.seqrf draws a WEIGHTED Tukey-hinge box (wtd.boxplot.tmr), NOT matplotlib's
    # unweighted np.percentile -- pinned per group as [whislo, q1, med, q3, whishi].
    box = load_golden("seqrf_uw_box")
    assert len(axR.patches) == 3
    # one box (PathPatch) per group; its x-extent == R's WEIGHTED [q1, q3] (the
    # current/old unweighted np.percentile box gives q1=2.5 for group 1, not 2.0).
    for patch, row in zip(axR.patches, box):
        xs = patch.get_path().vertices[:, 0]
        np.testing.assert_allclose(
            [xs.min(), xs.max()], [row[1], row[3]], rtol=0, atol=R_GOLDEN_ATOL
        )


# --- TraMineRextras seqplot.rf wrapper --------------------------------------


def _seqplot_rf(*args, **kwargs):
    return getattr(representatives, "seqplot_rf")(*args, **kwargs)


def test_seqplot_rf_both_matches_direct_unweighted_first_delegate():
    x = _seqrf_fixture()
    D = load_golden("seqrf_input_diss")
    sortv = load_golden("seqrf_input_sortv").ravel()
    wrapper_fig, (wrapper_left, wrapper_right) = _seqplot_rf(x, D, k=3, sortv=sortv)
    direct_fig, (direct_left, direct_right) = seqrfplot(
        x, D, k=3, sortv=sortv, weighted=False, grp_meth="first", yaxis=False
    )

    np.testing.assert_allclose(
        [band[0] for band in _band_heights(wrapper_left)],
        [band[0] for band in _band_heights(direct_left)],
    )
    assert len(wrapper_right.patches) == len(direct_right.patches)
    plt.close(wrapper_fig)
    plt.close(direct_fig)


@pytest.mark.parametrize("which_plot", ["medoids", "diss.to.med"])
def test_seqplot_rf_single_modes_draw_only_requested_panel(which_plot):
    x = _seqrf_fixture()
    D = load_golden("seqrf_input_diss")
    sortv = load_golden("seqrf_input_sortv").ravel()
    fig, ax = _seqplot_rf(x, D, k=3, sortv=sortv, which_plot=which_plot)

    if which_plot == "medoids":
        assert _band_heights(ax)
        assert "Sequences medoids" in ax.get_title()
    else:
        assert len(ax.patches) == 3
        assert "Dissimilarities to medoids" in ax.get_title()
    plt.close(fig)


def test_seqplot_rf_default_none_sort_maps_to_delegate_mds(monkeypatch):
    # seqplot_rf fits seqrf once itself (grp_meth/random_state go to the fit)
    # and hands the fitted result to seqrfplot via rf=, so no analytics travel
    # through figure attributes (cycle-8 audit F-47).
    fit_captured = {}
    plot_captured = {}
    real_seqrf = representatives.seqrf

    def fake_seqrf(*args, **kwargs):
        fit_captured.update(kwargs)
        return real_seqrf(*args, **kwargs)

    def fake_delegate(*args, **kwargs):
        plot_captured.update(kwargs)
        fig, axes = plt.subplots(1, 2)
        return fig, tuple(axes)

    monkeypatch.setattr(representatives, "seqrf", fake_seqrf)
    monkeypatch.setattr(representatives, "seqrfplot", fake_delegate)
    x = _seqrf_fixture()
    fig, _ = _seqplot_rf(x, load_golden("seqrf_input_diss"), k=3, sortv=None)
    assert fit_captured["sortv"] == "mds"
    assert fit_captured["weighted"] is False
    assert fit_captured["grp_meth"] == "first"
    assert plot_captured["sortv"] == "mds"
    assert plot_captured["weighted"] is False
    assert plot_captured["rf"] is not None
    plt.close(fig)


@pytest.mark.parametrize("which_plot", ["unknown", "distance", ""])
def test_seqplot_rf_rejects_unknown_plot_modes(which_plot):
    with pytest.raises(ValueError, match="which_plot.*both.*medoids.*diss.to.med"):
        _seqplot_rf(
            _seqrf_fixture(),
            load_golden("seqrf_input_diss"),
            k=3,
            which_plot=which_plot,
        )


def test_seqplot_rf_rejects_unknown_group_method():
    with pytest.raises(ValueError, match="grp_meth.*first.*random"):
        _seqplot_rf(_seqrf_fixture(), load_golden("seqrf_input_diss"), k=3, grp_meth="prop")


@pytest.mark.parametrize(
    "diss",
    [
        np.ones(4),
        np.ones((3, 4)),
        np.ones((2, 2)),
        np.full((12, 12), np.nan),
    ],
)
def test_seqplot_rf_validates_distance_shape_size_and_finiteness(diss):
    with pytest.raises(ValueError, match="diss"):
        _seqplot_rf(_seqrf_fixture(), diss, k=3)


@pytest.mark.parametrize("k", [0, -1, 2.5, True, 99])
def test_seqplot_rf_validates_group_count(k):
    with pytest.raises((TypeError, ValueError), match="k"):
        _seqplot_rf(_seqrf_fixture(), load_golden("seqrf_input_diss"), k=k)


def test_seqplot_rf_seeded_random_is_reproducible_and_global_rng_safe(monkeypatch):
    x = _seqrf_fixture()
    D = load_golden("seqrf_input_diss")
    sortv = load_golden("seqrf_input_sortv").ravel()
    fitted: list = []
    real_seqrf = representatives.seqrf

    def recording_seqrf(*args, **kwargs):
        result = real_seqrf(*args, **kwargs)
        fitted.append(result)
        return result

    monkeypatch.setattr(representatives, "seqrf", recording_seqrf)
    np.random.seed(991)
    expected_global = np.random.random(3)
    np.random.seed(991)
    first_fig, _ = _seqplot_rf(x, D, k=5, sortv=sortv, grp_meth="random", random_state=72)
    second_fig, _ = _seqplot_rf(x, D, k=5, sortv=sortv, grp_meth="random", random_state=72)

    np.testing.assert_array_equal(fitted[0]["rf"]["mdsk"], fitted[1]["rf"]["mdsk"])
    np.testing.assert_array_equal(np.random.random(3), expected_global)
    plt.close(first_fig)
    plt.close(second_fig)


def test_seqplot_rf_axes_titles_labels_and_messages(capsys):
    x = _seqrf_fixture()
    D = load_golden("seqrf_input_diss")
    sortv = load_golden("seqrf_input_sortv").ravel()
    fig, supplied = plt.subplots()
    returned_fig, ax = _seqplot_rf(
        x,
        D,
        k=3,
        sortv=sortv,
        which_plot="medoids",
        main="Careers",
        ylab="Frequency groups",
        yaxis=True,
        ax=supplied,
        border=False,
    )
    captured = capsys.readouterr()

    assert returned_fig is fig and ax is supplied
    assert ax.get_title() == "Careers: Sequences medoids"
    assert ax.get_ylabel() == "Frequency groups"
    assert [tick.get_text() for tick in ax.get_yticklabels()] == ["g1", "g2", "g3"]
    assert "Using k=3 frequency groups" in captured.err
    assert "Pseudo/median-based-R2" in captured.err
    assert "Pseudo/median-based-F statistic" in captured.err
    assert all(patch.get_linewidth() == 0 for patch in ax.patches)
    plt.close(fig)


def test_seqplot_rf_custom_palette_and_artist_style_are_forwarded():
    x = _seqrf_fixture()
    D = load_golden("seqrf_input_diss")
    sortv = load_golden("seqrf_input_sortv").ravel()
    palette = ["#CC0000", "#0044CC"]

    fig, ax = _seqplot_rf(
        x,
        D,
        k=3,
        sortv=sortv,
        which_plot="medoids",
        cpal=palette,
        alpha=0.45,
        linewidth=1.75,
    )

    first_code = int(load_golden("seqrf_uw_repcodes")[0, 0])
    expected_color = palette[first_code % len(palette)]
    assert ax.patches[0].get_facecolor() == pytest.approx(to_rgba(expected_color, 0.45))
    assert all(patch.get_alpha() == pytest.approx(0.45) for patch in ax.patches)
    assert all(patch.get_linewidth() == pytest.approx(1.75) for patch in ax.patches)
    plt.close(fig)


def test_seqplot_rf_does_not_mutate_inputs_and_is_top_level_identity_export():
    import traminex

    x = _seqrf_fixture()
    D = load_golden("seqrf_input_diss")
    sortv = load_golden("seqrf_input_sortv").ravel()
    original_data = x.data.copy()
    original_weights = x.weights.copy()
    original_diss = D.copy()
    original_sort = sortv.copy()
    fig, _ = _seqplot_rf(x, D, k=3, sortv=sortv)

    np.testing.assert_array_equal(x.data, original_data)
    np.testing.assert_array_equal(x.weights, original_weights)
    np.testing.assert_array_equal(D, original_diss)
    np.testing.assert_array_equal(sortv, original_sort)
    assert getattr(traminex, "seqplot_rf") is getattr(representatives, "seqplot_rf")
    plt.close(fig)
