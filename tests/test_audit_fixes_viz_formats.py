"""Regressions for the 2026-09-02 audit's visualization/format/surface findings."""

from __future__ import annotations

import warnings

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

import traminex
from traminex.sequences import STSList


def test_state_color_has_a_single_shared_implementation():
    from traminex.visualization import palette
    from traminex.visualization.plots import (
        alignment,
        frequency,
        representatives,
        summary,
        survival,
    )

    for module in (summary, frequency, representatives, survival, alignment):
        assert module._state_color is palette._state_color, module.__name__


def test_ensure_cmap_warns_when_a_declared_palette_is_too_short():
    from traminex.visualization.palette import ensure_cmap

    x = STSList(np.array([[0, 1, 2]]), alphabet=["A", "B", "C"], cpal=["#ff0000", "#0000ff"])
    with pytest.warns(UserWarning, match="cpal"):
        ensure_cmap(x, None)
    full = STSList(np.array([[0, 1]]), alphabet=["A", "B"], cpal=["#ff0000", "#0000ff"])
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        ensure_cmap(full, None)


def test_seqrfplot_returns_results_without_figure_attribute_smuggling():
    rng = np.random.default_rng(5)
    data = rng.integers(0, 3, size=(12, 6))
    x = STSList(data, alphabet=["A", "B", "C"])
    diss = traminex.seqdist(x, method="LCS")
    fig, _axes = traminex.seqrfplot(x, diss=diss, k=3)
    assert not hasattr(fig, "_traminex_seqrf")
    plt.close(fig)
    result = traminex.seqplot_rf(x, diss=diss, k=3)
    plt.close(result[0])


def test_compress_dss_decodes_external_markers():
    from traminex.seqformat import _compress_dss

    x = STSList(
        np.array([[0, 3, 1, 4]]),
        alphabet=["A", "B"],
        nr="*",
        void="%",
        nr_code=3,
        void_code=4,
    )
    out = _compress_dss(x)
    assert out.to_strings() == ["A-*-B"]


def test_materialize_sts_repair_is_code_keyed_for_both_marker_regimes():
    import warnings as w

    from traminex.seqformat import seqformat

    ext = STSList(
        np.array([[0, 3, 1, 4], [1, 1, 4, 4]]),
        alphabet=["A", "B"],
        labels=["Alpha", "Beta"],
        cpal=["#ff0000", "#0000ff"],
        nr="*",
        void="%",
        nr_code=3,
        void_code=4,
        missing_color="#00ff00",
    )
    with w.catch_warnings():
        w.simplefilter("ignore")
        out = seqformat(ext, from_="STS", to="DSS")
    assert out.labels == ["Alpha", "Beta", "*", "%"]
    assert out.cpal == ["#ff0000", "#0000ff", "#00ff00", "white"]


def test_stats_facade_docstring_names_all_source_modules():
    import traminex.stats as stats

    doc = stats.__doc__ or ""
    for module in ("stats_basic", "stats_complexity", "stats_utils", "stats_p1", "stats_p2"):
        assert module in doc, module


def test_dissrf_is_routed_from_its_home_module():
    import traminex as tx

    assert tx._LAZY_EXPORTS.get("repseq_rf") is not None
    assert "dissrf" in tx._LAZY_EXPORTS["repseq_rf"]
    from traminex.repseq_rf import dissrf as home

    assert tx.dissrf is home
