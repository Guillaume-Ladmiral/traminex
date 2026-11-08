"""R4 parity tests for TraMineRextras ``plot.seqimplic``."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from corestats_fixtures import fixture_a, fixture_b
from r_goldens_utils import R_GOLDEN_ATOL, R_GOLDENS_DIR
from scipy.stats import norm

import traminex
from traminex.stats_p2 import SeqImplic, seqimplic


def _plot(*args, **kwargs):
    return getattr(traminex, "plot_seqimplic")(*args, **kwargs)


def _series(case: str, level: str, state: str) -> pd.DataFrame:
    records = pd.read_csv(
        R_GOLDENS_DIR / "seqimplic_plot_r4_series.csv",
        na_values=["NA"],
        keep_default_na=True,
    )
    return records.loc[
        (records["case"] == case)
        & (records["level"].astype(str) == str(level))
        & (records["state"] == state)
    ].sort_values("position_index_zero")


def _metadata(case: str, key: str) -> str:
    records = pd.read_csv(
        R_GOLDENS_DIR / "seqimplic_plot_r4_metadata.csv", dtype=str, keep_default_na=False
    )
    return records.loc[(records["case"] == case) & (records["key"] == key), "value"].item()


def _condition(case: str) -> str:
    records = pd.read_csv(
        R_GOLDENS_DIR / "seqimplic_plot_r4_conditions.csv",
        dtype=str,
        keep_default_na=False,
    )
    return records.loc[records["case"] == case, "value"].item()


def _panel_axes(result):
    return list(np.asarray(result[1], dtype=object).reshape(-1))


def _state_lines(axis, state_count: int):
    return axis.lines[:state_count]


def test_plotted_series_match_r_sign_reversal_and_negative_clipping():
    result = seqimplic(fixture_a(), [1, 1, 2, 2])
    assert result.levels == _condition("default_levels").split("\x1e")
    fig, axes = _plot(result, conf_level=None)
    panels = _panel_axes((fig, axes))

    for axis, level in zip(panels, result.levels):
        for line, state in zip(_state_lines(axis, len(result.alphabet)), result.alphabet):
            expected = _series("A", level, state)
            np.testing.assert_allclose(
                line.get_ydata(),
                expected["plotted"],
                rtol=0,
                atol=R_GOLDEN_ATOL,
                equal_nan=True,
            )
            np.testing.assert_array_equal(line.get_xdata(), np.arange(1, 7))
    plt.close(fig)


def test_seqimplic_additive_plot_metadata_matches_source_stslist():
    x = fixture_a()
    x.cpal = ["#110000", "#220000", "#330000", "#440000"]
    x.xtstep = 2
    x.tick_last = True
    result = seqimplic(x, [1, 1, 2, 2])

    assert result.type == "indice"
    assert result.cpal == x.cpal
    assert result.xtstep == 2
    assert result.tick_last is True
    assert result.positions == [1, 2, 3, 4, 5, 6]


def test_seqimplic_preserves_custom_time_labels_from_r_metadata():
    x = fixture_a()
    expected = _metadata("custom_time_labels", "positions").split("\x1e")
    x.xtlab = expected

    result = seqimplic(x, [1, 1, 2, 2])
    assert result.positions == expected

    fig, axes = _plot(result, conf_level=None)
    axis = _panel_axes((fig, axes))[0]
    assert [tick.get_text() for tick in axis.get_xticklabels()] == expected
    plt.close(fig)


def test_default_titles_limits_colors_ticks_and_legend():
    result = seqimplic(fixture_a(), [1, 1, 2, 2])
    result.cpal = ["red", "green", "blue", "black"]
    result.xtstep = 2
    result.tick_last = True
    fig, axes = _plot(result, conf_level=None)
    panels = _panel_axes((fig, axes))

    assert [axis.get_title() for axis in panels] == result.levels
    expected_ylim = [float(value) for value in _metadata("A", "default_ylim").split("\x1e")]
    np.testing.assert_allclose(panels[0].get_ylim(), expected_ylim)
    assert [line.get_color() for line in _state_lines(panels[0], 4)] == result.cpal
    np.testing.assert_array_equal(panels[0].get_xticks(), [0.5, 2.5, 4.5, 5.5])
    assert [tick.get_text() for tick in panels[0].get_xticklabels()] == ["1", "3", "5", "6"]
    assert [text.get_text() for text in fig.legends[0].get_texts()] == result.labels
    plt.close(fig)


def test_missing_state_keeps_declared_state_label_and_palette_order():
    result = seqimplic(fixture_b(), [1, 1, 2, 2], with_missing=True)
    fig, axes = _plot(result, conf_level=None)
    axis = _panel_axes((fig, axes))[0]

    assert result.alphabet == _metadata("B_missing", "alphabet").split("\x1e")
    assert result.labels[-1] == "Missing"
    assert len(_state_lines(axis, len(result.alphabet))) == 4
    plt.close(fig)


@pytest.mark.parametrize("levels", [["2"], ["2", "1"], ["2", "2", "1"]])
def test_only_levels_selects_once_in_declared_level_order(levels):
    assert "at most once" in _condition("duplicate_levels")
    result = seqimplic(fixture_a(), [1, 1, 2, 2])
    fig, axes = _plot(result, only_levels=levels, conf_level=None)
    selected = [level for level in result.levels if level in levels]
    assert [axis.get_title() for axis in _panel_axes((fig, axes))] == selected
    plt.close(fig)


def test_unknown_and_empty_level_selection_are_actionable():
    assert "no panel" in _condition("unknown_levels")
    result = seqimplic(fixture_a(), [1, 1, 2, 2])
    for levels in (["unknown"], []):
        with pytest.raises(ValueError, match="only_levels"):
            _plot(result, only_levels=levels)


def test_indice_confidence_lines_use_normal_quantiles():
    assert "qnorm" in _condition("confidence_bounds")
    result = seqimplic(fixture_a(), [1, 1, 2, 2])
    fig, axes = _plot(result, conf_level=[0.95, 0.99], lwd=2)
    axis = _panel_axes((fig, axes))[0]
    confidence_lines = axis.lines[len(result.alphabet) :]

    np.testing.assert_allclose(
        [line.get_ydata()[0] for line in confidence_lines],
        [norm.ppf(0.95), norm.ppf(0.99)],
    )
    assert all(line.get_linestyle() == ":" for line in confidence_lines)
    assert all(line.get_linewidth() == 2 for line in confidence_lines)
    plt.close(fig)


def test_intensity_compatibility_uses_unit_limits_and_raw_confidence():
    result = SeqImplic(
        indices=np.array([[[0.2, -0.4, 0.8]]]),
        levels=["group"],
        alphabet=["A"],
        positions=["t1", "t2", "t3"],
        labels=["Alpha"],
        type="intensity",
        cpal=["purple"],
    )
    fig, axis = _plot(result, conf_level=0.8)

    assert axis.get_ylim() == (0.0, 1.0)
    assert axis.lines[1].get_ydata()[0] == pytest.approx(0.8)
    plt.close(fig)


def test_all_nan_indice_falls_back_to_unit_limits():
    result = SeqImplic(
        indices=np.full((1, 1, 3), np.nan),
        levels=["g"],
        alphabet=["A"],
        positions=[1, 2, 3],
    )
    fig, axis = _plot(result, conf_level=None)
    assert axis.get_ylim() == (0.0, 1.0)
    assert np.isnan(axis.lines[0].get_ydata()).all()
    plt.close(fig)


@pytest.mark.parametrize("conf_level", [0, 1, -0.1, 1.1, np.nan, "bad"])
def test_invalid_confidence_levels_are_rejected(conf_level):
    with pytest.raises((TypeError, ValueError), match="conf_level"):
        _plot(seqimplic(fixture_a(), [1, 1, 2, 2]), conf_level=conf_level)


def test_layout_rows_cols_bottom_axes_and_supplied_axes():
    base = seqimplic(fixture_a(), [1, 1, 2, 2])
    result = SeqImplic(
        indices=np.concatenate([base.indices, base.indices[:1]], axis=0),
        levels=["one", "two", "three"],
        alphabet=base.alphabet,
        positions=base.positions,
        labels=base.labels,
        cpal=base.cpal,
    )
    fig, supplied = plt.subplots(2, 2)
    supplied_flat = list(supplied.reshape(-1)[:3])
    returned_fig, axes = _plot(
        result,
        rows=2,
        cols=2,
        axes="bottom",
        ax=supplied_flat,
        conf_level=None,
    )
    panels = _panel_axes((returned_fig, axes))

    assert returned_fig is fig
    assert panels == supplied_flat
    assert panels[0].get_xticks().size == 0
    assert panels[1].get_xticks().size == 0
    assert panels[2].get_xticks().size > 0
    plt.close(fig)


@pytest.mark.parametrize(
    ("case", "proportion"),
    [("layout_right_0.2", 0.2), ("layout_right_0.4", 0.4)],
)
def test_legend_prop_reserves_r_plot_fraction(case, proportion):
    result = seqimplic(fixture_a(), [1, 1, 2, 2])
    fig, supplied = plt.subplots(1, 2)
    before = [axis.get_position().bounds for axis in supplied]

    _plot(
        result,
        ax=list(supplied),
        with_legend="right",
        legend_prop=proportion,
        conf_level=None,
    )

    plot_fraction = float(_metadata(case, "plot_fraction"))
    for original, axis in zip(before, supplied):
        positioned = axis.get_position().bounds
        assert positioned[0] == pytest.approx(original[0] * plot_fraction)
        assert positioned[2] == pytest.approx(original[2] * plot_fraction)
    assert fig.get_layout_engine() is None
    plt.close(fig)


def test_axis_title_style_and_legend_controls_forward():
    result = seqimplic(fixture_a(), [1, 1, 2, 2])
    fig, axes = _plot(
        result,
        main="Typical",
        ylim=(0, 4),
        xaxis=False,
        yaxis=False,
        ylab=None,
        with_legend=False,
        lwd=3,
        alpha=0.4,
        conf_level=None,
    )
    panels = _panel_axes((fig, axes))

    assert panels[0].get_title() == "Typical - 1"
    assert panels[0].get_ylim() == (0.0, 4.0)
    assert panels[0].get_xticks().size == 0
    assert panels[0].get_yticks().size == 0
    assert panels[0].get_ylabel() == ""
    assert all(line.get_linewidth() == 3 for line in _state_lines(panels[0], 4))
    assert all(line.get_alpha() == 0.4 for line in _state_lines(panels[0], 4))
    assert not fig.legends
    plt.close(fig)


def test_invalid_result_shapes_metadata_and_layout_are_rejected():
    valid = seqimplic(fixture_a(), [1, 1, 2, 2])
    malformed = [
        SeqImplic(np.zeros((2, 4)), valid.levels, valid.alphabet, valid.positions),
        SeqImplic(np.zeros((1, 4, 6)), valid.levels, valid.alphabet, valid.positions),
        SeqImplic(np.zeros((2, 3, 6)), valid.levels, valid.alphabet, valid.positions),
        SeqImplic(np.zeros((2, 4, 5)), valid.levels, valid.alphabet, valid.positions),
    ]
    with pytest.raises(TypeError, match="SeqImplic"):
        _plot(object())
    for value in malformed:
        with pytest.raises(ValueError):
            _plot(value)
    with pytest.raises(ValueError, match="rows.*cols"):
        _plot(valid, rows=1, cols=1)


def test_renderer_does_not_mutate_result_and_is_top_level_identity_export():
    result = seqimplic(fixture_a(), [1, 1, 2, 2])
    original_indices = result.indices.copy()
    original_levels = list(result.levels)
    original_cpal = None if result.cpal is None else list(result.cpal)
    fig, _ = _plot(result)

    np.testing.assert_array_equal(result.indices, original_indices)
    assert result.levels == original_levels
    assert result.cpal == original_cpal
    from traminex.visualization.plots import implicative

    assert traminex.plot_seqimplic is implicative.plot_seqimplic
    assert "SeqImplic" not in traminex.__all__
    plt.close(fig)
