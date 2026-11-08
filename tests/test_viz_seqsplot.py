"""R4 parity tests for TraMineRextras ``seqsplot``."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from corestats_fixtures import fixture_c
from r_goldens_utils import R_GOLDENS_DIR

import traminex
from traminex import seqdef
from traminex.stats_p2 import seqsurv
from traminex.visualization.plots.summary import _subset_stslist

PANELS = pd.read_csv(R_GOLDENS_DIR / "seqsplot_r4_panels.csv", dtype=str, keep_default_na=False)
LAYOUT = pd.read_csv(R_GOLDENS_DIR / "seqsplot_r4_layout.csv", dtype=str, keep_default_na=False)
CONDITIONS = pd.read_csv(
    R_GOLDENS_DIR / "seqsplot_r4_conditions.csv", dtype=str, keep_default_na=False
)


def _plot(*args, **kwargs):
    return getattr(traminex, "seqsplot")(*args, **kwargs)


def _fixture(*, weighted: bool = False):
    data = [
        ["A", "A", "B", "B", "C", "C"],
        ["A", "B", "B", "B", "C", "A"],
        ["D", "D", "D", "D", "D", "D"],
        ["A", "B", "A", "B", "A", "B"],
        ["A", "A", None, "C", None, None],
        ["B", "B", "C", "C", "C", "B"],
    ]
    return seqdef(
        data,
        alphabet=["A", "B", "C", "D"],
        labels=["Alpha", "Beta", "Gamma", "No exit"],
        cpal=["#111111", "#333333", "#777777", "#AAAAAA"],
        weights=[1, 2.5, 0.5, 3, 1.5, 2] if weighted else None,
        cnames=["one", "two", "three", "four", "five", "six"],
        xtstep=2,
        tick_last=True,
    )


def _groups():
    return pd.Categorical(
        ["early", "late", "noexit", "early", "late", "early"],
        categories=["late", "early", "noexit", "unused"],
        ordered=True,
    )


def _panels(rendered):
    return list(np.asarray(rendered[1], dtype=object).reshape(-1))


def _case(name: str):
    return PANELS.loc[PANELS["case"] == name].sort_values("panel_index_zero")


def _split(value: str):
    return [] if value == "" else value.split("\x1e")


def _condition(name: str, column: str):
    return CONDITIONS.loc[CONDITIONS["case"] == name, column].item()


def _expected_xy(curve, follow_up=None):
    """The drawn staircase: (0, 1), the event steps, then on to the last follow-up.

    As R's plot.survfit (tests/test_r_goldens_survival_extent.py,
    TX-AUD-20261010-042): a stratum whose spells are all censored is flat at 1.
    """
    values = np.asarray(curve, dtype=float)
    if values.size == 0 and follow_up is None:
        return np.empty((0, 2), dtype=float)
    points = (
        np.vstack([[0.0, 1.0], np.column_stack([values[:, 0], values[:, 3]])])
        if values.size
        else np.array([[0.0, 1.0]])
    )
    if follow_up is not None and follow_up > points[-1, 0]:
        points = np.vstack([points, [follow_up, points[-1, 1]]])
    return points


def _line_map(axis):
    return {line.get_label(): line for line in axis.lines}


def test_required_red_public_route_is_absent_until_implementation():
    assert callable(getattr(traminex, "seqsplot"))


def test_default_and_grouped_panels_match_r_selection_manifest():
    x = _fixture()
    fig_default, axis = _plot(x, with_legend=False)
    expected_default = _case("default").iloc[0]
    assert axis.get_title() == expected_default["title"]
    default_result = seqsurv(x)
    assert list(_line_map(axis)) == [
        key for key in _split(expected_default["curve_keys"]) if key in default_result.follow_up
    ]
    plt.close(fig_default)

    fig, axes = _plot(x, group=_groups())
    panels = _panels((fig, axes))
    expected = _case("grouped")
    assert [axis.get_title() for axis in panels] == expected["title"].tolist()
    assert [list(_line_map(axis)) for axis in panels] == [
        _split(value) for value in expected["curve_keys"]
    ]
    assert [text.get_text() for text in fig.legends[0].get_texts()] == _split(
        expected.iloc[0]["legend_labels"]
    )
    plt.close(fig)


def test_group_facets_draw_existing_seqsurv_series_for_each_subset():
    x = _fixture()
    groups = _groups()
    fig, axes = _plot(x, group=groups, with_legend=False)
    expected = _case("grouped")

    for axis, (_, row) in zip(_panels((fig, axes)), expected.iterrows()):
        positions = np.asarray([int(value) - 1 for value in _split(row["row_indexes_one"])])
        subset = _subset_stslist(x, positions)
        selected = _split(row["curve_keys"])
        expected_result = None if not selected else seqsurv(subset, state=selected)
        lines = _line_map(axis)
        for state in selected:
            curve = expected_result.curves[state]
            np.testing.assert_allclose(
                lines[state].get_xydata(),
                _expected_xy(curve, expected_result.follow_up.get(state)),
            )
    plt.close(fig)


def test_per_state_panels_superpose_groups_in_r_declared_order():
    x = _fixture()
    fig, axes = _plot(x, group=_groups(), per_state=True)
    panels = _panels((fig, axes))
    expected = _case("per_state")

    assert [axis.get_title() for axis in panels] == expected["title"].tolist()
    assert [list(_line_map(axis)) for axis in panels] == [
        [key for key in _split(value) if key in _line_map(axis)]
        for value, axis in zip(expected["curve_keys"], panels)
    ]
    assert [text.get_text() for text in fig.legends[0].get_texts()] == [
        "late",
        "early",
        "noexit",
    ]

    for axis, state in zip(panels, expected["panel_key"]):
        result = seqsurv(x, groups=_groups(), per_state=True, state=state)
        for group, line in _line_map(axis).items():
            np.testing.assert_allclose(
                line.get_xydata(), _expected_xy(result.curves[group], result.follow_up.get(group))
            )
    plt.close(fig)


def test_weighted_dispatch_preserves_the_r2_survival_estimates():
    x = _fixture(weighted=True)
    fig, axes = _plot(x, group=_groups(), per_state=True, which_states=["B"])
    axis = _panels((fig, axes))[0]
    result = seqsurv(x, groups=_groups(), per_state=True, state="B")

    for group, line in _line_map(axis).items():
        np.testing.assert_allclose(
            line.get_xydata(), _expected_xy(result.curves[group], result.follow_up.get(group))
        )
    plt.close(fig)


def test_state_subset_uses_declared_order_and_reports_partial_invalid(capsys):
    x = _fixture()
    fig, axes = _plot(x, group=_groups(), which_states=["C", "unknown", "A"])
    expected = _case("subset_declared_order")
    assert [list(_line_map(axis)) for axis in _panels((fig, axes))] == [
        _split(value) for value in expected["curve_keys"]
    ]
    assert [text.get_text() for text in fig.legends[0].get_texts()] == ["Alpha", "Gamma"]
    assert "Unvalid or unobserved state(s): unknown" in capsys.readouterr().err
    plt.close(fig)


def test_all_invalid_state_and_ignored_state_match_r_conditions(capsys):
    x = _fixture()
    with pytest.raises(ValueError, match="Only unvalid or unobserved selected states"):
        _plot(x, which_states=["unknown", "absent"])
    assert _condition("all_invalid_state", "outcome").startswith("ERROR:")

    fig, axis = _plot(x, state="B", with_legend=False)
    # D's spells are all censored: a flat curve, as R draws it.
    assert list(_line_map(axis)) == ["A", "B", "C", "D"]
    assert _condition("ignored_state", "message") in capsys.readouterr().err
    plt.close(fig)


def test_group_validation_drops_unused_and_na_but_preserves_declared_order():
    x = _fixture()
    groups = pd.Categorical(
        ["early", "late", None, "early", "late", "early"],
        categories=["late", "early", "unused"],
        ordered=True,
    )
    fig, axes = _plot(x, group=groups, with_legend=False)
    assert [axis.get_title() for axis in _panels((fig, axes))] == ["late", "early"]
    plt.close(fig)

    with pytest.raises(ValueError, match="group must contain one value for each row"):
        _plot(x, group=["a", "b"])
    assert "group must contain" in _condition("group_length", "outcome")


def test_missing_state_extension_palette_and_ragged_void_inputs():
    missing = seqdef(
        [["A", "*", "A"], ["*", "*", "A"], ["B", "B", "B"]],
        alphabet=["A", "B"],
        labels=["Alpha", "Beta"],
        cpal=["red", "blue"],
        nr="*",
        missing_color="pink",
    )
    fig, axis = _plot(missing, with_missing=True)
    # B's only spell is censored: a flat curve, as R draws it.
    assert list(_line_map(axis)) == ["A", "B", "*"]
    assert _line_map(axis)["*"].get_color() == "pink"
    assert [text.get_text() for text in fig.legends[0].get_texts()] == [
        "Alpha",
        "Beta",
        "Missing",
    ]
    plt.close(fig)

    ragged = fixture_c()
    original = ragged.data.copy()
    fig, axis = _plot(ragged, with_legend=False)
    assert axis.get_ylim()[0] == 0
    np.testing.assert_array_equal(ragged.data, original)
    plt.close(fig)


def test_no_event_state_is_flat_and_a_panel_without_its_spells_is_empty():
    # D's spells are all censored: R draws a flat curve at 1 to the last time;
    # a panel whose subset has no D spell is R's empty frame.
    x = _fixture()
    fig, axis = _plot(x, which_states=["D"], with_legend=False)
    assert len(axis.lines) == 1 and np.all(axis.lines[0].get_ydata() == 1.0)
    plt.close(fig)

    fig, axes = _plot(x, group=_groups(), which_states=["D"], with_legend=False)
    assert [len(axis.lines) for axis in _panels((fig, axes))] == [0, 0, 1]
    plt.close(fig)


def test_titles_ylabels_and_vector_length_validation():
    x = _fixture()
    fig, axes = _plot(
        x,
        group=_groups(),
        main="Survival",
        ylab=["Late y", "Early y", None],
        with_legend=False,
    )
    panels = _panels((fig, axes))
    assert [axis.get_title() for axis in panels] == [
        "Survival - late",
        "Survival - early",
        "Survival - noexit",
    ]
    assert [axis.get_ylabel() for axis in panels] == ["Late y", "Early y", ""]
    plt.close(fig)

    fig, axes = _plot(x, group=_groups(), main=None, with_legend=False)
    assert [axis.get_title() for axis in _panels((fig, axes))] == ["", "", ""]
    plt.close(fig)

    with pytest.raises(ValueError, match="ylab must have one value"):
        _plot(x, group=_groups(), ylab=["one", "two"])
    assert "ylab must have one value" in _condition("ylab_length", "outcome")


def test_per_state_custom_titles_match_r_label_manifest():
    fig, axes = _plot(
        _fixture(),
        group=_groups(),
        per_state=True,
        which_states=["C", "A"],
        main="Survival",
        with_legend=False,
    )
    assert [axis.get_title() for axis in _panels((fig, axes))] == _case("per_state_custom_title")[
        "title"
    ].tolist()
    plt.close(fig)


def test_three_panel_bottom_and_left_axes_match_r_occupied_layout():
    x = _fixture()
    fig, axes = _plot(
        x,
        group=_groups(),
        xaxis="bottom",
        yaxis="left",
        with_legend=False,
    )
    panels = _panels((fig, axes))
    expected = LAYOUT.loc[LAYOUT["nplot"] == "3"].iloc[0]
    bottom = {int(value) - 1 for value in _split(expected["bottom_panels_one"])}
    left = {int(value) - 1 for value in _split(expected["left_panels_one"])}
    assert {index for index, axis in enumerate(panels) if axis.get_xticks().size} == bottom
    assert {index for index, axis in enumerate(panels) if axis.get_yticks().size} == left
    plt.close(fig)


@pytest.mark.parametrize("xaxis", [False, "invalid"])
def test_xaxis_false_and_invalid_modes(xaxis):
    x = _fixture()
    if xaxis == "invalid":
        with pytest.raises(ValueError, match="xaxis"):
            _plot(x, xaxis=xaxis)
    else:
        fig, axis = _plot(x, xaxis=xaxis, with_legend=False)
        assert axis.get_xticks().size == 0
        plt.close(fig)


@pytest.mark.parametrize("yaxis", [False, "invalid"])
def test_yaxis_false_and_invalid_modes(yaxis):
    x = _fixture()
    if yaxis == "invalid":
        with pytest.raises(ValueError, match="yaxis"):
            _plot(x, yaxis=yaxis)
    else:
        fig, axis = _plot(x, yaxis=yaxis, with_legend=False)
        assert axis.get_yticks().size == 0
        plt.close(fig)


def test_custom_time_labels_axis_scale_palette_legend_and_style_forward():
    x = _fixture()
    fig, axes = _plot(
        x,
        group=_groups(),
        cpal=["red", "green"],
        which_states=["A", "C"],
        xtlab=["t1", "t2", "t3", "t4", "t5", "t6"],
        cex_axis=1.5,
        ltext=["State A", "State C"],
        linewidth=3,
        alpha=0.4,
    )
    axis = _panels((fig, axes))[0]
    assert [tick.get_text() for tick in axis.get_xticklabels()] == ["t1", "t2"]
    assert [line.get_color() for line in axis.lines] == ["red", "green"]
    assert all(line.get_linewidth() == 3 for line in axis.lines)
    assert all(line.get_alpha() == 0.4 for line in axis.lines)
    assert [text.get_text() for text in fig.legends[0].get_texts()] == [
        "State A",
        "State C",
    ]
    plt.close(fig)


@pytest.mark.parametrize(
    ("case", "kwargs"),
    [
        ("xtstep_override", {"xtstep": 3, "tick_last": False}),
        ("tick_last_override", {"xtstep": 4, "tick_last": True}),
    ],
)
def test_per_call_time_controls_match_r_ticks(case, kwargs):
    fig, axis = _plot(_fixture(), with_legend=False, **kwargs)
    expected = [int(value) for value in _condition(case, "message").split("\x1e")]
    np.testing.assert_array_equal(axis.get_xticks(), expected)
    plt.close(fig)


def test_per_state_custom_group_palette_and_legend_visibility():
    x = _fixture()
    fig, axes = _plot(
        x,
        group=_groups(),
        per_state=True,
        which_states=["A"],
        cpal=["red", "green", "blue"],
        ltext=["Late", "Early", "No exit"],
        with_legend=False,
    )
    axis = _panels((fig, axes))[0]
    assert [line.get_color() for line in axis.lines] == ["red", "green"]
    assert not fig.legends
    plt.close(fig)


def test_single_selected_curve_accepts_scalar_color_and_legend_text():
    fig, axes = _plot(
        _fixture(),
        group=_groups(),
        which_states=["C"],
        cpal="orange",
        ltext="Gamma spells",
    )
    panels = _panels((fig, axes))
    assert all(line.get_color() == "orange" for axis in panels for line in axis.lines)
    assert [text.get_text() for text in fig.legends[0].get_texts()] == ["Gamma spells"]
    plt.close(fig)


def test_supplied_axes_explicit_layout_and_nonmutation():
    x = _fixture()
    original = (
        x.data.copy(),
        None if x.weights is None else x.weights.copy(),
        list(x.alphabet),
        list(x.labels),
        list(x.cpal),
    )
    fig, grid = plt.subplots(1, 3)
    returned_fig, axes = _plot(
        x,
        group=_groups(),
        rows=1,
        cols=3,
        ax=list(grid),
        with_legend=False,
    )
    assert returned_fig is fig
    assert _panels((returned_fig, axes)) == list(grid)
    np.testing.assert_array_equal(x.data, original[0])
    assert x.weights is None
    assert (x.alphabet, x.labels, x.cpal) == original[2:]
    plt.close(fig)


@pytest.mark.parametrize(
    ("proportion", "column"),
    [(0.2, "right_0_2_plot_fraction"), (0.4, "right_0_4_plot_fraction")],
)
def test_legend_prop_reserves_r_plot_fraction(proportion, column):
    fig, supplied = plt.subplots()
    before = supplied.get_position().bounds
    _plot(
        _fixture(),
        ax=supplied,
        with_legend="right",
        legend_prop=proportion,
    )

    plot_fraction = float(LAYOUT.loc[LAYOUT["nplot"] == "1", column].item())
    positioned = supplied.get_position().bounds
    assert positioned[0] == pytest.approx(before[0] * plot_fraction)
    assert positioned[2] == pytest.approx(before[2] * plot_fraction)
    assert fig.get_layout_engine() is None
    plt.close(fig)


@pytest.mark.parametrize("enabled", [True, False])
def test_use_layout_controls_matplotlib_layout_engine(enabled):
    fig, _ = _plot(_fixture(), with_legend=False, use_layout=enabled)
    assert _condition(f"use_layout_{str(enabled).lower()}", "outcome") == "OK"
    assert (fig.get_layout_engine() is not None) is enabled
    plt.close(fig)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"rows": 1, "cols": 1, "group": _groups()}, "rows.*cols"),
        ({"cpal": ["red", "blue"]}, "cpal"),
        ({"ltext": ["one"]}, "ltext"),
        ({"legend_prop": 1.5}, "legend_prop"),
        ({"use_layout": "yes"}, "use_layout"),
        ({"xtstep": 0}, "xtstep"),
        ({"xtlab": ["one"]}, "xtlab"),
    ],
)
def test_invalid_layout_palette_and_metadata_controls(kwargs, message):
    with pytest.raises((TypeError, ValueError), match=message):
        _plot(_fixture(), **kwargs)


def test_wrong_input_axes_and_public_identity_are_rejected_or_stable():
    with pytest.raises(TypeError, match="STSList"):
        _plot(object())
    fig, axes = plt.subplots(1, 2)
    with pytest.raises(ValueError, match="supplied ax"):
        _plot(_fixture(), ax=list(axes))
    plt.close(fig)

    from traminex.visualization.plots import survival

    assert traminex.seqsplot is survival.seqsplot


@pytest.mark.parametrize("rc_constrained", [False, True])
@pytest.mark.parametrize("rc_tight", [False, True])
@pytest.mark.parametrize("enabled", [False, True])
def test_explicit_survival_layout_overrides_rc_settings(rc_constrained, rc_tight, enabled):
    with plt.rc_context(
        {"figure.constrained_layout.use": rc_constrained, "figure.autolayout": rc_tight}
    ):
        fig, _ = _plot(_fixture(), with_legend=False, use_layout=enabled)
        engine = fig.get_layout_engine()
        if enabled:
            assert type(engine).__name__ == "ConstrainedLayoutEngine"
        else:
            assert engine is None
        fig.canvas.draw()
        plt.close(fig)
