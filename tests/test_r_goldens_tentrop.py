"""R4 parity tests for TraMineRextras transversal-entropy plots."""

from __future__ import annotations

from collections import OrderedDict

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDEN_ATOL, R_GOLDENS_DIR

import traminex
from traminex.sequences import STSList
from traminex.visualization.plots import summary


def _fixture(*, weighted: bool = True) -> STSList:
    return STSList(
        data=np.array(
            [
                [0, 0, 1, 1, 2],
                [0, 1, 1, 2, 2],
                [1, 1, 0, 0, 0],
                [2, 2, 2, 2, 2],
                [0, 3, 1, 4, 4],
                [1, 0, 0, 0, 4],
            ],
            dtype=int,
        ),
        alphabet=["A", "B", "C"],
        labels=["A", "B", "C"],
        cpal=["#7FC97F", "#BEAED4", "#FDC086"],
        xtlab=["t1", "t2", "t3", "t4", "t5"],
        xtstep=2,
        tick_last=True,
        nr_code=3,
        void_code=4,
        weights=np.array([1.0, 2.0, 1.0, 3.0, 2.0, 1.0]) if weighted else None,
        start=15,
        ids=[f"person-{index}" for index in range(6)],
    )


def _group() -> pd.Categorical:
    return pd.Categorical(
        ["beta", "alpha", "beta", "gamma", None, "alpha"],
        categories=["unused", "gamma", "beta", "alpha"],
    )


def _golden_series(case: str, curve: str) -> pd.DataFrame:
    records = pd.read_csv(R_GOLDENS_DIR / "tentrop_r4_series.csv")
    return records.loc[(records["case"] == case) & (records["curve"] == curve)]


def _metadata(case: str, key: str) -> str:
    records = pd.read_csv(R_GOLDENS_DIR / "tentrop_r4_metadata.csv", dtype=str)
    return records.loc[(records["case"] == case) & (records["key"] == key), "value"].item()


def _condition(case: str) -> tuple[str, str]:
    records = pd.read_csv(
        R_GOLDENS_DIR / "tentrop_r4_conditions.csv", dtype=str, keep_default_na=False
    )
    kind, message = records.loc[records["case"] == case, "condition"].item().split("\x1e", 1)
    return kind, message


def _tentrop(*args, **kwargs):
    return getattr(summary, "seqplot_tentrop")(*args, **kwargs)


def _tentrop_m(*args, **kwargs):
    return getattr(summary, "seqplot_tentrop_m")(*args, **kwargs)


def test_grouped_entropy_lines_match_r_goldens_and_factor_order():
    fig, ax = _tentrop(_fixture(), _group())
    expected_names = _metadata("grouped", "curve_names").split("\x1e")

    assert [line.get_label() for line in ax.lines] == expected_names
    for line, name in zip(ax.lines, expected_names):
        expected = _golden_series("grouped", name)
        np.testing.assert_array_equal(line.get_xdata(), expected["position_one"])
        np.testing.assert_allclose(
            line.get_ydata(), expected["entropy"], rtol=0, atol=R_GOLDEN_ATOL
        )
    plt.close(fig)


def test_grouped_default_style_title_limits_ticks_and_legend_match_r():
    fig, ax = _tentrop(_fixture(), _group())

    assert ax.get_title() == _metadata("contract", "auto_title")
    np.testing.assert_allclose(
        ax.get_ylim(),
        [float(value) for value in _metadata("grouped", "ylim").split("\x1e")],
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )
    assert [line.get_color().upper() for line in ax.lines] == _metadata("grouped", "colors").split(
        "\x1e"
    )
    assert [line.get_linestyle() for line in ax.lines] == ["-", "--", ":"]
    assert [line.get_linewidth() for line in ax.lines] == [3.5, 3.5, 3.5]
    np.testing.assert_array_equal(ax.get_xticks(), [1, 3, 5])
    assert [tick.get_text() for tick in ax.get_xticklabels()] == ["t1", "t3", "t5"]
    assert [text.get_text() for text in ax.get_legend().get_texts()] == [
        "gamma",
        "beta",
        "alpha",
    ]
    plt.close(fig)


def test_multiple_mapping_lines_match_r_goldens_and_mapping_names():
    x = _fixture()
    curves = OrderedDict(
        [
            ("full", x),
            (
                "short",
                STSList(
                    data=x.data[:, :3],
                    alphabet=list(x.alphabet),
                    nr_code=x.nr_code,
                    void_code=x.void_code,
                    weights=x.weights,
                ),
            ),
            (
                "first_three",
                STSList(
                    data=x.data[:3],
                    alphabet=list(x.alphabet),
                    nr_code=x.nr_code,
                    void_code=x.void_code,
                    weights=x.weights[:3],
                ),
            ),
        ]
    )
    with pytest.warns(UserWarning, match="not all of same length"):
        fig, ax = _tentrop_m(curves)

    assert [line.get_label() for line in ax.lines] == list(curves)
    for line, name in zip(ax.lines, curves):
        expected = _golden_series("multiple", name)
        np.testing.assert_allclose(
            line.get_ydata(), expected["entropy"], rtol=0, atol=R_GOLDEN_ATOL
        )
    assert [line.get_color().upper() for line in ax.lines] == _metadata("multiple", "colors").split(
        "\x1e"
    )
    plt.close(fig)


def test_multiple_unequal_widths_warns_and_uses_longest_width():
    kind, oracle_message = _condition("multiple_unequal_widths")
    assert kind == "warning" and "not all of same length" in oracle_message
    x = _fixture()
    short = STSList(
        data=x.data[:, :2],
        alphabet=list(x.alphabet),
        nr_code=x.nr_code,
        void_code=x.void_code,
        weights=x.weights,
    )
    with pytest.warns(UserWarning, match="not all of same length"):
        fig, ax = _tentrop_m([x, short], with_legend=False)
    assert ax.get_xlim() == (1.0, 5.0)
    assert [len(line.get_xdata()) for line in ax.lines] == [5, 2]
    plt.close(fig)


@pytest.mark.parametrize(
    "group",
    [
        np.array([2, 1, 2, 3, 1, 3]),
        np.array(["z", "a", "z", "m", "a", "m"], dtype=object),
    ],
)
def test_plain_numeric_and_string_groups_have_factor_sorted_labels(group):
    fig, ax = _tentrop(_fixture(), group)
    assert [line.get_label() for line in ax.lines] == [str(value) for value in np.unique(group)]
    plt.close(fig)


def test_na_groups_are_omitted_and_all_na_is_rejected():
    fig, ax = _tentrop(_fixture(), ["a", None, "b", "a", None, "b"])
    assert [line.get_label() for line in ax.lines] == ["a", "b"]
    plt.close(fig)
    with pytest.raises(ValueError, match="at least one non-missing group"):
        _tentrop(_fixture(), [None] * 6)


def test_group_length_mismatch_is_actionable():
    kind, oracle_message = _condition("group_length_mismatch")
    assert kind == "error" and "same length" in oracle_message
    with pytest.raises(ValueError, match="group.*sequence"):
        _tentrop(_fixture(), ["a", "b"])


@pytest.mark.parametrize("value", [_fixture(), [_fixture()], (), [object(), _fixture()]])
def test_multiple_route_rejects_non_list_short_or_bad_inputs(value):
    if isinstance(value, list) and len(value) == 1:
        kind, oracle_message = _condition("multiple_fewer_than_two")
        assert kind == "error" and "at least 2" in oracle_message
    with pytest.raises((TypeError, ValueError), match="at least 2|STSList"):
        _tentrop_m(value)


def test_list_and_tuple_defaults_use_seq_number_labels():
    x = _fixture(weighted=False)
    for values in ([x, x.copy()], (x, x.copy())):
        fig, ax = _tentrop_m(values)
        assert [line.get_label() for line in ax.lines] == ["seq1", "seq2"]
        plt.close(fig)


def test_custom_labels_colors_types_widths_cycle_and_style_forwards():
    fig, ax = _tentrop(
        _fixture(),
        _group(),
        main="Custom",
        col=["red"],
        lty=["dashed", "dotted"],
        lwd=[1.0, 2.0],
        glabels=["G", "B", "A"],
        alpha=0.4,
    )
    assert ax.get_title() == "Custom"
    assert [line.get_color() for line in ax.lines] == ["red"] * 3
    assert [line.get_linestyle() for line in ax.lines] == ["--", ":", "--"]
    assert [line.get_linewidth() for line in ax.lines] == [1.0, 2.0, 1.0]
    assert [line.get_alpha() for line in ax.lines] == [0.4] * 3
    assert [text.get_text() for text in ax.get_legend().get_texts()] == ["G", "B", "A"]
    plt.close(fig)


@pytest.mark.parametrize(("keyword", "value"), [("col", []), ("lty", []), ("lwd", [])])
def test_empty_style_vectors_are_rejected(keyword, value):
    with pytest.raises(ValueError, match=keyword):
        _tentrop(_fixture(), _group(), **{keyword: value})


def test_bad_explicit_label_count_is_rejected():
    with pytest.raises(ValueError, match="glabels"):
        _tentrop(_fixture(), _group(), glabels=["only one"])


def test_axis_overrides_hidden_legend_none_title_horizontal_legend_and_supplied_axes():
    fig, supplied = plt.subplots()
    returned_fig, ax = _tentrop(
        _fixture(),
        _group(),
        main=None,
        ylim=(-0.2, 1.2),
        xtlab=["A", "B", "C", "D", "E"],
        xtstep=3,
        tick_last=True,
        with_legend=False,
        horiz=True,
        cex_legend=1.2,
        ax=supplied,
    )
    assert returned_fig is fig and ax is supplied
    assert ax.get_title() == ""
    assert ax.get_ylim() == (-0.2, 1.2)
    np.testing.assert_array_equal(ax.get_xticks(), [1, 4, 5])
    assert [tick.get_text() for tick in ax.get_xticklabels()] == ["A", "D", "E"]
    assert ax.get_legend() is None
    plt.close(fig)


def test_horizontal_legend_uses_one_row():
    fig, ax = _tentrop(_fixture(), _group(), horiz=True, legend_pos="lower right")
    assert ax.get_legend()._ncols == 3
    plt.close(fig)


def test_empty_and_all_missing_analytics_render_without_crashing():
    empty = STSList(data=np.empty((0, 3), dtype=int), alphabet=["A"])
    all_missing = STSList(data=np.ones((2, 3), dtype=int), alphabet=["A", "*"], nr_code=1)
    for values in ([empty, empty.copy()], [all_missing, all_missing.copy()]):
        fig, ax = _tentrop_m(values, with_legend=False)
        assert len(ax.lines) == 2
        assert all(np.isnan(line.get_ydata()).all() for line in ax.lines)
        plt.close(fig)


def test_callers_are_not_mutated_and_public_routes_are_identity_exports():
    x = _fixture()
    data = x.data.copy()
    weights = x.weights.copy()
    group = _group()
    original_group = group.copy()
    fig, _ = _tentrop(x, group)
    plt.close(fig)

    np.testing.assert_array_equal(x.data, data)
    np.testing.assert_array_equal(x.weights, weights)
    np.testing.assert_array_equal(group.codes, original_group.codes)
    assert list(group.categories) == list(original_group.categories)
    assert getattr(traminex, "seqplot_tentrop") is getattr(summary, "seqplot_tentrop")
    assert getattr(traminex, "seqplot_tentrop_m") is getattr(summary, "seqplot_tentrop_m")
