import warnings

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib import MatplotlibDeprecationWarning

from traminex import seqdef
from traminex.plotting import duration_boxplot, seqfrequency_plot


def _duration_fixture():
    matrix = np.array([[0, 1, 2, 2], [0, 1, 1, 2], [2, 2, 1, 0]], dtype=int)
    return seqdef(matrix)


def test_seqfrequency_and_duration_plots_smoke():
    x = _duration_fixture()
    fig_frequency, ax = seqfrequency_plot(x, top=2)
    assert ax is not None
    fig_duration, ax = duration_boxplot(x, labels=np.array([0, 0, 1]), by_state=False)
    assert ax is not None
    plt.close(fig_frequency)
    plt.close(fig_duration)


def test_duration_boxplot_has_no_matplotlib_deprecation_warning():
    x = _duration_fixture()

    with warnings.catch_warnings():
        warnings.simplefilter("error", MatplotlibDeprecationWarning)
        fig, ax = duration_boxplot(x, labels=np.array([0, 0, 1]), by_state=False)

    assert [tick.get_text() for tick in ax.get_xticklabels()] == ["0", "1"]
    assert len(ax.patches) == 2
    plt.close(fig)


def test_duration_boxplot_preserves_string_group_tick_labels():
    labels = np.array(["control", "control", "treated"])
    labels_before = labels.copy()

    fig, ax = duration_boxplot(_duration_fixture(), labels=labels, by_state=False, seaborn=False)

    assert [tick.get_text() for tick in ax.get_xticklabels()] == ["control", "treated"]
    np.testing.assert_array_equal(labels, labels_before)
    plt.close(fig)


def test_duration_boxplot_preserves_categorical_like_legend_labels():
    labels = np.array(["control", "control", "treated"], dtype=np.str_)
    labels_before = labels.copy()

    fig, ax = duration_boxplot(_duration_fixture(), labels=labels, by_state=True, seaborn=True)

    assert ax.get_legend() is not None
    assert [text.get_text() for text in ax.get_legend().get_texts()] == ["control", "treated"]
    np.testing.assert_array_equal(labels, labels_before)
    plt.close(fig)


@pytest.mark.parametrize("by_state,seaborn", [(False, False), (False, True), (True, False)])
def test_duration_boxplot_preserves_integer_compatibility_and_inputs(by_state, seaborn):
    labels = np.array([0, 0, 1], dtype=np.int64)
    labels_before = labels.copy()

    fig, ax = duration_boxplot(
        _duration_fixture(), labels=labels, by_state=by_state, seaborn=seaborn
    )

    np.testing.assert_array_equal(labels, labels_before)
    if not by_state:
        assert [tick.get_text() for tick in ax.get_xticklabels()] == ["0", "1"]
    plt.close(fig)


@pytest.mark.parametrize("by_state", [False, True])
@pytest.mark.parametrize("seaborn", [False, True])
@pytest.mark.parametrize(
    "labels",
    [
        np.array(["control", "treated"]),
        np.array(["control", "treated", "treated", "extra"]),
        np.array("control"),
        np.array([["control", "treated", "treated"]]),
    ],
    ids=["short", "long", "scalar", "two-dimensional"],
)
def test_duration_boxplot_rejects_invalid_label_shape(labels, by_state, seaborn):
    with pytest.raises(
        ValueError,
        match="labels must be one-dimensional with exactly one label per sequence",
    ):
        duration_boxplot(_duration_fixture(), labels=labels, by_state=by_state, seaborn=seaborn)
