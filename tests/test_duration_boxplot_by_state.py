"""duration_boxplot honours ``by_state`` on its matplotlib (no seaborn) path.

What was wrong (2026-09-22 audit S6-04): with ``by_state=True`` and
``seaborn=False`` the per-state split and ``labels`` were silently ignored; the
plot showed one box pooling every (sequence, state) cell, including the zero
cells of the missing and void marker columns.

Source of truth: the function's own contract -- the seaborn path draws one
distribution of per-sequence time spent per state (``seqistatd``), split by
cluster when ``labels`` is given. The matplotlib path must draw the same
groups. Structural void is never a substantive state
(``docs/VISUALIZATION_SERIES_CONTRACT.md``), and ``seqistatd`` is called
without missing values, so the marker columns are not states to plot.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from traminex import seqdef, seqistatd
from traminex.plotting import duration_boxplot

ROWS = [
    ["A", "A", "B", "B", "C", "C"],
    ["A", "B", "B", "C", "C", "C"],
    ["B", "B", None, "C", "D", "D"],
    ["A", "A", "A", "A", "B", "B"],
    ["C", "C", "C", "D", "D", None],
    ["A", None, "B", "B", "C", "C"],
    ["D", "D", "D", "D", "D", "D"],
    ["A", "A", "B", "B", "C", "C"],
]
CLUSTERS = np.array([1, 1, 2, 2, 1, 2, 1, 2])


def _sequences():
    return seqdef(ROWS, alphabet=["A", "B", "C", "D"])


def _medians(ax):
    # Matplotlib boxplot draws, per box, 2 whiskers, 2 caps, 1 median and
    # 1 flier line; the median is the fifth line of each group of six.
    return [float(np.mean(line.get_ydata())) for line in ax.lines[4::6]]


def test_by_state_without_seaborn_draws_one_box_per_state():
    x = _sequences()
    durations = seqistatd(x, include_missing=False)

    fig, ax = duration_boxplot(x, by_state=True, seaborn=False)

    assert [tick.get_text() for tick in ax.get_xticklabels()] == ["A", "B", "C", "D"]
    assert len(ax.patches) == 4
    assert _medians(ax) == pytest.approx([np.median(durations[:, k]) for k in range(4)])
    assert ax.get_title() == "Durations by state"
    plt.close(fig)


def test_by_state_with_labels_without_seaborn_splits_every_state_by_cluster():
    x = _sequences()
    durations = seqistatd(x, include_missing=False)

    fig, ax = duration_boxplot(x, CLUSTERS, by_state=True, seaborn=False)

    expected_ticks = [f"{state}\n{cluster}" for state in "ABCD" for cluster in (1, 2)]
    assert [tick.get_text() for tick in ax.get_xticklabels()] == expected_ticks
    assert len(ax.patches) == 8
    expected = [
        np.median(durations[CLUSTERS == cluster, k]) for k in range(4) for cluster in (1, 2)
    ]
    assert _medians(ax) == pytest.approx(expected)
    assert ax.get_title() == "Durations by state per cluster"
    plt.close(fig)


def test_seaborn_by_state_draws_the_same_states():
    pytest.importorskip("seaborn")
    x = _sequences()

    fig, ax = duration_boxplot(x, by_state=True, seaborn=True)

    assert [tick.get_text() for tick in ax.get_yticklabels()] == ["A", "B", "C", "D"]
    plt.close(fig)


def test_total_duration_path_is_unchanged():
    x = _sequences()
    durations = seqistatd(x, include_missing=False)

    fig, ax = duration_boxplot(x, seaborn=False)

    assert len(ax.patches) == 1
    assert _medians(ax) == pytest.approx([np.median(durations.sum(axis=1))])
    plt.close(fig)
