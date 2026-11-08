"""chronogram(seaborn=True) puts each time label at the centre of its column.

What was wrong (2026-09-22 audit S6-21): the seaborn heatmap draws time
position ``j`` as the cell ``[j, j + 1]``, centred at ``j + 0.5``, but the time
ticks were placed at ``0 .. L-1`` -- on the cell edges -- so every label sat half
a cell to the left of its column. The matplotlib path (``imshow``, cells
centred at ``j``) was right.

Source of truth: seaborn's heatmap geometry (the state ticks on the same axes
are already placed at ``k + 0.5``) and the matplotlib path's labels.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from traminex import seqdef
from traminex.plotting import chronogram

ROWS = [
    ["A", "A", "B", "B", "C", "C", "C", "C", "B", "A", "A", "B"],
    ["A", "B", "B", "C", "C", "C", "A", "A", "B", "B", "C", "C"],
    ["C", "C", "B", "A", "A", "A", "B", "B", "C", "A", "B", "C"],
]


def _ticks(ax):
    return np.asarray(ax.get_xticks(), dtype=float), [t.get_text() for t in ax.get_xticklabels()]


def test_seaborn_time_ticks_sit_on_the_cell_centres():
    pytest.importorskip("seaborn")
    x = seqdef(ROWS, alphabet=["A", "B", "C"])

    fig, ax = chronogram(x, seaborn=True)
    positions, labels = _ticks(ax)
    plt.close(fig)
    fig, ax = chronogram(x, seaborn=False)
    reference_positions, reference_labels = _ticks(ax)
    plt.close(fig)

    np.testing.assert_array_equal(positions, reference_positions + 0.5)
    assert labels == reference_labels
