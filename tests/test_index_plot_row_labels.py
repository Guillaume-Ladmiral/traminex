"""Index plots label a bounded number of rows (2026-09-22 audit S6-16).

What was wrong: ``sequence_index_plot`` gave every displayed row its own y-axis
tick label. Matplotlib lays out every label, so 5,000 rows took about 32 s to
render and save against 0.14 s without the labels (audit ``probe14b``), and
the labels were unreadably overlapped anyway.

Source of truth: R's base ``axis()`` in ``plot.stslist.R`` is given every row
label but draws only those that do not overlap. traminex now labels at most 50
evenly spaced rows, always including the first and last, with each label still
naming its own row; plots of up to 50 rows keep every label.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from traminex.plotting import sequence_index_plot
from traminex.sequences import STSList


def _sequences(n: int) -> STSList:
    data = np.random.default_rng(3).integers(0, 3, size=(n, 5))
    return STSList(data, alphabet=["A", "B", "C"], ids=[f"id{i}" for i in range(n)])


def test_many_rows_get_at_most_fifty_labels_naming_their_rows():
    x = _sequences(500)

    fig, ax = sequence_index_plot(x)

    ticks = [int(tick) for tick in ax.get_yticks()]
    labels = [label.get_text() for label in ax.get_yticklabels()]
    assert len(ticks) <= 50
    assert ticks[0] == 0 and ticks[-1] == 499
    assert labels == [f"id{row}" for row in ticks]
    plt.close(fig)


def test_few_rows_keep_every_label():
    x = _sequences(12)

    fig, ax = sequence_index_plot(x)

    assert [label.get_text() for label in ax.get_yticklabels()] == [f"id{i}" for i in range(12)]
    plt.close(fig)
