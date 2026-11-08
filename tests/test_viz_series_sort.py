"""RP7a foundation: the plot-series sort resolver + sequence_index_plot sortv.

resolve_sort_order mirrors TraMineR plot.stslist's ``sortv`` block
(SOURCES/TraMineR/R/plot.stslist.R): "from.start" = lexicographic with the first
position primary, "from.end" = last position primary, an array = stable order.
Structure-tested: the index plot's imshow rows == x.data reordered by the same
permutation.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import numpy as np
from corestats_fixtures import fixture_a

from traminex.visualization.plots.index import sequence_index_plot
from traminex.visualization.series import resolve_sort_order


def test_resolve_sort_order_none_is_identity():
    x = fixture_a()
    assert list(resolve_sort_order(x, None)) == [0, 1, 2, 3]


def test_resolve_sort_order_from_start_primary_is_first_column():
    x = fixture_a()
    # data first column: [0(A), 0(A), 3(D), 0(A)] -> the D-row (idx 2) sorts last
    order = resolve_sort_order(x, "from.start")
    # the three A-starting sequences come first (stable among equal first cols by
    # subsequent columns), the D-starting one last
    assert order[-1] == 2
    # verify it is a valid lexicographic order keyed on the whole row
    rows = x.data[order]
    assert rows[-1, 0] == 3  # 'D'


def test_resolve_sort_order_from_end_primary_is_last_column():
    x = fixture_a()
    order = resolve_sort_order(x, "from.end")
    # last columns: seq0=C(2), seq1=A(0), seq2=D(3), seq3=B(1)
    # ascending by last col -> [seq1(A), seq3(B), seq0(C), seq2(D)]
    assert list(order) == [1, 3, 0, 2]


def test_resolve_sort_order_numeric_is_stable_argsort():
    x = fixture_a()
    sortv = np.array([2.0, 0.0, 2.0, 1.0])
    # stable argsort: 0.0(idx1), 1.0(idx3), 2.0(idx0), 2.0(idx2)
    assert list(resolve_sort_order(x, sortv)) == [1, 3, 0, 2]


def test_index_plot_sortv_reorders_imshow_rows():
    x = fixture_a()
    order = resolve_sort_order(x, "from.end")
    fig, ax = sequence_index_plot(x, sortv="from.end")
    arr = np.asarray(ax.images[0].get_array())
    assert np.array_equal(arr, x.data[order])
