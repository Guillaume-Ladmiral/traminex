import matplotlib
import numpy as np

matplotlib.use("Agg")
from traminex import cluster_chronograms, overlay_representatives_on_index_plot, seqdef
from traminex.distances import pairwise


def test_grouped_and_overlay_plots_smoke():
    mat = np.array(
        [
            [0, 1, 2, 2],
            [0, 1, 1, 2],
            [2, 2, 1, 0],
            [2, 2, 2, 0],
        ],
        dtype=int,
    )
    x = seqdef(mat)
    labels = np.array([0, 0, 1, 1])
    D = pairwise(x.data, method="HAM")
    fig, ax = overlay_representatives_on_index_plot(x, labels, D)
    assert ax is not None
    fig, axes = cluster_chronograms(x, labels, include_missing=False, cols=2)
    assert axes is not None
