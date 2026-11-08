import matplotlib
import numpy as np

matplotlib.use("Agg")
from traminex import cluster_chronograms_with_representatives, seqdef
from traminex.distances import pairwise


def test_cluster_chronograms_with_representatives_smoke():
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
    fig, axes = cluster_chronograms_with_representatives(
        x, labels, D, include_missing=False, cols=2
    )
    assert axes is not None
