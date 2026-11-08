import numpy as np

from traminex import seqdef
from traminex.cluster import cluster_agglomerative, cluster_representatives
from traminex.distances import pairwise


def test_cluster_representatives_smoke():
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
    D = pairwise(x.data, method="HAM")
    res = cluster_agglomerative(D, n_clusters=2)
    reps = cluster_representatives(x, res.labels_, D)
    assert set(["cluster", "rep_index", "rep_sequence"]).issubset(reps.columns)
