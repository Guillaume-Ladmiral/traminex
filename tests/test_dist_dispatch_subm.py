import numpy as np

from traminex.distances import pairwise as pairwise_auto


def test_dispatch_with_subm_matrix_fallback():
    # simple two sequences, 3 states
    mat = np.array([[0, 1, 2], [2, 1, 0]], dtype=int)
    # subm: 0 on diag, 2 off-diag
    subm = np.full((3, 3), 2.0)
    np.fill_diagonal(subm, 0.0)
    D = pairwise_auto(mat, method="OM", indel=1.0, subm=subm)
    assert D.shape == (2, 2)
    assert D[0, 1] >= 0.0
