import numpy as np

from traminex.distances.costs import indel_from_subm, trate_substitution_matrix
from traminex.distances.python_kernels import pairwise as py_pairwise


def test_trate_and_indel_utils():
    mat = np.array(
        [
            [0, 1, 1, 2],
            [1, 0, 1, 2],
            [1, 1, 2, 2],
        ],
        dtype=int,
    )
    subm = trate_substitution_matrix(mat)
    K = int(mat.max()) + 1
    assert subm.shape == (K, K)
    assert np.all(np.diag(subm) == 0.0)
    indel = indel_from_subm(subm)
    assert indel > 0


def test_dhd_distance_basic():
    mat = np.array(
        [
            [0, 1, 2, 2],
            [0, 1, 1, 2],
        ],
        dtype=int,
    )
    D = py_pairwise(mat, method="DHD")
    assert D.shape == (2, 2)
    assert D[0, 1] >= 0.0
