import numpy as np

from traminex import seqdef
from traminex.distances import pairwise
from traminex.distances.costs import trate_substitution_matrix


def test_om_trate_convenience_matches_explicit():
    mat = np.array(
        [
            [0, 1, 2, 2],
            [0, 1, 1, 2],
            [0, 2, 2, 2],
        ],
        dtype=int,
    )
    x = seqdef(mat)
    D1 = pairwise(x.data, method="OM", subm="TRATE", trate_smoothing=1e-6, trate_normalize="none")
    sm = trate_substitution_matrix(x.data, smoothing=1e-6, mode="symmetric", normalize="none")
    D2 = pairwise(x.data, method="OM", subm=sm)
    assert np.allclose(D1, D2)
