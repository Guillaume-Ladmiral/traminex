import numpy as np

from traminex import seqdef
from traminex.distances.costs import local_indel_by_pos


def test_local_indel_entropy_mode():
    mat = np.array(
        [
            [0, 1, 2, 2],
            [0, 1, 1, 2],
            [0, 1, 2, 1],
        ],
        dtype=int,
    )
    x = seqdef(mat)
    indel_entropy = local_indel_by_pos(
        x.data, method="entropy", base=1.0, scale=0.5, smooth_window=3
    )
    assert indel_entropy.shape[0] == x.seqlen
    assert np.all(indel_entropy > 0)
