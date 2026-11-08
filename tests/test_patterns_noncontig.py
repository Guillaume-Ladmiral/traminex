import numpy as np

from traminex import seqdef
from traminex.patterns import frequent_subsequences


def test_frequent_subsequences_noncontiguous():
    mat = np.array(
        [
            [0, 1, 2, 2],
            [0, 1, 1, 2],
        ],
        dtype=int,
    )
    x = seqdef(mat)
    cont = frequent_subsequences(x, max_len=2, contiguous=True, min_support=1)
    noncont = frequent_subsequences(x, max_len=2, contiguous=False, min_support=1)
    assert len(noncont) >= len(cont)
