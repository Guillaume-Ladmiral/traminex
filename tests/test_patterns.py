import numpy as np

from traminex import seqdef
from traminex.patterns import frequent_subsequences


def test_frequent_subsequences_contiguous():
    mat = np.array(
        [
            [0, 1, 2, 2],
            [0, 1, 1, 2],
            [0, 1, 2, 1],
        ],
        dtype=int,
    )
    x = seqdef(mat)
    items = frequent_subsequences(x, max_len=3, contiguous=True, min_support=2)
    assert isinstance(items, list)
    assert len(items) > 0
