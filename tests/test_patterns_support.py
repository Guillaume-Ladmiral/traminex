import numpy as np

from traminex import seqdef
from traminex.patterns import frequent_subsequences


def test_frequent_subsequences_sequences_support():
    mat = np.array(
        [
            [0, 1, 2, 2],
            [0, 1, 1, 2],
            [0, 1, 2, 1],
        ],
        dtype=int,
    )
    x = seqdef(mat)
    items_seq = frequent_subsequences(
        x, max_len=4, contiguous=False, min_support=1, support_type="sequences"
    )
    assert isinstance(items_seq, list)
    assert len(items_seq) > 0
