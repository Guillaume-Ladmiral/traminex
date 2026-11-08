import numpy as np

from traminex import seqdef
from traminex.events import frequent_ngrams


def test_frequent_ngrams():
    mat = np.array(
        [
            [0, 1, 2, 2],
            [0, 1, 1, 2],
        ],
        dtype=int,
    )
    x = seqdef(mat)
    items = frequent_ngrams(x, n=2, min_support=1)
    # Should include ('0','1') mapped to alphabet strings, but alphabet inferred from ints as strings
    assert isinstance(items, list)
    assert len(items) > 0
