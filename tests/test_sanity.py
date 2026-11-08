import numpy as np

from traminex import STSList, seqdef
from traminex.distances.python_kernels import pairwise
from traminex.stats import seqistatd


def test_basic_seqdef_and_stats():
    mat = np.array([["A", "B", "C"], ["A", "B", "B"], ["C", "B", "A"]])
    x = seqdef(mat)
    assert isinstance(x, STSList)
    assert x.nseq == 3 and x.seqlen == 3
    st = seqistatd(x)
    assert st.shape == (3, len(x.alphabet))


def test_pairwise_om():
    mat = np.array([[0, 1, 2], [0, 1, 1]], dtype=int)
    D = pairwise(mat, method="OM", indel=1.0, subm=2.0)
    assert D.shape == (2, 2)
    assert D[0, 1] > 0
