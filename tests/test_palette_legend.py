import numpy as np

from traminex import seqdef
from traminex.plotting import sequence_index_plot, state_legend


def test_state_legend_smoke():
    mat = np.array([[0, 1, 2], [2, 1, 0]], dtype=int)
    x = seqdef(mat)
    # assign a simple palette
    x.cpal = ["#1f77b4", "#ff7f0e", "#2ca02c"]
    fig, ax = sequence_index_plot(x)
    patches, labels = state_legend(x, ax=ax)
    assert len(patches) == len(x.alphabet)
