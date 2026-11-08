import numpy as np

from traminex import (
    chronogram,
    seqdef,
    sequence_index_plot,
    state_distribution_plot,
    transition_heatmap,
)


def test_plotting_smoke(tmp_path):
    mat = np.array([[0, 1, 2, 2], [0, 1, 1, 2], [2, 2, 1, 0]], dtype=int)
    x = seqdef(mat)
    # Basic plots; we only test that they run and return axes
    fig, ax = sequence_index_plot(x, highlight=[1])
    assert ax is not None
    fig, ax = chronogram(x)
    assert ax is not None
    fig, ax = state_distribution_plot(x)
    assert ax is not None
    fig, ax = transition_heatmap(x, normalize=True)
    assert ax is not None
