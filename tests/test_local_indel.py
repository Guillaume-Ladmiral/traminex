import numpy as np
import pytest

from traminex import seqdef
from traminex.distances.costs import local_indel_by_pos
from traminex.distances.python_kernels import pairwise as py_pairwise


def test_local_indel_by_pos_reduces_cost_on_change_positions():
    # High change rate in the middle positions; local indel lowers cost there
    mat = np.array(
        [
            [0, 1, 2, 2, 2],
            [0, 2, 1, 2, 2],
            [0, 2, 2, 1, 2],
            [0, 1, 2, 1, 2],
        ],
        dtype=int,
    )
    x = seqdef(mat)
    D_base = py_pairwise(x.data, method="OM", indel=1.0, subm=2.0)
    indel_by_pos = local_indel_by_pos(x.data, base=1.0, scale=0.8, smooth_window=3)
    D_loc = py_pairwise(x.data, method="OM", indel=1.0, subm=2.0, indel_by_pos=indel_by_pos)
    # On average, local indel should not increase the distance matrix dramatically
    assert D_loc.mean() <= D_base.mean() + 1e-6


def test_even_smoothing_window_preserves_the_sequence_width():
    seqmat = np.array([[0, 1, 1, 2], [0, 0, 1, 2], [2, 1, 0, 0]], dtype=int)

    profile = local_indel_by_pos(seqmat, smooth_window=2)

    assert profile.shape == (seqmat.shape[1],)
    assert np.isfinite(profile).all()


@pytest.mark.parametrize("smooth_window", [True, 1.5])
def test_local_indel_rejects_non_integer_smoothing_windows(smooth_window):
    seqmat = np.array([[0, 1, 0], [1, 1, 0]], dtype=int)

    with pytest.raises(TypeError, match="smooth_window must be an integer"):
        local_indel_by_pos(seqmat, smooth_window=smooth_window)
