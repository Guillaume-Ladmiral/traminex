import numpy as np
import pytest

from traminex.distances.costs import trate_substitution_matrix_timevarying


def test_timevarying_trate_modes_normalizations_and_weights():
    seqmat = np.array(
        [[0, 1, 2, 2], [0, 1, 1, 2], [2, 2, 1, 0]],
        dtype=np.int32,
    )
    weights = np.array([1.0, 3.0, 2.0])

    row = trate_substitution_matrix_timevarying(
        seqmat,
        smoothing=0.25,
        mode="row",
        normalize="mean_offdiag",
        weights=weights,
    )
    col = trate_substitution_matrix_timevarying(
        seqmat,
        smoothing=0.25,
        mode="col",
        normalize="mean_offdiag",
        weights=weights,
    )
    clipped = trate_substitution_matrix_timevarying(
        seqmat,
        mode="symmetric",
        normalize="clip_0_2",
    )

    assert row.shape == col.shape == clipped.shape == (4, 3, 3)
    np.testing.assert_allclose(col, row.transpose(0, 2, 1), rtol=0, atol=1e-15)
    assert np.all((clipped >= 0) & (clipped <= 2))
    np.testing.assert_array_equal(np.diagonal(clipped, axis1=1, axis2=2), 0)


def test_timevarying_trate_empty_and_refusal_boundaries():
    empty = np.empty((2, 0), dtype=np.int32)
    np.testing.assert_array_equal(
        trate_substitution_matrix_timevarying(empty),
        np.zeros((0, 0, 0)),
    )

    seqmat = np.array([[0, 1], [1, 0]], dtype=np.int32)
    with pytest.raises(ValueError, match="2D"):
        trate_substitution_matrix_timevarying(seqmat[0])
    with pytest.raises(ValueError, match="lag must be 1"):
        trate_substitution_matrix_timevarying(seqmat, lag=2)
    with pytest.raises(ValueError, match="Unsupported TRATE mode"):
        trate_substitution_matrix_timevarying(seqmat, mode="invalid")
    with pytest.raises(ValueError, match="Unsupported normalize"):
        trate_substitution_matrix_timevarying(seqmat, normalize="invalid")
