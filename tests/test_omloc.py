"""Tests for OMloc (localized OM with position-dependent costs)."""

import numpy as np
import pytest

from traminex.distances import pairwise
from traminex.distances.costs import (
    local_indel_by_pos,
    trate_substitution_matrix,
)
from traminex.distances.omloc import omloc


@pytest.fixture
def seqmat():
    return np.array(
        [[0, 0, 1, 1, 2], [0, 1, 1, 2, 2], [2, 2, 1, 0, 0]],
        dtype=np.int32,
    )


def test_omloc_basic():
    """Test basic OMloc computation with default settings."""
    seqmat = np.array([[0, 0, 1, 1, 2], [0, 1, 1, 2, 2], [2, 2, 1, 0, 0]], dtype=np.int32)

    D = omloc(seqmat)

    assert D.shape == (3, 3)
    assert np.allclose(D, D.T), "Distance matrix should be symmetric"
    assert np.allclose(np.diag(D), 0), "Diagonal should be zero"
    assert np.all(D >= 0), "Distances should be non-negative"


def test_omloc_timevarying_subm():
    """Test OMloc with time-varying substitution matrix."""
    seqmat = np.array([[0, 0, 1, 1, 2], [0, 1, 1, 2, 2], [2, 2, 1, 0, 0]], dtype=np.int32)

    D_global = omloc(seqmat, subm_timevarying=False)
    D_timevar = omloc(seqmat, subm_timevarying=True)

    assert D_global.shape == D_timevar.shape
    # Time-varying should produce different results
    # (though they might be similar for small data)
    assert D_global.shape == (3, 3)
    assert D_timevar.shape == (3, 3)


def test_omloc_indel_methods():
    """Test different indel localization methods."""
    seqmat = np.array([[0, 0, 1, 1, 2], [0, 1, 1, 2, 2], [2, 2, 1, 0, 0]], dtype=np.int32)

    D_change = omloc(seqmat, indel_method="change")
    D_entropy = omloc(seqmat, indel_method="entropy")

    assert D_change.shape == (3, 3)
    assert D_entropy.shape == (3, 3)
    # Different methods should generally produce different results
    # but both should be valid


def test_omloc_subm_constant():
    """Test OMloc with constant substitution matrix."""
    seqmat = np.array([[0, 0, 1, 1, 2], [0, 1, 1, 2, 2], [2, 2, 1, 0, 0]], dtype=np.int32)

    D = omloc(seqmat, subm_base="CONSTANT")

    assert D.shape == (3, 3)
    assert np.all(D >= 0)


def test_omloc_with_scalar_subm():
    """Test OMloc with scalar substitution cost."""
    seqmat = np.array([[0, 0, 1, 1, 2], [0, 1, 1, 2, 2], [2, 2, 1, 0, 0]], dtype=np.int32)

    D = omloc(seqmat, subm_base=2.5)

    assert D.shape == (3, 3)
    assert np.all(D >= 0)


def test_omloc_with_custom_matrix():
    """Test OMloc with user-provided substitution matrix."""
    seqmat = np.array([[0, 0, 1, 1, 2], [0, 1, 1, 2, 2], [2, 2, 1, 0, 0]], dtype=np.int32)

    K = 3  # number of states
    custom_subm = np.ones((K, K), dtype=float) * 1.5
    np.fill_diagonal(custom_subm, 0.0)

    D = omloc(seqmat, subm_base=custom_subm)

    assert D.shape == (3, 3)
    assert np.all(D >= 0)


def test_omloc_normalization():
    """Test OMloc with different normalization options."""
    seqmat = np.array([[0, 0, 1, 1, 2], [0, 1, 1, 2, 2], [2, 2, 1, 0, 0]], dtype=np.int32)

    D_none = omloc(seqmat, norm=None)
    D_length = omloc(seqmat, norm="length")
    D_max = omloc(seqmat, norm="max")

    assert D_none.shape == D_length.shape == D_max.shape
    # Normalized distances should generally be smaller
    assert np.mean(D_length) <= np.mean(D_none) or np.mean(D_none) == 0


def test_omloc_indel_smoothing():
    """Test OMloc with indel cost smoothing."""
    seqmat = np.array([[0, 0, 1, 1, 2], [0, 1, 1, 2, 2], [2, 2, 1, 0, 0]], dtype=np.int32)

    D_nosmooth = omloc(seqmat, indel_smooth=None)
    D_smooth = omloc(seqmat, indel_smooth=3)

    assert D_nosmooth.shape == D_smooth.shape
    # Smoothing may affect results but both should be valid


def test_omloc_larger_dataset():
    """Test OMloc on a slightly larger dataset."""
    np.random.seed(42)
    seqmat = np.random.randint(0, 4, size=(10, 20), dtype=np.int32)

    D = omloc(seqmat)

    assert D.shape == (10, 10)
    assert np.allclose(D, D.T)
    assert np.allclose(np.diag(D), 0)
    assert np.all(D >= 0)


def test_omloc_condensed_output():
    """Test OMloc with condensed output format."""
    seqmat = np.array([[0, 0, 1, 1, 2], [0, 1, 1, 2, 2], [2, 2, 1, 0, 0]], dtype=np.int32)

    D_full = omloc(seqmat, condensed=False)
    D_cond = omloc(seqmat, condensed=True)

    assert D_full.shape == (3, 3)
    assert D_cond.shape == (3,)  # condensed upper triangle for n=3: 3*(3-1)/2 = 3
    # Verify condensed values match upper triangle
    expected = D_full[np.triu_indices(3, k=1)]
    assert np.allclose(D_cond, expected)


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(r"ignore:\s*Non\-symmetric\ substitution\ cost\ matrix!:UserWarning")
def test_omloc_weighted_global_trate_matches_explicit_regular_om(seqmat):
    """The wrapper consumes TRATE weights without leaking them to regular OM."""
    weights = np.array([1.0, 3.0, 2.0])
    indels = local_indel_by_pos(seqmat, method="entropy", base=1.2, scale=0.7)
    subm = trate_substitution_matrix(
        seqmat,
        smoothing=0.25,
        mode="row",
        normalize="mean_offdiag",
        weights=weights,
    )
    expected = pairwise(
        seqmat,
        method="OM",
        indel=float(indels.mean()),
        indel_by_pos=indels,
        subm=subm,
        norm=None,
    )

    actual = omloc(
        seqmat,
        indel_method="entropy",
        indel_base=1.2,
        indel_scale=0.7,
        subm_base="TRATE",
        trate_smoothing=0.25,
        trate_mode="row",
        trate_normalize="mean_offdiag",
        trate_weights=weights,
        norm=None,
    )

    np.testing.assert_allclose(actual, expected, rtol=0, atol=0)


@pytest.mark.parametrize(
    ("subm_base", "cost"),
    [
        ("CONSTANT", 2.0),
        (1.5, 1.5),
        (
            np.array(
                [[0.0, 1.5, 1.5], [1.5, 0.0, 1.5], [1.5, 1.5, 0.0]],
                dtype=float,
            ),
            1.5,
        ),
    ],
)
def test_omloc_timevarying_constant_paths_match_explicit_regular_om(seqmat, subm_base, cost):
    """String, scalar, and 2-D time-varying paths construct the same slices."""
    length = seqmat.shape[1]
    states = int(seqmat.max()) + 1
    subm_by_pos = np.full((length, states, states), cost, dtype=float)
    for matrix in subm_by_pos:
        np.fill_diagonal(matrix, 0.0)
    indels = local_indel_by_pos(seqmat, method="change", base=1.0, scale=0.8)
    expected = pairwise(
        seqmat,
        method="OM",
        indel=float(indels.mean()),
        indel_by_pos=indels,
        subm_by_pos=subm_by_pos,
        norm=None,
    )

    actual = omloc(
        seqmat,
        indel_scale=0.8,
        subm_timevarying=True,
        subm_base=subm_base,
        norm=None,
    )

    np.testing.assert_allclose(actual, expected, rtol=0, atol=0)


def test_omloc_timevarying_tensor_is_forwarded_without_reconstruction(seqmat):
    """A correctly shaped 3-D tensor is the explicit lower-level construction."""
    length = seqmat.shape[1]
    states = int(seqmat.max()) + 1
    subm_by_pos = np.empty((length, states, states), dtype=float)
    for position, matrix in enumerate(subm_by_pos):
        matrix.fill(1.0 + position / 4)
        np.fill_diagonal(matrix, 0.0)
    indels = local_indel_by_pos(seqmat)
    expected = pairwise(
        seqmat,
        method="OM",
        indel=float(indels.mean()),
        indel_by_pos=indels,
        subm_by_pos=subm_by_pos,
        norm="length",
    )

    actual = omloc(
        seqmat,
        subm_timevarying=True,
        subm_base=subm_by_pos,
        norm="length",
    )

    np.testing.assert_allclose(actual, expected, rtol=0, atol=0)


@pytest.mark.parametrize(
    ("bad", "match"),
    [
        (np.array([0, 1, 2], dtype=np.int32), "2D"),
        (np.array([[0.0, 1.0]], dtype=float), "integer"),
        (np.array([[False, True]], dtype=bool), "integer"),
        (np.array([[0, -1]], dtype=np.int32), "non-negative"),
        (np.empty((0, 3), dtype=np.int32), "at least one sequence"),
        (np.empty((2, 0), dtype=np.int32), "positive length"),
    ],
)
def test_omloc_rejects_invalid_sequence_matrices(bad, match):
    with pytest.raises(ValueError, match=match):
        omloc(bad)


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"indel_base": 0.0}, "indel_base.*positive"),
        ({"indel_base": np.inf}, "indel_base.*finite"),
        ({"indel_scale": -0.1}, "indel_scale.*non-negative"),
        ({"indel_scale": np.nan}, "indel_scale.*finite"),
        ({"indel_smooth": 0}, "indel_smooth.*positive integer"),
        ({"indel_smooth": 1.5}, "indel_smooth.*positive integer"),
        ({"subm_timevarying": 1}, "subm_timevarying.*bool"),
        ({"trate_smoothing": -0.1}, "trate_smoothing.*non-negative"),
        ({"trate_smoothing": np.inf}, "trate_smoothing.*finite"),
    ],
)
def test_omloc_rejects_invalid_parameter_domains(seqmat, kwargs, match):
    with pytest.raises((TypeError, ValueError), match=match):
        omloc(seqmat, **kwargs)


@pytest.mark.parametrize(
    ("subm_timevarying", "subm_base", "match"),
    [
        (False, "unknown", "Unknown subm_base"),
        (True, "unknown", "Unknown subm_base"),
        (False, np.ones(3), "2D"),
        (True, np.ones(3), "2D.*3D"),
        (False, np.zeros((3, 3, 3)), "subm_timevarying=True"),
        (False, np.zeros((2, 2)), "shape"),
        (True, np.zeros((4, 3, 3)), "shape"),
        (True, np.zeros((5, 2, 2)), "shape"),
        (False, np.full((3, 3), np.nan), "finite"),
        (False, np.full((3, 3), -1.0), "non-negative"),
        (False, np.ones((3, 3)), "zero diagonal"),
        (False, np.inf, "finite"),
        (True, np.nan, "finite"),
        (False, -1.0, "non-negative"),
        (True, -1.0, "non-negative"),
    ],
)
def test_omloc_rejects_invalid_substitution_specs(seqmat, subm_timevarying, subm_base, match):
    with pytest.raises(ValueError, match=match):
        omloc(seqmat, subm_timevarying=subm_timevarying, subm_base=subm_base)


@pytest.mark.parametrize(
    ("weights", "match"),
    [
        (np.ones(2), "length"),
        (np.ones((3, 1)), "one-dimensional"),
        (np.array([1.0, -1.0, 1.0]), "non-negative"),
        (np.array([1.0, np.nan, 1.0]), "finite"),
        (np.zeros(3), "positive total"),
    ],
)
def test_omloc_rejects_invalid_trate_weights(seqmat, weights, match):
    with pytest.raises(ValueError, match=match):
        omloc(seqmat, trate_weights=weights)


def test_omloc_rejects_trate_weights_for_non_trate_substitution(seqmat):
    with pytest.raises(ValueError, match="only valid with subm_base='TRATE'"):
        omloc(seqmat, subm_base="CONSTANT", trate_weights=np.ones(len(seqmat)))


@pytest.mark.parametrize(
    "small",
    [
        np.array([[0]], dtype=np.int8),
        np.array([[0], [0]], dtype=np.uint8),
        np.array([[0, 3], [3, 0]], dtype=np.int64),
    ],
)
def test_omloc_handles_single_state_length_and_sparse_code_boundaries(small):
    full = omloc(small, subm_base="CONSTANT", norm=None)

    assert full.shape == (len(small), len(small))
    assert np.isfinite(full).all()
    np.testing.assert_array_equal(np.diag(full), 0)


def test_omloc_accepts_array_like_integer_input(seqmat):
    expected = omloc(seqmat, subm_base=2.0)
    actual = omloc(seqmat.tolist(), subm_base=2.0)

    np.testing.assert_allclose(actual, expected, rtol=0, atol=0)
