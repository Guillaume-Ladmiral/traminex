import numpy as np
import pytest

import traminex.distances.python_kernels as python_kernels
from traminex import seqdef
from traminex.distances import pairwise
from traminex.distances.python_kernels import pairwise as py_pairwise


def test_duration_based_distances():
    mat = np.array([[0, 1, 2, 2], [0, 1, 1, 2], [2, 2, 1, 0]], dtype=int)
    x = seqdef(mat)
    for metric in ["l1", "cos", "chi2"]:
        D = py_pairwise(x.data, method="DUR", metric=metric)
        assert D.shape == (3, 3)


def test_jsd_transition_distance():
    mat = np.array([[0, 1, 2, 2], [0, 1, 1, 2]], dtype=int)
    x = seqdef(mat)
    D = py_pairwise(x.data, method="JSD_TRANS")
    assert D.shape == (2, 2)


def test_public_pairwise_extension_defaults_and_refusals_are_actionable():
    mat = np.array([[0, 1, 2, 2], [0, 1, 1, 2], [2, 2, 1, 0]], dtype=int)

    wham = pairwise(mat, method="WHAM")
    expected = np.mean(mat[0] != mat[1])
    assert wham[0, 1] == expected

    with pytest.raises(ValueError, match="length equal to sequence length"):
        pairwise(mat, method="WHAM", weights=np.ones(mat.shape[1] - 1))
    with pytest.raises(ValueError, match="Unsupported DUR metric"):
        pairwise(mat, method="DUR", metric="not-a-metric")
    with pytest.raises(ValueError, match="Unknown distance method"):
        pairwise(mat, method="not-a-method")


@pytest.mark.parametrize(
    "weights",
    [
        np.array([1.0, -1.0, 1.0, 1.0]),
        np.array([1.0, np.nan, 1.0, 1.0]),
        np.array([1.0, np.inf, 1.0, 1.0]),
        np.zeros(4),
    ],
)
def test_wham_rejects_weights_that_cannot_define_a_weighted_mean(weights):
    mat = np.array([[0, 1, 2, 2], [0, 1, 1, 2]], dtype=int)

    with pytest.raises(ValueError, match="finite, non-negative.*positive sum"):
        pairwise(mat, method="WHAM", weights=weights)


def test_jsd_transition_python_fallback_honors_eps(monkeypatch):
    mat = np.array([[0, 0, 0, 0], [0, 1, 1, 1]], dtype=int)
    eps = 0.25
    expected = python_kernels.jsd_transition_distance(mat[0], mat[1], 2, eps=eps)
    monkeypatch.setattr(python_kernels, "_HAVE_NUMBA", False)

    actual = pairwise(mat, method="JSD_TRANS", eps=eps)

    assert actual[0, 1] == pytest.approx(expected)


@pytest.mark.parametrize("eps", [-1.0, np.nan, np.inf])
def test_jsd_transition_rejects_invalid_smoothing(eps):
    mat = np.array([[0, 0, 0, 0], [0, 1, 1, 1]], dtype=int)

    with pytest.raises(ValueError, match="eps must be finite and non-negative"):
        pairwise(mat, method="JSD_TRANS", eps=eps)
