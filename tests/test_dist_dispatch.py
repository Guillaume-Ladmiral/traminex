import numpy as np
import pytest

from traminex import distances
from traminex.distances import pairwise as pairwise_auto
from traminex.distances.python_kernels import pairwise as pairwise_py


def test_dispatch_fallback_matches_python():
    mat = np.array([[0, 1, 2], [0, 1, 1]], dtype=int)
    D_auto = pairwise_auto(mat, method="OM", indel=1.0, subm=2.0)
    D_py = pairwise_py(mat, method="OM", indel=1.0, subm=2.0)
    assert D_auto.shape == D_py.shape
    assert np.allclose(D_auto, D_py)


def test_backend_selection_distinguishes_auto_fallback_and_explicit_request(monkeypatch):
    monkeypatch.setattr(distances, "_pb", None)

    monkeypatch.setattr(distances, "_backend", "auto")
    assert distances._impl_name() == "python"

    monkeypatch.setattr(distances, "_backend", "python")
    assert distances._impl_name() == "python"

    monkeypatch.setattr(distances, "_backend", "pybind11")
    with pytest.raises(RuntimeError, match="pybind11 backend was requested but is unavailable"):
        distances._impl_name()


def test_backend_selection_rejects_unknown_binding(monkeypatch):
    monkeypatch.setattr(distances, "_backend", "unknown")
    with pytest.raises(ValueError, match="TRAMINEX_BINDING must be one of"):
        distances._impl_name()


def test_dispatch_validation_rejects_invalid_norm_and_trate_options():
    mat = np.array([[0, 1], [1, 0]], dtype=int)

    with pytest.raises(ValueError, match="norm must be one of"):
        pairwise_auto(mat, method="OM", norm="gmean")
    with pytest.raises(ValueError, match="not supported for method LCS"):
        pairwise_auto(mat, method="LCS", norm="max")
    with pytest.raises(TypeError, match="string 'TRATE'"):
        pairwise_auto(mat, method="OM", subm="INVALID")
    with pytest.raises(TypeError, match="require subm='TRATE'"):
        pairwise_auto(mat, method="OM", trate_mode="row")
