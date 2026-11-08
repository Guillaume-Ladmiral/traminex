"""RP9 byte-identity differential: the compiled accel kernels must equal the
pure-Python reference **exactly** (np.array_equal, i.e. atol=0).

Two layers per method:
  * direct-kernel: `_pb_kernels.pairwise_*` vs `python_kernels` raw distance -- always
    runs when the compiled module is built, independent of TRAMINEX_BINDING.
  * dispatcher: `traminex.distances.pairwise(...)` vs `python_kernels.pairwise(...)` --
    exercises the full seqdist routing + normalisation; runs when the active backend is
    the compiled one (the default once the kernel is built, or TRAMINEX_BINDING=pybind11).

The whole module skips cleanly when the kernel is not built, so the default pure-Python
install stays green.
"""

from __future__ import annotations

import numpy as np
import pytest

from traminex import distances as D
from traminex.distances import python_kernels as PY

try:
    from traminex import _pb_kernels as PB
except ImportError:  # pragma: no cover - exercised only on pure-Python installs
    PB = None

pytestmark = pytest.mark.skipif(PB is None, reason="compiled _pb_kernels not built")

COMPILED = D._impl_name() == "pybind11"


def _corpora():
    """Randomised corpus + edge cases, as int32 [n, L] matrices."""
    rng = np.random.default_rng(20260625)
    mats = []
    for n, L, K in [(40, 12, 3), (60, 25, 6), (25, 40, 9), (80, 8, 2), (12, 50, 12)]:
        mats.append(rng.integers(0, K, size=(n, L), dtype=np.int32))
    # edge cases
    mats.append(np.zeros((1, 10), dtype=np.int32))  # single sequence
    mats.append(rng.integers(0, 5, size=(6, 0), dtype=np.int32))  # empty length L=0
    mats.append(np.full((5, 7), 3, dtype=np.int32))  # all identical
    mats.append(np.array([[0, 1, 2, 3, 4]], dtype=np.int32))  # single, distinct states
    # a "missing"-like high code mixed in (TraMineR treats nr as an ordinary state)
    m = rng.integers(0, 4, size=(20, 15), dtype=np.int32)
    m[m == 0] = 99
    mats.append(m)
    return mats


CORPORA = _corpora()


# --------------------------------------------------------------------------- HAM


@pytest.mark.parametrize("mat", CORPORA)
def test_ham_unit_kernel_byte_identical(mat):
    """Direct kernel: raw unit-cost Hamming == python_kernels raw HAM (mismatch count)."""
    compiled = PB.pairwise_ham_unit(mat)
    reference = PY.pairwise(mat, method="HAM")  # subm=None, norm=None -> raw mismatch count
    assert np.array_equal(np.asarray(compiled, float), np.asarray(reference, float))


@pytest.mark.skipif(not COMPILED, reason="active backend is not the compiled kernel")
@pytest.mark.parametrize("mat", CORPORA)
@pytest.mark.parametrize("norm", [None, "none", "maxlength", "auto"])
def test_ham_dispatcher_byte_identical(mat, norm):
    """Dispatcher: seqdist HAM (raw + length-normalised) == pure-Python, exactly."""
    compiled = D.pairwise(mat, method="HAM", norm=norm)
    reference = PY.pairwise(mat, method="HAM", norm=norm)
    assert np.array_equal(np.asarray(compiled, float), np.asarray(reference, float))


# --------------------------------------------------------------------------- LCS


@pytest.mark.parametrize("mat", CORPORA)
def test_lcs_kernel_byte_identical(mat):
    """Direct kernel: raw LCS distance (2L - 2*lcs) == python_kernels raw LCS."""
    compiled = PB.pairwise_lcs(mat)
    reference = PY.pairwise(mat, method="LCS")  # norm=None -> raw L1+L2-2*lcs
    assert np.array_equal(np.asarray(compiled, float), np.asarray(reference, float))


@pytest.mark.skipif(not COMPILED, reason="active backend is not the compiled kernel")
@pytest.mark.parametrize("mat", CORPORA)
@pytest.mark.parametrize("norm", [None, "length"])
def test_lcs_dispatcher_byte_identical(mat, norm):
    """Dispatcher: seqdist LCS (raw + length-normalised) == pure-Python, exactly."""
    compiled = D.pairwise(mat, method="LCS", norm=norm)
    reference = PY.pairwise(mat, method="LCS", norm=norm)
    assert np.array_equal(np.asarray(compiled, float), np.asarray(reference, float))


# ---------------------------------------------------------------------------- OM


def _subm_matrix(K, *, seed, symmetric):
    rng = np.random.default_rng(seed)
    M = rng.uniform(0.5, 4.0, size=(K, K))
    if symmetric:
        M = (M + M.T) / 2.0
    np.fill_diagonal(M, 0.0)  # OM substitution matrices are zero-diagonal (TraMineR)
    return np.ascontiguousarray(M, dtype=np.float64)


# (indel, subm) scalar variants -- include a NON-INTEGER indel (1.3) to stress the
# repeated-addition border accumulation order (where i*indel != indel+indel+...).
OM_SCALAR = [(1.0, 2.0), (1.0, 3.5), (1.3, 2.7), (0.7, 1.9)]


@pytest.mark.parametrize("mat", CORPORA)
@pytest.mark.parametrize("indel,subm", OM_SCALAR)
def test_om_scalar_kernel_byte_identical(mat, indel, subm):
    """Direct kernel: raw scalar OM == python_kernels raw OM (no normalisation)."""
    compiled = PB.pairwise_om(mat, indel, subm)
    reference = PY.pairwise(mat, method="OM", indel=indel, subm=subm)
    assert np.array_equal(np.asarray(compiled, float), np.asarray(reference, float))


@pytest.mark.parametrize("mat", CORPORA)
@pytest.mark.parametrize("symmetric", [True, False])
def test_om_matrix_kernel_byte_identical(mat, symmetric):
    """Direct kernel: raw matrix OM == python_kernels raw OM, zero-diagonal subm."""
    K = int(np.asarray(mat).max()) + 1 if np.asarray(mat).size else 1
    K = max(K, 2)
    M = _subm_matrix(K, seed=7, symmetric=symmetric)
    for indel in (1.0, 1.3):
        compiled = PB.pairwise_om_mat(mat, indel, M)
        reference = PY.pairwise(mat, method="OM", indel=indel, subm=M)
        assert np.array_equal(np.asarray(compiled, float), np.asarray(reference, float))


@pytest.mark.skipif(not COMPILED, reason="active backend is not the compiled kernel")
@pytest.mark.parametrize("mat", CORPORA)
@pytest.mark.parametrize("norm", [None, "length", "max"])
def test_om_dispatcher_scalar_byte_identical(mat, norm):
    """Dispatcher: seqdist OM scalar (all norms) == pure-Python, exactly."""
    for indel, subm in OM_SCALAR:
        compiled = D.pairwise(mat, method="OM", indel=indel, subm=subm, norm=norm)
        reference = PY.pairwise(mat, method="OM", indel=indel, subm=subm, norm=norm)
        assert np.array_equal(np.asarray(compiled, float), np.asarray(reference, float))


@pytest.mark.skipif(not COMPILED, reason="active backend is not the compiled kernel")
@pytest.mark.parametrize("mat", CORPORA)
@pytest.mark.parametrize("norm", [None, "length", "max"])
def test_om_dispatcher_matrix_byte_identical(mat, norm):
    """Dispatcher: seqdist OM matrix subm (all norms) == pure-Python, exactly."""
    arr = np.asarray(mat)
    K = max(int(arr.max()) + 1 if arr.size else 2, 2)
    M = _subm_matrix(K, seed=11, symmetric=True)
    compiled = D.pairwise(mat, method="OM", indel=1.3, subm=M, norm=norm)
    reference = PY.pairwise(mat, method="OM", indel=1.3, subm=M, norm=norm)
    assert np.array_equal(np.asarray(compiled, float), np.asarray(reference, float))


# ----------------------------------------------------------- direct input safety


@pytest.mark.parametrize(
    "mat",
    [
        np.array([[-1, 0]], dtype=np.int32),
        np.array([[0, 2]], dtype=np.int32),
    ],
)
def test_om_matrix_kernel_rejects_state_codes_outside_substitution_matrix(mat):
    with pytest.raises(ValueError, match=r"state codes.*\[0, 2\)"):
        PB.pairwise_om_mat(mat, 1.0, np.eye(2))


@pytest.mark.parametrize(
    "mat",
    [
        np.array([[-1, 0]], dtype=np.int32),
        np.array([[0, 2]], dtype=np.int32),
    ],
)
def test_dhd_kernel_rejects_state_codes_outside_frequency_matrix(mat):
    with pytest.raises(ValueError, match=r"state codes.*\[0, 2\)"):
        PB.pairwise_dhd(mat, np.ones((2, 2)), 0)


@pytest.mark.parametrize(
    "mat",
    [
        np.array([[-1, -1], [-1, -1]], dtype=np.int32),
        np.array([[-1, 0], [0, -1]], dtype=np.int32),
        np.array([[0, 1000], [1000, 0]], dtype=np.int32),
    ],
)
def test_lcs_kernel_remaps_negative_or_pathologically_sparse_state_codes(mat):
    # The stable Python dispatcher rejects negative categorical state codes.
    # Coordinate compression preserves exactly the equality relation that LCS
    # consumes while providing a valid independent dispatcher oracle.
    _, inverse = np.unique(mat, return_inverse=True)
    expected = PY.pairwise(inverse.reshape(mat.shape), method="LCS")
    actual = PB.pairwise_lcs(mat)
    assert np.array_equal(np.asarray(actual, float), np.asarray(expected, float))


@pytest.mark.parametrize("variant", [-1, 2])
def test_dhd_kernel_rejects_unknown_variant(variant):
    with pytest.raises(ValueError, match="variant"):
        PB.pairwise_dhd(np.array([[0, 1]], dtype=np.int32), np.ones((2, 2)), variant)


# ------------------------------------------------------------------ no surface change


def test_accel_adds_no_public_symbol():
    """RP9 is internal acceleration: the public API must not grow (gate 2)."""
    import traminex
    from traminex.api_manifest import PUBLIC_API

    assert traminex.__all__ == PUBLIC_API
    # the compiled-module diagnostics/kernels live on _pb_kernels, never on the public API
    for name in ("accel_max_threads", "accel_openmp_enabled", "pairwise_ham_unit"):
        assert name not in traminex.__all__


def test_impl_name_is_valid():
    """_impl_name() reports a known backend; with the kernel built it is the compiled one."""
    assert D._impl_name() in {"python", "pybind11"}
    assert hasattr(PB, "accel_openmp_enabled")
