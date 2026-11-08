"""Direct native kernels must validate state dtype and range before narrowing.

Private LCS still accepts arbitrary signed int32 labels because it uses their
equality relation, while matrix-cost kernels retain their alphabet bounds.
"""

import numpy as np
import pytest

PB = pytest.importorskip("traminex._pb_kernels", reason="compiled _pb_kernels not built")
METHODS = ("ham_unit", "om", "om_mat", "lcs", "lcp", "rlcp", "dhd")


def call_native(method, codes):
    kernel = getattr(PB, "pairwise_" + method)
    if method == "om_mat":
        return kernel(codes, 1.0, np.array([[0.0, 2.0], [2.0, 0.0]]))
    if method == "dhd":
        return kernel(codes, np.full((codes.shape[1], 2), 0.5), 0)
    return kernel(codes)


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("dtype", [np.int64, np.uint64], ids=["signed", "unsigned"])
def test_native_rejects_oversized_integer_before_int32_wrap(method, dtype):
    # 2**32 silently became zero under forcecast, changing categorical identity
    # and evading even the matrix-cost kernels' post-conversion bounds checks.
    codes = np.array([[0, 1], [2**32, 1]], dtype=dtype)
    with pytest.raises((ValueError, OverflowError), match="code|int32|range"):
        call_native(method, codes)


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("value", [0.0, 0.5], ids=["integral-float", "fractional-float"])
def test_native_requires_integer_state_dtype(method, value):
    codes = np.array([[0, 1], [value, 1]], dtype=np.float64)
    with pytest.raises((TypeError, ValueError), match="integer|dtype|code"):
        call_native(method, codes)


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("dtype", [np.int32, np.int64, np.uint32, np.uint64])
@pytest.mark.parametrize("layout", ["contiguous", "readonly", "strided"])
def test_native_preserves_valid_integer_codes_for_supported_array_layouts(method, dtype, layout):
    codes = np.array([[0, 1, 0], [1, 1, 0]], dtype=dtype)
    if layout == "readonly":
        codes.setflags(write=False)
    elif layout == "strided":
        storage = np.full((2, 6), 1, dtype=dtype)
        storage[:, ::2] = codes
        codes = storage[:, ::2]
        assert not codes.flags.c_contiguous
    # Hand-derived raw distances for one mismatch at the first position:
    # prefix length 0, suffix length 2, LCS length 2, DHD mean=(1-.5)/3.
    distance = {"ham_unit": 1, "om": 2, "om_mat": 2, "lcs": 2, "lcp": 6, "rlcp": 2, "dhd": 1 / 6}[
        method
    ]
    np.testing.assert_array_equal(call_native(method, codes), [[0, distance], [distance, 0]])


@pytest.mark.parametrize(
    "codes,want",
    [
        ([[-1, -1], [-1, -1]], 0),
        ([[-1, 0], [0, -1]], 2),
        ([[0, 1000], [1000, 0]], 2),
        ([[-(2**31), 2**31 - 1], [2**31 - 1, -(2**31)]], 2),
    ],
    ids=["identical-negative", "negative", "sparse", "int32-extremes"],
)
def test_private_lcs_keeps_valid_negative_and_sparse_int32_labels(codes, want):
    np.testing.assert_array_equal(
        PB.pairwise_lcs(np.array(codes, dtype=np.int32)), [[0, want], [want, 0]]
    )
