"""Distinct-state and warning controls for the October remediation."""

import warnings
from types import SimpleNamespace

import numpy as np
import pytest

from traminex import seqalign, seqdef, seqdist
from traminex.distances import pairwise, pairwise_between_sets
from traminex.distances import python_kernels as kernels


@pytest.mark.parametrize("dtype", [np.int64, np.uint64])
@pytest.mark.parametrize("side", [0, 1])
def test_between_sets_rejects_codes_before_narrowing(dtype, side):
    sets = [[np.array([0], dtype=dtype)], [np.array([2**32, 1], dtype=dtype)]]
    if side:
        sets.reverse()
    with pytest.raises(ValueError, match="32-bit"):
        pairwise_between_sets(*sets, method="LCS")
    assert pairwise_between_sets([[0]], [[2, 1]], method="LCS")[0, 0] == 3


def test_void_dispatch_rejects_overflow_before_deleting_padding():
    with pytest.raises(ValueError, match="32-bit"):
        pairwise(np.array([[0, 9], [2**32, 1]]), method="LCS", void=9)


@pytest.mark.parametrize("method", ["OM", "LCS", "LCP", "RLCP"])
@pytest.mark.parametrize("mode", ["gmean", "auto", "none", "maxlength"])
def test_low_level_norm_refusal_does_not_depend_on_padding(method, mode):
    for matrix, options in [(np.array([[0], [1]]), {}), (np.array([[0, 9], [1, 9]]), {"void": 9})]:
        with pytest.raises(ValueError, match="low-level"):
            pairwise(matrix, method=method, norm=mode, **options)
    with pytest.raises(ValueError, match="low-level"):
        pairwise_between_sets([[0]], [[1, 1]], method=method, norm=mode)


@pytest.mark.parametrize("mode", [None, "length"])
def test_padding_and_readonly_stride_preserve_valid_lcs(mode):
    padded = np.array([[0, 9, 9, 9], [1, 9, 9, 9]])[:, ::2]
    padded.flags.writeable = False
    np.testing.assert_array_equal(
        pairwise(padded, method="LCS", void=9, norm=mode),
        pairwise(np.array([[0], [1]]), method="LCS", norm=mode),
    )


@pytest.mark.parametrize("rows", [[np.array([-1])], [np.array([0.5])], [np.array([[0]])]])
def test_ragged_invalid_shape_dtype_and_negative_codes(rows):
    with pytest.raises((TypeError, ValueError)):
        pairwise_between_sets(rows, [[0, 1]], method="LCS")


def test_row_engine_rejects_wrapping_codes():
    with pytest.raises(ValueError, match="32-bit"):
        kernels._pairwise_rows([np.array([0]), np.array([2**32, 1])], "LCS")


def test_empty_alignment_can_be_displayed():
    result = seqalign(seqdef(np.empty((2, 0), dtype=object), alphabet=["A"]), (0, 1), indel=1)
    assert result.total_cost == 0
    assert "0" in repr(result)


def test_jsd_zero_smoothing_is_finite_and_warning_free():
    x = seqdef([["A", "A"], ["B", "B"], ["A", "A"]])
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        actual = seqdist(x, method="JSD_TRANS", eps=0)
    expected = np.array([[0, np.log(2), 0], [np.log(2), 0, np.log(2)], [0, np.log(2), 0]])
    np.testing.assert_allclose(actual, expected, rtol=1e-14, atol=0)


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ At\ least\ one\ substitution\ cost\ greater\ than\ twice\ the\ indel\ cost\.:traminex.distances.python_kernels.CostDiagnosticMessage"
)
def test_triangle_diagnostic_handles_finite_addition_overflow():
    matrix = np.full((3, 3), 1e308)
    np.fill_diagonal(matrix, 0)
    x = seqdef([["A"], ["B"]], alphabet=["A", "B", "C"])
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        actual = seqdist(x, method="OM", subm=matrix, indel=1)
    assert actual[0, 1] == 2


@pytest.mark.parametrize("dtype", [np.int64, np.uint64])
def test_reverse_prefix_wrapper_validates_before_int32_narrowing(dtype):
    from traminex import seqRLCP

    with pytest.raises(ValueError, match="32-bit"):
        seqRLCP(np.array([[0, 1], [2**32, 1]], dtype=dtype))
    np.testing.assert_array_equal(
        seqRLCP(np.array([[0, 1], [2, 1]], dtype=dtype)), [[0, 2], [2, 0]]
    )


@pytest.mark.parametrize("function", ["seqLLCS", "seqLLCP"])
def test_scalar_common_helpers_refuse_fractional_state_codes(function):
    import traminex as tx

    with pytest.raises(TypeError, match="integer"):
        getattr(tx, function)(np.array([0.1]), np.array([0.9]))


@pytest.mark.parametrize("function", ["seqLLCS", "seqLLCP"])
def test_scalar_common_helpers_preserve_large_unsigned_symbol_identity(function):
    import traminex as tx

    # These pure scalar comparisons have no int32 kernel boundary. Retain
    # valid integer symbols in their source dtype rather than narrowing them.
    a = np.array([2**63 + 3], dtype=np.uint64)
    b = np.array([3], dtype=np.uint64)
    assert getattr(tx, function)(a, b) == 0
    assert getattr(tx, function)(a, a) == 1


@pytest.mark.parametrize("deduplicate", [False, True])
@pytest.mark.parametrize("invalid", [np.array([0.5]), np.array([False])], ids=["float", "bool"])
def test_rectangular_ragged_dedup_cannot_discard_invalid_dtype(deduplicate, invalid):
    # The valid left row establishes key (0,) first. Coercing the invalid
    # right row to int for the key must not remove it before dtype validation.
    with pytest.raises(TypeError, match="integer"):
        pairwise_between_sets(
            np.array([[0]]),
            [invalid, np.array([0, 1])],
            method="LCS",
            deduplicate=deduplicate,
        )


@pytest.mark.parametrize("deduplicate", [False, True])
@pytest.mark.parametrize(
    "left,right,void",
    [
        (np.array([[0, 1]]), np.array([[0.5, 1, 9]]), 9),
        (np.array([[0]]), np.array([[False, True]]), 1),
    ],
    ids=["float", "bool"],
)
def test_rectangular_void_dedup_retains_original_dtype_validation(deduplicate, left, right, void):
    with pytest.raises(TypeError, match="integer"):
        pairwise_between_sets(left, right, method="LCS", void=void, deduplicate=deduplicate)


@pytest.mark.parametrize("deduplicate", [False, True])
@pytest.mark.parametrize("dtype", [np.int64, np.uint64])
def test_rectangular_void_code_cannot_hide_integer_overflow(deduplicate, dtype):
    with pytest.raises(ValueError, match="32-bit"):
        pairwise_between_sets(
            np.array([[0]], dtype=dtype),
            np.array([[0, 2**32]], dtype=dtype),
            method="LCS",
            void=2**32,
            deduplicate=deduplicate,
        )


@pytest.mark.parametrize("deduplicate", [False, True])
def test_rectangular_valid_ragged_duplicates_keep_expected_distance(deduplicate):
    np.testing.assert_array_equal(
        pairwise_between_sets(
            np.array([[0]]),
            [np.array([0]), np.array([0, 1])],
            method="LCS",
            deduplicate=deduplicate,
        ),
        [[0, 1]],
    )


@pytest.mark.parametrize(
    "dtype,value,error,match",
    [
        (np.int64, 2**32, ValueError, "32-bit"),
        (np.uint64, 2**32, ValueError, "32-bit"),
        (np.float64, 0.5, TypeError, "integer"),
        (np.bool_, False, TypeError, "integer"),
    ],
    ids=["signed-overflow", "unsigned-overflow", "fractional", "boolean"],
)
def test_external_generic_reference_validates_before_narrowing(dtype, value, error, match):
    reference = SimpleNamespace(data=np.array([[value]], dtype=dtype))
    with pytest.raises(error, match=match):
        seqdist(np.array([[0]]), method="OM", subm=2, indel=1, refseq=reference)


@pytest.mark.parametrize("dtype", [np.int64, np.uint64])
def test_external_generic_reference_preserves_valid_integer_symbol(dtype):
    reference = SimpleNamespace(data=np.array([[1]], dtype=dtype))
    np.testing.assert_array_equal(
        seqdist(np.array([[0]]), method="OM", subm=2, indel=1, refseq=reference), [2]
    )


@pytest.mark.parametrize("deduplicate", [False, True])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("dtype", [np.int8, np.uint8, np.int32, np.int64, np.uint64])
def test_dense_between_sets_validates_original_integer_types(dtype, reverse, deduplicate):
    a = np.array([[0, 9, 1, 9], [0, 9, 1, 9]], dtype=dtype)[:, ::2]
    b = np.array([[1, 0]], dtype=np.int64)
    a.flags.writeable = False
    left, right = (b, a) if reverse else (a, b)
    expected = np.full((len(left), len(right)), 2.0)
    np.testing.assert_array_equal(
        pairwise_between_sets(left, right, method="LCS", deduplicate=deduplicate), expected
    )
    np.testing.assert_array_equal(a, [[0, 1], [0, 1]])


@pytest.mark.parametrize("deduplicate", [False, True])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize(
    "invalid",
    [
        np.array([[True, False]]),
        np.array([[0.0, 1.0]]),
        np.array([[0, 2**32]], dtype=np.uint64),
        np.array([[0, -1]]),
    ],
)
def test_dense_between_sets_rejects_original_invalid_codes(invalid, reverse, deduplicate):
    valid = np.array([[0, 1]], dtype=np.int32)
    left, right = (valid, invalid) if reverse else (invalid, valid)
    with pytest.raises((TypeError, ValueError)):
        pairwise_between_sets(left, right, method="LCS", deduplicate=deduplicate)


@pytest.mark.parametrize(
    "values,n", [([7], 3), ([1, 2], 3), ([1, 2, 3, 4], 3), ([[1, 2, 3]], 3), ([[1]], None)]
)
def test_condensed_conversion_rejects_shape_and_cardinality(values, n):
    from traminex.distances import from_condensed

    with pytest.raises(ValueError):
        from_condensed(values, n=n)


@pytest.mark.parametrize("n", [True, np.bool_(False), -1, 2.5, 3.0, "3"])
def test_condensed_conversion_rejects_invalid_dimension(n):
    from traminex.distances import from_condensed

    with pytest.raises((TypeError, ValueError)):
        from_condensed([1, 2, 3], n=n)


@pytest.mark.parametrize("n", [0, 1, 2, 3, np.int64(4)])
def test_condensed_conversion_valid_roundtrip_and_order(n):
    from traminex.distances import from_condensed, to_condensed

    values = np.arange(n * (n - 1) // 2, dtype=float) + 1
    matrix = from_condensed(values, n=n)
    assert matrix.shape == (n, n)
    np.testing.assert_array_equal(to_condensed(matrix), values)
    np.testing.assert_array_equal(matrix, matrix.T)
    np.testing.assert_array_equal(np.diag(matrix), np.zeros(n))
    if n == 3:
        np.testing.assert_array_equal(matrix, [[0, 1, 2], [1, 0, 3], [2, 3, 0]])


@pytest.mark.parametrize(
    "a,b",
    [
        (1e308, 1e308),
        (1e308, 8e307),
        (1e308, -1e308),
        (1.0, np.nextafter(1.0, 2.0)),
        (3.0, 7.0),
        (0.0, 0.0),
        (np.nextafter(0.0, 1.0), np.nextafter(0.0, 1.0)),
        (np.nextafter(0.0, 1.0), 2 * np.nextafter(0.0, 1.0)),
        (np.finfo(float).max, np.nextafter(np.finfo(float).max, 0.0)),
    ],
)
def test_sanitization_uses_finite_correctly_rounded_symmetric_average(a, b):
    from fractions import Fraction

    from traminex.distances import sanitize_distance_matrix

    matrix = np.array([[9.0, a], [b, 9.0]])
    original = matrix.copy()
    expected = float((Fraction(float(a)) + Fraction(float(b))) / 2)
    with np.errstate(over="raise", invalid="raise"):
        actual = sanitize_distance_matrix(matrix)
    np.testing.assert_array_equal(actual, [[0, expected], [expected, 0]])
    np.testing.assert_array_equal(matrix, original)


def test_sanitization_preserves_nonfinite_replacement_contract():
    from traminex.distances import sanitize_distance_matrix

    np.testing.assert_array_equal(
        sanitize_distance_matrix([[0, np.inf], [4, np.nan]]), [[0, 4], [4, 0]]
    )
    np.testing.assert_array_equal(
        sanitize_distance_matrix(np.full((2, 2), np.nan)), np.zeros((2, 2))
    )


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ Empty\ sequence\(s\)\ \ \d+,\ \d+\ \ may\ cause\ inconsistent\ results!:UserWarning"
)
@pytest.mark.parametrize("subm", [2, np.array([[0, 2], [2, 0]]), "TRATE", "CONSTANT"])
@pytest.mark.parametrize("external", [False, True])
def test_empty_om_references_preserve_cost_forms_and_output_shape(subm, external):
    from traminex import STSList

    x = STSList(np.empty((2, 0), dtype=np.int32), ["A", "B"])
    reference = STSList(np.empty((1, 0), dtype=np.int32), ["A", "B"]) if external else 1
    if isinstance(subm, str) and subm == "TRATE":
        # Width 0 has no transition: TRATE costs are refused, as seqcost and
        # seqtrate refuse them (TX-AUD-20261011-040).
        with pytest.raises(ValueError, match="absolute lag must be smaller"):
            seqdist(x, "OM", subm=subm, indel=1, refseq=reference)
        return
    actual = seqdist(x, "OM", subm=subm, indel=1, refseq=reference)
    np.testing.assert_array_equal(actual, [0, 0])


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ Empty\ sequence\(s\)\ \ \d+\ \ may\ cause\ inconsistent\ results!:UserWarning"
)
@pytest.mark.parametrize("subm", [2, np.array([[0, 2], [2, 0]])])
def test_empty_om_external_reference_keeps_nonempty_indel_distance(subm):
    from traminex import STSList

    x = STSList(np.array([[0, 1], [9, 9]]), ["A", "B"], void_code=9)
    reference = STSList(np.empty((1, 0), dtype=np.int32), ["A", "B"])
    np.testing.assert_array_equal(seqdist(x, "OM", subm=subm, indel=1, refseq=reference), [2, 0])


@pytest.mark.parametrize(
    "norm", [None, "none", "auto", "maxlength", "length", "gmean", "maxdist", "max", "YujianBo"]
)
@pytest.mark.parametrize("width", [0, 2])
def test_empty_public_hamming_accepts_supported_normalizations(norm, width):
    from traminex import STSList

    x = STSList(np.empty((0, width), dtype=np.int32), ["A", "B"])
    assert seqdist(x, "HAM", norm=norm).shape == (0, 0)
    assert seqdist(x, "HAM", norm=norm, condensed=True).shape == (0,)
    # Empty reference index sets are explicitly unsupported, even for no rows.
    with pytest.raises(ValueError, match="indices"):
        seqdist(x, "HAM", norm=norm, ref_sets=([], []))


@pytest.mark.parametrize(
    "options",
    [{"norm": "invalid"}, {"subm": [1, 2]}, {"subm": [[0, -1], [-1, 0]]}, {"bad_option": True}],
)
def test_empty_public_hamming_still_validates_options_and_costs(options):
    from traminex import STSList

    x = STSList(np.empty((0, 2), dtype=np.int32), ["A", "B"])
    with pytest.raises((TypeError, ValueError)):
        seqdist(x, "HAM", **options)


@pytest.mark.parametrize("mode", ["symmetric", "row", "col"])
@pytest.mark.parametrize("normalization", ["none", "mean_offdiag", "clip_0_2"])
@pytest.mark.parametrize("smoothing", [0, 0.5])
@pytest.mark.parametrize("weights", [None, np.array([0.0, 0.0]), np.array([2.0, 1.0])])
def test_raw_singleton_trate_uses_no_transition_formula(mode, normalization, smoothing, weights):
    from traminex.distances.costs import trate_substitution_matrix

    x = np.array([[0], [1]], dtype=np.int32)
    expected_cost = 2 if mode == "symmetric" and normalization != "mean_offdiag" else 1
    expected = [[0, expected_cost], [expected_cost, 0]]
    np.testing.assert_array_equal(
        trate_substitution_matrix(
            x, mode=mode, normalize=normalization, smoothing=smoothing, weights=weights
        ),
        expected,
    )
    np.testing.assert_array_equal(trate_substitution_matrix(np.array([[0], [0]])), [[0]])


def test_raw_singleton_trate_public_builder_and_generated_om_agree_with_formula():
    from traminex import seqsubm

    x = np.array([[0], [1]], dtype=np.int32)
    np.testing.assert_array_equal(seqsubm(x, "TRATE"), [[0, 2], [2, 0]])
    np.testing.assert_array_equal(pairwise(x, "OM", subm="TRATE"), [[0, 2], [2, 0]])


@pytest.mark.parametrize("matrix", [np.empty((0, 0), dtype=int), np.array([[0], [1]])])
@pytest.mark.parametrize("options", [{"mode": "invalid"}, {"normalize": "invalid"}])
def test_raw_no_transition_trate_does_not_bypass_option_validation(matrix, options):
    from traminex.distances.costs import trate_substitution_matrix

    with pytest.raises(ValueError):
        trate_substitution_matrix(matrix, **options)


@pytest.mark.parametrize("width", [0, 1, 3])
def test_empty_timevarying_trate_preserves_position_dimension(width):
    from traminex.distances.costs import trate_substitution_matrix_timevarying

    actual = trate_substitution_matrix_timevarying(np.empty((0, width), dtype=np.int32))
    assert actual.shape == (width, 0, 0)


@pytest.mark.parametrize("kwargs", [{"lag": 2}, {"mode": "unknown"}, {"normalize": "unknown"}])
@pytest.mark.parametrize("width", [0, 3])
def test_empty_timevarying_trate_refuses_unsupported_options(kwargs, width):
    from traminex.distances.costs import trate_substitution_matrix_timevarying

    with pytest.raises(ValueError):
        trate_substitution_matrix_timevarying(np.empty((0, width), dtype=np.int32), **kwargs)


@pytest.mark.parametrize("smoothing", [-1.0, np.nan, np.inf, -np.inf])
@pytest.mark.parametrize("empty", [False, True])
@pytest.mark.parametrize("builder", ["stationary", "timevarying", "config"])
def test_trate_refuses_invalid_smoothing_before_transition_arithmetic(smoothing, empty, builder):
    from traminex.config_dataclasses import TRATEConfig
    from traminex.distances.costs import (
        trate_substitution_matrix,
        trate_substitution_matrix_timevarying,
    )

    rows = np.empty((0, 2), dtype=np.int32) if empty else np.array([[0, 1], [1, 0]], dtype=np.int32)
    with pytest.raises(ValueError, match="smoothing.*finite|smoothing.*non-negative"):
        if builder == "config":
            TRATEConfig(smoothing=smoothing)
        elif builder == "stationary":
            trate_substitution_matrix(rows, smoothing=smoothing)
        else:
            trate_substitution_matrix_timevarying(rows, smoothing=smoothing)
