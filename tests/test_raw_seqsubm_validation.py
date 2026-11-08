"""Raw substitution-cost entry points enforce the integer-state contract."""

import numpy as np
import pytest

import traminex as tx
from traminex.distances import costs

METHODS = ("TRATE", "CONSTANT", "INDELS", "INDELSLOG", "FEATURES", "FUTURE")
BUILDERS = (
    costs.trate_substitution_matrix,
    costs.constant_substitution_matrix,
    costs.indels_costs,
    costs.indelslog_costs,
    costs.future_costs,
    costs.trate_substitution_matrix_timevarying,
    costs.local_indel_by_pos,
)


@pytest.mark.parametrize("method", METHODS)
def test_raw_seqsubm_rejects_fractional_state_codes(method):
    kwargs = {"state_features": [[0.0], [1.0]]} if method == "FEATURES" else {}
    with pytest.raises(TypeError):
        tx.seqsubm(np.array([[0.0, 1.5]]), method=method, **kwargs)


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("matrix", [np.array([[-1, 1]]), np.array([0, 1])])
def test_raw_seqsubm_rejects_negative_codes_and_nonmatrix_input(method, matrix):
    kwargs = {"state_features": [[0.0], [1.0]]} if method == "FEATURES" else {}
    with pytest.raises(ValueError):
        tx.seqsubm(matrix, method=method, **kwargs)


@pytest.mark.parametrize("builder", BUILDERS)
def test_direct_cost_builder_rejects_fractional_state_codes(builder):
    with pytest.raises(TypeError):
        builder(np.array([[0.0, 1.5]]))


@pytest.mark.parametrize("builder", BUILDERS)
def test_direct_cost_builder_rejects_negative_state_codes(builder):
    with pytest.raises(ValueError):
        builder(np.array([[-1, 1]]))


@pytest.mark.parametrize("direct", [False, True])
@pytest.mark.parametrize("weights", [None, [1.0, 3.0]])
def test_transition_costs_widen_small_integer_dtypes_before_counting(direct, weights):
    matrix = np.array([[0, 16], [16, 0]], dtype=np.uint8)
    before = matrix.copy()
    if direct:
        actual = costs.trate_substitution_matrix(matrix, weights=weights)
    else:
        actual = tx.seqsubm(matrix, "TRATE", weights=weights)
    # Each observed state always transitions to the other, so both rates are
    # 1 and cost(0,16) = 2 - 1 - 1 = 0. Computing 16*17 in uint8 wraps to 16,
    # inventing a second 0->16 transition and losing the 16->0 transition.
    expected = np.full((17, 17), 2.0)
    np.fill_diagonal(expected, 0.0)
    expected[0, 16] = expected[16, 0] = 0.0
    np.testing.assert_array_equal(actual, expected)
    np.testing.assert_array_equal(matrix, before)


@pytest.mark.parametrize("weights", [None, [1.0, 3.0]])
def test_timevarying_transition_costs_widen_small_integer_dtypes(weights):
    matrix = np.array([[0, 16], [16, 0]], dtype=np.uint8)
    actual = costs.trate_substitution_matrix_timevarying(matrix, weights=weights)
    expected = np.full((2, 17, 17), 2.0)
    for position in range(2):
        np.fill_diagonal(expected[position], 0.0)
        expected[position, 0, 16] = expected[position, 16, 0] = 0.0
    np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize("method", ["INDELS", "INDELSLOG"])
@pytest.mark.parametrize("direct", [False, True])
def test_frequency_costs_accept_valid_uint64_states_on_supported_numpy_versions(method, direct):
    matrix = np.array([[0, 1], [1, 0]], dtype=np.uint64)
    if direct:
        builder = costs.indels_costs if method == "INDELS" else costs.indelslog_costs
        actual = builder(matrix)[1]
    else:
        actual = tx.seqsubm(matrix, method)
    # Each state has frequency 1/2: INDELS adds 2+2; INDELSLOG adds two
    # log(2/(1+1/2)) terms. Older NumPy rejects uint64 bincount input even
    # when all values are small enough to fit the established state range.
    off_diagonal = 4.0 if method == "INDELS" else 2 * np.log(4 / 3)
    np.testing.assert_allclose(actual, [[0, off_diagonal], [off_diagonal, 0]], rtol=0, atol=1e-15)
