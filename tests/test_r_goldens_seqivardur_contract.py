"""R-golden coverage for the complete ``seqivardur`` result object.

Generator: ``scripts/r_goldens/make_goldens_seqivardur_contract.R`` against
TraMineR 2.2-14.
"""

from __future__ import annotations

import numpy as np
import pytest
from r_goldens_utils import R_GOLDEN_ATOL, load_golden, load_golden_lines

from traminex import seqdef, seqivardur

ALPHABET_A = ["A", "B", "C", "D"]
SEQUENCES_A = [
    ["A", "A", "B", "C", "C", "C"],
    ["A", "B", "B", "B", "C", "A"],
    ["D", "D", "D", "D", "D", "D"],
    ["A", "B", "A", "B", "A", "B"],
]
ALPHABET_C = ["A", "B", "C"]
SEQUENCES_C = [
    ["A", "A", "B", "C", "C", "C"],
    ["A", "B", "B", "%", "%", "%"],
    ["A", "*", "A", "B", "%", "%"],
    ["C", "C", "%", "%", "%", "%"],
]


def _fixture_a():
    return seqdef(SEQUENCES_A, alphabet=ALPHABET_A)


def _fixture_c():
    return seqdef(SEQUENCES_C, alphabet=ALPHABET_C)


def _assert_contract(result, golden_name: str) -> None:
    golden = load_golden(golden_name)
    assert golden.shape == (3, 4)
    np.testing.assert_allclose(np.asarray(result), golden[0], rtol=0, atol=R_GOLDEN_ATOL)
    np.testing.assert_allclose(result.vmax, golden[1], rtol=0, atol=R_GOLDEN_ATOL)
    np.testing.assert_allclose(result.meand, golden[2], rtol=0, atol=R_GOLDEN_ATOL)


def test_seqivardur_type1_preserves_values_vmax_and_meand():
    _assert_contract(seqivardur(_fixture_a(), type=1), "seqivardur_contract_A_type1")


def test_seqivardur_type2_preserves_values_vmax_and_meand():
    _assert_contract(seqivardur(_fixture_a(), type=2), "seqivardur_contract_A_type2")


def test_seqivardur_ragged_with_missing_preserves_full_contract():
    _assert_contract(
        seqivardur(_fixture_c(), type=1, include_missing=True),
        "seqivardur_contract_C_wm1",
    )


def test_seqivardur_result_remains_ndarray_compatible():
    result = seqivardur(_fixture_a(), type=2)
    values = np.asarray(result)

    square_root = np.sqrt(result)
    arithmetic = result + 2.0
    np.testing.assert_allclose(square_root, np.sqrt(values), rtol=0, atol=0)
    np.testing.assert_allclose(arithmetic, values + 2.0, rtol=0, atol=0)
    assert type(square_root) is np.ndarray
    assert type(arithmetic) is np.ndarray
    assert not hasattr(arithmetic, "vmax")

    reshaped = result.reshape(2, 2)
    broadcast = result + np.ones((2, 4))
    sorted_values = np.sort(result)
    taken = result.take([2, 0])
    deleted = np.delete(result, 1)
    assert type(reshaped) is np.ndarray
    assert type(broadcast) is np.ndarray
    assert type(sorted_values) is np.ndarray
    assert type(taken) is np.ndarray
    assert type(deleted) is np.ndarray
    assert reshaped.shape == (2, 2)
    assert broadcast.shape == (2, 4)
    assert np.sum(result) == np.sum(values)
    with pytest.raises(AttributeError):
        result.shape = (2, 2)
    with pytest.raises(TypeError, match="in-place operations"):
        np.add(result, 1.0, out=result)
    with pytest.raises(TypeError, match="mutating ndarray method"):
        result.sort()
    with pytest.raises(TypeError, match="mutating ndarray method"):
        result.resize((2, 2))
    with pytest.raises(TypeError, match="mutating ndarray method"):
        result.setflags(write=True)
    with pytest.raises(ValueError):
        np.asarray(result).setflags(write=True)
    with pytest.raises(ValueError):
        result.flags.writeable = True

    sliced = result[[0, 2]]
    np.testing.assert_allclose(np.asarray(sliced), values[[0, 2]], rtol=0, atol=0)
    np.testing.assert_allclose(sliced.vmax, result.vmax[[0, 2]], rtol=0, atol=0)
    np.testing.assert_allclose(sliced.meand, result.meand[[0, 2]], rtol=0, atol=0)

    original_vmax = result.vmax.copy()
    sliced.vmax[0] = -1.0
    np.testing.assert_array_equal(result.vmax, original_vmax)
    copied = result.copy()
    copied.meand[0] = -1.0
    assert result.meand[0] != -1.0


def test_seqivardur_all_void_object_matches_oracle_refusal_boundary():
    oracle = dict(line.split("=", 1) for line in load_golden_lines("seqivardur_contract_errors"))
    assert "error" in oracle["class"].lower()
    assert "subscript out of bounds" in oracle["message"]

    all_void = seqdef([["%", "%", "%", "%"]], alphabet=["A", "B"])
    with pytest.raises(ValueError, match="non-void spell"):
        seqivardur(all_void)


def test_seqivardur_mixed_all_void_row_remains_supported_with_r_metadata():
    oracle = dict(line.split("=", 1) for line in load_golden_lines("seqivardur_contract_errors"))
    mixed = seqdef(
        [["A", "B", "A", "B"], ["%", "%", "%", "%"]],
        alphabet=["A", "B"],
    )

    for type_ in (1, 2):
        result = seqivardur(mixed, type=type_)
        expected = np.array(
            [
                float(oracle[f"mixed_type{type_}_variance"]),
                float(oracle[f"mixed_type{type_}_vmax"]),
                float(oracle[f"mixed_type{type_}_meand"]),
            ]
        )
        actual = np.array([result[1], result.vmax[1], result.meand[1]])
        np.testing.assert_allclose(actual, expected, rtol=0, atol=0, equal_nan=True)
