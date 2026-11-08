from __future__ import annotations

import inspect

import numpy as np
import pytest

import traminex.distances as distances

EXPECTED_PUBLIC_EXPORTS = (
    "om_trate_pairwise",
    "pairwise",
    "pairwise_between_sets",
    "seqdist",
    "to_condensed",
    "from_condensed",
    "sanitize_distance_matrix",
    "get_lcp_length",
    "get_rlcp_length",
    "seqsubm",
    "seqLLCS",
    "seqLLCP",
    "seqRLCP",
)


def test_distances_submodule_has_an_explicit_public_export_inventory():
    assert distances.__all__ == EXPECTED_PUBLIC_EXPORTS
    assert all(callable(getattr(distances, name)) for name in EXPECTED_PUBLIC_EXPORTS)
    assert "pairwise" not in __import__("traminex").__all__
    assert "pairwise_between_sets" not in __import__("traminex").__all__
    assert "python_kernels" not in distances.__all__
    assert "validate_pairwise_kwargs" not in distances.__all__


def test_pairwise_submodule_signatures_are_stable():
    pairwise = inspect.signature(distances.pairwise).parameters
    assert list(pairwise) == ["seqmat", "method", "kwargs"]
    assert pairwise["seqmat"].kind is inspect.Parameter.POSITIONAL_OR_KEYWORD
    assert pairwise["method"].default == "OM"
    assert pairwise["kwargs"].kind is inspect.Parameter.VAR_KEYWORD

    between = inspect.signature(distances.pairwise_between_sets).parameters
    assert list(between) == ["A", "B", "method", "deduplicate", "kwargs"]
    assert between["A"].kind is inspect.Parameter.POSITIONAL_OR_KEYWORD
    assert between["B"].kind is inspect.Parameter.POSITIONAL_OR_KEYWORD
    assert between["method"].kind is inspect.Parameter.KEYWORD_ONLY
    assert between["method"].default == "OM"
    assert between["deduplicate"].kind is inspect.Parameter.KEYWORD_ONLY
    assert between["deduplicate"].default is False
    assert between["kwargs"].kind is inspect.Parameter.VAR_KEYWORD


# A corpus strictly larger than A | B.  Deriving the expected block from
# ``vstack([left, right])`` instead would be circular for any method whose
# parameters are estimated from the data: that is the corpus a re-fitting
# implementation estimates from, so the comparison could not detect the re-fit
# (F-001).  The configurations below are all pair-local, so the cross-block is
# independent of the surrounding rows and the larger corpus is the honest
# baseline.  The dataset-fitted configurations are covered, with ``subm`` left at
# its default, in ``tests/test_distances_ref_sets_dataset_fit.py``.
_CROSS_BLOCK_CORPUS = np.array(
    [
        [0, 0, 1, 1],
        [0, 1, 1, 0],
        [1, 1, 0, 0],
        [0, 0, 1, 0],
        [1, 0, 0, 1],
        [1, 1, 1, 0],
        [0, 1, 0, 1],
    ],
    dtype=np.int32,
)
_CROSS_BLOCK_A = np.array([0, 1])
_CROSS_BLOCK_B = np.array([2, 3])


@pytest.mark.parametrize(
    ("method", "kwargs"),
    [
        ("DHD", {"subm": np.array([[0.0, 1.0], [1.0, 0.0]])}),
        ("OM", {"indel": 1.0, "subm": 2.0, "norm": "length"}),
        ("LCS", {"norm": "length"}),
        (
            "OMSPELL",
            {
                "indel": 1.0,
                "subm": 2.0,
                "expcost": 0.5,
                "tpow": 1.0,
                "norm": "maxlength",
            },
        ),
    ],
)
def test_pairwise_between_sets_is_the_rectangular_cross_block(method, kwargs):
    left = _CROSS_BLOCK_CORPUS[_CROSS_BLOCK_A]
    right = _CROSS_BLOCK_CORPUS[_CROSS_BLOCK_B]

    square = distances.pairwise(_CROSS_BLOCK_CORPUS.copy(), method=method, **kwargs)
    rectangular = distances.pairwise_between_sets(left, right, method=method, **kwargs)

    assert rectangular.shape == (len(left), len(right))
    np.testing.assert_array_equal(
        rectangular, np.asarray(square)[np.ix_(_CROSS_BLOCK_A, _CROSS_BLOCK_B)]
    )


def test_pairwise_between_sets_refuses_the_dataset_fitted_dhd_default():
    """The pinned ``subm`` above selects the pair-local DHD; the default is refused."""
    left = _CROSS_BLOCK_CORPUS[_CROSS_BLOCK_A]
    right = _CROSS_BLOCK_CORPUS[_CROSS_BLOCK_B]

    with pytest.raises(ValueError, match="estimated from the whole corpus"):
        distances.pairwise_between_sets(left, right, method="DHD")


def test_pairwise_normalization_and_deduplication_contracts():
    rows = np.array(
        [[0, 0, 1, 1], [0, 0, 1, 1], [0, 1, 1, 0], [1, 1, 0, 0]],
        dtype=np.int32,
    )

    raw = distances.pairwise(rows, method="LCS")
    normalized = distances.pairwise(rows, method="LCS", norm="length")
    np.testing.assert_allclose(normalized, raw / 8.0, rtol=0, atol=0)
    np.testing.assert_array_equal(
        distances.pairwise(rows, method="LCS", deduplicate=True),
        raw,
    )
    np.testing.assert_array_equal(
        distances.pairwise_between_sets(
            rows[:2],
            rows[2:],
            method="LCS",
            deduplicate=True,
        ),
        raw[:2, 2:],
    )


def test_pairwise_between_sets_ragged_rows_keep_rectangular_semantics():
    left = [np.array([0, 1, 2]), np.array([0, 1])]
    right = [np.array([0, 2]), np.array([2])]

    actual = distances.pairwise_between_sets(
        left,
        right,
        method="LCS",
        deduplicate=True,
    )

    np.testing.assert_array_equal(actual, np.array([[1.0, 2.0], [2.0, 3.0]]))
