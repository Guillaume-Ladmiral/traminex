from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pytest

import traminex as tx
from traminex.distances import om_trate_pairwise
from traminex.distances.costs import (
    future_costs,
    trate_substitution_matrix,
    trate_substitution_matrix_timevarying,
)
from traminex.sequences import STSList

_INVALID_WEIGHTS = (
    pytest.param(np.array([-1.0, 2.0]), "non-negative", id="negative"),
    pytest.param(np.array([np.nan, 1.0]), "finite", id="nan"),
    pytest.param(np.array([np.inf, 1.0]), "finite", id="infinite"),
)


@pytest.fixture
def sts() -> STSList:
    return tx.seqdef([["A", "B"], ["B", "A"]])


@pytest.mark.parametrize(("weights", "message"), _INVALID_WEIGHTS)
def test_stslist_rejects_invalid_weights_on_construction_and_assignment(
    weights: np.ndarray,
    message: str,
    sts: STSList,
) -> None:
    with pytest.raises(ValueError, match=message):
        STSList(np.array([[0], [1]], dtype=int), ["A", "B"], weights=weights)

    with pytest.raises(ValueError, match=message):
        sts.weights = weights


@pytest.mark.parametrize(
    ("data", "error", "message"),
    (
        pytest.param([[0]], TypeError, "data must be a numpy array", id="non-array"),
        pytest.param(
            np.array(0),
            ValueError,
            "data must be 2D",
            id="zero-dimensional",
        ),
    ),
)
def test_stslist_data_validation_precedes_weight_validation(
    data: object,
    error: type[Exception],
    message: str,
) -> None:
    with pytest.raises(error, match=message):
        STSList(data, ["A"], weights=np.array([1.0]))  # type: ignore[arg-type]


@pytest.mark.parametrize(("weights", "message"), _INVALID_WEIGHTS)
@pytest.mark.parametrize(
    "operation",
    (
        pytest.param(lambda x, w: tx.seqstatd(x, weights=w), id="seqstatd"),
        pytest.param(lambda x, w: tx.seqtab(x, idxs=0, weights=w), id="seqtab"),
    ),
)
def test_statistics_reject_invalid_explicit_weights(
    operation: Callable[[STSList, np.ndarray], object],
    weights: np.ndarray,
    message: str,
    sts: STSList,
) -> None:
    with pytest.raises(ValueError, match=message):
        operation(sts, weights)


@pytest.mark.parametrize(("weights", "message"), _INVALID_WEIGHTS)
@pytest.mark.parametrize(
    "operation",
    (
        pytest.param(
            lambda x, w: trate_substitution_matrix(x, weights=w),
            id="trate_substitution_matrix",
        ),
        pytest.param(lambda x, w: future_costs(x, weights=w), id="future_costs"),
        pytest.param(
            lambda x, w: trate_substitution_matrix_timevarying(x, weights=w),
            id="trate_substitution_matrix_timevarying",
        ),
        pytest.param(
            lambda x, w: om_trate_pairwise(x, trate_weights=w, norm=None),
            id="om_trate_pairwise",
        ),
        pytest.param(
            lambda x, w: tx.seqdist(
                x,
                method="OM",
                subm="TRATE",
                trate_weights=w,
                norm=None,
            ),
            id="seqdist",
        ),
    ),
)
def test_trate_paths_reject_invalid_explicit_weights(
    operation: Callable[[np.ndarray, np.ndarray], object],
    weights: np.ndarray,
    message: str,
) -> None:
    seqmat = np.array([[0, 1, 0], [1, 0, 1]], dtype=int)

    with pytest.raises(ValueError, match=message):
        operation(seqmat, weights)


def test_consumers_reject_in_place_corruption_of_attached_weights(sts: STSList) -> None:
    sts.weights = np.ones(sts.nseq)
    assert sts.weights is not None
    sts.weights[0] = -1.0

    with pytest.raises(ValueError, match="non-negative"):
        tx.seqstatd(sts)
    with pytest.raises(ValueError, match="non-negative"):
        tx.seqdist(sts, method="OM", subm="TRATE", norm=None)


def test_zero_sequence_weights_remain_supported(sts: STSList) -> None:
    zeros = np.zeros(sts.nseq)

    attached = tx.seqdef([["A", "B"], ["B", "A"]], weights=zeros)
    explicit = trate_substitution_matrix(sts.data, weights=zeros)

    np.testing.assert_array_equal(attached.weights, zeros)
    assert explicit.shape == (2, 2)


@pytest.mark.parametrize("dtype", [np.float64, np.int64])
def test_sequence_objects_keep_their_own_copy_of_the_weights(dtype: type) -> None:
    """Audit S2-13: a float64 weight array was stored by reference.

    ``np.asarray`` returned the caller's float64 array itself, so a later change
    to the caller's array silently changed the sequence object (an integer array
    was converted, hence copied). R copies on assignment; the object now always
    owns its weights.
    """
    caller = np.array([1, 2], dtype=dtype)

    defined = tx.seqdef([["A", "B"], ["B", "A"]], weights=caller)
    constructed = STSList(np.array([[0], [1]], dtype=int), ["A", "B"], weights=caller)
    assigned = tx.seqdef([["A", "B"], ["B", "A"]])
    assigned.weights = caller
    caller[0] = 99

    for x in (defined, constructed, assigned):
        assert x.weights is not None
        assert x.weights.tolist() == [1.0, 2.0]
        assert not np.shares_memory(x.weights, caller)
