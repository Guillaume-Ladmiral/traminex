"""A NumPy scalar substitution cost behaves as the number it is (TX-AUD-20261019-028).

``isinstance(subm, (int, float))`` missed NumPy scalars, which then fell into the matrix path: ``np.float32(2)``
and ``np.int64(2)`` raised IndexError ("too many indices for array: array is 0-dimensional") on 14 of 16 method
and data-shape cases. A NumPy scalar or a 0-d array now gives the matrix of the Python float; a bool is still
refused.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

import traminex as tx

RAGGED = tx.seqdef(
    [
        list("ABC") + [None],
        list("BA") + [None, None],
        list("ABC") + [None],
        list("CCBA"),
        list("BA") + [None, None],
    ],
    alphabet=list("ABC"),
    right="DEL",
)
DENSE = tx.seqdef(
    [list("ABC"), list("BAA"), list("ABC"), list("CCB"), list("BAA")],
    alphabet=list("ABC"),
)
METHODS = {
    "OM": dict(indel=1.0),
    "OMslen": dict(indel=1.0),
    "TWED": dict(nu=0.5, h=0.5),
    "OMspell": dict(indel=1.0),
}


def _quiet(*args, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return np.asarray(tx.seqdist(*args, **kwargs), dtype=float)


@pytest.mark.parametrize("data", ["ragged", "dense"])
@pytest.mark.parametrize("method", sorted(METHODS))
@pytest.mark.parametrize(
    "scalar", [np.float32(2), np.int64(2), np.array(2.0)], ids=["float32", "int64", "0d"]
)
def test_numpy_scalar_cost_equals_the_python_float(scalar, method: str, data: str) -> None:
    x = RAGGED if data == "ragged" else DENSE
    np.testing.assert_array_equal(
        _quiet(x, method=method, subm=scalar, **METHODS[method]),
        _quiet(x, method=method, subm=2.0, **METHODS[method]),
    )


def test_a_boolean_cost_is_still_refused() -> None:
    with pytest.raises((TypeError, ValueError)):
        _quiet(DENSE, method="OM", subm=True, indel=1.0)
