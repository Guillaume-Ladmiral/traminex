"""Raw-array INDELS/INDELSLOG costs agree with the R-validated object path.

``seqsubm`` on a plain integer matrix is a traminex extension; its docstring
promises R's ``indels <- 1/p`` with a non-finite value clamped to the largest
finite float. It floored a zero probability to ``1e-12`` first, so an absent
state cost ``1e12`` (and ``INDELSLOG`` was nudged off ``log 2``), where the
object path -- pinned to TraMineR by the ``seqcost`` goldens -- gives the
largest double (2026-09-22 audit S1a-10).
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

import traminex as tx

# State 1 never occurs, so its INDELS cost is 1/0.
X = np.array([[0, 2, 2, 0], [2, 0, 0, 2], [0, 0, 2, 2]])


@pytest.mark.parametrize("method", ["INDELS", "INDELSLOG"])
def test_an_absent_state_costs_what_the_object_path_gives(method):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        from_object = np.asarray(
            tx.seqcost(tx.seqdef(X, alphabet=[0, 1, 2]), method)["sm"], dtype=float
        )
        from_array = np.asarray(tx.seqsubm(X, method), dtype=float)

    np.testing.assert_array_equal(from_array, from_object)
