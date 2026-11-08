"""seqsubm's ``method`` default says what it does (TX-AUD-20261010-061).

R's ``seqsubm(...)`` is ``seqcost(...)$sm`` (seqcost.R:296-298) and ``seqcost``'s
``method`` has no default, so R stops without one. traminex refused too, but its
signature and stub showed ``method = "TRATE"``: a sentinel whose ``repr`` printed
``"TRATE"``, which only a raw integer array (a traminex extension) received. The
default is now ``None``: required for a sequence object, ``"TRATE"`` for a raw array.
"""

from __future__ import annotations

import inspect

import numpy as np
import pytest

import traminex as tx

SEQS = tx.seqdef(np.array([list("AABB"), list("ABAB")], dtype=object))
RAW = np.array([[0, 0, 1, 1], [0, 1, 0, 1]])


def test_the_method_default_is_none() -> None:
    assert inspect.signature(tx.seqsubm).parameters["method"].default is None


@pytest.mark.parametrize("kwargs", [{}, {"method": None}])
def test_a_sequence_object_requires_a_method(kwargs) -> None:
    with pytest.raises(TypeError, match="missing required argument: 'method'"):
        tx.seqsubm(SEQS, **kwargs)


def test_a_raw_array_defaults_to_trate() -> None:
    np.testing.assert_array_equal(tx.seqsubm(RAW), tx.seqsubm(RAW, method="TRATE"))
