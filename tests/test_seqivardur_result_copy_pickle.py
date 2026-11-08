"""seqivardur's result copies and pickles (TX-AUD-20261010-053).

SeqIVarDurResult delegates unknown attributes to its value array through
``__getattr__``; copy and pickle look up attributes such as ``__setstate__``
on an instance whose ``_values`` is not set yet, so the delegation called
itself until RecursionError.
"""

from __future__ import annotations

import copy
import pickle

import numpy as np
import pytest

import traminex as tx

SEQS = tx.seqdef(np.array([list("ABCCAB"), list("ABAACC"), list("AAABBB")], dtype=object))


def _same(left, right) -> None:
    assert type(left) is type(right)
    np.testing.assert_array_equal(np.asarray(left), np.asarray(right))
    np.testing.assert_array_equal(np.asarray(left.vmax), np.asarray(right.vmax))
    np.testing.assert_array_equal(np.asarray(left.meand), np.asarray(right.meand))


@pytest.mark.parametrize("operation", ["copy", "deepcopy", "pickle"])
def test_the_result_survives(operation: str) -> None:
    result = tx.seqivardur(SEQS)
    clone = {
        "copy": copy.copy,
        "deepcopy": copy.deepcopy,
        "pickle": lambda value: pickle.loads(pickle.dumps(value)),
    }[operation](result)
    _same(clone, result)


def test_private_names_are_not_delegated() -> None:
    with pytest.raises(AttributeError):
        tx.seqivardur(SEQS)._no_such_attribute
