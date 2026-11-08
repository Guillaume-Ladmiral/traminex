"""Canonicalizing keeps a real state spelled like an external marker (TX-AUD-20261010-034).

``STSList._canonicalized`` folds an external nr or void code (outside the
alphabet) into the alphabet. When the marker's text already named an alphabet
slot it merged into that slot, even when that slot was a real state: with
alphabet ``["A", "*"]`` and an external nr ``"*"``, the real ``"*"`` cells then
read as missing (lengths [2, 3] became [1, 1] after an identity
``seqgen_missing(p_cases=0)``). traminex's marker model gives a marker a
collision-free spelling when a real state uses its text (as seqdef does).
"""

from __future__ import annotations

import numpy as np

import traminex as tx
from traminex.sequences import STSList


def _object() -> STSList:
    # rows: A, real *, nr | real *, real *, A
    return STSList(data=np.array([[0, 1, 2], [1, 1, 0]]), alphabet=["A", "*"], nr="*", nr_code=2)


def test_identity_operations_keep_the_real_star_cells() -> None:
    x = _object()
    assert list(tx.seqlength(x, with_missing=False)) == [2, 3]
    assert list(tx.seqlength(tx.seqgen_missing(x, p_cases=0), with_missing=False)) == [2, 3]


def test_canonical_form_keeps_real_and_missing_cells_apart() -> None:
    canonical = _object()._canonicalized()
    assert canonical.data.tolist()[1] == [1, 1, 0]  # the real "*" cells keep code 1
    nr_code = int(canonical.nr_code)
    assert nr_code != 1 and canonical.data[0, 2] == nr_code
    assert canonical.alphabet[nr_code] == canonical.nr != "*"


def test_an_unused_slot_named_like_the_marker_is_still_merged() -> None:
    # No cell uses slot 1, so the external nr folds into it (the group_compare case).
    x = STSList(data=np.array([[0, 2], [0, 0]]), alphabet=["A", "*"], nr="*", nr_code=2)
    canonical = x._canonicalized()
    assert canonical.alphabet == ["A", "*"]
    assert canonical.nr_code == 1
    assert canonical.data.tolist() == [[0, 1], [0, 0]]
