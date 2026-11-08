"""seqtab(include_missing=False) leaves out missing cells, not void padding (TX-AUD-20261010-027).

``include_missing`` is a traminex extension of seqtab (R's seqtab has no such
option). With ``include_missing=False`` it dropped every sequence holding a
void cell, i.e. every shorter sequence of a ragged object, so ex1 (all ragged,
no missing cell) gave an empty table. Void is padding, not an observation; only
sequences with a missing cell are left out.
"""

from __future__ import annotations

import numpy as np

import traminex as tx


def _patterns(x: tx.STSList, **kwargs) -> list[str]:
    return tx.seqtab(x, idxs=0, **kwargs)["sequence"].tolist()


def test_void_padded_sequences_are_kept() -> None:
    ragged = tx.seqdef(np.array([["A", "B", "B"], ["A", "B", None], ["A", "A", "A"]], dtype=object))
    assert sorted(_patterns(ragged, include_missing=False)) == ["A/1-B/1", "A/1-B/2", "A/3"]


def test_an_all_ragged_object_keeps_every_sequence() -> None:
    ragged = tx.seqdef(
        np.array([["A", None, None], ["B", "B", None], ["A", "B", None]], dtype=object)
    )
    assert sorted(_patterns(ragged, include_missing=False)) == ["A/1", "A/1-B/1", "B/2"]


def test_a_sequence_with_a_missing_cell_is_left_out() -> None:
    gapped = tx.seqdef(
        np.array([["A", None, "B"], ["A", "B", "B"], ["A", "A", None]], dtype=object)
    )
    assert sorted(_patterns(gapped)) == ["A/1-*/1-B/1", "A/1-B/2", "A/2"]
    assert sorted(_patterns(gapped, include_missing=False)) == ["A/1-B/2", "A/2"]
