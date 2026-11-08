"""Per-sequence indicators return a 1-D array (TX-AUD-20261010-029).

docs/API_CONSISTENCY.md section 1c records traminex's convention: R returns
most per-sequence indicators as an ``n x 1`` matrix, traminex as a length-n
array. ``seqsubsn`` was the exception, returning shape ``(n, 1)`` (R's own
``seqsubsn`` is ``n x 1``, column "Subseq.").
"""

from __future__ import annotations

import numpy as np
import pytest

import traminex as tx

SEQS = tx.seqdef(
    np.array([list("AABCA"), list("ABBCC"), list("CCABA"), list("BBBBA")], dtype=object),
    alphabet=["A", "B", "C"],
)
INDICATORS = ["seqient", "seqtransn", "seqlength", "seqici", "seqST", "seqsubsn", "seqivolatility"]


@pytest.mark.parametrize("name", INDICATORS)
def test_a_per_sequence_indicator_is_one_dimensional(name: str) -> None:
    result = np.asarray(getattr(tx, name)(SEQS))
    assert result.shape == (SEQS.nseq,)
