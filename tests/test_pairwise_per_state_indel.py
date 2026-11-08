"""Low-level OM takes a per-state indel without void (TX-AUD-20261011-050).

The low-level layer documents "pairwise accepts a number or a per-state array of
numbers" (distances/_validation.py), and its void route and seqdist compute such
an indel, but the dense route (no void) converted it with ``float()`` and raised
``TypeError: only length-1 arrays can be converted to Python scalars``. The
expected values come from the optimal-matching recurrence written out here
(insert/delete a state s at cost indel[s], substitute at the constant cost).
"""

from __future__ import annotations

import numpy as np
import pytest

import traminex as tx
from traminex.distances import python_kernels as pk

CODES = np.array([[0, 1, 0], [1, 1, 0], [0, 0, 1]], dtype=np.int64)
INDEL = np.array([1.0, 2.0])
SUBM = 2.0


def _om(a: np.ndarray, b: np.ndarray) -> float:
    """Optimal matching by its recurrence, per-state indel costs."""
    table = np.zeros((len(a) + 1, len(b) + 1))
    for i in range(1, len(a) + 1):
        table[i, 0] = table[i - 1, 0] + INDEL[a[i - 1]]
    for j in range(1, len(b) + 1):
        table[0, j] = table[0, j - 1] + INDEL[b[j - 1]]
    for i in range(1, len(a) + 1):
        for j in range(1, len(b) + 1):
            substitution = 0.0 if a[i - 1] == b[j - 1] else SUBM
            table[i, j] = min(
                table[i - 1, j] + INDEL[a[i - 1]],
                table[i, j - 1] + INDEL[b[j - 1]],
                table[i - 1, j - 1] + substitution,
            )
    return float(table[-1, -1])


EXPECTED = np.array([[_om(a, b) for b in CODES] for a in CODES])


@pytest.mark.parametrize(
    "pairwise", [pk.pairwise, tx.distances.pairwise], ids=["python_kernels", "dispatcher"]
)
def test_the_dense_route_takes_a_per_state_indel(pairwise) -> None:
    np.testing.assert_allclose(
        pairwise(CODES, method="OM", indel=INDEL, subm=SUBM), EXPECTED, atol=1e-12
    )


def test_the_void_route_and_seqdist_agree() -> None:
    with_void = np.array([[0, 1, 0, 2], [1, 1, 0, 2], [0, 0, 1, 2]], dtype=np.int64)
    np.testing.assert_allclose(
        tx.distances.pairwise(with_void, method="OM", indel=INDEL, subm=SUBM, void=2),
        EXPECTED,
        atol=1e-12,
    )
    seqs = tx.seqdef(np.array([list("ABA"), list("BBA"), list("AAB")], dtype=object))
    np.testing.assert_allclose(
        tx.seqdist(seqs, method="OM", indel=list(INDEL), subm=SUBM), EXPECTED, atol=1e-12
    )


def test_a_vector_of_the_wrong_length_is_refused() -> None:
    with pytest.raises(ValueError, match="indel"):
        tx.distances.pairwise(CODES, method="OM", indel=np.array([1.0, 2.0, 3.0]), subm=SUBM)
