"""The low-level ``pairwise``'s legacy ``norm="length"`` keeps one meaning.

``docs/DISTANCES_NORMALIZATION.md`` defines the low-level vocabulary: its
``'length'`` divides LCS/LCP/RLCP distances by ``l1 + l2``, and for OM divides
the raw distance by the longer length. The dense path does that; the
variable-length row path -- reached as soon as any row contains the void
symbol, and always by ``pairwise_between_sets`` with a short row -- mapped
``'length'`` to TraMineR's ``maxlength`` instead. So the same equal-length pair
got LCS ``0.667`` or ``1.333``, LCP/RLCP ``1`` or ``2``, and OM with
``indel=2`` ``1.333`` or ``0.667`` depending on an unrelated row (2026-09-22
audit S1b-07). ``seqdist`` resolves the aliases first and is not affected.
"""

from __future__ import annotations

import numpy as np
import pytest

from traminex.distances import pairwise, pairwise_between_sets

DENSE = np.array([[0, 1, 2], [1, 1, 0], [2, 2, 1]])
VOID = 9
WITH_SHORT_ROW = np.vstack([DENSE, [[0, 1, VOID]]])
CASES = [
    ("LCS", {}),
    ("LCP", {}),
    ("RLCP", {}),
    ("OM", {"subm": 2.0, "indel": 2.0}),
]


@pytest.mark.parametrize("method,extra", CASES, ids=[m for m, _ in CASES])
def test_a_void_row_does_not_change_the_meaning_of_length(method, extra):
    dense = pairwise(DENSE, method=method, norm="length", **extra)
    with_void = pairwise(WITH_SHORT_ROW, method=method, norm="length", void=VOID, **extra)

    np.testing.assert_allclose(with_void[:3, :3], dense, rtol=0, atol=1e-12)


@pytest.mark.parametrize("method,extra", CASES, ids=[m for m, _ in CASES])
def test_the_ragged_values_follow_the_documented_denominators(method, extra):
    raw = pairwise(WITH_SHORT_ROW, method=method, void=VOID, **extra)
    normed = pairwise(WITH_SHORT_ROW, method=method, norm="length", void=VOID, **extra)
    lengths = np.array([3.0, 3.0, 3.0, 2.0])
    if method == "OM":
        denominator = np.maximum.outer(lengths, lengths)
    else:
        denominator = np.add.outer(lengths, lengths)

    np.testing.assert_allclose(normed, raw / denominator, rtol=0, atol=1e-12)


@pytest.mark.parametrize("method,extra", CASES, ids=[m for m, _ in CASES])
def test_pairwise_between_sets_keeps_the_same_meaning(method, extra):
    dense_cross = pairwise_between_sets(DENSE[:2], DENSE[2:], method=method, norm="length", **extra)
    ragged_cross = pairwise_between_sets(
        DENSE[:2], [DENSE[2], np.array([0, 1])], method=method, norm="length", **extra
    )

    np.testing.assert_allclose(ragged_cross[:, :1], dense_cross, rtol=0, atol=1e-12)
