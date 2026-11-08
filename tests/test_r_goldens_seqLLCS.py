"""R-golden parity for seqLLCS / seqLLCP (T12) vs TraMineR 2.2-14.

BREAKING re-point: the previous seqLLCS/seqLLCP returned an n x n distance
matrix; live R's seqLLCS(seq1, seq2) / seqLLCP(seq1, seq2) return the *scalar*
length of the longest common subsequence / prefix between two sequences,
compared over each sequence's non-void length (trailing void dropped, nr "*"
kept as a state). The LCS/LCP distance matrices now live on seqdist(method=...).
Goldens: scripts/r_goldens/make_goldens_seqLLCS.R (full pairwise matrices).
"""

from __future__ import annotations

import numpy as np
from corestats_fixtures import fixture_a, fixture_b, fixture_c
from r_goldens_utils import R_GOLDEN_ATOL, load_golden

from traminex import seqLLCP, seqLLCS


def _row(x, i):
    """Single-sequence STSList view of row i (mirrors R's seq[i,])."""
    r = x.copy()
    r.data = x.data[i : i + 1]
    return r


def _pairwise(fn, x):
    n = x.nseq
    out = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            out[i, j] = fn(_row(x, i), _row(x, j))
    return out


def _check(tag, fixture):
    x = fixture()
    np.testing.assert_allclose(
        _pairwise(seqLLCS, x), load_golden(f"seqLLCS_{tag}"), rtol=0, atol=R_GOLDEN_ATOL
    )
    np.testing.assert_allclose(
        _pairwise(seqLLCP, x), load_golden(f"seqLLCP_{tag}"), rtol=0, atol=R_GOLDEN_ATOL
    )


def test_seqLLCS_LLCP_a_voidfree():
    _check("A", fixture_a)


def test_seqLLCS_LLCP_b_interior_missing():
    _check("B", fixture_b)


def test_seqLLCS_LLCP_c_ragged_void():
    _check("C", fixture_c)


def test_seqLLCS_scalar_return_and_symmetry():
    x = fixture_a()
    v = seqLLCS(_row(x, 0), _row(x, 1))
    assert isinstance(v, int)
    assert seqLLCS(_row(x, 0), _row(x, 1)) == seqLLCS(_row(x, 1), _row(x, 0))
    # self-LCS equals the non-void length of the sequence
    assert seqLLCS(_row(x, 0), _row(x, 0)) == x.seqlen


def test_seqLLCS_rejects_mismatched_alphabet():
    import pytest

    x = fixture_a()
    y = x.copy()
    y.alphabet = list(x.alphabet) + ["Z"]
    with pytest.raises(ValueError, match="alphabet"):
        seqLLCS(_row(x, 0), _row(y, 0))
