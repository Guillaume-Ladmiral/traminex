"""R-golden parity for seqalign (T17 VERIFY) vs TraMineR 2.2-14.

seqalign(seqdata, indices, indel, sm, with.missing) is the optimal-matching
alignment of two sequences: operation list (E/S/I/D), aligned seq1/seq2 (with
"-" for gaps) and per-operation cost (sum = OM distance). R aligns over
seqlength() positions -- trailing void dropped, nr "*" kept as a state. Goldens:
scripts/r_goldens/make_goldens_seqalign.R (4-line text golden per case:
operations / seq1c / seq2c / costs).
"""

from __future__ import annotations

import numpy as np
from corestats_fixtures import fixture_a, fixture_c
from r_goldens_utils import R_GOLDEN_ATOL, load_golden, load_golden_lines

from traminex import seqalign


def _const_sm(k: int) -> np.ndarray:
    sm = np.full((k, k), 2.0, dtype=float)
    np.fill_diagonal(sm, 0.0)
    return sm


def _check(name, x, indices, indel, sm, *, include_missing=False):
    lines = load_golden_lines(f"{name}")
    exp_op = lines[0].split()
    exp_seq1 = lines[1].split()
    exp_seq2 = lines[2].split()
    exp_cost = [float(t) for t in lines[3].split()]

    sa = seqalign(x, indices, indel=indel, sm=sm, include_missing=include_missing)
    assert list(sa.operation) == exp_op, name
    assert list(sa.seq1) == exp_seq1, name
    assert list(sa.seq2) == exp_seq2, name
    np.testing.assert_allclose(sa.cost, exp_cost, rtol=0, atol=R_GOLDEN_ATOL)
    np.testing.assert_allclose(sa.total_cost, sum(exp_cost), rtol=0, atol=R_GOLDEN_ATOL)
    _check_matrices(name, sa)


def _check_matrices(name, sa):
    """The two DP matrices are part of R's public return value (cycle-6 section 9).

    They were unpinned, which hid a difference at the origin cell: R fills the
    first column with "D" and *then* the first row with "I"
    (``SOURCES/TraMineR/R/seqalign.R:22-23``), so ``opmatrix[1, 1]`` ends up
    "I", while traminex explicitly assigned "" there.
    """
    lines = load_golden_lines(f"{name}_matrices")
    rows, cols = (int(token) for token in lines[0].split())
    expected_operations = [
        ["" if token == "." else token for token in line.split()] for line in lines[1 : 1 + rows]
    ]
    expected_costs = [
        [float(token) for token in line.split()] for line in lines[1 + rows : 1 + 2 * rows]
    ]

    operations = np.asarray(sa.operation_matrix)
    assert operations.shape == (rows, cols), name
    assert operations.tolist() == expected_operations, name
    np.testing.assert_allclose(
        np.asarray(sa.cost_matrix, dtype=float),
        np.asarray(expected_costs, dtype=float),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


def test_seqalign_clean_constant_cost():
    a = fixture_a()
    k = len(a.alphabet)
    _check("seqalign_A12_const", a, (0, 1), 1.0, _const_sm(k))
    _check("seqalign_A24_const", a, (1, 3), 1.0, _const_sm(k))


def test_seqalign_clean_trate_cost():
    a = fixture_a()
    sm = load_golden("seqalign_smA_trate")  # the exact R seqsubm(TRATE) matrix
    _check("seqalign_A12_trate", a, (0, 1), 1.0, sm)


def test_seqalign_ragged_void_with_missing():
    # Pins the void-stripping path: C is ragged (trailing void) with interior nr "*".
    # with.missing=TRUE -> nr is a state, sm is (K+1)x(K+1), trailing void dropped.
    c = fixture_c()
    # cost alphabet = real states (A,B,C) + nr; void is not a state.
    n_states = sum(1 for k in range(len(c.alphabet)) if k not in (c.nr_code, c.void_code))
    k = n_states + 1  # + nr (with.missing=TRUE)
    _check("seqalign_C13_const_wm", c, (0, 2), 1.0, _const_sm(k), include_missing=True)
    _check("seqalign_C12_const_wm", c, (0, 1), 1.0, _const_sm(k), include_missing=True)
