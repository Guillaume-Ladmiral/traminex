"""RP7a: seqalignplot (plot.seqalign) structure test.

Three stacked coloured rows (seq1 / operations / seq2), one column per alignment
position. The operation row's cell colours == the op-type palette for sa.operation,
the seq rows' cells == the state palette for sa.seq1/seq2. Anchored to the
seqalign_A12_const golden (4-line: ops / seq1 / seq2 / costs).
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
from corestats_fixtures import fixture_a
from matplotlib.colors import to_rgba
from r_goldens_utils import load_golden_lines

from traminex import STSList
from traminex.visualization.palette import ensure_cmap
from traminex.visualization.plots.alignment import GAP_COLOR, OP_COLORS, _state_color, seqalignplot


def _rows(ax):
    """Rectangles grouped by y0 (0=seq2, 1=op, 2=seq1), each ordered by x."""
    groups = defaultdict(list)
    for p in ax.patches:
        groups[round(float(p.get_y()), 6)].append(p)
    return {
        y0: [tuple(np.round(p.get_facecolor(), 6)) for p in sorted(g, key=lambda r: r.get_x())]
        for y0, g in groups.items()
    }


def test_seqalignplot_rows_match_seqalign_golden():
    x = fixture_a()
    k = len(x.alphabet)
    sm = np.full((k, k), 2.0)
    np.fill_diagonal(sm, 0.0)
    fig, ax = seqalignplot(x, (0, 1), indel=1.0, sm=sm)

    lines = load_golden_lines("seqalign_A12_const")
    exp_op, exp_seq1, exp_seq2 = lines[0].split(), lines[1].split(), lines[2].split()
    n = len(exp_op)
    cm = ensure_cmap(x, None)
    alphabet = list(x.alphabet)

    def cell(sym):
        return to_rgba(GAP_COLOR) if sym not in alphabet else _state_color(cm, alphabet.index(sym))

    rows = _rows(ax)
    assert len(rows[1.0]) == n
    # operation row colours == op-type palette
    assert rows[1.0] == [tuple(np.round(to_rgba(OP_COLORS[o]), 6)) for o in exp_op]
    # seq1 (top, y=2) and seq2 (bottom, y=0) rows == state palette
    assert rows[2.0] == [tuple(np.round(to_rgba(cell(s)), 6)) for s in exp_seq1]
    assert rows[0.0] == [tuple(np.round(to_rgba(cell(s)), 6)) for s in exp_seq2]


def test_alignment_distinguishes_real_dash_state_from_insertion_gap():
    data = np.array([[0, 2], [1, 0]], dtype=np.uint8)
    data.setflags(write=False)
    source = STSList(data, ["-", "A", "%"], void_code=2, cpal=["red", "blue", "white"])
    fig, ax = seqalignplot(source, (0, 1))
    try:
        rows = _rows(ax)

        def rounded(color):
            return tuple(np.round(to_rgba(color), 6))

        # Align '-' to 'A -': an insertion followed by equality, cost 1.
        assert rows[1.0] == [rounded(OP_COLORS["I"]), rounded(OP_COLORS["E"])]
        assert rows[2.0] == [rounded(GAP_COLOR), rounded("red")]
        assert rows[0.0] == [rounded("blue"), rounded("red")]
        assert ax.get_title() == "Alignment (cost=1)"
        np.testing.assert_array_equal(source.data, data)
    finally:
        plt.close(fig)


def test_alignment_distinguishes_real_dash_state_from_deletion_gap():
    data = np.array([[1, 0], [0, 2]], dtype=np.uint8)
    data.setflags(write=False)
    source = STSList(data, ["-", "A", "%"], void_code=2, cpal=["red", "blue", "white"])
    fig, ax = seqalignplot(source, (0, 1))
    try:
        rows = _rows(ax)

        def rounded(color):
            return tuple(np.round(to_rgba(color), 6))

        assert rows[1.0] == [rounded(OP_COLORS["D"]), rounded(OP_COLORS["E"])]
        assert rows[2.0] == [rounded("blue"), rounded("red")]
        assert rows[0.0] == [rounded(GAP_COLOR), rounded("red")]
        assert ax.get_title() == "Alignment (cost=1)"
        np.testing.assert_array_equal(source.data, data)
    finally:
        plt.close(fig)
