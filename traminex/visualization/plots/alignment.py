"""Pairwise-alignment plot over the R-golden'd seqalign series (BUILD-from-R-design).

seqalignplot mirrors TraMineR's plot.seqalign: the optimal alignment of two
sequences drawn as three stacked coloured rows (seq1 / operations / seq2), one
column per alignment position (SOURCES/TraMineR/R/plot.seqalign.R). The operation
list / aligned sequences / per-op costs come from ``seqalign`` (R-golden'd in
``seqalign_*``). Original matplotlib; provenance docstring-only.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from ...align_p2 import seqalign
from ...sequences import _require_stslist
from ..palette import _state_color, ensure_cmap
from ..utils import _close_new_figures_on_error

# Operation -> colour: Equal / Substitution / Insertion / Deletion.
OP_COLORS = {"E": "#a6d854", "S": "#fc8d62", "I": "#8da0cb", "D": "#e78ac3"}
GAP_COLOR = "#dddddd"


@_close_new_figures_on_error
@_require_stslist
def seqalignplot(
    x, indices, *, indel=1.0, sm=None, include_missing=False, cpal=None, figsize=(8, 3), ax=None
):
    """Optimal alignment of two sequences -- re-expresses TraMineR ``plot.seqalign``.

    Three stacked rows, one column per alignment position:
    - top    = aligned ``seq1`` cells, coloured by state ('-' gaps grey)
    - middle = the operation at each position (E/S/I/D), coloured by op type
    - bottom = aligned ``seq2`` cells, coloured by state
    The operation list / aligned sequences come from ``seqalign`` (R-golden'd in
    ``seqalign_*``).

    Parameters
    ----------
    x : STSList
    indices : tuple[int, int]
        The zero-based rows of the two sequences to align (as :func:`seqalign`).
    indel, sm, include_missing : passed through to ``seqalign``.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
    """
    from matplotlib.patches import Rectangle

    # Both marker representations are equivalent inputs: with a declared
    # external nr code the marker text is not in ``x.alphabet``, so its cells
    # were painted in the gap grey (2026-09-22 audit S6-12).
    x = x._canonicalized()
    sa = seqalign(x, indices, indel=indel, sm=sm, include_missing=include_missing)
    ops = list(sa.operation)
    seq1 = list(sa.seq1)
    seq2 = list(sa.seq2)
    n = len(ops)
    alphabet = list(x.alphabet)
    cm = ensure_cmap(x, cpal)

    def _cell_color(symbol, *, gap: bool):
        if gap or symbol not in alphabet:
            return GAP_COLOR
        return _state_color(cm, alphabet.index(symbol))

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    else:
        fig = ax.figure

    # rows: y=2 seq1 (top), y=1 operations, y=0 seq2 (bottom)
    for i in range(n):
        ax.add_patch(
            Rectangle(
                (i, 2),
                1,
                1,
                facecolor=_cell_color(seq1[i], gap=ops[i] == "I"),
                edgecolor="white",
            )
        )
        ax.add_patch(
            Rectangle((i, 1), 1, 1, facecolor=OP_COLORS.get(ops[i], "#cccccc"), edgecolor="white")
        )
        ax.add_patch(
            Rectangle(
                (i, 0),
                1,
                1,
                facecolor=_cell_color(seq2[i], gap=ops[i] == "D"),
                edgecolor="white",
            )
        )
    ax.set_xlim(0, n)
    ax.set_ylim(0, 3)
    ax.set_yticks([0.5, 1.5, 2.5])
    ax.set_yticklabels(["seq2", "op", "seq1"])
    ax.set_xticks(np.arange(n) + 0.5)
    ax.set_xticklabels([str(i + 1) for i in range(n)])
    ax.set_xlabel("Alignment position")
    ax.set_title(f"Alignment (cost={sa.total_cost:g})")
    return fig, ax
