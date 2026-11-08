from __future__ import annotations

import matplotlib.pyplot as plt

from ..._rformat import r_as_character, r_round
from ...sequences import _require_stslist
from ...stats import seqtab
from ...stats_basic import _slice_seqtab_result
from ..palette import _state_color, ensure_cmap
from ..ticks import _apply_time_axis_on_cells
from ..utils import _close_new_figures_on_error, _import_seaborn


@_close_new_figures_on_error
@_require_stslist
def seqfplot(
    x,
    top: int | None = 10,
    *,
    group=None,
    main="auto",
    cpal=None,
    include_missing: bool = True,
    border: bool = True,
    figsize=(8, 4),
    ax=None,
):
    """Sequence frequency plot -- re-expresses TraMineR ``seqfplot`` / ``plot.stslist.freq``.

    Draws the ``top`` most frequent distinct sequences as horizontal coloured
    state-segment bands: one band per distinct sequence, each time cell coloured by
    the state at that position, the band thickness proportional to the sequence's
    (weighted) frequency. This mirrors R's
    ``barplot(apply(x, 1, seqgbar), col=cpal, width=table$Percent, horiz=TRUE)``
    (SOURCES/TraMineR/R/plot.stslist.freq.R). The bands are stacked bottom-up in
    descending-frequency order (most frequent at the bottom, as in R). As R's
    default axis, the y axis is labelled 0% at the bottom and the cumulative
    percentage of the shown sequences at the top, titled "Cum. % freq. (n=N)"
    ("weighted n=N" for weighted sequences, N rounded to 2 decimals).

    The frequencies/order come straight from ``seqtab`` (R-golden'd); this function
    only reshapes them into coloured bands.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
    """
    if group is not None or main != "auto":
        from .dispatch import seqplot

        # Called directly the renderer draws no legend; main and group keep it so.
        return seqplot(
            x,
            type="f",
            with_legend=False,
            group=group,
            main=main,
            top=top,
            cpal=cpal,
            include_missing=include_missing,
            border=border,
            figsize=figsize,
            ax=ax,
        )
    from matplotlib.patches import Rectangle

    # Both marker representations are equivalent inputs (identity for
    # in-alphabet objects; cycle-9 audit N-04/N-05).
    x = x._canonicalized()
    tab = seqtab(x, idxs=0, normalize=True, include_missing=include_missing)
    if top is not None:
        tab = _slice_seqtab_result(tab, slice(None, top))
    cm = ensure_cmap(x, cpal)
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    else:
        fig = ax.figure
    L = x.data.shape[1]
    props = tab["prop"].to_numpy(dtype=float)
    representatives = tab.attrs["representatives"]

    y0 = 0.0
    for prop, codes in zip(props, representatives.data, strict=True):
        for t in range(L):
            ax.add_patch(
                Rectangle(
                    (t, y0),
                    1.0,
                    prop,
                    facecolor=_state_color(cm, int(codes[t])),
                    edgecolor="white" if border else "none",
                    linewidth=0.3 if border else 0.0,
                )
            )
        y0 += prop
    ax.set_xlim(0, L)
    ax.set_ylim(0, max(y0, 1e-9))
    # R's default y axis: 0% at the bottom of the first band and the cumulative
    # percentage of the shown sequences at the top of the last one
    # (plot.stslist.freq.R:56-61, 102-110; TX-AUD-20261011-067).
    shown = r_round(float(tab.attrs["freq"]["Percent"].sum()), 1)
    ax.set_yticks([0.0, y0])
    ax.set_yticklabels(["0%", f"{r_as_character(shown)}%"])
    _apply_time_axis_on_cells(ax, x, L)
    ax.set_xlabel("Time")
    weighted = "weighted " if tab.attrs["weighted"] else ""
    n = r_as_character(r_round(float(tab.attrs["nbseq"]), 2))
    ax.set_ylabel(f"Cum. % freq. ({weighted}n={n})")
    ax.set_title("Sequence frequency plot")
    return fig, ax


@_close_new_figures_on_error
def seqfrequency_plot(x, top: int | None = 10, *, seaborn: bool = True, figsize=(8, 4), ax=None):
    """Plot the most frequent sequences from seqtab."""
    tab = seqtab(x, idxs=0, normalize=True)
    if top is not None:
        tab = _slice_seqtab_result(tab, slice(None, top))
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    else:
        fig = ax.figure
    if seaborn:
        sns = _import_seaborn("seqfrequency_plot")
        sns.barplot(x="prop", y="sequence", data=tab, ax=ax, orient="h", color="#1f77b4")
        ax.set_xlabel("Proportion")
        ax.set_ylabel("Sequence")
        ax.set_title("Top sequence frequencies")
        return fig, ax
    ax.barh(tab["sequence"], tab["prop"], color="#1f77b4")
    ax.set_xlabel("Proportion")
    ax.set_ylabel("Sequence")
    ax.set_title("Top sequence frequencies")
    return fig, ax
