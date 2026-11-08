from __future__ import annotations

from typing import Literal, Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

from ...sequences import STSList, _require_stslist
from ..legend import _draw_state_legend
from ..palette import _state_code_norm, ensure_cmap
from ..series import select_index_rows
from ..ticks import apply_time_axis
from ..utils import _close_new_figures_on_error, compute_even_ticks

# Row labels drawn by sequence_index_plot at most (the first and last included).
_MAX_ROW_LABELS = 50


@_close_new_figures_on_error
@_require_stslist
def sequence_index_plot(
    x: STSList,
    *,
    cmap: Optional[Sequence[str]] = None,
    figsize=(8, 4),
    highlight: Optional[Sequence[int]] = None,
    legend: bool = False,
    ax=None,
    sortv=None,
    sampling: int | None = None,
    sample_meth: Literal["systematic", "random"] = "systematic",
    weighted: bool = True,
    random_state: int | np.random.Generator | None = None,
    missing_alpha: float = 1.0,
    void_alpha: float = 0.2,
    apply_hatching: bool = False,
):
    """Plot sequence index plot (basic). Optionally highlight given sequence indices and add legend.

    The y axis labels each displayed row with its identifier (or 1-based row
    number); with more than 50 rows, 50 evenly spaced rows are labelled,
    including the first and the last.

    Parameters
    ----------
    x : STSList
        State sequence object.
    cmap : Optional[Sequence[str]]
        Custom color map for states.
    figsize : tuple, default=(8, 4)
        Figure size.
    highlight : Optional[Sequence[int]]
        Sequence indices to highlight with red borders (zero-based positions in
        the *plotted* order, i.e. after any ``sortv`` reordering).
    legend : bool, default=False
        Whether to add a legend.
    ax : matplotlib.axes.Axes, optional
        Axes to plot into. If None, creates new figure.
    sortv : None | "from.start" | "from.end" | array-like, optional
        Row sort order, faithful to TraMineR ``plot.stslist``/``seqIplot`` ``sortv``
        ("from.start"/"from.end" lexicographic, or one value per sequence). None
        keeps the input order.
    sampling : int, optional
        Maximum number of displayed sequences. ``None`` displays all rows.
    sample_meth : {"systematic", "random"}, default="systematic"
        Row selection method, applied after sorting; a unique prefix (``"sys"``)
        is accepted, as R's ``match.arg`` takes it.
    weighted : bool, default=True
        Draw each row with a height equal to its sequence's weight, as R's
        ``barplot(width = weights)`` (the rows touch, from 0 up to the total
        weight), and use the weights as random-selection probabilities. Rows of
        equal weight, or ``weighted=False``, are drawn with equal heights.
    random_state : int or numpy.random.Generator, optional
        Local random stream for reproducible random sampling.
    missing_alpha : float, default=1.0
        Alpha transparency for missing states (0=transparent, 1=opaque). The
        default draws them in the missing colour, opaque, as R does.
    void_alpha : float, default=0.2
        Alpha transparency for void states (0=transparent, 1=opaque).
    apply_hatching : bool, default=False
        Whether to apply hatching pattern to missing/void states.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
        The figure and axes objects.
    """
    # Both marker representations are equivalent inputs (identity for
    # in-alphabet objects; cycle-9 audit N-04/N-05).
    x = x._canonicalized()
    selection = select_index_rows(
        x,
        sortv=sortv,
        sampling=sampling,
        sample_meth=sample_meth,
        weighted=weighted,
        random_state=random_state,
    )
    data = selection.display_rows
    n, L = data.shape
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    else:
        fig = ax.figure

    cm = ensure_cmap(x, cmap)
    im = None
    # R draws each row with a height equal to its case's weight: barplot(width =
    # weights[idxs[sel]], horiz = TRUE) with space 0, the first row at the bottom
    # (plot.stslist.R:149; TX-AUD-20261010-059, TX-AUD-20261019-048). Rows of equal
    # weight keep the image route, centred on 0, 1, ..., n - 1.
    heights = np.asarray(selection.weights, dtype=float)
    weighted_rows = n > 0 and L > 0 and not bool(np.all(heights == 1.0))
    if weighted_rows:
        y_edges = np.concatenate(([0.0], np.cumsum(heights)))
    else:
        y_edges = np.arange(n + 1, dtype=float) - 0.5
    x_edges = np.arange(L + 1, dtype=float) - 0.5
    if n > 0 and L > 0:
        if weighted_rows:
            im = ax.pcolormesh(
                x_edges, y_edges, data, cmap=cm, norm=_state_code_norm(x), shading="flat"
            )
            ax.set_xlim(x_edges[0], x_edges[-1])
            ax.set_ylim(y_edges[0], y_edges[-1])
        else:
            im = ax.imshow(
                data,
                aspect="auto",
                interpolation="nearest",
                cmap=cm,
                norm=_state_code_norm(x),
                origin="lower",
            )
    else:
        ax.set_xlim(-0.5, max(L - 0.5, 0.5))
        ax.set_ylim(-0.5, max(n - 0.5, 0.5))

    def white_overlay(mask: np.ndarray, alpha: float) -> None:
        # A white layer of opacity 1 - alpha over the masked cells.
        if weighted_rows:
            ax.pcolormesh(
                x_edges,
                y_edges,
                np.ma.masked_where(~mask, np.ones(mask.shape)),
                cmap=ListedColormap(["white"]),
                vmin=0.0,
                vmax=1.0,
                alpha=1.0 - alpha,
                shading="flat",
            )
            return
        overlay = np.ones((n, L, 4))
        overlay[:, :, 3] = np.where(mask, 1.0 - alpha, 0.0)
        ax.imshow(
            overlay,
            aspect="auto",
            interpolation="nearest",
            extent=(x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]),
            origin="lower",
        )

    def hatch(mask: np.ndarray, edgecolor: str, pattern: str) -> None:
        for i, j in zip(*np.nonzero(mask)):
            ax.add_patch(
                plt.Rectangle(
                    (x_edges[j], y_edges[i]),
                    1,
                    y_edges[i + 1] - y_edges[i],
                    fill=False,
                    edgecolor=edgecolor,
                    hatch=pattern,
                    linewidth=0,
                )
            )

    # Apply missing/void styling if present
    for code, alpha, edgecolor, pattern in (
        (x.nr_code, missing_alpha, "gray", "///"),
        (x.void_code, void_alpha, "lightgray", "..."),
    ):
        if im is None or code is None:
            continue
        mask = data == code
        if np.any(mask):
            white_overlay(mask, alpha)
            if apply_hatching:
                hatch(mask, edgecolor, pattern)

    apply_time_axis(ax, x, L)
    # At most _MAX_ROW_LABELS evenly spaced row labels: laying out one label
    # per row made large plots take tens of seconds, and R's axis() drops
    # overlapping labels anyway (2026-09-22 audit S6-16).
    rows = compute_even_ticks(n, max_labels=_MAX_ROW_LABELS)
    centres = (y_edges[:-1] + y_edges[1:]) / 2.0
    ax.set_yticks([float(centres[row]) for row in rows])
    ax.set_yticklabels([str(selection.labels[row]) for row in rows])
    if highlight:
        for idx in highlight:
            if 0 <= idx < n:
                ax.axhline(y_edges[idx], color="red", linewidth=1.5)
                ax.axhline(y_edges[idx + 1], color="red", linewidth=1.5)
    ax.set_xlabel("Time")
    ax.set_ylabel("Sequence index")
    ax.set_title("Sequence index plot")
    if legend:
        _draw_state_legend(x, ax, cmap)
    return fig, ax


@_close_new_figures_on_error
@_require_stslist
def sequence_index_plot_grouped(
    x: STSList,
    labels: np.ndarray,
    *,
    order: Optional[Sequence[int] | np.ndarray] = None,
    figsize=(8, 4),
    cmap: Optional[Sequence[str]] = None,
    show_separators: bool = True,
    legend: bool = True,
):
    """Sequence index plot with sequences grouped by cluster labels.

    ``labels`` holds one label per sequence; ``order``, if given, is a
    permutation of the rows (each of 0 to number of sequences - 1 once),
    else the rows are grouped by label. Either rule broken raises ValueError.

    Returns (fig, ax, new_order) where new_order is the permutation applied to rows.
    """
    # Both marker representations are equivalent inputs (identity for
    # in-alphabet objects; cycle-9 audit N-04/N-05).
    x = x._canonicalized()
    labels = np.asarray(labels)
    if labels.shape[0] != x.nseq:
        raise ValueError("labels length must equal number of sequences")
    # determine order: group by label, then by original index unless a custom order is provided
    if order is None:
        resolved_order = np.argsort(labels, kind="stable")
    else:
        # A permutation of the rows, so each sequence is drawn once (TX-AUD-20261011-029).
        resolved_order = np.asarray(order)
        if resolved_order.dtype.kind not in "iu" or not np.array_equal(
            np.sort(resolved_order), np.arange(x.nseq)
        ):
            raise ValueError(
                "order must be a permutation of the sequence rows (0 to number of sequences - 1)"
            )
    data = x.data[resolved_order]
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    cm = ensure_cmap(x, cmap)
    ax.imshow(
        data,
        aspect="auto",
        interpolation="nearest",
        cmap=cm,
        norm=_state_code_norm(x),
        origin="lower",
    )
    apply_time_axis(ax, x, data.shape[1])
    if show_separators:
        # draw lines between clusters
        sorted_labels = labels[resolved_order]
        # Use equality comparison on adjacent elements to detect boundaries; works for numeric and string labels
        change_idx = np.where(sorted_labels[1:] != sorted_labels[:-1])[0]
        for i in change_idx:
            ax.axhline(i + 0.5, color="white", linewidth=1.0)
    ax.set_xlabel("Position")
    ax.set_ylabel("Sequence index (grouped)")
    ax.set_title("Sequence index plot (grouped by cluster)")
    if legend:
        _draw_state_legend(x, ax, cmap)
    return fig, ax, resolved_order.copy()


@_close_new_figures_on_error
@_require_stslist
def overlay_representatives_on_index_plot(
    x: STSList, labels: np.ndarray, D: Optional[np.ndarray] = None, *, figsize=(8, 4)
):
    """Convenience: draw grouped index plot and highlight cluster medoids if D provided."""
    from ...cluster import medoid_index

    labels = np.asarray(labels)
    order = np.argsort(labels, kind="stable")
    fig, ax, order = sequence_index_plot_grouped(x, labels, order=order, figsize=figsize)
    if D is not None:
        for c in np.unique(labels):
            idx = np.where(labels == c)[0]
            if idx.size == 0:
                continue
            sub = D[np.ix_(idx, idx)]
            mlocal = medoid_index(sub)
            med_global = int(idx[mlocal])
            # find the row in the reordered matrix
            pos = int(np.where(order == med_global)[0][0])
            ax.axhline(pos - 0.5, color="red", linewidth=1.5)
            ax.axhline(pos + 0.5, color="red", linewidth=1.5)
    return fig, ax
