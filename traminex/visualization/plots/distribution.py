from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from ...sequences import STSList, _real_state_codes, _require_stslist
from ...stats import seqstatd
from ..palette import ensure_cmap
from ..ticks import apply_time_axis
from ..utils import _close_new_figures_on_error, _legend_with_labels


@_close_new_figures_on_error
@_require_stslist
def state_distribution_plot(
    x: STSList,
    *,
    include_missing: bool = False,
    weights: np.ndarray | None = None,
    figsize=(8, 3),
    ax=None,
    mode: str = "d",  # 'd' (bars), 'Ht' (entropy line), 'dH' (bars + entropy)
    cmap=None,
    missing_alpha: float = 1.0,
    void_alpha: float = 0.4,
):
    """State distribution and entropy overlays.

    Parameters
    ----------
    x : STSList
        State sequence object.
    include_missing : bool, default=False
        Whether to include missing states in frequency computation.
    weights : np.ndarray | None, optional
        Per-sequence weights for frequency computation.
    figsize : tuple, default=(8, 3)
        Figure size.
    ax : matplotlib.axes.Axes, optional
        Axes to plot into. If None, creates new figure.
    mode : str, default="d"
        Plot mode:
        - 'd': per-position stacked bars of state proportions
        - 'Ht': entropy line only
        - 'dH': stacked bars with entropy line overlaid
    cmap : Matplotlib colormap or color sequence, optional
        State colors. Defaults to the palette stored on ``x`` when available.
    missing_alpha : float, default=1.0
        Alpha transparency for missing state bars (0=transparent, 1=opaque).
        The default draws them in the missing colour, opaque, as R does.
    void_alpha : float, default=0.4
        Alpha transparency for void state bars (0=transparent, 1=opaque).

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
        The figure and axes objects.
    """
    x = x._canonicalized()
    st = seqstatd(x, include_missing=include_missing, weights=weights)
    freq = st["Frequencies"]
    ent = st["Entropy"]
    L, K = freq.shape
    _ax_was_none = ax is None
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    xs = np.arange(1, L + 1)
    if mode in ("d", "dH"):
        bottom = np.zeros(L)
        colors = ensure_cmap(x, cmap)
        missing_code = x.nr_code
        void_code = x.void_code
        # R's legend (seqplot.R): state labels, "missing" only when missing
        # values are plotted, never void (2026-09-22 audit S6-11).
        names = x.labels if x.labels is not None else x.alphabet
        real_codes = set(_real_state_codes(x))
        show_missing = (
            include_missing
            and missing_code is not None
            and bool(np.any(np.asarray(x.data) == missing_code))
        )
        # Explicit handles: the automatic legend would drop a state whose label
        # starts with "_" (TX-REP-20261017-003).
        handles: list = []
        labels: list[str] = []
        for k in range(K):
            if k in real_codes:
                legend_label = str(names[k])
            elif k == missing_code and show_missing:
                legend_label = "missing"
            else:
                legend_label = None
            # Determine alpha for this state
            alpha = 1.0
            if missing_code is not None and k == missing_code:
                alpha = missing_alpha
            elif void_code is not None and k == void_code:
                alpha = void_alpha

            bars = ax.bar(
                xs,
                freq[:, k],
                bottom=bottom,
                width=0.8,
                label="_nolegend_" if legend_label is None else legend_label,
                color=colors(k),
                alpha=alpha,
            )
            if legend_label is not None:
                handles.append(bars)
                labels.append(legend_label)
            bottom += freq[:, k]
        ax.set_ylabel("Proportion")
        _legend_with_labels(ax, handles, labels, ncol=min(K, 6), fontsize="small")
    if mode in ("Ht", "dH"):
        ax.plot(xs, ent, color="black", linewidth=1.5, label="Entropy")
    apply_time_axis(ax, x, L, one_based=True)
    ax.set_xlabel("Time")
    if mode == "d":
        ax.set_title("State distribution plot")
    elif mode == "Ht":
        ax.set_title("Entropy over time (Ht)")
        ax.set_ylabel("Entropy")
    else:
        ax.set_title("State distribution with entropy (dH)")
    ax.set_xlim(0.5, L + 0.5)
    if _ax_was_none:
        plt.tight_layout()
    return fig, ax
