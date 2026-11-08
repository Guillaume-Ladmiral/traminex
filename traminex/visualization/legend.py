from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from ..sequences import _real_state_codes, _require_stslist
from .palette import _state_color, ensure_cmap
from .utils import _close_new_figures_on_error, _legend_with_labels, legend_figure


def _legend_entries(x, cmap=None) -> tuple[list[str], list]:
    """R ``seqlegend``'s entries: state labels and colours, then "missing".

    As in ``seqlegend.R`` (``with.missing = "auto"``), each real state is named
    by its label, a ``"missing"`` entry follows only when the data hold a
    missing value, and structural void is never listed (2026-09-22 audit
    S6-11). Colours come from the palette every renderer draws with, so a
    malformed ``cpal`` falls back exactly as in the plots (audit S6-14), and a
    plot drawn with a user ``cmap`` passes it, as R's legend takes the user
    ``cpal`` (TX-REP-20261021-002).
    """
    x = x._canonicalized()
    cmap = ensure_cmap(x, cmap)
    names = x.labels if x.labels is not None else x.alphabet
    codes = _real_state_codes(x)
    labels = [str(names[code]) for code in codes]
    colors = [_state_color(cmap, code) for code in codes]
    if x.nr_code is not None and np.any(np.asarray(x.data) == x.nr_code):
        labels.append("missing")
        colors.append(_state_color(cmap, x.nr_code))
    return labels, colors


@_close_new_figures_on_error
@_require_stslist
def state_legend(x, ax=None, *, ncol: int = 4):
    """Add a legend of the state colours to ``ax`` -- R ``seqlegend``.

    One entry per state, named by its label, plus ``"missing"`` when the data
    hold missing values; void is not a state and is not listed.

    Returns (handles, labels).
    """
    if ax is None:
        ax = plt.gca()
    return _draw_state_legend(x, ax, ncol=ncol)


def _draw_state_legend(x, ax, cmap=None, *, ncol: int = 4):
    """Draw :func:`state_legend`'s entries on ``ax`` in the palette ``cmap`` a plot used."""
    import matplotlib.patches as mpatches

    labels, colors = _legend_entries(x, cmap)
    patches = [mpatches.Patch(color=color, label=label) for label, color in zip(labels, colors)]
    _legend_with_labels(
        ax, patches, labels, ncol=min(ncol, len(patches)), fontsize="small", title="States"
    )
    return patches, labels


@_close_new_figures_on_error
@_require_stslist
def state_legend_figure(x, *, ncol: int = 5, figsize=(8, 1), fontsize: float | str = 10):
    """Create a standalone legend figure for the state's color palette and return (fig, ax).

    The entries are those of :func:`state_legend` (R ``seqlegend``).
    """
    labels, colors = _legend_entries(x)
    fig, ax = legend_figure(labels, colors, ncol=ncol, figsize=figsize, fontsize=fontsize)
    return fig, ax


__all__ = ["state_legend", "state_legend_figure", "legend_figure"]
