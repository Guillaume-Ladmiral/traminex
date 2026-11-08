"""Multidomain sequence plot over the per-domain raw series (BUILD-from-R-design).

seqplotMD mirrors TraMineR's seqplotMD: side-by-side per-domain panels, each a
``seqplot`` of one channel (SOURCES/TraMineR/R/seqplotMD.R). For ``type="i"`` each
panel is an index plot of that domain's state matrix; for ``type="d"`` a state
distribution. Original matplotlib (no Sequenzo borrow); provenance docstring-only.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from ...sequences import STSList
from ..series import resolve_sort_order
from ..utils import _close_new_figures_on_error
from .distribution import state_distribution_plot
from .index import sequence_index_plot


def _domain_items(value, *, name: str, count: int):
    if (
        isinstance(value, (str, bytes))
        or not hasattr(value, "__len__")
        or not hasattr(value, "__iter__")
    ):
        if name == "cpal":
            raise ValueError("cpal must be a collection with exactly one palette per domain")
        raise ValueError("dnames must contain exactly one name per domain")
    items = list(value)
    if len(items) != count:
        if name == "cpal":
            raise ValueError("cpal must be a collection with exactly one palette per domain")
        raise ValueError("dnames must contain exactly one name per domain")
    return items


@_close_new_figures_on_error
def seqplotMD(channels, *, type="i", cpal=None, dnames=None, figsize=(10, 4), **kwargs):
    """Multidomain sequence plot -- re-expresses TraMineR ``seqplotMD``.

    Parameters
    ----------
    channels : list of STSList
        One paired state-sequence object per domain (all with the same number of
        sequences). At least two are required.
    type : {"i", "d"}, default "i"
        Per-domain panel type: "i" = index plot, "d" = state distribution.
    cpal : collection, optional
        Ordered collection containing exactly one Matplotlib palette per domain.
    dnames : list of str, optional
        Domain titles (default "Dom1", "Dom2", ...).
    **kwargs
        Passed to every panel. With ``type="i"``, ``sortv="from.start"`` or
        ``"from.end"`` sorts the multidomain sequences
        (``seqMD(channels, what="MDseq")``) once and draws every panel in that
        order, as TraMineR's seqplotMD does with its default ``dom.crit=0``,
        so each row shows the same person in every domain.

    Returns
    -------
    fig, axes : matplotlib.figure.Figure, list[matplotlib.axes.Axes]
    """
    if not isinstance(channels, (list, tuple)) or len(channels) < 2:
        raise ValueError("channels should be a list of at least two state sequence objects")
    if any(not isinstance(channel, STSList) for channel in channels):
        raise TypeError("every channel must be an STSList state sequence object")
    ncase = channels[0].data.shape[0]
    for ch in channels:
        if ch.data.shape[0] != ncase:
            raise ValueError("domains with different number of sequences!")
    ndom = len(channels)
    if type not in ("i", "d"):
        raise ValueError(f"unsupported seqplotMD type {type!r} (use 'i' or 'd')")
    palettes = [None] * ndom if cpal is None else _domain_items(cpal, name="cpal", count=ndom)
    names = (
        [f"Dom{i + 1}" for i in range(ndom)]
        if dnames is None
        else _domain_items(dnames, name="dnames", count=ndom)
    )

    sortv = kwargs.get("sortv")
    if type == "i" and isinstance(sortv, str) and sortv in ("from.start", "from.end"):
        # seqplotMD.R:205-258: one order from the multidomain sequences,
        # handed to every panel as a rank vector.
        from ...multidomain_cat import seqMD

        order = resolve_sort_order(seqMD(list(channels), what="MDseq"), sortv)
        ranks = np.empty(ncase, dtype=float)
        ranks[order] = np.arange(ncase, dtype=float)
        kwargs = {**kwargs, "sortv": ranks}

    fig, axes = plt.subplots(1, ndom, figsize=figsize, constrained_layout=True)
    axes = list(axes) if ndom > 1 else [axes]
    for ax, ch, palette, name in zip(axes, channels, palettes, names):
        panel_kwargs = dict(kwargs)
        if palette is not None:
            panel_kwargs["cmap"] = palette
        if type == "i":
            sequence_index_plot(ch, ax=ax, **panel_kwargs)
        else:
            state_distribution_plot(ch, ax=ax, **panel_kwargs)
        ax.set_title(name)
    return fig, axes
