"""Event-duration plot ``seqedplot`` over the R-golden'd producing series (BUILD-from-R).

``seqedplot`` mirrors TraMineRextras' ``seqedplot`` (``seqedplot.R``): one coloured
curve per non-ignored event, of one of two kinds chosen by ``type``:

- ``type="survival"`` (default): the Kaplan-Meier survival step curve of the time
  to each event's first occurrence (non-occurrence censored at ``seqelength``);
- ``type="hazard"``: the weighted mean number of occurrences per equal-width age
  period.

The two series come from :mod:`traminex.events._duration`, which is pinned to live
R (``test_r_goldens_seqedplot.py``); this renderer is **structure-tested** against
those series (one curve per event, survival vertices, hazard abscissae, one panel
per group), not against pixel baselines (the RP7a rule). Original matplotlib (no
Sequenzo / ggseqplot borrow); provenance docstring-only.
"""

from __future__ import annotations

from typing import List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np

from ...events._duration import (
    _kept_event_codes,
    _seqedplot_ages_range,
    _seqedplot_follow_up,
    _seqedplot_hazard,
    _seqedplot_survival,
    _subset_seqelist,
)
from ..series import _factor_panels
from ..utils import _close_new_figures_on_error, _legend_with_labels


def _event_palette(nevent: int, cpal):
    """Per-event colours (R seqedplot: brewer Dark2 / Set3, else a long cycle).

    A given palette shorter than the events is recycled over them, as a short
    state palette is (docs/PLOTTING.md). R indexes ``cpal[event]``, so its later
    events get colour ``NA`` and are drawn in no colour (seqedplot.R:192-227;
    DIVERGENCE, TX-AUD-20261019-050).
    """
    if cpal is not None:
        # A string is one colour, as R's cpal = "red" (TX-REP-20261023-002).
        colours = [cpal] if isinstance(cpal, str) else list(cpal)
        if not colours:
            raise ValueError("cpal must contain at least one colour")
        return [colours[index % len(colours)] for index in range(nevent)]
    if nevent <= 8:
        cmap = plt.get_cmap("Dark2")
        return [cmap(index) for index in range(nevent)]
    if nevent <= 12:
        cmap = plt.get_cmap("Set3")
        return [cmap(index) for index in range(nevent)]
    cmap = plt.get_cmap("tab20")
    return [cmap(i % 20) for i in range(nevent)]


@_close_new_figures_on_error
def seqedplot(
    eseq,
    *,
    type: str = "survival",
    group: Optional[Sequence] = None,
    ignore: Optional[Sequence[str]] = None,
    breaks: int = 20,
    ages: Optional[Tuple[float, float]] = None,
    cpal=None,
    figsize=(8, 4),
    axes=None,
):
    """Event-duration plot -- re-expresses TraMineRextras ``seqedplot``.

    One coloured curve per non-ignored event. With ``type="survival"`` each curve
    is the unweighted Kaplan-Meier survival of the time to the event's first
    occurrence (a step starting at ``(0, 1)`` and dropping only at occurrence
    times); with ``type="hazard"`` each curve is the weighted mean number of
    occurrences per equal-width age period (a line over the period upper edges).
    ``group`` draws one panel per group level (sharing the hazard cut-points across
    panels); ``ignore`` drops events from the plot; ``ages`` clips the age window;
    ``breaks`` sets the number of hazard periods. The producing series are the
    live-R-golden'd :func:`_seqedplot_survival` / :func:`_seqedplot_hazard`.

    Parameters
    ----------
    eseq : SeqElist
        Event-sequence object (R ``seqelist``; see :func:`seqecreate`).
    type : {"survival", "hazard"}, default "survival"
    group : sequence, optional
        One value per event sequence; one panel per level of R's
        ``factor(group)`` -- numbers in numeric order, text in sorted order,
        categorical levels in their declared order. Sequences with a missing
        group value appear in no panel.
    ignore : sequence of str, optional
        Event labels to leave out (R ``ignore``).
    breaks : int, default 20
        Number of equal-width age periods for ``type="hazard"``.
    ages : (float, float), optional
        Age window; defaults to the global occurrence-age range.
    cpal : sequence, optional
        One colour per kept event (defaults to a qualitative palette); a single
        string is one colour. A shorter palette is recycled over the events, so
        each event keeps a visible colour; R draws the events beyond its palette
        in no colour.
    figsize : tuple, default (8, 4)
    axes : list of matplotlib.axes.Axes, optional
        Reuse these axes (exactly one per group panel) instead of creating a
        figure.

    Returns
    -------
    fig, axes : matplotlib.figure.Figure, list[matplotlib.axes.Axes]
        ``axes`` has one entry per group panel (length 1 when ``group`` is None).

    Raises
    ------
    TypeError
        If ``eseq`` is not a ``SeqElist``.
    ValueError
        If ``type`` is not "survival"/"hazard", if no event remains after
        ``ignore`` (R ``seqedplot.R:123``), if ``group`` does not hold one value
        per sequence, if ``axes`` does not hold one axis per panel, or -- for
        ``type="survival"`` -- if no observation time is set (all
        ``seqelength == -1``, R ``seqedplot.R:16-18``).
    """
    from ...events._object import SeqElist

    if not isinstance(eseq, SeqElist):
        raise TypeError("eseq must be a SeqElist event-sequence object")
    if type not in ("survival", "hazard"):
        raise ValueError("type should be one of 'survival' or 'hazard'.")

    kept_codes = _kept_event_codes(eseq, ignore)
    if len(kept_codes) < 1:
        raise ValueError("At least one event should remain aside the ignore list!")
    kept_labels = [eseq.dictionary[c - 1] for c in kept_codes]
    nevent = len(kept_codes)
    colors = _event_palette(nevent, cpal)

    # Group panels: R seqedplot.R:53 ``factor(group)`` -- numeric levels in
    # numeric order, missing values in no level (2026-09-22 audit S6-06).
    # Shared cut-points use the full-set range.
    if group is None:
        gindex: List[List[int]] = [list(range(len(eseq)))]
        titles: List[Optional[str]] = [None]
    else:
        if len(group) != len(eseq):
            raise ValueError("group must contain one value for each event sequence")
        gindex, level_titles = _factor_panels(group)
        titles = list(level_titles)
    nplot = len(gindex)

    effective_ages = ages if ages is not None else _seqedplot_ages_range(eseq, ignore)

    if axes is None:
        fig, axarr = plt.subplots(1, nplot, figsize=figsize, squeeze=False, constrained_layout=True)
        ax_list = list(axarr[0])
    else:
        ax_list = list(axes)
        if len(ax_list) != nplot:
            raise ValueError(f"axes must contain one axis for each of the {nplot} group panels")
        fig = ax_list[0].figure

    color_by_label = dict(zip(kept_labels, colors))
    # The first panel's curves, kept for the legend: explicit handles keep an
    # event whose label starts with "_" (batch-4 review M-4).
    legend_lines: list = []
    for ax, idx, title in zip(ax_list, gindex, titles):
        drawn = legend_lines if ax is ax_list[0] else []
        sub = _subset_seqelist(eseq, idx)
        if type == "survival":
            # R's plot.survfit runs each curve on to its last time, censored or
            # not, and stops at the right end of xlim = ages (seqedplot.R:217):
            # TX-AUD-20261011-026.
            follow_up = _seqedplot_follow_up(sub, ignore, ages)
            for label, km in _seqedplot_survival(sub, ignore, ages):
                times = km[:, 0]
                surv = km[:, 3]
                xs = np.concatenate([[0.0], times])
                ys = np.concatenate([[1.0], surv])
                last = min(follow_up[label], float(effective_ages[1]))
                if last > xs[-1]:
                    xs = np.append(xs, last)
                    ys = np.append(ys, ys[-1])
                drawn += ax.step(xs, ys, where="post", color=color_by_label[label], label=label)
            ax.set_ylim(0, 1.02)
            ax.set_ylabel("survival probability")
        else:
            onx, series = _seqedplot_hazard(sub, ignore, effective_ages, breaks)
            for label, ony in series:
                drawn += ax.plot(onx, ony, color=color_by_label[label], label=label)
            ax.set_ylabel("mean number of events")
        ax.set_xlabel("time")
        if title is not None:
            ax.set_title(title)

    if legend_lines:
        _legend_with_labels(
            fig,
            legend_lines,
            [line.get_label() for line in legend_lines],
            loc="upper right",
            fontsize="small",
            title="event",
        )
    return fig, ax_list
