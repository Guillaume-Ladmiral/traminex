"""Event-sequence plots over the R-golden'd event-mining series (BUILD-from-R-design).

These render the outputs of traminex's faithful event-mining engine (seqefsub ->
SubSeqElist, seqecmpgroup -> SubSeqElistChisq; events.py / events/_mining.py) the
way TraMineR's plot.subseqelist / plot.subseqelistchisq do (SOURCES/TraMineR/R/
subseqelist.R, seqecmpgroup.R), plus the parallel-coordinate event plot (seqpcplot)
over an Eseq/SeqElist. Original matplotlib; provenance docstring-only.
"""

from __future__ import annotations

import warnings
from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import is_color_like, to_rgba

from ..._deprecated import _caller_stacklevel
from ...events._mining import SubSeqElist, SubSeqElistChisq
from ...events._object import Eseq, SeqElist
from ...sequences import STSList
from ..utils import _close_new_figures_on_error


def _trajectory_colors(cpal):
    if cpal is None:
        return None
    if is_color_like(cpal):
        values = [cpal]
    else:
        try:
            values = list(cpal)
        except TypeError:
            values = [cpal]
    if not values:
        raise ValueError("cpal must contain at least one Matplotlib color")
    return [to_rgba(color) for color in values]


@_close_new_figures_on_error
def plot_subseqelist(x, *, freq: str = "Support", top=None, figsize=(8, 4), ax=None):
    """Frequent-subsequence bars -- re-expresses TraMineR ``plot.subseqelist``.

    One horizontal bar per mined subsequence, length = its ``Support`` (or
    ``Count``); the subsequence labels / frequencies come from ``seqefsub``'s
    ``SubSeqElist`` (R-golden'd in ``seqefsub_*``). Most frequent on top.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
    """
    if not isinstance(x, SubSeqElist):
        raise TypeError("x must be a SubSeqElist event-subsequence object")
    labels = list(x.subseq)
    vals = x.data[freq].to_numpy(dtype=float)
    if top is not None:
        labels, vals = labels[:top], vals[:top]
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    else:
        fig = ax.figure
    ax.barh(np.arange(len(labels)), vals, color="#4C72B0")
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel(freq)
    ax.set_title("Frequent subsequences")
    return fig, ax


@_close_new_figures_on_error
def plot_subseqelistchisq(x, *, top=None, figsize=(10, 4)):
    """Discriminating-subsequence bars -- re-expresses TraMineR ``plot.subseqelistchisq``.

    One panel per group, titled by the group label: a horizontal bar per
    subsequence, length = that group's relative frequency ``Freq.<group>``,
    coloured as R colours it: by the bin of the group's Pearson residual
    ``Resid.<group>`` cut at ``+/- -qnorm(c(0.05, 0.01))`` (Bonferroni-adjusted
    when seqecmpgroup used it), with the five Brewer "RdBu" colours, so a positive
    residual is blue (seqecmpgroup.R:105-142; TX-AUD-20261010-045). Every panel
    has the same frequency range, from ``min(0, smallest frequency)`` to the
    largest frequency of any group over the shown rows (the first ``top``), as
    R's ``ylim = "uniform"`` on ``x[1:top]``. The
    frequencies / residuals / group labels come from ``seqecmpgroup``'s
    ``SubSeqElistChisq`` (R-golden'd in ``seqecmpgroup_*``).

    Returns
    -------
    fig, axes : matplotlib.figure.Figure, list[matplotlib.axes.Axes]
    """
    if not isinstance(x, SubSeqElistChisq):
        raise TypeError("x must be a SubSeqElistChisq event-subsequence object")

    from scipy.stats import norm

    labels = list(x.subseq)
    groups = list(x.labels)
    levels = np.asarray([0.05, 0.01])
    bonferroni = getattr(x, "bonferroni", None) or {}
    if bonferroni.get("used"):
        levels = 1.0 - (1.0 - levels) ** (1.0 / float(bonferroni["ntest"]))
    bins = np.sort(-norm.ppf(levels))
    n = len(labels) if top is None else min(top, len(labels))
    # One frequency range for every panel, from all groups and the shown rows: R's
    # ylim = "uniform" (seqecmpgroup.R:136-138; TX-AUD-20261019-049) on x[1:top],
    # whose data are those rows (subseqelist.R:28-32; batch-13 review M-3).
    all_freq = x.data[[f"Freq.{g}" for g in groups]].to_numpy(dtype=float)[:n]
    fig, axes = plt.subplots(1, len(groups), figsize=figsize, sharey=True, constrained_layout=True)
    axes = list(axes) if len(groups) > 1 else [axes]
    for ax, g in zip(axes, groups):
        freq = x.data[f"Freq.{g}"].to_numpy(dtype=float)[:n]
        resid = x.data[f"Resid.{g}"].to_numpy(dtype=float)[:n]
        colors = [_residual_colour(value, bins) for value in resid]
        ax.barh(np.arange(n), freq, color=colors)
        ax.set_yticks(np.arange(n))
        ax.set_yticklabels(labels[:n])
        ax.invert_yaxis()
        if all_freq.size:
            ax.set_xlim(min(float(all_freq.min()), 0.0), float(all_freq.max()))
        ax.set_title(str(g))
        ax.set_xlabel("Frequency")
    return fig, axes


# brewer.pal(5, "RdBu"): strong negative ... strong positive.
_RESIDUAL_PALETTE = ("#CA0020", "#F4A582", "#F7F7F7", "#92C5DE", "#0571B0")


def _residual_colour(value: float, bins: np.ndarray) -> str:
    """R's ``cut(residual, c(-Inf, -l2, -l1, l1, l2, Inf), labels = cpal)``.

    The intervals are right-closed, so ``-Inf`` and ``NaN`` fall in none (R's
    ``NA``, a bar without fill).
    """
    low, high = float(bins[0]), float(bins[1])
    if np.isnan(value) or value == -np.inf:
        return "none"
    for index, upper in enumerate((-high, -low, low, high)):
        if value <= upper:
            return _RESIDUAL_PALETTE[index]
    return _RESIDUAL_PALETTE[4]


def _eseq_sequences(x):
    """Normalise an Eseq / SeqElist input into (list[Eseq], dictionary)."""
    if hasattr(x, "sequences"):  # SeqElist
        return list(x.sequences), list(x.dictionary)
    if hasattr(x, "codes"):  # single Eseq
        return [x], list(x.dictionary)
    seqs = list(x)
    return seqs, list(seqs[0].dictionary)


def _trajectory_key(sequence) -> tuple[tuple[int, ...], ...]:
    """R ``seqpcplot``'s trajectory identity: event codes grouped by distinct time.

    ``order.align = "first"`` ranks the distinct timestamps of a sequence, so
    simultaneous events share a rank (``seqpcplot.R:350-364``).
    """
    times = np.asarray(sequence.times, dtype=float)
    codes = np.asarray(sequence.codes, dtype=int)
    return tuple(
        tuple(sorted(int(code) for code in codes[times == time])) for time in np.unique(times)
    )


def _r_trajectory_text(key: tuple[tuple[int, ...], ...]) -> str:
    """The ``createstring`` text R sorts trajectories by (``factor(traj.string)``)."""
    return "-".join(
        "(" + ",".join(str(code) for code in group) + f")^{rank}"
        for rank, group in enumerate(key, start=1)
    )


def _time_ranked_events(sequence) -> tuple[np.ndarray, np.ndarray]:
    """R's points: x = rank of the event's distinct time, events of one time by level."""
    times = np.asarray(sequence.times, dtype=float)
    codes = np.asarray(sequence.codes, dtype=float)
    order = np.lexsort((codes, times))
    ranks = np.searchsorted(np.unique(times), times[order]) + 1
    return ranks.astype(float), codes[order]


def _state_trajectories(x: STSList):
    """A state sequence object as R's seqpcplot reads it: its DSS states at positions 1..k.

    R's ``seqdss`` drops missing cells before it collapses runs (``with.missing =
    FALSE``), so no missing level remains, and the y levels are the state labels
    (seqpcplot.R:133-157).
    """
    x = x._canonicalized()
    markers = {code for code in (x.nr_code, x.void_code) if code is not None}
    levels = [code for code in range(len(x.alphabet)) if code not in markers]
    level_of = {code: index + 1 for index, code in enumerate(levels)}
    weights = x.weights if x.weights is not None else np.ones(x.nseq)
    trajectories = []
    for row, weight in zip(x.data, weights):
        kept = [int(code) for code in row if int(code) in level_of]
        dss = [code for index, code in enumerate(kept) if index == 0 or code != kept[index - 1]]
        trajectories.append(
            SimpleNamespace(
                codes=np.array([level_of[code] for code in dss], dtype=int),
                times=np.arange(1, len(dss) + 1, dtype=float),
                weight=float(weight),
            )
        )
    labels = x.labels if x.labels is not None else x.alphabet
    return trajectories, [str(labels[code]) for code in levels]


def _filtered_sequences(seqs, filt) -> np.ndarray:
    """Per-sequence keep mask: ``filt`` evaluated on distinct trajectories.

    Mirrors ``seqpcplot.R:515-517`` and ``633-639``: the filter sees one
    weighted share per distinct trajectory, in R's trajectory order, followed
    by the share of sequences without events (``1 - sum``). A sequence is kept
    when its trajectory's filter value is positive (2026-09-22 audit S6-02).
    """
    keys = [_trajectory_key(sequence) for sequence in seqs]
    weights = np.array([getattr(sequence, "weight", 1.0) for sequence in seqs], dtype=float)
    total = weights.sum() if weights.sum() else 1.0
    trajectories = sorted({key for key in keys if key}, key=_r_trajectory_text)
    shares = [
        sum(weight for key, weight in zip(keys, weights) if key == trajectory) / total
        for trajectory in trajectories
    ]
    frequencies = np.array([*shares, 1.0 - sum(shares)], dtype=float)
    values = np.asarray(filt.apply(frequencies), dtype=float)
    kept = {trajectory for trajectory, value in zip(trajectories, values) if value > 0}
    return np.array([key in kept for key in keys], dtype=bool)


@_close_new_figures_on_error
def seqpcplot(x, *, filt=None, cpal=None, figsize=(8, 4), ax=None):
    """Parallel-coordinate event plot -- TraMineR ``seqpcplot`` / ``plot.seqpcplot``.

    Each event sequence is a polyline through its events in order: x = the rank
    of the event's distinct time within the sequence, so simultaneous events
    share an x (R's ``order.align = "first"``, seqpcplot.R:209-211), y = the event
    category (1-based code into the event dictionary). A state sequence object
    is drawn through its distinct successive states, as R does (seqpcplot.R:131-140;
    TX-AUD-20261011-027): x = the state's position in the DSS, y = the state.
    As R's ``seqdss`` does, missing cells are dropped before runs are
    collapsed, so no missing state is drawn. y-tick labels are the event
    dictionary or the state labels. When ``filt`` is a
    ``seqpcfilter`` object, it is evaluated as in R on the (weighted) share of
    each *distinct* trajectory -- events grouped by timestamp, so simultaneous
    events form one step -- followed by the share of sequences without events,
    and every sequence of a trajectory that keeps a positive filter value is
    drawn. Tied shares are ordered by R's trajectory text in code-point order;
    R uses the session collation, which can only matter for ties between
    trajectories that differ by simultaneous events. The event sequences come
    from ``seqecreate`` (faithful Eseq / SeqElist). NB: R aggregates distinct
    trajectories with frequency-weighted line widths; this draws one polyline
    per sequence. A supplied ``cpal`` is recycled over those polylines after
    filtering.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
    """
    if isinstance(x, STSList):
        seqs, dictionary = _state_trajectories(x)
    elif isinstance(x, (Eseq, SeqElist)):
        seqs, dictionary = _eseq_sequences(x)
    else:
        raise TypeError(
            "x must be an event-sequence (Eseq, SeqElist) or state-sequence (STSList) object"
        )
    colors = _trajectory_colors(cpal)
    keep = np.ones(len(seqs), dtype=bool)
    if filt is not None:
        keep = _filtered_sequences(seqs, filt)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    else:
        fig = ax.figure
    plotted_index = 0
    for s, k in zip(seqs, keep):
        if not k:
            continue
        positions, codes = _time_ranked_events(s)
        color = None if colors is None else colors[plotted_index % len(colors)]
        ax.plot(
            positions,
            codes,
            marker="o",
            alpha=0.6,
            color=color,
        )
        plotted_index += 1
    ax.set_yticks(np.arange(1, len(dictionary) + 1))
    ax.set_yticklabels(dictionary)
    ax.set_xlabel("Event rank" if not isinstance(x, STSList) else "Position")
    ax.set_ylabel("Event" if not isinstance(x, STSList) else "State")
    ax.set_title("Event-sequence parallel-coordinate plot")
    return fig, ax


@_close_new_figures_on_error
def plot_eseq(x, **kwargs):
    """Plot a single event sequence -- R ``plot.eseq`` (delegates to ``seqpcplot``)."""
    if not isinstance(x, Eseq):
        raise TypeError("x must be an Eseq event-sequence object")
    return seqpcplot(x, **kwargs)


@_close_new_figures_on_error
def plot_seqelist(x, **kwargs):
    """Plot a collection of event sequences -- R ``plot.seqelist`` (delegates to ``seqpcplot``)."""
    if not isinstance(x, SeqElist):
        raise TypeError("x must be a SeqElist event-sequence object")
    return seqpcplot(x, **kwargs)


@_close_new_figures_on_error
def seqeordplot(*args, **kwargs):
    """Deprecated event-ordering plot -- superseded by :func:`seqpcplot`.

    ``seqeordplot`` was **suppressed upstream in TraMineRextras 0.4.3
    (2018-08-22)** and superseded by TraMineR's ``seqpcplot`` (``NEWS.Rd``:
    "Function ``seqeordplot()`` suppressed. It is superseded by the function
    ``seqpcplot()`` of the TraMineR package."). traminex provides no analytical
    port; this thin forwarder emits a :class:`DeprecationWarning` and delegates to
    :func:`seqpcplot` so existing calls keep working. Prefer :func:`seqpcplot`.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
        Whatever :func:`seqpcplot` returns.
    """
    warnings.warn(
        "seqeordplot is deprecated: it was superseded by seqpcplot in "
        "TraMineRextras 0.4.3 (2018-08-22). Use seqpcplot instead.",
        DeprecationWarning,
        stacklevel=_caller_stacklevel(),
    )
    return seqpcplot(*args, **kwargs)
