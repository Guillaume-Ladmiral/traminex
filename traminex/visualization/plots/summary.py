"""State-summary plots over the R-golden'd series (BUILD-from-R-design).

seqmtplot / seqmsplot / seqHtplot / seqdHplot mirror TraMineR's
plot.stslist.meant / .modst / .statd drawing routines. They are built directly
from the R source (SOURCES/TraMineR/R/plot.stslist.{meant,modst,statd}.R) over
traminex's already-pinned series (seqmeant / seqmodst / seqstatd); the matplotlib
code is original (no Sequenzo borrow), so provenance is docstring-only.
"""

from __future__ import annotations

import warnings
from collections.abc import Mapping

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import t as _student_t

from ...sequences import _require_stslist
from ...stats import seqmeant, seqstatd
from ...stats_p1 import seqmodst
from ..palette import _state_color, ensure_cmap
from ..series import _factor_panels, grouped_transversal_entropy, transversal_entropy
from ..ticks import _apply_time_axis_on_cells, apply_time_axis
from ..utils import _close_new_figures_on_error, _legend_with_labels
from .distribution import state_distribution_plot


@_close_new_figures_on_error
@_require_stslist
def seqmtplot(
    x,
    *,
    group=None,
    main="auto",
    weighted: bool = True,
    include_missing: bool = False,
    cpal=None,
    errbar: bool = False,
    figsize=(6, 4),
    ax=None,
):
    """Mean time in each state -- re-expresses TraMineR ``seqmtplot`` / ``plot.stslist.meant``.

    One vertical bar per row of R's ``seqmeant``: the states of the alphabet,
    plus the missing state with ``include_missing`` (R's ``with.missing``, a bar
    of height 0 when nothing is missing); height = mean time spent in the state,
    coloured by the state palette (the missing state by ``missing_color``).
    With ``errbar=True`` (R's ``serr = TRUE``; R's default draws none,
    seqmeant.R:5) 95% confidence whiskers ``Mean +/- qt(.975, n-1) * SE`` are
    overlaid when ``df = nbseq - 1 >= 1`` (SOURCES/TraMineR/R/plot.stslist.meant.R;
    TX-AUD-20261010-043, -058). y-axis spans ``0..L``.

    The Mean / SE come straight from ``seqmeant`` (R-golden'd in ``seqmeant_*``).

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
    """
    if group is not None or main != "auto":
        from .dispatch import seqplot

        # Called directly the renderer draws no legend; main and group keep it so.
        return seqplot(
            x,
            type="mt",
            with_legend=False,
            group=group,
            main=main,
            weighted=weighted,
            include_missing=include_missing,
            cpal=cpal,
            errbar=errbar,
            figsize=figsize,
            ax=ax,
        )
    # Both marker representations are equivalent inputs (identity for
    # in-alphabet objects; cycle-9 audit N-04/N-05).
    x = x._canonicalized()
    df = seqmeant(x, weighted=weighted, include_missing=include_missing, prop=False, serr=True)
    # R's rows: the real states, then the missing state with with.missing
    # (seqmeant pads its result with both markers; TX-AUD-20261010-058).
    markers = {code for code in (x.nr_code, x.void_code) if code is not None}
    rows = [str(state) for code, state in enumerate(x.alphabet) if code not in markers]
    if include_missing:
        rows.append(str(x.nr))
    df = df.rename(index=str).reindex(rows, fill_value=0.0)
    means = df["Mean"].to_numpy(dtype=float)
    se = df["SE"].to_numpy(dtype=float)
    states = list(df.index)
    K = len(states)
    L = x.data.shape[1]
    # R plot.stslist.meant: df = nbseq - 1 where nbseq is seqmeant's weighted
    # case count (sum of weights), not the raw number of rows.
    if weighted and x.weights is not None:
        n = float(np.sum(x.weights))
    else:
        n = float(x.nseq)

    cm = ensure_cmap(x, cpal)
    colors = [
        (
            x.missing_color
            if include_missing and i == K - 1
            else _state_color(cm, x.alphabet.index(s) if s in x.alphabet else i)
        )
        for i, s in enumerate(states)
    ]

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    else:
        fig = ax.figure

    xs = np.arange(K)
    yerr = None
    if errbar and (n - 1) >= 1:
        q = float(_student_t.ppf(0.975, n - 1))
        yerr = q * se
    ax.bar(xs, means, color=colors, yerr=yerr, capsize=4 if yerr is not None else 0)
    ax.set_xticks(xs)
    ax.set_xticklabels([str(s) for s in states])
    ax.set_ylim(0, L)
    ax.set_xlabel("State")
    ax.set_ylabel("Mean time")
    ax.set_title("Mean time in each state")
    return fig, ax


@_close_new_figures_on_error
@_require_stslist
def seqmsplot(
    x,
    *,
    group=None,
    main="auto",
    weighted: bool = True,
    include_missing: bool = False,
    cpal=None,
    info: bool = True,
    figsize=(8, 3),
    ax=None,
):
    """Modal state sequence -- re-expresses TraMineR ``seqmsplot`` / ``plot.stslist.modst``.

    One bar per time position: height = the frequency of the modal (most frequent)
    state at that position, coloured by that modal state. Mirrors R's
    ``barplot(prof.freq, space=0, col=cpal, ylim=c(0, 1.2))`` where
    ``prof.freq[modal_state, t] = Frequencies[t]``
    (SOURCES/TraMineR/R/plot.stslist.modst.R). When ``info`` is set, the modal
    sequence's occurrence count and its frequency, relative to the weighted
    number of sequences as in R, are annotated.

    The modal codes / frequencies / occurrences come from ``seqmodst`` (R-golden'd
    in ``seqmodst_*``).

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
    """
    if group is not None or main != "auto":
        from .dispatch import seqplot

        # Called directly the renderer draws no legend; main and group keep it so.
        return seqplot(
            x,
            type="ms",
            with_legend=False,
            group=group,
            main=main,
            weighted=weighted,
            include_missing=include_missing,
            cpal=cpal,
            info=info,
            figsize=figsize,
            ax=ax,
        )
    # Both marker representations are equivalent inputs (identity for
    # in-alphabet objects; cycle-9 audit N-04/N-05).
    x = x._canonicalized()
    m = seqmodst(x, weighted=weighted, include_missing=include_missing)
    modal_codes = np.asarray(m.data[0], dtype=int)
    freqs = np.asarray(getattr(m, "frequencies_"), dtype=float)
    L = len(freqs)
    cm = ensure_cmap(x, cpal)
    colors = [_state_color(cm, int(c)) for c in modal_codes]

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    else:
        fig = ax.figure

    xs = np.arange(L)
    ax.bar(xs, freqs, width=1.0, align="edge", color=colors)
    ax.set_ylim(0, 1.2)
    ax.set_xlim(0, L)
    _apply_time_axis_on_cells(ax, x, L)
    ax.set_xlabel("Time")
    ax.set_ylabel("State freq.")
    ax.set_title("Modal state sequence")
    if info:
        # R plot.stslist.modst: n <- attr(x, "nbseq"), the sum of the weights
        # (seqmodst.R), and paste() prints round(., 1) without a trailing ".0"
        # (2026-09-22 audit S6-05).
        n = float(getattr(m, "nbseq_"))
        nbrep = int(getattr(m, "occurrences_"))
        ctfreq = round((nbrep / n) * 100, 1) if n else 0.0
        ax.text(
            L / 2.0,
            1.1,
            f"Modal state sequence ({nbrep} occurrences, freq={ctfreq:.15g}%)",
            ha="center",
            fontsize="small",
        )
    return fig, ax


def _subset_stslist(x, rows):
    """A new STSList holding only ``rows`` of ``x`` (for grouped series)."""
    from ...sequences import STSList, _remap_issues

    row_indices = np.atleast_1d(np.arange(x.nseq)[rows])
    return STSList(
        data=x.data[rows],
        alphabet=list(x.alphabet),
        labels=None if x.labels is None else list(x.labels),
        cpal=None if x.cpal is None else list(x.cpal),
        xtlab=None if x.xtlab is None else list(x.xtlab),
        xtstep=x.xtstep,
        tick_last=x.tick_last,
        nr=x.nr,
        void=x.void,
        nr_code=x.nr_code,
        void_code=x.void_code,
        weights=None if x.weights is None else x.weights[rows],
        start=x.start,
        ids=None if x.ids is None else [x.ids[int(i)] for i in row_indices],
        missing_color=x.missing_color,
        issues=_remap_issues(x.issues, row_indices),
    )


@_close_new_figures_on_error
@_require_stslist
def seqHtplot(x, *, group=None, include_missing: bool = False, figsize=(8, 3), ax=None):
    """Transversal entropy over time -- re-expresses TraMineR ``seqHtplot`` (type="Ht").

    A line of the per-position cross-sectional (Shannon) entropy of the state
    distribution, ``seqstatd$Entropy`` (SOURCES/TraMineR/R/plot.stslist.statd.R,
    ``type="Ht"``). The entropy series is R-golden'd in ``seqstatd_*_ve``.

    When ``group`` (one value per sequence) is given, one entropy line is drawn per
    group level -- the grouped transversal-entropy plot of TraMineRextras
    ``seqplot.tentrop`` (folded in here; each line is the ``seqstatd`` entropy of
    that group's subset). R's ``seqHtplot(group=)`` draws one panel per level
    instead of overlaying lines. Levels follow R ``factor`` order (text in R's
    collation), and sequences with a missing group value are left out, as in
    R. R's ``seqplot`` builds its groups with the internal ``group()``, which
    orders numeric groups as text (``1, 10, 2``) and titles every level with a
    leading space; traminex keeps ``factor``'s numeric order and plain titles,
    a recorded correction (TX-REP-20261017-002).

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
    """
    L = x.data.shape[1]
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    else:
        fig = ax.figure
    xs = np.arange(1, L + 1)
    ymax = 1.0
    if group is None:
        ent = np.asarray(seqstatd(x, include_missing=include_missing)["Entropy"], dtype=float)
        ax.plot(xs, ent, color="black", linewidth=1.5)
        ymax = max(ymax, float(ent.max()) if L else 1.0)
    else:
        if len(group) != x.nseq:
            raise ValueError("group must contain one value for each sequence")
        # Missing group values belong to no group, as with R's group()
        # (2026-09-22 audit S6-07).
        handles, levels = [], []
        for rows, level in zip(*_factor_panels(group)):
            sub = _subset_stslist(x, np.asarray(rows, dtype=int))
            ent = np.asarray(seqstatd(sub, include_missing=include_missing)["Entropy"], dtype=float)
            (line,) = ax.plot(xs, ent, linewidth=1.5, label=level)
            handles.append(line)
            levels.append(level)
            ymax = max(ymax, float(ent.max()) if L else 1.0)
        # Explicit labels through the helper: matplotlib drops a label that
        # starts with "_", such as a group named "_c".
        _legend_with_labels(ax, handles, levels, fontsize="small", title="Group")
    ax.set_ylim(0, ymax)
    ax.set_xlim(0.5, L + 0.5)
    apply_time_axis(ax, x, L, one_based=True)
    ax.set_xlabel("Time")
    ax.set_ylabel("Entropy")
    ax.set_title("Transversal entropy (Ht)")
    return fig, ax


_DARK3 = {
    1: ("#E16A86",),
    2: ("#E16A86", "#00AD9A"),
    3: ("#E16A86", "#50A315", "#009ADE"),
    4: ("#E16A86", "#909800", "#00AD9A", "#9183E6"),
    5: ("#E16A86", "#AA9000", "#00AA5A", "#00A6CA", "#B675E0"),
    6: ("#E16A86", "#B88A00", "#50A315", "#00AD9A", "#009ADE", "#C86DD7"),
    7: ("#E16A86", "#C18500", "#799D00", "#00AB6E", "#00A9BE", "#6C8EE6", "#D169D0"),
    8: (
        "#E16A86",
        "#C7821C",
        "#909800",
        "#00A846",
        "#00AD9A",
        "#00A2D3",
        "#9183E6",
        "#D766C9",
    ),
}
_SET1 = (
    "#E41A1C",
    "#377EB8",
    "#4DAF4A",
    "#984EA3",
    "#FF7F00",
    "#FFFF33",
    "#A65628",
    "#F781BF",
    "#999999",
)
_LINE_STYLES = {
    "solid": "-",
    "dashed": "--",
    "dotted": ":",
    "dotdash": "-.",
}


def _cycle_plot_values(value, default, count: int, name: str):
    if value is None:
        values = list(default)
    elif isinstance(value, str) or np.isscalar(value):
        values = [value]
    else:
        values = list(value)
    if not values:
        raise ValueError(f"{name} must not be empty")
    return [values[index % len(values)] for index in range(count)]


def _dark3_colors(count: int) -> list[str]:
    if count in _DARK3:
        return list(_DARK3[count])
    base = _DARK3[max(_DARK3)]
    return [base[index % len(base)] for index in range(count)]


def _legend_location(value: str) -> str:
    compact = value.replace(" ", "").lower()
    aliases = {
        "topright": "upper right",
        "topleft": "upper left",
        "bottomright": "lower right",
        "bottomleft": "lower left",
    }
    return aliases.get(compact, value)


def _draw_tentrop(
    curves,
    labels,
    *,
    npos,
    axis_source,
    list_route,
    main,
    col,
    lty,
    lwd,
    ylim,
    xtlab,
    xtstep,
    tick_last,
    with_legend,
    glabels,
    legend_pos,
    horiz,
    cex_legend,
    figsize,
    ax,
    style,
):
    count = len(curves)
    default_colors = _SET1 if list_route else _dark3_colors(count)
    colors = _cycle_plot_values(col, default_colors, count, "col")
    default_lty = ["solid", "dashed", "dotted", "solid", "dashed"]
    line_types = _cycle_plot_values(lty, default_lty, count, "lty")
    line_widths = _cycle_plot_values(lwd, [3.5], count, "lwd")
    if glabels is None:
        legend_labels = [str(label) for label in labels]
    else:
        legend_labels = [str(label) for label in glabels]
        if len(legend_labels) != count:
            raise ValueError("glabels must contain one label for each entropy curve")

    forbidden = {"color", "linestyle", "linewidth", "label"}.intersection(style)
    if forbidden:
        raise ValueError(f"use col/lty/lwd/glabels instead of style keys {sorted(forbidden)}")
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    else:
        fig = ax.figure

    handles = []
    for index, curve in enumerate(curves):
        linestyle = _LINE_STYLES.get(line_types[index], line_types[index])
        (line,) = ax.plot(
            np.arange(1, len(curve) + 1),
            curve,
            color=colors[index],
            linestyle=linestyle,
            linewidth=float(line_widths[index]),
            label=legend_labels[index],
            **style,
        )
        handles.append(line)

    finite = np.concatenate([curve[np.isfinite(curve)] for curve in curves])
    if ylim is None:
        if finite.size:
            lower = np.floor(10 * finite.min()) / 10
            upper = np.ceil(10 * finite.max()) / 10
            if lower == upper:
                upper = lower + 0.1
            ylim = (float(lower), float(upper))
        else:
            ylim = (0.0, 1.0)
    if len(ylim) != 2 or not np.all(np.isfinite(ylim)) or ylim[0] >= ylim[1]:
        raise ValueError("ylim must contain two finite increasing values")
    ax.set_ylim(*ylim)
    ax.set_xlim((1, npos) if npos > 1 else (0.5, 1.5))

    if xtlab is None:
        if list_route:
            labels_at_time = [f"t{position}" for position in range(1, npos + 1)]
        elif axis_source.xtlab is not None:
            labels_at_time = [str(value) for value in axis_source.xtlab]
        else:
            labels_at_time = [str(axis_source.start + index) for index in range(npos)]
    else:
        labels_at_time = [str(value) for value in xtlab]
    if len(labels_at_time) < npos:
        raise ValueError("xtlab must contain at least one label for each time position")
    step = axis_source.xtstep if xtstep is None else xtstep
    if isinstance(step, bool) or not isinstance(step, (int, np.integer)) or step <= 0:
        raise ValueError("xtstep must be a positive integer")
    use_tick_last = axis_source.tick_last if tick_last is None else bool(tick_last)
    tick_positions = list(range(1, npos + 1, int(step)))
    if npos and use_tick_last and npos not in tick_positions:
        tick_positions.append(npos)
    ax.set_xticks(tick_positions)
    ax.set_xticklabels([labels_at_time[position - 1] for position in tick_positions])
    ax.set_xlabel("")
    ax.set_ylabel("Entropy")
    ax.set_title("Tranversal Entropies" if main == "auto" else ("" if main is None else main))
    if with_legend:
        fontsize = float(plt.rcParams["font.size"]) * float(cex_legend)
        # Explicit labels through the helper: matplotlib drops a label that
        # starts with "_", such as a group named "_c".
        _legend_with_labels(
            ax,
            handles,
            legend_labels,
            loc=_legend_location(legend_pos),
            ncols=count if horiz else 1,
            fontsize=fontsize,
        )
    return fig, ax


@_close_new_figures_on_error
@_require_stslist
def seqplot_tentrop(
    seqdata,
    group,
    *,
    main: str | None = "auto",
    col=None,
    lty=None,
    lwd=3.5,
    ylim=None,
    xtlab=None,
    xtstep: int | None = None,
    tick_last: bool | None = None,
    with_legend: bool = True,
    glabels=None,
    legend_pos: str = "upper right",
    horiz: bool = False,
    cex_legend: float = 1,
    figsize=(8, 3),
    ax=None,
    **style,
):
    """Plot superposed ``seqstatd`` entropy curves by factor group."""
    labels, curves = grouped_transversal_entropy(seqdata, group)
    return _draw_tentrop(
        curves,
        labels,
        npos=seqdata.seqlen,
        axis_source=seqdata,
        list_route=False,
        main=main,
        col=col,
        lty=lty,
        lwd=lwd,
        ylim=ylim,
        xtlab=xtlab,
        xtstep=xtstep,
        tick_last=tick_last,
        with_legend=with_legend,
        glabels=glabels,
        legend_pos=legend_pos,
        horiz=horiz,
        cex_legend=cex_legend,
        figsize=figsize,
        ax=ax,
        style=style,
    )


@_close_new_figures_on_error
def seqplot_tentrop_m(
    seqdata_list,
    *,
    main: str | None = "auto",
    col=None,
    lty=None,
    lwd=3.5,
    ylim=None,
    xtlab=None,
    xtstep: int | None = None,
    tick_last: bool | None = None,
    with_legend: bool = True,
    glabels=None,
    legend_pos: str = "upper right",
    horiz: bool = False,
    cex_legend: float = 1,
    figsize=(8, 3),
    ax=None,
    **style,
):
    """Plot entropy curves for a list, tuple, or mapping of sequence objects."""
    from ...sequences import STSList

    if isinstance(seqdata_list, STSList):
        raise TypeError("seqdata_list must be a list of at least 2 STSList objects")
    if isinstance(seqdata_list, Mapping):
        values = list(seqdata_list.values())
        default_labels = list(seqdata_list.keys())
    elif isinstance(seqdata_list, (list, tuple)):
        values = list(seqdata_list)
        default_labels = [f"seq{index}" for index in range(1, len(values) + 1)]
    else:
        raise TypeError("seqdata_list must be a list of at least 2 STSList objects")
    if len(values) < 2:
        raise ValueError("seqdata_list must contain at least 2 STSList objects")
    if any(not isinstance(value, STSList) for value in values):
        raise TypeError("every seqdata_list member must be an STSList")
    widths = [value.seqlen for value in values]
    if len(set(widths)) > 1:
        warnings.warn(
            "sequence objects in seqdata_list are not all of same length",
            UserWarning,
            stacklevel=2,
        )
    curves = tuple(transversal_entropy(value) for value in values)
    return _draw_tentrop(
        curves,
        default_labels,
        npos=max(widths),
        axis_source=values[0],
        list_route=True,
        main=main,
        col=col,
        lty=lty,
        lwd=lwd,
        ylim=ylim,
        xtlab=xtlab,
        xtstep=xtstep,
        tick_last=tick_last,
        with_legend=with_legend,
        glabels=glabels,
        legend_pos=legend_pos,
        horiz=horiz,
        cex_legend=cex_legend,
        figsize=figsize,
        ax=ax,
        style=style,
    )


@_close_new_figures_on_error
@_require_stslist
def seqdHplot(
    x,
    *,
    group=None,
    main="auto",
    include_missing: bool = False,
    weights=None,
    figsize=(8, 3),
    ax=None,
):
    """State distribution with an entropy overlay -- TraMineR ``seqdHplot`` (type="dH").

    Stacked per-position state-distribution bars (``seqstatd$Frequencies``) with the
    transversal-entropy line (``seqstatd$Entropy``) overlaid -- R's
    ``plot.stslist.statd(type="dH")``. Delegates to ``state_distribution_plot`` in
    its ``mode="dH"`` (already structure-VERIFY'd to draw the seqstatd series).

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
    """
    if group is not None or main != "auto":
        from .dispatch import seqplot

        return seqplot(
            x,
            type="dH",
            group=group,
            main=main,
            include_missing=include_missing,
            weights=weights,
            figsize=figsize,
            ax=ax,
        )
    fig, ax = state_distribution_plot(
        x, include_missing=include_missing, weights=weights, figsize=figsize, ax=ax, mode="dH"
    )
    return fig, ax
