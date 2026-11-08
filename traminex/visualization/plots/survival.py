"""Spell-survival plot over the R-golden'd seqsurv series (BUILD-from-R-design).

seqsurvplot mirrors TraMineRextras' plot.stslist.surv: per-state Kaplan-Meier
survival step curves (SOURCES/TraMineRextras/R/plot.stslist.surv.R) over
traminex's already-pinned seqsurv series. Original matplotlib (no Sequenzo
borrow); provenance docstring-only.
"""

from __future__ import annotations

import sys
from math import ceil
from numbers import Integral, Real
from typing import Literal, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.lines import Line2D

from ..._deprecated import MISSING, resolve_deprecated_alias
from ...sequences import STSList, _require_stslist
from ...stats_p2 import _seqsurv_factor, _seqsurv_group_palette, seqsurv
from ..palette import _state_color, ensure_cmap
from ..utils import _close_new_figures_on_error, _legend_with_labels, _reserve_legend_space
from .summary import _subset_stslist


def _draw_survival_curves(ax, result, colors, style):
    """Draw the established ``SeqSurv`` step series on one axis.

    As R's plot.survfit, each curve starts at (0, 1), steps at its event times
    and runs on to the stratum's last follow-up time, censored or not; a stratum
    without events is a flat line at 1 (TX-AUD-20261010-042, TX-AUD-20261011-026).
    """
    lines = []
    follow_up = getattr(result, "follow_up", {}) or {}
    for key in result.strata:
        arr = np.asarray(result.curves.get(key), dtype=float)
        last = follow_up.get(key)
        if arr.size == 0 and last is None:
            continue
        times = arr[:, 0] if arr.size else np.empty(0)
        survival = arr[:, 3] if arr.size else np.empty(0)
        xs = np.concatenate([[0.0], times])
        ys = np.concatenate([[1.0], survival])
        if last is not None and last > xs[-1]:
            xs = np.append(xs, last)
            ys = np.append(ys, ys[-1])
        options = dict(style)
        where = options.pop("where", "post")
        color = None if colors is None else colors.get(key)
        (line,) = ax.step(
            xs,
            ys,
            where=where,
            color=color,
            label=str(key),
            **options,
        )
        lines.append(line)
    return lines


def _layout(count: int, rows, cols) -> tuple[int, int]:
    def positive(value, name):
        if value is None:
            return None
        if isinstance(value, bool) or not isinstance(value, Integral) or value <= 0:
            raise ValueError(f"{name} must be a positive integer")
        return int(value)

    nrows = positive(rows, "rows")
    ncols = positive(cols, "cols")
    if nrows is None and ncols is None:
        ncols = min(count, 2)
        nrows = ceil(count / ncols)
    elif nrows is None:
        nrows = ceil(count / ncols)
    elif ncols is None:
        ncols = ceil(count / nrows)
    if nrows * ncols < count:
        raise ValueError("rows * cols must provide one panel per plot")
    return nrows, ncols


def _resolve_axes(count: int, rows: int, cols: int, figsize, supplied, use_layout: bool = True):
    if supplied is None:
        use_figsize = figsize if figsize is not None else (5 * cols, 3.5 * rows)
        fig, grid = plt.subplots(
            rows,
            cols,
            figsize=use_figsize,
            squeeze=False,
            layout="constrained" if use_layout else "none",
        )
        flat = list(np.asarray(grid, dtype=object).reshape(-1))
        panels = flat[:count]
        for unused in flat[count:]:
            unused.remove()
        return fig, panels

    panels = (
        [supplied]
        if isinstance(supplied, Axes)
        else list(np.asarray(supplied, dtype=object).reshape(-1))
    )
    if len(panels) != count:
        raise ValueError(f"supplied ax must contain exactly {count} axes")
    if any(not isinstance(axis, Axes) for axis in panels):
        raise TypeError("supplied ax must contain matplotlib Axes objects")
    fig = panels[0].figure
    if any(axis.figure is not fig for axis in panels):
        raise ValueError("all supplied axes must belong to the same figure")
    return fig, panels


def _observed_states(x: STSList, with_missing: bool):
    observed_codes = {int(code) for code in np.asarray(x.data).reshape(-1)}
    marker_codes = {code for code in (x.nr_code, x.void_code) if code is not None}
    codes = [
        code
        for code in range(len(x.alphabet))
        if code in observed_codes and code not in marker_codes
    ]
    if with_missing and x.nr_code is not None and x.nr_code in observed_codes:
        codes.append(int(x.nr_code))
    return codes, [str(x.alphabet[code]) for code in codes]


def _select_states(x: STSList, which_states, with_missing: bool):
    available_codes, available = _observed_states(x, with_missing)
    if which_states is None:
        return available_codes, available
    if isinstance(which_states, str):
        requested = [which_states]
    else:
        try:
            requested = [str(value) for value in which_states]
        except TypeError as exc:
            raise TypeError("which_states must be a sequence of state names") from exc
    requested_set = set(requested)
    selected_codes = [
        code for code, state in zip(available_codes, available) if state in requested_set
    ]
    selected = [str(x.alphabet[code]) for code in selected_codes]
    bad = [state for state in requested if state not in available]
    if not selected:
        raise ValueError("Only unvalid or unobserved selected states!")
    if bad:
        print(
            f"[!!] Unvalid or unobserved state(s): {' '.join(bad)}",
            file=sys.stderr,
        )
    return selected_codes, selected


def _state_labels(x: STSList, codes: list[int]) -> list[str]:
    labels = []
    for code in codes:
        if code == x.nr_code:
            labels.append("Missing")
        elif x.labels is not None:
            labels.append(str(x.labels[code]))
        else:
            labels.append(str(x.alphabet[code]))
    return labels


def _normalize_axis_mode(value, name: str, secondary: str):
    if isinstance(value, (bool, np.bool_)):
        return "all" if bool(value) else False
    if not isinstance(value, str) or value not in {"all", secondary}:
        raise ValueError(f"{name} must be a boolean, 'all', or {secondary!r}")
    return value


def _panel_text(value, labels: list[str], name: str, default: str | None = None):
    if value is None:
        return [default if default is not None else ""] * len(labels)
    if isinstance(value, str):
        if name == "main" and value == "auto":
            return list(labels)
        if name == "main":
            return [f"{value} - {label}" for label in labels]
        return [value] * len(labels)
    try:
        values = list(value)
    except TypeError as exc:
        raise TypeError(f"{name} must be a string or a sequence") from exc
    if len(values) != len(labels):
        raise ValueError(f"{name} must have one value per group (or state) level")
    return ["" if item is None else str(item) for item in values]


def _palette(values, expected: int, name: str):
    if values is None:
        return None
    if isinstance(values, str):
        palette = [values]
    else:
        try:
            palette = list(values)
        except TypeError as exc:
            raise TypeError(f"{name} must be a sequence of colors") from exc
    if len(palette) != expected:
        raise ValueError(f"{name} must provide one color per plotted curve")
    return palette


def _time_ticks(x: STSList, result, xtlab, xtstep: int, tick_last: bool):
    if xtlab is None:
        labels = None
    else:
        labels = list(xtlab)
        if len(labels) != x.seqlen:
            raise ValueError("xtlab must provide one label per sequence position")
    maximum = int(getattr(result, "max_time", 0)) if result is not None else 0
    if result is not None and maximum == 0:
        for curve in result.curves.values():
            values = np.asarray(curve, dtype=float)
            if values.size:
                maximum = max(maximum, int(np.nanmax(values[:, 0])))
    if maximum == 0:
        maximum = x.seqlen
    ticks = list(range(1, maximum + 1, xtstep))
    if maximum and tick_last and maximum not in ticks:
        ticks.append(maximum)
    tick_labels = [str(tick) if labels is None else str(labels[tick - 1]) for tick in ticks]
    return maximum, ticks, tick_labels


@_close_new_figures_on_error
@_require_stslist
def seqsurvplot(
    x,
    *,
    groups=None,
    with_missing: bool = False,
    cpal=None,
    figsize=(8, 4),
    ax=None,
):
    """Spell-survival curves -- re-expresses TraMineRextras ``seqsurvplot`` / ``plot.stslist.surv``.

    One Kaplan-Meier survival step curve per state: x = time since the start of a
    spell, y = the survival probability (share of un-terminated spells), coloured
    by state. The (time, surv) steps come from ``seqsurv`` (R-golden'd in
    ``seqsurv_*``). The step starts at ``(0, 1)``, drops at each event time and
    runs on to the state's last follow-up time; a state without events is a flat
    line at 1, as R's ``plot.stslist.surv`` draws one curve per state
    (TX-AUD-20261010-042, TX-AUD-20261019-051). R draws no legend; the legend
    here names the states by their labels, as ``seqsplot``'s does.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
    """
    x = x._canonicalized()
    res = seqsurv(x, groups=groups, with_missing=with_missing)
    cm = ensure_cmap(x, cpal)
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    else:
        fig = ax.figure
    alphabet = list(x.alphabet)
    colors = {
        state: _state_color(cm, alphabet.index(state) if state in alphabet else 0)
        for state in res.states
    }
    lines = _draw_survival_curves(ax, res, colors, {})
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Time since start of the spell")
    ax.set_ylabel("Percentage of unterminated spells")
    ax.set_title("Spell survival")
    # Explicit handles and the helper keep a state whose label starts with "_"
    # (TX-REP-20261017-003, TX-AUD-20261019-060). The states are named by their
    # labels, as seqsplot's legend names them (seqsplot.R:91-92; TX-AUD-20261011-067).
    _legend_with_labels(
        ax,
        lines,
        _state_labels(x, [alphabet.index(line.get_label()) for line in lines]),
        fontsize="small",
        ncol=min(len(res.states), 6),
    )
    return fig, ax


@_close_new_figures_on_error
def seqsplot(
    seqdata: STSList,
    group: Sequence[object] | None = None,
    *,
    main: str | Sequence[str] | None = "auto",
    cpal: Sequence[object] | None = None,
    missing_color: object = None,
    ylab: str | Sequence[str] | None = None,
    yaxis: Literal["all", "left"] | bool = "all",
    xaxis: Literal["all", "bottom"] | bool = "all",
    xtlab: Sequence[str] | None = None,
    xtstep: int | None = None,
    tick_last: bool | None = None,
    cex_axis: float = 1,
    with_legend: Literal["auto", "right"] | bool = "auto",
    ltext: Sequence[str] | None = None,
    cex_legend: float = 1,
    use_layout: bool | None = None,
    legend_prop: float | None = None,
    rows: int | None = None,
    cols: int | None = None,
    which_states: Sequence[str] | None = None,
    per_state: bool = False,
    state: object = None,
    with_missing: bool = False,
    figsize: tuple[float, float] | None = None,
    ax: object = None,
    title=MISSING,
    cex_plot=MISSING,
    withlegend=MISSING,
    axes=MISSING,
    **style: object,
):
    """Plot spell-survival curves by group or by state.

    Every survival estimate is produced by :func:`traminex.seqsurv`. Group mode
    creates one panel per used group level; ``per_state=True`` creates one panel
    per selected observed state and superposes group curves.
    """
    main = resolve_deprecated_alias(
        main,
        title,
        canonical_name="main",
        alias_name="title",
        function_name="seqsplot",
        default="auto",
    )
    cex_axis = resolve_deprecated_alias(
        cex_axis,
        cex_plot,
        canonical_name="cex_axis",
        alias_name="cex_plot",
        function_name="seqsplot",
        default=1,
    )
    with_legend = resolve_deprecated_alias(
        with_legend,
        withlegend,
        canonical_name="with_legend",
        alias_name="withlegend",
        function_name="seqsplot",
        default="auto",
    )
    xaxis = resolve_deprecated_alias(
        xaxis,
        axes,
        canonical_name="xaxis",
        alias_name="axes",
        function_name="seqsplot",
        default="all",
    )
    if not isinstance(seqdata, STSList):
        raise TypeError("seqdata must be an STSList state sequence object")
    # Both marker representations are equivalent inputs (identity for
    # in-alphabet objects; cycle-9 audit N-04) — seqsurvplot already
    # canonicalizes; seqsplot must too before _observed_states/_state_labels.
    seqdata = seqdata._canonicalized()
    if not isinstance(per_state, (bool, np.bool_)):
        raise TypeError("per_state must be a boolean")
    use_xtstep = seqdata.xtstep if xtstep is None else xtstep
    if isinstance(use_xtstep, bool) or not isinstance(use_xtstep, Integral) or use_xtstep <= 0:
        raise ValueError("xtstep must be a positive integer")
    if tick_last is not None and not isinstance(tick_last, (bool, np.bool_)):
        raise TypeError("tick_last must be a boolean or None")
    use_tick_last = seqdata.tick_last if tick_last is None else bool(tick_last)
    if state is not None:
        print("'state' argument will be ignored", file=sys.stderr)
    if use_layout is not None and not isinstance(use_layout, (bool, np.bool_)):
        raise TypeError("use_layout must be a boolean or None")
    if not isinstance(with_legend, (str, bool, np.bool_)) or with_legend not in {
        "auto",
        "right",
        True,
        False,
    }:
        raise ValueError("with_legend must be 'auto', 'right', True, or False")
    xaxis_mode = _normalize_axis_mode(xaxis, "xaxis", "bottom")
    yaxis_mode = _normalize_axis_mode(yaxis, "yaxis", "left")
    try:
        axis_scale = float(cex_axis)
        legend_scale = float(cex_legend)
    except (TypeError, ValueError) as exc:
        raise ValueError("cex_axis and cex_legend must be positive numbers") from exc
    if not all(np.isfinite(value) and value > 0 for value in (axis_scale, legend_scale)):
        raise ValueError("cex_axis and cex_legend must be positive numbers")
    if legend_prop is not None and (
        not isinstance(legend_prop, Real) or not 0 < float(legend_prop) < 1
    ):
        raise ValueError("legend_prop must be strictly between zero and one")
    layout_enabled = (
        group is not None or with_legend is not False if use_layout is None else bool(use_layout)
    )
    if "color" in style or "label" in style:
        raise ValueError("color and label must be controlled by cpal and ltext")
    legend_ncol = style.pop("ncol", None)

    selected_codes, selected_states = _select_states(seqdata, which_states, bool(with_missing))
    selected_labels = _state_labels(seqdata, selected_codes)
    if group is None:
        group_values = np.full(seqdata.nseq, "Unique group", dtype=object)
        group_levels = ["Unique group"]
    else:
        try:
            group_values, group_levels = _seqsurv_factor(group, seqdata.nseq)
        except ValueError as exc:
            if "length" in str(exc):
                raise ValueError(
                    "group must contain one value for each row in the sequence object"
                ) from exc
            raise
    if not group_levels:
        raise ValueError("group must contain at least one non-missing value")
    grouped = pd.Categorical(group_values, categories=group_levels, ordered=True)

    if per_state:
        panel_keys = list(selected_states)
        panel_labels = list(selected_labels)
        curve_keys = list(group_levels)
        curve_labels = list(group_levels)
        default_palette = _seqsurv_group_palette(len(group_levels))
        if default_palette is None and cpal is None:
            print(
                " [!] too many groups (> 8), no automatic color palette assignation",
                file=sys.stderr,
            )
        curve_palette = _palette(cpal, len(group_levels), "cpal")
        if curve_palette is None:
            curve_palette = default_palette
        if curve_palette is None:
            # Explicit cycle colours, so the legend shows the curves' colours
            # (2026-09-22 audit S6-10).
            curve_palette = [f"C{index % 10}" for index in range(len(group_levels))]
        results = [
            seqsurv(
                seqdata,
                groups=grouped,
                per_state=True,
                state=panel_state,
                with_missing=bool(with_missing),
            )
            for panel_state in panel_keys
        ]
    else:
        panel_keys = list(group_levels)
        panel_labels = list(group_levels)
        curve_keys = list(selected_states)
        curve_labels = list(selected_labels)
        default_palette = []
        state_cmap = None
        # A palette shorter than the alphabet falls back as a whole, as in
        # every other renderer (2026-09-22 audit S6-14).
        full_palette = seqdata.cpal is not None and len(seqdata.cpal) >= len(seqdata.alphabet)
        for code in selected_codes:
            if code == seqdata.nr_code and missing_color is not None:
                default_palette.append(missing_color)
            elif full_palette:
                default_palette.append(seqdata.cpal[code])
            elif code == seqdata.nr_code:
                default_palette.append(seqdata.missing_color)
            else:
                # Without a palette entry, the shared default state palette:
                # the legend then shows the curve's colour and a state keeps
                # its colour in every panel (2026-09-22 audit S6-10).
                if state_cmap is None:
                    state_cmap = ensure_cmap(seqdata, None)
                default_palette.append(_state_color(state_cmap, code))
        curve_palette = _palette(cpal, len(selected_states), "cpal")
        if curve_palette is None:
            curve_palette = default_palette
        results = []
        for level in group_levels:
            positions = np.flatnonzero(group_values == level)
            subset = _subset_stslist(seqdata, positions)
            _, local_available = _observed_states(subset, bool(with_missing))
            local_states = [state for state in selected_states if state in local_available]
            results.append(
                None
                if not local_states
                else seqsurv(
                    subset,
                    state=local_states,
                    with_missing=bool(with_missing),
                )
            )

    if ltext is not None:
        curve_labels = [ltext] if isinstance(ltext, str) else [str(value) for value in ltext]
        if len(curve_labels) != len(curve_keys):
            raise ValueError("ltext must provide one value per plotted curve")
    titles = _panel_text(main, panel_labels, "main")
    ylabels = _panel_text(
        ylab,
        panel_labels,
        "ylab",
        default="Percentage of unterminated spells",
    )
    nrows, ncols = _layout(len(panel_keys), rows, cols)
    fig, panels = _resolve_axes(
        len(panel_keys),
        nrows,
        ncols,
        figsize,
        ax,
        use_layout=layout_enabled,
    )
    color_by_key = {key: color for key, color in zip(curve_keys, curve_palette)}
    bottom = {
        max(index for index in range(len(panels)) if index % ncols == column)
        for column in range(ncols)
        if any(index % ncols == column for index in range(len(panels)))
    }
    left = {index for index in range(len(panels)) if index % ncols == 0}
    font_size = float(plt.rcParams["font.size"]) * axis_scale

    for index, (axis, result, title, ylabel) in enumerate(zip(panels, results, titles, ylabels)):
        if result is not None:
            _draw_survival_curves(axis, result, color_by_key, style)
        maximum, ticks, tick_labels = _time_ticks(
            seqdata,
            result,
            xtlab,
            int(use_xtstep),
            bool(use_tick_last),
        )
        axis.set_ylim(0, 1.02)
        axis.set_xlim(0, max(1, maximum))
        axis.set_title(title)
        axis.set_xlabel("Time since start of the spell")
        axis.set_ylabel(ylabel)

        show_x = xaxis_mode == "all" or (xaxis_mode == "bottom" and index in bottom)
        if show_x:
            axis.set_xticks(ticks, tick_labels)
        else:
            axis.set_xticks([])
        show_y = yaxis_mode == "all" or (yaxis_mode == "left" and index in left)
        if not show_y:
            axis.set_yticks([])
        axis.tick_params(axis="both", labelsize=font_size)

    if with_legend is True or with_legend in {"auto", "right"}:
        handles = [
            Line2D(
                [],
                [],
                color=curve_palette[index],
                label=label,
                **{key: value for key, value in style.items() if key in {"alpha", "linewidth"}},
            )
            for index, label in enumerate(curve_labels)
        ]
        legend_options = {}
        if legend_prop is not None:
            position = "right" if with_legend == "right" else "bottom"
            legend_options["bbox_to_anchor"] = _reserve_legend_space(
                fig,
                panels,
                float(legend_prop),
                position,
            )
            legend_options["bbox_transform"] = fig.transFigure
            location = "center"
        else:
            location = "center right" if with_legend == "right" else "lower center"
        ncol = min(len(handles), 6) if legend_ncol is None else max(1, int(legend_ncol))
        _legend_with_labels(
            fig,
            handles,
            curve_labels,
            loc=location,
            fontsize=font_size * legend_scale,
            ncol=ncol,
            **legend_options,
        )

    returned_axes = panels[0] if len(panels) == 1 else np.asarray(panels, dtype=object)
    return fig, returned_axes


__all__ = ["seqsplot", "seqsurvplot"]
