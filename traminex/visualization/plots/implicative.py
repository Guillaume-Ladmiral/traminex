"""Reference renderer for TraMineRextras ``plot.seqimplic``."""

from __future__ import annotations

from math import ceil, sqrt
from numbers import Integral, Real
from typing import Literal, Sequence

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.lines import Line2D
from scipy.stats import norm

from ...stats_p2 import SeqImplic
from ..palette import _state_color, ensure_cmap
from ..utils import _close_new_figures_on_error, _legend_with_labels, _reserve_legend_space


def _validate_result(x: SeqImplic) -> np.ndarray:
    if not isinstance(x, SeqImplic):
        raise TypeError("x must be a SeqImplic result")
    try:
        values = np.asarray(x.indices, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("SeqImplic indices must be numeric") from exc
    if values.ndim != 3:
        raise ValueError("SeqImplic indices must be a three-dimensional array")
    expected = (len(x.levels), len(x.alphabet), len(x.positions))
    if values.shape != expected:
        raise ValueError("SeqImplic indices dimensions must match levels, alphabet, and positions")
    if x.type not in {"indice", "intensity"}:
        raise ValueError("SeqImplic type must be 'indice' or 'intensity'")
    if x.labels is not None and len(x.labels) != len(x.alphabet):
        raise ValueError("SeqImplic labels must match alphabet length")
    if x.cpal is not None and len(x.cpal) != len(x.alphabet):
        raise ValueError("SeqImplic cpal must match alphabet length")
    return values


def _select_levels(x: SeqImplic, only_levels) -> list[int]:
    if only_levels is None:
        return list(range(len(x.levels)))
    if isinstance(only_levels, (str, bytes)):
        requested = [only_levels]
    else:
        try:
            requested = list(only_levels)
        except TypeError as exc:
            raise TypeError("only_levels must be a sequence of level names") from exc
    selected = [index for index, level in enumerate(x.levels) if level in requested]
    if not selected:
        raise ValueError("only_levels must select at least one declared level")
    return selected


def _positive_integer(value, name: str) -> int | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, Integral) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def _layout(count: int, rows, cols) -> tuple[int, int]:
    nrows = _positive_integer(rows, "rows")
    ncols = _positive_integer(cols, "cols")
    if nrows is None and ncols is None:
        ncols = max(1, ceil(sqrt(count)))
        nrows = ceil(count / ncols)
    elif nrows is None:
        assert ncols is not None
        nrows = ceil(count / ncols)
    elif ncols is None:
        ncols = ceil(count / nrows)
    if nrows * ncols < count:
        raise ValueError("rows * cols must provide one panel per selected level")
    return nrows, ncols


def _resolve_axes(count: int, rows: int, cols: int, figsize, supplied):
    if supplied is None:
        fig, grid = plt.subplots(
            rows,
            cols,
            figsize=figsize,
            squeeze=False,
            constrained_layout=True,
        )
        flat = list(np.asarray(grid, dtype=object).reshape(-1))
        panels = flat[:count]
        for unused in flat[count:]:
            unused.remove()
        return fig, panels

    if isinstance(supplied, Axes):
        panels = [supplied]
    else:
        try:
            panels = list(np.asarray(supplied, dtype=object).reshape(-1))
        except (TypeError, ValueError) as exc:
            raise TypeError("ax must be a matplotlib Axes or a sequence of Axes") from exc
    if len(panels) != count:
        raise ValueError(f"supplied ax must contain exactly {count} axes")
    if any(not isinstance(axis, Axes) for axis in panels):
        raise TypeError("supplied ax must contain matplotlib Axes objects")
    fig = panels[0].figure
    if any(axis.figure is not fig for axis in panels):
        raise ValueError("all supplied axes must belong to the same figure")
    return fig, panels


def _confidence_levels(conf_level, result_type: str) -> list[float]:
    if conf_level is None:
        return []
    if isinstance(conf_level, (str, bytes)):
        raise TypeError("conf_level must contain probabilities")
    try:
        raw = np.asarray([conf_level] if np.isscalar(conf_level) else conf_level, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("conf_level must contain numeric probabilities") from exc
    values = raw.reshape(-1)
    if np.any(~np.isfinite(values)) or np.any((values <= 0) | (values >= 1)):
        raise ValueError("conf_level probabilities must be strictly between zero and one")
    if result_type == "indice":
        return [float(norm.ppf(value)) for value in values]
    return [float(value) for value in values]


def _resolve_ylim(values: np.ndarray, result_type: str, ylim) -> tuple[float, float]:
    if ylim is not None:
        try:
            limits = np.asarray(ylim, dtype=float).reshape(-1)
        except (TypeError, ValueError) as exc:
            raise ValueError("ylim must contain two finite increasing values") from exc
        if limits.size != 2 or not np.all(np.isfinite(limits)) or limits[0] >= limits[1]:
            raise ValueError("ylim must contain two finite increasing values")
        return float(limits[0]), float(limits[1])
    if result_type == "intensity":
        return 0.0, 1.0
    finite = values[np.isfinite(values)]
    maximum = float(finite.max()) if finite.size else 0.0
    return (0.0, maximum) if maximum > 0 else (0.0, 1.0)


def _resolve_titles(main, levels: list[str], all_levels: list[str]) -> list[str]:
    if main is None:
        return list(levels)
    if isinstance(main, str):
        return [f"{main} - {level}" for level in levels]
    try:
        titles = [str(value) for value in main]
    except TypeError as exc:
        raise TypeError("main must be a string or a sequence of titles") from exc
    if len(titles) == len(levels):
        return titles
    if len(titles) == len(all_levels):
        by_level = dict(zip(all_levels, titles))
        return [by_level[level] for level in levels]
    raise ValueError("main must provide one title per plotted or declared level")


def _tick_spec(labels: Sequence, step: int, tick_last: bool):
    indexes = list(range(0, len(labels), step))
    if tick_last and labels and len(labels) - 1 not in indexes:
        indexes.append(len(labels) - 1)
    return [index + 0.5 for index in indexes], [str(labels[index]) for index in indexes]


@_close_new_figures_on_error
def plot_seqimplic(
    x: SeqImplic,
    *,
    main: str | Sequence[str] | None = None,
    ylim: tuple[float, float] | None = None,
    xaxis: bool = True,
    ylab: str | None = "Implication",
    yaxis: bool = True,
    axes: Literal["all", "bottom"] | bool = "all",
    xtlab: Sequence | None = None,
    xtstep: int | None = None,
    tick_last: bool | None = None,
    cex_axis: float = 1,
    with_legend: Literal["auto", "right"] | bool = "auto",
    ltext: Sequence[str] | None = None,
    cex_legend: float = 1,
    legend_prop: float | None = None,
    rows: int | None = None,
    cols: int | None = None,
    conf_level: float | Sequence[float] | None = 0.95,
    lwd: float = 1,
    only_levels: Sequence[str] | None = None,
    figsize=None,
    ax=None,
    **style,
):
    """Plot implicative statistics by group level and sequence position.

    TraMineRextras stores signed residuals and negates them for plotting. This
    renderer applies that transformation to a local copy, replacing negative
    plotted values with ``NaN`` exactly as the R method does.
    """
    stored = _validate_result(x)
    selected = _select_levels(x, only_levels)
    nrows, ncols = _layout(len(selected), rows, cols)
    fig, panels = _resolve_axes(len(selected), nrows, ncols, figsize, ax)

    if axes not in {"all", "bottom", True, False}:
        raise ValueError("axes must be 'all', 'bottom', True, or False")
    if not isinstance(xaxis, bool) or not isinstance(yaxis, bool):
        raise TypeError("xaxis and yaxis must be booleans")
    use_step = x.xtstep if xtstep is None else xtstep
    if isinstance(use_step, bool) or not isinstance(use_step, Integral) or use_step <= 0:
        raise ValueError("xtstep must be a positive integer")
    use_tick_last = x.tick_last if tick_last is None else bool(tick_last)
    tick_labels = list(x.positions if xtlab is None else xtlab)
    if len(tick_labels) != len(x.positions):
        raise ValueError("xtlab must provide one label per sequence position")
    tick_locations, visible_tick_labels = _tick_spec(tick_labels, int(use_step), use_tick_last)

    try:
        axis_scale = float(cex_axis)
        legend_scale = float(cex_legend)
        line_width = float(lwd)
    except (TypeError, ValueError) as exc:
        raise ValueError("cex_axis, cex_legend, and lwd must be positive numbers") from exc
    if not all(
        np.isfinite(value) and value > 0 for value in (axis_scale, legend_scale, line_width)
    ):
        raise ValueError("cex_axis, cex_legend, and lwd must be positive numbers")
    if legend_prop is not None:
        if not isinstance(legend_prop, Real) or not 0 < float(legend_prop) < 1:
            raise ValueError("legend_prop must be strictly between zero and one")

    labels = list(x.alphabet if x.labels is None else x.labels)
    if ltext is not None:
        labels = [str(value) for value in ltext]
        if len(labels) != len(x.alphabet):
            raise ValueError("ltext must provide one label per alphabet state")
    colors = x.cpal
    if colors is not None and len(colors) != len(x.alphabet):
        raise ValueError("SeqImplic cpal must provide one color per alphabet state")
    if colors is None:
        # The shared default state palette, so the legend shows the lines'
        # colours (2026-09-22 audit S6-10).
        state_cmap = ensure_cmap(x, None)
        colors = [_state_color(state_cmap, code) for code in range(len(x.alphabet))]
    if with_legend not in {"auto", "right", True, False}:
        raise ValueError("with_legend must be 'auto', 'right', True, or False")
    for reserved in ("color", "label", "linewidth"):
        if reserved in style:
            raise ValueError(f"{reserved} must be controlled by cpal, ltext, or lwd")

    plotted = -stored.copy()
    plotted[plotted < 0] = np.nan
    limits = _resolve_ylim(plotted[selected], x.type, ylim)
    thresholds = _confidence_levels(conf_level, x.type)
    selected_levels = [x.levels[index] for index in selected]
    titles = _resolve_titles(main, selected_levels, x.levels)
    x_values = np.arange(1, len(x.positions) + 1)
    font_size = float(plt.rcParams["font.size"]) * axis_scale

    for panel_index, (axis, level_index, title) in enumerate(zip(panels, selected, titles)):
        for state_index, label in enumerate(labels):
            line_options = dict(style)
            line_options["color"] = colors[state_index]
            axis.plot(
                x_values,
                plotted[level_index, state_index],
                label=label,
                linewidth=line_width,
                **line_options,
            )
        for threshold in thresholds:
            axis.axhline(
                threshold,
                color="0.45",
                linestyle=":",
                linewidth=line_width,
            )
        axis.set_title(title)
        axis.set_ylim(limits)
        axis.set_xlim(0.5, len(x.positions) + 0.5)
        axis.set_ylabel("" if ylab is None else ylab)

        show_x = xaxis and axes is not False
        if axes == "bottom":
            show_x = show_x and panel_index // ncols == nrows - 1
        if show_x:
            axis.set_xticks(tick_locations, visible_tick_labels)
        else:
            axis.set_xticks([])
        if not yaxis:
            axis.set_yticks([])
        axis.tick_params(axis="both", labelsize=font_size)

    show_legend = with_legend is True or with_legend in {"auto", "right"}
    if show_legend:
        handles = [
            Line2D(
                [],
                [],
                color=colors[index],
                linewidth=line_width,
                label=label,
            )
            for index, label in enumerate(labels)
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
            location = "center right" if with_legend == "right" else "upper right"
        _legend_with_labels(
            fig,
            handles,
            labels,
            loc=location,
            fontsize=font_size * legend_scale,
            **legend_options,
        )

    returned_axes = panels[0] if len(panels) == 1 else np.asarray(panels, dtype=object)
    return fig, returned_axes


__all__ = ["plot_seqimplic"]
