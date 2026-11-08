"""The seqplot(type=...) dispatcher -- mirrors TraMineR's seqplot.R type routing.

A single entry point that routes a state-sequence object to the per-type renderer,
exactly as R's ``seqplot(seqdata, type=...)`` dispatches to plot.stslist[.*]
(SOURCES/TraMineR/R/seqplot.R). BUILD-from-R-design; provenance docstring-only.
"""

from __future__ import annotations

import functools
import inspect
import math
from typing import Any, Callable, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ...sequences import _require_stslist
from ..legend import _draw_state_legend
from ..series import _factor_panels
from ..utils import _close_new_figures_on_error
from .distribution import state_distribution_plot
from .events_viz import seqpcplot
from .frequency import seqfplot
from .index import sequence_index_plot
from .representatives import seqrfplot, seqrplot
from .summary import seqdHplot, seqHtplot, seqmsplot, seqmtplot

# type -> renderer, as R's seqplot routes them (pc included: seqplot.R:91).
_DISPATCH: dict[str, Callable[..., Any]] = {
    "i": sequence_index_plot,
    # Deliberate divergence: R's type="i" defaults to the first ten sequences
    # (idxs=1:10) while type="I" forces idxs=0 (all); traminex renders all
    # sequences for both spellings (sampling/subsetting is explicit via the
    # renderer's own options). docs/PLOTTING.md records this.
    "I": sequence_index_plot,
    "d": state_distribution_plot,
    "Ht": seqHtplot,
    "dH": seqdHplot,
    "f": seqfplot,
    "mt": seqmtplot,
    "ms": seqmsplot,
    "r": seqrplot,
    "rf": seqrfplot,
    # R routes type="pc" to seqpcplot, which draws a state sequence object through
    # its DSS (seqplot.R:91; TX-AUD-20261011-027).
    "pc": seqpcplot,
}

# Types whose plot carries no state legend in R's seqplot (seqplot.R:210).
_NO_STATE_LEGEND = {"Ht", "pc"}

_UNSUPPORTED_SOURCE_ALIASES = {
    "title",
    "cex_plot",
    "withlegend",
    "axes",
    "tlim",
    "ylas",
}


def _panel_titles(levels: Sequence[str], main: str | Sequence[str] | None) -> list[str]:
    """R's seqplot titles (seqplot.R:193-198) over factor()'s levels."""
    if main is None:
        return ["" for _ in levels]
    if isinstance(main, str):
        return list(levels) if main == "auto" else [f"{main} - {level}" for level in levels]
    titles = [str(title) for title in main]
    if len(titles) != len(levels):
        raise ValueError("main must be one title or one title per group level")
    return titles


def _palette_renderer(type_, renderer, kwargs):
    """Hand R's ``cpal`` to the renderer's own palette argument (seqplot.R:282-301).

    The index and distribution renderers take ``cmap`` (still accepted); the
    distribution-and-entropy plot is drawn by ``state_distribution_plot``'s
    ``"dH"`` mode, which ``seqdHplot`` delegates to, so it can take the palette;
    R drops ``cpal`` for the entropy plot (Ht). Returns the renderer to call.
    """
    if "cpal" not in kwargs or "cpal" in inspect.signature(renderer).parameters:
        return renderer
    cpal = kwargs.pop("cpal")
    if type_ == "Ht":
        return renderer
    if type_ == "dH":
        renderer = functools.partial(state_distribution_plot, mode="dH")
    if kwargs.get("cmap") is not None:
        raise TypeError("cpal and cmap cannot be given together")
    kwargs["cmap"] = cpal
    return renderer


def _panel_labels(ylab, count: int) -> list:
    """R's ``ylab``: one label for every panel, or one per group level (seqplot.R:173-176)."""
    if ylab is None or isinstance(ylab, str):
        return [ylab] * count
    labels = list(ylab)
    if len(labels) != count:
        raise ValueError("If a sequence, ylab must have one value per group level")
    return labels


def _legends(fig) -> list:
    return list(fig.legends) + [ax.get_legend() for ax in fig.axes if ax.get_legend()]


def _state_legend(
    seqdata, fig, axes, with_legend, palette, *, grouped: bool, before: list, own_figure: bool
) -> None:
    """R's seqplot legend: one state legend for the plot, or none (seqplot.R:426-456).

    Only the legends this call drew are looked at: ``before`` lists those a
    caller-given figure already had, which are left alone, and a figure-level
    legend is drawn only in a figure seqplot made (batch-13 review I-2).
    """
    known = {id(legend) for legend in before}
    drawn = [legend for legend in _legends(fig) if id(legend) not in known]
    if not with_legend or grouped:
        for legend in drawn:
            legend.remove()
        drawn = []
    if with_legend and not drawn:
        target = fig if own_figure and (grouped or len(axes) > 1) else axes[0]
        _draw_state_legend(seqdata, target, palette)


def _whole_set_rep_options(seqdata, kwargs):
    """R's type "r" preamble on the whole set, before the panels (seqplot.R:96-124).

    The distances are computed from the distance keywords when ``diss`` is not
    given, and ``dmax`` is ``max(diss)`` unless given, so that every panel uses
    the whole set's distances and maximum distance.
    """
    from ..._deprecated import MISSING, resolve_deprecated_alias
    from ...cluster import seqrep
    from ...distances import seqdist
    from .representatives import seqrplot

    options = dict(kwargs)
    options["diss"] = resolve_deprecated_alias(
        options.get("diss"),
        options.pop("dist_matrix", MISSING),
        canonical_name="diss",
        alias_name="dist_matrix",
        function_name="seqrep",
        default=None,
    )
    if options.get("diss") is None:
        named = set(inspect.signature(seqrplot).parameters) | set(
            inspect.signature(seqrep).parameters
        )
        distance = {key: options.pop(key) for key in list(options) if key not in named}
        options.pop("diss", None)
        options["diss"] = seqdist(seqdata, **distance)
    options["diss"] = np.asarray(options["diss"], dtype=float)
    if options.get("dmax") is None:
        options["dmax"] = float(np.max(options["diss"]))
    return options


def _whole_set_rf_ranks(seqdata, kwargs) -> np.ndarray:
    """R's type "rf" sort variable: the one-dimensional MDS of the whole set (seqplot.R:127-156).

    seqrf only orders the cases by ``sortv`` (dissrf.R:94), so the panels get the
    ranks of the whole-set axis: the same order, with the missing coordinate of
    a zero-weight case last, as R's ``order()`` puts ``wcmdscale``'s ``NA``.
    """
    from ...repseq_rf import _normalize_weights, _resolve_axis, _validate_dissimilarities

    matrix = _validate_dissimilarities(kwargs["diss"])
    weights = kwargs.get("weights")
    if weights is None:
        weights = seqdata.weights
    if not kwargs.get("weighted", True):
        weights = None
    resolved, weighted = _normalize_weights(matrix.shape[0], weights)
    axis = _resolve_axis(
        matrix,
        "mds",
        weights=resolved,
        weighted=weighted,
        squared=bool(kwargs.get("squared", False)),
    )
    ranks = np.empty(matrix.shape[0], dtype=float)
    ranks[np.argsort(axis, kind="stable")] = np.arange(matrix.shape[0], dtype=float)
    return ranks


def _group_panels(renderer, type_, seqdata, group, main, kwargs):
    """One panel per level of ``factor(group)``, as R's seqplot draws them."""
    if kwargs.get("ax") is not None:
        raise TypeError("seqplot draws its own panels with group; ax cannot be given")
    kwargs.pop("ax", None)
    if len(group) != seqdata.nseq:
        raise ValueError("group must contain one value for each row in the sequence object")
    rows, levels = _factor_panels(group)
    titles = _panel_titles(levels, main)
    # TraMineR.setlayout: min(nplot, 2) columns, rows filled by row.
    columns = max(1, min(len(levels), 2))
    lines = max(1, math.ceil(len(levels) / columns))
    figsize = kwargs.pop("figsize", None) or (6.0 * columns, 3.5 * lines)
    if type_ == "r":
        if any(len(index) < 2 for index in rows):
            raise ValueError(
                "For type = 'r', each group must have 2 or more cases. "
                "At least one group has only 1."
            )
        kwargs = _whole_set_rep_options(seqdata, kwargs)
    if type_ == "rf":
        # With group, R draws each level's medoid panel (plot.seqrf.R:4) and stops
        # for which.plot = "both" (seqplot.R:53-56).
        if kwargs.get("which_plot") == "both":
            raise ValueError('which_plot="both" applies only when group is None')
        if kwargs.get("which_plot") is None:
            kwargs["which_plot"] = "medoids"
        sortv = kwargs.get("sortv", "mds")
        if isinstance(sortv, str) and sortv == "mds" and kwargs.get("diss") is not None:
            kwargs["sortv"] = _whole_set_rf_ranks(seqdata, kwargs)
    fig, grid = plt.subplots(
        lines, columns, figsize=figsize, squeeze=False, constrained_layout=True
    )
    axes = list(grid.ravel())
    for unused in axes[len(levels) :]:
        unused.set_visible(False)
    from .summary import _subset_stslist

    # A ylab still here is rf's (seqrfplot labels its medoid panel itself); R gives
    # each panel its own ylab[np] (seqplot.R:173-176, 283; batch-13 review M-1).
    ylab = kwargs.get("ylab")
    panel_ylabs = None if ylab is None else _panel_labels(ylab, len(levels))
    for panel, (ax, index, title) in enumerate(zip(axes, rows, titles)):
        options = dict(kwargs)
        if panel_ylabs is not None:
            options["ylab"] = panel_ylabs[panel]
        # R subsets the dissimilarities and a sortv vector per panel (seqplot.R:317-369);
        # a weights vector is split the same way.
        if options.get("diss") is not None:
            options["diss"] = np.asarray(options["diss"], dtype=float)[np.ix_(index, index)]
        for name in ("sortv", "weights"):
            values = options.get(name)
            if values is not None and not isinstance(values, str) and len(values) == seqdata.nseq:
                if name == "sortv" and isinstance(
                    getattr(values, "dtype", None), pd.CategoricalDtype
                ):
                    # A categorical sortv keeps its levels in each panel, so the panel is
                    # ordered by level codes as R's order() orders a factor
                    # (TX-AUD-20261019-019).
                    options[name] = (
                        values.iloc[index] if isinstance(values, pd.Series) else values[index]
                    )
                else:
                    options[name] = np.asarray(values)[index]
        renderer(_subset_stslist(seqdata, index), ax=ax, **options)
        ax.set_title(title)
    return fig, axes[: len(levels)]


@_close_new_figures_on_error
@_require_stslist
def seqplot(seqdata, *, type="i", group=None, main="auto", **kwargs):
    """Dispatch a state-sequence plot by ``type`` -- re-expresses TraMineR ``seqplot``.

    Routes to the per-type renderer, mirroring ``seqplot.R``:
    ``i``/``I`` index, ``d`` distribution, ``Ht`` entropy, ``dH`` distribution +
    entropy, ``f`` frequency (coloured bands), ``mt`` mean time, ``ms`` modal state,
    ``r`` representatives, ``rf`` relative-frequency groups, ``pc`` the
    parallel-coordinate plot of the state sequences' DSS (``seqpcplot``). Extra
    keyword arguments pass through to the renderer (e.g. ``diss=`` for
    ``type="r"``/``"rf"``). As in R, ``weighted`` is honoured by
    every type: ``weighted=False`` plots the unweighted series.

    With ``group``, one panel per level of ``factor(group)`` (missing values in no
    panel), each the plot of that level's sequences, titled as R's seqplot does:
    by the level (``main="auto"``), by ``"<main> - <level>"`` for one ``main``, or
    by ``main`` itself when it gives one title per level (seqplot.R:163-205;
    TX-AUD-20261010-044). R's internal group() puts a space before each title and
    orders numeric groups as text; here the levels are factor()'s, as for
    seqHtplot. ``type="Ht"`` keeps its own grouped rendering (one line per
    level). Each panel gets its subset of ``diss`` and of a ``sortv`` vector.
    As in R (seqplot.R:96-156), types ``"r"`` and ``"rf"`` first compute on the
    whole set what the panels share: for ``"r"`` the distances (from the
    distance keywords when ``diss`` is not given) and ``dmax = max(diss)``, and
    for ``"rf"`` the default ``sortv="mds"`` axis; ``"r"`` refuses a level with
    fewer than two sequences. A per-sequence ``weights`` vector is split per
    panel. Without ``group``, a ``main`` other than ``"auto"`` titles the plot.
    ``seqfplot``, ``seqmtplot``, ``seqmsplot``, ``seqdHplot``, ``seqrplot`` and
    ``seqrfplot`` take ``group`` and ``main`` the same way.

    As R's seqplot (seqplot.R:6-11, 210, 282-301), every type takes ``cpal`` (one
    colour per state, handed to the renderer's palette argument; ignored for
    ``"Ht"``, as R ignores it), ``ylab`` (one label for every panel, or one per
    group level; for ``"rf"`` the medoid panel's) and ``with_legend``: ``"auto"``
    (default) or True draws one state legend for the figure, False none; types
    ``"Ht"`` and ``"pc"`` get no state legend (TX-AUD-20261011-067).

    Returns
    -------
    Whatever the routed renderer returns (typically ``(fig, ax)``); with
    ``group``, ``(fig, axes)`` with one axes per level.
    """
    unsupported = sorted(_UNSUPPORTED_SOURCE_ALIASES.intersection(kwargs))
    if unsupported:
        raise TypeError(
            "seqplot does not support source alias(es) "
            f"{unsupported}; its type-specific Python renderers have no shared "
            "canonical target for those plotting controls"
        )
    try:
        renderer = _DISPATCH[type]
    except KeyError:
        raise ValueError(
            f"unknown seqplot type {type!r}; expected one of {sorted(_DISPATCH)} "
            "(event sequences: call seqpcplot directly)"
        ) from None
    if "weighted" in kwargs and "weighted" not in inspect.signature(renderer).parameters:
        # R's seqplot hands `weighted` to seqstatd (d, Ht, dH) and seqtab (f),
        # whose renderers here read the attached weights: weighted=False plots
        # an unweighted copy (2026-09-22 audit S6-20).
        if not kwargs.pop("weighted"):
            seqdata = seqdata.copy()
            seqdata.weights = None
    # R's options for every type (seqplot.R:6-11): cpal, ylab and with.legend.
    with_legend = kwargs.pop("with_legend", "auto")
    if with_legend not in (True, False, "auto"):
        raise ValueError("with_legend must be True, False or 'auto'")
    renderer = _palette_renderer(type, renderer, kwargs)
    palette = kwargs.get("cpal", kwargs.get("cmap"))
    # seqrfplot labels its medoid panel itself and leaves the other unlabelled, as R.
    ylab = None if type == "rf" else kwargs.pop("ylab", None)
    grouped = group is not None and type != "Ht"
    supplied = kwargs.get("ax")
    if isinstance(supplied, (list, tuple)):
        supplied = supplied[0] if supplied else None
    before = _legends(supplied.figure) if supplied is not None else []
    if grouped:
        result = _group_panels(renderer, type, seqdata, group, main, kwargs)
    else:
        if group is not None:
            kwargs["group"] = group
        result = renderer(seqdata, **kwargs)
    if not isinstance(result, tuple):
        return result
    fig, axes = result[0], result[1]
    axes_list = list(axes) if isinstance(axes, (list, tuple)) else [axes]
    if not grouped and isinstance(main, str) and main != "auto":
        for ax in axes_list:
            ax.set_title(main)
    if ylab is not None:
        for ax, label in zip(axes_list, _panel_labels(ylab, len(axes_list))):
            ax.set_ylabel(label)
    if type not in _NO_STATE_LEGEND:
        _state_legend(
            seqdata,
            fig,
            axes_list,
            with_legend,
            palette,
            grouped=grouped,
            before=before,
            own_figure=supplied is None,
        )
    return result
