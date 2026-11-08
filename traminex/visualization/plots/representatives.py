"""Representative-sequence plots over the R-golden'd series (BUILD-from-R-design).

seqrplot / seqrfplot mirror TraMineR's plot.stslist.rep and plot.seqrf drawing
routines (SOURCES/TraMineR/R/plot.stslist.rep.R, plot.seqrf.R;
SOURCES/TraMineRextras/R/seqplot-rf.R) over traminex's already-pinned series
(seqrep / dissrep, seqrf / dissrf). The matplotlib code is original (no Sequenzo
borrow), so provenance is docstring-only.
"""

from __future__ import annotations

import sys
import warnings
from collections.abc import Mapping
from numbers import Integral
from typing import Any, Literal

import matplotlib.pyplot as plt
import numpy as np

from ..._rformat import r_as_character, r_number_text, r_round
from ...cluster import seqrep
from ...diss import disscenter
from ...repseq_rf import seqrf
from ...sequences import _require_stslist
from ..palette import _state_color, ensure_cmap
from ..ticks import _apply_time_axis_on_cells
from ..utils import _close_new_figures_on_error


def _draw_state_bands(ax, code_rows, heights, cm, *, border=True, artist_style=None):
    """Stack coloured state-segment bands bottom-up; band i has thickness heights[i].

    Returns the y-centres of the bands (for tick placement)."""
    from matplotlib.patches import Rectangle

    if artist_style is not None and not isinstance(artist_style, Mapping):
        raise TypeError("artist_style must be a mapping or None")
    patch_style = dict(artist_style or {})
    forbidden = {"color", "facecolor"}.intersection(patch_style)
    if forbidden:
        raise ValueError(f"use cpal instead of artist style keys {sorted(forbidden)}")
    patch_style.setdefault("edgecolor", "white" if border else "none")
    patch_style.setdefault("linewidth", 0.3 if border else 0.0)

    y0 = 0.0
    centres = []
    L = code_rows.shape[1] if len(code_rows) else 0
    for codes, h in zip(code_rows, heights):
        for t in range(L):
            ax.add_patch(
                Rectangle(
                    (t, y0),
                    1.0,
                    h,
                    facecolor=_state_color(cm, int(codes[t])),
                    **patch_style,
                )
            )
        centres.append(y0 + h / 2.0)
        y0 += h
    return centres, y0


def _wtd_fivenum(x, weights):
    """Weighted Tukey five-number summary -- port of TraMineR ``wtd.fivenum.tmr``.

    SOURCES/TraMineR/R/TraMineR-wtd-stats.R (TraMineR's vendored copy of ENmisc's
    weighted boxplot, by Erich Neuwirth). Returns ``[min, Q1, median, Q3, max]`` at
    the weighted hinge positions; for unequal weights the hinges use the
    ``interpolatedindex`` weighted-rank rule, exactly as R draws ``plot.seqrf``'s
    right panel. GPL port; provenance docstring-only.
    """
    x = np.asarray(x, dtype=float)
    weights = np.ones_like(x) if weights is None else np.asarray(weights, dtype=float)
    equalweights = bool(np.all(np.diff(weights) == 0)) if x.size > 1 else True
    keep = ~(np.isnan(x) | (weights == 0))
    x, weights = x[keep], weights[keep]
    if x.size == 0:
        return np.full(5, np.nan)
    order = np.argsort(x, kind="stable")
    x, weights = x[order], weights[order]
    n = float(weights.sum())
    if n == 0:
        return np.full(5, np.nan)
    if equalweights:
        half = 0.5 * np.floor(0.5 * (n + 3))
        d = np.array([1.0, half, 0.5 * (n + 1), n + 1 - half, n])
    elif x.size > 1:

        def interpolatedindex(myval):
            indices = np.arange(1, weights.size + 1)
            weightsbelow = np.concatenate([[0.0], np.cumsum(weights)[:-1]])
            weightsabove = n - weightsbelow - weights
            lowcands = weightsbelow < myval
            highcands = weightsabove < (n - myval)
            lo = int(indices[lowcands].max()) if lowcands.any() else 1
            hi = int(indices[highcands].min()) if highcands.any() else x.size
            return (lo + hi) / 2.0

        d = np.array(
            [1.0] + [interpolatedindex(q) for q in (0.25 * n, 0.5 * n, 0.75 * n)] + [float(x.size)]
        )
    else:
        d = np.ones(5)
    lo_idx = np.clip(np.floor(d).astype(int) - 1, 0, x.size - 1)
    hi_idx = np.clip(np.ceil(d).astype(int) - 1, 0, x.size - 1)
    return 0.5 * (x[lo_idx] + x[hi_idx])


def _wtd_boxplot_stats(x, weights, coef: float = 1.5):
    """Weighted boxplot statistics -- port of TraMineR ``wtd.boxplot.stats.tmr``.

    Returns the matplotlib ``bxp`` dict (``whislo``/``q1``/``med``/``q3``/``whishi``/
    ``fliers``) for one group, with the box hinges from :func:`_wtd_fivenum` and the
    whiskers shrunk to the ``coef*IQR`` non-outlier range -- the box R's
    ``plot.seqrf`` draws via ``wtd.boxplot.tmr``.
    """
    x = np.asarray(x, dtype=float)
    weights = None if weights is None else np.asarray(weights, dtype=float)
    nna = ~np.isnan(x)
    if x.size > 1:
        stats = _wtd_fivenum(x, weights)
    elif x.size == 1:
        stats = np.repeat(x, 5).astype(float)
    else:
        stats = np.full(5, np.nan)
    iqr = stats[3] - stats[1]
    out_mask = np.zeros(x.size, dtype=bool)
    if coef != 0 and np.isfinite(iqr):
        out_mask = (x < (stats[1] - coef * iqr)) | (x > (stats[3] + coef * iqr))
        if out_mask[nna].any():
            good = x[~out_mask]
            stats[0], stats[4] = good.min(), good.max()
    return {
        "label": "",
        "whislo": float(stats[0]),
        "q1": float(stats[1]),
        "med": float(stats[2]),
        "q3": float(stats[3]),
        "whishi": float(stats[4]),
        "fliers": [float(v) for v in x[out_mask & nna]],
    }


@_close_new_figures_on_error
@_require_stslist
def seqrplot(
    x,
    *,
    diss=None,
    criterion: Literal["density", "dist", "prob", "mscore"] = "density",
    coverage: float = 0.25,
    nrep=None,
    weighted: bool = True,
    cpal=None,
    border: bool = True,
    figsize=(8, 4),
    ax=None,
    group=None,
    main="auto",
    **dist_kwargs,
):
    """Representative sequences -- re-expresses TraMineR ``seqrplot`` / ``plot.stslist.rep``.

    Draws the selected representative sequences as horizontal coloured
    state-segment bands, each band's thickness proportional to the share of
    sequences it represents (``na(%)`` from ``seqrep``/``dissrep``) -- R's
    ``barplot(seqgbar(rep), width = Statistics$"na(%)"/100, horiz=TRUE)``
    (SOURCES/TraMineR/R/plot.stslist.rep.R). The representative set, the
    representative sequences and the coverage statistics come from ``seqrep``
    (R-golden'd in ``dissrep_*`` / ``seqrep_*``). R's optional mean-distance /
    discrepancy symbol overlay (``stats=TRUE``) is not drawn. As R, the bands carry
    no y labels and the y axis is titled "N representative(s) (n=N)" ("weighted
    n=N" for weighted sequences, N rounded to 2 decimals).

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
    """
    if group is not None or main != "auto":
        from .dispatch import seqplot

        # Called directly the renderer draws no legend; main and group keep it so.
        return seqplot(
            x,
            type="r",
            with_legend=False,
            group=group,
            main=main,
            diss=diss,
            criterion=criterion,
            coverage=coverage,
            nrep=nrep,
            weighted=weighted,
            cpal=cpal,
            border=border,
            figsize=figsize,
            ax=ax,
            **dist_kwargs,
        )
    if "cex_plot" in dist_kwargs:
        raise TypeError(
            "seqrplot does not support source alias 'cex_plot'; the Python renderer "
            "has no cex_with_axis plotting control"
        )
    if x.nseq < 2:
        # R's seqplot counts the whole set as one group (seqplot.R:97-101).
        raise ValueError(
            "For type = 'r', each group must have 2 or more cases. At least one group has only 1."
        )
    # Both marker representations are equivalent inputs (identity for
    # in-alphabet objects; cycle-9 audit N-04/N-05).
    x = x._canonicalized()
    rep = seqrep(
        x,
        diss=diss,
        criterion=criterion,
        coverage=coverage,
        nrep=nrep,
        weighted=weighted,
        **dist_kwargs,
    )
    codes = np.asarray(rep["representatives"], dtype=int)
    stats = rep["statistics"]
    na_pct = stats["na(%)"].to_numpy(dtype=float)[: len(codes)]
    heights = na_pct / 100.0

    cm = ensure_cmap(x, cpal)
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    else:
        fig = ax.figure
    _, ytop = _draw_state_bands(ax, codes, heights, cm, border=border)
    L = x.data.shape[1]
    ax.set_xlim(0, L)
    ax.set_ylim(0, max(ytop, 1e-9))
    # R's y axis: no label on the bands, titled "N representative(s) (n=N)"
    # ("weighted n=N"; plot.stslist.rep.R:109-110, 132-136; TX-REP-20261023-001).
    ax.set_yticks([])
    _apply_time_axis_on_cells(ax, x, L)
    ax.set_xlabel("Time")
    metadata = rep["metadata"]
    weighted_text = "weighted " if metadata.get("weighted") else ""
    n = r_as_character(r_round(float(metadata["weights_sum"]), 2))
    ax.set_ylabel(f"{len(codes)} representative(s) ({weighted_text}n={n})")
    ax.set_title("Representative sequences")
    return fig, ax


@_close_new_figures_on_error
@_require_stslist
def seqrfplot(
    x,
    diss,
    *,
    k=None,
    sortv="mds",
    weighted: bool = True,
    cpal=None,
    border: bool = True,
    figsize=(10, 4),
    which_plot: Literal["both", "medoids", "diss.to.med"] | None = None,
    ylab: str | None = None,
    yaxis: bool = True,
    main: str | None = None,
    ax=None,
    artist_style: Mapping[str, object] | None = None,
    rf: Mapping[str, Any] | None = None,
    group=None,
    **rf_kwargs,
):
    """Relative-frequency representative plot -- re-expresses TraMineR ``seqrfplot`` / ``plot.seqrf``.

    Two panels (R plot.seqrf / TraMineRextras seqplot.rf):
    - LEFT: the medoid sequence of each relative-frequency group, drawn as coloured
      state-segment bands (an index plot of ``seqrf$seqtoplot``), band thickness
      proportional to the group ``heights``.
    - RIGHT: a horizontal **weighted** Tukey-hinge boxplot per group of the
      within-group distances to the group medoid (``dissrf$dist_list`` with
      ``weights_list``), box widths proportional to ``heights`` -- the same
      ``wtd.boxplot.tmr`` box R draws (hinges/median/whiskers reproduce R's weighted
      five-number summary, not matplotlib's unweighted ``np.percentile``).
    The grouping / medoids / heights / residual distances come from ``seqrf``
    (R-golden'd in ``seqrf_*`` / ``dissrf_*``). The figure subtitle reports the
    representation quality ``R2`` / ``F`` from ``dissrf``. ``artist_style``
    applies common matplotlib patch properties to the rendered bands and boxes;
    those properties never enter the analytical ``seqrf`` call. ``which_plot``
    left at ``None`` draws both panels; with ``group``, each level gets its
    medoid panel, as R draws grouped plots (plot.seqrf.R:4), and ``"both"`` is
    refused as in R (seqplot.R:53-56).

    Returns
    -------
    fig, (ax_left, ax_right)
    """
    if group is not None:
        from .dispatch import seqplot

        # Called directly the renderer draws no legend; group keeps it so.
        return seqplot(
            x,
            type="rf",
            with_legend=False,
            group=group,
            main="auto" if main is None else main,
            diss=diss,
            k=k,
            sortv=sortv,
            weighted=weighted,
            cpal=cpal,
            border=border,
            figsize=figsize,
            which_plot=which_plot,
            ylab=ylab,
            yaxis=yaxis,
            ax=ax,
            artist_style=artist_style,
            rf=rf,
            **rf_kwargs,
        )
    if which_plot is None:
        which_plot = "both"
    if which_plot not in {"both", "medoids", "diss.to.med"}:
        raise ValueError("which_plot must be one of 'both', 'medoids', or 'diss.to.med'")
    # Both marker representations are equivalent inputs (identity for
    # in-alphabet objects; cycle-9 audit N-04/N-05).
    x = x._canonicalized()
    if artist_style is not None and not isinstance(artist_style, Mapping):
        raise TypeError("artist_style must be a mapping or None")
    forbidden = {"color", "facecolor"}.intersection(artist_style or {})
    if forbidden:
        raise ValueError(f"use cpal instead of artist style keys {sorted(forbidden)}")

    if rf is None:
        # As in R's plot.seqrf, an already-fitted seqrf result may be supplied
        # (rf=...) instead of refitting here.
        rf = seqrf(x, diss, k=k, sortv=sortv, weighted=weighted, **rf_kwargs)
    codes = np.asarray(rf["seqtoplot"], dtype=int)
    heights = np.asarray(rf["heights"], dtype=float)
    dist_list = rf["rf"]["dist_list"]
    weights_list = rf["rf"]["weights_list"]

    if which_plot == "both":
        if ax is None:
            fig, axes = plt.subplots(1, 2, figsize=figsize, constrained_layout=True)
            axL, axR = axes
        else:
            try:
                axL, axR = tuple(ax)
            except (TypeError, ValueError):
                raise ValueError("ax must contain two axes when which_plot='both'") from None
            if axL.figure is not axR.figure:
                raise ValueError("both supplied axes must belong to the same figure")
            fig = axL.figure
    else:
        if ax is None:
            fig, single_ax = plt.subplots(figsize=figsize, constrained_layout=True)
        else:
            single_ax = ax
            fig = ax.figure
        axL = single_ax if which_plot == "medoids" else None
        axR = single_ax if which_plot == "diss.to.med" else None

    cm = ensure_cmap(x, cpal)
    centres = (np.cumsum(heights) - heights / 2.0).tolist()
    ytop = float(heights.sum())
    if axL is not None:
        centres, ytop = _draw_state_bands(
            axL,
            codes,
            heights,
            cm,
            border=border,
            artist_style=artist_style,
        )
        L = x.data.shape[1]
        axL.set_xlim(0, L)
        axL.set_ylim(0, max(ytop, 1e-9))
        if yaxis:
            axL.set_yticks(centres)
            axL.set_yticklabels([f"g{i + 1}" for i in range(len(codes))])
        else:
            axL.set_yticks([])
        _apply_time_axis_on_cells(axL, x, L)
        axL.set_xlabel("Time")
        if ylab is not None:
            axL.set_ylabel(ylab)
        axL.set_title(
            "Group medoids"
            if which_plot == "both"
            else ("Sequences medoids" if main is None else f"{main}: Sequences medoids")
        )

    if axR is not None:
        # One weighted Tukey-hinge box per group at its band centre. R draws the boxes on
        # a 0..k axis with widths boxwex * h / max(h), boxwex = 0.8 (plot.seqrf.R:79-90;
        # graphics::bxp); on this 0..sum(h) axis that is 0.8 * h / max(h) * sum(h) / k,
        # inside each band (TX-AUD-20261019-018).
        box_widths = (
            0.8 * heights / heights.max() * heights.sum() / len(heights)
            if heights.size and heights.max() > 0
            else heights
        )
        bxpstats = [
            _wtd_boxplot_stats(np.asarray(d, dtype=float), np.asarray(w, dtype=float))
            for d, w in zip(dist_list, weights_list)
        ]
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="vert: bool was deprecated.*")
            axR.bxp(
                bxpstats,
                positions=list(centres),
                widths=list(box_widths),
                vert=False,
                patch_artist=True,
                boxprops=dict(artist_style or {}),
            )
        axR.set_ylim(0, max(ytop, 1e-9))
        if not yaxis:
            axR.set_yticks([])
        else:
            axR.set_yticks(centres)
            axR.set_yticklabels([f"g{i + 1}" for i in range(len(codes))])
        if ylab is not None:
            axR.set_ylabel(ylab)
        axR.set_xlabel("Distance to medoid")
        axR.set_title(
            "Distances to medoids"
            if which_plot == "both"
            else (
                "Dissimilarities to medoids"
                if main is None
                else f"{main}: Dissimilarities to medoids"
            )
        )

    r2 = float(rf["rf"]["R2"])
    fstat = float(rf["rf"]["Fstat"])

    def _fixed(value: float) -> str:  # R's sprintf("%0.2f") spells NaN and Inf
        return f"{value:.2f}" if np.isfinite(value) else r_number_text(value)

    quality = f"Representation quality: R2={_fixed(r2)} and F={_fixed(fstat)}"
    if which_plot == "both":
        fig.suptitle(quality if main is None else f"{main}\n{quality}")
    return (fig, (axL, axR)) if which_plot == "both" else (fig, axL if axL is not None else axR)


@_close_new_figures_on_error
@_require_stslist
def seqplot_rf(
    seqdata,
    diss,
    *,
    k: int | None = None,
    sortv: object = None,
    ylab: str | None = None,
    yaxis: bool = False,
    main: str | None = None,
    which_plot: Literal["both", "medoids", "diss.to.med"] = "both",
    grp_meth: Literal["first", "random"] = "first",
    random_state: int | np.random.Generator | None = None,
    **style,
):
    """Explicit ``seqplot.rf`` route with rendering style isolated from ``seqrf``."""
    if which_plot not in {"both", "medoids", "diss.to.med"}:
        raise ValueError("which_plot must be one of 'both', 'medoids', or 'diss.to.med'")
    if grp_meth not in {"first", "random"}:
        raise ValueError("grp_meth must be one of 'first' or 'random'")

    # Both marker representations are equivalent inputs (identity for
    # in-alphabet objects; cycle-9 audit N-05).
    seqdata = seqdata._canonicalized()
    matrix = np.asarray(diss, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("diss must be a square dissimilarity matrix")
    if matrix.shape[0] != seqdata.nseq:
        raise ValueError("diss size must equal the number of sequences")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("diss must contain only finite values")

    if k is None:
        k = seqdata.nseq // 10
    if isinstance(k, bool) or not isinstance(k, Integral):
        raise TypeError("k must be an integer")
    k = int(k)
    if not 1 <= k <= seqdata.nseq:
        raise ValueError("k must be between 1 and the number of sequences")

    print(f" [>] Using k={k} frequency groups", file=sys.stderr)
    plot_options = {
        name: style.pop(name) for name in ("cpal", "border", "figsize", "ax") if name in style
    }
    if "artist_style" in style:
        raise ValueError("pass artist properties directly, not through artist_style")
    fitted = seqrf(
        seqdata,
        matrix,
        k=k,
        sortv="mds" if sortv is None else sortv,
        weighted=False,
        grp_meth=grp_meth,
        random_state=random_state,
    )
    # seqplot.rf's own statistic, from squared distances to the group medoids and
    # to the global medoid (seqplot-rf.R:22-24, 139-146), not seqrf's
    # (TX-AUD-20261010-041); it also goes into the figure's quality line.
    to_global = matrix[:, int(disscenter(matrix, medoids_index="first"))]
    if "mdsk" in fitted["rf"]:
        groups = len(np.unique(np.asarray(fitted["rf"]["mdsk"])))
        to_group = np.asarray(fitted["rf"]["kmedoid_dist"], dtype=float)
    else:  # k = 1: the one group's medoid is the global medoid, as in R
        groups, to_group = 1, to_global
    with np.errstate(divide="ignore", invalid="ignore"):
        # IEEE division as R's: one group or one sequence per group gives NaN/Inf.
        r2_value = 1.0 - np.float64(np.sum(to_group**2)) / np.float64(np.sum(to_global**2))
        explained = r2_value / np.float64(groups - 1)
        unexplained = (1.0 - r2_value) / np.float64(seqdata.nseq - groups)
        r2, fstat = float(r2_value), float(explained / unexplained)
    shown = {**fitted, "rf": {**fitted["rf"], "R2": r2, "Fstat": fstat}}
    result = seqrfplot(
        seqdata,
        matrix,
        k=k,
        sortv="mds" if sortv is None else sortv,
        weighted=False,
        which_plot=which_plot,
        ylab=ylab,
        yaxis=yaxis,
        main=main,
        artist_style=style,
        rf=shown,
        **plot_options,
    )
    print(f" [>] Pseudo/median-based-R2: {r_number_text(r2, 7)}", file=sys.stderr)
    print(f" [>] Pseudo/median-based-F statistic: {r_number_text(fstat, 7)}", file=sys.stderr)
    return result
