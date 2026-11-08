from __future__ import annotations

from typing import Optional

import numpy as np
from matplotlib import rcParams
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import text_to_path
from matplotlib.ticker import Formatter


class _RTimeLabels(Formatter):
    """Time-axis tick labels drawn as R's ``axis()`` draws them (TX-AUD-20261011-067).

    R gives every tick its label, then draws the labels left to right and omits
    one that would come closer than the width of an "m" to the last label drawn
    (R's ``?axis``: ``gap.axis = 1`` for labels parallel to the axis). The same
    rule runs each time the axis is drawn, on its current size, so the labels
    shown follow the figure's size as R's follow the device's.
    """

    def __init__(self, labels: dict[float, str], fontsize=None):
        self.labels = labels
        self.fontsize = fontsize
        self._widths: dict[tuple[str, float], float] = {}

    def __call__(self, x, pos=None):
        return self.labels.get(float(x), "")

    def _width(self, text: str, size: float) -> float:
        # Width in points of ``text`` in the tick-label font.
        key = (text, size)
        if key not in self._widths:
            prop = FontProperties(size=size)
            self._widths[key] = text_to_path.get_text_width_height_descent(text, prop, False)[0]
        return self._widths[key]

    def format_ticks(self, values):
        texts = [self(value) for value in values]
        axes = getattr(self.axis, "axes", None)
        if axes is None or not texts:
            return texts
        size = FontProperties(
            size=self.fontsize if self.fontsize is not None else rcParams["xtick.labelsize"]
        ).get_size_in_points()
        scale = axes.figure.dpi / 72.0
        gap = self._width("m", size) * scale
        points = np.column_stack([np.asarray(values, dtype=float), np.zeros(len(values))])
        pixels = axes.transData.transform(points)[:, 0]
        shown = [""] * len(texts)
        last = -np.inf
        for index in np.argsort(pixels, kind="stable"):
            text = texts[index]
            if not text:
                continue
            half = self._width(text, size) * scale / 2.0
            if pixels[index] - half - last >= gap:
                shown[index] = text
                last = pixels[index] + half
        return shown


def _time_axis(
    ax,
    x,
    L: int,
    *,
    offset: float = 0.0,
    label_color: Optional[str] = None,
    label_fontsize: Optional[float] = None,
) -> None:
    """R's time ticks and labels, the tick of position t (0-based) at ``t + offset``."""
    step = max(int(getattr(x, "xtstep", 1)), 1)
    start = int(getattr(x, "start", 1))

    # R's positions seq(1, L, by = xtstep), plus L when tick.last (0-based here)
    ticks = list(range(0, max(L, 0), step))
    if ticks and getattr(x, "tick_last", False) and ticks[-1] != L - 1:
        ticks.append(L - 1)

    # Build full label vector for all positions 0..L-1
    xtlab = getattr(x, "xtlab", None)
    if isinstance(xtlab, (list, tuple)) and len(xtlab) >= L:
        labels = [str(xtlab[i]) for i in range(L)]
    else:
        labels = [str(start + i) for i in range(L)]

    positions = [float(t) + offset for t in ticks]
    ax.set_xticks(positions)
    ax.xaxis.set_major_formatter(
        _RTimeLabels({p: labels[t] for p, t in zip(positions, ticks)}, label_fontsize)
    )
    style: dict[str, object] = {}
    if label_color is not None:
        style["labelcolor"] = label_color
    if label_fontsize is not None:
        style["labelsize"] = label_fontsize
    if style:
        ax.tick_params(axis="x", **style)


def apply_time_axis(
    ax,
    x,
    L: int,
    *,
    one_based: bool = False,
    max_labels: int = 10,
    label_color: Optional[str] = None,
    label_fontsize: Optional[float] = None,
):
    """Apply R's time ticks and labels to an axis.

    As TraMineR's plots (plot.stslist.R:160-161): a tick at every ``x.xtstep``-th
    position from the first, plus the last when ``x.tick_last``, labelled with
    ``x.xtlab`` (else the positions numbered from ``x.start``). As R's ``axis()``,
    the labels are drawn left to right and one that would come closer than the
    width of an "m" to the last label drawn is omitted, on the axis' size when
    it is drawn (TX-AUD-20261011-067). ``max_labels`` is kept so that existing
    calls work; it no longer limits the labels. If ``one_based`` is True, tick
    positions are shifted by +1 (useful for bar plots at 1..L).
    """
    _time_axis(
        ax,
        x,
        L,
        offset=1.0 if one_based else 0.0,
        label_color=label_color,
        label_fontsize=label_fontsize,
    )


def _apply_time_axis_on_cells(ax, x, L: int) -> None:
    """:func:`apply_time_axis` for plots that draw position t (0-based) on the cell [t, t + 1].

    The ticks sit at the cell centres t + 0.5, as R puts the tick of position t
    at t - 0.5, the centre of its cell [t - 1, t] (plot.stslist.freq.R:91-92,
    plot.stslist.modst.R:94, plot.stslist.rep.R:144; TX-AUD-20261011-067).
    """
    _time_axis(ax, x, L, offset=0.5)
