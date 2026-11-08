from __future__ import annotations

from collections.abc import Callable
from functools import wraps
from io import BytesIO
from math import ceil
from typing import Optional, ParamSpec, Sequence, Tuple, TypeVar

import matplotlib.pyplot as plt
import numpy as np

_P = ParamSpec("_P")
_R = TypeVar("_R")


def _close_new_figures_on_error(function: Callable[_P, _R]) -> Callable[_P, _R]:
    """Close the figures a plotting call opened when that call raises.

    Renderers create their figure before every argument or caller style has
    been applied, so a failing call used to leave an open figure in pyplot's
    registry (2026-09-22 audit S6-09). Figures that existed before the call,
    including caller-supplied ones, are left alone.
    """

    @wraps(function)
    def guarded(*args: _P.args, **kwargs: _P.kwargs) -> _R:
        before = set(plt.get_fignums())
        try:
            return function(*args, **kwargs)
        except BaseException:
            for number in set(plt.get_fignums()) - before:
                plt.close(number)
            raise

    return guarded


def _import_seaborn(feature: str):
    """Import seaborn for ``feature`` or raise an ``ImportError`` naming the extra.

    seaborn is imported inside the plotting functions that use it, so a missing
    seaborn never reaches the package facade's extras message (audit S7-16). The
    wording matches that message and adds the ``seaborn=False`` opt-out.
    """
    try:
        import seaborn
    except ImportError as exc:
        from .. import _OPTIONAL_DEPENDENCY_EXTRAS

        extra = _OPTIONAL_DEPENDENCY_EXTRAS["seaborn"]
        raise ImportError(
            f"seaborn is required for {feature}; install the '{extra}' extra "
            f"(pip install 'traminex[{extra}]') or pass seaborn=False"
        ) from exc
    return seaborn


def _reserve_legend_space(fig, axes, proportion: float, position: str):
    """Reserve an R ``TraMineR.setlayout``-style figure fraction for a legend."""
    if position not in {"bottom", "right"}:
        raise ValueError("legend position must be 'bottom' or 'right'")
    plot_fraction = 1.0 - float(proportion)
    fig.set_layout_engine(None)
    for axis in axes:
        x0, y0, width, height = axis.get_position().bounds
        if position == "right":
            axis.set_position((x0 * plot_fraction, y0, width * plot_fraction, height))
        else:
            axis.set_position((x0, proportion + y0 * plot_fraction, width, height * plot_fraction))
    if position == "right":
        return (plot_fraction, 0.0, float(proportion), 1.0)
    return (0.0, 0.0, 1.0, float(proportion))


def compute_even_ticks(
    length: int, *, max_labels: int = 10, include_endpoints: bool = True
) -> list[int]:
    """Return up to max_labels evenly spaced tick indices over [0, length-1].

    - Guarantees endpoints (0 and length-1) if include_endpoints is True and length>0.
    - For short sequences (<= max_labels), returns all indices.
    """
    L = int(max(0, length))
    if L == 0:
        return []
    if max_labels is None or L <= int(max_labels):
        return list(range(L))
    count = max(2, int(max_labels))
    gap = (L - 1) / float(count - 1)
    ticks = [int(round(i * gap)) for i in range(count)]
    ticks = sorted(set(max(0, min(L - 1, t)) for t in ticks))
    if include_endpoints:
        if ticks[0] != 0:
            ticks[0] = 0
        if ticks[-1] != L - 1:
            ticks[-1] = L - 1
    return ticks


def determine_layout(
    num_items: int,
    *,
    layout: str = "column",
    nrows: Optional[int] = None,
    ncols: Optional[int] = None,
    default_cols: int = 3,
) -> Tuple[int, int]:
    """Decide (nrows, ncols) for a grid of subplots.

    - If both nrows and ncols are provided, they are validated and returned.
    - layout="column" uses a fixed number of columns (default_cols) and enough rows.
    - layout="grid" approximates a square grid via sqrt.
    """
    num = int(max(0, num_items))
    if num == 0:
        return 0, 0

    if (nrows is not None) and (ncols is not None):
        if int(nrows) * int(ncols) < num:
            raise ValueError(f"Provided layout ({nrows}x{ncols}) is too small for {num} items.")
        return int(nrows), int(ncols)

    if layout == "column":
        c = int(max(1, default_cols))
        r = int(ceil(num / float(c)))
        return r, c
    if layout == "grid":
        c = int(max(1, np.ceil(np.sqrt(num))))
        r = int(ceil(num / float(c)))
        return r, c
    raise ValueError(f"Unsupported layout style: {layout!r}")


def legend_figure(
    labels: Sequence[str],
    colors: Sequence,  # any matplotlib-compatible color spec
    *,
    ncol: int = 5,
    figsize: Tuple[float, float] = (8, 1),
    fontsize: float | str = 10,
):
    """Create a standalone legend figure and return (fig, ax).

    - labels and colors must be same length.
    - No frame; suitable for saving separately.
    """
    import matplotlib.patches as mpatches

    if len(labels) != len(colors):
        raise ValueError("labels and colors must have the same length")

    fig, ax = plt.subplots(figsize=figsize)
    ax.axis("off")

    patches = [mpatches.Patch(color=colors[i], label=str(labels[i])) for i in range(len(labels))]
    _legend_with_labels(
        ax,
        patches,
        [str(label) for label in labels],
        ncol=max(1, min(int(ncol), len(patches))),
        fontsize=fontsize,
        loc="center",
        frameon=False,
    )
    return fig, ax


def _legend_with_labels(target, handles, labels, **options):
    """``target.legend(handles, labels, **options)``, keeping a label that starts with "_".

    matplotlib 3.6.0, the declared floor, drops such a label even with explicit
    handles (3.11 keeps it), so the legend is made with stand-in labels and then
    given the real text, from which it is laid out when drawn
    (TX-AUD-20261019-060). Every legend of the package is made here.
    """
    labels = [str(label) for label in labels]
    stand_ins = ["x" * len(label) if label.startswith("_") else label for label in labels]
    legend = target.legend(handles, stand_ins, **options)
    for text, label in zip(legend.get_texts(), labels, strict=True):
        text.set_text(label)
    return legend


def figure_to_buffer(
    fig, *, fmt: str = "png", dpi: int = 200, bbox_inches: str = "tight"
) -> BytesIO:
    """Serialize a matplotlib figure to a BytesIO buffer and return it.

    Note: The figure is closed to free memory.
    """
    buf = BytesIO()
    fig.savefig(buf, format=fmt, dpi=dpi, bbox_inches=bbox_inches)
    plt.close(fig)
    buf.seek(0)
    return buf


def save_plot(
    fig_or_tuple, filepath: str, *, dpi: int = 300, bbox_inches: str | None = "tight", **kwargs
) -> None:
    """Unified helper to save plots to file.

    Parameters
    ----------
    fig_or_tuple : matplotlib.figure.Figure or tuple
        Either a Figure object or a tuple (fig, ax) as returned by plot functions.
    filepath : str
        Output file path. Format is inferred from extension (e.g., .png, .pdf, .svg).
    dpi : int, default=300
        Resolution for raster formats (PNG, JPG).
    bbox_inches : str or None, default="tight"
        Bounding box style: "tight" crops whitespace; "standard" (or None)
        keeps matplotlib's default margins, i.e. the full figure canvas.
    **kwargs
        Additional keyword arguments passed to fig.savefig().

    Examples
    --------
    >>> import numpy as np
    >>> from traminex import seqdef
    >>> from traminex.visualization import sequence_index_plot
    >>> sequences = seqdef(np.array([list("AABB"), list("ABBB"), list("BBAA")], dtype=object))
    >>> fig, ax = sequence_index_plot(sequences)
    >>> save_plot(fig, "index_plot.png", dpi=300)
    >>> # or
    >>> save_plot((fig, ax), "index_plot.pdf")
    """
    from pathlib import Path

    # Extract figure from tuple if needed
    if isinstance(fig_or_tuple, tuple):
        fig = fig_or_tuple[0]
    else:
        fig = fig_or_tuple

    # Validate filepath
    filepath = str(filepath)
    if not filepath:
        raise ValueError("filepath cannot be empty")

    # Create parent directory if needed
    path = Path(filepath)
    path.parent.mkdir(parents=True, exist_ok=True)

    import logging

    # "standard" is this helper's documented name for matplotlib's default
    # (None); savefig itself rejects it (2026-09-22 audit S6-19).
    if bbox_inches == "standard":
        bbox_inches = None

    # Save figure
    fig.savefig(filepath, dpi=dpi, bbox_inches=bbox_inches, **kwargs)
    logging.getLogger(__name__).info("Plot saved to: %s", filepath)
