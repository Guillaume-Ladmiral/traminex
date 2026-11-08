"""Position-wise discrepancy plot over the R-golden'd seqdiff series (BUILD-from-R-design).

seqdiffplot mirrors TraMineR's plot.seqdiff: a line plot of a per-position
between-group statistic (default "Pseudo R2") or the per-group discrepancy lines
over the sequence positions (SOURCES/TraMineR/R/seqdiff.R). The per-position
statistics / discrepancies come from ``seqdiff`` (R-golden'd in ``seqdiff_*``).
Original matplotlib; provenance docstring-only.
"""

from __future__ import annotations

from typing import Any, Callable, cast

import matplotlib.pyplot as plt

# NB: import via the package (not ``from ...seqdiff import seqdiff``). The seqdiff
# *module* name collides with the seqdiff *function* export; a direct submodule
# import would bind ``traminex.seqdiff`` to the module and shadow the function.
# Going through ``traminex.__getattr__`` returns the function and self-heals the
# binding (see traminex/__init__.py).
from traminex import seqdiff

from ..utils import _close_new_figures_on_error, _legend_with_labels


@_close_new_figures_on_error
def seqdiffplot(
    x,
    group,
    *,
    method: str = "LCS",
    norm="auto",
    weighted: bool = True,
    stat: str = "Pseudo R2",
    figsize=(8, 4),
    ax=None,
    **kwargs,
):
    """Position-wise discrepancy -- re-expresses TraMineR ``plot.seqdiff``.

    With ``stat`` a single statistic name (default "Pseudo R2"), draws one line of
    that per-position statistic. With ``stat="discrepancy"``, draws one line per
    group level (plus "Total") of the per-position discrepancy. The series come
    from ``seqdiff`` (R-golden'd in ``seqdiff_*``); x = sequence positions.

    The defaults are ``seqdiff``'s (LCS distances with ``norm="auto"``), as R's
    ``plot.seqdiff`` draws the object ``seqdiff`` returns. ``seqdist_args=``
    (passed through ``**kwargs``) replaces ``method`` and ``norm``, as in
    ``seqdiff``.

    Returns
    -------
    fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
    """
    compute_seqdiff = cast(Callable[..., Any], seqdiff)
    out = compute_seqdiff(x, group, method=method, norm=norm, weighted=weighted, **kwargs)
    if stat == "discrepancy":
        df = out["discrepancy"]
        ylabel = "Discrepancy"
    else:
        df = out["stat"][[stat]]
        ylabel = stat
    positions = list(df.index)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    else:
        fig = ax.figure
    lines = []
    for col in df.columns:
        (line,) = ax.plot(positions, df[col].to_numpy(dtype=float), marker="o", label=str(col))
        lines.append(line)
    ax.set_xlabel("Position")
    ax.set_ylabel(ylabel)
    ax.set_title("Position-wise group discrepancy")
    if len(df.columns) > 1:
        # Explicit handles and the helper keep a group whose label starts with "_"
        # (TX-REP-20261017-003, TX-AUD-20261019-060).
        _legend_with_labels(ax, lines, [str(col) for col in df.columns], fontsize="small")
    return fig, ax
