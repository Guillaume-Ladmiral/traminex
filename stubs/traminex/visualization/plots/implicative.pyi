from typing import Literal, Sequence

from ...stats_p2 import SeqImplic

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
): ...
