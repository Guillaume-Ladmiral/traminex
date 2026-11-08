from typing import Literal, Sequence

from ...sequences import STSList
from ...stats_p2 import seqsurv as seqsurv
from ..palette import ensure_cmap as ensure_cmap

def seqsurvplot(
    x,
    *,
    groups=None,
    with_missing: bool = False,
    cpal=None,
    figsize=(8, 4),
    ax=None,
): ...
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
    title: str | Sequence[str] | None = ...,
    cex_plot: float = ...,
    withlegend: Literal["auto", "right"] | bool = ...,
    axes: Literal["all", "bottom"] | bool = ...,
    **style: object,
): ...
