from ...stats import seqmeant as seqmeant
from ...stats import seqstatd as seqstatd
from ...stats_p1 import seqmodst as seqmodst
from ..palette import ensure_cmap as ensure_cmap
from ..ticks import apply_time_axis as apply_time_axis
from .distribution import state_distribution_plot as state_distribution_plot

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
): ...
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
): ...
def seqHtplot(x, *, group=None, include_missing: bool = False, figsize=(8, 3), ax=None): ...
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
): ...
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
): ...
def seqdHplot(
    x,
    *,
    group=None,
    main="auto",
    include_missing: bool = False,
    weights=None,
    figsize=(8, 3),
    ax=None,
): ...
