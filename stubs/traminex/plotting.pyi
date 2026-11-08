from collections.abc import Mapping
from typing import Any, Literal, Optional, Sequence

import numpy as np

from .samm import SAMMResult
from .sequences import STSList

def save_plot(
    fig_or_tuple,
    filepath: str,
    *,
    dpi: int = 300,
    bbox_inches: str | None = "tight",
    **kwargs,
) -> None: ...
def sequence_index_plot(
    x: STSList,
    *,
    cmap: Optional[Sequence[str]] = None,
    figsize=(8, 4),
    highlight: Optional[Sequence[int]] = None,
    legend: bool = False,
    ax=None,
    sortv=None,
    sampling: int | None = None,
    sample_meth: Literal["systematic", "random"] = "systematic",
    weighted: bool = True,
    random_state: int | np.random.Generator | None = None,
    missing_alpha: float = 1.0,
    void_alpha: float = 0.2,
    apply_hatching: bool = False,
): ...
def plot_SAMM(x: SAMMResult, type: str = "d", **kwargs: object): ...
def plot_seqimplic(
    x,
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
def sequence_index_plot_grouped(
    x: STSList,
    labels: np.ndarray,
    *,
    order: Optional[Sequence[int]] = None,
    figsize=(8, 4),
    cmap: Optional[Sequence[str]] = None,
    show_separators: bool = True,
    legend: bool = True,
): ...
def overlay_representatives_on_index_plot(
    x: STSList,
    labels: np.ndarray,
    D: np.ndarray | None = None,
    *,
    figsize=(8, 4),
): ...
def state_distribution_plot(
    x: STSList,
    *,
    include_missing: bool = False,
    weights: np.ndarray | None = None,
    figsize=(8, 3),
    ax=None,
    mode: str = "d",
    cmap=None,
    missing_alpha: float = 1.0,
    void_alpha: float = 0.4,
): ...
def chronogram(
    x: STSList, *, include_missing: bool = False, figsize=(8, 3), seaborn: bool = False, ax=None
): ...
def transition_heatmap(
    x: STSList,
    *,
    normalize: bool = True,
    weighted: bool = True,
    figsize=(5, 4),
    seaborn: bool = False,
    ax=None,
): ...
def cluster_chronograms(
    x: STSList,
    labels: np.ndarray,
    *,
    include_missing: bool = False,
    cols: int = 3,
    figsize=(8, 3),
    shared_scale: bool = True,
    layout: str = "column",
    nrows: int | None = None,
    ncols: int | None = None,
): ...
def cluster_chronograms_with_representatives(
    x: STSList,
    labels: np.ndarray,
    D: np.ndarray,
    *,
    include_missing: bool = False,
    cols: int = 3,
    figsize=(8, 3),
    shared_scale: bool = True,
    marker_size: int = 3,
    layout: str = "column",
    nrows: int | None = None,
    ncols: int | None = None,
): ...
def seqfplot(
    x,
    top: int | None = 10,
    *,
    group=None,
    main="auto",
    cpal=None,
    include_missing: bool = True,
    border: bool = True,
    figsize=(8, 4),
    ax=None,
): ...
def seqfrequency_plot(
    x, top: int | None = 10, *, seaborn: bool = True, figsize=(8, 4), ax=None
): ...
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
def seqrplot(
    x,
    *,
    diss=None,
    criterion: str = "density",
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
): ...
def seqrfplot(
    x,
    diss,
    *,
    k=None,
    sortv: str = "mds",
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
): ...
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
): ...
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
def seqedplot(
    eseq,
    *,
    type: str = "survival",
    group: Optional[Sequence] = None,
    ignore: Optional[Sequence[str]] = None,
    breaks: int = 20,
    ages: Optional[tuple[float, float]] = None,
    cpal=None,
    figsize=(8, 4),
    axes=None,
): ...
def seqeordplot(*args, **kwargs): ...
def seqalignplot(
    x,
    indices,
    *,
    indel=1.0,
    sm=None,
    include_missing: bool = False,
    cpal=None,
    figsize=(8, 3),
    ax=None,
): ...
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
): ...
def seqplotMD(channels, *, type: str = "i", cpal=None, dnames=None, figsize=(10, 4), **kwargs): ...
def plot_subseqelist(x, *, freq: str = "Support", top=None, figsize=(8, 4), ax=None): ...
def plot_subseqelistchisq(x, *, top=None, figsize=(10, 4)): ...
def seqpcplot(x, *, filt=None, cpal=None, figsize=(8, 4), ax=None): ...
def plot_eseq(x, **kwargs): ...
def plot_seqelist(x, **kwargs): ...
def seqplot(seqdata, *, type=..., group=None, main="auto", **kwargs): ...
def duration_boxplot(
    x,
    labels: np.ndarray | None = None,
    *,
    by_state: bool = False,
    seaborn: bool = True,
    figsize=(8, 4),
): ...
def state_legend(x, ax=None, *, ncol: int = 4): ...
def state_legend_figure(x, *, ncol: int = 5, figsize=(8, 1), fontsize: float | str = 10): ...
