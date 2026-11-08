from collections.abc import Mapping
from typing import Any, Literal

import numpy as np

from ...cluster import seqrep as seqrep
from ...repseq_rf import seqrf as seqrf
from ..palette import ensure_cmap as ensure_cmap

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
