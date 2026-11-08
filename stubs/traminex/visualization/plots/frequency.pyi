from ...stats import seqtab as seqtab
from ..palette import ensure_cmap as ensure_cmap

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
