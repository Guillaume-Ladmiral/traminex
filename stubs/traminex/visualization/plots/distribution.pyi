import numpy as np

from ...sequences import STSList as STSList
from ...stats import seqstatd as seqstatd
from ..ticks import apply_time_axis as apply_time_axis

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
