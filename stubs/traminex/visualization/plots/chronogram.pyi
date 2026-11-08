import numpy as np

from ...sequences import STSList as STSList
from ...stats import seqstatd as seqstatd
from ..ticks import apply_time_axis as apply_time_axis
from ..utils import determine_layout as determine_layout

def chronogram(
    x: STSList, *, include_missing: bool = False, figsize=(8, 3), seaborn: bool = False, ax=None
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
