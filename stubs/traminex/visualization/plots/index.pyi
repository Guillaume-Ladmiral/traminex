from typing import Literal, Sequence

import numpy as np

from ...sequences import STSList as STSList
from ..palette import ensure_cmap as ensure_cmap
from ..series import select_index_rows as select_index_rows
from ..ticks import apply_time_axis as apply_time_axis

def sequence_index_plot(
    x: STSList,
    *,
    cmap: Sequence[str] | None = None,
    figsize=(8, 4),
    highlight: Sequence[int] | None = None,
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
def sequence_index_plot_grouped(
    x: STSList,
    labels: np.ndarray,
    *,
    order: Sequence[int] | None = None,
    figsize=(8, 4),
    cmap: Sequence[str] | None = None,
    show_separators: bool = True,
    legend: bool = True,
): ...
def overlay_representatives_on_index_plot(
    x: STSList, labels: np.ndarray, D: np.ndarray | None = None, *, figsize=(8, 4)
): ...
