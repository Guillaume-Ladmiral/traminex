from dataclasses import dataclass
from typing import Literal

import numpy as np

from ..sequences import STSList

__all__ = ["IndexSelection", "resolve_sort_order", "select_index_rows"]

@dataclass(frozen=True)
class IndexSelection:
    source_rows: np.ndarray
    display_rows: np.ndarray
    weights: np.ndarray
    labels: tuple[object, ...]

def resolve_sort_order(x, sortv=None, *, weighted: bool = True) -> np.ndarray: ...
def select_index_rows(
    x: STSList,
    *,
    sortv: object = None,
    sampling: int | None = None,
    sample_meth: Literal["systematic", "random"] = "systematic",
    weighted: bool = True,
    random_state: int | np.random.Generator | None = None,
) -> IndexSelection: ...
def transversal_entropy(x: STSList) -> np.ndarray: ...
def grouped_transversal_entropy(
    x: STSList, group
) -> tuple[tuple[object, ...], tuple[np.ndarray, ...]]: ...
