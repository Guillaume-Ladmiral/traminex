from typing import Literal

import numpy as np
from _typeshed import Incomplete

from .sequences import STSList as STSList

class SeqAlign:
    operation: Incomplete
    seq1: Incomplete
    seq2: Incomplete
    cost: Incomplete
    total_cost: Incomplete
    cost_matrix: Incomplete
    operation_matrix: Incomplete
    def __init__(
        self,
        operation: np.ndarray,
        seq1: np.ndarray,
        seq2: np.ndarray,
        cost: np.ndarray,
        cost_matrix: np.ndarray,
        operation_matrix: np.ndarray,
    ) -> None: ...

def seqalign(
    x: STSList,
    indices: tuple[int, int],
    *,
    indel: float = 1.0,
    sm: np.ndarray | None = None,
    include_missing: bool = False,
) -> SeqAlign: ...
def dissrf(
    diss: np.ndarray,
    *,
    k: int | None = None,
    sortv=...,
    weights: np.ndarray | None = None,
    grp_meth: Literal["prop", "first", "random"] = "prop",
    squared: bool = False,
    pow: int | None = None,
    random_state: int | np.random.Generator | None = None,
) -> dict: ...
