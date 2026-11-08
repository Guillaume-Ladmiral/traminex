from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd

from ..sequences import STSList
from .costs import (
    SeqCostMatrix,
)
from .costs import (
    indel_from_subm as indel_from_subm,
)
from .costs import (
    local_indel_by_pos as local_indel_by_pos,
)
from .costs import (
    seqcost as seqcost,
)
from .costs import (
    trate_substitution_matrix as trate_substitution_matrix,
)
from .costs import (
    trate_substitution_matrix_timevarying as trate_substitution_matrix_timevarying,
)
from .omloc import omloc as omloc

__all__: tuple[str, ...]

def om_trate_pairwise(
    seqmat: np.ndarray,
    *,
    indel: float | None = None,
    trate_smoothing: float = 0.0,
    trate_mode: str = "symmetric",
    trate_normalize: str = "none",
    trate_weights: np.ndarray | None = None,
    norm: str | None = "length",
    deduplicate: bool = False,
    condensed: bool = False,
) -> np.ndarray: ...
def pairwise(seqmat: np.ndarray, method: str = "OM", **kwargs) -> np.ndarray: ...
def pairwise_between_sets(
    A: np.ndarray | Sequence[np.ndarray],
    B: np.ndarray | Sequence[np.ndarray],
    *,
    method: str = "OM",
    deduplicate: bool = False,
    **kwargs,
) -> np.ndarray: ...
def seqdist(x, method: str = "OM", **kwargs) -> np.ndarray: ...
def to_condensed(D: np.ndarray) -> np.ndarray: ...
def from_condensed(d: np.ndarray, n: int | None = None) -> np.ndarray: ...
def sanitize_distance_matrix(D: np.ndarray) -> np.ndarray: ...
def get_lcp_length(a: np.ndarray, b: np.ndarray) -> int: ...
def get_rlcp_length(a: np.ndarray, b: np.ndarray) -> int: ...
def seqsubm(
    seqmat: STSList | np.ndarray,
    method: str | None = None,
    cval: float | None = None,
    **kwargs: Any,
) -> np.ndarray | pd.DataFrame | SeqCostMatrix | None: ...
def seqLLCS(seq1, seq2) -> int: ...
def seqLLCP(seq1, seq2) -> int: ...
def seqRLCP(x, *, norm: str | None = None) -> np.ndarray: ...
