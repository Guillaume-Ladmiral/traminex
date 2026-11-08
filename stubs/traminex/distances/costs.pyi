from collections.abc import Mapping, Sequence
from typing import Any, Literal, TypedDict, overload

import numpy as np
import pandas as pd

from ..sequences import STSList

class SeqCostMatrix:
    values: np.ndarray
    from_states: tuple[str, ...]
    to_states: tuple[str, ...]
    time_labels: tuple[str, ...]
    def __init__(
        self,
        values: np.ndarray,
        from_states: tuple[str, ...],
        to_states: tuple[str, ...],
        time_labels: tuple[str, ...],
    ) -> None: ...
    @property
    def shape(self) -> tuple[int, int, int]: ...
    @property
    def ndim(self) -> int: ...
    def __array__(self, dtype=..., copy=...) -> np.ndarray: ...
    def __getitem__(self, key): ...

class _SeqCostResult(TypedDict):
    indel: float | pd.Series | pd.DataFrame
    sm: pd.DataFrame | SeqCostMatrix

class _SeqCostProximityResult(TypedDict):
    indel: float | pd.Series | pd.DataFrame
    prox: pd.DataFrame | SeqCostMatrix

@overload
def seqcost(
    seqdata: STSList,
    method: str,
    cval: float | None = None,
    *,
    with_missing: bool = False,
    miss_cost: float | Sequence[float] | None = None,
    time_varying: bool = False,
    weighted: bool = True,
    transition: str = "both",
    lag: int = 1,
    miss_cost_fixed: bool | None = None,
    state_features: pd.DataFrame | np.ndarray | None = None,
    feature_weights: float | Sequence[float] | None = None,
    feature_type: Mapping[str, Any] | None = None,
    proximities: Literal[False] = False,
    check_costs: bool = True,
) -> _SeqCostResult: ...
@overload
def seqcost(
    seqdata: STSList,
    method: str,
    cval: float | None = None,
    *,
    with_missing: bool = False,
    miss_cost: float | Sequence[float] | None = None,
    time_varying: bool = False,
    weighted: bool = True,
    transition: str = "both",
    lag: int = 1,
    miss_cost_fixed: bool | None = None,
    state_features: pd.DataFrame | np.ndarray | None = None,
    feature_weights: float | Sequence[float] | None = None,
    feature_type: Mapping[str, Any] | None = None,
    proximities: Literal[True],
    check_costs: bool = True,
) -> _SeqCostProximityResult: ...
@overload
def seqcost(
    seqdata: STSList,
    method: str,
    cval: float | None = None,
    *,
    with_missing: bool = False,
    miss_cost: float | Sequence[float] | None = None,
    time_varying: bool = False,
    weighted: bool = True,
    transition: str = "both",
    lag: int = 1,
    miss_cost_fixed: bool | None = None,
    state_features: pd.DataFrame | np.ndarray | None = None,
    feature_weights: float | Sequence[float] | None = None,
    feature_type: Mapping[str, Any] | None = None,
    proximities: bool = False,
    check_costs: bool = True,
) -> _SeqCostResult | _SeqCostProximityResult: ...
def trate_substitution_matrix(
    seqmat: np.ndarray,
    *,
    smoothing: float = 0.0,
    mode: str = "symmetric",
    normalize: str = "none",
    weights: np.ndarray | None = None,
) -> np.ndarray: ...
def constant_substitution_matrix(seqmat: np.ndarray, cval: float = 2.0) -> np.ndarray: ...
def indelslog_costs(
    seqmat: np.ndarray,
    *,
    include_missing: bool = False,
    missing_codes: tuple[int, int] | None = None,
) -> tuple[float, np.ndarray]: ...
def indels_costs(
    seqmat: np.ndarray,
    *,
    include_missing: bool = False,
    missing_codes: tuple[int, int] | None = None,
) -> tuple[float, np.ndarray]: ...
def future_costs(
    seqmat: np.ndarray, *, weights: np.ndarray | None = None
) -> tuple[float, np.ndarray]: ...
def gower_distance(features, *, weights=None, types=None, ordinal_levels=None) -> np.ndarray: ...
def trate_substitution_matrix_timevarying(
    seqmat: np.ndarray,
    *,
    smoothing: float = 0.0,
    mode: str = "symmetric",
    normalize: str = "none",
    weights: np.ndarray | None = None,
    lag: int = 1,
) -> np.ndarray: ...
def indel_from_subm(subm: np.ndarray, *, strategy: str = "mean_offdiag_half") -> float: ...
def local_indel_by_pos(
    seqmat: np.ndarray,
    *,
    method: str = "change",
    base: float = 1.0,
    scale: float = 1.0,
    smooth_window: int | None = None,
) -> np.ndarray: ...
