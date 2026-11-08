from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd

from .sequences import STSList as STSList

def group_p(group: Sequence[object], weights: Sequence[float] | None = None) -> pd.Categorical: ...
def sorti(x: STSList, start: str = "end", sort_index: bool = True) -> np.ndarray: ...
def sortv(x: STSList, start: str = "end") -> np.ndarray: ...
def rowmode(v: Any, except_: Any = None) -> str | None | np.ndarray: ...
def seqtabstocc(
    x: STSList,
    include_missing: bool = False,
    sep: str = "-",
    *,
    with_missing: bool | None = ...,
    idxs: int | Sequence[int] = ...,
    weighted: bool = True,
    tlim: int | Sequence[int] | None = ...,
    format: str | None = ...,
) -> pd.DataFrame: ...
def seqindic(
    x: STSList,
    *,
    indic: str | Sequence[str] = ...,
    with_missing: bool = False,
    ipos_args: Mapping[str, Any] | None = None,
    prec_args: Mapping[str, Any] | None = None,
    w: float = 0.5,
) -> pd.DataFrame: ...
def seqindic_dyn(
    x: STSList,
    *,
    indic: str | Sequence[str] = "cplx",
    window_size: float = 0.2,
    sliding: bool = True,
    wstep: int = 1,
    with_missing: bool = False,
    endmiss_as_void: bool = False,
    silent_indic: bool = True,
    **kwargs,
) -> np.ndarray: ...
def seqasnum(x: STSList, *, with_missing: bool = False) -> np.ndarray: ...
def seqnum(x: STSList, *, with_missing: bool = False) -> STSList: ...

class SeqTrateResult:
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

def seqtrate(
    x: STSList,
    *,
    sel_states: Sequence[str] | None = None,
    time_varying: bool = False,
    weighted: bool = True,
    lag: int = 1,
    include_missing: bool = False,
    count: bool = False,
    statl: Sequence[str] | None = None,
    as_dataframe: bool = True,
) -> pd.DataFrame | np.ndarray | SeqTrateResult: ...
def seqhasmiss(x: STSList, per_sequence: bool = False): ...
