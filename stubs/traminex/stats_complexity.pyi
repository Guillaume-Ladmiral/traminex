from collections.abc import Sequence

import numpy as np
from numpy.lib.mixins import NDArrayOperatorsMixin

from .sequences import STSList as STSList
from .stats_basic import seqistatd as seqistatd
from .stats_p2 import seqpm as seqpm

class SeqIVarDurResult(NDArrayOperatorsMixin):
    vmax: np.ndarray
    meand: np.ndarray
    def __init__(self, values: np.ndarray, *, vmax: np.ndarray, meand: np.ndarray) -> None: ...
    @property
    def shape(self) -> tuple[int, ...]: ...
    @property
    def dtype(self) -> np.dtype: ...
    @property
    def ndim(self) -> int: ...
    @property
    def size(self) -> int: ...
    def __array__(self, dtype=..., copy=...) -> np.ndarray: ...
    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs): ...
    def __getitem__(self, key): ...
    def __len__(self) -> int: ...
    def __iter__(self): ...
    def copy(self) -> SeqIVarDurResult: ...

def seqtransn_count(x: STSList, *, normalized: bool = True) -> np.ndarray: ...
def turbulence_proxy(x: STSList, *, normalized: bool = True) -> float: ...
def seqST(
    x: STSList,
    *,
    norm: bool = False,
    include_missing: bool = False,
    type: int = 1,
    silent: bool = True,
) -> np.ndarray: ...
def seqdss(x: STSList, *, include_missing: bool = False) -> STSList: ...
def seqdur(x: STSList, *, include_missing: bool = False) -> np.ndarray: ...
def seqivolatility(
    x: STSList, *, w: float = 0.5, include_missing: bool = False, adjust: bool = True
) -> np.ndarray: ...
def turbulence_exact(
    x: STSList, *, norm: bool = True, include_missing: bool = False, type: int = 1
) -> float: ...
def seqivardur(x: STSList, *, type: int = 1, include_missing: bool = False) -> SeqIVarDurResult: ...
def seqlogp(
    x: STSList,
    *,
    prob: str | np.ndarray = "trate",
    time_varying: bool = True,
    begin: str | np.ndarray = "freq",
    weighted: bool = True,
    include_missing: bool = False,
) -> np.ndarray: ...
def seqmaintokens(
    x: STSList,
    *,
    k: int = 8,
    mint: float | None = None,
    weighted: bool = True,
    include_missing: bool = False,
) -> np.ndarray: ...
def seqprecarity(
    x: STSList,
    *,
    correction: float | Sequence[float] | np.ndarray | None = None,
    otto: float = 0.2,
    a: float = 1.0,
    b: float = 1.2,
    method: str = "TRATEDSS",
    state_order: list[str] | None = None,
    state_equiv: list[list[str]] | None = None,
    stprec: list[float] | None = None,
    include_missing: bool = False,
    penalized: str | bool = "BOTH",
    weight_type: str = "ADD",
    pow: float = 1.0,
    border_effect: float = 10.0,
) -> np.ndarray: ...
def seqidegrad(
    x: STSList,
    *,
    state_order: list[str] | None = None,
    state_equiv: list[list[str]] | None = None,
    stprec: list[float] | None = None,
    with_missing: bool = False,
    penalized: str | bool = "BOTH",
    method: str = "RANK",
    weight_type: str = "ADD",
    pow: float | bool = 1.0,
    spell_integr: bool | None = None,
    border_effect: float = 10.0,
) -> np.ndarray: ...
def seqinsecurity(
    x: STSList,
    *,
    correction: float | Sequence[float] | np.ndarray | None = None,
    state_order: list[str] | None = None,
    state_equiv: list[list[str]] | None = None,
    stprec: list[float] | None = None,
    with_missing: bool = False,
    pow: float = 1.0,
    spow: float | None = None,
    bound: bool = False,
    method: str = "RANK",
    penalized: str | bool = "BOTH",
    weight_type: str = "ADD",
    border_effect: float = 10.0,
) -> np.ndarray: ...
def seqprecstart(
    x: STSList,
    *,
    state_order: list[str] | None = None,
    state_equiv: list[list[str]] | None = None,
    stprec: list[float] | None = None,
    include_missing: bool = False,
) -> np.ndarray: ...
def seqibad(
    x: STSList,
    *,
    pow: float = 1.0,
    state_order: list[str] | None = None,
    state_equiv: list[list[str]] | None = None,
    stprec: list[float] | None = None,
    include_missing: bool = False,
) -> np.ndarray: ...
def seqintegr(
    x: STSList, *, state: str | None = None, pow: float = 1.0, include_missing: bool = False
) -> np.ndarray: ...
def seqipos(
    x: STSList,
    *,
    pos_states: str | Sequence[str] | None = None,
    neg_states: str | Sequence[str] | None = None,
    index: str = "share",
    pow: float = 1.0,
    w: float = 0.5,
    dss: bool | None = None,
    include_missing: bool = False,
) -> np.ndarray: ...
