import math
from dataclasses import dataclass
from typing import Any, Literal, Mapping, Sequence, overload

import numpy as np
import pandas as pd

from .sequences import STSList as STSList

def seqient(
    x: STSList,
    *,
    norm: bool = True,
    base: float = math.e,
    include_missing: bool = False,
    silent: bool = True,
) -> np.ndarray: ...
def seqici(x: STSList, *, include_missing: bool = False, silent: bool = True) -> np.ndarray: ...
def seqgen(
    n: int,
    length: int,
    alphabet: Sequence[Any] | str | int | float,
    *,
    p: Sequence[float] | None = None,
    seed: int | None = None,
) -> np.ndarray: ...
def seqrecode(
    x: STSList,
    recodes: dict[str, list[str]],
    *,
    otherwise: str | None = None,
    labels: Mapping[str, str] | Sequence[str] | None = None,
    cpal: Mapping[str, str] | Sequence[str] | None = None,
) -> STSList: ...
def seqconc(
    data,
    var: Sequence[int | str] | int | str | None = None,
    sep: str = "-",
    vname: str = "Sequence",
    void=...,
) -> np.ndarray: ...
@overload
def seqsep(seqdata: list[str], sl: int = 1, sep: str = "-") -> list[str]: ...
@overload
def seqsep(seqdata: np.ndarray, sl: int = 1, sep: str = "-") -> np.ndarray: ...
@overload
def seqsep(seqdata: pd.Series, sl: int = 1, sep: str = "-") -> pd.Series: ...
def seqdecomp(
    data,
    var: Sequence[int | str] | int | str | None = None,
    sep: str = "-",
    miss: str | Sequence[str] = "NA",
    vnames: Sequence[str] | str | None = None,
) -> np.ndarray: ...
def seqdim(data) -> np.ndarray: ...
def seqmpos(seq1, seq2, *, with_missing: bool = False) -> int: ...
def seqcomp(x, y) -> bool: ...
def seqsubsn(x: STSList, *, dss: bool = True, with_missing: bool = False) -> np.ndarray: ...
def seqpm(x: STSList, pattern: str, *, sep: str = ""): ...
@dataclass
class SeqImplic:
    indices: np.ndarray
    levels: list[str]
    alphabet: list[str]
    positions: list
    labels: list[str] | None = None
    type: str = "indice"
    cpal: list[str] | None = None
    xtstep: int = 1
    tick_last: bool = False

def seqimplic(
    x: STSList,
    group: Sequence,
    *,
    with_missing: bool = False,
    weighted: bool = True,
    na_rm: bool = True,
) -> SeqImplic: ...
def seqsha(
    x: STSList,
    time: Sequence[int],
    event: Sequence[bool],
    *,
    include_present: bool = False,
    align_end: bool = False,
    covar=None,
): ...
def seqgranularity(x: STSList, *, tspan: int = 3, method: str = "last") -> STSList: ...
@dataclass
class SeqPcFilter:
    type: str
    value: str
    level: float | None = None
    def apply(self, x) -> np.ndarray: ...

def seqpcfilter(method: str = "minfreq", level: float | None = 0.05) -> SeqPcFilter: ...
def seqauto(x: STSList, *, order: int = 1, measure: str = "cv") -> np.ndarray: ...
def seqcta(
    x: STSList,
    *,
    subseq: int = 5,
    time=None,
    event=None,
    initial_state: str | None = None,
    covar=None,
): ...

class SeqSurv:
    curves: dict[str, np.ndarray]
    states: list[str]
    curve_axis: str = "state"
    strata: list[str]
    groups: list[str]
    ltext: list[str]
    cpal: list[str] | None
    xtstep: int
    tick_last: bool
    max_time: int
    follow_up: dict[str, float]
    def __init__(
        self,
        curves: dict[str, np.ndarray],
        states: list[str],
        curve_axis: str = "state",
        strata: list[str] = ...,
        groups: list[str] = ...,
        ltext: list[str] = ...,
        cpal: list[str] | None = None,
        xtstep: int = 1,
        tick_last: bool = False,
        max_time: int = 0,
        follow_up: dict[str, float] = ...,
    ) -> None: ...
    def summary(self): ...

def seqsurv(
    x: STSList,
    *,
    groups: Sequence[object] | None = None,
    per_state: bool = False,
    state: Sequence[str] | str | None = None,
    with_missing: bool = False,
) -> SeqSurv: ...
def seqpolyads(
    seqlist: Sequence[STSList],
    *,
    a: int = 1,
    method: str = "HAM",
    pair_weights: Sequence[float] | None = None,
    s: int = 36963,
    T: int = 1000,
    core: int | Literal["auto"] = 1,
    replace: bool = True,
    weighted: bool = True,
    with_missing: bool = False,
    rand_weight_type: int = 1,
    role_weights: Sequence[float] | None = None,
    show_time: bool = False,
    **distance_kwargs: object,
) -> dict[str, np.ndarray]: ...
