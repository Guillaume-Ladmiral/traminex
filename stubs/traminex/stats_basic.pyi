from collections.abc import Sequence

import numpy as np
import pandas as pd

from .sequences import STSList as STSList

class SeqStatDResult(dict[str, np.ndarray]):
    nbseq: float
    weighted: bool
    norm: bool
    xtlab: list[str]
    xtstep: int
    tick_last: bool
    cpal: list[str] | None
    def __init__(
        self,
        frequencies: np.ndarray,
        valid_states: np.ndarray,
        entropy: np.ndarray,
        *,
        nbseq: float,
        weighted: bool,
        norm: bool,
        xtlab: list[str],
        xtstep: int,
        tick_last: bool,
        cpal: list[str] | None,
    ) -> None: ...

def seqlength(
    x: STSList, with_missing: bool = True, *, as_dataframe: bool = False
) -> np.ndarray | pd.DataFrame: ...
def seqistatd(x: STSList, *, include_missing: bool = False, prop: bool = False) -> np.ndarray: ...
def seqstatd(
    x: STSList,
    *,
    weighted: bool = True,
    include_missing: bool = False,
    norm: bool = True,
    weights: np.ndarray | None = None,
) -> SeqStatDResult: ...
def seqtab(
    x: STSList,
    idxs: int | Sequence[int] = ...,
    weighted: bool = True,
    format: str = "SPS",
    tlim: int | Sequence[int] | None = ...,
    *,
    sep: str = "-",
    weights: np.ndarray | None = None,
    sort: bool = True,
    normalize: bool = True,
    include_missing: bool = True,
) -> pd.DataFrame: ...
def seqtransn(
    x: STSList,
    *,
    with_missing: bool = False,
    norm: bool = False,
    pweight: bool = False,
) -> np.ndarray: ...
def _transition_count_matrix(x: STSList, *, normalize: bool = False) -> np.ndarray: ...
def entropy(x: STSList, *, normalized: bool = True) -> float: ...
