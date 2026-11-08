from typing import Sequence

import numpy as np
import pandas as pd

from .sequences import STSList as STSList

def seqmodst(x: STSList, *, weighted: bool = True, include_missing: bool = False) -> STSList: ...
def seqmeant(
    x: STSList,
    *,
    weighted: bool = True,
    include_missing: bool = False,
    prop: bool = False,
    serr: bool = False,
) -> pd.DataFrame: ...
def seqstatf(
    x: STSList, *, weighted: bool = True, include_missing: bool = False
) -> pd.DataFrame: ...
def seqstatl(
    data,
    var: Sequence[int | str] | int | str | None = None,
    format: str = "STS",
) -> list[str] | list[float]: ...
def seqfind(ref: STSList, target: STSList) -> np.ndarray: ...
def seqfpos(
    x: STSList,
    state: str | None = None,
    *,
    states: Sequence[str] | None = None,
    include_missing: bool = False,
) -> np.ndarray | pd.DataFrame: ...
def seqfposend(
    x: STSList,
    state: str | None = None,
    *,
    states: Sequence[str] | None = None,
    include_missing: bool = False,
    lead: int = 0,
    from_seq_start: bool = True,
) -> np.ndarray | pd.DataFrame: ...
