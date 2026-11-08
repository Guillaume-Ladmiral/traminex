from collections.abc import Sequence
from typing import Literal

import numpy as np
import pandas as pd

from .sequences import STSList

def dissCompare(
    diss: np.ndarray,
    group: Sequence[object],
    *,
    set: Sequence[object] | None = None,
    s: int = 100,
    seed: int | np.random.Generator | None = 36963,
    stat: Literal["all", "LRT", "BIC"] = "all",
    squared: bool | Literal["LRTonly"] = "LRTonly",
    weighted: bool | Literal["by.group"] = True,
    weights: Sequence[float] | None = None,
    BFopt: Literal[1, 2] | None = None,
    inc_df: bool = False,
) -> pd.DataFrame: ...
def seqCompare(
    seqdata: STSList | Sequence[STSList],
    seqdata2: STSList | Sequence[STSList] | None = None,
    *,
    group: Sequence[object] | None = None,
    set: Sequence[object] | None = None,
    s: int = 100,
    seed: int | np.random.Generator | None = 36963,
    stat: Literal["all", "LRT", "BIC"] = "all",
    squared: bool | Literal["LRTonly"] = "LRTonly",
    weighted: bool | Literal["by.group"] = True,
    opt: Literal[1, 2] | None = None,
    BFopt: Literal[1, 2] | None = None,
    inc_df: bool = False,
    method: str = "LCS",
    **distance_kwargs: object,
) -> pd.DataFrame: ...
def seqLRT(
    seqdata: STSList | Sequence[STSList],
    seqdata2: STSList | Sequence[STSList] | None = None,
    *,
    group: Sequence[object] | None = None,
    set: Sequence[object] | None = None,
    s: int = 100,
    seed: int | np.random.Generator | None = 36963,
    squared: bool | Literal["LRTonly"] = "LRTonly",
    weighted: bool | Literal["by.group"] = True,
    opt: Literal[1, 2] | None = None,
    BFopt: Literal[1, 2] | None = None,
    inc_df: bool = False,
    method: str = "LCS",
    **distance_kwargs: object,
) -> pd.DataFrame: ...
def seqBIC(
    seqdata: STSList | Sequence[STSList],
    seqdata2: STSList | Sequence[STSList] | None = None,
    *,
    group: Sequence[object] | None = None,
    set: Sequence[object] | None = None,
    s: int = 100,
    seed: int | np.random.Generator | None = 36963,
    squared: bool | Literal["LRTonly"] = "LRTonly",
    weighted: bool | Literal["by.group"] = True,
    opt: Literal[1, 2] | None = None,
    BFopt: Literal[1, 2] | None = None,
    inc_df: bool = False,
    method: str = "LCS",
    **distance_kwargs: object,
) -> pd.DataFrame: ...
