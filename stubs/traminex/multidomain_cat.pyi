from typing import Any, Sequence

import numpy as np

from .sequences import STSList as STSList

def seqMD(
    channels: Sequence[STSList],
    method: str | None = None,
    norm: str = "none",
    indel: Any = "auto",
    sm: Any = None,
    with_missing: Any = None,
    full_matrix: bool = True,
    link: str = "sum",
    cval: float = 2.0,
    miss_cost: float = 2.0,
    cweight: Sequence[float] | None = None,
    what: str = "MDseq",
    ch_sep: str = "+",
    fill_with_miss: bool = True,
    check_costs: bool = True,
) -> STSList | dict | np.ndarray: ...
def seqdistmc(
    channels: Sequence[STSList],
    what: str = "diss",
    ch_sep: str = "@@@@TraMineRSep@@@@",
    **kwargs: Any,
) -> STSList | dict | np.ndarray: ...
