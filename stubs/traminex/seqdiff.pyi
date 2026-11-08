from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd

from .sequences import STSList

def seqdiff(
    x: STSList,
    group: np.ndarray | pd.Series | STSList,
    *,
    cmprange: Sequence[int] = ...,
    method: str = "LCS",
    norm: str | None = "auto",
    seqdist_args: Mapping[str, Any] | None = None,
    with_missing: bool = False,
    weighted: bool = True,
    squared: bool = False,
    seqdist_arg: Mapping[str, Any] | None = None,
) -> dict[str, pd.DataFrame]: ...
