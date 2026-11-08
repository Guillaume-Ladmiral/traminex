from typing import Mapping, Sequence

import numpy as np
import pandas as pd

from .sequences import STSList

def seqdomassoc(
    seqdata_dom: Sequence[STSList],
    *,
    assoc: str | Sequence[str] = ...,
    rep_method: str = "overall",
    diss_dom: Sequence[np.ndarray] | None = None,
    wrange: Sequence[float] | None = None,
    cross_table: bool = False,
    seqrep_args: Mapping[str, object] | None = None,
    seqrf_args: Mapping[str, object] | None = None,
    with_missing: bool = False,
    weighted: bool = True,
    struct_zero: bool = True,
    p_value: bool = True,
    dnames: Sequence[str] | None = None,
) -> pd.DataFrame: ...
