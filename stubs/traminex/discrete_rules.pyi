from collections.abc import Mapping, Sequence
from typing import Literal

import pandas as pd

from .events import SubSeqElist

def createdatadiscrete(
    ids: Sequence[object],
    data: pd.DataFrame,
    vars: Sequence[str],
    agemin: pd.DataFrame,
    agemax: pd.DataFrame,
    supvar: Sequence[str] | None = None,
) -> dict[str, object]: ...
def seqerulesdisc(
    fsubseq: SubSeqElist,
    datadiscr: Mapping[str, object],
    tsef: pd.DataFrame,
    *,
    pvalue: float = 0.1,
    supvars: Sequence[str] | None = None,
    adjust: bool = True,
    topt: bool = False,
    link: Literal["cloglog", "logit"] = "cloglog",
    dep: str | None = None,
) -> pd.DataFrame: ...
