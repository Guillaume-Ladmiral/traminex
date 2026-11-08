from typing import Any, Sequence

import numpy as np
import pandas as pd

from .seqformat import seqdef_from_long as seqdef_from_long
from .sequences import STSList as STSList

def seqdef(
    data: STSList | pd.DataFrame | np.ndarray | Sequence[Sequence[str | int | float | None]],
    var: Sequence[int | str] | int | str | None = None,
    *,
    informat: str = "STS",
    stsep: str | None = None,
    alphabet: str | Sequence[str] | None = None,
    states: str | Sequence[str] | None = None,
    id: str | Sequence[object] | None = None,
    weights: Sequence[float] | None = None,
    start: int | None = 1,
    left: object = None,
    right: object = "DEL",
    gaps: object = None,
    missing: object = None,
    void: str = "%",
    nr: str = "*",
    cnames: str | Sequence[str] | None = None,
    xtstep: int = 1,
    tick_last: bool = False,
    cpal: str | Sequence[str] | None = None,
    missing_color: str = "darkgrey",
    labels: str | Sequence[str] | None = None,
    **seqformat_kwargs: Any,
) -> STSList: ...
def FCE_to_TSE(
    seqdata: pd.DataFrame,
    *,
    id: int | str | Sequence[object] | None = None,
    cols: Sequence[int | str] | int | str = (),
    eventlist: Sequence[object] | None = None,
    firstEvent: object | None = None,
    time_col: str = "time",
) -> pd.DataFrame: ...
def TSE_to_STS(
    df: pd.DataFrame,
    *,
    id_col: int | str = "id",
    time_col: int | str = "time",
    state_col: int | str | None = "state",
    event_col: int | str = "event",
    stm: pd.DataFrame | None = None,
    tmin: int = 1,
    tmax: int | None = None,
    first_state: str = "None",
    **kwargs: Any,
) -> STSList: ...
def STS_to_SPELL(
    x: STSList,
    *,
    with_missing: bool = True,
    right: object = None,
    pdata: pd.DataFrame | None = None,
    pvar: Sequence[int | str] | None = None,
) -> pd.DataFrame: ...
def STS_to_SRS(x: STSList, *, nrep: int = 0) -> pd.DataFrame: ...
def HSPELL_to_STS(
    df: pd.DataFrame,
    *,
    id_col: str | None = None,
    start_col: str | None = None,
    end_col: str | None = None,
    state_col: str | None = None,
    begin: Sequence[str | int] | None = None,
    end: Sequence[str | int] | None = None,
    status: Sequence[str | int] | None = None,
    fixed_status: str | Sequence[str] | None = None,
    pvar: str | None = None,
    overwrite: bool = True,
    overlap_priority: str = "later",
    fillblanks: str | None = None,
    tmin: int | None = None,
    tmax: int | None = None,
    endObs: str | float | Sequence[float] | None = None,
    time_grid: Sequence[int] | None = None,
    fill: str | None = None,
    interval_closed: str = "both",
    tmin_col: str | None = None,
    tmax_col: str | None = None,
    **kwargs,
) -> STSList: ...
