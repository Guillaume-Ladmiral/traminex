from typing import Any, Mapping, Sequence

import pandas as pd

from .sequences import STSList as STSList

def seqformat(
    data: Any,
    *,
    from_: str | None = None,
    to: str | None = None,
    id_col: str | None = None,
    time_col: str | None = None,
    state_col: str | None = None,
    var: Sequence[int | str] | int | str | None = None,
    compress: bool = False,
    nrep: int | None = None,
    tevent: Any = None,
    stsep: str | None = None,
    covar: Sequence[int | str] | None = None,
    sps_in: Mapping[str, str] | None = None,
    sps_out: Mapping[str, str] | None = None,
    id: int | str | Sequence[object] | None = 1,
    begin: int | str = 2,
    end: int | str = 3,
    status: int | str = 4,
    process: bool = True,
    pdata: pd.DataFrame | str | None = None,
    pvar: Sequence[int | str] | None = None,
    limit: int | float = 100,
    overwrite: bool = True,
    fillblanks: str | None = None,
    tmin: int | None = None,
    tmax: int | None = None,
    missing: str | Sequence[str] = "*",
    with_missing: bool = True,
    right: str | None = "DEL",
    compressed: bool | None = None,
    nr: str | None = None,
) -> pd.DataFrame | STSList | list[str]: ...
def seqdef_from_long(
    data: pd.DataFrame,
    *,
    id_col: str = "id",
    time_col: str = "time",
    state_col: str = "state",
    time_grid: Sequence[object] | None = None,
    **kwargs,
) -> STSList: ...
def seqstart(
    x: STSList,
    data_start: int | Sequence[int],
    new_start: int | Sequence[int],
    *,
    tmin: int | None = None,
    tmax: int | None = None,
    missing: str | None = None,
) -> STSList: ...
def toPersonPeriod(x: STSList, *, time_col: str = "timestamp") -> pd.DataFrame: ...
def STS_to_TSE(
    x: STSList, *, state_changes_only: bool = True, time_col: str = "time"
) -> pd.DataFrame: ...
def to_SPS(
    x: STSList, *, state_sep: str = "-", sdsep: str = "/", compress: bool = True
) -> list[str]: ...
def to_DSS(x: STSList, *, state_sep: str = "-") -> list[str]: ...
