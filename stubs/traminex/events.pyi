from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Self, overload

import numpy as np
import pandas as pd

from .sequences import STSList

def seqetm(
    x: STSList = ...,
    method: str = "transition",
    *,
    use_labels: bool = True,
    sep: str = ">",
    bp: str = "",
    ep: str = "end",
    seq: STSList = ...,
) -> pd.DataFrame: ...
def seqe2stm(
    events: list[str],
    *,
    drop_matrix: pd.DataFrame | None = None,
    drop_list: dict[str, object] | None = None,
    first_state: str = "None",
) -> pd.DataFrame: ...
@dataclass
class Eseq:
    id_code: int
    times: np.ndarray
    codes: np.ndarray
    dictionary: list[str]
    weight: float = ...
    obs_time: float = ...
    def as_character(self) -> str: ...

class SeqElist:
    sequences: list[Eseq]
    dictionary: list[str]
    def __init__(
        self,
        sequences: list[Eseq],
        dictionary: list[str],
        *,
        flat_id: list[str],
        flat_time: np.ndarray,
        flat_label: list[str],
    ) -> None: ...
    @property
    def id(self) -> list[str]: ...
    @property
    def time(self) -> np.ndarray: ...
    @property
    def label(self) -> list[str]: ...
    def __len__(self) -> int: ...
    def __getitem__(self, idx): ...
    def __iter__(self): ...
    @property
    def weights(self) -> np.ndarray: ...
    @weights.setter
    def weights(self, value) -> None: ...
    @property
    def lengths(self) -> np.ndarray: ...
    @lengths.setter
    def lengths(self, value) -> None: ...
    def as_character(self) -> list[str]: ...

def seqeid(eseq=..., *, s=...) -> np.ndarray | int: ...
def seqeage(eseq, event_list) -> np.ndarray: ...
def seqeweight(eseq=..., *, s=...) -> np.ndarray | float: ...
def seqecontain(
    eseq=...,
    event_list=...,
    *,
    unknown_exclude: bool = False,
    seq=...,
    eventList=...,
    exclude: bool = ...,
) -> np.ndarray | bool: ...
def seqecreate(
    data: pd.DataFrame | STSList | None = None,
    *,
    id: Sequence[object] | None = None,
    timestamp: Sequence[float] | None = None,
    event: Sequence[object] | None = None,
    end_event: str | float | None = None,
    tevent: str | np.ndarray = "transition",
    use_labels: bool = True,
    weighted: bool = True,
    id_col: str = "id",
    time_col: str = "timestamp",
    event_col: str = "event",
    endEvent: str | float | None = ...,
) -> SeqElist: ...
def is_eseq(x: object = ..., *, s: object = ...) -> bool: ...
def is_seqelist(x: object = ..., *, s: object = ...) -> bool: ...
def is_stslist(x: object) -> bool: ...
@dataclass
class SeqEConstraint:
    max_gap: float = ...
    window_size: float = ...
    age_min: float = ...
    age_max: float = ...
    age_max_end: float = ...
    count_method: int = ...

def seqeconstraint(
    max_gap: float = -1,
    window_size: float = -1,
    age_min: float = -1,
    age_max: float = -1,
    age_max_end: float = -1,
    count_method: int | str = 1,
    *,
    maxGap: float = ...,
    windowSize: float = ...,
    ageMin: float = ...,
    ageMax: float = ...,
    ageMaxEnd: float = ...,
    countMethod: int | str = ...,
) -> SeqEConstraint: ...

class SubSeqElist:
    eseq: SeqElist
    constraint: SeqEConstraint
    subseq: list[str]
    data: pd.DataFrame
    type: str
    count_total: float | None
    def __init__(
        self,
        eseq,
        constraint,
        subseq_objs,
        data,
        type: str = "frequent",
        count_total: float | None = None,
    ) -> None: ...
    def __len__(self) -> int: ...
    def as_character(self) -> list[str]: ...
    @overload
    def __getitem__(self, idx: int | np.integer) -> str: ...
    @overload
    def __getitem__(self, idx: slice | Sequence[int] | Sequence[bool] | np.ndarray) -> Self: ...

def is_subseqelist(x: object) -> bool: ...
def seqefsub(
    eseq: SeqElist = ...,
    *,
    str_subseq: list[str] | None = None,
    min_support: float | None = None,
    pmin_support: float | None = None,
    constraint: SeqEConstraint | None = None,
    max_k: int = -1,
    weighted: bool = True,
    seq: SeqElist = ...,
    strsubseq: list[str] | None = ...,
    minSupport: float | None = ...,
    pMinSupport: float | None = ...,
    maxK: int = ...,
) -> SubSeqElist: ...
def seqeapplysub(
    subseq: SubSeqElist,
    method: int | str | None = None,
    constraint: SeqEConstraint | None = None,
    rules: bool = False,
) -> pd.DataFrame: ...
def seqerules(
    fsubseq: SubSeqElist, sortv: str | None = None, decreasing: bool = False
) -> pd.DataFrame: ...

class SubSeqElistChisq(SubSeqElist):
    labels: list[str]
    bonferroni: dict
    def __init__(
        self,
        eseq,
        constraint,
        subseq_objs,
        data,
        type: str = "chisq",
        labels=None,
        bonferroni=None,
    ) -> None: ...

def seqecmpgroup(
    subseq: SubSeqElist,
    group,
    method: str = "chisq",
    pvalue_limit: float | None = None,
    weighted: bool = True,
) -> SubSeqElistChisq: ...
def seqedist(
    eseq,
    idcost,
    vparam: float,
    interval: str | int | bool = "No",
    norm: str | int | bool = "YujianBo",
) -> np.ndarray: ...
def frequent_events(
    events: EventSequence, min_support: int = 1
) -> list[tuple[str, int, float]]: ...
def events_transition_matrix(
    events: EventSequence, normalize: bool = False, id_weights: Mapping[str, float] | None = None
) -> tuple[list[str], np.ndarray]: ...
@dataclass
class EventSequence:
    id: list[str]
    time: np.ndarray
    label: list[str]

def seqelength(eseq=..., *, s=...) -> np.ndarray | float: ...
def seqentrans(fsubseq, *, avg_occ: bool = False): ...
def events_to_state_sequences(
    events: EventSequence,
    *,
    states: list[str] | None = None,
    time_grid: np.ndarray | None = None,
    fill: str = "_",
    align_missing: bool = False,
) -> STSList: ...
def frequent_ngrams(
    x: STSList, n: int = 2, min_support: int = 1
) -> list[tuple[tuple[str, ...], int, float]]: ...
def frequent_event_pairs(
    events: EventSequence, min_support: int = 1
) -> list[tuple[tuple[str, str], int, float]]: ...
def events_transition_degrees(
    events: EventSequence, normalize: bool = False, id_weights: Mapping[str, float] | None = None
) -> tuple[list[str], np.ndarray, np.ndarray]: ...
def frequent_event_pairs_timegap(
    events: EventSequence, max_gap: int, min_support: int = 1
) -> list[tuple[tuple[str, str], int, float]]: ...
