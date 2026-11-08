from collections.abc import Sequence

import numpy as np
import pandas as pd

from .sequences import STSList

class SAMMResult:
    frame: pd.DataFrame
    sequences: STSList
    spell: tuple[str, ...]
    typology: np.ndarray
    state_columns: tuple[str, ...]
    sublength: int
    def __init__(
        self,
        frame: pd.DataFrame,
        sequences: STSList,
        spell: tuple[str, ...],
        typology: np.ndarray,
        state_columns: tuple[str, ...],
        sublength: int,
    ) -> None: ...

def seqsamm(
    x: STSList,
    sublength: int,
    covar: pd.DataFrame | None = None,
) -> SAMMResult: ...
def seqsammseq(samm: SAMMResult, spell: str) -> STSList: ...
def seqsammeha(
    samm: SAMMResult,
    spell: str,
    typology: Sequence[object],
    persper: bool = True,
) -> pd.DataFrame: ...
