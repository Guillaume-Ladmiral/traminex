from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np
import pandas as pd

@dataclass
class DissMultifactor:
    variables: list[str]
    pseudo_F: np.ndarray
    pseudo_R2: np.ndarray
    p_value: np.ndarray
    SS: np.ndarray
    mfac: pd.DataFrame

def dissmfacw(
    diss: np.ndarray,
    factors: Mapping[str, Sequence],
    weights: Sequence[float] | None = None,
    *,
    nperm: int = 1000,
    seed: int | None = None,
    gower: bool = False,
    squared: bool = False,
) -> DissMultifactor: ...
