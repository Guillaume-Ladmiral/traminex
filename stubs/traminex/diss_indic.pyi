from typing import Optional

import numpy as np
import pandas as pd

def dissindic(
    diss: np.ndarray,
    group: np.ndarray,
    *,
    gower: bool = False,
    squared: bool = False,
    weights: Optional[np.ndarray] = None,
) -> pd.DataFrame: ...
