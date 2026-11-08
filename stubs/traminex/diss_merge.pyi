from typing import Any, Optional

import numpy as np

def dissmergegroups(
    diss: np.ndarray,
    group: np.ndarray,
    *,
    weights: Optional[np.ndarray] = None,
    measure: str = "ASW",
    crit: float = 0.2,
    ref: str = "max",
    min_group: int = 4,
    small: float = 0.05,
    squared: bool = False,
    silent: bool = True,
) -> dict[str, Any]: ...
