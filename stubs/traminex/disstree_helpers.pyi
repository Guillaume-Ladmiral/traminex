from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd

def disstreeleaf(tree: Any, label: bool = False, collapse: str = ", ") -> np.ndarray: ...
def disstree_get_rules(tree: Any, collapse: str = "; ") -> list[str | None]: ...
def disstree_assign(
    rules: Sequence[str | None],
    profile: pd.DataFrame | pd.Series | dict,
    covar: Sequence[str] | None = None,
) -> np.ndarray: ...
