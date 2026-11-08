from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

@dataclass
class DissTreeResult:
    model: Any
    coords: np.ndarray
    labels: np.ndarray
    eigvecs: np.ndarray
    eigvals: np.ndarray
    row_mean_d2: np.ndarray
    grand_mean_d2: float

def disstree_fit(
    D: np.ndarray,
    y,
    *,
    n_components: int = 5,
    max_depth: int | None = None,
    min_samples_leaf: int = 2,
    random_state: int = 0,
) -> DissTreeResult: ...
def disstree_plot(result: DissTreeResult, ax=None): ...
def disstree_leaf_summary(result: DissTreeResult) -> pd.DataFrame: ...
def disstree_feature_importance(result: DissTreeResult) -> pd.DataFrame: ...
def disstree_predict(result: DissTreeResult, D_new: np.ndarray) -> np.ndarray: ...
def disstree2dot_explore(result: DissTreeResult, out_file: str) -> None: ...
