from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from .sequences import STSList as STSList
from .stats import seqistatd as seqistatd

@dataclass
class SeqTreeResult:
    model: object
    feature_names: list[str]
    labels: np.ndarray
    leaf_index: np.ndarray

def seqtree_fit(
    x: STSList,
    y,
    *,
    features: str = "durations+transitions+summary",
    max_depth: int | None = None,
    min_samples_leaf: int = 2,
    random_state: int = 0,
) -> SeqTreeResult: ...
def seqtree_leaf_summary(
    x: STSList | SeqTreeResult | None = None,
    result: SeqTreeResult | None = None,
) -> pd.DataFrame: ...
def seqtree_plot(result: SeqTreeResult, ax=None): ...
def seqtree_cv_fit(
    x: STSList,
    y,
    *,
    features: str = "durations+transitions+summary",
    param_grid: dict[str, Any] | None = None,
    cv: int = 5,
    scoring: str | None = None,
    random_state: int = 0,
): ...
def seqtree_predict(
    result: SeqTreeResult, x: STSList, *, features: str = "durations+transitions+summary"
) -> np.ndarray: ...
