from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from .disstree_explore import (
    DissTreeResult as DissTreeResult,
)
from .disstree_explore import (
    disstree2dot_explore as disstree2dot_explore,
)
from .disstree_explore import (
    disstree_feature_importance as disstree_feature_importance,
)
from .disstree_explore import (
    disstree_fit as disstree_fit,
)
from .disstree_explore import (
    disstree_leaf_summary as disstree_leaf_summary,
)
from .disstree_explore import (
    disstree_plot as disstree_plot,
)
from .disstree_explore import (
    disstree_predict as disstree_predict,
)
from .disstree_viz import (
    disstree2dot as disstree2dot,
)
from .disstree_viz import (
    disstree2dotp as disstree2dotp,
)
from .disstree_viz import (
    disstreedisplay as disstreedisplay,
)

@dataclass
class DissTreeSplit:
    varindex: int
    prob: np.ndarray
    info: dict[str, Any]
    index: np.ndarray | None = ...
    breaks: float | None = ...
    na_group: int | None = ...
    labels: list[Any] | None = ...

@dataclass
class DissTreeNode:
    id: int
    info: dict[str, Any]
    split: DissTreeSplit | None = ...
    kids: list[DissTreeNode] | None = ...

@dataclass
class DissTree:
    root: DissTreeNode
    fitted: np.ndarray
    info: dict[str, Any]
    data: pd.DataFrame
    weights: np.ndarray
    adjustment: Any = ...

def disstreeleaf(
    tree: DissTree | DissTreeNode, label: bool = False, collapse: str = ", "
) -> np.ndarray: ...
def disstree_get_rules(tree: DissTree, collapse: str = "; ") -> list[str | None]: ...
def disstree_assign(
    rules: list[str | None], profile: pd.DataFrame, covar: list[str] | None = None
) -> np.ndarray: ...
def disstree(
    diss: np.ndarray,
    predictors: pd.DataFrame,
    weights: np.ndarray | None = None,
    min_size: float = 0.05,
    max_depth: int = 5,
    R: int = 1000,
    pval: float = 0.01,
    weight_permutation: str = "replicate",
    squared: bool = False,
    first: str | None = None,
    random_state: int | None = None,
    *,
    minSize: float = ...,
    maxdepth: int = ...,
) -> DissTree: ...
