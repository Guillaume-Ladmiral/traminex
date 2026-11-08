from typing import Any

import numpy as np

from .disstree import DissTree as DissTree
from .disstree_viz import seqtree2dot as seqtree2dot
from .disstree_viz import seqtreedisplay as seqtreedisplay
from .seqtree_explore import (
    SeqTreeResult as SeqTreeResult,
)
from .seqtree_explore import (
    seqtree_cv_fit as seqtree_cv_fit,
)
from .seqtree_explore import (
    seqtree_fit as seqtree_fit,
)
from .seqtree_explore import (
    seqtree_leaf_summary as seqtree_leaf_summary,
)
from .seqtree_explore import (
    seqtree_plot as seqtree_plot,
)
from .seqtree_explore import (
    seqtree_predict as seqtree_predict,
)
from .sequences import STSList

def seqtree(
    seqdata: STSList,
    predictors: Any,
    *,
    diss: np.ndarray | None = None,
    weighted: bool = True,
    seqdist_args: dict[str, Any] | None = None,
    min_size: float = 0.05,
    max_depth: int = 5,
    R: int = 1000,
    pval: float = 0.01,
    weight_permutation: str = "replicate",
    squared: bool = False,
    first: str | None = None,
    random_state: int | None = None,
    minSize: float = ...,
    maxdepth: int = ...,
    seqdist_arg: dict[str, Any] | None = ...,
) -> DissTree: ...
