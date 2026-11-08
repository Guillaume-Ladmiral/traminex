from dataclasses import dataclass
from typing import Literal

import numpy as np
import pandas as pd

@dataclass
class ClusterResult:
    labels_: np.ndarray
    medoid_index: int | None = ...

def medoid_index(D: np.ndarray) -> int: ...
def cluster_agglomerative(D: np.ndarray, n_clusters: int) -> ClusterResult: ...
def cluster_summary(x, labels: np.ndarray, D: np.ndarray | None = None) -> pd.DataFrame: ...
def cluster_representatives(
    x, labels: np.ndarray, D: np.ndarray, *, sep: str = "-"
) -> pd.DataFrame: ...
def seqrep(
    x,
    *,
    diss: np.ndarray | None = None,
    criterion: Literal["density", "freq", "dist", "prob", "mscore", "random"] = "density",
    score: np.ndarray | None = None,
    decreasing: bool = True,
    coverage: float = 0.25,
    nrep: int | None = None,
    pradius: float = 0.10,
    dmax: float | None = None,
    weighted: bool = True,
    trep: float = ...,
    tsim: float = ...,
    dist_matrix: np.ndarray | None = ...,
    random_state: int | np.random.Generator | None = None,
    **dist_kwargs,
) -> dict: ...
def seqrep_grp(
    x,
    group=None,
    *,
    diss: np.ndarray | None = None,
    ret: str = "stat",
    criterion: str = "density",
    coverage: float = 0.25,
    weighted: bool = True,
    with_missing: bool = False,
    mdis: np.ndarray | None = ...,
    **dist_kwargs,
): ...
