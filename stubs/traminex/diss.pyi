from dataclasses import dataclass
from typing import Any, Literal, Sequence

import numpy as np

def double_center_matrix(D: np.ndarray) -> np.ndarray: ...
def disscenter_double_center(D: np.ndarray) -> np.ndarray: ...
def disscenter(
    diss: np.ndarray,
    group: np.ndarray | None = None,
    medoids_index: str | bool | None = None,
    allcenter: bool = False,
    weights: np.ndarray | None = None,
    squared: bool = False,
) -> Any: ...
def gower_matrix(
    D: np.ndarray, *, squared: bool = ..., weights: np.ndarray | None = ...
) -> np.ndarray: ...
def dissvar_grp(
    diss: np.ndarray,
    group: np.ndarray | None = None,
    *,
    weights: np.ndarray | None = None,
    squared: bool = False,
) -> dict[str, float]: ...
def weighted_pairwise_sum(
    diss: np.ndarray,
    indices: np.ndarray | None = None,
    weights: np.ndarray | None = None,
) -> float: ...
def weighted_inertia(
    diss: np.ndarray,
    indices: np.ndarray | None = None,
    weights: np.ndarray | None = None,
) -> float: ...
def dissvar(
    diss: np.ndarray, weights: np.ndarray | None = None, squared: bool = False
) -> float: ...
def dissvar_summary(D: np.ndarray, labels: np.ndarray | None = None) -> dict[str, Any]: ...
@dataclass
class DissAssocResult:
    groups: Any
    anova_table: Any
    stat: Any
    R: int
    weight_permutation: str
    perms: dict[str, Any] = ...

def dissassoc(
    diss: np.ndarray,
    group: np.ndarray,
    weights: np.ndarray | None = None,
    R: int = 1000,
    weight_permutation: Literal[
        "none",
        "diss",
        "replicate",
        "rounded-replicate",
        "group",
        "permute",
        "random-sampling",
    ] = "replicate",
    squared: bool = False,
    random_state: int | np.random.Generator | None = None,
) -> DissAssocResult: ...
def dissassoc_effect(D: np.ndarray, labels: np.ndarray) -> dict[str, Any]: ...
def dissassoc_test(
    D: np.ndarray, labels: np.ndarray, *, n_perm: int = 999, random_state: int | None = None
) -> dict[str, Any]: ...
def dissdomassoc(
    domain_dists: list[np.ndarray],
    *,
    jointdiss: np.ndarray | None = None,
    what: Sequence[str] = ...,
    weights: np.ndarray | None = None,
    dnames: Sequence[str] | None = None,
    w_rank: bool = False,
) -> dict[str, Any]: ...
