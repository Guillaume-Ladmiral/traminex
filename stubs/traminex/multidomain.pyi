import numpy as np

from .sequences import SeqMD as SeqMD

def seqMD_gridsearch_weights(
    obj: SeqMD,
    groups: list[int] | np.ndarray,
    *,
    method: str = "LCS",
    weight_grid: list[list[float]] | None = None,
    weight_step: float = 0.1,
    metric: str = "silhouette",
    per_channel_kwargs: list[dict] | None = None,
    return_all: bool = False,
) -> list[float] | dict: ...
def seqMD_optimize_weights(
    obj: SeqMD,
    groups: list[int] | np.ndarray,
    *,
    method: str = "LCS",
    metric: str = "silhouette",
    per_channel_kwargs: list[dict] | None = None,
    init_weights: list[float] | None = None,
    max_iter: int = 50,
    tol: float = 0.001,
) -> dict: ...
