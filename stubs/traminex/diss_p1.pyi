from typing import Literal

import numpy as np

def dissrep(
    diss: np.ndarray,
    *,
    criterion: Literal["density", "freq", "dist", "random"] = "density",
    score: np.ndarray | None = None,
    decreasing: bool = True,
    coverage: float = 0.25,
    nrep: int | None = None,
    pradius: float = 0.1,
    dmax: float | None = None,
    weights: np.ndarray | None = None,
    random_state: int | np.random.Generator | None = None,
    trep: float = ...,
    tsim: float = ...,
) -> dict: ...
