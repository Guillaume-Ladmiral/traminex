from typing import Any, Literal, Optional

import numpy as np

def dissrf(
    diss: np.ndarray,
    *,
    k: int | None = None,
    sortv=...,
    weights: np.ndarray | None = None,
    grp_meth: Literal["prop", "first", "random"] = "prop",
    squared: bool = False,
    pow: int | None = None,
    random_state: int | np.random.Generator | None = None,
) -> dict[str, Any]: ...
def seqrf(
    x,
    diss: np.ndarray,
    *,
    k: Optional[int] = None,
    sortv=...,
    weights: np.ndarray | None = None,
    weighted: bool = True,
    grp_meth: Literal["prop", "first", "random"] = "prop",
    squared: bool = False,
    pow: int | None = None,
    random_state: int | np.random.Generator | None = None,
) -> dict[str, Any]: ...
