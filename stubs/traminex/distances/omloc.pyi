import numpy as np

from .costs import (
    local_indel_by_pos as local_indel_by_pos,
)
from .costs import (
    trate_substitution_matrix_timevarying as trate_substitution_matrix_timevarying,
)

def omloc(
    seqmat: np.ndarray,
    *,
    indel_method: str = "change",
    indel_base: float = 1.0,
    indel_scale: float = 1.0,
    indel_smooth: int | None = None,
    subm_timevarying: bool = False,
    subm_base: str | np.ndarray | float = "TRATE",
    trate_smoothing: float = 0.0,
    trate_mode: str = "symmetric",
    trate_normalize: str = "none",
    norm: str | None = "length",
    **kwargs,
) -> np.ndarray: ...
