"""OMloc (Localized OM) - position-dependent indel and/or substitution costs."""

from __future__ import annotations

import numpy as np

from .costs import local_indel_by_pos, trate_substitution_matrix_timevarying


def _validate_substitution_array(
    value: np.ndarray,
    *,
    length: int,
    state_count: int,
    timevarying: bool,
) -> np.ndarray:
    """Return a finite localized-OM substitution array with its contract checked."""
    array = np.asarray(value, dtype=float)
    if timevarying:
        if array.ndim not in {2, 3}:
            raise ValueError("subm_base array must be 2D [K,K] or 3D [L,K,K]")
        expected = (
            (state_count, state_count) if array.ndim == 2 else (length, state_count, state_count)
        )
        if array.shape != expected:
            raise ValueError(f"subm_base array has shape {array.shape}; expected {expected}")
    else:
        if array.ndim == 3:
            raise ValueError("a 3D subm_base requires subm_timevarying=True")
        if array.ndim != 2:
            raise ValueError("subm_base array must be 2D [K,K]")
        expected = (state_count, state_count)
        if array.shape != expected:
            raise ValueError(f"subm_base array has shape {array.shape}; expected {expected}")

    if not np.all(np.isfinite(array)):
        raise ValueError("subm_base costs must be finite")
    if np.any(array < 0):
        raise ValueError("subm_base costs must be non-negative")
    if not np.all(np.diagonal(array, axis1=-2, axis2=-1) == 0):
        raise ValueError("subm_base must have a zero diagonal")
    return array


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
) -> np.ndarray:
    """OMloc: OM with localized (position-dependent) indel and/or substitution costs.

    This is a convenience wrapper that:
    1. Computes per-position indel costs using `local_indel_by_pos()`
    2. Optionally computes time-varying substitution matrix [L,K,K]
    3. Calls `pairwise()` with `indel_by_pos` and/or `subm_by_pos`

    Parameters
    ----------
    seqmat : np.ndarray
        Non-empty integer-coded sequence matrix of shape ``(n, L)`` with
        non-negative state codes. Array-like integer input is normalized with
        ``np.asarray``.

    indel_method : str, default "change"
        Method for computing local indel costs:
        - "change": based on transition probability at each position
        - "entropy": based on normalized entropy at each position

    indel_base : float, default 1.0
        Base indel cost before localization

    indel_scale : float, default 1.0
        Scaling factor for the localization signal

    indel_smooth : int | None, default None
        Smoothing window size for indel profile (moving average)

    subm_timevarying : bool, default False
        If True, compute position-specific substitution matrix [L,K,K]
        If False, use a single global substitution matrix

    subm_base : str | np.ndarray | float, default "TRATE"
        Base substitution specification:
        - "TRATE": compute from transition rates (time-varying if subm_timevarying=True)
        - "CONSTANT": constant off-diagonal costs
        - np.ndarray: use ``(K, K)`` globally, or ``(K, K)`` / ``(L, K, K)``
          when ``subm_timevarying=True``
        - float: scalar substitution cost

    trate_smoothing : float, default 0.0 (R-faithful; no smoothing)
        Smoothing parameter for TRATE computation

    trate_mode : str, default "symmetric"
        TRATE mode: "symmetric", "row", "col"

    trate_normalize : str, default "none"
        TRATE normalization: "none", "mean_offdiag", "clip_0_2"

    norm : str | None, default "length"
        Distance normalization: None, "length", "max"

    **kwargs
        Additional arguments passed to `pairwise()` (e.g., deduplicate,
        condensed). ``trate_weights`` is consumed by this wrapper when
        ``subm_base="TRATE"`` and is not forwarded to regular OM.

    Returns
    -------
    np.ndarray
        Float distance output:
        - square matrix of shape ``(n, n)``
        - or condensed vector of shape ``(n*(n-1)/2,)`` when ``condensed=True``

    Examples
    --------
    >>> import numpy as np
    >>> from traminex.distances.omloc import omloc
    >>> seqmat = np.array([[0,0,1,1,2], [0,1,1,2,2], [2,2,1,0,0]])
    >>> # Position-dependent indel only
    >>> D1 = omloc(seqmat, subm_timevarying=False)
    >>> # Both position-dependent indel and substitution
    >>> D2 = omloc(seqmat, subm_timevarying=True)

    Notes
    -----
    - The `indel_method="change"` is recommended for typical life-course data
      where transitions are non-uniform across time.
    - Time-varying substitution matrices significantly increase computation time
      but can improve distance accuracy when transition rates vary by position.
    - For sequences of length L with K states, time-varying substitution requires
      storing L matrices of size K×K.
    """
    seqmat = np.asarray(seqmat)
    if seqmat.ndim != 2:
        raise ValueError("seqmat must be 2D [n,L]")
    if np.issubdtype(seqmat.dtype, np.bool_) or not np.issubdtype(seqmat.dtype, np.integer):
        raise ValueError("seqmat must use an integer dtype")
    if seqmat.shape[0] == 0:
        raise ValueError("seqmat must contain at least one sequence")
    if seqmat.shape[1] == 0:
        raise ValueError("seqmat sequences must have positive length")
    if np.any(seqmat < 0):
        raise ValueError("seqmat state codes must be non-negative")

    if not np.isfinite(indel_base):
        raise ValueError("indel_base must be finite")
    if indel_base <= 0:
        raise ValueError("indel_base must be strictly positive")
    if not np.isfinite(indel_scale):
        raise ValueError("indel_scale must be finite")
    if indel_scale < 0:
        raise ValueError("indel_scale must be non-negative")
    if indel_smooth is not None:
        if (
            isinstance(indel_smooth, (bool, np.bool_))
            or not isinstance(indel_smooth, (int, np.integer))
            or indel_smooth <= 0
        ):
            raise ValueError("indel_smooth must be a positive integer or None")
        indel_smooth = int(indel_smooth)
    if not isinstance(subm_timevarying, (bool, np.bool_)):
        raise TypeError("subm_timevarying must be bool")
    subm_timevarying = bool(subm_timevarying)
    if not np.isfinite(trate_smoothing):
        raise ValueError("trate_smoothing must be finite")
    if trate_smoothing < 0:
        raise ValueError("trate_smoothing must be non-negative")

    length = seqmat.shape[1]
    state_count = int(seqmat.max()) + 1
    trate_weights_arg = kwargs.pop("trate_weights", None)
    is_trate = isinstance(subm_base, str) and subm_base.upper() == "TRATE"
    trate_weights: np.ndarray | None = None
    if trate_weights_arg is not None:
        if not is_trate:
            raise ValueError("trate_weights is only valid with subm_base='TRATE'")
        trate_weights = np.asarray(trate_weights_arg, dtype=float)
        if trate_weights.ndim != 1:
            raise ValueError("trate_weights must be one-dimensional")
        if trate_weights.shape[0] != seqmat.shape[0]:
            raise ValueError("trate_weights length must equal number of sequences")
        if not np.all(np.isfinite(trate_weights)):
            raise ValueError("trate_weights must be finite")
        if np.any(trate_weights < 0):
            raise ValueError("trate_weights must be non-negative")
        if float(trate_weights.sum()) <= 0:
            raise ValueError("trate_weights must have a positive total")

    # Compute position-dependent indel costs
    indel_by_pos = local_indel_by_pos(
        seqmat,
        method=indel_method,
        base=indel_base,
        scale=indel_scale,
        smooth_window=indel_smooth,
    )

    # Determine substitution matrix
    subm_by_pos: np.ndarray | None = None
    subm: float | np.ndarray | None = None

    if subm_timevarying:
        # Compute time-varying substitution
        if isinstance(subm_base, str):
            if subm_base.upper() == "TRATE":
                subm_by_pos = trate_substitution_matrix_timevarying(
                    seqmat,
                    smoothing=trate_smoothing,
                    mode=trate_mode,
                    normalize=trate_normalize,
                    weights=trate_weights,
                )
            elif subm_base.upper() == "CONSTANT":
                # Time-varying constant matrix (same value at all positions)
                subm_by_pos = np.full((length, state_count, state_count), 2.0, dtype=float)
                for matrix in subm_by_pos:
                    np.fill_diagonal(matrix, 0.0)
            else:
                raise ValueError(f"Unknown subm_base: {subm_base}")
        elif isinstance(subm_base, np.ndarray):
            arr = _validate_substitution_array(
                subm_base,
                length=length,
                state_count=state_count,
                timevarying=True,
            )
            if arr.ndim == 3:  # Already [L,K,K]
                subm_by_pos = arr
            else:  # [K,K] - replicate for all positions
                subm_by_pos = np.tile(arr[None, :, :], (length, 1, 1))
        else:
            # Scalar - replicate as constant matrix for all positions
            scalar = float(subm_base)
            if not np.isfinite(scalar):
                raise ValueError("subm_base scalar cost must be finite")
            if scalar < 0:
                raise ValueError("subm_base scalar cost must be non-negative")
            subm_by_pos = np.full((length, state_count, state_count), scalar, dtype=float)
            for matrix in subm_by_pos:
                np.fill_diagonal(matrix, 0.0)
    else:
        # Global substitution matrix (not time-varying)
        if isinstance(subm_base, str):
            if subm_base.upper() == "TRATE":
                from .costs import trate_substitution_matrix

                subm = trate_substitution_matrix(
                    seqmat,
                    smoothing=trate_smoothing,
                    mode=trate_mode,
                    normalize=trate_normalize,
                    weights=trate_weights,
                )
            elif subm_base.upper() == "CONSTANT":
                subm = np.full((state_count, state_count), 2.0, dtype=float)
                np.fill_diagonal(subm, 0.0)
            else:
                raise ValueError(f"Unknown subm_base: {subm_base}")
        elif isinstance(subm_base, np.ndarray):
            subm = _validate_substitution_array(
                subm_base,
                length=length,
                state_count=state_count,
                timevarying=False,
            )
        else:
            subm = float(subm_base)
            if not np.isfinite(subm):
                raise ValueError("subm_base scalar cost must be finite")
            if subm < 0:
                raise ValueError("subm_base scalar cost must be non-negative")

    # Call pairwise with localized parameters
    # Import here to avoid circular import
    from . import pairwise as pairwise_func

    # Use mean of indel_by_pos as the scalar indel (for border initialization)
    # Even though we use indel_by_pos, the pairwise function still needs a scalar value
    indel_scalar = float(np.mean(indel_by_pos))

    return pairwise_func(
        seqmat,
        method="OM",
        indel=indel_scalar,
        indel_by_pos=indel_by_pos,
        subm=subm,
        subm_by_pos=subm_by_pos,
        norm=norm,
        **kwargs,
    )
