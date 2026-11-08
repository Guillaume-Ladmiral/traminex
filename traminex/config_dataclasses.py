"""Configuration dataclasses for supported distance parameters.

Provides structured, type-safe configuration for TRATE substitution matrices,
the implemented DHD variants, and other runtime-supported parameter sets.
"""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral, Real
from typing import Optional, Sequence

import numpy as np

#: seqdist's cost methods for ``subm`` (seqdist.R:151).
_SUBM_METHODS = ("TRATE", "CONSTANT", "INDELS", "INDELSLOG")


def _check_seqdist_norm(norm: object) -> None:
    """Refuse a ``norm`` that seqdist refuses for OM, LCS, LCP and RLCP."""
    if norm is None:
        return
    from .distances import _CLASSIC_NORM_DISPLAY, _CLASSIC_NORM_SPELLINGS

    if not isinstance(norm, str) or norm.strip().lower() not in _CLASSIC_NORM_SPELLINGS:
        accepted = ", ".join(_CLASSIC_NORM_DISPLAY["OM"])
        raise ValueError(f"norm must be one of {accepted}, got {norm!r}")


@dataclass
class TRATEConfig:
    """Configuration for TRATE substitution matrix computation.

    Parameters
    ----------
    smoothing : float, default=0.0
        Smoothing factor added to transition rates. The R-faithful default is
        0.0 (TraMineR's TRATE route applies no smoothing); pass a positive
        value only to deliberately regularize sparse transition tables.
    mode : str, default="symmetric"
        Transition-cost construction mode:
        - "symmetric": 2 - (T[i,j] + T[j,i])
        - "row": 1 - T[i,j]
        - "col": 1 - T[j,i]
    normalize : str, default="none"
        Normalization strategy:
        - "none": No normalization (raw costs)
        - "mean_offdiag": Scale by mean off-diagonal cost
        - "clip_0_2": Clip costs to [0, 2]
    weights : Optional[np.ndarray], default=None
        Per-sequence weights for computing transition rates.

    Examples
    --------
    >>> import numpy as np
    >>> from traminex import seqdef, seqdist
    >>> sequences = seqdef(np.array([list("AABBC"), list("ABBCC"), list("CCBBA")], dtype=object))
    >>> from traminex.distances import trate_substitution_matrix
    >>> config = TRATEConfig(smoothing=1e-6, mode="symmetric", normalize="mean_offdiag")
    >>> sm = trate_substitution_matrix(sequences.data, **config.as_dict())
    """

    smoothing: float = 0.0
    mode: str = "symmetric"
    normalize: str = "none"
    weights: Optional[np.ndarray] = None

    def as_dict(self) -> dict:
        """Convert to kwargs dictionary for trate_substitution_matrix."""
        return {
            "smoothing": self.smoothing,
            "mode": self.mode,
            "normalize": self.normalize,
            "weights": self.weights,
        }

    def __post_init__(self):
        """Validate configuration parameters."""
        if isinstance(self.smoothing, (bool, np.bool_)) or not isinstance(
            self.smoothing, (Real, np.number)
        ):
            raise TypeError(f"smoothing must be a number, got {self.smoothing!r}")
        if not np.isfinite(self.smoothing) or self.smoothing < 0:
            raise ValueError(f"smoothing must be non-negative and finite, got {self.smoothing}")

        valid_modes = {"symmetric", "row", "col"}
        if self.mode not in valid_modes:
            raise ValueError(f"mode must be one of {valid_modes}, got {self.mode!r}")

        valid_norms = {"none", "mean_offdiag", "clip_0_2"}
        if self.normalize not in valid_norms:
            raise ValueError(f"normalize must be one of {valid_norms}, got {self.normalize!r}")


@dataclass
class DHDConfig:
    """Configuration for Dynamic Hamming Distance computation.

    ``DHDConfig()`` selects ``seqdist(method="DHD")``'s default, TraMineR's
    transition-rate DHD. Setting ``variant``, ``smooth_window`` or
    ``profile="calibrated"`` opts into traminex's legacy marginal-frequency
    metric, which is not TraMineR's DHD.

    Parameters
    ----------
    variant : str or None, default=None
        ``None`` (default) leaves the DHD metric to ``seqdist``: TraMineR's.
        A legacy marginal-frequency variant is selected by:
        - "mean": cost_t(x,y) = 1 - 0.5 * (p_t[x] + p_t[y])
        - "sum":  cost_t(x,y) = 2 - (p_t[x] + p_t[y]); equals 2 × "mean"
    smooth_window : Optional[int], default=None
        Moving average window size for smoothing per-position frequencies.
        If None, no smoothing is applied.
    profile : str, default="default"
        Configuration profile to use:
        - "default": Use specified variant and smooth_window
        - "calibrated": Load from explicit calibration path(s) passed to seqdist
          (``calibrated_params_path``/``calibrated_params_paths``) or from
          env var ``TRAMINEX_CALIBRATED_PARAMS``; with ``variant`` (or
          ``smooth_window``) left unset, the calibration file's value is used.
          A value that is set overrides the file's. When a value must come
          from the file, ``seqdist`` raises ``ValueError`` if no file loads,
          if it has no ``"dhd"`` section or if that section holds a key other
          than ``variant`` and ``smooth_window``.

    Examples
    --------
    >>> import numpy as np
    >>> from traminex import seqdef, seqdist
    >>> sequences = seqdef(np.array([list("AABBC"), list("ABBCC"), list("CCBBA")], dtype=object))
    >>> D = seqdist(sequences, method="DHD", **DHDConfig().as_dict())  # TraMineR DHD
    >>> legacy = DHDConfig(variant="mean", smooth_window=5)
    >>> D_legacy = seqdist(sequences, method="DHD", **legacy.as_dict())
    """

    variant: Optional[str] = None
    smooth_window: Optional[int] = None
    profile: str = "default"

    def as_dict(self) -> dict:
        """Convert to kwargs dictionary for DHD distance computation."""
        # The dispatcher selects the legacy metric whenever ``variant`` is
        # present, so it is emitted only when set.
        result: dict[str, object] = {"dhd_profile": self.profile}
        if self.variant is not None:
            result["variant"] = self.variant
        if self.smooth_window is not None:
            result["smooth_window"] = self.smooth_window
        return result

    def __post_init__(self):
        """Validate configuration parameters."""
        valid_variants = {"mean", "sum"}
        if self.variant is not None:
            # The spellings seqdist accepts (normalize_dhd_variant; TX-AUD-20261019-026).
            spelled = str(self.variant).lower()
            if spelled not in valid_variants:
                raise ValueError(f"variant must be one of {valid_variants}, got {self.variant!r}")
            self.variant = spelled

        if self.smooth_window is not None:
            if isinstance(self.smooth_window, bool) or not isinstance(self.smooth_window, Integral):
                raise TypeError("smooth_window must be an integer")
            if self.smooth_window < 1:
                raise ValueError(f"smooth_window must be positive, got {self.smooth_window}")

        valid_profiles = {"default", "calibrated"}
        if self.profile not in valid_profiles:
            raise ValueError(f"profile must be one of {valid_profiles}, got {self.profile!r}")


@dataclass
class OMConfig:
    """Configuration for Optimal Matching distance computation.

    Parameters
    ----------
    indel : float | sequence of float | str, default="auto"
        Insertion/deletion cost, as ``seqdist`` takes it: a number, one cost per
        state, or ``"auto"`` (``seqdist``'s default: half the largest
        substitution cost, or the cost-method's own indel; TX-AUD-20261011-043).
    subm : float | np.ndarray | str, default=2.0
        Substitution cost:
        - float: Uniform cost for all mismatches
        - np.ndarray: K×K cost matrix (finite, non-negative, zero diagonal)
        - "TRATE", "CONSTANT", "INDELS" or "INDELSLOG": ``seqdist``'s cost
          methods ("TRATE" with ``trate_config``)
    norm : Optional[str], default=None
        Normalization, any ``seqdist`` OM spelling: "none", "auto",
        "maxlength", "gmean", "maxdist", "YujianBo", or the aliases "length"
        and "max".
    indel_by_pos : Optional[np.ndarray], default=None
        Position-dependent indel costs (array of length L).
    subm_by_pos : Optional[np.ndarray], default=None
        Position-dependent substitution costs (array of shape [L, K, K]).
    trate_config : Optional[TRATEConfig], default=None
        Configuration for TRATE substitution matrix if subm="TRATE".

    Examples
    --------
    >>> import numpy as np
    >>> from traminex import seqdef, seqdist
    >>> sequences = seqdef(np.array([list("AABBC"), list("ABBCC"), list("CCBBA")], dtype=object))
    >>> # Simple OM with uniform costs
    >>> config = OMConfig(indel=1.0, subm=2.0, norm="length")
    >>> D = seqdist(sequences, method="OM", **config.as_dict())

    >>> # OM with TRATE substitution matrix
    >>> trate = TRATEConfig(smoothing=1e-6, mode="symmetric")
    >>> config = OMConfig(subm="TRATE", trate_config=trate)
    >>> D = seqdist(sequences, method="OM", **config.as_dict())
    """

    indel: float | Sequence[float] | np.ndarray | str = "auto"
    subm: float | np.ndarray | str = 2.0
    norm: Optional[str] = None
    indel_by_pos: Optional[np.ndarray] = None
    subm_by_pos: Optional[np.ndarray] = None
    trate_config: Optional[TRATEConfig] = None

    def as_dict(self) -> dict:
        """Convert to kwargs dictionary for OM distance computation."""
        indel = self.indel
        if not isinstance(indel, str) and np.ndim(indel) == 1:
            indel = np.asarray(indel, dtype=float)
        result = {
            "indel": indel,
            "subm": self.subm,
        }
        if self.norm is not None:
            result["norm"] = self.norm
        if self.indel_by_pos is not None:
            result["indel_by_pos"] = self.indel_by_pos
        if self.subm_by_pos is not None:
            result["subm_by_pos"] = self.subm_by_pos

        # Add TRATE config if subm="TRATE"
        if (
            isinstance(self.subm, str)
            and self.subm.lower() == "trate"
            and self.trate_config is not None
        ):
            trate_kwargs = self.trate_config.as_dict()
            result.update({f"trate_{k}": v for k, v in trate_kwargs.items() if k != "weights"})
            if trate_kwargs.get("weights") is not None:
                result["trate_weights"] = trate_kwargs["weights"]

        return result

    def __post_init__(self):
        """Validate configuration parameters."""
        # A boolean is an int to Python and was silently used as 1 or 0
        # (audit S2-17).
        # The values seqdist takes, refused as seqdist refuses them (TX-AUD-20261011-043).
        if isinstance(self.indel, (bool, np.bool_)):
            raise TypeError("indel must be a number, not a boolean")
        if isinstance(self.subm, (bool, np.bool_)):
            raise TypeError("subm must be a number, an array or 'TRATE', not a boolean")
        if isinstance(self.indel, str):
            if self.indel.lower() != "auto":
                raise ValueError(
                    "indel must be a number, a per-state vector of numbers or 'auto', "
                    f"got {self.indel!r}"
                )
        else:
            try:
                indel_values = np.asarray(self.indel, dtype=float)
            except (TypeError, ValueError):
                raise ValueError(
                    f"indel must be a number or a per-state vector of numbers, got {self.indel!r}"
                ) from None
            if indel_values.ndim > 1:
                raise ValueError("indel must be a number or a one-dimensional per-state vector")
            if (
                indel_values.size == 0
                or not np.all(np.isfinite(indel_values))
                or np.any(indel_values <= 0)
            ):
                raise ValueError(f"indel must be strictly positive, got {self.indel}")

        if self.trate_config is not None and not (
            isinstance(self.subm, str) and self.subm.upper() == "TRATE"
        ):
            # as_dict would drop it, as seqdist refuses TRATE-specific keywords without
            # subm="TRATE" (TX-AUD-20261019-026).
            raise TypeError("trate_config requires subm='TRATE'")
        if isinstance(self.subm, str):
            if self.subm.upper() not in _SUBM_METHODS:
                raise ValueError(
                    "subm must be a number, a K x K matrix or one of "
                    f"{', '.join(_SUBM_METHODS)}, got {self.subm!r}"
                )
        elif np.ndim(self.subm) == 0:
            try:
                cost = float(self.subm)
            except (TypeError, ValueError):
                raise TypeError(
                    "subm must be a number, a K x K matrix or one of "
                    f"{', '.join(_SUBM_METHODS)}, got {self.subm!r}"
                ) from None
            if not np.isfinite(cost) or cost < 0:
                raise ValueError(f"subm must be non-negative, got {self.subm}")
        else:
            try:
                matrix = np.asarray(self.subm, dtype=float)
            except (TypeError, ValueError):
                raise ValueError("subm must be a K x K matrix of numbers") from None
            if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
                raise ValueError("subm must be a square K x K matrix")
            if not np.all(np.isfinite(matrix)) or np.any(matrix < 0):
                raise ValueError("subm costs must be finite and non-negative")
            if np.any(np.diag(matrix) != 0):
                raise ValueError("subm must have a zero diagonal")

        _check_seqdist_norm(self.norm)


@dataclass
class LCSLCPConfig:
    """Configuration for LCS/LCP distance computation.

    Parameters
    ----------
    norm : Optional[str], default=None
        Normalization, any ``seqdist`` spelling for these methods: "none",
        "auto", "maxlength", "gmean", "maxdist", "YujianBo", "length" or "max".

    Examples
    --------
    >>> import numpy as np
    >>> from traminex import seqdef, seqdist
    >>> sequences = seqdef(np.array([list("AABBC"), list("ABBCC"), list("CCBBA")], dtype=object))
    >>> config = LCSLCPConfig(norm="length")
    >>> D_lcs = seqdist(sequences, method="LCS", **config.as_dict())
    >>> D_lcp = seqdist(sequences, method="LCP", **config.as_dict())
    """

    norm: Optional[str] = None

    def as_dict(self) -> dict:
        """Convert to kwargs dictionary for LCS/LCP distance computation."""
        result = {}
        if self.norm is not None:
            result["norm"] = self.norm
        return result

    def __post_init__(self):
        """Validate configuration parameters."""
        _check_seqdist_norm(self.norm)


@dataclass
class DURConfig:
    """Configuration for duration-based distance computation.

    Parameters
    ----------
    metric : str, default="l1"
        Distance metric on duration distributions:
        - "l1": Manhattan distance
        - "cos": Cosine distance
        - "chi2": Chi-square distance

    Examples
    --------
    >>> import numpy as np
    >>> from traminex import seqdef, seqdist
    >>> sequences = seqdef(np.array([list("AABBC"), list("ABBCC"), list("CCBBA")], dtype=object))
    >>> config = DURConfig(metric="chi2")
    >>> D = seqdist(sequences, method="DUR", **config.as_dict())
    """

    metric: str = "l1"

    def as_dict(self) -> dict:
        """Convert to kwargs dictionary for DUR distance computation."""
        return {"metric": self.metric}

    def __post_init__(self):
        """Validate configuration parameters."""
        valid_metrics = {"l1", "cos", "chi2"}
        if self.metric not in valid_metrics:
            raise ValueError(f"metric must be one of {valid_metrics}, got {self.metric!r}")


# Convenience function to create default configs
def default_configs() -> dict:
    """Return dictionary of default configurations for all distance methods.

    Returns
    -------
    dict
        Dictionary mapping method names to their default configurations.

    Examples
    --------
    >>> import numpy as np
    >>> from traminex import seqdef, seqdist
    >>> sequences = seqdef(np.array([list("AABBC"), list("ABBCC"), list("CCBBA")], dtype=object))
    >>> configs = default_configs()
    >>> om_config = configs["OM"]
    >>> D = seqdist(sequences, method="OM", **om_config.as_dict())
    """
    return {
        "OM": OMConfig(),
        "LCS": LCSLCPConfig(),
        "LCP": LCSLCPConfig(),
        "RLCP": LCSLCPConfig(),
        "DHD": DHDConfig(),
        "DUR": DURConfig(),
        "TRATE": TRATEConfig(),
    }
