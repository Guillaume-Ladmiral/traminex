from __future__ import annotations

import os
import warnings
from math import isqrt
from numbers import Integral
from typing import Any, Optional

import numpy as np
import pandas as pd

from .._deprecated import _caller_stacklevel
from ..sequences import STSList, _real_state_codes, _validate_sequence_weights
from . import python_kernels as _py
from ._normalization import _NORM_ALIASES, _maxscost, normalize_value, resolve_norm_code
from ._validation import (
    TRATE_OM_KWARGS,
    check_costs_requested,
    is_refused_cost_scalar,
    refuse_unread_calibration_path,
    resolve_r_seqdist_arguments,
    validate_pairwise_kwargs,
)
from .costs import (
    SeqCostMatrix,
    _trate_substitution_matrix_rows,
    trate_substitution_matrix,
)
from .costs import (
    # Re-exported for callers that want the historical mean-off-diagonal rule
    # explicitly; ``om_trate_pairwise`` itself now resolves its automatic indel
    # the way R and ``seqdist`` do (T-06). Declared in
    # ``stubs/traminex/distances/__init__.pyi``.
    indel_from_subm as indel_from_subm,
)
from .costs import (
    local_indel_by_pos as local_indel_by_pos,
)
from .costs import (
    seqcost as seqcost,
)
from .costs import (
    trate_substitution_matrix_timevarying as trate_substitution_matrix_timevarying,
)
from .omloc import omloc as omloc

__all__ = (
    "om_trate_pairwise",
    "pairwise",
    "pairwise_between_sets",
    "seqdist",
    "to_condensed",
    "from_condensed",
    "sanitize_distance_matrix",
    "get_lcp_length",
    "get_rlcp_length",
    "seqsubm",
    "seqLLCS",
    "seqLLCP",
    "seqRLCP",
)

# Try to import compiled backends if present
# Spaces stripped and an empty value read as "auto", as TRAMINEX_NUMBA is read
# (python_kernels._resolve_numba_mode; TX-AUD-20261019-031).
_backend = os.environ.get("TRAMINEX_BINDING", "auto").strip().lower() or "auto"
if _backend not in {"auto", "python", "pybind11"}:
    # Validated once at import, exactly like TRAMINEX_NUMBA (docs/DISTANCES.md):
    # an unknown value must not surface only at the first dispatch call.
    raise ValueError("TRAMINEX_BINDING must be one of: auto, python, pybind11")

_pb = None

if _backend in ("auto", "pybind11"):
    try:
        from .. import _pb_kernels as _pb  # type: ignore
    except ImportError:
        _pb = None


def _impl_name() -> str:
    if _backend not in {"auto", "python", "pybind11"}:
        raise ValueError("TRAMINEX_BINDING must be one of: auto, python, pybind11")
    if _backend == "python":
        return "python"
    if _backend == "pybind11":
        if _pb is None:
            raise RuntimeError("pybind11 backend was requested but is unavailable")
        return "pybind11"
    if _pb is not None:
        return "pybind11"
    return "python"


#: ``seqdist``'s classic OM/LCS/LCP/RLCP normalization vocabulary. R's own six
#: names (``seqdist.R``: auto, none, maxlength, gmean, maxdist, YujianBo) plus the
#: traminex-friendly aliases ``'length'`` (== ``maxlength``) and ``'max'``
#: (== ``maxdist``), accepted by ``seqdist`` for all four methods. The low-level
#: ``_map_norm`` below is narrower: ``'max'`` for OM only.
_CLASSIC_NORM_METHODS = frozenset({"OM", "LCS", "LCP", "RLCP"})
_CLASSIC_NORM_SPELLINGS = frozenset(
    {"length", "max", "none", "maxlength", "gmean", "maxdist", "yujianbo", "auto"}
)
_CLASSIC_NORM_DISPLAY = {
    "OM": ("length", "max", "none", "maxlength", "gmean", "maxdist", "YujianBo", "auto"),
    "LCS": ("length", "max", "none", "maxlength", "gmean", "maxdist", "YujianBo", "auto"),
    "LCP": ("length", "max", "none", "maxlength", "gmean", "maxdist", "YujianBo", "auto"),
    "RLCP": ("length", "max", "none", "maxlength", "gmean", "maxdist", "YujianBo", "auto"),
}


def _map_norm(method: str, norm: Optional[str]) -> Optional[str]:
    """Validate and normalize the `norm` parameter.

    This is the **low-level** ``pairwise`` vocabulary, which is deliberately
    narrower than ``seqdist``'s: only ``None``, ``"length"`` (OM/LCS/LCP/RLCP) and
    ``"max"`` (OM only). ``seqdist`` validates the full TraMineR vocabulary itself
    before dispatch, so a user never sees this message for a spelling ``seqdist``
    accepts.
    """
    if norm is None:
        return None
    n = str(norm).lower()
    m = method.upper()
    if n not in {"length", "max"}:
        raise ValueError(
            "norm must be one of: None, 'length', 'max' for the low-level pairwise() "
            "API; seqdist() additionally accepts 'none', 'maxlength', 'gmean', "
            "'maxdist', 'YujianBo' and 'auto'"
        )
    if m in {"LCS", "LCP", "RLCP"} and n != "length":
        raise ValueError(f"norm='{n}' is not supported for method {m}; use 'length' or None")
    if m not in {"OM", "LCS", "LCP", "RLCP"}:
        raise ValueError(f"norm is not supported for method {m}")
    return n


def _normalize_om_output(
    D: np.ndarray,
    *,
    seqmat: np.ndarray,
    subm: float | np.ndarray,
    indel: float,
    norm: Optional[str],
) -> np.ndarray:
    """Apply OM post-normalization consistently across backends."""
    if norm == "length" and seqmat.shape[1] > 0:
        return D / float(seqmat.shape[1])
    if norm == "max":
        # The bound takes the whole off-diagonal maximum, capped at an indel pair as
        # R caps it (fmin2(maxscost, 2*indel)); R's upper-triangle scan understates it
        # for asymmetric matrices (TX-AUD-20261019-008).
        max_sub = _py._subm_max_bound(subm)
        denom = float(seqmat.shape[1]) * min(max_sub, 2.0 * indel)
        if denom > 0:
            return D / denom
    return D


def _normalize_square_output(
    distances: np.ndarray,
    *,
    maxdist: float,
    length_1: float,
    length_2: float,
    code: int,
) -> np.ndarray:
    """Apply one pair-independent TraMineR bound to a square raw matrix.

    Each cell is normalized on its own: with asymmetric costs a reference
    distance fills only the (sequence, reference) cell, which can lie below the
    diagonal (TX-REP-20261017-001; batch-6 review C-1).
    """
    if code == 0:
        return distances
    out = np.asarray(distances, dtype=float).copy()
    for i in range(out.shape[0]):
        for j in range(i + 1, out.shape[1]):
            upper = normalize_value(out[i, j], maxdist, length_1, length_2, code)
            lower = (
                upper
                if out[j, i] == out[i, j]
                else normalize_value(out[j, i], maxdist, length_1, length_2, code)
            )
            out[i, j], out[j, i] = upper, lower
    return out


def om_trate_pairwise(
    seqmat: np.ndarray,
    *,
    indel: float | None = None,
    trate_smoothing: float = 0.0,
    trate_mode: str = "symmetric",
    trate_normalize: str = "none",
    trate_weights: np.ndarray | None = None,
    norm: Optional[str] = "length",
    deduplicate: bool = False,
    condensed: bool = False,
) -> np.ndarray:
    """Compute OM pairwise distances with a TRATE-derived substitution matrix.

    Parameters
    ----------
    seqmat : np.ndarray
        Integer-coded sequence matrix of shape ``(n, L)`` where:
        - ``n`` is the number of sequences
        - ``L`` is sequence length
        Recommended dtype is integer (for example ``int32``/``int64``).
    indel : float | None, default=None
        Scalar insertion/deletion cost. If ``None``, resolved exactly as
        ``seqdist(subm="TRATE", indel="auto")`` and R's ``seqdist.R:359-365``
        do it: ``max(sm) / 2`` over the whole substitution matrix. Until the
        cycle-10 follow-on (T-06) this helper instead used
        ``0.5 * mean(off-diagonal)``, so the package's two TRATE entry points
        returned different numbers for the same corpus; pass ``indel=``
        explicitly for any other rule.
    trate_smoothing : float, default=0.0 (R-faithful; no smoothing)
        Smoothing used when building the TRATE substitution matrix.
    trate_mode : str, default="symmetric"
        TRATE construction mode.
    trate_normalize : str, default="none"
        TRATE matrix normalization mode.
    trate_weights : np.ndarray | None, default=None
        Optional per-sequence weights of shape ``(n,)``.
    norm : None | "length" | "max", default="length"
        Optional distance normalization mode, in the low-level ``pairwise``
        vocabulary: ``"length"`` divides each distance by the longer sequence
        length, without the indel factor TraMineR's ``maxlength`` applies, so
        with an indel other than 1 it differs from
        ``seqdist(norm="maxlength")`` (docs/DISTANCES_NORMALIZATION.md).
    deduplicate : bool, default=False
        If True, computes on unique rows then expands back.
    condensed : bool, default=False
        If True, returns condensed upper-triangular vector.

    Returns
    -------
    np.ndarray
        - If ``condensed=False``: square float distance matrix of shape ``(n, n)``.
        - If ``condensed=True``: float vector of shape ``(n*(n-1)/2,)``.
    """
    sm = trate_substitution_matrix(
        seqmat,
        smoothing=float(trate_smoothing),
        mode=str(trate_mode),
        normalize=str(trate_normalize),
        weights=trate_weights,
    )
    if indel is None:
        # R's rule, shared with the seqdist TRATE route (seqdist.R:359-365).
        indel = 0.5 * float(_py._subm_max_all(sm))
    D = pairwise(
        seqmat,
        method="OM",
        indel=float(indel),
        subm=sm,
        norm=norm,
        deduplicate=deduplicate,
        condensed=False,  # compute full then optionally condense
    )
    return to_condensed(D) if condensed else D


# Order-based methods whose per-pair kernels admit variable-length sequences, so the
# structural ``void`` symbol of ragged sequences can be deleted before the distance.
_VOID_DELETE_METHODS = (
    "LCS",
    "OM",
    "LCP",
    "RLCP",
    "CHI2",
    "EUCLID",
    "OMSLEN",
    "OMSPELL",
    "TWED",
    "NMS",
    "NMSMST",
    "SVRSPELL",
    "OMSTRAN",
    "OMLOC",
)

# HAM and DHD are admitted on void-bearing input too, deliberately. Void deletion
# may leave every sequence the same length, in which case they succeed; when it
# does not, they reach the more specific "not defined for sequences of different
# length" error, which is more useful than the generic refusal below.
_EQUAL_LENGTH_VOID_METHODS = ("HAM", "DHD")

#: The complete set the void guard admits. The refusal message is built from this
#: same tuple so the two cannot state different sets (cycle-6 F-006).
_VOID_ADMITTED_METHODS = _VOID_DELETE_METHODS + _EQUAL_LENGTH_VOID_METHODS


def _validated_reference_indices(values, *, n_rows: int, one_based: bool, name: str) -> np.ndarray:
    """Validate an index vector without silently truncating fractional values."""
    raw = np.atleast_1d(np.asarray(values))
    if raw.size == 0 or raw.dtype.kind not in {"i", "u", "f"}:
        raise ValueError(f"{name} must contain positive integer indices")
    numeric = raw.astype(float)
    if not np.all(np.isfinite(numeric)) or not np.all(numeric == np.floor(numeric)):
        raise ValueError(f"{name} must contain positive integer indices")
    indices = numeric.astype(int)
    lower = 1 if one_based else 0
    if np.any(indices < lower) or np.any(indices >= n_rows + int(one_based)):
        raise ValueError(f"{name} contains an index outside the sequence data")
    return indices - 1 if one_based else indices


def _trate_profile_settings(kwargs: dict) -> tuple[float, str, str]:
    """TRATE smoothing, mode and normalisation for ``subm="TRATE"`` with a profile.

    ``trate_profile="calibrated"`` reads the calibration file through
    :func:`traminex.config._calibrated_profile_values`, which fails loudly when no
    file loads and lets an explicit ``trate_*`` argument override the file's
    value (TX-AUD-20261010-032); ``"default"`` uses the arguments alone.
    """
    profile = str(kwargs.get("trate_profile", "default")).lower()
    if profile not in {"default", "calibrated"}:
        raise ValueError("trate_profile must be 'default' or 'calibrated'")
    if profile == "default":
        refuse_unread_calibration_path(kwargs)
    explicit = {
        key: kwargs[f"trate_{key}"]
        for key in ("smoothing", "mode", "normalize")
        if f"trate_{key}" in kwargs
    }
    if profile == "calibrated":
        from ..config import _calibrated_profile_values

        explicit = _calibrated_profile_values(
            "trate",
            explicit,
            kwargs.get("calibrated_params_paths", kwargs.get("calibrated_params_path")),
        )
    return (
        float(explicit.get("smoothing", 0.0)),
        str(explicit.get("mode", "symmetric")),
        str(explicit.get("normalize", "none")),
    )


def _low_level_rows(rows: list[np.ndarray], method: str, options: dict, **row_kwargs) -> np.ndarray:
    """Row-engine call on behalf of the low-level ``pairwise`` vocabulary.

    The row engine speaks TraMineR's vocabulary, where ``'length'`` is
    ``maxlength``. The low-level API documents its own ``'length'``: the raw
    distance over ``l1 + l2`` for LCS/LCP/RLCP and over the longer length for OM
    -- what its dense path computes (docs/DISTANCES_NORMALIZATION.md). Without
    this, a void row anywhere in the input changed every normalized value
    (2026-09-22 audit S1b-07).
    """
    m = method.upper()
    norm = options.get("norm", None)
    if m in _CLASSIC_NORM_METHODS:
        norm = _map_norm(m, norm)
    legacy_length = (
        isinstance(norm, str)
        and norm.strip().lower() == "length"
        and m in {"OM", "LCS", "LCP", "RLCP"}
    )
    if not legacy_length:
        return _py._pairwise_rows(rows, method, options, **row_kwargs)
    raw_options = {key: value for key, value in options.items() if key != "norm"}
    raw = _py._pairwise_rows(rows, method, raw_options, **row_kwargs)
    lengths = np.array([len(row) for row in rows], dtype=float)
    denominator = (
        np.maximum.outer(lengths, lengths) if m == "OM" else np.add.outer(lengths, lengths)
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(denominator > 0, raw / denominator, 0.0)


def _pairwise_void_delete(seqmat: np.ndarray, void: int, method: str, **kwargs) -> np.ndarray:
    """R-faithful order-based distances on ragged sequences (TraMineR ``seqdist``).

    TraMineR deletes the structural ``void`` symbol (the right padding of ragged
    sequences) before computing an order-based distance, while keeping the missing
    state ``nr`` as a regular state (``with.missing=TRUE``). traminex's compiled
    kernels operate on the rectangular matrix and would treat ``void`` as an ordinary
    state, inflating common subsequences. When ``void`` cells are present this routes
    the per-sequence ``void``-deleted codes to the row engine
    (``python_kernels._pairwise_rows``), which evaluates pairs in vectorized batches
    where the method has a batched kernel.

    Only the order-based methods in :data:`_VOID_DELETE_METHODS` admit
    variable-length sequences. Normalization is computed per pair from the two
    canonical (void-deleted) row lengths.
    """
    m = method.upper()
    if m not in _VOID_ADMITTED_METHODS:
        raise NotImplementedError(
            f"{m!r} does not accept sequences containing the structural void "
            "symbol. Void-bearing input is supported for the order-based methods "
            f"{list(_VOID_DELETE_METHODS)}, which admit variable-length "
            f"sequences, and for {list(_EQUAL_LENGTH_VOID_METHODS)} when void "
            "deletion leaves every sequence the same length; drop the void "
            "positions or use one of those methods"
        )
    ref_sets = kwargs.pop("ref_sets", None)
    condensed = bool(kwargs.pop("condensed", False))
    kwargs.pop("deduplicate", None)
    kwargs.pop("_no_dedup", None)
    n_states = kwargs.pop("_n_states", None)
    missing_code = kwargs.pop("_missing_code", None)

    seqmat = np.asarray(seqmat)
    seqs = [row[row != void] for row in seqmat]

    def rows_distances() -> np.ndarray:
        return _low_level_rows(
            seqs,
            m,
            kwargs,
            n_states=n_states,
            missing_code=missing_code,
            time_width=int(seqmat.shape[1]),
        )

    if ref_sets is not None:
        idx_a, idx_b = ref_sets
        n_rows = len(seqs)
        idx_a = _validated_reference_indices(idx_a, n_rows=n_rows, one_based=False, name="ref_sets")
        idx_b = _validated_reference_indices(idx_b, n_rows=n_rows, one_based=False, name="ref_sets")
        # The pairs (a, b), a first: R's d(set-1 element, set-2 element); the square
        # matrix holds d(earlier row, later row) (TX-AUD-20261019-009).
        with _py._only_pairs(n_rows, idx_a, idx_b):
            return rows_distances()[np.ix_(idx_a, idx_b)]
    D = rows_distances()
    return to_condensed(D) if condensed else D


def _has_asymmetric_costs(subm: Any) -> bool:
    """Whether ``subm`` is a substitution-cost array whose direction shows."""
    if subm is None or isinstance(subm, str) or np.ndim(subm) < 2:
        return False
    return bool(_py._cost_array_is_asymmetric(np.asarray(subm, dtype=float)))


def _configuration_is_dataset_fitted(method: str, kwargs: dict) -> bool:
    """Return whether the method's parameters are estimated from the corpus.

    TraMineR fits substitution costs and state frequencies from the whole
    ``seqdata`` *before* selecting a reference set, so restricting the corpus to
    the two reference sets changes the answer for these configurations.

    Finer than :func:`_deduplication_changes_dataset_fit`, which reports ``DHD``
    unconditionally: DHD with a caller-supplied ``subm`` routes to
    ``_pairwise_dhd_explicit_subm`` and is pair-local.  That distinction is
    load-bearing -- ``weightedcluster._seqdist`` reaches this module's
    cross-block entry point with an explicit DHD ``subm``.
    """
    mu = method.upper()
    if mu in {"CHI2", "EUCLID", "OMSTRAN"}:
        return True
    if mu == "DHD":
        return kwargs.get("subm", None) is None
    subm = kwargs.get("subm", None)
    return mu == "OM" and isinstance(subm, str) and subm.lower() == "trate"


def _deduplication_changes_dataset_fit(method: str, kwargs: dict) -> bool:
    """Return whether duplicate removal would change fitted method metadata.

    Derived from :func:`_configuration_is_dataset_fitted` so the two cannot drift
    apart: a method added to that set is automatically covered here. The only
    difference is ``DHD``, which is treated as corpus-sensitive even with an
    explicit ``subm`` -- deduplication is merely an optimisation, so being
    conservative there costs nothing, whereas the reference-set route must
    distinguish the pair-local explicit-``subm`` DHD that
    ``weightedcluster._seqdist`` relies on.

    Verified equivalent to the previous explicit form over all 19 dispatch keys
    (plus lower-case and unknown names) x 7 ``subm`` shapes: 336 combinations, 0
    differences.
    """
    return method.upper() == "DHD" or _configuration_is_dataset_fitted(method, kwargs)


def _resolve_om_trate(seqmat: np.ndarray, kwargs: dict, void: int | None) -> dict:
    """OM's documented low-level ``subm="TRATE"`` as a cost matrix, for every route.

    The rates come from the rows without their void cells, so the void route and a
    per-state indel vector get the costs seqdist resolves for the same data
    (TX-AUD-20261019-024); the dense rows are counted as before.
    """
    subm = kwargs.get("subm", None)
    if isinstance(subm, str) and subm.lower() == "trate":
        smoothing, mode, normalize = _trate_profile_settings(kwargs)
        weights = kwargs.get("trate_weights", None)
        matrix = np.asarray(seqmat)
        if void is not None and np.any(matrix == void):
            rows = [row[row != void] for row in matrix]
            sm = _trate_substitution_matrix_rows(
                rows,
                n_states=max((int(row.max()) + 1 for row in rows if row.size), default=0),
                smoothing=smoothing,
                mode=mode,
                normalize=normalize,
                weights=weights,
            )
        else:
            sm = trate_substitution_matrix(
                matrix, smoothing=smoothing, mode=mode, normalize=normalize, weights=weights
            )
        kwargs = {k: v for k, v in kwargs.items() if k not in TRATE_OM_KWARGS}
        kwargs["subm"] = sm
    elif isinstance(subm, str):
        raise TypeError("OM subm accepts float, ndarray, or string 'TRATE'")
    elif any(k in kwargs for k in TRATE_OM_KWARGS):
        raise TypeError("TRATE-specific kwargs require subm='TRATE'")
    return kwargs


def pairwise(seqmat: np.ndarray, method: str = "OM", **kwargs) -> np.ndarray:
    """Compute pairwise sequence distances using the fastest available backend.

    Parameters
    ----------
    seqmat : np.ndarray
        Integer-coded sequence matrix of shape ``(n, L)`` with integer dtype.
        Each row is one sequence and each column one aligned position.
    method : str, default="OM"
        Distance method name (for example ``OM``, ``HAM``, ``LCS``, ``LCP``, ``RLCP``,
        ``DHD``, ``DUR``, ``JSD_TRANS``, ``CHI2``).
    **kwargs
        Method-specific options plus dispatcher options:
        - ``deduplicate`` (bool): compute on unique sequences then expand
        - ``ref_sets`` (two sequences of zero-based row positions): compute the
          rectangular A×B distances, each with the A row first as R's
          ``refseq=list(a, b)`` (``seqdist``'s ``refseq`` takes R's one-based
          indices instead)
        - ``condensed`` (bool): return condensed vector for square outputs
        - ``norm``: ``None | "length" | "max"`` (method-dependent)

    Returns
    -------
    np.ndarray
        Float distance output:
        - Square matrix ``(n, n)`` by default
        - Rectangular matrix ``(|A|, |B|)`` when ``ref_sets`` is used
        - Condensed vector ``(n*(n-1)/2,)`` when ``condensed=True``
    """
    # Ragged sequences: the structural ``void`` symbol is deleted before order-based
    # distances (TraMineR ``seqdist``). ``void`` is opt-in (passed by ``seqdist`` from a
    # state-sequence object's ``void_code``); the dense path below is byte-unchanged when
    # no void cells are present.
    seqmat = _py._validate_dense_sequence_matrix(seqmat)
    # One matrix conversion before routing (TX-AUD-20261019-003).
    if "subm" in kwargs:
        kwargs["subm"] = _py._as_cost_matrix(kwargs["subm"])
    if method.upper() in _CLASSIC_NORM_METHODS:
        mapped_norm = _map_norm(method, kwargs.get("norm", None))
    void = kwargs.pop("void", None)
    if method.upper() == "OM":
        # Before any route: the void and indel-vector routes branch off below
        # (TX-AUD-20261019-024).
        kwargs = _resolve_om_trate(seqmat, kwargs, void)
    if void is not None:
        _sm = np.asarray(seqmat)
        if np.any(_sm == void):
            return _pairwise_void_delete(_sm, int(void), method, **kwargs)
    if method.upper() == "OM" and np.ndim(kwargs.get("indel", 1.0)) == 1:
        # The dense kernels take one indel cost; a per-state vector goes through
        # the void route's row engine with no void code (no state is -1), as it
        # does when void is present and as seqdist does (TX-AUD-20261011-050).
        return _pairwise_void_delete(np.asarray(seqmat), -1, method, **kwargs)

    # Extract QoL kwargs (and avoid passing them to backends)
    deduplicate = bool(kwargs.pop("deduplicate", False))
    _no_dedup = bool(kwargs.pop("_no_dedup", False))
    _allow_zero_indel = bool(kwargs.pop("_allow_zero_indel", False))
    condensed = bool(kwargs.pop("condensed", False))
    ref_sets = kwargs.pop("ref_sets", None)

    # Dataset-fitted methods must be estimated from the whole corpus, exactly as
    # TraMineR does before it selects a reference set, so the cross-block equals the
    # corresponding block of the full-corpus matrix.  Computing the square matrix here
    # and slicing it is the same route ``seqdist`` takes for ``refseq``/``ref_sets``
    # and that ``_pairwise_void_delete`` takes for ragged input; it is correct by
    # construction rather than by threading each method's fitted parameters through a
    # narrower corpus.  Pair-local methods keep the efficient rectangular path below.
    if ref_sets is not None and _configuration_is_dataset_fitted(method, kwargs):
        idx_a, idx_b = ref_sets
        full_options = dict(kwargs)
        if _allow_zero_indel:
            full_options["_allow_zero_indel"] = True
        full = np.asarray(pairwise(seqmat, method=method, deduplicate=deduplicate, **full_options))
        n_rows = len(seqmat)
        return full[
            np.ix_(
                _validated_reference_indices(
                    idx_a, n_rows=n_rows, one_based=False, name="ref_sets"
                ),
                _validated_reference_indices(
                    idx_b, n_rows=n_rows, one_based=False, name="ref_sets"
                ),
            )
        ]

    validate_pairwise_kwargs(method, kwargs, dispatcher=True)
    # The classic OM/LCS/LCP/RLCP family maps norm here ("length"/"max") for the compiled
    # backends. CHI2/EUCLID use a within-method boolean norm; the OM-family spell/warp
    # methods (TWED/OMslen/OMspell/OMstran) carry the TraMineR norm vocabulary
    # (maxlength/gmean/maxdist/YujianBo/auto) straight through to their kernels.
    if method.upper() in ("OM", "LCS", "LCP", "RLCP"):
        if mapped_norm is None:
            kwargs.pop("norm", None)
        else:
            kwargs["norm"] = mapped_norm

    # Between-sets distances if requested
    if ref_sets is not None:
        idxA, idxB = ref_sets
        n_rows = len(seqmat)
        A = np.asarray(seqmat)[
            _validated_reference_indices(idxA, n_rows=n_rows, one_based=False, name="ref_sets")
        ]
        B = np.asarray(seqmat)[
            _validated_reference_indices(idxB, n_rows=n_rows, one_based=False, name="ref_sets")
        ]
        between_options = dict(kwargs)
        if _allow_zero_indel:
            between_options["_allow_zero_indel"] = True
        return pairwise_between_sets(
            A,
            B,
            method=method,
            deduplicate=deduplicate,
            **between_options,
        )

    # Deduplicate path (square only)
    if deduplicate and not _no_dedup and not _deduplication_changes_dataset_fit(method, kwargs):
        X = np.asarray(seqmat)
        unique_indices: list[int] = []
        lookup: dict[tuple[int, ...], int] = {}
        inv = np.empty(X.shape[0], dtype=int)
        for index, row in enumerate(X):
            key = tuple(int(state) for state in row)
            unique_index = lookup.get(key)
            if unique_index is None:
                unique_index = len(unique_indices)
                lookup[key] = unique_index
                unique_indices.append(index)
            inv[index] = unique_index
        uniq = X[unique_indices]
        if uniq.shape[0] < X.shape[0]:
            unique_options = dict(kwargs)
            if _allow_zero_indel:
                unique_options["_allow_zero_indel"] = True
            Duniq = pairwise(
                uniq,
                method=method,
                _no_dedup=True,
                condensed=False,
                **unique_options,
            )
            Dfull = Duniq[inv][:, inv]
            return to_condensed(Dfull) if condensed else Dfull
        # else fall through to normal path

    m = method.upper()
    if m == "OM":
        indel_value = kwargs.get("indel", 1.0)
        if np.isscalar(indel_value):
            resolved_indel = float(np.asarray(indel_value).item())
            if (
                np.isnan(resolved_indel)
                or resolved_indel < 0
                or (resolved_indel == 0 and not _allow_zero_indel)
            ):
                raise ValueError("OM indel cost must be strictly positive")
            kwargs["indel"] = resolved_indel
        fixed_subm = kwargs.get("subm", 2.0)
        if fixed_subm is None:
            fixed_subm = 2.0
        kwargs["subm"] = _py._validate_substitution_costs(
            fixed_subm,
            seqmat=seqmat,
            method="OM",
            allow_scalar=True,
        )
        if isinstance(kwargs["subm"], np.ndarray):
            _py._warn_asymmetric_costs(kwargs["subm"], stacklevel=3)
        indel_by_pos, subm_by_pos = _py._validate_om_position_costs(
            seqmat,
            indel_by_pos=kwargs.get("indel_by_pos", None),
            subm_by_pos=kwargs.get("subm_by_pos", None),
        )
        if indel_by_pos is not None:
            kwargs["indel_by_pos"] = indel_by_pos
        if subm_by_pos is not None:
            kwargs["subm_by_pos"] = subm_by_pos
            _py._warn_asymmetric_costs(subm_by_pos, stacklevel=3)
    elif m == "HAM" and kwargs.get("subm", None) is not None:
        kwargs["subm"] = _py._validate_substitution_costs(
            kwargs["subm"],
            seqmat=seqmat,
            method="HAM",
            allow_scalar=False,
        )
        _py._warn_asymmetric_costs(kwargs["subm"], stacklevel=3)
    if _allow_zero_indel:
        kwargs["_allow_zero_indel"] = True

    impl = _impl_name()

    # Fast backends. Unit-cost HAM uses the compiled integer-count kernel for its
    # raw and supported length-normalized forms. Explicit substitution costs and
    # unsupported normalization requests route through the R-faithful Python
    # kernel so their result or refusal remains identical across backends.
    if impl == "pybind11" and _pb is not None:
        if m == "HAM" and hasattr(_pb, "pairwise_ham_unit"):
            # Compiled fast path for the R-faithful unit-cost HAM (subm=None): the raw
            # per-position mismatch count, optionally length-normalised. Integer counts
            # → byte-identical (atol=0) to the pure-Python HAM for any n, L. Explicit-subm
            # HAM (float sum, numpy accumulation order) stays on the Python kernel.
            subm = kwargs.get("subm", None)
            norm = kwargs.get("norm", None)
            if subm is None and norm in (None, "none", "maxlength", "auto", "length", "max"):
                Dham = _pb.pairwise_ham_unit(np.asarray(seqmat, dtype=np.int32))
                if norm not in (None, "none"):
                    Lham = int(np.asarray(seqmat).shape[1])
                    denom = float(Lham)  # maxcost == 1.0 for unit-cost HAM
                    if denom > 0:
                        Dham = Dham / denom
                return to_condensed(Dham) if condensed else Dham
            # subm given or an unsupported norm → the Python kernel (identical result/error)
            D = _py.pairwise(seqmat, method=method, **kwargs)
            return to_condensed(D) if condensed else D
        if m == "OM":
            if kwargs.get("indel_by_pos") is not None or kwargs.get("subm_by_pos") is not None:
                D = _py.pairwise(seqmat, method=method, **kwargs)
                return to_condensed(D) if condensed else D
            indel = float(kwargs.get("indel", 1.0))
            subm = kwargs.get("subm", 2.0)
            norm = kwargs.get("norm", None)
            if isinstance(subm, np.ndarray) and _py._cost_array_is_asymmetric(subm):
                D = _py.pairwise(seqmat, method=method, **kwargs)
                return to_condensed(D) if condensed else D
            if np.isposinf(indel):
                D = _py.pairwise(seqmat, method=method, **kwargs)
                return to_condensed(D) if condensed else D
            if isinstance(subm, np.ndarray):
                D = _pb.pairwise_om_mat(seqmat, indel, np.asarray(subm, dtype=np.float64))
            else:
                D = _pb.pairwise_om(seqmat, indel, float(subm))
            D = _normalize_om_output(D, seqmat=seqmat, subm=subm, indel=indel, norm=norm)
            return to_condensed(D) if condensed else D
        if m == "LCS" and hasattr(_pb, "pairwise_lcs"):
            D = _pb.pairwise_lcs(seqmat)
            if kwargs.get("norm", None) == "length":
                L = seqmat.shape[1]
                if L > 0:
                    D = D / float(2 * L)
            return to_condensed(D) if condensed else D
        if m == "LCP" and hasattr(_pb, "pairwise_lcp"):
            D = _pb.pairwise_lcp(seqmat)
            if kwargs.get("norm", None) == "length":
                L = seqmat.shape[1]
                if L > 0:
                    D = D / float(2 * L)
            return to_condensed(D) if condensed else D
        if m == "RLCP":
            if hasattr(_pb, "pairwise_rlcp"):
                D = _pb.pairwise_rlcp(seqmat)
                if kwargs.get("norm", None) == "length":
                    L = seqmat.shape[1]
                    if L > 0:
                        D = D / float(2 * L)
            else:
                D = _py.pairwise(seqmat, method=method, **kwargs)
            return to_condensed(D) if condensed else D
        if m == "DHD":
            # Compiled kernel only implements the legacy marginal metric; the default
            # TraMineR-faithful DHD and the explicit-subm path go through Python (_py.pairwise).
            if (
                _py._is_legacy_dhd(kwargs)
                and kwargs.get("subm") is None
                and hasattr(_pb, "pairwise_dhd")
            ):
                freq, var_code, _variant = _py._resolve_legacy_dhd_settings(seqmat, kwargs)
                D = _pb.pairwise_dhd(seqmat, np.asarray(freq, dtype=np.float64), int(var_code))
            else:
                D = _py.pairwise(seqmat, method=method, **kwargs)
            return to_condensed(D) if condensed else D

    # Fallback to Python implementations (includes HAM/OM/LCS/LCP/RLCP/DHD)
    D = _py.pairwise(seqmat, method=method, **kwargs)
    return to_condensed(D) if condensed else D


def pairwise_between_sets(
    A,
    B,
    *,
    method: str = "OM",
    deduplicate: bool = False,
    **kwargs,
) -> np.ndarray:
    """Compute rectangular distances between two sequence sets.

    Parameters
    ----------
    A : array-like
        Integer-coded matrix of shape ``(nA, L)`` or an iterable of canonical
        one-dimensional rows with variable lengths.
    B : array-like
        Integer-coded matrix of shape ``(nB, L)`` or an iterable of canonical
        one-dimensional rows with variable lengths.
    method : str, default="OM"
        Distance method name passed to :func:`pairwise`.
    deduplicate : bool, default=False
        If True, enable deduplication in the internal computation.
    **kwargs
        Additional arguments passed through to :func:`pairwise`.

    Returns
    -------
    np.ndarray
        Float distance matrix of shape ``(nA, nB)``.

    Raises
    ------
    ValueError
        If ``method`` and ``kwargs`` select a configuration whose parameters are
        estimated from the corpus (``CHI2``, ``EUCLID``, ``OMSTRAN``, ``DHD``
        without an explicit ``subm``, or ``OM`` with ``subm="TRATE"``).  Given
        only ``A`` and ``B`` this function cannot know the corpus those
        parameters should be fitted from, and TraMineR fits them from the whole
        ``seqdata`` before selecting a reference set.

    Notes
    -----
    Equal-width dense matrices use the compiled stacked-data fast path.  When
    either set contains variable-length rows, the structural ``void`` symbol is
    removed (when supplied) and the R-facing row engine computes the
    cross-block.  Cell ``(a, b)`` is the distance with ``A``'s row first, as
    TraMineR's ``refseq=list(a, b)`` gives d(set-1 element, set-2 element); with
    asymmetric costs this holds also when ``B`` repeats a row of ``A``
    (TX-AUD-20261019-009).
    """
    options = dict(kwargs)
    options.pop("condensed", None)
    options.pop("ref_sets", None)
    # One matrix conversion before routing: the ragged branch below calls the row
    # engine directly (TX-AUD-20261019-003; batch-7 review I-3).
    if "subm" in options:
        options["subm"] = _py._as_cost_matrix(options["subm"])

    if _configuration_is_dataset_fitted(method, options):
        raise ValueError(
            f"{method.upper()} parameters are estimated from the whole corpus, so a "
            "standalone two-set call cannot reproduce the full-corpus cross-block. "
            "Pass the whole corpus and select the two sets instead: "
            "pairwise(corpus, method=..., ref_sets=(a, b)), or use "
            "seqdist(x, ..., refseq=[a, b]) for the R-facing route."
        )

    def dense_matrix(value):
        try:
            array = np.asarray(value)
        except ValueError:
            return None
        return array if array.ndim == 2 else None

    dense_a = dense_matrix(A)
    dense_b = dense_matrix(B)
    if dense_a is not None and dense_b is not None and dense_a.shape[1] == dense_b.shape[1]:
        # Validate original dtypes before promotion (bool/int and uint64/int64
        # otherwise change their acceptance). Narrow only after range checks.
        _py._validate_dense_sequence_matrix(dense_a)
        _py._validate_dense_sequence_matrix(dense_b)
        combined = np.vstack([dense_a.astype(np.int32), dense_b.astype(np.int32)])
        nA = dense_a.shape[0]
        # Only the pairs (a, b), a from A first, are evaluated: with asymmetric costs that
        # is R's d(set-1 element, set-2 element), also when B repeats a row of A, and the
        # deduplicated square matrix would take the first occurrence's direction instead
        # (TX-AUD-20261019-009).
        with _py._only_pairs(combined.shape[0], range(nA), range(nA, combined.shape[0])):
            Dfull = pairwise(
                combined,
                method=method,
                deduplicate=deduplicate and not _has_asymmetric_costs(options.get("subm")),
                condensed=False,
                **options,
            )
        return Dfull[:nA, nA:]

    void = options.pop("void", None)
    n_states = options.pop("_n_states", None)
    missing_code = options.pop("_missing_code", None)

    def canonical_rows(value, name: str) -> list[np.ndarray]:
        dense = dense_matrix(value)
        if dense is not None:
            candidates = [dense[i] for i in range(dense.shape[0])]
        else:
            try:
                candidates = list(value)
            except TypeError as exc:
                raise ValueError(f"{name} must be a 2D matrix or iterable of rows") from exc
        rows: list[np.ndarray] = []
        for row in candidates:
            array = np.asarray(row)
            if array.ndim != 1:
                raise ValueError(f"{name} must contain one-dimensional sequence rows")
            _py._validate_dense_sequence_matrix(array[None, :])
            if void is not None:
                array = array[array != void]
            rows.append(array)
        return rows

    rows_a = canonical_rows(A, "A")
    rows_b = canonical_rows(B, "B")
    rows = rows_a + rows_b

    # Dataset-derived methods must see every original row (or aggregated
    # weights), so keep their preprocessing untouched.  For pair-local
    # methods, deduplicate canonical rows exactly as the dense dispatcher does
    # and expand the cross-block afterward.
    if (
        deduplicate
        and not _deduplication_changes_dataset_fit(method, options)
        and not _has_asymmetric_costs(options.get("subm"))
    ):
        unique_rows: list[np.ndarray] = []
        lookup: dict[tuple[int, ...], int] = {}
        inverse = np.empty(len(rows), dtype=int)
        for index, row in enumerate(rows):
            key = tuple(int(state) for state in row)
            unique_index = lookup.get(key)
            if unique_index is None:
                unique_index = len(unique_rows)
                lookup[key] = unique_index
                unique_rows.append(row)
            inverse[index] = unique_index
        if len(unique_rows) < len(rows):
            unique_distances = _low_level_rows(
                unique_rows,
                method,
                options,
                n_states=n_states,
                missing_code=missing_code,
            )
            Dfull = unique_distances[inverse][:, inverse]
            nA = len(rows_a)
            return Dfull[:nA, nA:]

    nA = len(rows_a)
    with _py._only_pairs(len(rows), range(nA), range(nA, len(rows))):
        Dfull = _low_level_rows(
            rows,
            method,
            options,
            n_states=n_states,
            missing_code=missing_code,
        )
    return Dfull[:nA, nA:]


def _require_transition_rate_data(n_states: int, width: int) -> None:
    """Refuse the data seqcost(TRATE) and seqtrate refuse, in their order.

    Transition-rate costs need at least two states and at least one transition
    (a sequence width above the lag of 1); R stops on both too
    (TX-AUD-20261011-039, -040).
    """
    if n_states < 2:
        raise ValueError("TRATE costs require at least two states in this mode")
    if width <= 1:
        raise ValueError("absolute lag must be smaller than the sequence width")


def _cost_state_labels(source, size: int) -> list[str] | None:
    """The state names R's checkcost prints for cost rows ``0..size-1``.

    A cost matrix is indexed by state code, so each row is named by decoding
    its code (``nr`` included); a plain array has no names.
    """
    if not isinstance(source, STSList):
        return None
    try:
        return [str(source._state_for_code(code)) for code in range(size)]
    except IndexError:
        return None


def _distance_codebook(x: STSList) -> STSList:
    """Put real states before nr/void so public costs address logical ranks.

    Manual STSList objects may place declared markers anywhere in the codebook.
    Public cost matrices follow the real alphabet, then missing, with void
    excluded. Ordinary factory objects already have this layout.
    """
    x = x._canonicalized()
    order = _real_state_codes(x) + [code for code in (x.nr_code, x.void_code) if code is not None]
    if order == list(range(len(x.alphabet))):
        return x
    inverse = np.empty(len(order), dtype=np.int64)
    inverse[order] = np.arange(len(order))
    result = x.copy()
    result.data = inverse[np.asarray(x.data, dtype=np.int64)]
    result.alphabet = [x.alphabet[code] for code in order]
    if x.labels is not None:
        result.labels = [x.labels[code] for code in order]
    if x.cpal is not None and len(x.cpal) >= len(x.alphabet):
        result.cpal = [x.cpal[code] for code in order]
    for attribute in ("nr_code", "void_code"):
        code = getattr(x, attribute)
        if code is not None:
            setattr(result, attribute, int(inverse[code]))
    return result


def seqdist(x, method: str = "OM", **kwargs) -> np.ndarray:
    """TraMineR-style distance wrapper for STS objects or raw integer matrices.

    Parameters
    ----------
    x : STSList | np.ndarray
        Either:
        - an STS-like object exposing ``.data`` of shape ``(n, L)``, or
        - a raw integer-coded matrix of shape ``(n, L)``.
    method : str, default="OM"
        Distance method name. ``method="HAM"`` returns the **raw** per-position
        substitution-cost sum by default (R-faithful, ``norm="none"``; unit cost per
        mismatch, or ``sum_t sm[a,b]`` with an explicit ``subm``); pass
        ``norm="maxlength"`` (== ``"auto"`` for fixed-length data) for the ``/L``
        proportion.
    **kwargs
        R-compatible distance options plus the lower-level :func:`pairwise`
        options. Common options include:
        - ``with_missing`` (bool): treat the missing state ``nr`` as an ordinary
          state, as R's ``seqdist(..., with.missing=TRUE)`` does. **Required**
          when the sequences contain missing values: otherwise ``seqdist``
          raises ``ValueError: with_missing must be true when sequences contain
          missing values``. Conversely, when neither the sequences nor an
          external ``refseq`` contains a missing value, it is reset to
          ``False`` (``seqdist.R:111-116``), for a plain integer matrix as
          well as a sequence object. It is accepted through ``**kwargs`` and
          so does not appear in the signature or the stub.
        - ``indel``, ``subm``, ``norm``; R-facing methods that consume
          substitution costs accept ``subm="TRATE"``, ``"CONSTANT"``,
          ``"INDELS"``, or ``"INDELSLOG"`` with method-specific behavior.
          ``norm`` accepts the full TraMineR vocabulary for the classic
          OM/LCS/LCP/RLCP family — ``None``, ``"none"``, ``"length"``,
          ``"max"``, ``"maxlength"``, ``"gmean"``, ``"maxdist"``,
          ``"YujianBo"``, ``"auto"`` — case-insensitively. The traminex-only
          aliases ``"length"`` (== ``"maxlength"``) and ``"max"``
          (== ``"maxdist"``) are accepted by every method that takes a
          TraMineR ``norm``, not only the classic four
        - ``deduplicate``, ``condensed``, ``full_matrix``
        - ``refseq=k`` for distances to one 1-based row (``0`` selects the
          first most-frequent sequence)
        - ``refseq=(rows, cols)`` for two positive 1-based index sets, or a
          one-row STS object for an external reference
        - ``ref_sets=(rows, cols)`` for the existing zero-based Python
          rectangular interface; it cannot be combined with ``refseq``
        - method-specific knobs (for example ``metric`` for ``DUR``)

    Returns
    -------
    np.ndarray
        Float distance output:
        - square ``(n, n)`` by default
        - a length-``n`` vector for one internal/external reference
        - rectangular ``(|A|, |B|)`` for two reference sets
        - condensed vector when ``condensed=True``
    """
    # seqdist.R:18-37 maps the deprecated methods and a logical norm (audit S1a-12).
    if isinstance(method, str) and method.upper() in {"OMOPT", "LCSOPT"}:
        replacement = "OM" if method.upper() == "OMOPT" else "LCS"
        warnings.warn(
            f"method {method!r} is deprecated in TraMineR; using the equivalent {replacement!r}",
            DeprecationWarning,
            stacklevel=_caller_stacklevel(),
        )
        method = replacement
    if isinstance(kwargs.get("norm", None), (bool, np.bool_)):
        mapped_norm = "auto" if kwargs["norm"] else "none"
        warnings.warn(
            f"norm={bool(kwargs['norm'])} is deprecated in TraMineR; using {mapped_norm!r}",
            DeprecationWarning,
            stacklevel=_caller_stacklevel(),
        )
        kwargs["norm"] = mapped_norm
    if isinstance(kwargs.get("refseq", None), (bool, np.bool_)):
        # R: is.a.positive.integer(TRUE) is FALSE, so refseq = TRUE is refused.
        raise ValueError(
            "refseq must be a non-negative integer index, two index sets, or one "
            "sequence object; a boolean is not an index (TraMineR: 'refseq' invalid type)"
        )
    with_missing = bool(kwargs.pop("with_missing", False))
    weighted = bool(kwargs.pop("weighted", True))
    # Set only by the index and index-set reference routes, whose second call on
    # the same data must not print again the input diagnostics the first call
    # printed: R prints each once per call (batch-6 review I-1).
    diagnostics_shown = bool(kwargs.pop("_diagnostics_shown", False))
    # Both marker representations are equivalent inputs (canonicalization is
    # the identity for in-alphabet objects). The kernels bound the maximum
    # code VALUE while the state-count bound below counts states, so a
    # declared external nr_code/void_code (a code outside
    # range(len(alphabet))) must be folded into the alphabet before codes are
    # extracted.
    if isinstance(x, STSList):
        original = x
        x = _distance_codebook(x)
        explicit_void = kwargs.get("void", None)
        if explicit_void is not None and explicit_void in original._codebook_codes():
            kwargs["void"] = x._code_for_state(original._state_for_code(int(explicit_void)))
    if isinstance(kwargs.get("refseq"), STSList):
        kwargs["refseq"] = _distance_codebook(kwargs["refseq"])
    # A matrix given as nested lists is read as the array it spells, on every route
    # (TX-AUD-20261019-003). A logical or complex sm is left for R's ordered checks,
    # which refuse it at the sm step, after norm and indel (batch-9 review M-1).
    if "subm" in kwargs and not is_refused_cost_scalar(kwargs["subm"]):
        kwargs["subm"] = _py._as_cost_matrix(kwargs["subm"])
    pending_refseq = kwargs.get("refseq", None)
    full_matrix = kwargs.pop("full_matrix", None)
    if full_matrix is not None:
        requested_condensed = not bool(full_matrix)
        if "condensed" in kwargs and bool(kwargs["condensed"]) != requested_condensed:
            raise ValueError("full_matrix and condensed request conflicting output forms")
        kwargs["condensed"] = requested_condensed

    nr_code = None
    pending_ref_nr_code = None
    active_state_count = None
    source_weights = None

    def is_data_wrapper(value) -> bool:
        # ndarray.data is a memoryview, whereas the historical generic STS-like
        # contract accepts any non-ndarray wrapper exposing a matrix-like .data.
        # A NumPy scalar has a memoryview .data too and is a number, such as an
        # index refseq=np.int64(2) (batch-6 review M-4).
        return not isinstance(value, (np.ndarray, np.generic)) and hasattr(value, "data")

    source_obj = x if is_data_wrapper(x) else None
    if source_obj is not None:
        seqmat = source_obj.data
        source_weights = getattr(x, "weights", None)
        nr_code = getattr(x, "nr_code", None)
        # A weighted object's own weights (selected by ``weighted=``) replace a
        # sequence-weight keyword, which was therefore dropped silently (audit
        # S1a-06). WHAM's ``weights=`` are position weights and stay allowed.
        sequence_weight_keywords = ["trate_weights"]
        if str(method).upper() in {"CHI2", "OMSTRAN"}:
            sequence_weight_keywords.append("weights")
        for keyword in sequence_weight_keywords:
            if source_weights is not None and kwargs.get(keyword) is not None:
                raise ValueError(
                    f"{keyword}= conflicts with the weights the sequence object "
                    "carries; set them on the object and choose them with weighted="
                )
    else:
        seqmat = np.asarray(x)
    # seqdist.R:101-109: a sequence with nothing to compare -- no element other
    # than void, and missing unless with.missing -- gets R's message, and OMloc
    # stops; R checks this before its with.missing rules (TX-AUD-20261011-042).
    if source_obj is not None:
        not_compared = {getattr(x, "void_code", None)}
        if not with_missing:
            not_compared.add(nr_code)
        not_compared.discard(None)
        empty_rows = [
            str(index)
            for index, row in enumerate(np.asarray(seqmat), start=1)
            if all(int(code) in not_compared for code in row)
        ]
        if empty_rows:
            if str(method).upper() == "OMLOC":
                raise ValueError(f"OMloc does not support empty sequences {', '.join(empty_rows)}")
            if not diagnostics_shown:
                warnings.warn(
                    f" [!!] Empty sequence(s)  {', '.join(empty_rows)}  may cause inconsistent results!",
                    UserWarning,
                    stacklevel=2,
                )
    # A plain integer array carries no missing-state code, so only an external
    # reference can bring nr into it; the rule is otherwise identical.
    has_missing = nr_code is not None and bool(np.any(np.asarray(seqmat) == nr_code))
    if is_data_wrapper(pending_refseq):
        ref_nr_code = getattr(pending_refseq, "nr_code", None)
        pending_ref_nr_code = ref_nr_code
        has_missing = has_missing or (
            ref_nr_code is not None
            and bool(np.any(np.asarray(getattr(pending_refseq, "data")) == ref_nr_code))
        )
    if has_missing and not with_missing:
        raise ValueError("with_missing must be true when sequences contain missing values")
    # TraMineR resets with.missing=FALSE when neither the source nor an
    # external reference actually materializes nr (seqdist.R:111-116), for
    # either input type.  A merely declared missing symbol therefore does not
    # widen DHD's active cost dimensions, and a plain array never widens the
    # costs seqcost generates for it.
    if not has_missing:
        if with_missing:
            # R says so (seqdist.R:114-117; TX-AUD-20261011-042).
            warnings.warn(
                " [!!] seqdist: 'with.missing' set as FALSE as 'seqdata' has no "
                "non-void missing values",
                UserWarning,
                stacklevel=2,
            )
        with_missing = False
    if source_obj is not None:
        # Ragged state-sequence objects carry a structural ``void`` code; pass it through
        # so the distance layer deletes void (TraMineR ``seqdist``) instead of treating it
        # as a state. Inert when no void cells are present.
        void_code = getattr(x, "void_code", None)
        if void_code is not None and "void" not in kwargs:
            kwargs["void"] = void_code
        if hasattr(x, "alphabet"):
            real_state_count = len(_real_state_codes(x))
            active_state_count = real_state_count + int(with_missing and nr_code is not None)
    seqmat = _py._validate_dense_sequence_matrix(seqmat)

    # ``pairwise`` predates the R-facing wrapper and retains two legacy aliases:
    # ``norm='length'`` and ``norm='max'``.  Consume R's explicit no-normalization
    # spelling here without changing those lower-level contracts.
    m = method.upper()
    # R's seqdist takes every one of its arguments for every method: run R's
    # checks and drop what the method does not use (TX-AUD-20261010-071).
    resolve_r_seqdist_arguments(m, kwargs, n_states=active_state_count)
    if source_weights is not None and m in {"CHI2", "OMSTRAN"}:
        kwargs["weights"] = (
            np.asarray(source_weights, dtype=float)
            if weighted
            else np.ones(len(seqmat), dtype=float)
        )

    # Resolve R's dataset-derived OM costs before either the ragged adapter or
    # an external reference appends a zero-weight row.  Transition rates are
    # fitted on the source sequences only; an observed-missing external
    # reference can nevertheless widen the active alphabet.
    om_subm = kwargs.get("subm", None)
    generated_method_subm = False
    if m == "OM" and isinstance(om_subm, str) and om_subm.lower() == "trate":
        generated_method_subm = True
        smoothing, trate_mode, trate_normalize = _trate_profile_settings(kwargs)

        dense_source = np.asarray(seqmat, dtype=np.int32)
        source_void = kwargs.get("void", None)
        trate_rows = [
            row[row != source_void] if source_void is not None else row for row in dense_source
        ]
        observed_states = max((int(row.max()) + 1 for row in trate_rows if len(row)), default=0)
        trate_state_count = max(
            observed_states,
            0 if active_state_count is None else int(active_state_count),
            0 if pending_ref_nr_code is None or not with_missing else int(pending_ref_nr_code) + 1,
        )
        _require_transition_rate_data(trate_state_count, int(dense_source.shape[1]))
        explicit_trate_weights = kwargs.get("trate_weights", None)
        trate_weights: np.ndarray | None
        if source_weights is not None:
            trate_weights = (
                _validate_sequence_weights(source_weights, len(trate_rows))
                if weighted
                else np.ones(len(trate_rows), dtype=float)
            )
        else:
            trate_weights = (
                None
                if explicit_trate_weights is None
                else _validate_sequence_weights(
                    explicit_trate_weights,
                    len(trate_rows),
                    name="trate_weights",
                )
            )
        resolved_subm = _trate_substitution_matrix_rows(
            trate_rows,
            n_states=trate_state_count,
            smoothing=smoothing,
            mode=trate_mode,
            normalize=trate_normalize,
            weights=trate_weights,
        )
        # seqcost.R derives its automatic indel before applying the fixed
        # missing-state cost.  A dataset with transitions only through missing
        # can therefore retain a zero indel even though the final missing row
        # and column are overwritten with cost 2.
        trate_auto_indel = 0.5 * _py._subm_max_all(resolved_subm)
        effective_missing_code = nr_code if nr_code is not None else pending_ref_nr_code
        if with_missing and effective_missing_code is not None:
            missing_index = int(effective_missing_code)
            resolved_subm[missing_index, :] = 2.0
            resolved_subm[:, missing_index] = 2.0
            resolved_subm[missing_index, missing_index] = 0.0

        indel_arg = kwargs.get("indel", None)
        if indel_arg is None or (isinstance(indel_arg, str) and indel_arg.lower() == "auto"):
            kwargs["indel"] = trate_auto_indel
            if trate_auto_indel == 0:
                kwargs["_allow_zero_indel"] = True
        kwargs = {key: value for key, value in kwargs.items() if key not in TRATE_OM_KWARGS}
        # R resolves sm = "TRATE" through seqcost, which checks the generated matrix
        # once, with seqcost's own indel (seqcost.R:286; TX-AUD-20261011-080).
        _py._warn_cost_quality(
            resolved_subm,
            indel=trate_auto_indel,
            labels=_cost_state_labels(source_obj, int(np.shape(resolved_subm)[-1])),
            stacklevel=3,
        )
        kwargs["subm"] = resolved_subm

    # TraMineR accepts four source-active substitution-cost method names and
    # resolves them through seqcost before method dispatch.  Keep OM/TRATE's
    # established profile-aware extension above; resolve every other faithful
    # method string here.  LCS validates the name and then deliberately ignores
    # the selected cost method because its costs are fixed by definition.
    cost_method_arg = kwargs.get("subm", None)
    string_cost_methods = {"TRATE", "CONSTANT", "INDELS", "INDELSLOG"}
    string_cost_consumers = {
        "OM",
        "OMLOC",
        "OMSLEN",
        "OMSPELL",
        "OMSTRAN",
        "HAM",
        "LCS",
        "TWED",
    }
    if isinstance(cost_method_arg, str) and m in string_cost_consumers:
        resolved_cost_method = cost_method_arg.upper()
        if resolved_cost_method not in string_cost_methods:
            raise ValueError("subm must be one of: TRATE, CONSTANT, INDELS, INDELSLOG")
        if m == "LCS":
            kwargs.pop("subm")
        else:
            if isinstance(source_obj, STSList):
                cost_source = source_obj
            else:
                state_count = int(seqmat.max()) + 1 if seqmat.size else 0
                cost_source = STSList(
                    data=np.asarray(seqmat, dtype=int).copy(),
                    alphabet=[str(code) for code in range(state_count)],
                    weights=(
                        None
                        if source_weights is None
                        else np.asarray(source_weights, dtype=float).copy()
                    ),
                )
            cval = 1.0 if m == "HAM" and resolved_cost_method == "CONSTANT" else 2.0
            generated_costs = seqcost(
                cost_source,
                resolved_cost_method,
                cval=cval,
                with_missing=with_missing,
                miss_cost=cval,
                weighted=weighted,
            )
            kwargs["subm"] = np.asarray(generated_costs["sm"], dtype=float)
            if m == "OMSTRAN" and not np.all(np.isfinite(kwargs["subm"])):
                # Batch-9 review I-1: the substitution costs are at fault, not the indel.
                raise ValueError(
                    f"OMstran needs finite substitution costs, but the {resolved_cost_method} "
                    "costs generated for these data are not finite"
                    + (
                        ": each state the data never show gets the largest double as its "
                        "indel, and the substitution cost between two such states overflows"
                        if resolved_cost_method == "INDELS"
                        else ""
                    )
                )
            generated_method_subm = True
            indel_arg = kwargs.get("indel", None)
            auto_indel = indel_arg is None or (
                isinstance(indel_arg, str) and indel_arg.lower() == "auto"
            )
            if resolved_cost_method == "INDELS":
                # seqcost has checked R's costs; the stand-in of a state no compared
                # sequence contains is kept out of the scale (TX-AUD-20261019-002).
                canonical_source = cost_source._canonicalized()
                present_codes = np.unique(canonical_source.data)
                if canonical_source.void_code is not None:
                    present_codes = present_codes[present_codes != canonical_source.void_code]
                if is_data_wrapper(pending_refseq):
                    # An external reference is compared too: its states are present.
                    reference_codes = np.asarray(getattr(pending_refseq, "data")).ravel()
                    reference_void = getattr(pending_refseq, "void_code", None)
                    if reference_void is not None:
                        reference_codes = reference_codes[reference_codes != reference_void]
                    present_codes = np.union1d(present_codes, reference_codes)
                kwargs["subm"], generated_indel = _py._without_unobserved_stand_ins(
                    kwargs["subm"], np.asarray(generated_costs["indel"], dtype=float), present_codes
                )
                generated_costs = {**generated_costs, "indel": generated_indel}
            if auto_indel and m not in {"OMLOC", "HAM", "TWED"}:
                kwargs["indel"] = generated_costs["indel"]

    # R generates LCS's costs, and HAM's when sm is not given, with
    # seqsubm(seqdata, "CONSTANT", ...) (seqdist.R:370-374, 452-455), and that
    # seqcost checks them once, with its indel of 1 (seqcost.R:286;
    # TX-AUD-20261011-080): a one-state alphabet gives "max <= 0".
    if (
        not diagnostics_shown
        and active_state_count is not None
        and (m == "LCS" or (m == "HAM" and kwargs.get("subm", None) is None))
    ):
        constant_costs = np.full(
            (active_state_count, active_state_count), 2.0 if m == "LCS" else 1.0, dtype=float
        )
        np.fill_diagonal(constant_costs, 0.0)
        _py._warn_cost_quality(
            constant_costs,
            indel=1.0,
            labels=_cost_state_labels(source_obj, active_state_count),
            stacklevel=3,
        )

    # seqdist.R resolves a generic automatic indel against an explicit
    # substitution matrix before dispatch.  TWED uses a distinct formula;
    # OM+TRATE was resolved above before its missing-cost overwrite.
    if m in {"OM", "OMSLEN", "OMSPELL", "OMSTRAN"}:
        explicit_subm = kwargs.get("subm", None)
        if (
            np.isscalar(explicit_subm)
            and not isinstance(explicit_subm, str)
            and float(np.asarray(explicit_subm, dtype=float).item()) < 0
        ):
            raise ValueError(f"{m} substitution costs cannot be negative")
        indel_arg = kwargs.get("indel", None)
        auto_indel = indel_arg is None or (
            isinstance(indel_arg, str) and indel_arg.lower() == "auto"
        )
        if explicit_subm is not None and not isinstance(explicit_subm, str) and auto_indel:
            # ``seqdist.R:359-365``: ``indel <- max(sm)/2`` over the WHOLE matrix.
            # That scope is observable and must not be narrowed to the observed
            # sub-matrix: in live R, changing only the cost between two unobserved
            # states changes the normalized distances.
            automatic_indel = 0.5 * _py._subm_max_all(explicit_subm)
            if not np.isfinite(automatic_indel):
                # An ``Inf`` anywhere in the matrix -- even in a cell no observed
                # state can reach -- would make every distance ``inf``. R refuses
                # this input too (``checkcost`` computes ``tol <- 1e-07 * max(sma)``
                # and then evaluates ``NaN > Inf``, giving "missing value where
                # TRUE/FALSE needed"). traminex need not reproduce R's internal
                # error, but must not return a silently meaningless matrix.
                raise ValueError(
                    f"{m} substitution costs must be finite to derive an automatic "
                    "indel: max(subm)/2 is not finite, so every distance would be "
                    "infinite. Pass an explicit finite indel, or remove the "
                    "non-finite substitution cost."
                )
            kwargs["indel"] = automatic_indel

    # TWED's automatic indel is its own formula, 2 * max(sm) + nu + h over the
    # whole matrix, for an explicit matrix and a generated cost method alike
    # (seqdist.R:359-364, 433-439). The kernel resolves it from ``indel=None``
    # (``_twed_indel_value``); the string itself is the pairwise layer's to
    # refuse. checkcost's indel diagnostics cannot fire for that value (it
    # always exceeds max(sm) / 2), so handing them None loses nothing.
    if m == "TWED":
        twed_indel = kwargs.get("indel", None)
        if isinstance(twed_indel, str) and twed_indel.lower() == "auto":
            kwargs["indel"] = None

    explicit_subm_value = kwargs.get("subm", None)
    if (
        m in {"OM", "OMLOC", "OMSLEN", "OMSPELL", "OMSTRAN", "HAM", "TWED"}
        and explicit_subm_value is not None
        and not isinstance(explicit_subm_value, str)
        and not generated_method_subm
    ):
        substitution_costs = np.asarray(explicit_subm_value, dtype=float)
        ignored_codes = (kwargs.get("void", None),)
        if not _py._observed_substitution_costs_are_finite(
            substitution_costs,
            np.asarray(seqmat),
            ignored_codes=ignored_codes,
        ):
            raise ValueError(f"{m} substitution costs must be finite")
        if (
            m != "HAM"
            and substitution_costs.ndim == 2
            and source_obj is not None
            and hasattr(source_obj, "alphabet")
        ):
            # checkcost (TraMineR-checkcost.R): one row and column per state of
            # the alphabet, plus nr under with.missing -- exactly (audit S1b-01).
            # A plain matrix declares no alphabet; the kernels require it to cover
            # the observed codes. R's HAM path never reaches checkcost.
            expected = len(_real_state_codes(source_obj)) + int(with_missing)
            if substitution_costs.shape != (expected, expected):
                raise ValueError(
                    f"{m} subm must be {expected}x{expected}: one row and column per "
                    "state of the alphabet" + (", plus the missing state" if with_missing else "")
                )
        if m != "HAM":
            # checkcost stops on these before its triangle-inequality diagnostic
            # runs, so refuse here rather than first in the kernels.
            void_code = kwargs.get("void", None)
            _py._refuse_invalid_fixed_costs(
                m,
                substitution_costs,
                [row[row != void_code] if void_code is not None else row for row in seqmat],
            )

    indel_value = kwargs.get("indel", None)
    indel_subm = kwargs.get("subm", None)
    explicit_indel_contract = (
        indel_subm is not None and not isinstance(indel_subm, str) and not generated_method_subm
    )
    if (
        m in {"OM", "OMSLEN", "OMSPELL", "OMSTRAN", "TWED"}
        and indel_value is not None
        and explicit_indel_contract
    ):
        if not isinstance(indel_value, str):
            indel_values = np.atleast_1d(np.asarray(indel_value, dtype=float))
            invalid = bool(np.any(np.isnan(indel_values)))
            if m == "OMSTRAN" and np.any(np.isposinf(indel_values)):
                raise ValueError("OMSTRAN indel cost must be finite and strictly positive")
            if not invalid:
                invalid = (
                    float(indel_values.max()) <= 0
                    if m == "TWED" and indel_values.size
                    else bool(np.any(indel_values <= 0))
                )
            if invalid:
                raise ValueError(f"{m} indel cost must be strictly positive")

    # R ``checkcost``'s non-fatal cost-quality diagnostics
    # (SOURCES/TraMineR/R/TraMineR-checkcost.R, invoked from seqdist.R:380-390)
    # for a *caller-supplied* cost matrix, ``sm.type == "matrix"``, or a DHD
    # array, ``sm.type == "array"``. ``indel`` is passed for the OM family and
    # TWED (seqdist.R:381). TraMineR 2.2-14 also runs checkcost inside seqcost
    # (seqcost.R:286) on the matrices generated from a method string; traminex
    # does so in seqcost and in its own OM/TRATE and DHD routes
    # (TX-AUD-20261011-080), so a generated matrix (``generated_method_subm``)
    # is not checked again here. For HAM R runs no explicit-matrix check at all: in ``if (method %in%
    # c(om.methods, "TWED")) if (check.costs) checkcost(...) else if (method ==
    # "HAM") ...`` the ``else`` binds to the inner ``if``, so the HAM branch is
    # never reached. traminex deliberately still refuses an invalid HAM matrix
    # and gives it the substitution triangle-inequality diagnostic
    # (docs/API_PARITY.md).
    if (
        m in {"OM", "OMLOC", "OMSLEN", "OMSPELL", "OMSTRAN", "TWED", "HAM"}
        and explicit_subm_value is not None
        and not isinstance(explicit_subm_value, str)
        and not generated_method_subm
        and not diagnostics_shown
    ):
        cost_matrix = np.asarray(explicit_subm_value, dtype=float)
        _py._warn_cost_quality(
            cost_matrix,
            indel=None if m == "HAM" else kwargs.get("indel", None),
            labels=_cost_state_labels(source_obj, cost_matrix.shape[-1] if cost_matrix.ndim else 0),
            triangle_only=m == "HAM",
            stacklevel=3,
        )
    if (
        m == "DHD"
        and explicit_subm_value is not None
        and not isinstance(explicit_subm_value, str)
        and not diagnostics_shown
    ):
        dhd_costs = np.asarray(explicit_subm_value, dtype=float)
        # Only an array the DHD validator will accept: R stops on the others
        # before its diagnostics, and so does traminex.
        if (
            dhd_costs.ndim == 3
            and dhd_costs.shape[0] == np.asarray(seqmat).shape[1]
            and dhd_costs.shape[1] == dhd_costs.shape[2]
            and (active_state_count is None or dhd_costs.shape[1] == active_state_count)
            and bool(np.all(np.isfinite(dhd_costs)))
            and bool(np.all(dhd_costs >= 0))
            and bool(np.all(np.diagonal(dhd_costs, axis1=1, axis2=2) == 0))
        ):
            _py._warn_cost_quality(
                dhd_costs,
                labels=_cost_state_labels(source_obj, dhd_costs.shape[-1]),
                stacklevel=3,
            )

    norm = kwargs.get("norm", None)
    if isinstance(norm, str) and norm.lower() == "none":
        if m in {"OM", "LCS", "LCP", "RLCP"}:
            kwargs["norm"] = None
        elif m in {"NMS", "NMSMST", "SVRSPELL"}:
            kwargs.pop("norm")

    void_for_rows = kwargs.get("void", None)
    if (
        active_state_count is not None
        and void_for_rows is not None
        and np.any(np.asarray(seqmat) == void_for_rows)
    ):
        kwargs["_n_states"] = active_state_count
        if with_missing and nr_code is not None:
            kwargs["_missing_code"] = int(nr_code)

    refseq = kwargs.pop("refseq", None)
    ref_sets = kwargs.pop("ref_sets", None)
    if refseq is not None and ref_sets is not None:
        raise ValueError("refseq and ref_sets cannot be used together")

    validation_kwargs = dict(kwargs)
    for transport_key in (
        "void",
        "_n_states",
        "_missing_code",
        "_allow_zero_indel",
        "deduplicate",
        "_no_dedup",
        "condensed",
    ):
        validation_kwargs.pop(transport_key, None)
    validate_pairwise_kwargs(m, validation_kwargs, dispatcher=True, layer="seqdist")

    if m in {"CHI2", "EUCLID"}:
        if kwargs.get("breaks", None) is not None and not diagnostics_shown:
            # seqdist.R:256-259 warns and ignores step/overlap when breaks is given, with
            # msg.warn.ign2's text (TraMineR-logging_helpers.R:98-100; batch-9 review M-2).
            for ignored in ("step", "overlap"):
                if ignored in kwargs:
                    warnings.warn(
                        f" [!!] '{ignored}' has been ignored because 'breaks' is specified",
                        UserWarning,
                        stacklevel=2,
                    )
        requested_pdotj = kwargs.get("global_pdotj", None)
        if isinstance(requested_pdotj, str):
            # seqdist-CHI2.R:54-60: the only accepted string is "obs", resolved as
            # seqmeant(seqdata, weighted=..., with.missing=..., prop=TRUE). Left to
            # itself this reached float("obs") and raised "could not convert string
            # to float: 'obs'", naming neither the option nor its accepted values.
            if requested_pdotj != "obs":
                raise ValueError(
                    "global_pdotj should be either 'obs' or a vector of proportions "
                    "summing up to 1"
                )
            from ..stats_p1 import seqmeant

            # seqmeant's rows follow the whole codebook (real states, then the
            # nr/void markers), but the state space is the alphabet plus nr under
            # with_missing, so keep its leading rows only (audit S1a-01). A plain
            # matrix's alphabet is 0..max, exactly as its string cost methods use.
            if source_obj is not None:
                mean_source = x
            else:
                plain_states = int(seqmat.max()) + 1 if seqmat.size else 0
                plain_weights = kwargs.get("weights", None)
                mean_source = STSList(
                    data=np.asarray(seqmat, dtype=int).copy(),
                    alphabet=[str(code) for code in range(plain_states)],
                    weights=(
                        None if plain_weights is None else np.asarray(plain_weights, float).copy()
                    ),
                )
            state_space = (
                active_state_count if active_state_count is not None else len(mean_source.alphabet)
            )
            # seqmeant spells R's with.missing as include_missing.
            kwargs["global_pdotj"] = np.asarray(
                seqmeant(mean_source, weighted=weighted, include_missing=with_missing, prop=True),
                dtype=float,
            ).ravel()[:state_space]

    # ``seqdist`` accepts the full TraMineR normalization vocabulary for the classic
    # OM/LCS/LCP/RLCP family, plus the traminex-friendly aliases ``'length'`` and
    # ``'max'``. The low-level ``_map_norm`` reached through ``pairwise`` knows only
    # the latter two, so an invalid spelling used to be rejected with a message
    # naming three of the nine accepted values.
    if m in _CLASSIC_NORM_METHODS:
        requested_norm = kwargs.get("norm", None)
        if (
            requested_norm is not None
            and str(requested_norm).lower() not in _CLASSIC_NORM_SPELLINGS
        ):
            accepted = ", ".join(f"'{name}'" for name in _CLASSIC_NORM_DISPLAY[m])
            raise ValueError(
                f"norm must be one of: None, {accepted} (case-insensitive) for method {m}"
            )
        if isinstance(requested_norm, str):
            # The friendly aliases are exactly TraMineR's maxlength/maxdist
            # (docs/DISTANCES_NORMALIZATION.md); resolving them here keeps
            # them out of the legacy pairwise vocabulary, whose 'length'
            # divides LCS/LCP/RLCP distances by l1+l2 instead.
            alias = requested_norm.strip().lower()
            if alias in _NORM_ALIASES:
                kwargs["norm"] = _NORM_ALIASES[alias]

    # R stops when 'sm' is missing for the OM family and TWED
    # (seqdist.R: msg.stop.miss("sm")); only HAM (CONSTANT 1), DHD (TRATE,
    # cval=4) and LCS (CONSTANT 2) auto-generate their costs.
    if m in {"OM", "OMLOC", "OMSLEN", "OMSPELL", "OMSTRAN", "TWED"} and "subm" not in kwargs:
        raise ValueError(
            f"subm must be specified for method {m} (TraMineR: \"'sm' must be "
            'specified"); pass a matrix, a scalar cost, or one of '
            "TRATE, CONSTANT, INDELS, INDELSLOG"
        )

    if refseq is not None or ref_sets is not None:
        ref_kwargs = dict(kwargs)
        ref_kwargs.pop("condensed", None)
        n_rows = len(seqmat)

        if is_data_wrapper(refseq):
            ext_data = _py._validate_dense_sequence_matrix(refseq.data).astype(np.int32, copy=False)
            if ext_data.ndim != 2 or ext_data.shape[0] != 1:
                raise ValueError("external refseq must contain exactly one sequence")
            if source_obj is not None:
                source_alphabet = getattr(source_obj, "alphabet", None)
                ref_alphabet = getattr(refseq, "alphabet", None)

                def logical_alphabet(obj, alphabet):
                    return [alphabet[code] for code in _real_state_codes(obj)]

                if (
                    source_alphabet is not None
                    and ref_alphabet is not None
                    and logical_alphabet(source_obj, source_alphabet)
                    != logical_alphabet(refseq, ref_alphabet)
                ):
                    raise ValueError("external refseq must use the same alphabet")
                if getattr(source_obj, "nr", None) != getattr(refseq, "nr", None):
                    raise ValueError("external refseq must use the same missing-value symbol")
            elif getattr(refseq, "alphabet", None) is not None:
                # A plain array declares no alphabet, so R's identical-alphabet
                # rule reduces to this: every code the array uses as a state must
                # be a real state of the reference.  Otherwise one code would
                # mean a state in the source and nr (or nothing) in the reference.
                marker_codes = {
                    getattr(refseq, "nr_code", None),
                    getattr(refseq, "void_code", None),
                }
                ref_state_codes = {
                    code for code in range(len(refseq.alphabet)) if code not in marker_codes
                }
                source_void = ref_kwargs.get("void", None)
                foreign = sorted(
                    int(code)
                    for code in np.unique(np.asarray(seqmat))
                    if code != source_void and int(code) not in ref_state_codes
                )
                if foreign:
                    raise ValueError(
                        "external refseq must use the same alphabet: the sequences use "
                        f"state code(s) {foreign}, which are not states of the reference"
                    )
            ref_kwargs.pop("_n_states", None)
            ref_kwargs.pop("_missing_code", None)
            ref_kwargs.pop("deduplicate", None)
            ref_kwargs.pop("_no_dedup", None)
            base_void = ref_kwargs.pop("void", None)
            ext_void = getattr(refseq, "void_code", None)
            base_rows = [
                row[row != base_void] if base_void is not None else row
                for row in np.asarray(seqmat)
            ]
            ext_row = ext_data[0]
            if ext_void is not None:
                ext_row = ext_row[ext_row != ext_void]
            if len(ext_row) > max((len(row) for row in base_rows), default=0):
                raise ValueError("external refseq cannot be longer than the longest sequence")
            combined = base_rows + [ext_row]
            method_weights = ref_kwargs.get("weights", None)
            dhd_indels_method = (
                m == "DHD"
                and isinstance(ref_kwargs.get("subm", None), str)
                and ref_kwargs["subm"].upper() in {"INDELS", "INDELSLOG"}
            )
            if source_weights is not None:
                base_weights = (
                    np.asarray(source_weights, dtype=float)
                    if weighted or dhd_indels_method
                    else np.ones(n_rows, dtype=float)
                )
            elif method_weights is not None and m in {"CHI2", "OMSTRAN"}:
                base_weights = np.asarray(method_weights, dtype=float)
            else:
                base_weights = np.ones(n_rows, dtype=float)
            if base_weights.ndim != 1 or base_weights.shape[0] != n_rows:
                raise ValueError("weights length must equal number of source sequences")
            combined_weights = np.concatenate([base_weights, np.array([0.0])])
            ref_nr_code = getattr(refseq, "nr_code", None)
            combined_nr_code = nr_code if nr_code is not None else ref_nr_code
            combined_state_count = active_state_count
            if with_missing and combined_nr_code is not None:
                combined_state_count = max(
                    0 if combined_state_count is None else combined_state_count,
                    int(combined_nr_code) + 1,
                )
            # An explicit substitution matrix was bounded against the source
            # rows only (``_validate_fixed_substitution_costs`` uses
            # ``seqmat.max()``), so a reference contributing a higher code --
            # its missing state -- used to reach the kernel and fail with a bare
            # ``IndexError`` naming nothing (cycle-10 follow-on T-02). Reuse the
            # dispatcher's own wording now that the combined width is known.
            explicit_ref_subm = ref_kwargs.get("subm", None)
            if explicit_ref_subm is not None and not isinstance(explicit_ref_subm, str):
                ref_costs = np.asarray(explicit_ref_subm, dtype=float)
                if ref_costs.ndim == 2:
                    needed = max(
                        (int(row.max()) + 1 for row in combined if row.size),
                        default=0 if combined_state_count is None else combined_state_count,
                    )
                    if ref_costs.shape[0] < needed:
                        raise ValueError(
                            f"{m} subm must cover every observed state code "
                            f"(need at least {needed}x{needed})"
                        )
            # Only the pairs (sequence, reference) are needed, with the sequence
            # first: R's d(x_i, reference) for asymmetric costs
            # (TX-REP-20261017-001).
            with _py._only_pairs(n_rows + 1, range(n_rows), [n_rows]):
                cross = _py._pairwise_rows(
                    combined,
                    m,
                    ref_kwargs,
                    n_states=combined_state_count,
                    missing_code=(
                        int(combined_nr_code)
                        if with_missing and combined_nr_code is not None
                        else None
                    ),
                    weights=combined_weights,
                    reference_index=n_rows,
                    time_width=int(np.asarray(seqmat).shape[1]),
                )
            return cross[:n_rows, n_rows]

        vector_message = (
            "refseq index sets must be two collections of 1-based indices, such as "
            "([1, 2], [3]); a vector of indices is refused as R refuses c(1, 2) "
            "(TraMineR: 'refseq' invalid type or unknown value); one index is an integer"
        )
        if isinstance(refseq, np.ndarray) and refseq.ndim >= 1:
            # R's vector c(1, 2): neither one index nor a list of two sets
            # (seqdist.R:67-99; TX-AUD-20261011-044).
            raise ValueError(vector_message)
        if refseq is not None and not isinstance(
            refseq, (list, tuple, int, float, np.integer, np.floating)
        ):
            # R stops on any other refseq (seqdist.R:94-96); every refseq refusal
            # is a ValueError (batch-6 review M-4).
            raise ValueError("refseq must be an integer, two index sets, or one STS sequence")

        recursive_kwargs = dict(ref_kwargs)
        recursive_kwargs.pop("_n_states", None)
        # Deduplication would change the second call's rows, so the pair filter
        # would not apply and asymmetric costs would lose R's direction; the
        # filter already limits the work (batch-6 review I-2).
        recursive_kwargs.pop("deduplicate", None)
        if source_weights is not None and m in {"CHI2", "OMSTRAN"}:
            # These weights were injected above from the object; the recursive
            # call on the same object injects them again according to
            # ``weighted``, and an explicit copy would trip the conflict guard.
            recursive_kwargs.pop("weights", None)

        def reference_matrix(left, right) -> np.ndarray:
            # Only the pairs joining the reference rows to the others are
            # evaluated by the row engine; the other cells stay 0 and are sliced
            # away (TX-AUD-20261010-022). Dataset-fitted costs still see every row.
            with _py._only_pairs(n_rows, left, right):
                return seqdist(
                    source_obj if source_obj is not None else np.asarray(seqmat, dtype=np.int32),
                    method=method,
                    with_missing=with_missing,
                    weighted=weighted,
                    full_matrix=True,
                    _diagnostics_shown=True,
                    **recursive_kwargs,
                )

        if ref_sets is not None:
            idx_a, idx_b = ref_sets
            idx_a = _validated_reference_indices(
                idx_a, n_rows=n_rows, one_based=False, name="ref_sets"
            )
            idx_b = _validated_reference_indices(
                idx_b, n_rows=n_rows, one_based=False, name="ref_sets"
            )
            return reference_matrix(idx_a, idx_b)[np.ix_(idx_a, idx_b)]

        if isinstance(refseq, (list, tuple)):
            if any(np.isscalar(part) for part in refseq):
                raise ValueError(vector_message)
            if len(refseq) > 2:
                # seqdist.R: "Only first two elements of the 'refseq' list are used!"
                warnings.warn(
                    "Only the first two elements of the refseq list are used",
                    UserWarning,
                    stacklevel=2,
                )
                refseq = list(refseq[:2])
            if len(refseq) != 2:
                raise ValueError("refseq list must contain two 1-based index sets")
            idx_a = _validated_reference_indices(
                refseq[0], n_rows=n_rows, one_based=True, name="refseq index sets"
            )
            idx_b = _validated_reference_indices(
                refseq[1], n_rows=n_rows, one_based=True, name="refseq index sets"
            )
            return reference_matrix(idx_a, idx_b)[np.ix_(idx_a, idx_b)]

        if isinstance(refseq, (int, float, np.integer, np.floating)):
            ref_value = float(refseq)
            if not np.isfinite(ref_value) or ref_value != np.floor(ref_value):
                raise ValueError("refseq must be a non-negative integer index")
            ref_index = int(ref_value)
            if ref_index < 0:
                raise ValueError("refseq must be a non-negative integer index")
            if ref_index > n_rows:
                # Every refseq refusal is a ValueError (TX-AUD-20261011-044).
                raise ValueError(
                    "refseq index is outside the sequence data: it must be at most the "
                    "number of sequences"
                )
            if ref_index == 0:
                void = ref_kwargs.get("void", None)
                modal_rows = [
                    tuple(row[row != void]) if void is not None else tuple(row)
                    for row in np.asarray(seqmat)
                ]
                row_weights = (
                    np.ones(n_rows, dtype=float)
                    if source_weights is None
                    else np.asarray(source_weights, dtype=float)
                )
                counts: dict[tuple, float] = {}
                for row, weight in zip(modal_rows, row_weights):
                    counts[row] = counts.get(row, 0.0) + float(weight)

                # R calls ``seqtab(seqnum(seqdata), idxs=1)`` (seqdist.R:514, 570-574).
                # seqnum codes the real states 0..k-1; seqtab writes each sequence as
                # its compressed SPS text, a missing cell as ``*`` whatever the nr
                # spelling and void dropped, groups the weights with ``tapply`` (factor
                # levels in R's collation) and sorts decreasingly, then ``seqfind``
                # returns the first matching input row (TX-AUD-20261019-001). The keys
                # hold only digits, ``/``, ``-`` and ``*``, and the characters that can
                # meet at one position (a digit against ``*``, ``/`` or ``-``) order alike
                # in code points and in R's collation, so plain string order is R's.
                markers = {code for code in (nr_code, void) if code is not None}
                alphabet_size = (
                    len(source_obj.alphabet)
                    if source_obj is not None and getattr(source_obj, "alphabet", None) is not None
                    else None
                )
                real_rank: dict[int, int] = (
                    {}
                    if alphabet_size is None
                    else {
                        code: rank
                        for rank, code in enumerate(
                            code for code in range(alphabet_size) if code not in markers
                        )
                    }
                )

                def numeric_sps_key(row: tuple) -> str:
                    tokens: list[str] = []
                    position = 0
                    while position < len(row):
                        state = int(row[position])
                        if nr_code is not None and state == int(nr_code):
                            symbol = "*"
                        else:
                            symbol = str(real_rank.get(state, state))
                        end = position + 1
                        while end < len(row) and row[end] == row[position]:
                            end += 1
                        tokens.append(f"{symbol}/{end - position}")
                        position = end
                    return "-".join(tokens)

                # seqtab drops the sequences whose every cell is the nr code
                # (seqtab.R:17); a void cell is not nr, so such a row stays.
                non_missing: dict[tuple, bool] = {}
                for source_row, key in zip(np.asarray(seqmat), modal_rows):
                    keeps = nr_code is None or bool(np.any(source_row != int(nr_code)))
                    non_missing[key] = non_missing.get(key, False) or keeps
                eligible = [row for row in counts if non_missing[row]]
                if not eligible:
                    eligible = list(counts)
                max_count = max(counts[row] for row in eligible)
                selected = min(
                    (row for row in eligible if counts[row] == max_count),
                    key=numeric_sps_key,
                )
                ref_index = modal_rows.index(selected)
            else:
                ref_index -= 1
            return reference_matrix(range(n_rows), [ref_index])[:, ref_index]

    r_norm = str(kwargs.get("norm", "")).strip().lower()
    r_norm_requested = r_norm in {
        "auto",
        "maxlength",
        "gmean",
        "maxdist",
        "yujianbo",
    }

    # ``pairwise`` keeps its legacy classic-distance vocabulary.  Dense
    # ``seqdist`` calls using R's vocabulary run the compiled raw kernel when
    # their rows are rectangular, then apply R's normalization here.  Only
    # genuinely variable rows (or vector indels) use the row reference engine.
    if m in {"OM", "LCS", "LCP", "RLCP"} and r_norm_requested:
        classic_kwargs = dict(kwargs)
        condensed = bool(classic_kwargs.pop("condensed", False))
        void = classic_kwargs.pop("void", None)
        classic_kwargs.pop("_n_states", None)
        dense_input = np.asarray(seqmat, dtype=np.int32)
        has_void = void is not None and bool(np.any(dense_input == void))
        vector_indel = m == "OM" and not np.isscalar(classic_kwargs.get("indel", 1.0))
        if has_void or vector_indel:
            classic_kwargs.pop("deduplicate", None)
            classic_rows: list[np.ndarray] = [
                np.asarray(row[row != void] if void is not None else row, dtype=np.int32)
                for row in dense_input
            ]
            D = _py._pairwise_rows_classic(
                classic_rows, m, classic_kwargs, n_states=active_state_count
            )
        else:
            norm_arg = classic_kwargs.pop("norm")
            raw_kwargs = dict(classic_kwargs)
            raw_kwargs["norm"] = None
            D = pairwise(
                dense_input,
                method=m,
                condensed=False,
                void=void,
                **raw_kwargs,
            )
            length = int(dense_input.shape[1])
            code = resolve_norm_code(m, norm_arg)
            if m == "OM":
                indel = float(classic_kwargs.get("indel", 1.0))
                subm = classic_kwargs.get("subm", 2.0)
                if isinstance(subm, str):
                    # Defence in depth (F-010): unreachable today, because seqdist
                    # resolves every accepted cost string to an ndarray before this
                    # point and rejects an unknown one earlier. Kept so the branch
                    # stays correct if that resolution is ever restructured -- with
                    # a message naming all four strings seqdist accepts, not the
                    # 'TRATE'-only set of the low-level pairwise layer.
                    if subm.lower() != "trate":
                        raise TypeError(
                            "seqdist subm accepts float, ndarray, or one of: "
                            "TRATE, CONSTANT, INDELS, INDELSLOG"
                        )
                    subm = trate_substitution_matrix(
                        dense_input,
                        smoothing=float(classic_kwargs.get("trate_smoothing", 0.0)),
                        mode=str(classic_kwargs.get("trate_mode", "symmetric")),
                        normalize=str(classic_kwargs.get("trate_normalize", "none")),
                        weights=classic_kwargs.get("trate_weights", None),
                    )
                max_substitution = _py._subm_max_bound(subm)
                # Position-dependent costs (a traminex extension) replace the
                # scalar ones in the recurrence, so the bounds take their maxima,
                # as docs/OMLOC.md states for 'max' (audit S1b-02).
                position_subm = classic_kwargs.get("subm_by_pos", None)
                position_indel = classic_kwargs.get("indel_by_pos", None)
                if position_subm is not None and np.size(position_subm):
                    max_substitution = float(np.max(position_subm))
                if position_indel is not None and np.size(position_indel):
                    indel = float(np.max(position_indel))
                maxscost = _maxscost(max_substitution, indel, code)
                # OMdistance.cpp: maxpossiblecost = abs(n-m)*indel + maxscost*min(m,n);
                # the first term is 0 * Inf = NaN for an infinite indel (S1b-12).
                D = _normalize_square_output(
                    D,
                    maxdist=abs(length - length) * indel + length * maxscost,
                    length_1=length * indel,
                    length_2=length * indel,
                    code=code,
                )
            else:
                D = _normalize_square_output(
                    D,
                    maxdist=2.0 * length,
                    length_1=float(length),
                    length_2=float(length),
                    code=code,
                )
        return to_condensed(D) if condensed else D

    # HAM shares DHD's fixed-position normalizer.  Apply the complete R
    # vocabulary only at the seqdist boundary; pairwise retains its existing
    # raw/maxlength behavior.
    if m == "HAM" and r_norm_requested:
        ham_kwargs = dict(kwargs)
        condensed = bool(ham_kwargs.pop("condensed", False))
        void = ham_kwargs.pop("void", None)
        ham_kwargs.pop("_n_states", None)
        ham_kwargs.pop("deduplicate", None)
        norm_arg = ham_kwargs.pop("norm")
        ham_rows: list[np.ndarray] = [
            np.asarray(row[row != void] if void is not None else row, dtype=np.int32)
            for row in np.asarray(seqmat, dtype=np.int32)
        ]
        if len({len(row) for row in ham_rows}) > 1:
            raise ValueError("HAM is not defined for sequences of different length")
        dense = (
            np.vstack(ham_rows).astype(np.int32, copy=False)
            if ham_rows
            else np.asarray(seqmat, dtype=np.int32)
        )
        subm = ham_kwargs.get("subm", None)
        D = pairwise(dense, method="HAM", norm=None, **ham_kwargs)
        maxcost = 1.0 if subm is None else _py._subm_max_all(subm)
        maxdist = float(dense.shape[1]) * maxcost
        D = _py._normalize_fixed_position_output(D, maxdist, "HAM", norm_arg)
        return to_condensed(D) if condensed else D

    # DHD's generated time-varying TRATE costs need STS metadata to distinguish
    # the missing state from an ordinary highest-numbered state.  Keep all R
    # normalization at this boundary; raw ``pairwise`` calls retain their
    # historical meaning.
    if m == "DHD" and not _py._is_legacy_dhd(kwargs):
        void = kwargs.pop("void", None)
        dense = np.asarray(seqmat, dtype=np.int32)
        canonical_dense = dense
        if void is not None and np.any(dense == void):
            rows = [row[row != void] for row in dense]
            if len({len(row) for row in rows}) > 1:
                raise ValueError("DHD is not defined for sequences of different length")
            canonical_dense = np.vstack(rows).astype(np.int32, copy=False)
        transported_n_states = kwargs.pop("_n_states", None)
        kwargs.pop("_missing_code", None)
        condensed = bool(kwargs.pop("condensed", False))
        norm_arg = kwargs.pop("norm", None)
        subm = kwargs.pop("subm", None)
        observed_state_count = int(canonical_dense.max()) + 1 if canonical_dense.size else 0
        declared_state_count = (
            int(active_state_count)
            if active_state_count is not None
            else (int(transported_n_states) if transported_n_states is not None else None)
        )
        state_count = observed_state_count if declared_state_count is None else declared_state_count
        if state_count < observed_state_count:
            raise ValueError("declared state count cannot be smaller than observed DHD states")
        cost_method = subm.upper() if isinstance(subm, str) else None
        if subm is None or cost_method == "TRATE":
            _require_transition_rate_data(state_count, int(canonical_dense.shape[1]))
            dhd_weights = (
                np.asarray(source_weights, dtype=float)
                if weighted and source_weights is not None
                else None
            )
            # R generates DHD's costs through seqcost (seqsubm(..., "TRATE",
            # time.varying = TRUE), seqdist.R:459-463), which checks them once
            # (seqcost.R:286; TX-AUD-20261011-080).
            generated_dhd_costs = _py._dhd_costs_traminer(
                dense,
                cval=4.0,
                weights=dhd_weights,
                missing_code=None if nr_code is None else int(nr_code),
                n_states=state_count,
                ignore_code=None if void is None else int(void),
            )
            _py._warn_cost_quality(
                generated_dhd_costs[: canonical_dense.shape[1]],
                labels=_cost_state_labels(source_obj, int(generated_dhd_costs.shape[-1])),
                stacklevel=3,
            )
            D = _py._pairwise_dhd_traminer(
                dense,
                cval=4.0,
                weights=dhd_weights,
                missing_code=None if nr_code is None else int(nr_code),
                n_states=state_count,
                void_code=None if void is None else int(void),
                norm=norm_arg,
                costs=generated_dhd_costs,
            )
        elif cost_method in {"INDELS", "INDELSLOG"}:
            costs = _py._dhd_indels_costs(
                dense,
                method=cost_method,
                weights=(
                    np.asarray(source_weights, dtype=float) if source_weights is not None else None
                ),
                n_states=state_count,
                ignore_code=None if void is None else int(void),
            )
            logical_costs = costs[: canonical_dense.shape[1]]
            _py._warn_cost_quality(
                logical_costs,
                labels=_cost_state_labels(source_obj, int(logical_costs.shape[-1])),
                stacklevel=3,
            )
            if cost_method == "INDELS":
                # TX-AUD-20261019-002: a state absent from a position keeps its
                # stand-in out of that position's bound (diagnostics above use R's costs).
                logical_costs = _py._dhd_without_unobserved_stand_ins(
                    logical_costs, canonical_dense
                )
            D = _py._pairwise_dhd_explicit_subm(canonical_dense, logical_costs)
            maxdist = _py._dhd_max_cost_sum(logical_costs)
            D = _py._normalize_fixed_position_output(D, maxdist, "DHD", norm_arg)
        elif isinstance(subm, str):
            raise ValueError(_py._dhd_cost_string_message(subm))
        else:
            costs = _py._validate_dhd_cost_array(
                subm,
                time_width=int(dense.shape[1]),
                n_states=state_count,
            )
            logical_costs = costs[: canonical_dense.shape[1]]
            D = _py._pairwise_dhd_explicit_subm(canonical_dense, logical_costs)
            maxdist = _py._dhd_max_cost_sum(logical_costs)
            D = _py._normalize_fixed_position_output(D, maxdist, "DHD", norm_arg)
        return to_condensed(D) if condensed else D

    if m in {"OM", "OMSLEN", "OMSPELL", "TWED"} and not np.isscalar(kwargs.get("indel", 1.0)):
        condensed = bool(kwargs.pop("condensed", False))
        void = kwargs.pop("void", None)
        kwargs.pop("_n_states", None)
        kwargs.pop("_missing_code", None)
        kwargs.pop("deduplicate", None)
        variable_rows: list[np.ndarray] = [
            np.asarray(row[row != void] if void is not None else row, dtype=np.int32)
            for row in np.asarray(seqmat, dtype=np.int32)
        ]
        D = _py._pairwise_rows(
            variable_rows,
            m,
            kwargs,
            n_states=active_state_count,
            missing_code=(int(nr_code) if with_missing and nr_code is not None else None),
        )
        return to_condensed(D) if condensed else D

    void_code_for_rows = kwargs.get("void", None)
    has_void_cells = void_code_for_rows is not None and bool(
        np.any(np.asarray(seqmat) == void_code_for_rows)
    )
    if (
        active_state_count is not None
        and not has_void_cells
        and m in {"CHI2", "EUCLID", "NMS", "NMSMST", "SVRSPELL"}
    ):
        # These size global_pdotj/prox on the state space, which is the alphabet
        # (plus nr under with_missing) -- not the codes a corpus happens to use, as
        # the dense kernels infer (audit S1b-05). The row engine takes the declared
        # count directly, exactly as ragged input already reaches it.
        condensed = bool(kwargs.pop("condensed", False))
        for transport_key in ("void", "_n_states", "_missing_code", "deduplicate", "_no_dedup"):
            kwargs.pop(transport_key, None)
        dense_rows = np.asarray(seqmat, dtype=np.int32)
        D = _py._pairwise_rows(
            [row for row in dense_rows],
            m,
            kwargs,
            n_states=active_state_count,
            missing_code=(int(nr_code) if with_missing and nr_code is not None else None),
            time_width=int(dense_rows.shape[1]),
        )
        return to_condensed(D) if condensed else D

    return pairwise(np.asarray(seqmat, dtype=np.int32), method=method, **kwargs)


def to_condensed(D: np.ndarray) -> np.ndarray:
    """Convert a square distance matrix to condensed upper-triangular form.

    Parameters
    ----------
    D : np.ndarray
        Square matrix of shape ``(n, n)``.

    Returns
    -------
    np.ndarray
        Vector of shape ``(n*(n-1)/2,)`` containing upper-triangular entries.
    """
    D = np.asarray(D)
    if D.ndim != 2 or D.shape[0] != D.shape[1]:
        raise ValueError("D must be a square matrix")
    n = D.shape[0]
    iu = np.triu_indices(n, k=1)
    return D[iu]


def from_condensed(d: np.ndarray, n: Optional[int] = None) -> np.ndarray:
    """Convert condensed distance form back to a full symmetric matrix.

    Parameters
    ----------
    d : np.ndarray
        Condensed vector of shape ``(n*(n-1)/2,)``.
    n : int | None, default=None
        Number of rows/columns in reconstructed square matrix. If ``None``, inferred
        from ``len(d)``.

    Returns
    -------
    np.ndarray
        Square symmetric float matrix of shape ``(n, n)`` with zero diagonal.
    """
    d = np.asarray(d)
    if d.ndim != 1:
        raise ValueError("d must be a one-dimensional condensed vector")
    m = d.size
    if n is None:
        # Integer arithmetic avoids rounding the triangular-size inverse.
        n = (1 + isqrt(1 + 8 * m)) // 2
    elif isinstance(n, (bool, np.bool_)) or not isinstance(n, Integral):
        raise TypeError("n must be a non-negative integer")
    n = int(n)
    if n < 0:
        raise ValueError("n must be a non-negative integer")
    if n * (n - 1) // 2 != m:
        raise ValueError("Length of d is not a valid condensed size for n")
    D = np.zeros((n, n), dtype=float)
    iu = np.triu_indices(n, k=1)
    D[iu] = d
    D[(iu[1], iu[0])] = d
    return D


def sanitize_distance_matrix(D: np.ndarray) -> np.ndarray:
    """Sanitize a distance matrix for downstream algorithms.

    Parameters
    ----------
    D : np.ndarray
        Square matrix of shape ``(n, n)``.

    Returns
    -------
    np.ndarray
        Float matrix of shape ``(n, n)`` with:
        - non-finite values replaced,
        - symmetry enforced,
        - zero diagonal.
    """
    X = np.asarray(D, dtype=float)
    if X.ndim != 2 or X.shape[0] != X.shape[1]:
        raise ValueError("D must be a square matrix")
    # replace non-finite with max finite value
    finite = np.isfinite(X)
    if not finite.all():
        max_valid = np.nanmax(X[finite]) if np.any(finite) else 0.0
        X = np.where(finite, X, max_valid)
    # Sum first at ordinary/subnormal magnitudes to retain representable tiny
    # values. Halve first only where addition could overflow: these large
    # operands halve exactly; a tiny partner cannot affect their rounding.
    half_max = np.finfo(float).max / 2
    if np.max(X, initial=0.0) <= half_max and np.min(X, initial=0.0) >= -half_max:
        # Ordinary matrices need no per-cell masks or indexed temporaries.
        X = (X + X.T) * 0.5
    else:
        large = (np.abs(X) > half_max) | (np.abs(X.T) > half_max)
        average = np.empty_like(X)
        np.add(X, X.T, out=average, where=~large)
        np.multiply(average, 0.5, out=average, where=~large)
        average[large] = X[large] * 0.5 + X.T[large] * 0.5
        X = average
    np.fill_diagonal(X, 0.0)
    return X


def get_lcp_length(a: np.ndarray, b: np.ndarray) -> int:
    """Return the length of the longest common prefix between two sequences a and b."""
    a = np.asarray(a)
    b = np.asarray(b)
    L = int(min(a.shape[0], b.shape[0]))
    k = 0
    while k < L and a[k] == b[k]:
        k += 1
    return k


def get_rlcp_length(a: np.ndarray, b: np.ndarray) -> int:
    """Return the length of the longest common prefix computed on reversed sequences (RLCP)."""
    a = np.asarray(a)
    b = np.asarray(b)
    a_rev = a[::-1]
    b_rev = b[::-1]
    L = int(min(a_rev.shape[0], b_rev.shape[0]))
    k = 0
    while k < L and a_rev[k] == b_rev[k]:
        k += 1
    return k


def seqsubm(
    seqmat: STSList | np.ndarray,
    method: str | None = None,
    cval: float | None = None,
    **kwargs: Any,
) -> np.ndarray | pd.DataFrame | SeqCostMatrix | None:
    """Return public ``seqcost(...).sm`` costs or use the raw-array extension.

    Parameters
    ----------
    seqmat : STSList | np.ndarray
        A public state-sequence object, or a two-dimensional raw integer-coded
        sequence matrix with non-negative state codes fitting a 32-bit integer.
    method : str | None
        For an :class:`STSList`, the required public TraMineR ``seqcost`` method:
        R's ``seqcost`` has no default for it (``seqcost.R:296-298``), so
        ``None`` raises. For a raw array (a Python extension), ``None`` means
        ``"TRATE"``.
    cval : float | None
        Public ``seqcost`` base or constant cost. It may be passed positionally.
    **kwargs
        An :class:`STSList` forwards the complete public ``seqcost`` keyword
        surface unchanged. Raw arrays retain their method-specific Python
        extension keywords, including ``feature_types`` for ``FEATURES``.

    Returns
    -------
    pandas.DataFrame | SeqCostMatrix | numpy.ndarray | None
        The labelled ``sm`` component for an :class:`STSList`, ``None`` when
        ``proximities=True`` (R ``NULL``), or an unlabelled NumPy matrix for a
        raw integer array.
    """
    if isinstance(seqmat, STSList):
        if method is None:
            raise TypeError("seqsubm() missing required argument: 'method'")
        result = seqcost(seqmat, method, cval, **kwargs)
        value = result.get("sm")
        if value is None:
            return None
        if not isinstance(value, (pd.DataFrame, SeqCostMatrix)):
            raise TypeError("seqcost returned an invalid substitution-cost object")
        return value

    raw_method = "TRATE" if method is None else method
    raw_kwargs = dict(kwargs)
    # R's check.costs, as seqcost takes it for a sequence object (TX-AUD-20261011-041).
    check_costs_requested(raw_kwargs.pop("check_costs", True))
    if cval is not None:
        raw_kwargs["cval"] = cval
    return _seqsubm_from_array(np.asarray(seqmat), method=raw_method, **raw_kwargs)


def _seqsubm_from_array(seqmat: np.ndarray, method: str = "TRATE", **kwargs: Any) -> np.ndarray:
    """Build an unlabelled substitution-cost matrix from integer-coded sequences.

    This is the retained raw-array Python extension. Supported methods are:
        - ``TRATE``: transition-rate based matrix
        - ``CONSTANT``: constant off-diagonal costs
        - ``INDELS``: reciprocal-frequency based costs (sm[i,j] = 1/p_i + 1/p_j)
        - ``INDELSLOG``: log-frequency based costs (sm[i,j] = log(2/(1+p_i)) + log(2/(1+p_j)))
        - ``FEATURES``: Gower dissimilarity over a per-state feature table (needs ``state_features``)
        - ``FUTURE``: chi-squared distance on common-future transition rates
    """
    seqmat = _py._validate_dense_sequence_matrix(seqmat)
    m = str(method).upper()
    if m == "TRATE":
        from .costs import trate_substitution_matrix

        allowed = {"smoothing", "mode", "normalize", "weights"}
        unknown = sorted(set(kwargs) - allowed)
        if unknown:
            raise TypeError(
                f"Unsupported kwargs for seqsubm(method='TRATE'): {unknown}. Allowed kwargs: {sorted(allowed)}"
            )
        return trate_substitution_matrix(seqmat, **kwargs)
    if m == "CONSTANT":
        from .costs import constant_substitution_matrix

        if "constant" in kwargs:
            raise TypeError("constant is not supported; use cval=<float>")
        allowed = {"cval"}
        unknown = sorted(set(kwargs) - allowed)
        if unknown:
            raise TypeError(
                f"Unsupported kwargs for seqsubm(method='CONSTANT'): {unknown}. Allowed kwargs: {sorted(allowed)}"
            )
        cval = float(kwargs.get("cval", 2.0))
        return constant_substitution_matrix(seqmat, cval=cval)
    if m == "FUTURE":
        from .costs import future_costs

        allowed = {"weights"}
        unknown = sorted(set(kwargs) - allowed)
        if unknown:
            raise TypeError(
                f"Unsupported kwargs for seqsubm(method='FUTURE'): {unknown}. Allowed kwargs: {sorted(allowed)}"
            )
        _indel, sm = future_costs(seqmat, weights=kwargs.get("weights", None))
        return sm
    if m == "FEATURES":
        from .costs import gower_distance

        allowed = {"state_features", "feature_weights", "feature_types", "ordinal_levels"}
        unknown = sorted(set(kwargs) - allowed)
        if unknown:
            raise TypeError(
                f"Unsupported kwargs for seqsubm(method='FEATURES'): {unknown}. Allowed kwargs: {sorted(allowed)}"
            )
        state_features = kwargs.get("state_features", None)
        if state_features is None:
            raise ValueError("seqsubm(method='FEATURES') requires state_features=<[K,F] table>")
        K = int(seqmat.max()) + 1 if seqmat.size else 0
        feats = np.asarray(state_features, dtype=object)
        if feats.ndim != 2 or feats.shape[0] != K:
            raise ValueError(
                f"state_features must have one row per state (expected {K} rows, got {feats.shape[0] if feats.ndim == 2 else 'non-2D'})"
            )
        return gower_distance(
            feats,
            weights=kwargs.get("feature_weights", None),
            types=kwargs.get("feature_types", None),
            ordinal_levels=kwargs.get("ordinal_levels", None),
        )
    if m == "INDELS":
        from .costs import indels_costs

        allowed = {"include_missing"}
        unknown = sorted(set(kwargs) - allowed)
        if unknown:
            raise TypeError(
                f"Unsupported kwargs for seqsubm(method='INDELS'): {unknown}. Allowed kwargs: {sorted(allowed)}"
            )
        _indel, sm = indels_costs(
            seqmat, include_missing=bool(kwargs.get("include_missing", False))
        )
        return sm
    if m == "INDELSLOG":
        from .costs import indelslog_costs

        allowed = {"include_missing"}
        unknown = sorted(set(kwargs) - allowed)
        if unknown:
            raise TypeError(
                f"Unsupported kwargs for seqsubm(method='INDELSLOG'): {unknown}. Allowed kwargs: {sorted(allowed)}"
            )
        _indel, sm = indelslog_costs(
            seqmat, include_missing=bool(kwargs.get("include_missing", False))
        )
        return sm
    raise ValueError(
        f"Unknown seqsubm method: {method!r}. "
        "Supported: TRATE, CONSTANT, INDELS, INDELSLOG, FEATURES, FUTURE"
    )


def _llcs_extract(seq):
    """Return ``(codes, alphabet)`` for one ``seqLLCS``/``seqLLCP`` argument.

    Every ``void`` cell is dropped, while the missing state ``nr`` ("*") is kept
    as an ordinary state. R compares the first ``seqlength()`` cells, which is the
    same when the voids trail and differs for an interior void, which seqrecode
    can make (a recorded divergence, TX-AUD-20261011-005). Accepts a **single-sequence**
    :class:`STSList` (R's ``seqLLCS(seq[i, ], seq[j, ])``) or a 1-D/2-D array of
    integer codes (``alphabet`` is then ``None`` and the same-alphabet check is
    skipped).

    A multi-sequence object is refused. ``SOURCES/TraMineR/R/seqLLCS.R`` and
    ``seqLLCP.R`` reach ``.C(C_cLCS, as.integer(seq1), as.integer(seq2), ...)``,
    and a multi-row ``stslist`` is a ``data.frame`` -- hence a list -- which
    ``as.integer`` cannot coerce, so R stops. The message matches the
    ``exactly_one`` guard in ``traminex/stats_p2.py::_sequence_codes`` that
    already governs the sibling pair ``seqmpos``/``seqcomp``.
    """
    if hasattr(seq, "data") and hasattr(seq, "alphabet"):
        if isinstance(seq, STSList):
            seq = seq._canonicalized()
        data = np.asarray(seq.data)
        if data.ndim == 2 and data.shape[0] != 1:
            raise ValueError("sequence objects must contain exactly one sequence")
        row = np.asarray(seq.data[0])
        vc = getattr(seq, "void_code", None)
        if vc is not None:
            row = row[row != vc]
        # Storage positions differ across external/in-alphabet marker layouts.
        # Compare ordinary states by their real-alphabet rank and reserve the
        # next logical code for missing, independent of its storage code.
        real_codes = _real_state_codes(seq)
        code_map = {code: rank for rank, code in enumerate(real_codes)}
        nr = getattr(seq, "nr_code", None)
        if nr is not None:
            code_map[nr] = len(real_codes)
        logical = np.array([code_map[int(code)] for code in row], dtype=np.int64)
        return logical, [seq.alphabet[code] for code in real_codes]
    arr = np.asarray(seq)
    if arr.ndim == 2:
        # As for a multi-row object: silently taking row 0 hid the caller's
        # mistake (audit S1a-09).
        if arr.shape[0] != 1:
            raise ValueError("sequence arrays must contain exactly one sequence")
        arr = arr[0]
    if arr.ndim != 1:
        raise ValueError("sequence arrays must be one-dimensional or contain exactly one row")
    if arr.dtype.kind not in {"i", "u"}:
        raise TypeError("sequence arrays must contain integer state codes")
    return arr, None


def _lcs_length(a: np.ndarray, b: np.ndarray) -> int:
    """Length of the longest common subsequence of two integer-coded rows."""
    la, lb = int(a.shape[0]), int(b.shape[0])
    if la == 0 or lb == 0:
        return 0
    prev = np.zeros(lb + 1, dtype=np.int64)
    cur = np.zeros(lb + 1, dtype=np.int64)
    for i in range(1, la + 1):
        ai = a[i - 1]
        cur[0] = 0
        for j in range(1, lb + 1):
            if ai == b[j - 1]:
                cur[j] = prev[j - 1] + 1
            elif prev[j] >= cur[j - 1]:
                cur[j] = prev[j]
            else:
                cur[j] = cur[j - 1]
        prev, cur = cur, prev
    return int(prev[lb])


def seqLLCS(seq1, seq2) -> int:
    """Length of the longest common subsequence between two sequences (R ``seqLLCS``).

    Port of TraMineR 2.2-14 ``seqLLCS`` (GPL): a 2-argument function
    returning the scalar length of the longest common subsequence of ``seq1`` and
    ``seq2``. The comparison runs over the non-void cells of each sequence (the
    missing state ``nr`` is treated as an ordinary state); R compares the first
    ``seqlength()`` cells, so an interior void gives another value (a recorded
    divergence, TX-AUD-20261011-005). Both sequences must share the same alphabet.

    Parameters
    ----------
    seq1, seq2 : STSList | np.ndarray
        Sequence objects containing exactly one sequence each, or 1-D arrays of
        integer state codes. Multirow objects raise ValueError, as in R.

    Returns
    -------
    int
        Length of the longest common subsequence.

    Notes
    -----
    This is a per-pair scalar, not a distance matrix. The full LCS-based
    *distance* matrix is available via ``seqdist(method="LCS")`` / ``pairwise``.
    """
    a, alph1 = _llcs_extract(seq1)
    b, alph2 = _llcs_extract(seq2)
    if alph1 is not None and alph2 is not None and alph1 != alph2:
        raise ValueError("The alphabet of both sequences have to be same.")
    return _lcs_length(a, b)


def seqLLCP(seq1, seq2) -> int:
    """Length of the longest common prefix between two sequences (R ``seqLLCP``).

    Port of TraMineR 2.2-14 ``seqLLCP`` (GPL): a 2-argument function
    returning the scalar length of the longest common prefix of ``seq1`` and
    ``seq2``, compared over the non-void cells of each sequence (``nr`` kept as
    an ordinary state); R compares the first ``seqlength()`` cells, so an
    interior void gives another value (a recorded divergence,
    TX-AUD-20261011-005). Both sequences must share the same alphabet.

    Parameters
    ----------
    seq1, seq2 : STSList | np.ndarray
        Sequence objects containing exactly one sequence each, or 1-D arrays of
        integer state codes. Multirow objects raise ValueError, as in R.

    Returns
    -------
    int
        Length of the longest common prefix.

    Notes
    -----
    Per-pair scalar, not a distance matrix; the LCP-based *distance* matrix is
    available via ``seqdist(method="LCP")`` / ``pairwise``.
    """
    a, alph1 = _llcs_extract(seq1)
    b, alph2 = _llcs_extract(seq2)
    if alph1 is not None and alph2 is not None and alph1 != alph2:
        raise ValueError("The alphabet of both sequences have to be same.")
    return get_lcp_length(a, b)


def seqRLCP(x, *, norm: Optional[str] = None) -> np.ndarray:
    """Wrapper for RLCP-based distance matrix (reverse LCP on sequences).

    .. note::

       **traminex-only extension with no TraMineR counterpart.** Despite sharing
       their prefix, this function does *not* follow the ``seqLLCS``/``seqLLCP``
       contract. Those two are ports of R functions that take **two single
       sequences** and return a **scalar**; ``seqRLCP`` takes **one whole
       sequence object** and returns the full ``(n, n)`` **distance matrix**.
       ``exists("seqRLCP", where = asNamespace("TraMineR"))`` is ``FALSE``, and
       the same holds for ``TraMineRextras``, so there is no R behaviour to be
       faithful to. The name is retained because the public API is frozen at 226
       names (``ECOSYSTEM.md``); the divergence is recorded here rather than
       resolved by a rename. For the scalar reverse-LCP length of two sequences
       use :func:`get_rlcp_length`; for the matrix, this function gives
       ``seqdist(method="RLCP", with_missing=True)``'s raw distances.

    Parameters
    ----------
    x : STSList | np.ndarray
        STS-like object or integer-coded matrix of shape ``(n, L)``.
    norm : None | "length"
        Optional normalization mode. When "length", divides by the sum of the two
        sequence lengths (``2 * L`` for equal lengths) -- the low-level
        ``pairwise`` meaning, not ``seqdist``'s ``maxlength`` alias.

    Returns
    -------
    np.ndarray
        Float square distance matrix of shape ``(n, n)``.
    """
    if norm is not None and str(norm).lower() != "length":
        raise ValueError("seqRLCP norm must be None or 'length'")
    if isinstance(x, STSList):
        # seqdist deletes void (a ragged object's padding is not a state) and
        # canonicalizes the markers; the missing state stays a compared state,
        # as it always was here (audit S1a-08).
        canonical = x._canonicalized()
        # with_missing only when there is something missing: otherwise seqdist prints
        # R's notice that it resets it, which R's helpers never print (TX-AUD-20261019-025).
        has_missing = canonical.nr_code is not None and bool(
            np.any(np.asarray(canonical.data) == canonical.nr_code)
        )
        raw = np.asarray(seqdist(canonical, method="RLCP", with_missing=has_missing), dtype=float)
        void_code = canonical.void_code
        lengths = np.array(
            [
                int(np.sum(row != void_code)) if void_code is not None else len(row)
                for row in np.asarray(canonical.data)
            ],
            dtype=float,
        )
        if norm is None:
            return raw
        denominator = np.add.outer(lengths, lengths)
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(denominator > 0, raw / denominator, 0.0)
    seqmat = np.asarray(x.data if hasattr(x, "data") else x)
    return pairwise(seqmat, method="RLCP", norm=norm)
