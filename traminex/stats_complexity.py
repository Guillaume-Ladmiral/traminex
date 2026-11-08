from __future__ import annotations

import sys
from collections.abc import Sequence
from typing import Tuple

import numpy as np
from numpy.lib.mixins import NDArrayOperatorsMixin

from ._deprecated import _caller_stacklevel
from ._state_order import StateOrderContract, _resolve_state_order
from .sequences import STSList, _alphabet_codes, _require_stslist
from .sequences import alphabet as sequence_alphabet
from .stats_basic import seqistatd, seqtransn
from .stats_p2 import seqpm as seqpm


class SeqIVarDurResult(NDArrayOperatorsMixin):
    """Array-compatible spell-duration variances with TraMineR metadata."""

    __array_priority__ = 1000

    _values: np.ndarray
    vmax: np.ndarray
    meand: np.ndarray

    def __init__(
        self,
        values: np.ndarray,
        *,
        vmax: np.ndarray,
        meand: np.ndarray,
    ) -> None:
        array = np.asarray(values, dtype=float)
        self._values = np.frombuffer(array.tobytes(order="C"), dtype=float).reshape(array.shape)
        for name, metadata in (("vmax", vmax), ("meand", meand)):
            array = np.asarray(metadata, dtype=float)
            if array.shape != self._values.shape:
                raise ValueError(f"{name} must match the variance shape")
            setattr(self, name, array.copy())

    @property
    def shape(self) -> tuple[int, ...]:
        return self._values.shape

    @property
    def dtype(self) -> np.dtype:
        return self._values.dtype

    @property
    def ndim(self) -> int:
        return self._values.ndim

    @property
    def size(self) -> int:
        return self._values.size

    def __array__(self, dtype=None, copy=None) -> np.ndarray:
        if copy is None:
            return np.asarray(self._values, dtype=dtype)
        return np.array(self._values, dtype=dtype, copy=copy)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        """Evaluate NumPy operations as plain arrays without stale metadata."""
        arrays = tuple(
            value._values if isinstance(value, SeqIVarDurResult) else value for value in inputs
        )
        if kwargs.get("out") is not None:
            if any(isinstance(value, SeqIVarDurResult) for value in kwargs["out"]):
                raise TypeError("in-place operations are not supported for SeqIVarDurResult")
            kwargs["out"] = tuple(
                value._values if isinstance(value, SeqIVarDurResult) else value
                for value in kwargs["out"]
            )
        return getattr(ufunc, method)(*arrays, **kwargs)

    def __getitem__(self, key):
        values = self._values[key]
        if np.ndim(values) == 0:
            return values
        return SeqIVarDurResult(values, vmax=self.vmax[key], meand=self.meand[key])

    def __len__(self) -> int:
        return len(self._values)

    def __iter__(self):
        return iter(self._values)

    def __getattr__(self, name: str):
        """Delegate ndarray conveniences as metadata-free operations."""
        if name.startswith("_"):
            # copy and pickle probe dunders (``__setstate__``, ``__deepcopy__``)
            # before ``_values`` exists; delegating them recursed forever or
            # copied only the array (TX-AUD-20261010-053).
            raise AttributeError(name)
        if name in {
            "fill",
            "partition",
            "put",
            "resize",
            "setfield",
            "setflags",
            "sort",
        }:
            raise TypeError(f"mutating ndarray method {name}() is not supported")
        return getattr(self._values, name)

    def __repr__(self) -> str:
        return f"SeqIVarDurResult({self._values!r}, vmax={self.vmax!r}, meand={self.meand!r})"

    def copy(self) -> SeqIVarDurResult:
        return SeqIVarDurResult(self._values, vmax=self.vmax, meand=self.meand)


@_require_stslist
def seqtransn_count(x: STSList, *, normalized: bool = True) -> np.ndarray:
    """Per-sequence number of transitions, as R's ``seqtransn`` counts them (NOT the Elzinga turbulence).

    Counts the state changes of each sequence with its missing cells left out of
    the spells, as R's ``seqtransn`` does without ``with.missing``; entering the
    void padding of a shorter sequence is not a change. If normalized=True
    (default), divides by the sequence's own length minus 1 (missing cells
    counted, void not), as ``seqtransn(norm = TRUE)``; a sequence of length 1 or
    less gives 0. Equal to ``seqtransn(x, norm=normalized)`` (TX-AUD-20261011-068).

    Note
    ----
    This is a simple transition-count proxy, not the canonical turbulence score.
    For the standard Elzinga turbulence index (as in TraMineR's seqST), use
    ``seqST`` or ``turbulence_exact`` instead.
    """
    # It counted changes into and out of the missing marker and divided by the
    # padded width; seqtransn is R's count and normalization.
    return np.asarray(seqtransn(x, norm=normalized), dtype=float)


@_require_stslist
def turbulence_proxy(x: STSList, *, normalized: bool = True) -> float:
    """Weighted mean of the transition-count proxy across all sequences.

    This is a dataset-level summary of ``seqtransn_count``. It is NOT the
    canonical Elzinga turbulence. For the true turbulence summary use
    ``turbulence_exact`` (wraps ``seqST``).
    """
    per = seqtransn_count(x, normalized=normalized)
    if x.weights is None:
        return float(per.mean()) if per.size > 0 else 0.0
    w = np.asarray(x.weights, dtype=float)
    s = w.sum()
    return float((per * w).sum() / s) if s > 0 else 0.0


def _dss_and_durations(
    row: np.ndarray,
    *,
    include_missing: bool,
    nr_code: int | None,
    void_code: int | None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return (dss_codes, run_lengths) for a single integer-coded sequence row."""
    if not include_missing:
        bad = []
        if nr_code is not None:
            bad.append(int(nr_code))
        if void_code is not None:
            bad.append(int(void_code))
        if bad:
            mask = np.ones(row.shape[0], dtype=bool)
            for b in bad:
                mask &= row != b
            row = row[mask]
    if row.size == 0:
        return np.array([], dtype=int), np.array([], dtype=int)
    dss = [int(row[0])]
    lens = []
    cur = int(row[0])
    cnt = 1
    for v in row[1:]:
        vv = int(v)
        if vv == cur:
            cnt += 1
        else:
            lens.append(cnt)
            dss.append(vv)
            cur = vv
            cnt = 1
    lens.append(cnt)
    return np.array(dss, dtype=int), np.array(lens, dtype=int)


def _strip_void(row: np.ndarray, void_code: int | None) -> np.ndarray:
    """Drop structural void (right-padding) positions from a coded row.

    R's seqlength/seqdss/seqdur always exclude ``void`` (it is padding, never a
    state). Removing it up front makes the seqdss/seqdur/seqivolatility ports
    faithful for ragged (unequal-length, void-padded) sequences as well as
    rectangular ones, without touching the shared ``_dss_and_durations`` helper.
    """
    if void_code is None:
        return row
    return row[row != void_code]


@_require_stslist
def seqdur(x: STSList, *, include_missing: bool = False) -> np.ndarray:
    """Durations of the successive distinct states (TraMineR seqdur).

    Returns a float array of shape (n_sequences, max_spells), right-padded with
    NaN. ``include_missing`` maps to R's ``with.missing``: when False, missing
    (nr) positions do not count toward a spell and do not break a run. Structural
    ``void`` (right-padding) is always dropped (R excludes it in both modes).
    """
    per: list[np.ndarray] = []
    max_sp = 0
    for i in range(x.nseq):
        row = _strip_void(x.data[i], x.void_code)
        _, lens = _dss_and_durations(
            row,
            include_missing=include_missing,
            nr_code=x.nr_code,
            void_code=None,
        )
        per.append(lens)
        if lens.size > max_sp:
            max_sp = lens.size
    out = np.full((x.nseq, max(max_sp, 1)), np.nan, dtype=float)
    for i, lens in enumerate(per):
        if lens.size:
            out[i, : lens.size] = lens
    return out


@_require_stslist
def seqdss(x: STSList, *, include_missing: bool = False) -> STSList:
    """Distinct successive states (TraMineR seqdss).

    Collapses each run of equal states to a single state, returning a new
    STSList over the same alphabet, right-padded with the void symbol. With
    ``include_missing=False``, missing (nr) positions do not break a run.
    Structural ``void`` (right-padding) is always dropped before collapsing.
    """
    x = x._canonicalized()
    per: list[np.ndarray] = []
    max_sp = 0
    for i in range(x.nseq):
        row = _strip_void(x.data[i], x.void_code)
        dss, _ = _dss_and_durations(
            row,
            include_missing=include_missing,
            nr_code=x.nr_code,
            void_code=None,
        )
        per.append(dss)
        if dss.size > max_sp:
            max_sp = dss.size

    out_alph = list(x.alphabet)
    void_symbol = x.void
    void_code = x.void_code
    if void_code is None:
        # An undeclared marker spelling may already name a real state.
        # Give new DSS padding its own collision-free identity.
        while void_symbol in out_alph or void_symbol == x.nr:
            void_symbol += "%"
        void_code = len(out_alph)
        out_alph.append(void_symbol)
    nr_code = x.nr_code

    data = np.full((x.nseq, max(max_sp, 1)), void_code, dtype=int)
    for i, dss in enumerate(per):
        if dss.size:
            data[i, : dss.size] = dss

    label_by_state = dict(zip(x.alphabet, x.labels)) if x.labels is not None else None
    color_by_state = dict(zip(x.alphabet, x.cpal)) if x.cpal is not None else None

    return STSList(
        data=data,
        alphabet=out_alph,
        labels=(
            None
            if label_by_state is None
            else [label_by_state.get(state, state) for state in out_alph]
        ),
        cpal=(
            None
            if color_by_state is None
            else [
                color_by_state.get(state, "white" if state == void_symbol else x.missing_color)
                for state in out_alph
            ]
        ),
        nr=x.nr,
        void=void_symbol,
        nr_code=nr_code,
        void_code=void_code,
        weights=None if x.weights is None else x.weights.copy(),
        start=x.start,
        xtlab=[f"ST{position}" for position in range(1, data.shape[1] + 1)],
        xtstep=x.xtstep,
        tick_last=x.tick_last,
        ids=None if x.ids is None else list(x.ids),
        missing_color=x.missing_color,
        issues=x.issues,
    )


@_require_stslist
def seqivolatility(
    x: STSList, *, w: float = 0.5, include_missing: bool = False, adjust: bool = True
) -> np.ndarray:
    """Within-sequence volatility (TraMineR seqivolatility).

    ``w * visit_term + (1 - w) * transp`` where ``transp`` is the per-sequence
    normalized transition count ((#DSS - 1)/(seqlength - 1)) and ``visit_term`` is
    ``(nvisit - 1)/(alph_size - 1)`` (0 if nvisit <= 1) when ``adjust=True``,
    else ``nvisit / alph_size``. ``include_missing`` maps to R's ``with.missing``.
    Structural ``void`` (right-padding) is excluded everywhere: ``seqlength`` is
    the per-row non-void length, the DSS drops void, and ``nvisit`` ignores it.
    ``adjust`` must be a boolean, as in R.
    """
    x = x._canonicalized()
    if not isinstance(adjust, (bool, np.bool_)):
        # seqivolatility.R:5 (audit S3-17).
        raise TypeError("adjust should be logical")
    if not (0.0 <= float(w) <= 1.0):
        raise ValueError("w should be in the range [0, 1]")

    missing_codes = {c for c in (x.nr_code, x.void_code) if c is not None}
    k_real = len(x.alphabet) - len(missing_codes)
    alph_size = k_real + (1 if include_missing else 0)

    sdist = seqistatd(x, include_missing=include_missing)
    # void is structural padding, never a visited state -> exclude it from nvisit
    if x.void_code is not None and 0 <= x.void_code < sdist.shape[1]:
        sdist = sdist.copy()
        sdist[:, x.void_code] = 0
    nvisit = (sdist > 0).sum(axis=1).astype(float)

    transp = np.zeros(x.nseq, dtype=float)
    for i in range(x.nseq):
        row = _strip_void(x.data[i], x.void_code)
        seql = row.shape[0]  # R seqlength (with.missing=TRUE default) = ncol - #void
        dss, _ = _dss_and_durations(
            row,
            include_missing=include_missing,
            nr_code=x.nr_code,
            void_code=None,
        )
        ntrans = float(max(dss.size - 1, 0))
        denom = float(seql - 1)
        transp[i] = (ntrans / denom) if denom > 0 else 0.0

    if adjust:
        if alph_size > 1:
            visit_term = np.where(nvisit - 1.0 <= 0.0, 0.0, (nvisit - 1.0) / (alph_size - 1.0))
        else:
            visit_term = np.zeros(x.nseq, dtype=float)
        ret = w * visit_term + (1.0 - w) * transp
    else:
        ret = w * (nvisit / alph_size) + (1.0 - w) * transp
    return np.asarray(ret, dtype=float)


def _num_distinct_subseq(dss: np.ndarray) -> int:
    """Number of distinct subsequences (including empty) of the DSS sequence."""
    N = 1
    last_val_to_N: dict[int, int] = {}
    for sym in dss:
        symbol = int(sym)
        prev_N = N
        if symbol in last_val_to_N:
            N = 2 * N - last_val_to_N[symbol]
        else:
            N = 2 * N
        last_val_to_N[symbol] = prev_N
    return int(N)


def _subsequence_count_as_float(count: int) -> float:
    """Round an exact nonnegative count to binary64, including overflow to Inf.

    The integer recurrence avoids R's intermediate Inf - Inf for repetitions.
    Its final public float result still overflows as IEEE arithmetic does.
    """
    try:
        return float(count)
    except OverflowError:
        return float("inf")


def _dur_stats(
    lens: np.ndarray, total_len: int, *, type: int = 1, K: int = 0, nvisit: int = 0
) -> Tuple[float, float, float]:
    """Return variance, maximum variance, and mean spell duration.

    ``type=1`` is the realvar of the spell durations with population denominator
    ``d`` (number of spells). ``type=2`` (R seqivardur type=2) additionally treats
    each non-visited alphabet state as a zero-duration spell: the denominator is
    ``d + nnvisit`` where ``nnvisit = K - nvisit`` (alphabet states never visited),
    and ``var.max`` uses ``maxnnv`` (``K-1`` when ``d==1`` else ``K-2``, 0 when
    ``K<2``). ``K`` is the effective alphabet size, ``nvisit`` the number of
    distinct visited states.
    """
    d = int(lens.size)
    if d == 0:
        # A row with no spell (every cell missing, include_missing=False). R's
        # type 1 is realvar/mean of an empty vector: NaN. Type 2 counts the K
        # unvisited states as zero-length spells, so variance and mean are 0
        # and var.max is 0/(K-2): 0 for K >= 3, NaN (0/0) for K <= 2
        # (seqivardur.R:25-40; audit S3-11).
        if type == 2 and K > 0:
            return 0.0, (0.0 if K >= 3 else float("nan")), 0.0
        return float("nan"), float("nan"), float("nan")

    if type == 2:
        nnvisit = max(int(K) - int(nvisit), 0)
        denom = float(d + nnvisit)
        meand = float(total_len) / denom if denom > 0 else 0.0
        s2 = (float(np.sum((lens.astype(float) - meand) ** 2)) + nnvisit * (meand**2)) / denom
        if K < 2:
            maxnnv = 0
        else:
            maxnnv = (K - 1) if d == 1 else (K - 2)
        denom_max = float(d + maxnnv)
        meand_max = float(total_len) / denom_max if denom_max > 0 else 0.0
        s2max = (
            (d - 1) * (1.0 - meand_max) ** 2
            + (float(total_len) - d + 1 - meand_max) ** 2
            + maxnnv * (meand_max**2)
        ) / denom_max
    else:
        meand = float(total_len) / float(d)
        s2 = float(np.sum((lens.astype(float) - meand) ** 2) / float(d))
        if d == 1:
            s2max = 0.0
        else:
            a = float(total_len - (d - 1))
            s2max = (((a - meand) ** 2) + (d - 1) * ((1.0 - meand) ** 2)) / float(d)

    return s2, float(s2max), meand


def _dur_var_and_max(
    lens: np.ndarray, total_len: int, *, type: int = 1, K: int = 0, nvisit: int = 0
) -> Tuple[float, float]:
    """Return the two duration statistics consumed by ``seqST``."""
    variance, vmax, _ = _dur_stats(lens, total_len, type=type, K=K, nvisit=nvisit)
    return variance, vmax


@_require_stslist
def seqST(
    x: STSList,
    *,
    norm: bool = False,
    include_missing: bool = False,
    type: int = 1,
    silent: bool = True,
) -> np.ndarray:
    """Exact Elzinga sequence turbulence per sequence (TraMineR seqST).

    ``type`` must be 1 or 2, as in R. ``silent=False`` emits TraMineR's
    extraction and turbulence progress messages.

    Exact integer subsequence counts are rounded to float64 before evaluating
    turbulence, including Inf and normalized NaN when the count overflows.
    This avoids the upstream repeated-state Inf - Inf / undefined turb.phi
    failure; that large-count correction is documented in the semantic ledger.
    """
    x = x._canonicalized()
    if not isinstance(silent, (bool, np.bool_)):
        raise TypeError("silent must be a boolean")
    if type not in (1, 2):
        # R forwards type to seqivardur, which stops (seqivardur.R:5; audit S3-17).
        raise ValueError("type should be 1 or 2")
    n = x.nseq
    if not silent:
        print(" [>] extracting symbols and durations ...", file=sys.stderr)
        print(
            f" [>] computing turbulence type {type} for {n} sequence(s) ...",
            file=sys.stderr,
        )
    out = np.zeros(n, dtype=float)

    # Effective alphabet size (R seqST): structural void is never a state, and the
    # missing state counts only when include_missing. Mirrors seqivolatility's
    # alph_size and R's `alph <- alphabet(seqdata); if (with.missing) c(alph, nr)`.
    missing_codes = {c for c in (x.nr_code, x.void_code) if c is not None}
    k_real = len(x.alphabet) - len(missing_codes)
    K = k_real + (1 if include_missing else 0)

    # Per-row void-free views: R's seqdss/seqdur/seqlength always drop void.
    rows = [_strip_void(x.data[i], x.void_code) for i in range(n)]

    maxT = None
    if norm:
        # R: maxlength = max(seqlength(seqdata)) (default with.missing=TRUE = the
        # max non-void length); the reference alphabet has K distinct states.
        Lmax = max((int(r.shape[0]) for r in rows), default=0)
        if Lmax > 0 and K >= 1:
            labels_idx = np.arange(K, dtype=int)
            reference_codes: list[int] = []
            while len(reference_codes) < Lmax:
                reference_codes.extend(labels_idx.tolist())
            reference_sequence = np.array(reference_codes[:Lmax], dtype=int)
            # Every reference symbol is a regular state (no nr/void to strip).
            dss_m, lens_m = _dss_and_durations(
                reference_sequence, include_missing=True, nr_code=None, void_code=None
            )
            phi_m = _subsequence_count_as_float(_num_distinct_subseq(dss_m))
            nvisit_m = int(np.unique(dss_m).size) if dss_m.size > 0 else 0
            s2_m, s2max_m = _dur_var_and_max(
                lens_m, int(lens_m.sum()), type=type, K=K, nvisit=nvisit_m
            )
            maxT = (
                float(np.log2(phi_m * ((s2max_m + 1.0) / (s2_m + 1.0))))
                if (s2_m + 1.0) > 0
                else 1.0
            )
        else:
            maxT = 1.0

    for i in range(n):
        row = rows[i]
        dss, lens = _dss_and_durations(
            row, include_missing=include_missing, nr_code=x.nr_code, void_code=None
        )
        total_len = int(lens.sum()) if lens.size > 0 else 0
        phi = _subsequence_count_as_float(_num_distinct_subseq(dss))
        nvisit = int(np.unique(dss).size) if dss.size > 0 else 0
        s2, s2max = _dur_var_and_max(lens, total_len, type=type, K=K, nvisit=nvisit)
        # s2 is a variance (>= 0) or NaN for a row without a spell; NaN propagates.
        T = float(np.log2(phi * ((s2max + 1.0) / (s2 + 1.0))))
        if norm and maxT is not None:
            # R: Tx < 1 -> 0, otherwise (Tx - 1) / (maxT - 1) in IEEE arithmetic, so a
            # degenerate maximum (maxT == 1: one-state alphabet or length-1 data) gives
            # NaN and a NaN turbulence stays NaN (seqST.R:93-99; audit S3-10).
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = np.float64(T - 1.0) / np.float64(maxT - 1.0)
            out[i] = 0.0 if T < 1.0 else float(ratio)
        else:
            out[i] = T
    return out


@_require_stslist
def turbulence_exact(
    x: STSList, *, norm: bool = True, include_missing: bool = False, type: int = 1
) -> float:
    """Weighted mean of exact turbulence (seqST) across sequences."""
    per = seqST(x, norm=norm, include_missing=include_missing, type=type)
    if x.weights is None:
        return float(per.mean()) if per.size > 0 else 0.0
    w = np.asarray(x.weights, dtype=float)
    s = w.sum()
    return float((per * w).sum() / s) if s > 0 else 0.0


def seqivardur(x: STSList, *, type: int = 1, include_missing: bool = False) -> SeqIVarDurResult:
    """Return spell-duration variance plus TraMineR's ``vmax`` and ``meand``."""
    if not isinstance(x, STSList):
        raise TypeError("x must be a state sequence object")
    x = x._canonicalized()
    if type not in (1, 2):
        raise ValueError("type should be 1 or 2")
    dur = seqdur(x, include_missing=include_missing)  # (n, maxsp), NaN-padded
    if not np.any(~np.isnan(dur)):
        # TraMineR refuses this object through a raw seqdss subscript error.
        # Keep the same boundary but expose the representation issue clearly.
        raise ValueError("seqivardur requires at least one non-void spell")
    missing_codes = {
        c for c in (x.nr_code, x.void_code) if c is not None and 0 <= c < len(x.alphabet)
    }
    k_real = len(x.alphabet) - len(missing_codes)
    K = k_real + (1 if include_missing else 0)
    sdist = seqistatd(x, include_missing=include_missing)
    if x.void_code is not None and 0 <= x.void_code < sdist.shape[1]:
        sdist = sdist.copy()
        sdist[:, x.void_code] = 0
    nvisit = (sdist > 0).sum(axis=1).astype(int)
    out = np.full(x.nseq, np.nan, dtype=float)
    vmax = np.full(x.nseq, np.nan, dtype=float)
    meand = np.full(x.nseq, np.nan, dtype=float)
    for i in range(x.nseq):
        lens = dur[i][~np.isnan(dur[i])].astype(int)
        # An empty row inside an otherwise supported object gets R's values from
        # _dur_stats: NaN for type 1; type 2 adds K zero-duration states.
        out[i], vmax[i], meand[i] = _dur_stats(
            lens,
            int(lens.sum()),
            type=type,
            K=K,
            nvisit=int(nvisit[i]),
        )
    return SeqIVarDurResult(out, vmax=vmax, meand=meand)


def seqlogp(
    x: STSList,
    *,
    prob: str | np.ndarray = "trate",
    time_varying: bool = True,
    begin: str | np.ndarray = "freq",
    weighted: bool = True,
    include_missing: bool = False,
) -> np.ndarray:
    """Negative log-probability under TraMineR's sequence models.

    ``prob`` accepts fitted transition rates, fitted state frequencies, a
    stationary K-by-K matrix, or a K-by-K-by-time array. ``begin`` accepts the
    first-position distribution, the global state distribution, or a supplied
    length-K vector. As in TraMineR, zero/negative/NaN probabilities contribute
    nothing and positive infinity contributes ``-Inf``; no smoothing is added.
    """
    if not isinstance(x, STSList):
        raise TypeError("x must be a state sequence object")
    x = x._canonicalized()
    if not include_missing and x.nr_code is not None and bool((x.data == x.nr_code).any()):
        raise ValueError("Non-void missing values in seqdata, include_missing should be True!")

    rows = [_strip_void(x.data[i], x.void_code) for i in range(x.nseq)]
    nonempty = [row for row in rows if row.size]
    observed = np.concatenate(nonempty) if nonempty else np.array([], dtype=int)
    if observed.size == 0:
        return np.zeros(x.nseq, dtype=float)
    nstates = int(observed.max()) + 1
    max_length = max((row.size for row in rows), default=0)
    ntransitions = max(max_length - 1, 0)

    if weighted and x.weights is not None:
        weights = np.asarray(x.weights, dtype=float)
        if weights.ndim != 1 or weights.size != x.nseq:
            raise ValueError(f"sequence weights must have length {x.nseq}")
    else:
        weights = np.ones(x.nseq, dtype=float)

    # Resolve the starting distribution before the probability model, matching
    # R's validation order when both supplied arrays are malformed.
    if isinstance(begin, str):
        if begin == "freq":
            from .stats_basic import seqstatd

            frequencies = seqstatd(
                x,
                weighted=weighted,
                include_missing=include_missing,
            )["Frequencies"]
            first_probabilities = (
                np.asarray(frequencies[0, :nstates], dtype=float)
                if frequencies.shape[0]
                else np.zeros(nstates, dtype=float)
            )
        elif begin == "global.freq":
            counts = np.asarray(seqistatd(x, include_missing=include_missing), dtype=float)
            totals = (counts * weights[:, None]).sum(axis=0)[:nstates]
            denominator = float(totals.sum())
            first_probabilities = totals / denominator if denominator > 0 else totals
        else:
            raise ValueError("begin must be 'freq', 'global.freq', or a numeric array")
    else:
        try:
            first_probabilities = np.asarray(begin, dtype=float).ravel(order="F")
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"begin must be 'freq', 'global.freq', or a numeric array of length {nstates}"
            ) from exc
        if first_probabilities.size != nstates:
            raise ValueError(f"begin numeric array must have length {nstates}")

    if isinstance(prob, str):
        if prob not in {"trate", "freq"}:
            raise ValueError("prob must be 'trate', 'freq', or a numeric matrix/array")
        if prob == "trate":
            from .stats_utils import _seqtrate_states, seqtrate

            fitted = np.asarray(
                seqtrate(
                    x,
                    time_varying=time_varying,
                    weighted=weighted,
                    include_missing=include_missing,
                    as_dataframe=False,
                ),
                dtype=float,
            )
            # seqtrate's axes are real states followed by nr; probabilities
            # below are addressed by the object's actual (possibly interleaved)
            # codes. Reindex only fitted axes; supplied arrays retain their API.
            _, fitted_codes = _seqtrate_states(x, None, include_missing)
            active = [
                (axis, code)
                for axis, code in enumerate(fitted_codes)
                if code is not None and code < nstates
            ]
            transitions = np.zeros((nstates, nstates, ntransitions), dtype=float)
            for source_axis, source_code in active:
                for target_axis, target_code in active:
                    transitions[source_code, target_code, :] = (
                        fitted[source_axis, target_axis, :ntransitions]
                        if time_varying
                        else fitted[source_axis, target_axis]
                    )
        elif time_varying:
            from .stats_basic import seqstatd

            frequencies = np.asarray(
                seqstatd(
                    x,
                    weighted=weighted,
                    include_missing=include_missing,
                )["Frequencies"],
                dtype=float,
            )
            transitions = np.zeros((nstates, nstates, ntransitions), dtype=float)
            for time in range(ntransitions):
                transitions[:, :, time] = frequencies[time + 1, :nstates][None, :]
        else:
            counts = np.asarray(seqistatd(x, include_missing=include_missing), dtype=float)
            totals = (counts * weights[:, None]).sum(axis=0)[:nstates]
            denominator = float(totals.sum())
            frequencies = totals / denominator if denominator > 0 else totals
            stationary = np.repeat(frequencies[None, :], nstates, axis=0)
            transitions = np.repeat(stationary[:, :, None], ntransitions, axis=2)
    else:
        try:
            supplied = np.asarray(prob, dtype=float)
        except (TypeError, ValueError) as exc:
            raise ValueError("prob must be 'trate', 'freq', or a numeric matrix/array") from exc
        if supplied.ndim == 2:
            matrix_shape = (nstates, nstates)
            if supplied.shape != matrix_shape:
                raise ValueError(f"probability matrix must have shape {nstates} x {nstates}")
            transitions = np.repeat(supplied[:, :, None], ntransitions, axis=2)
        elif supplied.ndim == 3:
            array_shape = (nstates, nstates, ntransitions)
            if supplied.shape != array_shape:
                raise ValueError(
                    "time-varying probability array must have shape "
                    f"{nstates} x {nstates} x {ntransitions}"
                )
            transitions = supplied
        else:
            raise ValueError(
                "prob must be a two-dimensional matrix or three-dimensional "
                f"array with state width {nstates}"
            )

    out = np.zeros(x.nseq, dtype=float)
    for i, row in enumerate(rows):
        if row.size == 0:
            continue
        initial = first_probabilities[int(row[0])]
        if not np.isnan(initial) and initial > 0:
            out[i] -= np.log(initial)
        for time in range(1, row.size):
            probability = transitions[int(row[time - 1]), int(row[time]), time - 1]
            if not np.isnan(probability) and probability > 0:
                out[i] -= np.log(probability)
    return out


@_require_stslist
def seqmaintokens(
    x: STSList,
    *,
    k: int = 8,
    mint: float | None = None,
    weighted: bool = True,
    include_missing: bool = False,
) -> np.ndarray:
    """Main tokens (states) by mean time spent (TraMineR seqmaintokens), 0-based.

    Returned indices address ``alphabet(x, with_missing=include_missing)`` — the
    real states, plus the missing marker ``nr`` when ``include_missing=True`` —
    exactly as R's 1-based indices address ``statl``.

    Note
    ----
    Faithful to TraMineR 2.2-14. When ``mint`` is given **and** more than ``k``
    states pass the threshold, R (``seqmaintokens.R``) returns positions *into the
    thresholded subset* rather than alphabet indices — i.e. ``sort(order(meant[which
    (meant>=mint)], decreasing=TRUE)[1:k])`` — which are not the actual top-``k``
    state indices. This behaviour (arguably an R quirk) is reproduced exactly for
    parity; pinned by ``seqmaintokens_J_mint_k2``.
    """
    x = x._canonicalized()
    if not k >= 1:
        raise ValueError("k must be a strictly positive integer")
    from .stats_p1 import seqmeant

    # R selects over `seqmeant`, whose rows are `statl`. traminex's `seqmeant` is
    # marker-padded to the internal codebook width (a recorded DIVERGENCE), so the
    # `nr`/`void` marker rows must be dropped here or a `k` above the real
    # alphabet size would return marker codes as "main tokens". The phantom `nr`
    # row that R appends even for missing-free data contributes a zero mean.
    padded = np.asarray(seqmeant(x, weighted=weighted, include_missing=include_missing)["Mean"])
    meant = np.array(
        [
            0.0 if code is None else float(padded[code])
            for code in _alphabet_codes(x, include_missing)
        ],
        dtype=float,
    )
    if mint is not None:
        main = np.where(meant >= mint)[0]
        if main.size > k:
            # R: order() indices into meant[main] (1-based subset positions), not the
            # alphabet indices -> reproduce as 0-based subset positions.
            order = np.argsort(-meant[main], kind="stable")[:k]
            main = np.sort(order)
    else:
        main = np.sort(np.argsort(-meant, kind="stable")[:k])
    return main.astype(int)


def _spell_integr_weights(x: STSList, *, pow: float, include_missing: bool) -> list[np.ndarray]:
    """Per-spell integration weights (R seqidegrad.private spell.integr branch).

    R builds an SPS sequence whose state at each spell is the spell index repeated
    ``duration`` times, then takes ``seqintegr(., pow=pow)``. With no missing in that
    SPS sequence this reduces to: for spell ``j`` occupying the contiguous block of
    positions ``[cum+1 .. cum+dur_j]`` in ``1..Ltot`` (``Ltot`` = total spell
    duration), the weight is ``sum(block**pow) / sum((1:Ltot)**pow)``. Computed
    directly from ``seqdur`` so it stays aligned with the DSS spells used below.
    """
    durs = seqdur(x, include_missing=include_missing)
    out: list[np.ndarray] = []
    for i in range(x.nseq):
        di = durs[i][~np.isnan(durs[i])].astype(int)
        ltot = int(di.sum())
        if ltot == 0:
            out.append(np.zeros(0, dtype=float))
            continue
        wpos = np.arange(1, ltot + 1, dtype=float) ** float(pow)
        base = wpos.sum()
        weights = np.empty(di.size, dtype=float)
        cum = 0
        for k, d in enumerate(di):
            weights[k] = wpos[cum : cum + d].sum() / base
            cum += d
        out.append(weights)
    return out


def _canonicalize_ranked_sequence(
    x: STSList,
    *,
    contract: StateOrderContract,
) -> tuple[STSList, list[int], list[list[int]]]:
    """Apply R's pre-transition incomparable/equivalent-state rewrites."""
    labels = list(contract.labels)
    transformed = x.copy()
    label_to_code = {state: code for code, state in enumerate(transformed.alphabet)}
    noncomp_positions = np.flatnonzero(contract.code_to_group < 0).astype(int).tolist()
    noncomp_labels = [labels[position] for position in noncomp_positions]
    noncomp_codes = {label_to_code[state] for state in noncomp_labels if state in label_to_code}

    if noncomp_codes:
        # R repeatedly copies each non-initial incomparable cell from its left
        # neighbour.  A left-to-right pass reaches the same fixed point while
        # deliberately leaving an initial incomparable state untouched.
        for column in range(1, transformed.seqlen):
            mask = np.isin(transformed.data[:, column], list(noncomp_codes))
            transformed.data[mask, column] = transformed.data[mask, column - 1]

    equiv_positions: list[list[int]] = []
    for canonical_group in contract.groups:
        positions = list(canonical_group)
        group = [labels[position] for position in positions]
        equiv_positions.append(positions)
        representative = label_to_code.get(group[0]) if group else None
        if group and representative is None and group[0] == transformed.nr:
            # R's stslist factor carries an unobserved nr level.  Materialize the
            # analogous code only in this private working copy when an
            # equivalence class selects that phantom state as its representative.
            representative = len(transformed.alphabet)
            transformed.alphabet.append(group[0])
            if transformed.labels is not None:
                transformed.labels.append(group[0])
            if transformed.cpal is not None:
                transformed.cpal.append(transformed.missing_color)
            transformed.nr_code = representative
            label_to_code[group[0]] = representative
        members = [label_to_code[state] for state in group if state in label_to_code]
        if representative is not None and members:
            transformed.data[np.isin(transformed.data, members)] = representative

    return transformed, noncomp_positions, equiv_positions


def _seqidegrad_correction(
    x: STSList,
    contract: StateOrderContract,
    *,
    method: str = "TRATEDSS",
    weight_type: str = "ADD",
    penalized: str | bool = "BOTH",
    include_missing: bool = False,
    border_effect: float = 10.0,
    spell_integr: bool = False,
    pow: float = 1.0,
) -> np.ndarray:
    """TraMineR seqidegrad.private correction term.

    Partial orders forward-fill non-initial incomparable states before DSS
    construction, and equivalence classes collapse to their first state.  Penalty
    signs and ignored rows/columns then follow the resolved TraMineR state order.

    With ``spell_integr=True`` (the public ``seqidegrad`` / type-2 ``seqinsecurity``
    default) each transition ``j-1 -> j`` is weighted by the integration potential of
    the destination spell ``j`` (``pow`` exponent), and the result is the raw
    penalized sum ``transpen`` (R: ``prop.transpen <- transpen``); with
    ``spell_integr=False`` it is the normalized ``transpen / transw``."""
    labels = list(contract.labels)
    costs = np.asarray(contract.precedence, dtype=float).reshape(-1)
    if costs.size != len(labels):
        raise ValueError("length(stprec) should equal length of alphabet")
    # seqidegrad.R:74-119: a logical penalized means BOTH/NO, "NO" returns before
    # any other check, then method, the rate-method options and penalized are
    # validated in this order.
    if isinstance(penalized, (bool, np.bool_)):
        penalized = "BOTH" if penalized else "NO"
    if penalized == "NO":
        return np.zeros(x.nseq, dtype=float)
    method_names = (
        "FREQ",
        "TRATE",
        "TRATEDSS",
        "RANK",
        "FREQ+",
        "TRATE+",
        "TRATEDSS+",
        "RANK+",
        "ONE",
    )
    if not isinstance(method, str) or method not in method_names:
        raise ValueError(f"'method' must be one of: {', '.join(method_names)}")
    use_mean_tr = method.endswith("+")
    method = method.rstrip("+")
    if method in ("FREQ", "TRATE", "TRATEDSS"):
        if not isinstance(weight_type, str) or weight_type not in ("ADD", "INV", "LOGINV"):
            raise ValueError("'weight_type' must be one of: ADD, INV, LOGINV")
        if not border_effect > 1:
            raise ValueError("border_effect should be strictly greater than one")
    if not isinstance(penalized, str) or penalized not in ("NEG", "POS", "BOTH"):
        raise ValueError("'penalized' must be one of: NEG, POS, BOTH, NO")

    transformed, noncomp_positions, equiv_positions = _canonicalize_ranked_sequence(
        x,
        contract=contract,
    )
    targets = _alph_targets(transformed, include_missing)
    K = len(labels)
    code_to_pos = {
        code: labels.index(transformed.alphabet[code])
        for code in targets
        if transformed.alphabet[code] in labels
    }
    dss = seqdss(transformed, include_missing=include_missing)

    # R's FREQ/TRATE/TRATEDSS go through seqtrate(), whose default is weighted=TRUE,
    # so each transition is counted with its sequence's weight. The i-th DSS / source
    # sequence is the i-th original sequence, so its weight applies in every loop.
    w = (
        np.ones(transformed.nseq, dtype=float)
        if transformed.weights is None
        else np.asarray(transformed.weights, dtype=float)
    )
    if method == "RANK":
        tr = np.abs(costs[:, None] - costs[None, :])
    elif method == "ONE":
        # seqidegrad.R:130-131: every transition weighs one (the diagonal is
        # reset below, as R does again at line 307).
        tr = np.ones((K, K), dtype=float)
    else:
        if method == "FREQ":
            tr = np.zeros((K, K), dtype=float)
            for i in range(dss.nseq):
                row = _strip_void(dss.data[i], dss.void_code)
                codes = [code_to_pos[c] for c in row if c in code_to_pos]
                for j in range(1, len(codes)):
                    tr[codes[j - 1], codes[j]] += w[i]
            s = tr.sum()
            tr = tr / s if s > 0 else tr
        else:
            # TRATEDSS / TRATE: row-normalized transition counts over the rows
            # of the DSS (TRATEDSS) or the original sequence (TRATE), in
            # alphabet(with.missing) order. As in R's seqtrate, only adjacent
            # cells form a pair, and a pair touching void (or, without
            # with.missing, a missing cell) is skipped rather than bridged; this
            # keeps void out of both the pairs and the row-sum denominator.
            src = dss if method == "TRATEDSS" else transformed
            counts = np.zeros((K, K), dtype=float)
            for i in range(src.nseq):
                row = src.data[i]
                for j in range(1, len(row)):
                    before, after = int(row[j - 1]), int(row[j])
                    if before in code_to_pos and after in code_to_pos:
                        counts[code_to_pos[before], code_to_pos[after]] += w[i]
            row_sums = counts.sum(axis=1, keepdims=True)
            tr = np.divide(counts, row_sums, out=np.zeros_like(counts), where=row_sums > 0)
        np.fill_diagonal(tr, 0.0)
        if np.any(tr > 1.0 - 0.1 / border_effect):
            tr = tr - tr / border_effect
        if weight_type == "ADD":
            tr = 1.0 - tr
        elif weight_type == "INV":
            eps = 1e-10
            tr = (1.0 + eps) / (tr + eps)
        elif weight_type == "LOGINV":
            eps = 1e-10
            tr = np.log((1.0 + eps) / (tr + eps))
        diag = np.diag(tr).copy()
        diag[diag == 0] = 1.0
        tr = tr / diag[:, None]

    # Penalty signs follow the resolved ranked states, followed by incomparable
    # states.  The latter rows/columns are zeroed below, matching R.
    rank = [
        int(group) if group >= 0 else len(contract.groups) + position
        for position, group in enumerate(contract.code_to_group)
    ]
    signs = np.zeros((K, K), dtype=float)
    for col_a in range(K):
        for col_b in range(K):
            if rank[col_a] < rank[col_b] and penalized in ("NEG", "BOTH"):
                signs[col_a, col_b] = 1.0
            elif rank[col_a] > rank[col_b] and penalized in ("POS", "BOTH"):
                signs[col_a, col_b] = -1.0
    for positions in equiv_positions:
        if not positions:
            continue
        ix = np.ix_(positions, positions)
        signs[ix] = 0.0
        tr[ix] = 0.0
        for position in positions[1:]:
            tr[:, position] = 0.0
            tr[position, :] = 0.0
    for position in noncomp_positions:
        signs[:, position] = 0.0
        signs[position, :] = 0.0
        tr[:, position] = 0.0
        tr[position, :] = 0.0
    np.fill_diagonal(tr, 0.0)

    integr = (
        _spell_integr_weights(transformed, pow=pow, include_missing=include_missing)
        if spell_integr
        else None
    )
    out = np.zeros(transformed.nseq, dtype=float)
    for i in range(transformed.nseq):
        row = _strip_void(dss.data[i], dss.void_code)
        codes = [code_to_pos[c] for c in row if c in code_to_pos]
        wi = integr[i] if integr is not None else None
        transw = transpen = 0.0
        for j in range(1, len(codes)):
            a, b = codes[j - 1], codes[j]
            wj = float(wi[j]) if wi is not None else 1.0
            transw += tr[a, b] * wj
            transpen += tr[a, b] * signs[a, b] * wj
        if spell_integr:
            out[i] = transpen if transw > 0 else 0.0
        else:
            out[i] = (transpen / transw) if transw > 0 else 0.0
        if use_mean_tr:
            # seqidegrad.R:334-337: scale by the mean transition weight per
            # spell; R's 0/0 for a sequence without spells yields NaN.
            out[i] *= transw / row.size if row.size else np.nan
    return out


def _coerce_precarity_correction(
    correction: float | Sequence[float] | np.ndarray,
    *,
    nseq: int,
) -> np.ndarray:
    """Normalize an explicit R-style correction to one value per sequence."""
    values = np.asarray(correction, dtype=float)
    if values.ndim == 0 or values.size == 1:
        return np.full(nseq, float(values.reshape(-1)[0]), dtype=float)
    values = values.reshape(-1)
    if values.size != nseq:
        raise ValueError(f"correction must have length 1 or {nseq}")
    return values


def _resolve_correction_contract(
    x: STSList,
    contract: StateOrderContract,
    *,
    state_equiv: list[list[str]] | None,
    with_missing: bool,
) -> StateOrderContract:
    """Apply seqidegrad.private's second seqprecstart normalization."""
    return _resolve_state_order(
        x,
        state_equiv=state_equiv,
        stprec=contract.precedence,
        with_missing=with_missing,
    )


def _seqindic_seqprecarity(
    x: STSList,
    *,
    correction: float | Sequence[float] | np.ndarray | None = None,
    otto: float = 0.2,
    a: float = 1.0,
    b: float = 1.2,
    method: str = "TRATEDSS",
    state_order: list[str] | None = None,
    state_equiv: list[list[str]] | None = None,
    stprec: list[float] | None = None,
    include_missing: bool = False,
    penalized: str | bool = "BOTH",
    weight_type: str = "ADD",
    pow: float = 1.0,
    border_effect: float = 10.0,
) -> np.ndarray:
    """R-faithful private precarity path used by ranked ``seqindic``.

    ``state_order`` may be full or partial, and ``state_equiv`` supplies R-style
    equivalence classes for both start costs and the degradation correction.
    """
    contract = _resolve_state_order(
        x,
        state_order=state_order,
        state_equiv=state_equiv,
        stprec=stprec,
        with_missing=include_missing,
    )
    sp = contract.precedence
    labels = list(sequence_alphabet(x, with_missing=include_missing))
    code_to_pos = {
        code: labels.index(x.alphabet[code])
        for code in _alph_targets(x, include_missing)
        if x.alphabet[code] in labels
    }
    sdss = seqdss(x, include_missing=include_missing)
    from .stats_p2 import seqici

    ici = np.asarray(seqici(x, include_missing=include_missing), float)
    if correction is None:
        # seqprecarity.private passes its first seqprecstart result into
        # seqidegrad.private, which resolves explicit costs a second time.  That
        # turns equal non-comparable costs into R equivalence classes.
        correction_contract = _resolve_correction_contract(
            x,
            contract,
            state_equiv=state_equiv,
            with_missing=include_missing,
        )
        correction_values = _seqidegrad_correction(
            x,
            correction_contract,
            method=method,
            weight_type=weight_type,
            penalized=penalized,
            include_missing=include_missing,
            border_effect=border_effect,
            spell_integr=False,
            pow=float(pow),
        )
    else:
        correction_values = _coerce_precarity_correction(correction, nseq=x.nseq)
    out = np.zeros(x.nseq, dtype=float)
    for i in range(x.nseq):
        row = _strip_void(sdss.data[i], sdss.void_code)
        first = code_to_pos.get(int(row[0]), None) if row.size else None
        # empty / all-void row: R matches the first state into the alphabet -> NA,
        # so stprec[NA] = NA propagates -> the whole index is NaN.
        sval = sp[first] if first is not None else np.nan
        out[i] = otto * sval + (1.0 - otto) * (ici[i] ** a) * ((1.0 + correction_values[i]) ** b)
    return out


def _seqindic_seqidegrad(
    x: STSList,
    *,
    state_order: list[str] | None = None,
    state_equiv: list[list[str]] | None = None,
    stprec: list[float] | None = None,
    with_missing: bool = False,
    penalized: str | bool = "BOTH",
    method: str = "RANK",
    weight_type: str = "ADD",
    pow: float = 1.0,
    spell_integr: bool = True,
    border_effect: float = 10.0,
) -> np.ndarray:
    """R-faithful private degradation path used by ranked ``seqindic``.

    Penalized count of the (spell-integration-weighted) downward transitions between
    successive states, ranked by their start cost (``method="RANK"``) or transition
    probability (``method`` in ``FREQ``/``TRATE``/``TRATEDSS``). With ``spell_integr``
    (R's default, since the R ``pow`` argument is numeric there) each transition is
    weighted by the destination spell's integration potential and the result is the
    raw penalized sum; otherwise it is the mean penalty per transition.

    Partial orders and equivalence classes follow TraMineR's incomparable-state
    forwarding and transition-collapse rules. Provenance: TraMineR 2.2-14
    ``seqidegrad.R``.
    """
    contract = _resolve_state_order(
        x,
        state_order=state_order,
        state_equiv=state_equiv,
        stprec=stprec,
        with_missing=with_missing,
    )
    return _seqidegrad_correction(
        x,
        contract,
        method=method,
        weight_type=weight_type,
        penalized=penalized,
        include_missing=with_missing,
        border_effect=border_effect,
        spell_integr=spell_integr,
        pow=float(pow),
    )


def _seqindic_seqinsecurity(
    x: STSList,
    *,
    correction: float | Sequence[float] | np.ndarray | None = None,
    state_order: list[str] | None = None,
    state_equiv: list[list[str]] | None = None,
    stprec: list[float] | None = None,
    with_missing: bool = False,
    pow: float = 1.0,
    spow: float | None = None,
    bound: bool = False,
    method: str = "RANK",
    penalized: str | bool = "BOTH",
    weight_type: str = "ADD",
    border_effect: float = 10.0,
) -> np.ndarray:
    """R-faithful private insecurity path used by ranked ``seqindic``.

    The type-2 branch of ``seqprecarity.private``: with start- and spell-integration
    both on, ``insec = stprec[first] * integr1 + (ici + correction)`` where ``integr1``
    is the integration potential of the first spell (``spow`` exponent), ``ici`` is the
    complexity index, and ``correction`` is the spell-integration-weighted degradation
    term (``seqidegrad`` with ``method`` -- ``"RANK"`` by default). With ``bound=True``
    the index is clamped to ``[min, max]`` of the start cost over the states the
    sequence actually visits.

    Partial orders and equivalence classes share the TraMineR ranked-state model.
    Provenance: TraMineR 2.2-14 ``seqprecarity.R`` (``seqprecarity.private`` type=2).
    """
    if spow is None:
        spow = pow
    contract = _resolve_state_order(
        x,
        state_order=state_order,
        state_equiv=state_equiv,
        stprec=stprec,
        with_missing=with_missing,
    )
    sp = contract.precedence
    labels = list(sequence_alphabet(x, with_missing=with_missing))
    code_to_pos = {
        code: labels.index(x.alphabet[code])
        for code in _alph_targets(x, with_missing)
        if x.alphabet[code] in labels
    }
    sdss = seqdss(x, include_missing=with_missing)
    from .stats_p2 import seqici

    ici = np.asarray(seqici(x, include_missing=with_missing), float)
    if correction is None:
        correction_contract = _resolve_correction_contract(
            x,
            contract,
            state_equiv=state_equiv,
            with_missing=with_missing,
        )
        correction_values = _seqidegrad_correction(
            x,
            correction_contract,
            method=method,
            weight_type=weight_type,
            penalized=penalized,
            include_missing=with_missing,
            border_effect=border_effect,
            spell_integr=True,
            pow=float(pow),
        )
    else:
        correction_values = _coerce_precarity_correction(correction, nseq=x.nseq)
    integr_w = _spell_integr_weights(x, pow=float(spow), include_missing=with_missing)
    out = np.zeros(x.nseq, dtype=float)
    for i in range(x.nseq):
        row = _strip_void(sdss.data[i], sdss.void_code)
        first = code_to_pos.get(int(row[0]), None) if row.size else None
        sval = sp[first] if first is not None else np.nan
        integr1 = float(integr_w[i][0]) if integr_w[i].size else np.nan
        prec = sval * integr1 + (ici[i] + correction_values[i])
        if bound and row.size:
            observed = [code_to_pos[int(c)] for c in row if int(c) in code_to_pos]
            if observed:
                lo = min(sp[p] for p in observed)
                hi = max(sp[p] for p in observed)
                prec = min(max(prec, lo), hi)
        out[i] = prec
    return out


@_require_stslist
def seqprecarity(
    x: STSList,
    *,
    correction: float | Sequence[float] | np.ndarray | None = None,
    otto: float = 0.2,
    a: float = 1.0,
    b: float = 1.2,
    method: str = "TRATEDSS",
    state_order: list[str] | None = None,
    state_equiv: list[list[str]] | None = None,
    stprec: list[float] | None = None,
    include_missing: bool = False,
    penalized: str | bool = "BOTH",
    weight_type: str = "ADD",
    pow: float = 1.0,
    border_effect: float = 10.0,
) -> np.ndarray:
    """Sequence precarity index, type=1 (TraMineR seqprecarity).

    Partial orders and equivalence classes use the same immutable state contract
    as ``seqibad``, ``seqidegrad``, and ranked ``seqindic``. ``method``,
    ``penalized``, ``weight_type`` and ``border_effect`` configure the
    degradation correction exactly as in :func:`seqidegrad`.
    """
    x = x._canonicalized()
    targets = _alph_targets(x, include_missing)
    contract = _resolve_state_order(
        x,
        state_order=state_order,
        state_equiv=state_equiv,
        stprec=stprec,
        with_missing=include_missing,
    )
    sp = contract.precedence
    code_to_pos = {c: p for p, c in enumerate(targets)}
    sdss = seqdss(x, include_missing=include_missing)
    from .stats_p2 import seqici

    ici = np.asarray(seqici(x, include_missing=include_missing), float)
    if correction is None:
        correction_contract = _resolve_correction_contract(
            x,
            contract,
            state_equiv=state_equiv,
            with_missing=include_missing,
        )
        correction_values = _seqidegrad_correction(
            x,
            correction_contract,
            method=method,
            weight_type=weight_type,
            penalized=penalized,
            include_missing=include_missing,
            border_effect=border_effect,
            spell_integr=False,
            pow=float(pow),
        )
    else:
        correction_values = _coerce_precarity_correction(correction, nseq=x.nseq)
    out = np.zeros(x.nseq, dtype=float)
    for i in range(x.nseq):
        row = _strip_void(sdss.data[i], sdss.void_code)
        first = code_to_pos.get(int(row[0]), None) if row.size else None
        sval = sp[first] if first is not None else np.nan
        out[i] = otto * sval + (1.0 - otto) * (ici[i] ** a) * ((1.0 + correction_values[i]) ** b)
    return out


def _r_logical_pow(pow: float | bool | np.bool_, spell_integr: bool | None) -> tuple[float, bool]:
    """Resolve ``pow`` and spell integration as R's public ``seqidegrad`` does.

    A logical ``pow`` is R's spell-integration switch and is then used as the
    number 0 or 1 (``seqidegrad.R:10-15``); a numeric ``pow`` keeps spell
    integration on. An explicit ``spell_integr`` (a traminex argument with no
    R counterpart) overrides that rule.
    """
    if isinstance(pow, (bool, np.bool_)):
        rule = bool(pow)
    else:
        rule = True
    return float(pow), rule if spell_integr is None else bool(spell_integr)


@_require_stslist
def seqidegrad(
    x: STSList,
    *,
    state_order: list[str] | None = None,
    state_equiv: list[list[str]] | None = None,
    stprec: list[float] | None = None,
    with_missing: bool = False,
    penalized: str | bool = "BOTH",
    method: str = "RANK",
    weight_type: str = "ADD",
    pow: float | bool = 1.0,
    spell_integr: bool | None = None,
    border_effect: float = 10.0,
) -> np.ndarray:
    """Degradation index per sequence (TraMineR ``seqidegrad``).

    ``method`` weighs each transition by the start-cost difference
    (``"RANK"``), the transition frequency (``"FREQ"``), the transition rate in
    the sequences (``"TRATE"``) or in their distinct-state sequences
    (``"TRATEDSS"``), or by one (``"ONE"``); the ``"+"`` variants of the first
    four multiply the result by the sequence's mean transition weight per
    spell. ``penalized`` is ``"NEG"``, ``"POS"``, ``"BOTH"`` or ``"NO"``, and a
    logical is read as R reads it (``True`` is ``"BOTH"``, ``False`` is
    ``"NO"``, which returns zeros). ``weight_type`` (``"ADD"``, ``"INV"``,
    ``"LOGINV"``) and ``border_effect`` (strictly above one) are validated for
    the rate methods only, as in R; other values raise ``ValueError``.

    As in R, a logical ``pow`` switches spell integration (``pow=False``
    computes the index without it) and is then used as 0 or 1; a numeric
    ``pow`` keeps spell integration. ``spell_integr`` (no R counterpart)
    overrides that rule when given.
    """
    pow_value, spell_integration = _r_logical_pow(pow, spell_integr)
    x = x._canonicalized()
    contract = _resolve_state_order(
        x,
        state_order=state_order,
        state_equiv=state_equiv,
        stprec=stprec,
        with_missing=with_missing,
    )
    return _seqidegrad_correction(
        x,
        contract,
        method=method,
        weight_type=weight_type,
        penalized=penalized,
        include_missing=with_missing,
        border_effect=border_effect,
        spell_integr=spell_integration,
        pow=pow_value,
    )


@_require_stslist
def seqinsecurity(
    x: STSList,
    *,
    correction: float | Sequence[float] | np.ndarray | None = None,
    state_order: list[str] | None = None,
    state_equiv: list[list[str]] | None = None,
    stprec: list[float] | None = None,
    with_missing: bool = False,
    pow: float = 1.0,
    spow: float | None = None,
    bound: bool = False,
    method: str = "RANK",
    penalized: str | bool = "BOTH",
    weight_type: str = "ADD",
    border_effect: float = 10.0,
) -> np.ndarray:
    """Sequence insecurity index per sequence (TraMineR ``seqinsecurity``).

    ``method``, ``penalized``, ``weight_type`` and ``border_effect`` configure
    the degradation correction exactly as in :func:`seqidegrad`.
    """
    x = x._canonicalized()
    if spow is None:
        spow = pow
    targets = _alph_targets(x, with_missing)
    contract = _resolve_state_order(
        x,
        state_order=state_order,
        state_equiv=state_equiv,
        stprec=stprec,
        with_missing=with_missing,
    )
    sp = contract.precedence
    code_to_pos = {c: p for p, c in enumerate(targets)}
    sdss = seqdss(x, include_missing=with_missing)
    from .stats_p2 import seqici

    ici = np.asarray(seqici(x, include_missing=with_missing), float)
    if correction is None:
        correction_contract = _resolve_correction_contract(
            x,
            contract,
            state_equiv=state_equiv,
            with_missing=with_missing,
        )
        correction_values = _seqidegrad_correction(
            x,
            correction_contract,
            method=method,
            weight_type=weight_type,
            penalized=penalized,
            include_missing=with_missing,
            border_effect=border_effect,
            spell_integr=True,
            pow=float(pow),
        )
    else:
        correction_values = _coerce_precarity_correction(correction, nseq=x.nseq)
    integr_w = _spell_integr_weights(x, pow=float(spow), include_missing=with_missing)
    out = np.zeros(x.nseq, dtype=float)
    for i in range(x.nseq):
        row = _strip_void(sdss.data[i], sdss.void_code)
        first = code_to_pos.get(int(row[0]), None) if row.size else None
        sval = sp[first] if first is not None else np.nan
        integr1 = float(integr_w[i][0]) if integr_w[i].size else np.nan
        prec = sval * integr1 + (ici[i] + correction_values[i])
        if bound and row.size:
            observed = [code_to_pos[int(c)] for c in row if int(c) in code_to_pos]
            if observed:
                lo = min(sp[p] for p in observed)
                hi = max(sp[p] for p in observed)
                prec = min(max(prec, lo), hi)
        out[i] = prec
    return out


def _alph_targets(x: STSList, include_missing: bool) -> list[int]:
    """Column states for the with.missing alphabet: real states, then nr if asked
    (void never included). Preserves alphabet order."""
    miss = {c for c in (x.nr_code, x.void_code) if c is not None}
    targets = [i for i in range(len(x.alphabet)) if i not in miss]
    if include_missing and x.nr_code is not None:
        targets.append(x.nr_code)
    return targets


@_require_stslist
def seqintegr(
    x: STSList, *, state: str | None = None, pow: float = 1.0, include_missing: bool = False
) -> np.ndarray:
    """Integrative potential / position-weighted state presence (TraMineR seqintegr).

    ``state``, when given, must be a single state of
    ``alphabet(x, with_missing=include_missing)``, as in R: the missing state is
    accepted only with ``include_missing=True`` (a zero column when no cell is
    missing) and the void marker never.
    """
    x = x._canonicalized()
    L = x.seqlen
    wpos = (np.arange(1, L + 1, dtype=float)) ** float(pow)
    base = float(wpos.sum())
    # The denominator sums the weights of the sequence's own positions: void cells are
    # never part of it, nonresponse only with include_missing=True. R counts the void
    # padding too when with.missing=TRUE (seqintegr.R:30-34), against seqintegr.Rd's
    # positions 1..L, so its value depended on the data's width (TX-AUD-20261019-012).
    miss = {c for c in (x.nr_code, x.void_code) if c is not None}
    excluded = {x.void_code} if include_missing else {x.nr_code, x.void_code}
    excluded_codes = {c for c in excluded if c is not None}
    suminteg = np.full(x.nseq, base, dtype=float)
    if excluded_codes:
        excluded_arr = np.array(sorted(excluded_codes), dtype=int)
        for i in range(x.nseq):
            suminteg[i] = base - float(wpos[np.isin(x.data[i], excluded_arr)].sum())
    if state is not None:
        # R validates against alphabet(seqdata, with.missing), never the
        # internal codebook: markers are refused unless the missing state is
        # requested, which R accepts even when no cell is missing
        # (seqintegr.R:9-20; audit S3-17).
        if isinstance(state, (list, tuple, np.ndarray)):
            if len(state) != 1:
                raise ValueError("When non null, 'state' must be a single state")
            state = state[0]
        if state not in sequence_alphabet(x, with_missing=include_missing):
            raise ValueError(f"state {state} not in the alphabet!")
        if include_missing and state == x.nr and (x.nr_code is not None or state not in x.alphabet):
            # -1 matches no cell: a missing state absent from the data scores 0.
            targets = [x.nr_code if x.nr_code is not None else -1]
        else:
            targets = [list(x.alphabet).index(state)]
    else:
        targets = _alph_targets(x, include_missing)
    out = np.zeros((x.nseq, len(targets)), dtype=float)
    for i in range(x.nseq):
        denom = suminteg[i]
        for j, code in enumerate(targets):
            num = float(wpos[x.data[i] == code].sum())
            # empty / all-void row (denom==0) -> R divides by zero -> NaN
            out[i, j] = (num / denom) if denom > 0 else np.nan
    if state is not None:
        return out.ravel()
    # Phantom nr column(s): include_missing requested but nr absent from the data ->
    # R's alphabet(with.missing) still appends an (all-zero) nr column.
    k_real = len([i for i in range(len(x.alphabet)) if i not in miss])
    n_phantom = (k_real + (1 if include_missing else 0)) - len(targets)
    if n_phantom > 0:
        out = np.hstack([out, np.zeros((x.nseq, n_phantom), dtype=float)])
    return out


def seqprecstart(
    x: STSList,
    *,
    state_order: list[str] | None = None,
    state_equiv: list[list[str]] | None = None,
    stprec: list[float] | None = None,
    include_missing: bool = False,
) -> np.ndarray:
    """Per-state precarity start-cost vector over the alphabet (TraMineR seqprecstart).

    The returned vector has length ``K_eff = k_real + (1 if include_missing else 0)``
    — the size of R's ``alphabet(with.missing)`` — so a phantom nr cost (= 1, the top
    of the ramp) is appended when ``include_missing`` is requested but nr is absent
    from the data.
    """
    contract = _resolve_state_order(
        x,
        state_order=state_order,
        state_equiv=state_equiv,
        stprec=stprec,
        with_missing=include_missing,
    )
    return contract.precedence.copy()


@_require_stslist
def seqibad(
    x: STSList,
    *,
    pow: float = 1.0,
    state_order: list[str] | None = None,
    state_equiv: list[list[str]] | None = None,
    stprec: list[float] | None = None,
    include_missing: bool = False,
) -> np.ndarray:
    """Sequence 'badness' = sum over states of stprec[state]*integr[:, state] (TraMineR seqibad)."""
    x = x._canonicalized()
    contract = _resolve_state_order(
        x,
        state_order=state_order,
        state_equiv=state_equiv,
        stprec=stprec,
        with_missing=include_missing,
    )
    integr = seqintegr(x, pow=pow, include_missing=include_missing)
    return np.asarray(integr @ contract.precedence, dtype=float)


def _seqindic_seqibad(
    x: STSList,
    *,
    pow: float = 1.0,
    state_order: list[str] | None = None,
    state_equiv: list[list[str]] | None = None,
    stprec: list[float] | None = None,
    include_missing: bool = False,
) -> np.ndarray:
    """R-faithful private badness path used by ranked ``seqindic``."""
    contract = _resolve_state_order(
        x,
        state_order=state_order,
        state_equiv=state_equiv,
        stprec=stprec,
        with_missing=include_missing,
    )
    integr = seqintegr(x, pow=pow, include_missing=include_missing)
    return np.asarray(integr @ contract.precedence, dtype=float)


@_require_stslist
def seqipos(
    x: STSList,
    *,
    pos_states: str | Sequence[str] | None = None,
    neg_states: str | Sequence[str] | None = None,
    index: str = "share",
    pow: float = 1.0,
    w: float = 0.5,
    dss: bool | None = None,
    include_missing: bool = False,
) -> np.ndarray:
    """Positive/negative spell share, integration, or volatility (TraMineR seqipos).

    R's obsolete ``index="integration"`` is accepted as ``"integr"`` with a
    ``DeprecationWarning``.
    """
    x = x._canonicalized()
    from .stats_p2 import seqrecode

    if index == "integration":
        # seqipos.R:11-14 maps the obsolete name, with a notice (audit S3-17).
        import warnings

        warnings.warn(
            "index='integration' obsolete, use 'integr' instead",
            DeprecationWarning,
            stacklevel=_caller_stacklevel(),
        )
        index = "integr"
    if index not in ("share", "integr", "volatility"):
        raise ValueError("index should be one of share, integr, or volatility")
    if pos_states is None and neg_states is None:
        raise ValueError("at least one of pos_states and neg_states must be given")
    if isinstance(pos_states, str):
        pos_states = [pos_states]
    elif pos_states is not None:
        pos_states = list(pos_states)
    if isinstance(neg_states, str):
        neg_states = [neg_states]
    elif neg_states is not None:
        neg_states = list(neg_states)
    alph = sequence_alphabet(x, with_missing=include_missing)
    if pos_states is not None:
        invalid = [state for state in pos_states if state not in alph]
        if invalid:
            raise ValueError(f"invalid values in pos_states: {', '.join(invalid)}")
    if neg_states is not None:
        invalid = [state for state in neg_states if state not in alph]
        if invalid:
            raise ValueError(f"invalid values in neg_states: {', '.join(invalid)}")
    if pos_states is None:
        pos_states = [s for s in alph if s not in (neg_states or [])]
    if neg_states is None:
        neg_states = [s for s in alph if s not in pos_states]
    if len(pos_states) != len(set(pos_states)):
        raise ValueError("Multiple occurrences of same state in pos_states")
    if len(neg_states) != len(set(neg_states)):
        raise ValueError("Multiple occurrences of same state in neg_states")
    if dss is None:
        dss = index == "share"
    s = seqdss(x, include_missing=include_missing) if dss else x
    sbin = seqrecode(s, {"p": list(pos_states), "n": list(neg_states)}, otherwise=s.void)
    if index == "share":
        pidx = list(sbin.alphabet).index("p") if "p" in sbin.alphabet else None
        nidx = list(sbin.alphabet).index("n") if "n" in sbin.alphabet else None
        pmask = (sbin.data == pidx) if pidx is not None else np.zeros_like(sbin.data, bool)
        nmask = (sbin.data == nidx) if nidx is not None else np.zeros_like(sbin.data, bool)
        npos = pmask.sum(axis=1).astype(float)
        nneg = nmask.sum(axis=1).astype(float)
        denom = npos + nneg
        return np.divide(
            npos,
            denom,
            out=np.full_like(npos, np.nan),
            where=denom > 0,
        )
    if index == "integr":
        return np.asarray(seqintegr(sbin, state="p", pow=pow, include_missing=include_missing))
    return np.asarray(seqivolatility(sbin, w=w, include_missing=include_missing))
