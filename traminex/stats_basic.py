from __future__ import annotations

import warnings
from collections.abc import Sequence
from numbers import Integral
from typing import Any

import numpy as np
import pandas as pd

from ._deprecated import _caller_stacklevel
from .sequences import STSList, _remap_issues, _require_stslist, _validate_sequence_weights


class SeqStatDResult(dict[str, np.ndarray]):
    """Mutable seqstatd mapping with TraMineR result metadata."""

    nbseq: float
    weighted: bool
    norm: bool
    xtlab: list[str]
    xtstep: int
    tick_last: bool
    cpal: list[str] | None

    def __init__(
        self,
        frequencies: np.ndarray,
        valid_states: np.ndarray,
        entropy: np.ndarray,
        *,
        nbseq: float,
        weighted: bool,
        norm: bool,
        xtlab: list[str],
        xtstep: int,
        tick_last: bool,
        cpal: list[str] | None,
    ) -> None:
        super().__init__(
            Frequencies=frequencies,
            ValidStates=valid_states,
            Entropy=entropy,
        )
        self.nbseq = float(nbseq)
        self.weighted = bool(weighted)
        self.norm = bool(norm)
        self.xtlab = list(xtlab)
        self.xtstep = int(xtstep)
        self.tick_last = bool(tick_last)
        self.cpal = None if cpal is None else list(cpal)


def seqlength(
    x: STSList,
    with_missing: bool = True,
    *,
    as_dataframe: bool = False,
) -> np.ndarray | pd.DataFrame:
    """Return each sequence's non-void length (R: ``seqlength``).

    - If as_dataframe=False (default), returns a numpy array (n,).
    - If as_dataframe=True, returns a DataFrame with columns [id, seqlength].
    - Missing states count by default; structural void never counts.
    """
    if not isinstance(x, STSList):
        raise TypeError("x must be a state sequence object")

    included = np.ones(x.data.shape, dtype=bool)
    if x.void_code is not None:
        included &= x.data != x.void_code
    if not with_missing and x.nr_code is not None:
        included &= x.data != x.nr_code
    arr = included.sum(axis=1, dtype=int)
    if not as_dataframe:
        return arr

    ids = getattr(x, "ids", None)
    if ids is None:
        ids = np.arange(1, x.nseq + 1, dtype=int)
    return pd.DataFrame({"id": ids, "seqlength": arr})


@_require_stslist
def seqistatd(x: STSList, *, include_missing: bool = False, prop: bool = False) -> np.ndarray:
    """State durations per sequence (counts or row proportions).

    Returns an array shape (n_sequences, K) with counts of each state.
    Structural void is never a state. Missing values are included only when
    requested. With ``prop=True``, included state counts are divided by their
    row total, preserving TraMineR's ``0 / 0 -> NaN`` for all-void rows.
    """
    # Both marker representations must give one shape: the per-state axis is
    # read from the codebook, which only the canonical form pads with the
    # marker columns the R goldens pin (2026-09-22 audit F-04).
    x = x._canonicalized()
    state_codes = x._codebook_codes(
        include_external_nr=include_missing,
        include_external_void=False,
    )
    K = len(state_codes)
    if x.nseq == 0 or K == 0:
        dtype = float if prop else int
        return np.zeros((x.nseq, K), dtype=dtype)
    # Vectorised: one np.bincount call over the entire flattened matrix.
    # linear_idx maps each element to its (row, col) position in the (n, K) output.
    flat = x.data.ravel()  # shape (n * L,)
    row_idx = np.repeat(np.arange(x.nseq), x.seqlen)  # which sequence each element belongs to
    code_to_column = {code: column for column, code in enumerate(state_codes)}
    columns = np.fromiter(
        (code_to_column.get(int(code), -1) for code in flat),
        dtype=int,
        count=flat.size,
    )
    represented = columns >= 0
    linear_idx = row_idx[represented] * K + columns[represented]
    counts = np.bincount(linear_idx, minlength=x.nseq * K).reshape(x.nseq, K).astype(int)

    excluded = {x.void_code}
    if not include_missing:
        excluded.add(x.nr_code)
    excluded.discard(None)
    excluded_columns = {
        code_to_column[code] for code in excluded if code is not None and code in code_to_column
    }
    for column in excluded_columns:
        counts[:, column] = 0

    if not prop:
        return counts

    proportions = np.zeros((x.nseq, K), dtype=float)
    included = [column for column in range(K) if column not in excluded_columns]
    if included:
        included_counts = counts[:, included]
        totals = included_counts.sum(axis=1, dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            proportions[:, included] = included_counts / totals[:, None]
    return proportions


@_require_stslist
def seqstatd(
    x: STSList,
    *,
    weighted: bool = True,
    include_missing: bool = False,
    norm: bool = True,
    weights: np.ndarray | None = None,
) -> SeqStatDResult:
    """Cross-sectional distributions, valid-state totals, and entropy.

    The mutable three-key mapping preserves the established Python ``L x K``
    frequency orientation. TraMineR result metadata is available as attributes.
    An explicit ``weights`` vector overrides attached weights when
    ``weighted=True``; ``weighted=False`` ignores both. With
    ``include_missing=True`` the missing state always counts in the normalised
    entropy's maximum and in ``cpal``, as in R, even when no cell is missing.
    ``cpal`` is R's palette for statl (+ nr), so it has one colour per R state;
    ``Frequencies`` has one column per codebook entry (``x.alphabet``, markers
    included), so the two lengths differ whenever the codebook differs from
    statl (TX-AUD-20261010-028).
    """
    # Both marker representations must give one shape: the per-state axis is
    # read from the codebook, which only the canonical form pads with the
    # marker columns the R goldens pin (2026-09-22 audit F-04).
    x = x._canonicalized()
    L = x.seqlen
    state_codes = x._codebook_codes(
        include_external_nr=include_missing,
        include_external_void=False,
    )
    K = len(state_codes)
    n = x.nseq
    if not weighted:
        w_seq = np.ones(n, dtype=float)
    elif weights is not None:
        w_seq = np.asarray(weights, dtype=float)
    elif x.weights is not None:
        w_seq = np.asarray(x.weights, dtype=float)
    else:
        w_seq = np.ones(n, dtype=float)
    w_seq = _validate_sequence_weights(w_seq, n)
    result_weighted = bool(weighted and not np.all(w_seq == 1.0))
    freq = np.zeros((L, K), dtype=float)
    valid = np.zeros(L, dtype=float)
    # Columns counted as states = R's statl: the real alphabet, plus nr when
    # with.missing. Structural void is NEVER a state (R never adds void to
    # statl), so it is excluded regardless of include_missing.
    include_mask = np.array(
        [code != x.void_code and (include_missing or code != x.nr_code) for code in state_codes],
        dtype=bool,
    )
    code_to_column = {code: column for column, code in enumerate(state_codes)}
    data_T = x.data.T  # shape (L, n) — iterate over positions in the outer loop
    for t in range(L):
        col = data_T[t]  # shape (n,)
        columns = np.fromiter(
            (code_to_column.get(int(code), -1) for code in col),
            dtype=int,
            count=n,
        )
        represented = columns >= 0
        freq[t] = np.bincount(
            columns[represented],
            weights=w_seq[represented],
            minlength=K,
        )
        s = float(freq[t, include_mask].sum())
        valid[t] = s
        if s > 0:
            freq[t, include_mask] /= s
        else:
            freq[t, include_mask] = np.nan
        freq[t, ~include_mask] = 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        entropy_terms = np.where(freq > 0.0, freq * np.log(freq), 0.0)
    ent = -np.sum(entropy_terms, axis=1)
    ent[valid == 0.0] = np.nan
    # Maximum entropy is the entropy of a uniform distribution over statl, i.e.
    # the number of state columns counted (real alphabet, + nr when with.missing).
    # R's statl gains nr under with.missing even when no cell is missing
    # (seqstatd.R:12-17), so an object without an nr code still counts it.
    K_eff = int(include_mask.sum())
    absent_nr = include_missing and (x.nr_code is None or x.nr_code not in state_codes)
    if absent_nr:
        K_eff += 1
    if norm:
        with np.errstate(divide="ignore", invalid="ignore"):
            ent = ent / np.log(K_eff)

    if x.xtlab is None:
        xtlab = [f"T{x.start + offset}" for offset in range(L)]
    else:
        xtlab = list(x.xtlab)
    if x.cpal is None:
        result_cpal = None
    else:
        result_cpal = [
            (
                x.missing_color
                if code == x.nr_code
                else (x.cpal[code] if code < len(x.cpal) else x.missing_color)
            )
            for code, included in zip(state_codes, include_mask)
            if included
        ]
        if absent_nr:
            result_cpal.append(x.missing_color)
    return SeqStatDResult(
        freq,
        valid,
        ent,
        nbseq=float(w_seq.sum()),
        weighted=result_weighted,
        norm=norm,
        xtlab=xtlab,
        xtstep=x.xtstep,
        tick_last=x.tick_last,
        cpal=result_cpal,
    )


_SEQTAB_DEFAULT_IDXS = tuple(range(1, 11))


class _OmittedSeqtabTlim:
    def __repr__(self) -> str:
        return "None"


_SEQTAB_TLIM_OMITTED: Any = _OmittedSeqtabTlim()


def _subset_seqtab_rows(
    x: STSList,
    row_indices: Sequence[int] | np.ndarray,
    *,
    weights: np.ndarray | None = None,
) -> STSList:
    rows = np.asarray(row_indices, dtype=int).reshape(-1)
    subset_weights = (
        np.asarray(weights, dtype=float).copy()
        if weights is not None
        else (None if x.weights is None else np.asarray(x.weights, dtype=float)[rows].copy())
    )
    return STSList(
        data=x.data[rows].copy(),
        alphabet=list(x.alphabet),
        labels=None if x.labels is None else list(x.labels),
        cpal=None if x.cpal is None else list(x.cpal),
        xtlab=None if x.xtlab is None else list(x.xtlab),
        xtstep=x.xtstep,
        tick_last=x.tick_last,
        nr=x.nr,
        void=x.void,
        nr_code=x.nr_code,
        void_code=x.void_code,
        weights=subset_weights,
        start=x.start,
        ids=None if x.ids is None else [x.ids[int(index)] for index in rows],
        missing_color=x.missing_color,
        issues=_remap_issues(x.issues, rows),
    )


def _seqtab_is_default_idxs(value: object) -> bool:
    if isinstance(value, (str, bytes)) or value is None:
        return False
    try:
        array = np.asarray(value)
    except (TypeError, ValueError):
        return False
    return bool(array.ndim == 1 and array.tolist() == list(_SEQTAB_DEFAULT_IDXS))


def _resolve_seqtab_idxs(value: object, unique_count: int) -> list[int]:
    if isinstance(value, (bool, np.bool_)) or value is None or isinstance(value, (str, bytes)):
        raise TypeError("idxs must be scalar 0 or a nonempty vector of positive integer ranks")

    if isinstance(value, (Integral, np.integer)):
        rank = int(value)
        if rank == 0:
            return list(range(1, unique_count + 1))
        if rank < 0:
            raise ValueError("idxs must be 0 or contain only strictly positive integer ranks")
        requested = [rank]
    else:
        try:
            array = np.asarray(value)
        except (TypeError, ValueError) as error:
            raise TypeError(
                "idxs must be scalar 0 or a nonempty vector of positive integer ranks"
            ) from error
        if array.ndim <= 1 and array.size == 1:
            # R has no scalars: c(0) is 0 and means every rank (seqtab.R:12, 50;
            # audit S3-21), so one value takes the scalar route.
            item = array.reshape(-1).tolist()[0]
            if isinstance(item, (bool, np.bool_)) or not isinstance(item, (Integral, np.integer)):
                raise TypeError("idxs ranks must be integers, not logical or fractional values")
            return _resolve_seqtab_idxs(int(item), unique_count)
        if array.ndim != 1 or array.size == 0:
            raise ValueError("idxs must be a nonempty one-dimensional rank vector")
        requested = []
        for item in array.tolist():
            if isinstance(item, (bool, np.bool_)) or not isinstance(item, (Integral, np.integer)):
                raise TypeError("idxs ranks must be integers, not logical or fractional values")
            rank = int(item)
            if rank <= 0:
                raise ValueError("idxs vectors must contain only strictly positive integer ranks")
            requested.append(rank)

    if unique_count == 0 or max(requested) > unique_count:
        return list(range(1, unique_count + 1))
    return requested


def _slice_seqtab_result(table: pd.DataFrame, rows: slice | Sequence[int]) -> pd.DataFrame:
    """Slice a seqtab frame while keeping representative metadata aligned."""
    positions = np.arange(len(table))[rows]
    positions = np.asarray(positions, dtype=int).reshape(-1)
    result = table.iloc[positions].copy()
    result["cumprop"] = np.cumsum(result["prop"].to_numpy(dtype=float))
    result.attrs["freq"] = table.attrs["freq"].iloc[positions].copy()
    for name in ("idxs", "idxf", "source_indices"):
        values = table.attrs[name]
        result.attrs[name] = [values[int(position)] for position in positions]
    representatives = table.attrs["representatives"]
    rep_weights = (
        None
        if representatives.weights is None
        else np.asarray(representatives.weights, dtype=float)[positions]
    )
    result.attrs["representatives"] = _subset_seqtab_rows(
        representatives,
        positions,
        weights=rep_weights,
    )
    return result


def seqtab(
    x: STSList,
    idxs: int | Sequence[int] = _SEQTAB_DEFAULT_IDXS,
    weighted: bool = True,
    format: str = "SPS",
    tlim: int | Sequence[int] | None = _SEQTAB_TLIM_OMITTED,
    *,
    sep: str = "-",
    weights: np.ndarray | None = None,
    sort: bool = True,
    normalize: bool = True,
    include_missing: bool = True,
) -> pd.DataFrame:
    """Frequency table with TraMineR ranks and a pandas result re-expression.

    The R-compatible defaults return the ten highest-frequency ranks in
    compressed ``SPS`` form. ``idxs=0`` (or ``[0]``, R's ``c(0)``) selects every distinct pattern.
    The DataFrame columns remain the established Python contract:
    ``sequence``, ``count``, ``prop``, and ``cumprop``. R result metadata and
    representative sequence rows are available under :attr:`DataFrame.attrs`:
    ``idxs`` holds R's one-based frequency ranks, ``idxf`` R's one-based row of
    each pattern's first occurrence among the sequences seqtab keeps (empty
    sequences dropped, as R does), and ``source_indices`` (traminex's own) the
    zero-based row of that occurrence in ``x``, for indexing ``x.data``.

    ``sep``, explicit ``weights``, ``sort``, ``normalize``, and
    ``include_missing`` are retained Python extensions; ``include_missing=False``
    leaves out the sequences that hold a missing cell (void padding does not
    count). ``weighted=False`` ignores both attached and explicit weights. Python deliberately rejects
    logical, fractional, empty, and malformed rank values rather than relying
    on R's permissive indexing coercions.
    """
    if not isinstance(x, STSList):
        raise TypeError("x must be an STSList state sequence object")
    if format not in {"SPS", "STS"}:
        raise ValueError("format must be one of: SPS or STS")

    if tlim is not _SEQTAB_TLIM_OMITTED:
        if not _seqtab_is_default_idxs(idxs):
            raise ValueError(
                "idxs and tlim cannot be specified together; tlim is deprecated, use idxs instead"
            )
        warnings.warn(
            "tlim is deprecated, use idxs instead",
            DeprecationWarning,
            stacklevel=_caller_stacklevel(),
        )
        idxs = tlim  # type: ignore[assignment]

    if not weighted:
        base_weights = np.ones(x.nseq, dtype=float)
    elif weights is not None:
        base_weights = np.asarray(weights, dtype=float)
    elif x.weights is not None:
        base_weights = np.asarray(x.weights, dtype=float)
    else:
        base_weights = np.ones(x.nseq, dtype=float)
    base_weights = _validate_sequence_weights(base_weights, x.nseq)

    eligible = np.ones(x.nseq, dtype=bool)
    if x.nr_code is not None:
        eligible &= ~np.all(x.data == int(x.nr_code), axis=1)
    if not include_missing and x.nr_code is not None:
        # Only a missing cell disqualifies a sequence; void is padding, so a
        # shorter sequence is kept (TX-AUD-20261010-027).
        eligible &= ~np.any(x.data == int(x.nr_code), axis=1)
    eligible_indices = np.flatnonzero(eligible)
    eligible_x = _subset_seqtab_rows(x, eligible_indices)
    effective_weights = base_weights[eligible_indices].copy()

    real_state_codes = {code for code in (x.nr_code, x.void_code) if code is not None}
    real_states = [
        str(state) for code, state in enumerate(x.alphabet) if code not in real_state_codes
    ]
    if any("-" in state for state in real_states):
        warnings.warn(
            "'-' character in state codes may cause invalid results",
            UserWarning,
            stacklevel=2,
        )

    if format == "SPS":
        from .seqformat import _sts_to_sps

        sequence_labels = np.asarray(
            _sts_to_sps(
                eligible_x,
                sps_out={"xfix": "", "sdsep": "/"},
                right="DEL",
            ),
            dtype=object,
        )
    else:
        from .stats_p2 import seqconc

        sequence_labels = seqconc(eligible_x, sep=sep)

    # return_index gives each label's first row directly (a stable sort), which
    # replaced a per-label scan that was quadratic in the number of distinct
    # sequences (audit S3-18).
    unique_labels, first_index, inverse = np.unique(
        sequence_labels, return_index=True, return_inverse=True
    )
    counts = np.bincount(
        inverse,
        weights=effective_weights,
        minlength=len(unique_labels),
    ).astype(float)
    if sort and len(unique_labels):
        from ._categorical_order import r_character_order

        label_order = {
            label: position
            for position, label in enumerate(r_character_order(unique_labels.astype(str).tolist()))
        }
        order = np.asarray(
            sorted(
                range(len(unique_labels)),
                key=lambda position: (
                    -counts[position],
                    label_order[str(unique_labels[position])],
                ),
            ),
            dtype=int,
        )
    else:
        order = np.arange(len(unique_labels))
    unique_labels = unique_labels[order]
    counts = counts[order]
    first_index = first_index[order]

    full_total = float(effective_weights.sum())
    if full_total != 0:
        percentages = counts / full_total * 100.0
    else:
        percentages = np.full_like(counts, np.nan)
    resolved_idxs = _resolve_seqtab_idxs(idxs, len(unique_labels))
    selected_positions = np.asarray([rank - 1 for rank in resolved_idxs], dtype=int)
    selected_labels = unique_labels[selected_positions]
    selected_counts = counts[selected_positions]
    selected_percentages = percentages[selected_positions]

    first_positions = np.asarray(first_index, dtype=int)
    representative_positions = first_positions[selected_positions]
    source_indices = eligible_indices[representative_positions]
    representatives = _subset_seqtab_rows(
        x,
        source_indices,
        weights=selected_counts,
    )

    if normalize:
        proportions = selected_percentages / 100.0
        cumulative = np.cumsum(proportions)
    else:
        proportions = np.zeros_like(selected_counts)
        cumulative = np.zeros_like(selected_counts)
    result = pd.DataFrame(
        {
            "sequence": selected_labels,
            "count": selected_counts,
            "prop": proportions,
            "cumprop": cumulative,
        }
    )
    result.attrs.update(
        {
            "freq": pd.DataFrame(
                {"Freq": selected_counts, "Percent": selected_percentages},
                index=pd.Index(selected_labels, name=None),
            ),
            "nbseq": full_total,
            "weighted": bool(not np.all(effective_weights == 1.0)),
            "idxs": list(resolved_idxs),
            # R's idxf: the one-based row of each pattern's first occurrence among the
            # sequences seqtab keeps (seqtab.R:55, 75; TX-AUD-20261019-034).
            "idxf": (representative_positions.astype(int) + 1).tolist(),
            "source_indices": source_indices.astype(int).tolist(),
            "format": format,
            "representatives": representatives,
        }
    )
    return result


@_require_stslist
def seqtransn(
    x: STSList,
    *,
    with_missing: bool = False,
    norm: bool = False,
    pweight: bool = False,
) -> np.ndarray:
    """Number of transitions (state changes) in each sequence (TraMineR ``seqtransn``).

    For each sequence this is ``length(seqdss) - 1``: the number of distinct
    successive states minus one. ``with_missing`` controls whether the missing
    code is a state in the DSS (R ``with.missing``). ``norm=True`` divides by
    ``seqlength(seqdata) - 1`` (the *non-void* length, 0 for length-1 sequences).
    ``pweight=True`` weights each transition by ``1 - p`` where ``p`` is its
    dataset transition rate (``seqtrate``), as in R; since that rate matrix
    has no missing state, a sequence with a transition to or from the missing
    state (``with_missing=True``) gets NaN, as R's ``NA``.

    Returns an array of shape ``(n_sequences,)``. (The previous dataset-level
    K x K transition *matrix* is available privately as
    :func:`_transition_count_matrix`; for transition rates use :func:`seqtrate`.)
    """
    x = x._canonicalized()
    from .stats_complexity import seqdss

    dss = seqdss(x, include_missing=with_missing)
    if dss.void_code is not None:
        dssl = (dss.data != dss.void_code).sum(axis=1)
    else:
        dssl = np.full(dss.nseq, dss.data.shape[1], dtype=int)
    trans = np.maximum(dssl - 1, 0).astype(float)

    if pweight:
        from .stats_utils import seqasnum, seqtrate

        tr = np.asarray(seqtrate(x, as_dataframe=False), dtype=float)
        dss_num = seqasnum(dss)  # 0-based codes over the real alphabet; NaN for void/nr
        trans = np.zeros(x.nseq, dtype=float)
        for i in range(x.nseq):
            for j in range(1, int(dssl[i])):
                a, b = dss_num[i, j - 1], dss_num[i, j]
                if np.isnan(a) or np.isnan(b):
                    # A missing DSS cell (with_missing=True) has no row or
                    # column in seqtrate's matrix; R's tr[NA, .] is NA and
                    # makes the whole count NA (seqtransn.R:16-24).
                    trans[i] = np.nan
                else:
                    trans[i] += 1.0 - tr[int(a), int(b)]

    if norm:
        if x.void_code is not None:
            seql = (x.data != x.void_code).sum(axis=1).astype(float)
        else:
            seql = np.full(x.nseq, x.seqlen, dtype=float)
        denom = seql - 1.0
        trans = np.divide(trans, denom, out=np.zeros_like(trans), where=denom > 0)
    return trans


def _transition_count_matrix(x: STSList, *, normalize: bool = False) -> np.ndarray:
    """Dataset-level transition count matrix (K×K); private (not R ``seqtransn``).

    Aggregates transitions across consecutive positions over all sequences,
    weighted by ``x.weights``. With ``normalize=True`` rows are row-normalized
    (equivalent to :func:`seqtrate` over the full alphabet). Kept for the
    transition heatmap; the public R name ``seqtransn`` is the per-sequence
    transition count. The void padding is not a state: pairs that involve it
    are not transitions, so its row and column are zero.
    """
    x = x._canonicalized()
    K = len(x.alphabet)
    if x.seqlen < 2:
        return np.zeros((K, K), dtype=float)
    # Linear state-pair indexes can exceed the storage dtype even when both
    # state codes fit (for example uint8 code 19 in a twenty-state alphabet).
    a = x.data[:, :-1].astype(np.int64).ravel()
    b = x.data[:, 1:].astype(np.int64).ravel()
    reps = x.seqlen - 1
    w_seq = np.asarray(x.weights, dtype=float) if x.weights is not None else np.ones(x.nseq)
    w_full = np.repeat(w_seq, reps)
    if x.void_code is not None:
        # audit S3-19: moving into (or out of) the padding is not a transition.
        real_pair = (a != x.void_code) & (b != x.void_code)
        a, b, w_full = a[real_pair], b[real_pair], w_full[real_pair]
    idx = a * K + b
    M = np.bincount(idx, weights=w_full, minlength=K * K).astype(float).reshape(K, K)
    if normalize:
        row_sums = M.sum(axis=1, keepdims=True)
        M = np.divide(M, row_sums, out=np.zeros_like(M), where=row_sums > 0)
    return M


@_require_stslist
def entropy(x: STSList, *, normalized: bool = True) -> float:
    """Global Shannon entropy over all positions and sequences.

    Uses optional sequence weights if present. If normalized=True, divide by log(K).
    The void padding of unequal-length sequences is not a state: its cells are
    not counted and it is not one of the ``K`` states. The missing state is: when
    the object's codebook has the missing marker, its cells are counted and it is
    one of the ``K`` states, as ``seqstatd(include_missing=True)`` counts it.
    """
    x = x._canonicalized()
    # audit S3-19: void is padding, never a state.
    K = len(x.alphabet) - (1 if x.void_code is not None else 0)
    if K <= 1:
        return 0.0
    vals = x.data.ravel()
    nseq = x.nseq
    reps = x.seqlen
    if x.weights is not None:
        w_seq = np.asarray(x.weights, dtype=float)
        w_full = np.repeat(w_seq, reps)
    else:
        w_full = np.ones(vals.shape[0], dtype=float)
    if x.void_code is not None:
        state_cell = vals != x.void_code
        vals, w_full = vals[state_cell], w_full[state_cell]
    counts = np.bincount(vals, weights=w_full, minlength=len(x.alphabet)).astype(float)
    total = counts.sum()
    if total <= 0:
        return 0.0
    # A zero frequency contributes 0, as R's entropy.R:5-7 drops it; adding a
    # constant inside the log biased uniform data below 1 (TX-AUD-20261010-095).
    p = counts[counts > 0] / total
    H = float(-np.sum(p * np.log(p)))
    if normalized:
        return H / float(np.log(K))
    return H
