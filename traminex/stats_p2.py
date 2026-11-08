"""Priority 2 parity functions: complexity, generation, recoding."""

from __future__ import annotations

import math
import re
import sys
import unicodedata
import warnings
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Iterable, List, Literal, Mapping, Optional, Sequence, cast

import numpy as np
import pandas as pd

from .seqformat import _cell_is_na, _collision_free_marker, _r_state_text, _table_to_matrix
from .sequences import STSList, _require_stslist
from .stats_basic import seqistatd


def _is_na(v) -> bool:
    if v is None:
        return True
    try:
        return bool(np.isnan(float(v)))
    except (TypeError, ValueError):
        return False


def _implicativestat_weighted(
    x: Sequence,
    y: Sequence,
    weights: np.ndarray,
    *,
    resid: str = "adjusted",
    continuity: bool = True,
):
    """TraMineRextras ``implicativestat`` (weighted, multi-level; used by ``seqimplic``).

    Unlike the ``TraMineR`` internal ``implicativestat`` ported for ``seqerules`` (which
    is hard-wired to 0/1 levels), the TraMineRextras variant keeps the actual factor
    levels of ``x``/``y`` and accepts case ``weights`` and a continuity correction. With
    ``resid="adjusted"`` and ``type="indice"`` (the ``seqimplic`` settings) the cell is
    ``(Nnbj - Nexpnbj) / sqrt(Nexpnbj * (W_condj/n) * (1 - W_condi/n))`` where
    ``Nnbj = W[condi & !condj] (+0.5 if continuity)``. Provenance: TraMineRextras 0.6.9
    ``seqimplic.R``.
    """
    xs = np.array([str(v) for v in x], dtype=object)
    ys = np.array([str(v) for v in y], dtype=object)
    w = np.asarray(weights, dtype=float)
    xgrp = sorted(set(xs.tolist()))
    ygrp = sorted(set(ys.tolist()))
    n = w.sum()
    result = np.zeros((len(xgrp), len(ygrp)), dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        for i, xl in enumerate(xgrp):
            condi = xs == xl
            for j, yl in enumerate(ygrp):
                condj = ys == yl
                nnbj = w[condi & ~condj].sum()
                if continuity:
                    nnbj += 0.5
                nexpnbj = w[condi].sum() * w[~condj].sum() / n
                if resid == "adjusted":
                    indice = (nnbj - nexpnbj) / np.sqrt(
                        nexpnbj * (w[condj].sum() / n) * (1 - w[condi].sum() / n)
                    )
                else:  # standard
                    indice = (nnbj - nexpnbj) / np.sqrt(nexpnbj)
                result[i, j] = indice
    return result, xgrp, ygrp


@dataclass
class SeqImplic:
    """Result of :func:`seqimplic` (R ``seqimplic`` object; the numeric table).

    ``indices`` is the ``(n_levels, n_alphabet, n_positions)`` array of implication
    statistics; cells for states / group levels absent at a position are ``NaN``. The
    plotting metadata mirrors the additive fields carried by TraMineRextras so
    :func:`traminex.plot_seqimplic` can render the result without the source
    sequence object.
    """

    indices: np.ndarray
    levels: List[str]
    alphabet: List[str]
    positions: List
    labels: Optional[List[str]] = None
    type: str = "indice"
    cpal: Optional[List[str]] = None
    xtstep: int = 1
    tick_last: bool = False


def _seqimplic_label(value: object) -> str:
    """R's ``as.character`` of one group value: ``2``, ``1e+05``, ``TRUE`` (audit S3-20)."""
    from ._rformat import r_as_character

    label = r_as_character(value)
    return "NA" if label is None else label


def _seqimplic_levels(group: object, kept: list) -> List[str]:
    """``levels(factor(group))`` for the kept rows, as ``seqimplic``'s str labels.

    A pandas categorical keeps its declared order (unused levels dropped, as
    ``factor()`` of a factor does), numbers sort numerically, logicals put
    False first and anything else follows R's character collation; NA is never
    a level (audit S3-15). Labels are R's ``as.character`` of the values, so
    values with the same label share a level (audit S3-20).
    """
    from numbers import Real

    from ._categorical_order import r_character_order

    observed = [value for value in kept if not _is_na(value)]
    present = {_seqimplic_label(value) for value in observed}
    categories = None
    if isinstance(group, pd.Categorical):
        categories = list(group.categories)
    elif isinstance(group, pd.Series) and isinstance(group.dtype, pd.CategoricalDtype):
        categories = list(group.cat.categories)
    if categories is not None:
        ordered = [_seqimplic_label(value) for value in categories]
    elif observed and all(isinstance(value, (bool, np.bool_)) for value in observed):
        ordered = [_seqimplic_label(value) for value in sorted(set(observed))]
    elif observed and all(
        isinstance(value, Real) and not isinstance(value, (bool, np.bool_)) for value in observed
    ):
        ordered = [_seqimplic_label(value) for value in sorted(set(observed))]
    else:
        ordered = r_character_order([_seqimplic_label(value) for value in observed])
    return [label for label in dict.fromkeys(ordered) if label in present]


@_require_stslist
def seqimplic(
    x: STSList,
    group: Sequence,
    *,
    with_missing: bool = False,
    weighted: bool = True,
    na_rm: bool = True,
) -> "SeqImplic":
    """Implicative statistics of a group membership over sequence positions.

    For each position, computes the adjusted-residual implication statistic (with a
    continuity correction) between the ``group`` factor and the state observed at that
    position, restricted to the non-void (and non-missing, unless ``with_missing``)
    cases. Returns a :class:`SeqImplic` whose ``indices`` array is levels x alphabet x
    positions. Provenance: TraMineRextras 0.6.9 ``seqimplic.R``.

    ``levels`` follow R's ``factor(group)``: a pandas categorical keeps its
    declared order, numbers sort numerically and text follows R's collation;
    labels are R's ``as.character`` (``2`` not ``2.0``, ``TRUE`` not ``True``).
    ``with_missing=True`` always appends the missing state (an all-NaN slice
    when no cell is missing), as R does. ``group`` gives one entry per sequence;
    another length is refused (R stops on a longer group and returns no finite
    value for a shorter one).
    """
    grp = list(group)
    n0 = x.nseq
    if len(grp) != n0:
        # seqimplic.Rd: one membership per sequence. R stops on a longer group and
        # recycles the row mask of a shorter one into a result without a finite
        # value (seqimplic.R:69-72; TX-AUD-20261011-073).
        raise ValueError(
            f"group must have one entry per sequence ({n0} sequences, {len(grp)} entries)"
        )
    idx = [i for i in range(n0) if not _is_na(grp[i])] if na_rm else list(range(n0))
    data = x.data[idx]
    grp_str = [_seqimplic_label(grp[i]) for i in idx]
    grp_na = [_is_na(grp[i]) for i in idx]
    if x.weights is None or not weighted:
        w = np.ones(len(idx), dtype=float)
    else:
        w = np.asarray(x.weights, dtype=float)[idx]

    levels = _seqimplic_levels(group, [grp[i] for i in idx])
    real_codes = [i for i in range(len(x.alphabet)) if i not in (x.void_code, x.nr_code)]
    alph = [x._state_for_code(code) for code in real_codes]
    discard_missing = False
    if with_missing:
        # R appends the missing state even when no cell is missing (audit S3-15).
        alph.append(x.nr if x.nr_code is None else x._state_for_code(x.nr_code))
    elif x.nr_code is not None and np.any(x.data == x.nr_code):
        discard_missing = True

    miss_codes = [c for c in (x.nr_code, x.void_code) if c is not None]
    L = data.shape[1]
    indices = np.full((len(levels), len(alph), L), np.nan, dtype=float)
    lev_pos = {lev: i for i, lev in enumerate(levels)}
    alph_pos = {a: i for i, a in enumerate(alph)}
    for i in range(L):
        col = data[:, i]
        if discard_missing and miss_codes:
            cond = ~np.isin(col, miss_codes)
        elif x.void_code is not None:
            cond = col != x.void_code
        else:
            cond = np.ones(len(col), dtype=bool)
        if not cond.any():
            continue
        if any(grp_na[j] for j in range(len(col)) if cond[j]):
            # na_rm=False: R's sums over an NA group value are NA, so every
            # cell at this position stays NA.
            continue
        xg = [grp_str[j] for j in range(len(col)) if cond[j]]
        yv = [x._state_for_code(int(col[j])) for j in range(len(col)) if cond[j]]
        res, xlev, ylev = _implicativestat_weighted(
            xg, yv, w[cond], resid="adjusted", continuity=True
        )
        for a, xl in enumerate(xlev):
            if xl not in lev_pos:
                continue
            for b, yl in enumerate(ylev):
                if yl in alph_pos:
                    indices[lev_pos[xl], alph_pos[yl], i] = res[a, b]

    if x.labels is None:
        plot_labels = list(alph)
    else:
        plot_labels = [
            x.labels[code] if code < len(x.labels) else x._state_for_code(code)
            for code in real_codes
        ]
        if with_missing:
            plot_labels.append(alph[-1])
    if with_missing:
        plot_labels[-1] = "Missing"

    plot_cpal = None
    if x.cpal is not None:
        plot_cpal = [x.cpal[code] if code < len(x.cpal) else x.missing_color for code in real_codes]
        if with_missing:
            plot_cpal.append(x.missing_color)

    return SeqImplic(
        indices=indices,
        levels=levels,
        alphabet=alph,
        positions=list(x.xtlab) if x.xtlab is not None else list(range(1, L + 1)),
        labels=plot_labels,
        type="indice",
        cpal=plot_cpal,
        xtstep=x.xtstep,
        tick_last=x.tick_last,
    )


@_require_stslist
def seqient(
    x: STSList,
    *,
    norm: bool = True,
    base: float = math.e,
    include_missing: bool = False,
    silent: bool = True,
) -> np.ndarray:
    """Within-sequence entropy (longitudinal entropy).

    Computes Shannon entropy for each sequence based on the distribution
    of states within that sequence. Measures the diversity/variability of
    states within individual trajectories.

    Parameters
    ----------
    x : STSList
        Input sequence object
    norm : bool, default True
        If True, normalize by log(K) where K is alphabet size. As in R, a
        one-state alphabet gives NaN (0/0).
    base : float, default e
        Logarithm base used for the entropy itself. TraMineR's normalization
        denominator remains the natural logarithm of the state count.
    include_missing : bool, default False
        If True, include missing states in entropy calculation. Structural void
        is never treated as a state.
    silent : bool, default True
        If False, emit TraMineR's entropy progress message on stderr.

    Returns
    -------
    np.ndarray
        Array of shape (n,) with entropy value for each sequence

    Examples
    --------
    >>> import numpy as np
    >>> from traminex import seqdef, seqient
    >>> data = np.array([
    ...     ['A', 'A', 'A', 'A'],  # Low entropy (all same)
    ...     ['A', 'B', 'C', 'D'],  # High entropy (all different)
    ...     ['A', 'A', 'B', 'B'],  # Medium entropy
    ... ])
    >>> sts = seqdef(data)
    >>> ent = seqient(sts)
    >>> bool(ent[0] < ent[2] < ent[1])
    True
    """
    x = x._canonicalized()
    try:
        base = float(base)
    except (TypeError, ValueError) as exc:
        raise ValueError("base must be finite, positive, and not equal to 1") from exc
    if not np.isfinite(base) or base <= 0 or base == 1:
        raise ValueError("base must be finite, positive, and not equal to 1")

    if not silent:
        print(
            f" [>] computing entropy for {x.nseq} sequences ...",
            file=sys.stderr,
        )

    # Get per-sequence state counts
    istatd = np.asarray(seqistatd(x, include_missing=include_missing), dtype=float)  # (n, K)
    # Distribution is over statl (real alphabet, + nr when include_missing); the
    # structural void is NEVER a state, so on ragged data its column must not
    # enter the per-sequence distribution (else the entropy is inflated).
    K = istatd.shape[1]
    statl = [j for j in range(K) if j != x.void_code and (include_missing or j != x.nr_code)]

    # Compute entropy for each sequence
    ient = np.zeros(x.nseq, dtype=float)
    log_base = np.log(base)

    for i in range(x.nseq):
        counts = istatd[i, statl]
        total = counts.sum()
        if total > 0:
            p = counts / total
            p = p[p > 0]
            ient[i] = -np.sum(p * np.log(p) / log_base)
        else:
            ient[i] = 0.0

    # Normalize by maximum entropy log(K_eff), where K_eff is the effective
    # with.missing alphabet size = k_real + (1 if include_missing). This always
    # counts the nr state under include_missing — even when nr is absent from the
    # data (R's alphabet(with.missing) appends a phantom nr) — and never counts void.
    if norm:
        miss = {c for c in (x.nr_code, x.void_code) if c is not None}
        k_real = len([i for i in range(len(x.alphabet)) if i not in miss])
        K_eff = k_real + (1 if include_missing else 0)
        # R divides unconditionally: with one state log(1) = 0 and the
        # normalised entropy is NaN (0/0), as in seqient.R:22-25 (audit S3-10).
        with np.errstate(divide="ignore", invalid="ignore"):
            ient = ient / np.log(K_eff)

    return ient


@_require_stslist
def seqici(
    x: STSList,
    *,
    include_missing: bool = False,
    silent: bool = True,
) -> np.ndarray:
    """Complexity index (longitudinal entropy × transition rate).

    Combines within-sequence entropy and normalized transition count to
    measure sequence complexity. Higher values indicate more complex
    trajectories with both high diversity and high turnover.

    Formula: C = sqrt(entropy * transitions)

    Parameters
    ----------
    x : STSList
        Input sequence object
    include_missing : bool, default False
        If True, include the missing state nr; structural void remains excluded
    silent : bool, default True
        Suppress TraMineR's computation progress message.

    Returns
    -------
    np.ndarray
        Array of shape (n,) with complexity index for each sequence

    Notes
    -----
    This is the TraMineR seqici implementation. The complexity index
    captures both the diversity of states (entropy) and the amount of
    change (transitions).

    Examples
    --------
    >>> import numpy as np
    >>> from traminex import seqdef, seqici
    >>> data = np.array([
    ...     ['A', 'A', 'A', 'A'],  # Low complexity (no diversity, no change)
    ...     ['A', 'B', 'C', 'D'],  # High complexity (high diversity, many changes)
    ...     ['A', 'A', 'A', 'B'],  # Low complexity (low diversity, one change)
    ... ])
    >>> sts = seqdef(data)
    >>> complexity = seqici(sts)
    >>> bool(complexity[0] < complexity[2] < complexity[1])
    True
    """
    if not isinstance(silent, (bool, np.bool_)):
        raise TypeError("silent must be a boolean")
    if not silent:
        print(
            f" [>] computing complexity index for {x.nseq} sequences ...",
            file=sys.stderr,
        )

    # Number of transitions (normalized) = R seqtransn(norm=TRUE) =
    # (DSS length - 1) / (non-void seqlength - 1). Computed over void-stripped rows
    # so structural void padding is never counted as a state/transition (R treats
    # void as NA); matches the per-sequence transition term used by seqivolatility.
    # Length-1 rows get 0 transitions from the loop below, and the entropy term
    # still decides the result: NaN for a one-state alphabet, as in R's
    # sqrt(0 * NaN) (audit S3-10).
    from .stats_complexity import _dss_and_durations, _strip_void

    trans_norm = np.zeros(x.nseq, dtype=float)
    for i in range(x.nseq):
        row = _strip_void(x.data[i], x.void_code)
        seql = row.shape[0]
        dss, _ = _dss_and_durations(
            row, include_missing=include_missing, nr_code=x.nr_code, void_code=None
        )
        ntrans = float(max(dss.size - 1, 0))
        trans_norm[i] = (ntrans / float(seql - 1)) if seql > 1 else 0.0

    # Longitudinal entropy (normalized)
    ient = seqient(x, norm=True, include_missing=include_missing)

    # Complexity index
    comp_index = np.sqrt(trans_norm * ient)

    return comp_index


def seqgen(
    n: int,
    length: int,
    alphabet: Sequence[Any] | str | int | float,
    *,
    p: Optional[Sequence[float]] = None,
    seed: Optional[int] = None,
) -> np.ndarray:
    """Generate random sequences.

    Creates n sequences of specified length by sampling states from the
    given alphabet with optional probabilities.

    Parameters
    ----------
    n : int
        Number of sequences to generate
    length : int
        Length of each sequence
    alphabet : sequence of states, or one state
        States to sample from: text, numbers or booleans. The draws keep the
        alphabet's own type, as R's ``sample()`` does. One value, even given as
        a scalar, is drawn as itself; R reads one number ``x >= 1`` as ``1:x``
        (``seqgen(.., alphabet = 5)`` draws 1 to 5), which traminex corrects
        (TX-AUD-20261011-001).
    p : sequence of float, optional
        Probability weights over the alphabet, one per state; normalized, as R's
        ``sample()`` does (``[5, 3, 2]`` draws as ``[0.5, 0.3, 0.2]``). A negative,
        missing or infinite weight, or no positive weight, is refused with
        ``sample()``'s message. If None, uniform distribution.
    seed : int, optional
        Random seed for reproducibility

    Returns
    -------
    np.ndarray
        Array of shape (n, length) of the alphabet's type (text, integer, float or
        boolean) with the randomly generated sequences.

    Important
    ---------
    The return value is a raw numpy array of states, **not** an STSList.
    To use the result with any traminex statistical or distance function you
    must convert it first::

        data = seqgen(10, 5, ['A', 'B', 'C'], seed=42)
        sts = seqdef(data)  # <- required before passing to seqdist, seqici, etc.

    Examples
    --------
    >>> from traminex import seqgen, seqdef
    >>> # Generate 10 sequences of length 5 from alphabet ['A', 'B', 'C']
    >>> data = seqgen(10, 5, ['A', 'B', 'C'], seed=42)
    >>> data.shape
    (10, 5)
    >>> # Convert to STSList before using with traminex functions
    >>> sts = seqdef(data, alphabet=['A', 'B', 'C'])
    >>> # With probabilities (70% A, 20% B, 10% C)
    >>> data2 = seqgen(10, 5, ['A', 'B', 'C'], p=[0.7, 0.2, 0.1], seed=42)
    >>> sts2 = seqdef(data2, alphabet=['A', 'B', 'C'])
    """
    rng = np.random.default_rng(seed)
    # The draws keep the alphabet's own type, as R's sample() does; one value is
    # drawn as itself, where R reads one number x >= 1 as 1:x (TX-AUD-20261011-001).
    alphabet_arr = np.atleast_1d(np.asarray(alphabet))
    if alphabet_arr.ndim != 1:
        raise ValueError("alphabet must be one state or a one-dimensional sequence of states")

    if p is not None:
        p_arr = np.array(p, dtype=float)
        if p_arr.ndim != 1 or len(p_arr) != len(alphabet_arr):
            raise ValueError("p must have same length as alphabet")
        # p holds probability weights, as R's sample() reads them: each checked in
        # turn, then normalized (TX-AUD-20261011-072).
        for weight in p_arr:
            if not np.isfinite(weight):
                raise ValueError("NA in probability vector")
            if weight < 0:
                raise ValueError("negative probability")
        if not np.any(p_arr > 0):
            raise ValueError("too few positive probabilities")
        p_arr = p_arr / p_arr.sum()
    else:
        p_arr = None

    # Generate random sequences
    # Numbers and booleans keep their dtype; text stays an object array of states.
    result = np.empty(
        (n, length), dtype=object if alphabet_arr.dtype.kind in "USO" else alphabet_arr.dtype
    )
    for i in range(n):
        result[i, :] = rng.choice(alphabet_arr, size=length, replace=True, p=p_arr)

    return result


@_require_stslist
def seqrecode(
    x: STSList,
    recodes: Dict[str, List[str]],
    *,
    otherwise: Optional[str] = None,
    labels: Optional[Mapping[str, str] | Sequence[str]] = None,
    cpal: Optional[Mapping[str, str] | Sequence[str]] = None,
) -> STSList:
    """Recode sequence alphabet with mapping rules.

    Groups states according to a recoding scheme, creating a new STSList
    with a reduced alphabet.

    Parameters
    ----------
    x : STSList
        Input sequence object
    recodes : dict[str, list[str]]
        Recoding rules: {new_state: [old_states_to_map]}
        Example: {'Work': ['FT', 'PT'], 'NotWork': ['UE', 'Inactive']}
    otherwise : str, optional
        Code for the states not mentioned in recodes. If None, they retain
        their original codes. The missing and void markers are not states:
        they stay markers unless recodes names them, so no sequence changes
        length (R also sends them to ``otherwise``; TX-AUD-20261019-033).
        An ``otherwise`` naming a marker itself (not a real state spelled
        like an unused one) sends every unlisted cell to that marker, as
        R's ``seqipos`` uses it.
    labels : dict[str, str] or sequence of str, optional
        Labels for the new states: a mapping, or, as in R, one label per new
        state in alphabet order. If None, use state codes.
    cpal : dict[str, str] or sequence of str, optional
        Colours for the new states, given the same way. If None, as R's
        seqrecode (seqrecode.R:88-104): a recoded state takes the colour of
        its first source state, the ``otherwise`` state that of the first
        original state no recode uses, a retained state keeps its own; a
        missing source gives ``missing_color`` (TX-AUD-20261011-024).

    Returns
    -------
    STSList
        New STSList with recoded alphabet

    Examples
    --------
    >>> import numpy as np
    >>> from traminex import seqdef, seqrecode
    >>> data = np.array([
    ...     ['FT', 'FT', 'PT', 'UE'],
    ...     ['PT', 'PT', 'FT', 'FT'],
    ... ])
    >>> sts = seqdef(data)
    >>> # Recode: FT/PT -> Work, UE -> NotWork
    >>> recoded = seqrecode(sts, {
    ...     'Work': ['FT', 'PT'],
    ...     'NotWork': ['UE']
    ... })
    >>> recoded.alphabet
    ['Work', 'NotWork']
    """
    # Both marker representations are equivalent inputs (the identity for
    # in-alphabet objects).  The mapping loops below iterate ``x.alphabet``, so
    # a declared *external* ``nr_code``/``void_code`` -- a code outside
    # ``range(len(alphabet))`` -- would never be visited, so its cells would be
    # handled unlike an in-alphabet marker's (cycle-10 N-15).
    x = x._canonicalized()

    # Build mapping from old to new codes
    old_to_new = {}
    new_alphabet = list(recodes.keys())

    # Map specified recodings
    for new_state, old_states in recodes.items():
        for old_state in old_states:
            if old_state in x.alphabet:
                old_to_new[old_state] = new_state

    # Handle otherwise: it recodes the real states recodes does not list. The
    # missing and void markers stay markers unless recodes names them, so no
    # sequence changes length (owner decision 2026-10-23, TX-AUD-20261019-033);
    # R's recodeFunc sends them to otherwise too (seqrecode.R:69-77). An
    # otherwise naming a marker itself (seqipos's otherwise = void, as R's
    # seqipos), not a real state spelled like one, still sends every unlisted
    # cell to that marker.
    marker_states = {x.alphabet[code] for code in (x.nr_code, x.void_code) if code is not None}
    real_states = {
        state for code, state in enumerate(x.alphabet) if code not in (x.nr_code, x.void_code)
    }
    if otherwise in (x.nr, x.void) and otherwise not in real_states:
        marker_states = set()
    if otherwise is not None:
        if otherwise not in new_alphabet:
            new_alphabet.append(otherwise)
        for state in x.alphabet:
            if state in old_to_new:
                continue
            if state in marker_states:
                old_to_new[state] = state
                # A recode named like the marker already lists it.
                if state not in new_alphabet:
                    new_alphabet.append(state)
            else:
                old_to_new[state] = otherwise
    else:
        # States not in recodes keep their original code
        for state in x.alphabet:
            if state not in old_to_new:
                old_to_new[state] = state
                if state not in new_alphabet:
                    new_alphabet.append(state)

    # Apply recoding to data
    new_data = np.empty_like(x.data, dtype=object)
    for i in range(x.nseq):
        for j in range(x.seqlen):
            old_state = x._state_for_code(int(x.data[i, j]))
            new_data[i, j] = old_to_new.get(old_state, old_state)

    # Build new STSList
    from .format import seqdef

    marker_codes = {code for code in (x.nr_code, x.void_code) if code is not None}
    marker_texts = {x.alphabet[code] for code in marker_codes}
    real_new_states = [state for state in new_alphabet if state not in marker_texts]

    def per_state(values: Mapping[str, str] | Sequence[str], name: str) -> dict[str, str]:
        if isinstance(values, Mapping):
            return dict(values)
        listed = list(values)
        if len(listed) != len(real_new_states):
            raise ValueError(
                f"{name} must give one entry per new state ({len(real_new_states)}: "
                f"{real_new_states}), or be a mapping"
            )
        return dict(zip(real_new_states, listed))

    # Prepare labels
    if labels is None:
        new_labels = None
    else:
        label_of = per_state(labels, "labels")
        new_labels = [label_of.get(s, s) for s in new_alphabet]

    def inherited_colour(source: str | None) -> str:
        if source is None or source not in x.alphabet:
            return "#808080"
        code = x.alphabet.index(source)
        if code == x.nr_code:
            return x.missing_color
        if code == x.void_code:
            warnings.warn(
                f"[!] no previous color defined for the 'void' ({source}) state, setting it to "
                "'white'",
                UserWarning,
                stacklevel=3,
            )
            return "white"
        return x.cpal[code] if x.cpal is not None and code < len(x.cpal) else "#808080"

    # Prepare colours as R's seqrecode does (seqrecode.R:88-104).
    if cpal is None and x.cpal is not None:
        used = {old for olds in recodes.values() for old in olds}
        others = [
            state
            for code, state in enumerate(x.alphabet)
            if code not in marker_codes and state not in used
        ]
        new_cpal_list = []
        for state in new_alphabet:
            if state in recodes:
                sources = list(recodes[state])
                new_cpal_list.append(inherited_colour(sources[0] if sources else None))
            elif otherwise is not None and state == otherwise:
                new_cpal_list.append(inherited_colour(others[0] if others else None))
            elif state in marker_texts:
                new_cpal_list.append(x.cpal[x.alphabet.index(state)])
            else:
                new_cpal_list.append(inherited_colour(state))
    elif cpal is not None:
        colour_of = per_state(cpal, "cpal")
        new_cpal_list = [
            (
                colour_of.get(s, "#808080")
                if s not in marker_texts or x.cpal is None
                else colour_of.get(s, x.cpal[x.alphabet.index(s)])
            )
            for s in new_alphabet
        ]
    else:
        new_cpal_list = None

    # Decoding loses the distinction between a real marker-shaped state and
    # a declared marker. Keep absent markers out of the recoded real alphabet
    # while the STS factory materializes it; restore their unused text below.
    factory_nr = (
        x.nr
        if x.nr_code is not None or x.nr not in x.alphabet
        else _collision_free_marker(x.nr, [*new_alphabet, x.void], "NR")
    )
    factory_void = (
        x.void
        if x.void_code is not None or x.void not in x.alphabet
        else _collision_free_marker(x.void, [*new_alphabet, factory_nr], "VOID")
    )
    result = seqdef(
        new_data,
        alphabet=new_alphabet,
        labels=new_labels,
        weights=cast(Optional[Sequence[float]], x.weights),
        cnames=None if x.xtlab is None else list(x.xtlab),
        xtstep=x.xtstep,
        tick_last=x.tick_last,
        start=x.start,
        id=None if x.ids is None else list(x.ids),
        nr=factory_nr,
        void=factory_void,
        missing_color=x.missing_color,
    )
    if result.nr_code is None:
        result.nr = x.nr
    if result.void_code is None:
        result.void = x.void

    # Colours go to the states by name: seqdef places the markers after the
    # real states, so the list built in recode order does not line up with the
    # result's alphabet when a recode is named like nr or void (batch-3 review I1).
    if new_cpal_list is not None:
        colour_by_state = dict(zip(new_alphabet, new_cpal_list))
        result.cpal = [
            (
                result.missing_color
                if code == result.nr_code
                else "white" if code == result.void_code else colour_by_state.get(state, "#808080")
            )
            for code, state in enumerate(result.alphabet)
        ]
    result.issues = (
        None
        if x.issues is None
        else {category: list(indexes) for category, indexes in x.issues.items()}
    )

    return result


# ---------------------------------------------------------------------------
# Sequence-operation indicators (RP6a; BUILD-from-R, TraMineR 2.2-14 GPL source)
# ---------------------------------------------------------------------------


def _sequence_codes(seq, *, exactly_one: bool = False):
    """Return the coded matrix and marker metadata of an STS sequence object."""
    if not isinstance(seq, STSList):
        raise TypeError("sequences must be sequence objects")
    if exactly_one and seq.nseq != 1:
        raise ValueError("sequence objects must contain exactly one sequence")
    return np.asarray(seq.data), seq.nr_code, seq.void_code


def _sequence_labels(seq: STSList, codes: np.ndarray) -> np.ndarray:
    """Decode ``seq`` codes to their state or marker texts (R's factor labels)."""
    labels = np.empty(codes.shape, dtype=object)
    for code in np.unique(codes):
        labels[codes == code] = str(seq._state_for_code(int(code)))
    return labels


class _SeqConcVoidDefault:
    def __repr__(self) -> str:
        return "NA"


_SEQCONC_VOID_DEFAULT = _SeqConcVoidDefault()


class _SeqConcResult(np.ndarray):
    """One-dimensional array with the observable R column metadata."""

    vname: str
    row_names: tuple[str, ...]

    def __new__(
        cls,
        values: Sequence[str],
        *,
        vname: str,
        row_names: Sequence[object],
    ) -> "_SeqConcResult":
        result = np.asarray(values, dtype=object).view(cls)
        result.vname = vname
        result.row_names = tuple(str(value) for value in row_names)
        return result

    def __array_finalize__(self, source) -> None:
        if source is None:
            return
        self.vname = getattr(source, "vname", "Sequence")
        self.row_names = getattr(source, "row_names", ())


def _seqconc_vname(vname: object) -> str:
    if isinstance(vname, (str, bytes)) or np.isscalar(vname):
        return str(vname)
    raise ValueError("vname must have length one")


def _seqconc_var(var: object) -> object:
    if var is None:
        return None
    if isinstance(var, np.ndarray):
        if var.ndim != 1:
            raise ValueError("'var' must be a vector of one-based integers or strings")
        return None if var.size == 0 else var.tolist()
    if np.isscalar(var) and _cell_is_na(var):
        return None
    if isinstance(var, Sequence) and not isinstance(var, (str, bytes, bytearray)) and len(var) == 0:
        return None
    return var


def _seqconc_raw_matrix(
    data,
    var: Sequence[int | str] | int | str | None,
) -> tuple[np.ndarray, list[object]]:
    normalized_var = _seqconc_var(var)
    source = data
    if not isinstance(data, (pd.DataFrame, pd.Series, np.ndarray)):
        source = np.asarray(data, dtype=object)
        if source.ndim == 0:
            source = source.reshape(1)

    if (
        isinstance(normalized_var, (int, np.integer))
        and not isinstance(normalized_var, (bool, np.bool_))
        and int(normalized_var) == 0
    ):
        if isinstance(source, pd.DataFrame):
            source = source.iloc[:, :0]
        elif isinstance(source, np.ndarray):
            matrix = source.reshape(-1, 1) if source.ndim == 1 else source
            if matrix.ndim != 2:
                raise ValueError("'data' must be a data frame, matrix, or vector")
            source = matrix[:, :0]
        else:
            raise ValueError("'var' must be a vector of one-based integers or strings")
        normalized_var = None

    matrix, ids, _ = _table_to_matrix(
        source,
        var=cast(Sequence[int | str] | int | str | None, normalized_var),
        allow_empty_columns=True,
    )
    if ids is not None:
        row_names: list[object] = ids
    elif isinstance(data, pd.DataFrame):
        row_names = list(range(1, matrix.shape[0] + 1))
    else:
        row_names = [f"[{index}]" for index in range(1, matrix.shape[0] + 1)]
    return matrix, row_names


def _seqconc_cell_text(value: object) -> str:
    return "NA" if _cell_is_na(value) else _r_state_text(value)


def seqconc(
    data,
    var: Sequence[int | str] | int | str | None = None,
    sep: str = "-",
    vname: str = "Sequence",
    void=_SEQCONC_VOID_DEFAULT,
) -> np.ndarray:
    """Concatenate sequence rows with the complete public R input surface.

    ``STSList`` inputs use their own structural-void marker, as TraMineR does;
    raw vectors, matrices, Series, and data frames support one-based/name
    selectors, literal separators, and missing or explicit void removal.

    The return remains a one-dimensional ndarray for established Python use.
    Its private ndarray subclass exposes R's observable column metadata as
    ``vname`` and ``row_names``. R's explicit ``void=NULL`` error is retained
    for raw data.
    """
    if not isinstance(sep, str):
        raise TypeError("sep must be a string")
    column_name = _seqconc_vname(vname)

    if isinstance(data, STSList):
        values = []
        for row in data.data:
            codes = row[row != data.void_code] if data.void_code is not None else row
            values.append(sep.join(data._state_for_code(int(code)) for code in codes))
        row_names = list(data.ids) if data.ids is not None else list(range(1, data.nseq + 1))
        return _SeqConcResult(values, vname=column_name, row_names=row_names)

    if void is None:
        raise ValueError("void=NULL has no value; argument is of length zero")
    matrix, row_names = _seqconc_raw_matrix(data, var)
    drop_missing = void is _SEQCONC_VOID_DEFAULT or _cell_is_na(void)
    if not drop_missing and not np.isscalar(void):
        raise ValueError("void must have length one")

    values = []
    for row in matrix:
        kept = []
        for value in row:
            if drop_missing:
                if _cell_is_na(value):
                    continue
            elif not _cell_is_na(value) and value == void:
                continue
            kept.append(_seqconc_cell_text(value))
        values.append(sep.join(kept))
    return _SeqConcResult(values, vname=column_name, row_names=row_names)


def seqsep(seqdata, sl: int = 1, sep: str = "-"):
    """Separate fixed-width states in character sequences (R ``seqsep``).

    This is an exact port of TraMineR 2.2-14 ``seqsep`` for character vectors,
    with container-preserving Python returns: a list, one-dimensional ndarray,
    or pandas Series produces the same container category. Series index/name and
    string dtype are retained; ndarray object/unicode dtype categories are retained
    without truncating expanded strings.

    TraMineR's source retains the previous separated value when it encounters an
    empty string after a non-empty element. If the first element is empty it aborts
    because that value is not yet defined; both observed behaviors are preserved.
    """
    if isinstance(sl, bool) or not isinstance(sl, (int, np.integer)) or sl <= 0:
        raise ValueError("sl must be a positive integer")
    if not isinstance(sep, str):
        raise TypeError("sep must be a string")
    state_length = int(sl)

    container = None
    result_dtype = None
    if isinstance(seqdata, list):
        arr = np.asarray(seqdata, dtype=object)
        if arr.ndim != 1:
            raise ValueError("seqdata must be one-dimensional")
        values = list(seqdata)
        container = "list"
    elif isinstance(seqdata, np.ndarray):
        if seqdata.ndim != 1:
            raise ValueError("seqdata must be a one-dimensional numpy.ndarray")
        values = seqdata.tolist()
        container = "ndarray"
        result_dtype = object if seqdata.dtype.kind == "O" else str
    else:
        import pandas as pd

        if not isinstance(seqdata, pd.Series):
            raise TypeError(
                "seqdata must be a list, one-dimensional numpy.ndarray, or pandas.Series"
            )
        values = seqdata.tolist()
        container = "series"

    if not values:
        raise ValueError("missing value where TRUE/FALSE needed")

    result = []
    separated = None
    error = "Number of characters does not match number of states*states length in sequence"
    for position, value in enumerate(values, start=1):
        if not isinstance(value, str):
            raise TypeError(f"seqdata element {position} must be a string")
        length = len(value)
        if length % state_length != 0:
            raise ValueError(f"{error}{position}")
        if length > 0:
            separated = sep.join(
                value[start : start + state_length] for start in range(0, length, state_length)
            )
        if separated is None:
            raise ValueError("object 'subseq' not found")
        result.append(separated)

    if container == "list":
        return result
    if container == "ndarray":
        return np.asarray(result, dtype=result_dtype)

    result_dtype = seqdata.dtype if isinstance(seqdata.dtype, pd.StringDtype) else object
    return pd.Series(result, index=seqdata.index, name=seqdata.name, dtype=result_dtype)


class _SeqDecompResult(np.ndarray):
    """Object matrix carrying the observable R dimnames."""

    row_names: tuple[str, ...]
    column_names: tuple[str, ...]

    def __new__(
        cls,
        values: np.ndarray,
        *,
        row_names: Sequence[object],
        column_names: Sequence[object],
    ) -> "_SeqDecompResult":
        result = np.asarray(values, dtype=object).view(cls)
        result.row_names = tuple(str(value) for value in row_names)
        result.column_names = tuple(str(value) for value in column_names)
        return result

    def __array_finalize__(self, source) -> None:
        if source is None:
            return
        self.row_names = getattr(source, "row_names", ())
        self.column_names = getattr(source, "column_names", ())


def _seqdecomp_named_var(data, var):
    if not isinstance(data, _SeqConcResult) or var is None:
        return var
    values = [var] if isinstance(var, str) else list(var)
    if not values or not all(isinstance(value, str) for value in values):
        return var
    if any(value != data.vname for value in values):
        raise ValueError("undefined columns selected")
    return [1] * len(values)


def _seqdecomp_explicit_row_names(data) -> list[object] | None:
    if isinstance(data, _SeqConcResult):
        return list(data.row_names)
    if isinstance(data, pd.DataFrame):
        return list(data.index)
    if isinstance(data, pd.Series) and not isinstance(data.index, pd.RangeIndex):
        return list(data.index)
    return None


def _seqdecomp_split(value: object, pattern: re.Pattern[str] | None) -> list[object]:
    if _cell_is_na(value):
        return [None]
    if not isinstance(value, str):
        raise TypeError("non-character argument")
    if value == "":
        return []
    if pattern is None:
        return list(value)

    parts: list[object] = []
    start = 0
    matched = False
    for match in pattern.finditer(value):
        matched = True
        parts.append(value[start : match.start()])
        start = match.end()
    if not matched:
        return [value]
    if start < len(value):
        parts.append(value[start:])
    return parts


def _seqdecomp_missing_values(miss: str | Sequence[str]) -> set[str]:
    if isinstance(miss, str):
        return {miss}
    try:
        values = list(miss)
    except TypeError as exc:
        raise ValueError("miss must be a string or sequence of strings") from exc
    if not all(isinstance(value, str) for value in values):
        raise ValueError("miss must be a string or sequence of strings")
    return set(values)


def _seqdecomp_column_names(vnames: object, width: int) -> list[object]:
    if vnames is None:
        return [f"[{index}]" for index in range(1, width + 1)]
    if isinstance(vnames, (str, bytes)) or np.isscalar(vnames):
        values: list[object] = [vnames]
    else:
        try:
            values = list(cast(Iterable[object], vnames))
        except TypeError as exc:
            raise ValueError("vnames must have one name per output column") from exc
    if len(values) != width:
        raise ValueError("vnames length does not match the output column count")
    return values


def seqdecomp(
    data,
    var: Sequence[int | str] | int | str | None = None,
    sep: str = "-",
    miss: str | Sequence[str] = "NA",
    vnames: Sequence[str] | str | None = None,
) -> np.ndarray:
    """Decompress state strings with TraMineR's public tokenization contract.

    Extraction applies ``var`` before R column-major flattening. ``sep`` uses
    regular-expression semantics (with ``""`` splitting single characters),
    each token listed in ``miss`` becomes ``None``, and rows are right-padded.
    The returned ndarray subclass exposes R row and column names through
    ``row_names`` and ``column_names`` while remaining an ordinary object
    array for established Python consumers.
    """
    if not isinstance(sep, str):
        raise TypeError("sep must be a string")
    try:
        pattern = None if sep == "" else re.compile(sep)
    except re.error as exc:
        raise ValueError(f"invalid separator regular expression: {exc}") from exc

    resolved_var = _seqdecomp_named_var(data, var)
    explicit_row_names = _seqdecomp_explicit_row_names(data)
    matrix, _ = _seqconc_raw_matrix(data, resolved_var)
    flat = matrix.ravel(order="F")
    if explicit_row_names is not None and len(explicit_row_names) != len(flat):
        raise ValueError("length of row dimnames does not equal the output extent")
    row_names = (
        explicit_row_names
        if explicit_row_names is not None
        else [f"[{index}]" for index in range(1, len(flat) + 1)]
    )

    missing_values = _seqdecomp_missing_values(miss)
    parts = [_seqdecomp_split(value, pattern) for value in flat]
    width = max((len(tokens) for tokens in parts), default=0)
    if width == 0:
        raise ValueError("empty decompression has zero columns and invalid dimnames")

    output = np.empty((len(parts), width), dtype=object)
    output[:] = None
    for row, tokens in enumerate(parts):
        for column, token in enumerate(tokens):
            output[row, column] = None if token is None or token in missing_values else token
    column_names = _seqdecomp_column_names(vnames, width)
    return _SeqDecompResult(
        output,
        row_names=row_names,
        column_names=column_names,
    )


def _is_missing_value(value: object) -> bool:
    """R's NA among Python values: None, a float NaN or pandas' NA."""
    if value is None or value is pd.NA:
        return True
    return isinstance(value, (float, np.floating)) and bool(np.isnan(value))


def seqdim(data) -> np.ndarray:
    """Dimensions ``[Rows, Columns]`` of a sequence object (R ``seqdim``).

    Faithful port of TraMineR 2.2-14 ``seqdim`` (GPL). For an :class:`STSList`
    returns ``[nseq, seqlen]``. Otherwise R's ``seqfcheck`` rule decides, over
    every element (seqdim.R:5-15, seqfcheck.R:3-23): numeric data is never
    decomposed; one-column character data (a vector, or a one-column matrix or
    data frame; numbers mixed with text count as text, as R's ``c()`` makes
    them) is compressed when any element holds ``"-"`` (else ``":"``) and
    is decomposed with :func:`seqdecomp`; data with several columns is never
    decomposed. A vector then counts as one sequence, ``[1, length]``; a matrix
    or data frame gives ``[nrow, ncol]`` (TX-AUD-20261011-071).
    """
    if hasattr(data, "data") and hasattr(data, "alphabet"):
        return np.array([data.nseq, data.seqlen], dtype=int)
    if isinstance(data, pd.DataFrame):
        matrix, vector = data.to_numpy(dtype=object), False
    elif isinstance(data, str) or np.ndim(data) == 0:
        matrix, vector = np.array([[data]], dtype=object), True
    else:
        arr = np.asarray(data, dtype=object)
        vector = arr.ndim == 1
        matrix = arr.reshape(-1, 1) if vector else arr
    present = [value for value in matrix.ravel() if not _is_missing_value(value)]
    numeric = bool(present) and all(
        isinstance(value, (int, float, np.integer, np.floating)) and not isinstance(value, bool)
        for value in present
    )
    if not numeric and matrix.shape[1] == 1:
        texts = [str(value) for value in present]
        # Data that is not all numeric is character data, as R's c() makes it.
        column = [value if _is_missing_value(value) else str(value) for value in matrix[:, 0]]
        for sep in ("-", ":"):
            if any(sep in text for text in texts):
                decomposed = np.asarray(seqdecomp(column, sep=sep))
                return np.array([decomposed.shape[0], decomposed.shape[1]], dtype=int)
    if vector:
        return np.array([1, matrix.shape[0]], dtype=int)
    return np.array([matrix.shape[0], matrix.shape[1]], dtype=int)


def seqmpos(seq1, seq2, *, with_missing: bool = False) -> int:
    """Number of matching positions between two sequences (R ``seqmpos``).

    Faithful port of TraMineR 2.2-14 ``seqmpos`` (GPL): counts positions where
    ``seq1`` and ``seq2`` carry the same state, excluding void positions always
    and missing (``nr``) positions unless ``with_missing=True``. States are
    compared by label, as R compares factor labels, so the two objects may use
    different alphabets or alphabet orders; each object's void and missing
    cells are identified by its own markers.
    """
    a, nr1, vc1 = _sequence_codes(seq1)
    b, nr2, vc2 = _sequence_codes(seq2)
    if a.shape != b.shape:
        raise ValueError("sequence objects must have matching dimensions")
    comp = np.ones(a.shape, dtype=bool)
    excluded = [(a, vc1), (b, vc2)]
    if not with_missing:
        excluded += [(a, nr1), (b, nr2)]
    for codes, marker in excluded:
        if marker is not None:
            comp &= codes != marker
    same = _sequence_labels(seq1, a) == _sequence_labels(seq2, b)
    return int(np.sum(same & comp))


def seqcomp(x, y) -> bool:
    """Whether two sequences are identical (R ``seqcomp``).

    Faithful port of TraMineR 2.2-14 ``seqcomp`` (GPL): ``True`` iff the two
    sequences have the same non-void length and equal states at every position.
    States are compared by label, as R compares factor labels, so the two
    objects may use different alphabets or alphabet orders.
    """
    a, _, vc1 = _sequence_codes(x, exactly_one=True)
    b, _, vc2 = _sequence_codes(y, exactly_one=True)
    a = a[0]
    b = b[0]
    if vc1 is not None:
        a = a[a != vc1]
    if vc2 is not None:
        b = b[b != vc2]
    if len(a) != len(b):
        return False
    return bool(np.all(_sequence_labels(x, a) == _sequence_labels(y, b)))


def _nsubs(codes: list) -> int:
    """Number of distinct subsequences (incl. empty) of a coded sequence.

    Port of TraMineR's ``nsubs`` dynamic program: ``N[i] = 2*N[i-1]`` minus the
    contribution double-counted at a state's previous occurrence.
    """
    sl = len(codes)
    if sl == 0:
        return 1
    N = [0] * (sl + 2)
    N[1] = 1
    last: Dict[int, int] = {}
    for i in range(2, sl + 2):
        N[i] = 2 * N[i - 1]
        c = int(codes[i - 2])
        lv = last.get(c, 0)
        if lv > 0:
            N[i] -= N[lv]
        last[c] = i - 1
    return N[sl + 1]


@_require_stslist
def seqsubsn(x: STSList, *, dss: bool = True, with_missing: bool = False) -> np.ndarray:
    """Number of distinct subsequences per sequence (R ``seqsubsn``).

    Faithful port of TraMineR 2.2-14 ``seqsubsn`` (GPL). With ``dss=True``
    (default) sequences are first collapsed to their distinct successive states
    (DSS); then the count of distinct subsequences (including the empty one) is
    computed over the void-stripped codes (missing ``nr`` is also stripped unless
    ``with_missing=True``). Returns a length-n float array: R returns an
    ``n x 1`` matrix (column "Subseq."), and traminex's per-sequence indicators
    are 1-D (docs/API_CONSISTENCY.md section 1c; TX-AUD-20261010-029).

    The recurrence uses exact integers before final float64 rounding (including
    overflow to Inf), deliberately avoiding R's intermediate Inf - Inf for very
    long repeated-state sequences.
    """
    x = x._canonicalized()
    from .stats_complexity import _subsequence_count_as_float, seqdss

    src = seqdss(x, include_missing=with_missing) if dss else x
    out = np.zeros(src.nseq, dtype=float)
    for i in range(src.nseq):
        row = src.data[i]
        if src.void_code is not None:
            row = row[row != src.void_code]
        if not with_missing and src.nr_code is not None:
            row = row[row != src.nr_code]
        out[i] = _subsequence_count_as_float(_nsubs(list(row)))
    return out


_POSIX_CLASSES: dict[str, Callable[[str], bool]] = {
    "alpha": str.isalpha,
    "digit": lambda c: c in "0123456789",
    "alnum": str.isalnum,
    "upper": str.isupper,
    "lower": str.islower,
    "space": str.isspace,
    "blank": lambda c: c in " \t",
    "punct": lambda c: c.isprintable() and not c.isalnum() and not c.isspace(),
    "cntrl": lambda c: unicodedata.category(c) == "Cc",
    "print": str.isprintable,
    "graph": lambda c: c.isprintable() and not c.isspace(),
    "xdigit": lambda c: c in "0123456789abcdefABCDEF",
}


def _translate_posix_classes(pattern: str, text: str) -> str:
    """Replace POSIX bracket classes (``[:alpha:]`` ...) for Python's ``re``.

    Each class inside a bracket expression becomes the escaped characters of
    ``text`` that satisfy it, so it matches in ``text`` exactly what R's
    locale-aware class matches there (TX-AUD-20261010-054). A class no
    character of ``text`` satisfies becomes one character absent from ``text``,
    which keeps the bracket valid and matches nothing.
    """
    out: list[str] = []
    i, n, in_bracket = 0, len(pattern), False
    while i < n:
        char = pattern[i]
        if not in_bracket:
            if char == "\\" and i + 1 < n:
                out.append(pattern[i : i + 2])
                i += 2
                continue
            out.append(char)
            i += 1
            if char == "[":
                in_bracket = True
                if i < n and pattern[i] == "^":
                    out.append("^")
                    i += 1
                if i < n and pattern[i] == "]":
                    out.append("\\]")
                    i += 1
            continue
        if pattern.startswith("[:", i):
            end = pattern.find(":]", i + 2)
            name = pattern[i + 2 : end] if end != -1 else ""
            if name in _POSIX_CLASSES:
                members = sorted({c for c in text if _POSIX_CLASSES[name](c)})
                if not members:
                    filler = 0xE000
                    while chr(filler) in text:
                        filler += 1
                    members = [chr(filler)]
                out.append("".join(re.escape(c) for c in members))
                i = end + 2
                continue
        if char == "[":
            out.append("\\[")
        elif char == "\\":
            # A backslash is an ordinary character inside an R bracket.
            out.append("\\\\")
        else:
            out.append(char)
            if char == "]":
                in_bracket = False
        i += 1
    return "".join(out)


def seqpm(x: STSList, pattern: str, *, sep: str = ""):
    """Find sequences matching a pattern (R ``seqpm``).

    Faithful port of TraMineR 2.2-14 ``seqpm`` (GPL): searches the regular
    expression ``pattern`` in the ``seqconc(x, sep=sep)`` strings. Returns a dict
    with ``MTab`` (a one-row DataFrame ``pattern``/``nbocc``) and ``MIndex`` (the
    1-based indices of matching sequences). POSIX bracket classes such as
    ``[[:upper:]]`` match as in R (Unicode classes, as R's under C.UTF-8), and a
    backslash inside a bracket is an ordinary character, as in R; the rest of
    the pattern is read by Python's ``re``.
    """
    import pandas as pd

    strs = [str(s) for s in seqconc(x, sep=sep)]
    compiled = re.compile(_translate_posix_classes(pattern, "".join(strs)))
    idx = [i + 1 for i, s in enumerate(strs) if compiled.search(s) is not None]
    return {
        "MTab": pd.DataFrame({"pattern": [pattern], "nbocc": [len(idx)]}),
        "MIndex": np.array(idx, dtype=int),
    }


@_require_stslist
def seqsha(
    x: STSList,
    time: Sequence[int],
    event: Sequence[bool],
    *,
    include_present: bool = False,
    align_end: bool = False,
    covar=None,
):
    """Sequence-history person-period converter (TraMineRextras ``seqsha``).

    Expands each sequence into one row per time step ``1..time[i]`` (the person-period
    form), flags ``event`` only at the last row of each id, and lays the *past* states
    into a state-history block. With ``align_end`` the history is right-aligned into
    ``Tm<k>`` columns; otherwise it is left-aligned into the original position columns.
    ``include_present`` toggles whether the current position is included (R's exact,
    slightly counter-intuitive, ``age > aa`` vs ``age >= aa`` convention is reproduced).
    Returns a DataFrame ``id, time, event, <history...>`` (optionally with ``covar``
    columns appended by id). Provenance: TraMineRextras 0.6.9 ``seqsha.R``.

    Note: in the start-aligned path with ``max(time) < ncol`` R aborts (it assigns
    ``ncol(seqdata)[1:ma]`` names to a wider matrix -- a length-mismatch error);
    traminex emits the full-width result instead (the de-bugged R intent). Documented
    divergence (audit-and-fix).

    As R, it stops on a missing time or event and on a time below 1 or not a
    whole number; ``event``
    is logical (booleans, or the numbers 0 and 1), as the Rd says, and text
    events are refused where R keeps them as a text column.
    """
    import pandas as pd

    n = x.nseq
    time_values = list(np.asarray(time, dtype=object).ravel())
    event_values = list(np.asarray(event, dtype=object).ravel())
    if len(time_values) != n or len(event_values) != n:
        raise ValueError("time and event must have one entry per sequence")
    # In R's order (TX-AUD-20261011-069): R stops on a missing value
    # (seqsha.R:6-8), then on a time below 1 (its rep() and [<- fail); the Rd
    # takes a logical event, where R would keep text events as a text column.
    if any(pd.isna(value) for value in time_values + event_values):
        raise ValueError(" [!] cannot handle missing data in the time or event variables")
    t_float = np.asarray(time_values, dtype=float)
    if np.any(t_float < 1):
        raise ValueError(
            "time must be at least 1 for every sequence (R stops on a time of 0 or below)"
        )
    # R also stops on a fractional or an infinite time (batch-15 review M-1).
    if not np.all(np.isfinite(t_float)) or np.any(t_float != np.floor(t_float)):
        raise ValueError(
            "time must be a whole number of periods for every sequence (R stops otherwise)"
        )
    if not all(
        isinstance(value, (bool, np.bool_))
        or (isinstance(value, (int, float, np.integer, np.floating)) and value in (0, 1))
        for value in event_values
    ):
        raise ValueError("event must be logical: booleans, or the numbers 0 and 1")
    t = t_float.astype(int)
    ev_in = np.asarray([bool(value) for value in event_values], dtype=bool)
    ids = np.repeat(np.arange(1, n + 1), t)
    periods = np.concatenate([np.arange(1, ti + 1) for ti in t]) if n else np.array([], dtype=int)
    nrows = ids.size
    ev = np.zeros(nrows, dtype=bool)
    if nrows:
        idmax = np.cumsum(t) - 1
        ev[idmax] = ev_in

    L = x.seqlen
    sdata = np.empty(x.data.shape, dtype=object)
    for index, code in np.ndenumerate(x.data):
        sdata[index] = x._state_for_code(int(code))
    age = periods
    ma = int(age.max()) if nrows else 0
    if ma > L:
        raise ValueError(
            " [!] Maximum time of event occurrence is higher than the longest sequence!"
        )

    past = np.full((nrows, L), None, dtype=object)
    if align_end:
        start = 1 if include_present else 2
        for aa in range(start, ma + 1):
            cond = age == aa
            ids_a = ids[cond] - 1
            if include_present:
                past[np.ix_(cond, np.arange(L - aa, L))] = sdata[ids_a, 0:aa]
            else:
                past[np.ix_(cond, np.arange(L - aa + 1, L))] = sdata[ids_a, 0 : (aa - 1)]
        colnames = [f"Tm{k}" for k in range(L, 0, -1)]
    else:
        for aa in range(1, ma + 1):
            cond = (age > aa) if include_present else (age >= aa)
            ids_a = ids[cond] - 1
            past[cond, aa - 1] = sdata[ids_a, aa - 1]
        colnames = list(x.xtlab) if x.xtlab else [f"T{k}" for k in range(1, L + 1)]

    df = pd.DataFrame({"id": ids, "time": periods, "event": ev})
    for j, cn in enumerate(colnames):
        df[cn] = past[:, j]
    if covar is not None:
        cov = pd.DataFrame(covar).reset_index(drop=True)
        df = pd.concat(
            [df.reset_index(drop=True), cov.iloc[ids - 1].reset_index(drop=True)], axis=1
        )
    return df


_GRAN_MISS = -1


def _granularity_block_valid(block: np.ndarray, real_codes: set, first: bool) -> np.ndarray:
    """First / last valid (real-state) code per row of a block (R ``firstvalid``)."""
    out = np.full(block.shape[0], _GRAN_MISS, dtype=int)
    for r in range(block.shape[0]):
        pos = [j for j in range(block.shape[1]) if int(block[r, j]) in real_codes]
        if pos:
            out[r] = int(block[r, pos[0] if first else pos[-1]])
    return out


def _granularity_block_mostfreq(block: np.ndarray, real_codes: list) -> np.ndarray:
    """Most-frequent real-state code per row (R ``seqistatd`` + ``which.max``)."""
    out = np.full(block.shape[0], _GRAN_MISS, dtype=int)
    rc = np.asarray(real_codes, dtype=int)
    for r in range(block.shape[0]):
        counts = np.array([(block[r] == c).sum() for c in rc])
        if counts.max() > 0:
            out[r] = int(rc[int(np.argmax(counts))])  # first max, alphabet order
    return out


@_require_stslist
def seqgranularity(x: STSList, *, tspan: int = 3, method: str = "last") -> STSList:
    """Coarsen the time granularity of a sequence object (TraMineRextras ``seqgranularity``).

    Aggregates every ``tspan`` consecutive positions into one, picking the
    ``first`` / ``last`` cell of each block, the ``first.valid`` / ``last.valid``
    (in-alphabet) cell, or the ``mostfreq`` state (a tie goes to the first of the tied
    states in the alphabet order, as R's ``which.max`` over the alphabet-ordered
    frequencies does, seqgranularity.R:56; without an explicit ``alphabet=`` that
    order is R's ``seqdef`` order). Returns a new ``STSList`` of length
    ``ceil(maxlen/tspan)``. Trailing missing collapse to void (R's internal re-``seqdef``
    with ``right="DEL"`` when the input has void), otherwise to ``nr``. Provenance:
    TraMineRextras 0.6.9 ``seqgranularity.R``.
    """
    import math

    if tspan == 1:
        return x.copy()
    if tspan < 1:
        raise ValueError("tspan must be greater than 1")
    methods = ("first", "first.valid", "last", "last.valid", "mostfreq")
    if method not in methods:
        raise ValueError(f"method must be one of: {' '.join(methods)}")

    n = x.nseq
    nr, void = x.nr_code, x.void_code
    if void is not None:
        lgth = int((x.data != void).sum(axis=1).max())
    else:
        lgth = x.seqlen
    real_codes = [i for i in range(len(x.alphabet)) if i not in (void, nr)]
    real_set = set(real_codes)

    new_lgth = math.ceil(lgth / tspan)
    new_lgth_f = lgth // tspan
    newseq = np.full((n, new_lgth), _GRAN_MISS, dtype=int)

    prev = (tspan - 1) if method == "first" else 0
    for i in range(1, new_lgth_f + 1):
        cols = slice(tspan * (i - 1), tspan * i)
        if method in ("first", "last"):
            newseq[:, i - 1] = x.data[:, tspan * i - prev - 1]
        elif method in ("first.valid", "last.valid"):
            newseq[:, i - 1] = _granularity_block_valid(
                x.data[:, cols], real_set, method == "first.valid"
            )
        else:  # mostfreq
            newseq[:, i - 1] = _granularity_block_mostfreq(x.data[:, cols], real_codes)

    if new_lgth > new_lgth_f:
        last = new_lgth - 1
        if lgth > tspan * new_lgth_f + 1:
            block = x.data[:, tspan * new_lgth_f : lgth]
            if method == "first":
                newseq[:, last] = x.data[:, tspan * new_lgth_f]
            elif method in ("first.valid", "last.valid"):
                newseq[:, last] = _granularity_block_valid(block, real_set, method == "first.valid")
            elif method == "mostfreq":
                newseq[:, last] = _granularity_block_mostfreq(block, real_codes)
            else:  # last
                newseq[:, last] = x.data[:, lgth - 1]
        else:
            newseq[:, last] = x.data[:, lgth - 1]

    # Original nr / void cells become missing (R: newseq[==nr|==void] <- missing.char)
    if nr is not None:
        newseq[newseq == nr] = _GRAN_MISS
    if void is not None:
        newseq[newseq == void] = _GRAN_MISS

    is_void = void is not None and bool((x.data == void).any())
    new_alph = list(x.alphabet)
    nr_code, void_code = nr, void
    new_nr = x.nr
    if (newseq == _GRAN_MISS).any():
        if nr_code is None:
            nr_code = len(new_alph)
            new_nr = _collision_free_marker(x.nr, [*new_alph, x.void], "NR")
            new_alph.append(new_nr)
        if is_void:
            if void_code is None:
                void_code = len(new_alph)
                new_alph.append(x.void)
            for r in range(n):
                c = new_lgth
                while c > 0 and newseq[r, c - 1] == _GRAN_MISS:
                    c -= 1
                newseq[r, c:] = void_code
                seg = newseq[r, :c]
                seg[seg == _GRAN_MISS] = nr_code
            newseq[newseq == _GRAN_MISS] = nr_code
        else:
            newseq[newseq == _GRAN_MISS] = nr_code

    label_by_state = dict(zip(x.alphabet, x.labels)) if x.labels is not None else None
    color_by_state = dict(zip(x.alphabet, x.cpal)) if x.cpal is not None else None
    new_xtlab: list[str] | None = None
    if x.xtlab is not None:
        xtlab_positions = list(range(0, lgth, tspan))
        if new_lgth > new_lgth_f and xtlab_positions:
            xtlab_positions[-1] = lgth - 1
        new_xtlab = [x.xtlab[position] for position in xtlab_positions]

    return STSList(
        data=newseq,
        alphabet=new_alph,
        labels=(
            None
            if label_by_state is None
            else [label_by_state.get(state, state) for state in new_alph]
        ),
        cpal=(
            None
            if color_by_state is None
            else [
                color_by_state.get(state, "white" if state == x.void else x.missing_color)
                for state in new_alph
            ]
        ),
        nr=new_nr,
        void=x.void,
        nr_code=nr_code,
        void_code=void_code,
        weights=None if x.weights is None else x.weights.copy(),
        start=x.start,
        xtlab=new_xtlab,
        xtstep=math.ceil(x.xtstep / tspan),
        tick_last=x.tick_last,
        ids=None if x.ids is None else list(x.ids),
        missing_color=x.missing_color,
        issues=x.issues,
    )


def _pc_minfreq(x: np.ndarray, level: float) -> np.ndarray:
    return (np.asarray(x, dtype=float) >= level).astype(float)


def _pc_cumfreq(x: np.ndarray, level: float) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    order = np.argsort(-x, kind="stable")  # decreasing (R order(x, decreasing=TRUE))
    cs = np.cumsum(x[order])
    reached = np.nonzero(cs >= level)[0]
    tmp = (int(reached[0]) + 1) if reached.size else x.size  # R which(...)[1]
    ret = np.zeros(x.size, dtype=float)
    ret[order[:tmp]] = 1.0
    return ret


def _pc_linear(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    return (x - x.min()) / (x.max() - x.min())


@dataclass
class SeqPcFilter:
    """A parallel-coordinate line filter (R ``seqpcfilter`` / ``construct.filter``).

    Wraps a colouring rule -- ``minfreq`` (colour lines with frequency >= ``level``),
    ``cumfreq`` (colour the most frequent lines up to cumulative ``level``) or
    ``linear`` (rescale to [0, 1]). Call :meth:`apply` to evaluate the rule on a vector
    of line frequencies. The ``seqpcplot`` rendering that consumes the filter is RP7.
    """

    type: str
    value: str
    level: Optional[float] = None

    def apply(self, x) -> np.ndarray:
        level = self.level
        if self.value == "minfreq":
            assert level is not None
            return _pc_minfreq(x, level)
        if self.value == "cumfreq":
            assert level is not None
            return _pc_cumfreq(x, level)
        return _pc_linear(x)


def seqpcfilter(method: str = "minfreq", level: Optional[float] = 0.05) -> SeqPcFilter:
    """Construct a parallel-coordinate line filter (TraMineR ``seqpcfilter``).

    ``method`` is one of ``"minfreq"``, ``"cumfreq"`` or ``"linear"``, or a unique
    prefix of one (``"cum"``), as R's ``match.arg`` takes it. Returns a
    :class:`SeqPcFilter`; the ``seqpcplot`` plot that applies it is deferred to RP7.
    Provenance: TraMineR 2.2-14 ``seqpcplot.R`` (``seqpcfilter`` / ``construct.filter``).

    Note: ``cumfreq`` with a ``level`` above the maximum cumulative weight (unreachable
    for real proportions in ``[0, 1]``) makes R abort (``which(...)[1]`` -> ``NA`` ->
    ``order[1:NA]``); traminex selects all lines instead. Documented divergence.
    """
    # R's match.arg(method) (seqpcplot.R:1100): an exact name, else the one name it
    # begins, case-sensitive. The level check reads the matched name, where R reads
    # the given one and lets ("cum", NULL) through (TX-AUD-20261011-002).
    from ._rformat import r_match_arg

    matched = r_match_arg(method, ("minfreq", "cumfreq", "linear"))
    if matched is None:
        raise ValueError("method must be one of minfreq, cumfreq, linear")
    method = matched
    if level is None and method in ("minfreq", "cumfreq"):
        raise ValueError("'seqpcfilter' requires a valid 'level' value.")
    return SeqPcFilter(type="function", value=method, level=level)


@_require_stslist
def seqauto(x: STSList, *, order: int = 1, measure: str = "cv") -> np.ndarray:
    """Sequence (state) autocorrelation by lag (TraMineRextras ``seqauto``).

    Builds the shifted-replication form (``STS_to_SRS``) and, for each lag
    ``i = 1..order``, cross-tabulates the current state against the state ``i`` steps
    back and measures their association with Cramer's V (``measure="cv"``):
    ``sqrt(chi2 / (N * (min(dim) - 1)))`` from ``chisq.test(correct=FALSE)``. Returns an
    ``order x 2`` array ``[auto-association, p-value]``. Provenance: TraMineRextras 0.6.9
    ``seqauto.R``.

    As in R, the table keeps every state either column takes (an all-zero row or
    column makes both values NaN), a one-row or one-column table gets
    ``chisq.test``'s goodness-of-fit p-value and an infinite (or NaN) Cramer's V,
    and a 1 x 1 table is refused.
    """
    import pandas as pd
    from scipy.stats import chi2 as chi2_distribution
    from scipy.stats import chi2_contingency

    from .format import STS_to_SRS

    if measure != "cv":
        raise ValueError(f"Measure {measure} not available, should be one of cv")
    srs = STS_to_SRS(x)
    out = np.full((order, 2), np.nan, dtype=float)
    for i in range(1, order + 1):
        current, lagged = srs["T"], srs[f"T-{i}"]
        # R's table(srs[, nc], srs[, nc - i]) keeps every level either column
        # takes, even one seen only beside an NA in the other column; such a
        # level is an all-zero row or column (audit S3-13).
        sub = srs[["T", f"T-{i}"]].dropna()
        ct = pd.crosstab(sub["T"], sub[f"T-{i}"]).reindex(
            index=pd.unique(current.dropna()), columns=pd.unique(lagged.dropna()), fill_value=0
        )
        table = ct.to_numpy(dtype=float)
        total = float(table.sum())
        if total == 0:
            raise ValueError("at least one entry of 'x' must be positive")
        if min(table.shape) == 1:
            # chisq.test turns a one-row or one-column table into a
            # goodness-of-fit test against equal probabilities.
            counts = table.ravel()
            if counts.size == 1:
                raise ValueError("'x' must at least have 2 elements")
            expected = total / counts.size
            chi2 = float(np.sum((counts - expected) ** 2 / expected))
            p = float(chi2_distribution.sf(chi2, counts.size - 1))
        elif np.any(table.sum(axis=0) == 0) or np.any(table.sum(axis=1) == 0):
            # A zero margin makes chisq.test's (x - E)^2 / E a 0/0: NaN.
            chi2 = p = float("nan")
        else:
            chi2, p, _, _ = chi2_contingency(table, correction=False)
        with np.errstate(divide="ignore", invalid="ignore"):
            out[i - 1, 0] = np.sqrt(np.float64(chi2) / np.float64(total * (min(table.shape) - 1)))
        out[i - 1, 1] = p
    return out


@_require_stslist
def seqcta(
    x: STSList,
    *,
    subseq: int = 5,
    time=None,
    event=None,
    initial_state: Optional[str] = None,
    covar=None,
):
    """Competing-trajectory-analysis converter (TraMineRextras ``seqcta``).

    Recodes the sequences (via the shifted-replication form ``STS_to_SRS``) into a
    person-period table suited to a discrete-time competing-risk model: for each id and
    each retained time the last ``subseq`` states form ``T1..Tsubseq`` (filled only at
    the event time), with an ``event`` flag and a ``lastobs`` flag. The transition time
    / event are taken from ``time``/``event`` or, if absent, derived from ``seqdur`` and
    ``initial_state``; the time is right-censored at ``seqlength - subseq``. Returns a
    DataFrame ``id, time, event, lastobs, T1..Tsubseq`` (optionally with ``covar`` by id).
    Provenance: TraMineRextras 0.6.9 ``seqcta.R``.
    """
    x = x._canonicalized()
    if not isinstance(subseq, (int, np.integer)) or isinstance(subseq, bool) or subseq < 1:
        raise ValueError("subseq must be a positive integer")
    if subseq > x.seqlen:
        # R crashes on the same input with an unrelated colnames-length error;
        # the recode needs subseq trailing positions per sequence.
        raise ValueError(f"subseq ({int(subseq)}) cannot exceed the sequence length ({x.seqlen})")
    import pandas as pd

    from .format import STS_to_SRS
    from .stats_complexity import seqdur

    n = x.nseq
    L = x.seqlen
    if x.void_code is not None:
        sl = (x.data != x.void_code).sum(axis=1).astype(float)
    else:
        sl = np.full(n, float(L))

    if time is not None:
        timevent = np.asarray(time, dtype=float).copy()
        ev = np.asarray(event, dtype=bool).copy() if event is not None else ~np.isnan(timevent)
    else:
        if initial_state is None or initial_state not in x.alphabet:
            raise ValueError(
                " [!] If time and event are not provided, initial_state should be one state "
                "of the alphabet of seqdata."
            )
        dur0 = seqdur(x)[:, 0]
        timevent = np.where(np.isnan(dur0), 0.0, dur0).astype(float)
        init_code = list(x.alphabet).index(initial_state)
        timevent[x.data[:, 0] != init_code] = 0.0
        ev = timevent < sl

    maxtime = sl - subseq
    censored = timevent > maxtime
    timevent[censored] = maxtime[censored]
    ev[censored] = False

    srs = STS_to_SRS(x)
    state_cols = [c for c in srs.columns if c not in ("id", "idx")]
    last_cols = state_cols[-subseq:]
    df = srs[["id", "idx"] + last_cols].copy()

    ids0 = df["id"].to_numpy(dtype=int) - 1
    idx_arr = df["idx"].to_numpy(dtype=float)
    timelimit = timevent[ids0] + subseq
    ev_row = ev[ids0] & (idx_arr == timelimit)
    keep = idx_arr <= timelimit
    df = df[keep].reset_index(drop=True)
    ev_row = ev_row[keep]

    tcols = [f"T{k}" for k in range(1, subseq + 1)]
    df.columns = ["id", "time"] + tcols
    df["event"] = ev_row
    df.loc[~df["event"], tcols] = np.nan
    df["time"] = df["time"] - subseq + 1
    df = df[df["time"] > 0].reset_index(drop=True)
    df["lastobs"] = df.groupby("id")["time"].transform("max") == df["time"]
    df = df[["id", "time", "event", "lastobs"] + tcols]
    if covar is not None:
        cov = pd.DataFrame(covar).reset_index(drop=True)
        df = pd.concat(
            [
                df.reset_index(drop=True),
                cov.iloc[df["id"].to_numpy(int) - 1].reset_index(drop=True),
            ],
            axis=1,
        )
    return df


def _weighted_km(
    dur: np.ndarray,
    status: np.ndarray,
    w: np.ndarray,
    *,
    confidence: bool = False,
    robust: bool | None = None,
):
    """Weighted Kaplan-Meier curve at event times (matches R ``survfit`` + ``summary``).

    Returns five-column summary rows, or seven columns including ``lower`` and
    ``upper`` when ``confidence=True``, only at event times.  ``robust`` can
    transport survfit's whole-model fractional-weight decision into a stratum.

    ``std.err`` follows R's ``survival:::survfitKM`` estimator selection
    (``has.rwt <- any(weights != floor(weights))``):

    - **integer weights** (incl. the unweighted all-ones case): the Greenwood
      log-scale estimator ``surv * sqrt(cumsum(d / (n*(n-d))))`` -- which yields ``NaN``
      at ``surv == 0`` (``n == d`` -> division by zero), exactly as R;
    - **non-integer weights**: the *infinitesimal-jackknife* (robust) SE, since R sets
      ``robust=TRUE`` for fractional weights. With ``lambda_k = d_k/n_k``,
      ``Y_i(t_k)=w_i[dur_i>=t_k]``, ``dN_i(t_k)=w_i[dur_i=t_k & event_i]``, the influence
      of observation ``i`` is ``U_i(t_m) = -S(t_m) * sum_{k<=m, n_k>d_k}
      (dN_i - Y_i*lambda_k)/(n_k - d_k)`` and ``std.err(t_m) = sqrt(sum_i U_i^2)`` -- which
      yields ``0`` at ``surv == 0`` (the ``n_k>d_k`` guard + ``S=0``), again as R.
    """
    rows: list[tuple[float, ...]] = []
    if not status.any():
        return rows
    has_rwt = bool(np.any(w != np.floor(w))) if robust is None else bool(robust)
    s = 1.0
    var = 0.0  # Greenwood cumulative
    infl = np.zeros(dur.shape[0], dtype=float)  # IJ per-observation cumulative term
    with np.errstate(divide="ignore", invalid="ignore"):
        for t in np.unique(dur[status]):
            at_risk_mask = dur >= t
            n = w[at_risk_mask].sum()
            d = w[(dur == t) & status].sum()
            if d == 0:
                # Only zero-weight spells end here: summary.survfit keeps rows
                # with n.event > 0 only, and the estimate is unchanged (audit S3-14).
                continue
            s *= 1.0 - d / n
            if has_rwt:
                if n > d:
                    dN = np.where((dur == t) & status, w, 0.0)
                    Y = np.where(at_risk_mask, w, 0.0)
                    infl += (dN - Y * (d / n)) / (n - d)
                se = float(s * np.sqrt(np.sum(infl * infl)))
            else:
                var += d / (n * (n - d))
                se = float(s * np.sqrt(var))
            if s == 0:
                # survival 3.8.9 stores Greenwood uncertainty on the log scale
                # (Inf at a terminal event), but robust fractional-weight fits
                # store a survival-scale zero.  survfit_confint therefore emits
                # NA/NA for the former and the exact boundary 0/0 for the latter.
                lower = upper = 0.0 if has_rwt and se == 0 else float("nan")
            elif se == 0:
                lower = upper = float(s)
            else:
                z95 = 1.959963984540054
                log_se = se / s
                lower = float(np.exp(np.log(s) - z95 * log_se))
                upper = float(min(np.exp(np.log(s) + z95 * log_se), 1.0))
            core = (float(t), float(n), float(d), float(s), se)
            rows.append(core + (lower, upper) if confidence else core)
    return rows


@dataclass
class SeqSurv:
    """Spell-survival curves (R ``seqsurv`` / ``stslist.surv``; the numeric table).

    ``curves`` is insertion ordered and maps a state or group stratum to an
    ``(m, 7)`` array of ``time, n.risk, n.event, surv, std.err, lower, upper``
    rows.  New metadata fields have compatibility defaults so historical direct
    ``SeqSurv(curves, states)`` construction remains valid.
    """

    curves: Dict[str, np.ndarray]
    states: List[str]
    curve_axis: str = "state"
    strata: List[str] = field(default_factory=list)
    groups: List[str] = field(default_factory=lambda: ["1"])
    ltext: List[str] = field(default_factory=list)
    cpal: List[str] | None = None
    xtstep: int = 1
    tick_last: bool = False
    max_time: int = 0
    # Last follow-up time of each stratum, censored or not: R's plot.survfit draws
    # a curve on to it (TX-AUD-20261010-042, TX-AUD-20261011-026).
    follow_up: Dict[str, float] = field(default_factory=dict)

    def summary(self):
        import pandas as pd

        columns = [
            "states",
            "time",
            "n.risk",
            "n.event",
            "surv",
            "std.err",
            "lower",
            "upper",
            "curve_axis",
            "state",
            "group",
        ]
        parts = []
        strata = self.strata or list(self.curves)
        for stratum in strata:
            arr = self.curves.get(stratum)
            if arr is None or len(arr) == 0:
                continue
            values = np.asarray(arr, dtype=float)
            if values.ndim != 2 or values.shape[1] not in {5, 7}:
                raise ValueError("SeqSurv curves must have five or seven columns")
            if values.shape[1] == 5:
                values = np.column_stack(
                    [values, np.full((values.shape[0], 2), np.nan, dtype=float)]
                )
            state_value = stratum if self.curve_axis == "state" else None
            group_value = (
                (self.groups[0] if self.groups else "1") if self.curve_axis == "state" else stratum
            )
            parts.append(
                pd.DataFrame(
                    {
                        "states": stratum,
                        "time": values[:, 0],
                        "n.risk": values[:, 1],
                        "n.event": values[:, 2],
                        "surv": values[:, 3],
                        "std.err": values[:, 4],
                        "lower": values[:, 5],
                        "upper": values[:, 6],
                        "curve_axis": self.curve_axis,
                        "state": state_value,
                        "group": group_value,
                    }
                )
            )
        return (
            pd.concat(parts, ignore_index=True).loc[:, columns]
            if parts
            else pd.DataFrame(columns=columns)
        )


_DARK2 = [
    "#1B9E77",
    "#D95F02",
    "#7570B3",
    "#E7298A",
    "#66A61E",
    "#E6AB02",
    "#A6761D",
    "#666666",
]


def _seqsurv_factor(groups, nseq: int) -> tuple[np.ndarray, List[str]]:
    """R ``factor(groups)`` values and used level order for seqsurv.

    Labels are R's ``as.character`` (``2``, ``TRUE``) and character levels
    follow R's collation (``a``, ``A``, ``b``), as ``factor()`` gives them
    (audit S3-22).
    """
    import pandas as pd

    from ._categorical_order import r_character_order
    from ._rformat import r_as_character

    if groups is None:
        return np.full(nseq, "1", dtype=object), ["1"]
    if len(groups) != nseq:
        raise ValueError(f"groups must have length {nseq}")

    def label(value: object) -> str:
        text = r_as_character(value)
        return "NA" if text is None else text

    series = pd.Series(groups)
    nonmissing = [value for value in series.tolist() if not pd.isna(value)]
    if isinstance(series.dtype, pd.CategoricalDtype):
        used = set(nonmissing)
        levels = [label(value) for value in series.cat.categories if value in used]
    else:
        unique = list(dict.fromkeys(nonmissing))
        if unique and all(
            isinstance(value, (int, float, np.integer, np.floating))
            and not isinstance(value, (bool, np.bool_))
            for value in unique
        ):
            levels = [label(value) for value in sorted(unique, key=float)]
        elif unique and all(isinstance(value, (bool, np.bool_)) for value in unique):
            levels = [label(value) for value in sorted(unique)]
        else:
            levels = r_character_order([label(value) for value in unique])

    levels = list(dict.fromkeys(levels))
    values = np.asarray(
        [None if pd.isna(value) else label(value) for value in series.tolist()],
        dtype=object,
    )
    return values, levels


def _seqsurv_group_palette(level_count: int) -> List[str] | None:
    if level_count == 1:
        return [_DARK2[2]]
    if level_count == 2:
        return _DARK2[:2]
    if 3 <= level_count <= 8:
        return _DARK2[:level_count]
    return None


def seqsurv(
    x: STSList,
    *,
    groups: Sequence[object] | None = None,
    per_state: bool = False,
    state: Sequence[str] | str | None = None,
    with_missing: bool = False,
) -> SeqSurv:
    """Weighted spell-survival curves by state or factor group.

    With ``per_state=False``, exactly one group is supported and one curve is
    fitted per selected observed state.  With ``per_state=True``, spells are
    filtered to the selected state set and one curve is fitted per group level.
    Curves reproduce survival 3.8.9's default 95% log confidence intervals.

    Unlike upstream TraMineRextras 0.6.9, ``with_missing=True`` is implemented:
    the declared missing state participates in spell and state routes without
    emitting upstream's misleading "isn't implemented yet" message.
    """
    from .format import STS_to_SPELL

    if not isinstance(x, STSList):
        raise TypeError("seqdata isn't a sequence object, use seqdef to create one")
    x = x._canonicalized()

    group_values, group_levels = _seqsurv_factor(groups, x.nseq)
    if not per_state and len(group_levels) != 1:
        raise ValueError(
            "With per_state=False, only a single group is supported. "
            "Consider using seqsplot() instead."
        )

    marker_codes = {code for code in (x.nr_code, x.void_code) if code is not None}
    observed_codes = {int(code) for code in np.asarray(x.data).ravel()}
    available: List[str] = [
        str(label)
        for code, label in enumerate(x.alphabet)
        if code in observed_codes and code not in marker_codes
    ]
    if with_missing and x.nr_code is not None:
        available.append(str(x.nr))

    if state is None:
        selected = list(available)
    else:
        requested = [state] if isinstance(state, str) else list(state)
        requested_set = {str(value) for value in requested}
        selected = [label for label in available if label in requested_set]
        if not selected:
            raise ValueError(" [!] state contains only unobserved states")

    # R's seqsurv builds its spells with seqformat(to="SPELL") (seqsurv.R:15), which
    # re-runs seqdef(right="DEL") first.
    spell = STS_to_SPELL(x, with_missing=with_missing, right="DEL")
    ids = spell["id"].to_numpy(dtype=int)
    begin = spell["begin"].to_numpy(dtype=int)
    end = spell["end"].to_numpy(dtype=int)
    dur = (end - begin + 1).astype(float)
    if x.void_code is not None:
        seqlen = (x.data != x.void_code).sum(axis=1)
    else:
        seqlen = np.full(x.nseq, x.seqlen, dtype=int)
    length = seqlen[ids - 1]
    status = end != length
    weights = (
        np.asarray(x.weights, dtype=float)[ids - 1]
        if x.weights is not None
        else np.ones(len(ids), dtype=float)
    )
    spell_states = spell["states"].astype(str).to_numpy(dtype=object)
    spell_groups = group_values[ids - 1]
    selected_mask = np.isin(spell_states, np.asarray(selected, dtype=object))
    model_weights = weights[selected_mask]
    robust = bool(np.any(model_weights != np.floor(model_weights)))

    curves: Dict[str, np.ndarray] = {}
    follow_up: Dict[str, float] = {}
    if per_state:
        for level in group_levels:
            mask = selected_mask & (spell_groups == level)
            if np.any(mask):
                follow_up[level] = float(np.max(dur[mask]))
            rows = _weighted_km(
                dur[mask],
                status[mask],
                weights[mask],
                confidence=True,
                robust=robust,
            )
            curves[level] = np.asarray(rows, dtype=float) if rows else np.empty((0, 7), dtype=float)
        curve_axis = "group"
        ltext = list(group_levels)
        palette = _seqsurv_group_palette(len(group_levels))
    else:
        for label in selected:
            mask = spell_states == label
            if np.any(mask):
                follow_up[label] = float(np.max(dur[mask]))
            rows = _weighted_km(
                dur[mask],
                status[mask],
                weights[mask],
                confidence=True,
                robust=robust,
            )
            curves[label] = np.asarray(rows, dtype=float) if rows else np.empty((0, 7), dtype=float)
        curve_axis = "state"
        code_by_label = {str(label): code for code, label in enumerate(x.alphabet)}
        ltext = [
            (
                str(x.labels[code_by_label[label]])
                if x.labels is not None and code_by_label[label] < len(x.labels)
                else label
            )
            for label in selected
        ]
        # A palette shorter than the alphabet (which STSList accepts) is left
        # to the renderers' shared warn-and-fall-back rule instead of raising
        # IndexError here (2026-09-22 audit S6-14).
        palette = (
            [str(x.cpal[code_by_label[label]]) for label in selected]
            if x.cpal is not None and len(x.cpal) >= len(x.alphabet)
            else None
        )

    return SeqSurv(
        curves=curves,
        states=selected,
        curve_axis=curve_axis,
        strata=list(curves),
        groups=group_levels,
        ltext=ltext,
        cpal=palette,
        xtstep=int(x.xtstep),
        tick_last=bool(x.tick_last),
        max_time=int(np.max(dur[selected_mask])) if np.any(selected_mask) else 0,
        follow_up=follow_up,
    )


def _stack_role_codes(
    roles: Sequence[STSList],
) -> tuple[np.ndarray, list[str], int | None, int | None]:
    """Recode the seqpolyads roles onto one codebook and stack their rows.

    The codebook is the first role's real states, then its nr and void texts when
    any role declares that marker. Every role's real states must be among the
    first role's (the caller checks R's alphabet comparison first).
    """
    from .sequences import alphabet as sequence_alphabet

    real = [str(state) for state in sequence_alphabet(roles[0])]
    use_nr = any(role.nr_code is not None for role in roles)
    use_void = any(role.void_code is not None for role in roles)
    alphabet = list(real)
    nr_code = void_code = None
    if use_nr:
        nr_code = len(alphabet)
        alphabet.append(roles[0].nr)
    if use_void:
        void_code = len(alphabet)
        alphabet.append(roles[0].void)
    blocks = []
    for role in roles:
        mapping = np.empty(len(role.alphabet), dtype=int)
        for code, state in enumerate(role.alphabet):
            if code == role.nr_code:
                mapping[code] = cast(int, nr_code)
            elif code == role.void_code:
                mapping[code] = cast(int, void_code)
            else:
                mapping[code] = real.index(str(state))
        blocks.append(mapping[np.asarray(role.data, dtype=int)])
    return np.vstack(blocks), alphabet, nr_code, void_code


def _per_state(role: STSList, alphabet: list[str], field: str) -> list[str] | None:
    """The first role's labels or colours re-indexed onto the stacked codebook."""
    values = getattr(role, field)
    if values is None:
        return None
    canonical = role._canonicalized()
    by_state = dict(zip(canonical.alphabet, getattr(canonical, field)))
    fallback = {"labels": None, "cpal": role.missing_color}[field]
    return [by_state.get(state, state if fallback is None else fallback) for state in alphabet]


def seqpolyads(
    seqlist: Sequence[STSList],
    *,
    a: int = 1,
    method: str = "HAM",
    pair_weights: Sequence[float] | None = None,
    s: int = 36963,
    T: int = 1000,
    core: int | Literal["auto"] = 1,
    replace: bool = True,
    weighted: bool = True,
    with_missing: bool = False,
    rand_weight_type: int = 1,
    role_weights: Sequence[float] | None = None,
    show_time: bool = False,
    **distance_kwargs: object,
) -> Dict[str, np.ndarray]:
    """Linked-polyad sequence analysis (TraMineRextras ``seqpolyads``).

    ``a=1`` independently samples one observed sequence for each role; ``a=2``
    additionally resamples positions within each selected sequence.  All random
    index plans are generated serially by one local NumPy generator before any
    distance work, so ``core`` is an execution-only control and cannot alter the
    stream.  Unlike upstream's overwritten-argument defect, ``s`` is the actual
    Python seed (an intentional ``DIVERGENCE/FULL`` contract).
    """
    import warnings
    from math import comb
    from time import perf_counter

    from .distances import seqdist
    from .format import seqdef
    from .sequences import alphabet as sequence_alphabet

    started = perf_counter()
    try:
        roles = list(seqlist)
    except TypeError:
        raise ValueError("seqlist must be a list of at least two stslist objects") from None
    P = len(roles)
    if P < 2:
        raise ValueError("seqlist must be a list of at least two stslist objects")
    if any(not isinstance(role, STSList) for role in roles):
        raise TypeError("At least one element of seqlist is not a stslist object")
    roles = [role._canonicalized() for role in roles]
    if a not in {1, 2}:
        raise ValueError("Bad 'a' value")
    if core != "auto" and (
        not isinstance(core, (int, np.integer)) or isinstance(core, (bool, np.bool_))
    ):
        raise ValueError('core must be numeric or "auto"')
    if not isinstance(T, (int, np.integer)) or isinstance(T, (bool, np.bool_)) or T < 2:
        raise ValueError("T must be an integer greater than one")
    if not isinstance(s, (int, np.integer)) or isinstance(s, (bool, np.bool_)):
        raise TypeError("s must be an integer seed")
    if rand_weight_type not in {1, 2}:
        raise ValueError("rand_weight_type must be 1 or 2")

    n = roles[0].nseq
    ncol = roles[0].seqlen
    alph0 = sequence_alphabet(roles[0], with_missing=with_missing)
    for p in range(1, P):
        if (
            roles[p].nseq != n
            or roles[p].seqlen != ncol
            or sequence_alphabet(roles[p], with_missing=with_missing) != alph0
        ):
            raise ValueError("All stslist objects must have same size and alphabet")

    d = comb(P, 2)
    pair_weight_array = (
        np.ones(d, dtype=float) if pair_weights is None else np.asarray(pair_weights, dtype=float)
    )
    if pair_weight_array.ndim != 1 or pair_weight_array.size != d:
        raise ValueError(f"pair_weights must have length {d}")
    if role_weights is None:
        role_weight_array = np.full(P, 1.0 / P, dtype=float)
    else:
        role_weight_array = np.asarray(role_weights, dtype=float)
        if role_weight_array.ndim != 1 or role_weight_array.size != P:
            raise ValueError("length of role_weights not equal to length of seqlist")

    # Stack the P role objects into one sequence object. R compares only the roles'
    # alphabet(seqdata, with.missing) (seqpolyads.R:52-57), so their codebooks may
    # differ (a nr code in one role only): recode every role onto one codebook.
    data, stacked_alphabet, stacked_nr, stacked_void = _stack_role_codes(roles)
    if all(role.weights is not None for role in roles):
        all_weights = np.concatenate([np.asarray(role.weights, dtype=float) for role in roles])
    else:
        all_weights = np.ones(n * P, dtype=float)
    seqall = STSList(
        data=data,
        alphabet=stacked_alphabet,
        labels=_per_state(roles[0], stacked_alphabet, "labels"),
        cpal=_per_state(roles[0], stacked_alphabet, "cpal"),
        nr=roles[0].nr,
        void=roles[0].void,
        nr_code=stacked_nr,
        void_code=stacked_void,
        weights=all_weights,
        ids=None,
        missing_color=roles[0].missing_color,
    )
    cj = np.arange(P) * n
    iu = np.triu_indices(P, k=1)

    weights = all_weights[:n] if weighted else np.ones(n, dtype=float)

    # Generate the complete Monte-Carlo plan before either distance call.  The
    # worker/execution control never owns an RNG and therefore cannot change it.
    rng = np.random.default_rng(int(s))
    probability = weights / weights.sum()
    sampled_local = rng.choice(n, size=(T, P), replace=True, p=probability)
    sampled_global = sampled_local + cj[np.newaxis, :]
    if a == 2:
        if replace:
            position_plan = rng.integers(0, ncol, size=(T, P, ncol))
        else:
            position_plan = np.empty((T, P, ncol), dtype=int)
            for iteration in range(T):
                for role in range(P):
                    position_plan[iteration, role] = rng.permutation(ncol)

    alldist = np.asarray(
        seqdist(
            seqall,
            method=method,
            with_missing=with_missing,
            weighted=weighted,
            **distance_kwargs,
        ),
        dtype=float,
    )

    observed = np.zeros(n, dtype=float)
    for j in range(n):
        idxs = cj + j
        sub = alldist[np.ix_(idxs, idxs)]
        observed[j] = float((sub[iu] * pair_weight_array).sum() / d)
    mean_obs = float((observed * weights).sum() / weights.sum())

    random_dist = np.empty(T, dtype=float)
    if a == 1:
        for iteration in range(T):
            idxs = sampled_global[iteration]
            sub = alldist[np.ix_(idxs, idxs)]
            random_dist[iteration] = float((sub[iu] * pair_weight_array).sum() / d)
    else:
        sampled_rows = data[sampled_global]
        randomized_codes = np.take_along_axis(sampled_rows, position_plan, axis=2)
        raw = np.empty(randomized_codes.shape, dtype=object)
        labels = np.asarray(seqall.alphabet, dtype=object)
        raw[:] = labels[randomized_codes]
        special = np.zeros(randomized_codes.shape, dtype=bool)
        if seqall.nr_code is not None:
            special |= randomized_codes == int(seqall.nr_code)
        has_void = seqall.void_code is not None and bool(
            np.any(seqall.data == int(seqall.void_code))
        )
        if seqall.void_code is not None:
            special |= randomized_codes == int(seqall.void_code)
        raw[special] = None
        if has_void:
            warnings.warn(
                "Found voids in input sequences: dropping all right missings "
                "in the resampled sequences!",
                UserWarning,
                stacklevel=2,
            )
        real_alphabet = sequence_alphabet(seqall)
        resampled_nr = _collision_free_marker(seqall.nr, [*real_alphabet, seqall.void], "NR")
        resampled_void = _collision_free_marker(seqall.void, [*real_alphabet, resampled_nr], "VOID")
        randomized = seqdef(
            raw.reshape(T * P, ncol),
            alphabet=real_alphabet,
            nr=resampled_nr,
            void=resampled_void,
            right="DEL" if has_void else None,
        )
        random_matrix = np.asarray(
            seqdist(
                randomized,
                method=method,
                with_missing=with_missing,
                **distance_kwargs,
            ),
            dtype=float,
        )
        role_offsets = np.arange(P)
        for iteration in range(T):
            idxs = iteration * P + role_offsets
            sub = random_matrix[np.ix_(idxs, idxs)]
            random_dist[iteration] = float((sub[iu] * pair_weight_array).sum() / d)

    l_weights = np.ones(T, dtype=float)
    if weighted:
        member_weights = all_weights[sampled_global]
        if rand_weight_type == 2:
            role_denominators = member_weights.sum(axis=0)
            member_weights = np.divide(
                member_weights,
                role_denominators[np.newaxis, :],
                out=np.zeros_like(member_weights),
                where=role_denominators[np.newaxis, :] != 0,
            )
        l_weights = np.sum(member_weights * role_weight_array[np.newaxis, :], axis=1)
    l_weights = l_weights / l_weights.sum()

    test_p = np.array([float(((observed[j] < random_dist) * l_weights).sum()) for j in range(n)])
    v95 = (test_p > 0.95).astype(float)
    mean_rand = float((random_dist * l_weights).sum())
    mean_u = mean_rand - observed
    from scipy.stats import t as _t

    se = random_dist.std(ddof=1) / np.sqrt(T)
    u_p = 2 * _t.sf(np.abs(mean_u / se), T - 1) if se > 0 else np.full(n, np.nan)

    result = {
        "mean.dist": np.array([mean_obs, mean_rand], dtype=float),
        "U": mean_u,
        "U.tp": np.asarray(u_p, dtype=float),
        "V": test_p,
        "V.95": v95,
        "observed.dist": observed,
        "random.dist": random_dist,
    }
    if show_time:
        print(f"seqpolyads elapsed: {perf_counter() - started:.6f} seconds")
    return result
