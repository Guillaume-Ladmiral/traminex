from __future__ import annotations

import decimal
import math
import os
import warnings
from collections.abc import Callable, Iterable, Iterator
from contextlib import contextmanager
from contextvars import ContextVar
from numbers import Real
from typing import Any, cast

import numpy as np

from .._rformat import r_number_text
from . import _batched
from ._normalization import _maxscost, normalize_value, omspell_bound, resolve_norm_code
from ._spells import spell_encode, spell_length_per_position, spell_length_power
from ._validation import (
    is_r_positive_integer,
    is_refused_cost_scalar,
    normalize_dhd_variant,
    validate_pairwise_kwargs,
    validate_smooth_window,
)

try:  # optional numba accelerators
    from numba import njit
    from numba import prange as _numba_prange

    # Numba's runtime ``prange`` object is callable but recent stubs model it as
    # a module.  Preserve the exact runtime object while giving static checkers
    # its effective return type.
    prange = cast(Callable[..., range], _numba_prange)

    _NUMBA_IMPORTED = True
except ImportError:  # pragma: no cover
    _NUMBA_IMPORTED = False

#: Accepted values of ``TRAMINEX_NUMBA``.
NUMBA_MODES = ("auto", "off")


def _resolve_numba_mode() -> str:
    """Read and validate ``TRAMINEX_NUMBA``.

    Numba is a **third** distance backend: it is active by default whenever
    ``numba`` imports, and it changes which kernel computes a result. Before this
    switch existed there was no supported way to select the reference kernels --
    ``TRAMINEX_BINDING`` chooses between the Python and compiled *implementations*
    and leaves Numba on under every one of its three values -- so the only way to
    run the suite against the reference path was to monkeypatch a private flag
    (cycle-4b F-008).

    ``TRAMINEX_NUMBA`` is deliberately a separate variable rather than a fourth
    ``TRAMINEX_BINDING`` value: the two select orthogonal things, and folding them
    into one three-valued domain would make ``python`` mean two different things --
    which is exactly the confusion F-015 found in the documentation.
    """
    mode = os.environ.get("TRAMINEX_NUMBA", "auto").strip().lower() or "auto"
    if mode not in NUMBA_MODES:
        raise ValueError(f"TRAMINEX_NUMBA must be one of: {', '.join(NUMBA_MODES)}")
    return mode


#: Whether the Numba kernels are *dispatched to*. Distinct from
#: :data:`_NUMBA_IMPORTED`, which records whether the ``@njit`` kernels could be
#: defined at all: the maintained benchmark CLI toggles this flag at runtime to
#: compare backends, so the kernels must exist even when dispatch is disabled.
#: The mode is read first, so an invalid value is refused whether or not numba
#: imports (TX-AUD-20261011-049).
_NUMBA_MODE = _resolve_numba_mode()
_HAVE_NUMBA = _NUMBA_IMPORTED and _NUMBA_MODE == "auto"


def _scalar_float(value) -> float:
    """Convert a value already established as scalar without widening array types."""
    return float(np.asarray(value).item())


def _validate_dense_sequence_matrix(seqmat) -> np.ndarray:
    """Validate the stable dense-dispatcher state-code contract."""
    matrix = np.asarray(seqmat)
    if matrix.ndim != 2:
        raise ValueError("seqmat must be a 2D integer-coded sequence matrix")
    if matrix.dtype.kind not in {"i", "u"}:
        raise TypeError("seqmat must contain integer state codes")
    if matrix.size and int(matrix.min()) < 0:
        raise ValueError("seqmat state codes must be non-negative")
    if matrix.size and int(matrix.max()) > np.iinfo(np.int32).max:
        # The kernels work in int32; a larger code silently wrapped (2**32 -> 0),
        # so different sequences compared equal (audit S1b-11).
        raise ValueError("seqmat state codes must fit a 32-bit integer (at most 2147483647)")
    return matrix


def _validate_substitution_costs(
    subm,
    *,
    seqmat: np.ndarray,
    method: str,
    allow_scalar: bool,
) -> float | np.ndarray:
    """Validate fixed substitution costs before any Python or native dispatch."""
    if np.isscalar(subm):
        if not allow_scalar:
            raise ValueError(f"{method} subm must be square")
        value = _scalar_float(subm)
        if not np.isfinite(value):
            raise ValueError(f"{method} substitution costs must be finite")
        if value < 0:
            raise ValueError(f"{method} substitution costs cannot be negative")
        return value

    costs = np.asarray(subm, dtype=float)
    if costs.ndim != 2 or costs.shape[0] != costs.shape[1]:
        raise ValueError(f"{method} subm must be square")
    observed_states = int(seqmat.max()) + 1 if seqmat.size else 0
    if costs.shape[0] < observed_states:
        raise ValueError(
            f"{method} subm must cover every observed state code "
            f"(need at least {observed_states}x{observed_states})"
        )
    if not _observed_substitution_costs_are_finite(costs, seqmat):
        raise ValueError(f"{method} substitution costs must be finite")
    if np.any(costs < 0):
        raise ValueError(f"{method} substitution costs cannot be negative")
    if costs.size and np.any(np.diag(costs) != 0):
        raise ValueError(f"{method} substitution-cost diagonal must be zero")
    return costs


def _observed_substitution_costs_are_finite(
    costs: np.ndarray,
    seqmat: np.ndarray,
    *,
    ignored_codes: tuple[int | None, ...] = (),
) -> bool:
    """Whether every cost reachable from the observed state codes is finite."""
    matrix = np.asarray(costs, dtype=float)
    if matrix.ndim == 0:
        return bool(np.isfinite(matrix.item()))

    states = np.unique(np.asarray(seqmat, dtype=int))
    ignored = {int(code) for code in ignored_codes if code is not None}
    states = states[
        np.asarray(
            [int(code) not in ignored and 0 <= int(code) < matrix.shape[0] for code in states],
            dtype=bool,
        )
    ]
    if states.size == 0:
        return True
    return bool(np.all(np.isfinite(matrix[np.ix_(states, states)])))


def _validate_om_position_costs(
    seqmat: np.ndarray,
    *,
    indel_by_pos,
    subm_by_pos,
) -> tuple[np.ndarray | None, np.ndarray | None]:
    """Validate the documented dense position-dependent OM extension."""
    length = seqmat.shape[1]
    observed_states = int(seqmat.max()) + 1 if seqmat.size else 0

    resolved_indel = None
    if indel_by_pos is not None:
        resolved_indel = np.asarray(indel_by_pos, dtype=float)
        if resolved_indel.ndim != 1 or resolved_indel.shape[0] != length:
            raise ValueError("indel_by_pos must be one-dimensional with sequence length")
        if not np.all(np.isfinite(resolved_indel)) or np.any(resolved_indel <= 0):
            raise ValueError("indel_by_pos costs must be finite and strictly positive")

    resolved_subm = None
    if subm_by_pos is not None:
        resolved_subm = np.asarray(subm_by_pos, dtype=float)
        expected = (length, observed_states, observed_states)
        if resolved_subm.shape != expected:
            raise ValueError(f"subm_by_pos must have shape {expected}")
        if not np.all(np.isfinite(resolved_subm)):
            raise ValueError("subm_by_pos substitution costs must be finite")
        if np.any(resolved_subm < 0):
            raise ValueError("subm_by_pos substitution costs cannot be negative")
        if resolved_subm.size and np.any(np.diagonal(resolved_subm, axis1=1, axis2=2) != 0):
            raise ValueError("subm_by_pos substitution-cost diagonals must be zero")

    return resolved_indel, resolved_subm


#: The row engine evaluates its dynamic programs for many pairs at once
#: (:mod:`traminex.distances._batched`, TX-AUD-20261010-022). ``False`` runs the
#: per-pair reference implementations instead; the tests compare the two.
_USE_BATCHED = True

#: ``(row count, left mask, right mask)``: boolean row masks of the only pairs a
#: caller needs, set by :func:`_only_pairs` around a reference (``refseq``)
#: computation.
_PAIR_FILTER: ContextVar[tuple[int, np.ndarray, np.ndarray] | None] = ContextVar(
    "traminex_pair_filter", default=None
)

#: Pairs are generated, evaluated and written in blocks of about this many, so
#: the bookkeeping stays bounded whatever the number of rows; listing every pair
#: as Python objects took about 14 times the result matrix (batch-2 review I-2).
_PAIR_BLOCK = 1 << 15


@contextmanager
def _only_pairs(n: int, left: Iterable[int], right: Iterable[int]) -> Iterator[None]:
    """Within the block, the row engine evaluates only pairs joining ``left`` and ``right``.

    A ``refseq`` needs one distance per sequence (or one block); the full matrix
    it is sliced from used to be computed in full (TX-AUD-20261010-022). Every
    other cell stays 0 and must be discarded by the caller.

    The filter applies to every row-engine call over exactly ``n`` rows inside
    the block, and reads that call's row positions as the caller's ``0..n-1``.
    A call over another number of rows (a deduplicated computation) evaluates
    every pair. A call over the same number of rows in another order would be
    filtered with the wrong positions, so the block must not contain one; the
    callers (``seqdist``'s ``refseq`` routes: an index or index sets over the
    data, an external reference over the data with the reference appended)
    keep the rows and their order.

    With asymmetric costs each pair is evaluated with its ``left`` member
    first and only the ``(left, right)`` cells are written, as R computes
    d(sequence, reference) and d(set-1 element, set-2 element)
    (TX-REP-20261017-001); symmetric costs fill both cells.
    """
    count = int(n)
    in_left = np.zeros(count, dtype=bool)
    in_right = np.zeros(count, dtype=bool)
    in_left[np.fromiter((int(i) for i in left), dtype=np.int64)] = True
    in_right[np.fromiter((int(j) for j in right), dtype=np.int64)] = True
    token = _PAIR_FILTER.set((count, in_left, in_right))
    try:
        yield
    finally:
        _PAIR_FILTER.reset(token)


def _pair_masks(n: int) -> tuple[np.ndarray, np.ndarray] | None:
    """The :func:`_only_pairs` row masks for a call over ``n`` rows, or None."""
    selected = _PAIR_FILTER.get()
    if selected is None or selected[0] != n:
        return None
    return selected[1], selected[2]


def _worth_batching(n: int, from_pairs: int) -> bool:
    """Whether a call over ``n`` rows evaluates at least ``from_pairs`` pairs.

    The batched sweep costs about as much for one pair as for a few dozen,
    while the per-pair reference costs one dynamic program per pair, so below a
    method's crossover the reference is faster (batch-2 review M-4); both give
    the same bits. Each caller passes its crossover, the ratio of the two times
    for one pair rounded up, measured at lengths 150 and 300 on Ubuntu 24.04
    with Python 3.12 and NumPy 2.3. The count honours :func:`_only_pairs` but
    not deduplication, so asymmetric costs over repeated rows may batch early.
    """
    masks = _pair_masks(n)
    if masks is None:
        return n * (n - 1) // 2 >= from_pairs
    in_left, in_right = masks
    count = 0
    for a in range(n - 1):
        b = np.arange(a + 1, n, dtype=np.int64)
        count += int(np.count_nonzero((in_left[a] & in_right[b]) | (in_right[a] & in_left[b])))
        if count >= from_pairs:
            return True
    return False


def _stable_groups(rows) -> tuple[list[list[int]], list[int]]:
    """Identical rows grouped in first-occurrence order, with each group's first row."""
    groups: list[list[int]] = []
    representatives: list[int] = []
    group_by_row: dict[tuple[int, ...], int] = {}
    for index, row in enumerate(rows):
        key = tuple(int(value) for value in np.asarray(row).tolist())
        group_index = group_by_row.get(key)
        if group_index is None:
            group_index = len(groups)
            group_by_row[key] = group_index
            groups.append([])
            representatives.append(index)
        groups[group_index].append(index)
    return groups, representatives


def _upper_pair_blocks(
    count: int, has_left: np.ndarray | None, has_right: np.ndarray | None
) -> Iterator[tuple[np.ndarray, np.ndarray]]:
    """Every ``a < b`` below ``count`` in row-major order, as int64 index-array blocks.

    With masks, only pairs with one side in ``has_left`` and the other in
    ``has_right`` are kept. Each block holds about :data:`_PAIR_BLOCK` pairs,
    so all ``count * (count - 1) / 2`` pairs are never held at once.
    """
    firsts: list[np.ndarray] = []
    seconds: list[np.ndarray] = []
    size = 0
    for a in range(count - 1):
        b = np.arange(a + 1, count, dtype=np.int64)
        if has_left is not None and has_right is not None:
            b = b[(has_left[a] & has_right[b]) | (has_right[a] & has_left[b])]
        if b.size == 0:
            continue
        firsts.append(np.full(b.size, a, dtype=np.int64))
        seconds.append(b)
        size += b.size
        if size >= _PAIR_BLOCK:
            yield np.concatenate(firsts), np.concatenate(seconds)
            firsts, seconds, size = [], [], 0
    if size:
        yield np.concatenate(firsts), np.concatenate(seconds)


def _pair_blocks(rows, asymmetric: bool):
    """The pairs R evaluates, in its order, filtered by :func:`_only_pairs`, in blocks.

    Yields ``(left, right, writer)``: int64 row indices of the pairs to evaluate
    and ``writer(distances, values)``, which stores their values. Symmetric
    costs: every ``i < j``. Asymmetric costs: the first rows of distinct groups
    in first-occurrence order (TraMineR deduplicates before its kernels), each
    value copied to every member pair of the two groups.
    """
    n = len(rows)
    masks = _pair_masks(n)
    if not asymmetric:
        in_left, in_right = masks if masks is not None else (None, None)
        for left, right in _upper_pair_blocks(n, in_left, in_right):
            yield left, right, _symmetric_writer(left, right)
        return
    groups, representatives = _stable_groups(rows)
    first_rows = np.asarray(representatives, dtype=np.int64)
    if masks is not None:
        # A reference distance is d(sequence, reference) and an index-set block
        # d(set-1 element, set-2 element), whichever comes first in the data; with
        # asymmetric costs the direction shows, so each pair is evaluated with its
        # ``left`` member first and only those cells are written
        # (TX-REP-20261017-001). Identical rows stay at 0.
        in_left, in_right = masks
        group_left = np.asarray([bool(in_left[members].any()) for members in groups])
        group_right = np.asarray([bool(in_right[members].any()) for members in groups])
        for a, b in _oriented_pair_blocks(len(groups), group_left, group_right):
            yield first_rows[a], first_rows[b], _oriented_group_writer(groups, masks, a, b)
        return
    has_left = has_right = None
    singletons = all(len(members) == 1 for members in groups)
    for a, b in _upper_pair_blocks(len(groups), has_left, has_right):
        left, right = first_rows[a], first_rows[b]
        if singletons:
            yield left, right, _symmetric_writer(left, right)
        else:
            yield left, right, _group_writer(groups, a, b)


def _oriented_pair_blocks(
    count: int, has_left: np.ndarray, has_right: np.ndarray
) -> Iterator[tuple[np.ndarray, np.ndarray]]:
    """Every ``(a, b)`` with ``a != b``, ``has_left[a]`` and ``has_right[b]``, in blocks."""
    firsts: list[np.ndarray] = []
    seconds: list[np.ndarray] = []
    size = 0
    candidates = np.flatnonzero(has_right).astype(np.int64)
    for a in np.flatnonzero(has_left).tolist():
        b = candidates[candidates != a]
        if b.size == 0:
            continue
        firsts.append(np.full(b.size, a, dtype=np.int64))
        seconds.append(b)
        size += b.size
        if size >= _PAIR_BLOCK:
            yield np.concatenate(firsts), np.concatenate(seconds)
            firsts, seconds, size = [], [], 0
    if size:
        yield np.concatenate(firsts), np.concatenate(seconds)


def _oriented_group_writer(groups, masks, a: np.ndarray, b: np.ndarray):
    """Store each value in the cells (left member of group a, right member of group b)."""
    in_left, in_right = masks

    def write(distances: np.ndarray, values: np.ndarray) -> None:
        for first, second, value in zip(a.tolist(), b.tolist(), values.tolist()):
            rows_in = [index for index in groups[first] if in_left[index]]
            cols_in = [index for index in groups[second] if in_right[index]]
            distances[np.ix_(rows_in, cols_in)] = value

    return write


def _symmetric_writer(left: np.ndarray, right: np.ndarray):
    def write(distances: np.ndarray, values: np.ndarray) -> None:
        distances[left, right] = values
        distances[right, left] = values

    return write


def _group_writer(groups: list[list[int]], a: np.ndarray, b: np.ndarray):
    def write(distances: np.ndarray, values: np.ndarray) -> None:
        for first, second, value in zip(a.tolist(), b.tolist(), values.tolist()):
            distances[np.ix_(groups[first], groups[second])] = value
            distances[np.ix_(groups[second], groups[first])] = value

    return write


def _pairwise_apply(n: int, fn):
    rows = range(n)  # distinct keys: plain i < j order
    distances = np.zeros((n, n), dtype=float)
    for left, right, write in _pair_blocks(rows, asymmetric=False):
        values = np.empty(left.size, dtype=float)
        for k, (i, j) in enumerate(zip(left.tolist(), right.tolist())):
            values[k] = fn(i, j)
        write(distances, values)
    return distances


def _pairwise_apply_stable_unique(rows, fn) -> np.ndarray:
    """Evaluate distinct rows in first-occurrence order and expand the result.

    TraMineR removes duplicated sequences before entering its distance kernels.
    First-occurrence order is observable only for accepted asymmetric costs.
    """
    n = len(rows)
    distances = np.zeros((n, n), dtype=float)
    for left, right, write in _pair_blocks(rows, asymmetric=True):
        values = np.empty(left.size, dtype=float)
        for k, (i, j) in enumerate(zip(left.tolist(), right.tolist())):
            values[k] = fn(i, j)
        write(distances, values)
    return distances


def _evaluate_raw_batches(left: np.ndarray, right: np.ndarray, raw_batch):
    """Raw distances of pairs ``(left[k], right[k])`` (one indirection for the tests)."""
    return raw_batch(left, right)


def _pairwise_apply_batched(
    rows,
    asymmetric: bool,
    raw_batch,
    finish: Callable[[int, int, Any], float],
    special: Callable[[int, int], float | None] | None = None,
) -> np.ndarray:
    """:func:`_pairwise_apply` (or the stable-unique order) with batched raw distances.

    ``special(i, j)`` answers pairs that need no dynamic program (or None);
    ``raw_batch(left, right)`` gives the raw distances of the others at once and
    ``finish(i, j, raw)`` turns each into the value the per-pair path returns.
    """
    n = len(rows)
    distances = np.zeros((n, n), dtype=float)
    for left, right, write in _pair_blocks(rows, asymmetric):
        firsts, seconds = left.tolist(), right.tolist()
        values = np.empty(left.size, dtype=float)
        pending = np.ones(left.size, dtype=bool)
        if special is not None:
            for k, (i, j) in enumerate(zip(firsts, seconds)):
                answer = special(i, j)
                if answer is not None:
                    values[k] = answer
                    pending[k] = False
        positions = np.flatnonzero(pending)
        if positions.size:
            raws = _evaluate_raw_batches(left[positions], right[positions], raw_batch)
            for k, raw in zip(positions.tolist(), raws):
                values[k] = finish(firsts[k], seconds[k], raw)
        write(distances, values)
    return distances


def _in_batches(raw, left: np.ndarray, right: np.ndarray, width: int, *, depth: int = 1):
    """Apply ``raw(left, right)`` to slices sized for ``width`` cells per pair."""
    size = _batched.batch_size(width, depth)
    parts = [
        raw(left[start : start + size], right[start : start + size])
        for start in range(0, len(left), size)
    ]
    if not parts:
        return np.empty(0, dtype=float)
    if isinstance(parts[0], list):
        return [vector for part in parts for vector in part]
    return np.concatenate(parts)


def _longest(rows) -> int:
    return max((len(row) for row in rows), default=0)


def _asymmetric_slice_flags(costs) -> np.ndarray:
    """Flag the square slices whose final two axes differ beyond R's tolerance.

    The tolerance is ``checkcost``'s ``1e-07 * max(sma)``, taken once over the
    whole matrix or array; a ``(K, K)`` matrix is one slice.
    """
    array = np.asarray(costs, dtype=float)
    if array.ndim < 2 or array.shape[-1] != array.shape[-2] or array.size == 0:
        return np.zeros(0, dtype=bool)
    finite_values = np.abs(array[np.isfinite(array)])
    maximum = float(np.max(finite_values)) if finite_values.size else 0.0
    tolerance = 1e-7 * maximum
    transposed = np.swapaxes(array, -1, -2)
    with np.errstate(invalid="ignore"):
        differences = np.abs(array - transposed)
    exceeded = (np.isfinite(array) != np.isfinite(transposed)) | (
        np.isfinite(differences) & (differences > tolerance)
    )
    size = array.shape[-1]
    return np.asarray(exceeded.reshape(-1, size, size).any(axis=(1, 2)), dtype=bool)


def _asymmetric_cost_slices(costs) -> int:
    """Count the slices :func:`_asymmetric_slice_flags` flags."""
    return int(np.count_nonzero(_asymmetric_slice_flags(costs)))


def _cost_array_is_asymmetric(costs) -> bool:
    """Return whether the final two axes differ beyond the existing R tolerance."""
    return _asymmetric_cost_slices(costs) > 0


def _warn_asymmetric_costs(costs, *, stacklevel: int) -> None:
    """Emit TraMineR's diagnostic for accepted asymmetric substitution costs.

    ``TraMineR-checkcost.R:115-116`` (2.2-14) warns once per asymmetric matrix,
    so once per asymmetric position of a time-varying array.
    """
    for _ in range(_asymmetric_cost_slices(costs)):
        warnings.warn("Non-symmetric substitution cost matrix!", UserWarning, stacklevel=stacklevel)


def _pairwise_apply_for_costs_batched(rows, subm, raw_batch, finish, special=None) -> np.ndarray:
    """:func:`_pairwise_apply_for_costs` with batched raw distances (same order, same warning)."""
    asymmetric = np.ndim(subm) == 2 and _cost_array_is_asymmetric(np.asarray(subm, dtype=float))
    if asymmetric:
        _warn_asymmetric_costs(subm, stacklevel=5)
    return _pairwise_apply_batched(rows, asymmetric, raw_batch, finish, special)


def _pairwise_apply_for_costs(rows, subm, fn) -> np.ndarray:
    """Evaluate pairs in the order R would for these costs.

    TraMineR deduplicates before its kernels, which only matters for asymmetric
    costs: then the first-occurrence order fixes each comparison's direction,
    so identical sequences get identical distances, and ``checkcost`` warns
    (2026-09-22 audit S1b-03). Symmetric costs keep the plain row order.
    """
    if np.ndim(subm) == 2 and _cost_array_is_asymmetric(np.asarray(subm, dtype=float)):
        _warn_asymmetric_costs(subm, stacklevel=5)
        return _pairwise_apply_stable_unique(rows, fn)
    return _pairwise_apply(len(rows), fn)


class CostDiagnosticMessage(UserWarning):
    """A cost diagnostic that R signals through ``message()`` rather than ``warning()``.

    Python has no separate catchable informational-condition channel, so R's
    ``message()`` diagnostics are surfaced as warnings.  Giving them their own
    category keeps R's channel distinction expressible: a caller -- or a composite
    function that mirrors R's ``suppressMessages(...)`` around an internal
    ``seqdist`` call, as ``seqdiff`` does -- can silence exactly the diagnostics R
    silences, without also silencing R's genuine ``warning()`` diagnostics.

    It subclasses ``UserWarning`` so existing ``pytest.warns(UserWarning)``
    expectations and default warning filters continue to see it.
    """


def _first_triangle_inequality_violation(
    matrix: np.ndarray, tol: float
) -> tuple[int, int, int] | None:
    """Port of TraMineR's ``checktriangleineq`` (``src/checktriangleineq.c``).

    Returns the first ``(i, j, z)`` for which going from ``i`` to ``j`` through an
    intermediate ``z`` costs at least ``tol`` less than the direct cost, or
    ``None`` when the matrix satisfies the inequality.  Indices are zero-based;
    the R original returns them one-based.

    The C loop is reproduced exactly, including its details: R stores matrices
    column-major, so ``matrix[j + i * n]`` is element ``[j, i]``; the diagonal of
    row ``i`` is forced to zero before that row is scanned; and ``z`` ranges over
    **every** index, including ``i`` and ``j`` themselves.
    """
    costs = np.array(matrix, dtype=float, copy=True)
    n = int(costs.shape[0])
    # ``Inf - Inf`` is a silent NaN in the R original's C loop, and a NaN
    # comparison there is simply false.  NumPy raises ``RuntimeWarning: invalid
    # value encountered in scalar subtract`` for the same arithmetic, which
    # reached the user ahead of the R-ordered diagnostic this function feeds
    # (cycle-10 N-18).  Suppressing it here restores the C semantics exactly:
    # the comparison below is false for NaN either way, so no branch changes.
    with np.errstate(invalid="ignore", over="ignore"):
        for i in range(n):
            costs[i, i] = 0.0
            for j in range(i + 1, n):
                direct = costs[j, i]
                for z in range(n):
                    if direct - (costs[z, i] + costs[z, j]) >= tol:
                        return i, j, z
    return None


def _checkcost_tolerance(costs: np.ndarray) -> float:
    """R ``checkcost``'s default tolerance: ``1e-07 * max(sma)``."""
    finite = costs[np.isfinite(costs)]
    return 1e-7 * float(np.max(finite)) if finite.size else 0.0


def _warn_cost_quality(
    costs,
    *,
    indel=None,
    labels=None,
    triangle_only: bool = False,
    symmetry: bool = False,
    stacklevel: int,
) -> None:
    """Emit R ``checkcost``'s non-fatal cost diagnostics (TraMineR 2.2-14).

    Faithful to ``SOURCES/TraMineR/R/TraMineR-checkcost.R`` for a ``(K, K)``
    matrix or a time-first ``(T, K, K)`` array, in R's order:

    1. a largest cost of 0 or below: reported alone, and nothing else is
       checked (lines 5-8);
    2. zero costs between distinct states, counted over R's ``as.dist`` lower
       triangle and reported once, with the first pair (lines 56-67);
    3. the first substitution triangle-inequality violation, reported once
       (lines 74-93);
    4. for a matrix only, **either** an indel-vector triangle-inequality
       violation, when ``indel`` has more than one element, **or** a cost
       greater than twice a scalar ``indel`` (lines 96-111).

    R signals all four with ``message()``, so they are
    :class:`CostDiagnosticMessage` warnings here, with R's text minus the
    newline ``message()`` appends; for an array each carries R's
    ``"At pos t,"`` prefix. R's fifth diagnostic, asymmetry, is a
    ``warning()`` the kernels emit (:func:`_warn_asymmetric_costs`); with
    ``symmetry=True`` it is emitted here instead, as a ``UserWarning`` at the
    end of each asymmetric matrix's checks, which is R's position (lines
    114-117; seqMD checks its channels this way, audit S10-11). R's
    zero-indel message (lines 39-41) cannot arise: traminex refuses such an
    indel (docs/API_PARITY.md).

    ``labels`` names the states as R's alphabet does (plus ``nr`` under
    ``with.missing``); without it the integer codes stand in. Numbers read
    as R's ``format()``. ``triangle_only`` runs only (3), after (1)'s
    short-circuit: traminex's documented choice for HAM, which R never checks.
    """
    array = np.asarray(costs, dtype=float)
    if array.ndim not in (2, 3) or array.shape[-1] != array.shape[-2] or array.size == 0:
        return
    time_varying = array.ndim == 3
    slices = array if time_varying else array[None, :, :]
    size = array.shape[-1]
    names = [str(label) for label in labels] if labels is not None else []
    if len(names) != size:
        names = [str(code) for code in range(size)]

    def emit(text: str) -> None:
        warnings.warn(text, CostDiagnosticMessage, stacklevel=stacklevel + 1)

    def number(value: float) -> str:
        return r_number_text(value, digits=7)

    if float(np.max(array)) <= 0:
        if not triangle_only:
            emit(" [!!] max substitution cost less or equal to 0!")
        return
    tol = _checkcost_tolerance(array)
    asymmetric = _asymmetric_slice_flags(slices) if symmetry else None
    report_zero = not triangle_only
    report_triangle = True
    for position, matrix in enumerate(slices, start=1):
        prefix = f"At pos {position}," if time_varying else ""
        if report_zero:
            zero_pairs = [
                (row, col)
                for col in range(size)
                for row in range(col + 1, size)
                if matrix[row, col] == 0
            ]
            if zero_pairs:
                verb = "is 0!" if len(zero_pairs) == 1 else "are 0!"
                emit(f" [!!] [!] {prefix} {len(zero_pairs)} substitution cost value(s) {verb}")
                row, col = zero_pairs[0]
                emit(
                    " [!!]      for example, cost is zero for the pair"
                    f" ( {names[col]} , {names[row]} )"
                )
                report_zero = False
        if report_triangle:
            violation = _first_triangle_inequality_violation(matrix, tol)
            if violation is not None:
                i, j, z = violation
                a, b, c = names[i], names[j], names[z]
                emit(
                    f" [!!] {prefix} at least, one substitution cost doesn't respect the"
                    " triangle inequality.\n"
                    f" [!!] replacing {a} with {c} (cost={number(matrix[i, z])}) and then"
                    f" {c} with {b} (cost={number(matrix[z, j])})\n"
                    f" [!] costs less than replacing directly {a} with {b}"
                    f" (cost={number(matrix[i, j])})\n"
                    f" [!] total difference ([{a}=>{c}] + [{c}=>{b}] - [{a}=>{b}]):"
                    f" {number(matrix[i, z] + matrix[z, j] - matrix[i, j])}"
                )
                report_triangle = False
        if not (time_varying or triangle_only or indel is None):
            indel_values = np.atleast_1d(np.asarray(indel, dtype=float))
            if indel_values.size > 1:
                # R builds cbind(sm, indel) / rbind(sm2, c(indel, 0)) and re-checks.
                if indel_values.size == size:
                    augmented = np.zeros((size + 1, size + 1), dtype=float)
                    augmented[:-1, :-1] = matrix
                    augmented[:-1, -1] = indel_values
                    augmented[-1, :-1] = indel_values
                    if _first_triangle_inequality_violation(augmented, tol) is not None:
                        emit(
                            " [!!] at least one indel cost does not respect the triangle"
                            " inequality."
                        )
            elif indel_values.size == 1 and np.any(matrix > 2.0 * indel_values[0]):
                emit(
                    " [!!] At least one substitution cost greater than twice the indel cost.\n"
                    "      Such substitution costs will never be used."
                )
        if asymmetric is not None and asymmetric[position - 1]:
            warnings.warn(
                # emit() adds a frame of its own; this call does not.
                "Non-symmetric substitution cost matrix!",
                UserWarning,
                stacklevel=stacklevel,
            )


def _subm_max_offdiag(subm: float | np.ndarray) -> float:
    """TraMineR OM-family maximum, scanning only the strict upper triangle.

    ``OMdistance::setParameters`` scans ``i < j``. That value enters raw distances
    (TWED's border value, OMloc's local costs) and the bounds R gets right (OMloc,
    OMstran), so it stays R's there; the normalization bounds of OM, OMslen, OMspell
    and TWED take :func:`_subm_max_bound` instead (TX-AUD-20261019-008), except that the
    bound of gmean of OMloc, OMspell and TWED keeps this value (owner decision of 2026-10-22,
    TX-REP-20261022-002).
    """
    if isinstance(subm, (int, float)):
        return float(subm)
    M = np.asarray(subm, dtype=float)
    if M.size == 0:
        return 0.0
    upper = M[np.triu_indices(M.shape[0], k=1)]
    return float(upper.max()) if upper.size else 0.0


def _subm_max_bound(subm: float | np.ndarray) -> float:
    """The largest off-diagonal substitution cost, for the normalization bound.

    R's bound reads ``maxscost`` from the strict upper triangle (:func:`_subm_max_offdiag`),
    so with a larger cost below the diagonal it is not the maximum possible distance that
    seqdist.Rd defines, and normalized OM, OMslen, OMspell and TWED distances exceed 1
    (TX-AUD-20261019-008, an R defect corrected). Every off-diagonal cell counts here; the
    raw-distance uses of ``maxscost`` (OMloc's local costs, TWED's border value) keep R's,
    and so do the bounds of OMspell and TWED except their gmean, which keeps R's value
    (:func:`_subm_max_offdiag`; TX-REP-20261022-002).
    """
    if np.ndim(subm) == 0:
        return float(subm)
    M = np.asarray(subm, dtype=float)
    if M.size == 0:
        return 0.0
    off = M[~np.eye(M.shape[0], dtype=bool)]
    return float(off.max()) if off.size else 0.0


def _subm_max_all(subm: float | np.ndarray) -> float:
    """Maximum over the complete matrix for R paths that call ``max(sm)``."""
    if isinstance(subm, (int, float)):
        return float(subm)
    M = np.asarray(subm, dtype=float)
    return float(M.max()) if M.size else 0.0


#: Methods whose caller-supplied costs R checks with ``checkcost``
#: (``seqdist.R:377-381``): the OM family and TWED.
_CHECKED_COST_METHODS = frozenset({"OM", "OMLOC", "OMSLEN", "OMSPELL", "OMSTRAN", "TWED"})


def _refuse_invalid_fixed_costs(method: str, subm: object, rows: list[np.ndarray]) -> None:
    """Refuse the substitution matrices R's ``checkcost`` stops on.

    Every row engine and dense branch calls this, so a negative cost, a non-zero
    diagonal, a matrix narrower than the observed codes or a non-finite cost between
    observed states is refused on every path, as it already was on the dense OM one
    (2026-09-22 audit S1b-01; TX-AUD-20261019-010). The
    exact-width check needs the declared alphabet and lives in ``seqdist``.
    """
    if subm is None:
        return
    if isinstance(subm, str):
        # Only seqdist resolves cost-method strings for these methods; the bare
        # float() failure named neither the option nor a remedy (audit S1a-14).
        raise ValueError(
            f"{method} subm={subm!r}: cost-method strings are resolved by seqdist; "
            "call seqdist(..., subm=...) or pass a substitution-cost matrix to pairwise"
        )
    costs = np.asarray(subm, dtype=float)
    if costs.ndim == 0:
        if costs < 0:
            raise ValueError(f"{method} substitution costs cannot be negative")
        return
    if costs.ndim != 2 or costs.shape[0] != costs.shape[1]:
        raise ValueError(f"{method} subm must be square")
    observed = max((int(row.max()) + 1 for row in rows if len(row)), default=0)
    if costs.shape[0] < observed:
        raise ValueError(
            f"{method} subm must cover every observed state code "
            f"(need at least {observed}x{observed})"
        )
    # A NaN or infinite cost between observed states is refused, as on the dense OM
    # route; one no row can reach stays accepted (docs/API_PARITY.md; TX-AUD-20261019-010).
    codes = np.concatenate([np.asarray(row, dtype=int).ravel() for row in rows]) if rows else []
    if not _observed_substitution_costs_are_finite(costs, np.asarray(codes, dtype=int)):
        raise ValueError(f"{method} substitution costs must be finite")
    if np.any(costs < 0):
        raise ValueError(f"{method} substitution costs cannot be negative")
    if costs.size and np.any(np.diag(costs) != 0):
        raise ValueError(f"{method} substitution-cost diagonal must be zero")


def _state_indel_spec(indel, *, n_states: int, method: str) -> tuple[float | np.ndarray, float]:
    """Validate a scalar/state-vector indel and return it plus its R normalization max."""
    if isinstance(indel, str):
        raise TypeError(f"{method} indel must be a scalar or state-cost vector")
    if np.isscalar(indel):
        value = _scalar_float(indel)
        if np.isnan(value) or value <= 0:
            raise ValueError(f"{method} indel cost must be strictly positive")
        return value, value
    values = np.asarray(indel, dtype=float)
    if values.ndim != 1 or values.shape[0] != n_states:
        raise ValueError(f"{method} indel vector length must equal the active state count")
    if np.any(np.isnan(values)) or np.any(values <= 0):
        raise ValueError(f"{method} indel cost must be strictly positive")
    return values, float(values.max()) if values.size else 0.0


def _twed_indel_value(indel, *, subm, nu: float, h: float, n_states: int) -> float:
    """Resolve TWED's distinct automatic indel and R's vector-to-max coercion."""
    if indel is None or (isinstance(indel, str) and indel.lower() == "auto"):
        return 2.0 * _subm_max_all(subm) + nu + h
    if isinstance(indel, str):
        raise TypeError("TWED indel must be a scalar or state-cost vector")
    if np.isscalar(indel):
        maximum = _scalar_float(indel)
    else:
        values = np.asarray(indel, dtype=float)
        if values.ndim != 1 or values.shape[0] != n_states:
            raise ValueError("TWED indel vector length must equal the active state count")
        maximum = float(values.max()) if values.size else 0.0
        # seqdist.R:301-305, with msg.warn's " [!!]" prefix (TX-REP-20261023-008).
        warnings.warn(
            " [!!] indel vector not supported by the chosen method, max(indel) used instead!",
            UserWarning,
            stacklevel=3,
        )
    if np.isnan(maximum) or maximum <= 0:
        raise ValueError("TWED indel cost must be strictly positive")
    return maximum


def ham_distance(a: np.ndarray, b: np.ndarray) -> float:
    if a.size == 0:
        return 0.0
    return float(np.mean(a != b))


def om_distance(
    a: np.ndarray,
    b: np.ndarray,
    indel: float = 1.0,
    subm: float | np.ndarray = 2.0,
    *,
    indel_by_pos: np.ndarray | None = None,
    subm_by_pos: np.ndarray | None = None,
    indel_by_state: np.ndarray | None = None,
) -> float:
    """Simple OM (edit distance) with optional position-dependent indel cost and uniform or matrix substitution cost.
    DP in O(L^2).

    ``indel_by_state`` (used by OMstran) gives a per-state indel cost: deleting/inserting a
    token of state ``s`` costs ``indel_by_state[s]``. It takes precedence over ``indel_by_pos``.
    """
    # TraMineR's OM kernel accepts a positive infinite scalar indel but makes
    # every non-identity comparison unreachable, including equal-length paths
    # that could theoretically use substitutions only.
    if indel_by_state is None and np.isposinf(indel):
        return 0.0 if np.array_equal(a, b) else float("inf")
    L1, L2 = len(a), len(b)
    use_subm_by_pos = subm_by_pos is not None
    if use_subm_by_pos:
        S = np.asarray(subm_by_pos, dtype=float)  # shape [L, K, K]

        # no scalar path when per-position provided
        def sc_pos(i, j, t):
            return 0.0 if i == j else float(S[t, i, j])

    else:
        if isinstance(subm, (int, float)):

            def sc(x, y):
                return 0.0 if x == y else float(subm)

        else:
            M = np.asarray(subm, dtype=float)

            def sc(x, y):
                return float(M[x, y])

    def indel_a(pos: int) -> float:
        if indel_by_state is not None:
            return float(indel_by_state[int(a[pos])])
        if indel_by_pos is not None:
            return float(indel_by_pos[pos])
        return indel

    def indel_b(pos: int) -> float:
        if indel_by_state is not None:
            return float(indel_by_state[int(b[pos])])
        if indel_by_pos is not None:
            return float(indel_by_pos[pos])
        return indel

    dp = np.zeros((L1 + 1, L2 + 1), dtype=float)
    # initialize borders with variable indel if provided
    for i in range(1, L1 + 1):
        dp[i, 0] = dp[i - 1, 0] + indel_a(i - 1)
    for j in range(1, L2 + 1):
        dp[0, j] = dp[0, j - 1] + indel_b(j - 1)
    for i in range(1, L1 + 1):
        for j in range(1, L2 + 1):
            c_i = indel_a(i - 1)
            c_j = indel_b(j - 1)
            if use_subm_by_pos:
                # choose time index for substitution: align by min(i-1, j-1)
                t = min(i - 1, j - 1)
                sub_cost = sc_pos(int(a[i - 1]), int(b[j - 1]), t)
            else:
                sub_cost = sc(a[i - 1], b[j - 1])
            dp[i, j] = min(
                dp[i - 1, j] + c_i,  # deletion from a
                dp[i, j - 1] + c_j,  # insertion into a
                dp[i - 1, j - 1] + sub_cost,
            )
    return float(dp[L1, L2])


def subseq_counts(
    a: np.ndarray,
    b: np.ndarray,
    ksize: int,
    *,
    prox: np.ndarray | None = None,
) -> np.ndarray:
    """Return TraMineR's hard- or soft-NMS common-subsequence count vector.

    Hard matching follows ``NMSdistance.cpp``.  A supplied ``prox`` follows the
    distinct multiplicative recurrence in ``NMSMSTSoftdistanceII.cpp``.
    """
    if prox is not None:
        return _soft_subseq_counts(a, b, ksize, prox)

    m, n = len(a), len(b)
    kvect = np.zeros(ksize, dtype=float)
    if m == 0 or n == 0 or ksize == 0:
        return kvect

    hmat = np.zeros((m, n), dtype=float)
    matches: list[tuple[int, int]] = []
    htot = 0.0
    for i in range(m):
        for j in range(n):
            w = 1.0 if a[i] == b[j] else 0.0
            if w != 0.0:
                hmat[i, j] = w
                matches.append((i, j))
                htot += w
    kvect[0] = htot
    if m > 1 and n > 1:
        k = 1
        while k < ksize and htot > 0:
            vmat = np.zeros((m, n), dtype=float)
            for j in range(n - 2, -1, -1):
                for i in range(m - 2, -1, -1):
                    vmat[i, j] = (
                        vmat[i + 1, j] + vmat[i, j + 1] - vmat[i + 1, j + 1] + hmat[i + 1, j + 1]
                    )
            if vmat[0, 0] == 0:
                kvect[k] = 0.0
                break
            htot = 0.0
            newh = np.zeros((m, n), dtype=float)
            for i, j in matches:
                newh[i, j] = vmat[i, j]
                htot += vmat[i, j]
            hmat = newh
            kvect[k] = htot
            k += 1
    return kvect


def _soft_subseq_counts(
    a: np.ndarray,
    b: np.ndarray,
    ksize: int,
    prox: np.ndarray,
) -> np.ndarray:
    """Port ``NMSMSTSoftdistanceII::computeattr`` for NMS proximity costs."""
    m, n = len(a), len(b)
    counts = np.zeros(ksize, dtype=float)
    if m == 0 or n == 0 or ksize == 0:
        return counts

    active_rows = m + 1
    active_columns = n + 1
    matches = np.zeros((active_rows, active_columns), dtype=float)
    recurrence = np.zeros_like(matches)
    total = 0.0
    for i in range(m):
        for j in range(n):
            weight = float(prox[int(a[i]), int(b[j])])
            matches[i, j] = weight
            recurrence[i, j] = weight
            total += weight

    counts[0] = total + 1.0
    if total == 0.0:
        return counts

    order = 0
    while active_rows != 0 and active_columns != 0:
        order += 1
        for i in range(active_rows):
            subtotal = 0.0
            for j in range(active_columns - 1, -1, -1):
                preceding = subtotal
                subtotal += recurrence[i, j]
                recurrence[i, j] = preceding

        total = 0.0
        for j in range(active_columns):
            subtotal = 0.0
            for i in range(active_rows - 1, -1, -1):
                preceding = subtotal
                subtotal += recurrence[i, j]
                recurrence[i, j] = preceding * matches[i, j]
                total += recurrence[i, j]

        if total == 0.0 or order >= ksize:
            return counts
        counts[order] = total
        active_rows -= 1
        active_columns -= 1
    return counts


def _validate_proximity_matrix(prox, *, n_states: int, method: str) -> np.ndarray:
    """Validate the TraMineR NMS/SVRspell proximity-matrix contract."""
    try:
        matrix = np.asarray(prox, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{method} prox must be of size {n_states} x {n_states}") from exc

    expected = (n_states, n_states)
    if matrix.ndim != 2 or matrix.shape != expected:
        raise ValueError(f"{method} prox must be of size {n_states} x {n_states}")
    if not np.all(np.isfinite(matrix)):
        raise ValueError(f"{method} prox must contain only finite values")

    tolerance = np.sqrt(np.finfo(float).eps)
    eigenvalues = np.linalg.eigvals(matrix)
    if np.any(eigenvalues.real < -tolerance):
        values = " ".join(f"{value.real:.3g}" for value in eigenvalues)
        warnings.warn(
            f"{method} 'prox' is not positive semi-definite. Eigen values: {values}",
            UserWarning,
            stacklevel=4,
        )
    return matrix


def _rescale_by_power_of_two(*arrays: np.ndarray) -> int:
    """Divide every array in place by the power of two that brings their largest
    magnitude into [0.5, 1); return its exponent. Exact: only the exponent changes."""
    peak = max(float(np.max(np.abs(array))) for array in arrays)
    if peak == 0.0 or not math.isfinite(peak):
        return 0
    scale = math.frexp(peak)[1]
    for array in arrays:
        np.ldexp(array, -scale, out=array)
    return scale


def nmsmst_counts(
    a: np.ndarray, b: np.ndarray, dur_a: np.ndarray, dur_b: np.ndarray, ksize: int
) -> np.ndarray:
    """NMSMST minimum-shared-time subsequence counts (TraMineR NMSMSTdistance.cpp).

    ``kvect[k]`` is the total minimum-shared-time over all common length-(k+1) subsequences,
    where a matched spell pair contributes ``min(dur_a[i], dur_b[j])``.

    With Numba dispatched, :func:`_nmsmst_counts_nb` runs the same loops with the same
    floating-point operations in the same order, so its counts equal these bit for bit
    (2026-10-25 study, finding F-06; ``tests/test_nmsmst_svrspell_parity.py``).
    """
    if _HAVE_NUMBA:
        return _nmsmst_counts_nb(*_spell_kernel_arrays(a, b, dur_a, dur_b), int(ksize))
    return _nmsmst_counts(a, b, dur_a, dur_b, ksize, None)


def _spell_kernel_arrays(
    a: np.ndarray, b: np.ndarray, dur_a: np.ndarray, dur_b: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """The spell rows and durations as the compiled spell kernels take them, one type each."""
    return (
        np.ascontiguousarray(a, dtype=np.int64),
        np.ascontiguousarray(b, dtype=np.int64),
        np.ascontiguousarray(dur_a, dtype=np.float64),
        np.ascontiguousarray(dur_b, dtype=np.float64),
    )


def _nmsmst_counts(
    a: np.ndarray,
    b: np.ndarray,
    dur_a: np.ndarray,
    dur_b: np.ndarray,
    ksize: int,
    exponents_out: list[int] | None,
) -> np.ndarray:
    """:func:`nmsmst_counts`; with ``exponents_out`` the state is divided by a power of two
    before each length and ``kvect[k] * 2**exponents_out[k]`` is the count
    (TX-AUD-20261019-029)."""
    m, n = len(a), len(b)
    kvect = np.zeros(ksize, dtype=float)
    if m == 0 or n == 0 or ksize == 0:
        return kvect
    R, C = m + 1, n + 1
    e1 = np.zeros((R, C), dtype=float)
    e = np.zeros((R, C), dtype=float)
    t1 = np.zeros((R, C), dtype=float)
    t = np.zeros((R, C), dtype=float)
    st = 0.0
    for i in range(m):
        for j in range(n):
            if a[i] == b[j]:
                e1[i, j] = 1.0
                e[i, j] = 1.0
                mst = min(float(dur_a[i]), float(dur_b[j]))
                t1[i, j] = mst
                t[i, j] = mst
                st += mst
    kvect[0] = st
    shift = 0
    if exponents_out is not None:
        exponents_out[:] = [0] * ksize
    if st == 0.0:
        return kvect
    mrows, ncols, k = R, C, 0
    while mrows != 0 and ncols != 0 and k + 1 < ksize:
        k += 1
        if exponents_out is not None:
            shift += _rescale_by_power_of_two(e, t)
        for irow in range(mrows):
            s = st_ = 0.0
            for jcol in range(ncols - 1, -1, -1):
                temp, tempt = s, st_
                s += e[irow, jcol]
                st_ += t[irow, jcol]
                e[irow, jcol] = temp
                t[irow, jcol] = tempt
        temps = temptt = 0.0
        for jcol in range(ncols):
            s = st_ = 0.0
            for irow in range(mrows - 1, -1, -1):
                temp, tempt = s, st_
                s += e[irow, jcol]
                st_ += t[irow, jcol]
                e[irow, jcol] = temp * e1[irow, jcol]
                t[irow, jcol] = e1[irow, jcol] * (e[irow, jcol] * t1[irow, jcol] + tempt)
                temps += e[irow, jcol]
                temptt += t[irow, jcol]
        if temps == 0.0:
            break
        kvect[k] = temptt
        if exponents_out is not None:
            exponents_out[k] = shift
        mrows -= 1
        ncols -= 1
    return kvect


def svrspell_counts(
    a: np.ndarray,
    b: np.ndarray,
    dur_a: np.ndarray,
    dur_b: np.ndarray,
    ksize: int,
    softmatch: np.ndarray,
) -> np.ndarray:
    """SVRspell soft duration-weighted subsequence counts (TraMineR NMSDURSoftdistance.cpp).

    Tracks existence (e) and duration moments (t_i, t_j, t_ij) with soft matching; ``kvect[k]``
    is the total ``t_ij`` over common length-(k+1) soft subsequences.

    With Numba dispatched, :func:`_svrspell_counts_nb` computes them bit for bit as
    :func:`nmsmst_counts` describes.
    """
    if _HAVE_NUMBA:
        soft = np.ascontiguousarray(softmatch, dtype=np.float64)
        return _svrspell_counts_nb(*_spell_kernel_arrays(a, b, dur_a, dur_b), int(ksize), soft)
    return _svrspell_counts(a, b, dur_a, dur_b, ksize, softmatch, None)


def _svrspell_counts(
    a: np.ndarray,
    b: np.ndarray,
    dur_a: np.ndarray,
    dur_b: np.ndarray,
    ksize: int,
    softmatch: np.ndarray,
    exponents_out: list[int] | None,
) -> np.ndarray:
    """:func:`svrspell_counts`; ``exponents_out`` works as in :func:`_nmsmst_counts`."""
    m, n = len(a), len(b)
    kvect = np.zeros(ksize, dtype=float)
    if m == 0 or n == 0 or ksize == 0:
        return kvect
    R, C = m + 1, n + 1
    e1 = np.zeros((R, C))
    e = np.zeros((R, C))
    t_i = np.zeros((R, C))
    t_j = np.zeros((R, C))
    t_ij = np.zeros((R, C))
    tot_t_ij = 0.0
    for i in range(m):
        for j in range(n):
            sf = float(softmatch[int(a[i]), int(b[j])])
            ti = float(dur_a[i])
            tj = float(dur_b[j])
            e1[i, j] = sf
            e[i, j] = sf
            t_i[i, j] = sf * ti
            t_j[i, j] = sf * tj
            t_ij[i, j] = sf * ti * (sf * tj)
            tot_t_ij += t_ij[i, j]
    kvect[0] = tot_t_ij
    shift = 0
    if exponents_out is not None:
        exponents_out[:] = [0] * ksize
    if tot_t_ij == 0.0:
        return kvect
    mrows, ncols, k = R, C, 0
    while mrows != 0 and ncols != 0 and k + 1 < ksize:
        k += 1
        if exponents_out is not None:
            shift += _rescale_by_power_of_two(e, t_i, t_j, t_ij)
        for irow in range(mrows):
            se = si = sj = sij = 0.0
            for jcol in range(ncols - 1, -1, -1):
                te, ti_, tj_, tij_ = se, si, sj, sij
                se += e[irow, jcol]
                si += t_i[irow, jcol]
                sj += t_j[irow, jcol]
                sij += t_ij[irow, jcol]
                e[irow, jcol] = te
                t_i[irow, jcol] = ti_
                t_j[irow, jcol] = tj_
                t_ij[irow, jcol] = tij_
        tot_e = 0.0
        for jcol in range(ncols):
            se = si = sj = sij = 0.0
            for irow in range(mrows - 1, -1, -1):
                te, ti_, tj_, tij_ = se, si, sj, sij
                se += e[irow, jcol]
                si += t_i[irow, jcol]
                sj += t_j[irow, jcol]
                sij += t_ij[irow, jcol]
                e[irow, jcol] = te
                t_i[irow, jcol] = ti_
                t_j[irow, jcol] = tj_
                t_ij[irow, jcol] = tij_
                tot_e += e[irow, jcol]
        if tot_e == 0.0:
            break
        tot_t_ij = 0.0
        for irow in range(mrows):
            ti = float(dur_a[irow]) if irow < m else 0.0
            for jcol in range(ncols):
                tj = float(dur_b[jcol]) if jcol < n else 0.0
                sf = e1[irow, jcol]
                e[irow, jcol] = sf * e[irow, jcol]
                t_ij[irow, jcol] = sf * (
                    t_ij[irow, jcol]
                    + ti * t_j[irow, jcol]
                    + tj * t_i[irow, jcol]
                    + e[irow, jcol] * ti * tj
                )
                t_i[irow, jcol] = sf * (t_i[irow, jcol] + ti * e[irow, jcol])
                t_j[irow, jcol] = sf * (t_j[irow, jcol] + tj * e[irow, jcol])
                tot_t_ij += t_ij[irow, jcol]
        kvect[k] = tot_t_ij
        if exponents_out is not None:
            exponents_out[k] = shift
        mrows -= 1
        ncols -= 1
    return kvect


def lcs_distance(a: np.ndarray, b: np.ndarray) -> float:
    """LCS-based distance: d = |a| + |b| - 2 * LCS(a,b)."""
    L1, L2 = len(a), len(b)
    dp = np.zeros((L1 + 1, L2 + 1), dtype=int)
    for i in range(L1):
        for j in range(L2):
            if a[i] == b[j]:
                dp[i + 1, j + 1] = dp[i, j] + 1
            else:
                dp[i + 1, j + 1] = max(dp[i, j + 1], dp[i + 1, j])
    lcs = dp[L1, L2]
    return float(L1 + L2 - 2 * lcs)


def lcp_distance(a: np.ndarray, b: np.ndarray) -> float:
    """LCP-based distance: d = |a| + |b| - 2 * LCP(a,b)."""
    L1 = min(len(a), len(b))
    k = 0
    while k < L1 and a[k] == b[k]:
        k += 1
    return float(len(a) + len(b) - 2 * k)


def rlcp_distance(a: np.ndarray, b: np.ndarray) -> float:
    """RLCP-based distance: common suffix length; d = |a|+|b| - 2*RLCP."""
    L1 = len(a)
    L2 = len(b)
    k = 0
    m = L1 if L1 < L2 else L2
    while k < m and a[L1 - 1 - k] == b[L2 - 1 - k]:
        k += 1
    return float(L1 + L2 - 2 * k)


def twed_distance(
    a: np.ndarray,
    b: np.ndarray,
    *,
    nu: float,
    lam: float,
    subm: float | np.ndarray = 2.0,
    indel: float | None = None,
) -> float:
    """Time-Warp Edit Distance (Marteau 2009), faithful to TraMineR TWEDdistance.cpp.

    Uses the substitution matrix ``subm`` as the local cost, integer timestamps (so the
    timestamp gap between adjacent positions is 1), stiffness ``nu`` and penalty ``lam``
    (``h`` in seqdist). State 0 is the dummy predecessor before the first position.
    """
    L1, L2 = len(a), len(b)
    if isinstance(subm, (int, float)):
        subm_f = float(subm)

        def sc(x: int, y: int) -> float:
            return 0.0 if x == y else subm_f

    else:
        M = np.asarray(subm, dtype=float)

        def sc(x: int, y: int) -> float:
            return float(M[x, y])

    if indel is None:
        # seqdist.R's automatic TWED indel uses max over the complete matrix,
        # while the border value below uses R's i < j maximum.
        indel = 2.0 * _subm_max_all(subm) + nu + lam
    if L1 == 0 or L2 == 0:
        # TWEDdistance.cpp explicitly overwrites the empty border with n*indel.
        return float((L1 + L2) * indel)
    infinite_indel = np.isposinf(indel)
    # TraMineR fills the impossible border cells with a finite ``1000 * maxscost``
    # sentinel. traminex keeps the TWED definition's infinite border: the two
    # agree whenever the sentinel exceeds every real warping cost, but with
    # costs small relative to nu/h R's sentinel binds and caps distances -- a
    # recorded divergence (docs/API_PARITY.md; tests/test_seqdist_numeric_edges.py).
    # An infinite indel initializes fmat[0,0] as 0*Inf=NaN and exposes the
    # exact finite sentinel value, which is reproduced here.
    big = 1000.0 * _subm_max_offdiag(subm) if infinite_indel else float("inf")
    fmat = np.zeros((L1 + 1, L2 + 1), dtype=float)
    if infinite_indel:
        fmat[0, 0] = np.nan
    for i in range(1, L1 + 1):
        i_state = int(a[i - 1])
        i_m1 = 0 if i == 1 else int(a[i - 2])
        for j in range(1, L2 + 1):
            j_state = int(b[j - 1])
            j_m1 = 0 if j == 1 else int(b[j - 2])
            if i_state == j_state and i_m1 == j_m1:
                cost = 0.0
            else:
                cost = sc(i_state, j_state) + sc(i_m1, j_m1)
            sub = fmat[i - 1, j - 1] + cost + 2.0 * nu * abs(i - j)
            if j - 1 == 0:
                i_warp = big
                if i - 1 > 0:
                    sub = big
            else:
                i_warp = fmat[i, j - 1] + sc(j_state, j_m1) + nu + lam
            if i - 1 == 0:
                j_warp = big
                if j - 1 > 0:
                    sub = big
            else:
                j_warp = fmat[i - 1, j] + sc(i_state, i_m1) + nu + lam
            fmat[i, j] = min(i_warp, j_warp, sub)
    return float(fmat[L1, L2])


def omslen_distance(
    a: np.ndarray,
    b: np.ndarray,
    *,
    indel: float | np.ndarray,
    subm: float | np.ndarray,
    link: str,
    dur_a: np.ndarray,
    dur_b: np.ndarray,
    h: float = 0.5,
) -> float:
    """OMslen distance (TraMineR OMvdistance.cpp): OM weighted by spell lengths.

    Each position carries its spell length; the cost weight is dur^(-h) (R prep:
    dur.mat ^ (-1*h), h defaults to 0.5). getIndel(state,t) = indel[state] * w[t];
    getSubCost uses link="mean" -> sm/2*(w_i+w_j) or
    link="gmean" -> sm*sqrt(w_i*w_j).
    """
    L1, L2 = len(a), len(b)
    w_a = spell_length_power(dur_a, -float(h))
    w_b = spell_length_power(dur_b, -float(h))
    if isinstance(subm, (int, float)):
        sval = float(subm)

        def sc(x: int, y: int) -> float:
            return 0.0 if x == y else sval

    else:
        M = np.asarray(subm, dtype=float)

        def sc(x: int, y: int) -> float:
            return float(M[x, y])

    if np.isscalar(indel):
        scalar_indel = _scalar_float(indel)

        def indel_cost(state: int) -> float:
            return scalar_indel

    else:
        indel_values = np.asarray(indel, dtype=float)

        def indel_cost(state: int) -> float:
            return float(indel_values[state])

    if link == "mean":

        def sub_cost(x: int, y: int, wi: float, wj: float) -> float:
            # seqdist.R:737 halves the cost first, so a cost at the largest double
            # does not overflow to Inf (audit S1a-05); halving is exact.
            return (sc(x, y) / 2.0) * (wi + wj)

    elif link == "gmean":

        def sub_cost(x: int, y: int, wi: float, wj: float) -> float:
            return sc(x, y) * float(np.sqrt(wi * wj))

    else:
        raise ValueError("OMslen link must be 'mean' or 'gmean'")

    dp = np.zeros((L1 + 1, L2 + 1), dtype=float)
    for i in range(1, L1 + 1):
        dp[i, 0] = dp[i - 1, 0] + indel_cost(int(a[i - 1])) * w_a[i - 1]
    for j in range(1, L2 + 1):
        dp[0, j] = dp[0, j - 1] + indel_cost(int(b[j - 1])) * w_b[j - 1]
    for i in range(1, L1 + 1):
        for j in range(1, L2 + 1):
            dp[i, j] = min(
                dp[i - 1, j] + indel_cost(int(a[i - 1])) * w_a[i - 1],
                dp[i, j - 1] + indel_cost(int(b[j - 1])) * w_b[j - 1],
                dp[i - 1, j - 1] + sub_cost(int(a[i - 1]), int(b[j - 1]), w_a[i - 1], w_b[j - 1]),
            )
    return float(dp[L1, L2])


def omspell_distance(
    a: np.ndarray,
    b: np.ndarray,
    *,
    indel: float | np.ndarray,
    subm: float | np.ndarray,
    timecost: float,
    tpow: float = 1.0,
) -> float:
    """OMspell distance (TraMineR OMPerdistance.cpp): OM on the DSS with duration costs.

    On the spell (DSS) sequences with durations d = dur^tpow - 1 (seqdist.R subtracts
    c=1 for OMspell): getIndel(state,d) = indel[state] + timecost*d;
    getSubCost = timecost*|d_i-d_j| if same state, else sm[i,j] + (d_i+d_j)*timecost.
    """
    sa, da = spell_encode(a)
    sb, db = spell_encode(b)
    da = spell_length_power(da, float(tpow)) - 1.0
    db = spell_length_power(db, float(tpow)) - 1.0
    m, n = len(sa), len(sb)
    if isinstance(subm, (int, float)):
        sval = float(subm)

        def sc(x: int, y: int) -> float:
            return 0.0 if x == y else sval

    else:
        M = np.asarray(subm, dtype=float)

        def sc(x: int, y: int) -> float:
            return float(M[x, y])

    if np.isscalar(indel):
        scalar_indel = _scalar_float(indel)

        def state_indel(state: int) -> float:
            return scalar_indel

    else:
        indel_values = np.asarray(indel, dtype=float)

        def state_indel(state: int) -> float:
            return float(indel_values[state])

    def get_indel(state: int, dur: float) -> float:
        return state_indel(state) + timecost * dur

    def sub_cost(xi: int, xj: int, di: float, dj: float) -> float:
        if xi == xj:
            return timecost * abs(di - dj)
        return sc(xi, xj) + (di + dj) * timecost

    dp = np.zeros((m + 1, n + 1), dtype=float)
    for i in range(1, m + 1):
        dp[i, 0] = dp[i - 1, 0] + get_indel(int(sa[i - 1]), da[i - 1])
    for j in range(1, n + 1):
        dp[0, j] = dp[0, j - 1] + get_indel(int(sb[j - 1]), db[j - 1])
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            dp[i, j] = min(
                dp[i - 1, j] + get_indel(int(sa[i - 1]), da[i - 1]),
                dp[i, j - 1] + get_indel(int(sb[j - 1]), db[j - 1]),
                dp[i - 1, j - 1] + sub_cost(int(sa[i - 1]), int(sb[j - 1]), da[i - 1], db[j - 1]),
            )
    return float(dp[m, n])


# Optional numba-accelerated helpers
# Defined whenever numba is importable, independently of whether dispatch is
# enabled: benchmarks/run_benchmarks.py toggles ``_HAVE_NUMBA`` at runtime to
# compare the two backends, so the kernels must exist even under
# ``TRAMINEX_NUMBA=off``. ``@njit`` compiles lazily, so defining them is cheap.
if _NUMBA_IMPORTED:

    @njit
    def _c_weighted_pair_sum_nb(
        D: np.ndarray, members: np.ndarray, w: np.ndarray
    ) -> tuple[float, float]:
        # traminex.diss._c_weighted_pair_sum, as weightedinertia.cpp:17-31 runs it:
        # i outer, each later j inner, result += w[j] * w[i] * D[j, i], one double.
        size = members.shape[0]
        result = 0.0
        totweights = 0.0
        for a in range(size):
            i = members[a]
            wi = w[i]
            for b in range(a + 1, size):
                j = members[b]
                result += w[j] * wi * D[j, i]
            totweights += wi
        return result, totweights

    @njit
    def _nmsmst_counts_nb(
        a: np.ndarray, b: np.ndarray, dur_a: np.ndarray, dur_b: np.ndarray, ksize: int
    ) -> np.ndarray:
        # _nmsmst_counts without rescaling, line for line: the same operations in the same order.
        m, n = a.shape[0], b.shape[0]
        kvect = np.zeros(ksize, dtype=np.float64)
        if m == 0 or n == 0 or ksize == 0:
            return kvect
        R, C = m + 1, n + 1
        e1 = np.zeros((R, C), dtype=np.float64)
        e = np.zeros((R, C), dtype=np.float64)
        t1 = np.zeros((R, C), dtype=np.float64)
        t = np.zeros((R, C), dtype=np.float64)
        st = 0.0
        for i in range(m):
            for j in range(n):
                if a[i] == b[j]:
                    e1[i, j] = 1.0
                    e[i, j] = 1.0
                    # Python's min(x, y): x unless y < x.
                    mst = dur_b[j] if dur_b[j] < dur_a[i] else dur_a[i]
                    t1[i, j] = mst
                    t[i, j] = mst
                    st += mst
        kvect[0] = st
        if st == 0.0:
            return kvect
        mrows, ncols, k = R, C, 0
        while mrows != 0 and ncols != 0 and k + 1 < ksize:
            k += 1
            for irow in range(mrows):
                s = 0.0
                st_ = 0.0
                for jcol in range(ncols - 1, -1, -1):
                    temp, tempt = s, st_
                    s += e[irow, jcol]
                    st_ += t[irow, jcol]
                    e[irow, jcol] = temp
                    t[irow, jcol] = tempt
            temps = 0.0
            temptt = 0.0
            for jcol in range(ncols):
                s = 0.0
                st_ = 0.0
                for irow in range(mrows - 1, -1, -1):
                    temp, tempt = s, st_
                    s += e[irow, jcol]
                    st_ += t[irow, jcol]
                    e[irow, jcol] = temp * e1[irow, jcol]
                    t[irow, jcol] = e1[irow, jcol] * (e[irow, jcol] * t1[irow, jcol] + tempt)
                    temps += e[irow, jcol]
                    temptt += t[irow, jcol]
            if temps == 0.0:
                break
            kvect[k] = temptt
            mrows -= 1
            ncols -= 1
        return kvect

    @njit
    def _svrspell_counts_nb(
        a: np.ndarray,
        b: np.ndarray,
        dur_a: np.ndarray,
        dur_b: np.ndarray,
        ksize: int,
        softmatch: np.ndarray,
    ) -> np.ndarray:
        # _svrspell_counts without rescaling, line for line: the same operations in the same order.
        m, n = a.shape[0], b.shape[0]
        kvect = np.zeros(ksize, dtype=np.float64)
        if m == 0 or n == 0 or ksize == 0:
            return kvect
        R, C = m + 1, n + 1
        e1 = np.zeros((R, C), dtype=np.float64)
        e = np.zeros((R, C), dtype=np.float64)
        t_i = np.zeros((R, C), dtype=np.float64)
        t_j = np.zeros((R, C), dtype=np.float64)
        t_ij = np.zeros((R, C), dtype=np.float64)
        tot_t_ij = 0.0
        for i in range(m):
            for j in range(n):
                sf = softmatch[a[i], b[j]]
                ti = dur_a[i]
                tj = dur_b[j]
                e1[i, j] = sf
                e[i, j] = sf
                t_i[i, j] = sf * ti
                t_j[i, j] = sf * tj
                t_ij[i, j] = sf * ti * (sf * tj)
                tot_t_ij += t_ij[i, j]
        kvect[0] = tot_t_ij
        if tot_t_ij == 0.0:
            return kvect
        mrows, ncols, k = R, C, 0
        while mrows != 0 and ncols != 0 and k + 1 < ksize:
            k += 1
            for irow in range(mrows):
                se = 0.0
                si = 0.0
                sj = 0.0
                sij = 0.0
                for jcol in range(ncols - 1, -1, -1):
                    te, ti_, tj_, tij_ = se, si, sj, sij
                    se += e[irow, jcol]
                    si += t_i[irow, jcol]
                    sj += t_j[irow, jcol]
                    sij += t_ij[irow, jcol]
                    e[irow, jcol] = te
                    t_i[irow, jcol] = ti_
                    t_j[irow, jcol] = tj_
                    t_ij[irow, jcol] = tij_
            tot_e = 0.0
            for jcol in range(ncols):
                se = 0.0
                si = 0.0
                sj = 0.0
                sij = 0.0
                for irow in range(mrows - 1, -1, -1):
                    te, ti_, tj_, tij_ = se, si, sj, sij
                    se += e[irow, jcol]
                    si += t_i[irow, jcol]
                    sj += t_j[irow, jcol]
                    sij += t_ij[irow, jcol]
                    e[irow, jcol] = te
                    t_i[irow, jcol] = ti_
                    t_j[irow, jcol] = tj_
                    t_ij[irow, jcol] = tij_
                    tot_e += e[irow, jcol]
            if tot_e == 0.0:
                break
            tot_t_ij = 0.0
            for irow in range(mrows):
                ti = dur_a[irow] if irow < m else 0.0
                for jcol in range(ncols):
                    tj = dur_b[jcol] if jcol < n else 0.0
                    sf = e1[irow, jcol]
                    e[irow, jcol] = sf * e[irow, jcol]
                    t_ij[irow, jcol] = sf * (
                        t_ij[irow, jcol]
                        + ti * t_j[irow, jcol]
                        + tj * t_i[irow, jcol]
                        + e[irow, jcol] * ti * tj
                    )
                    t_i[irow, jcol] = sf * (t_i[irow, jcol] + ti * e[irow, jcol])
                    t_j[irow, jcol] = sf * (t_j[irow, jcol] + tj * e[irow, jcol])
                    tot_t_ij += t_ij[irow, jcol]
            kvect[k] = tot_t_ij
            mrows -= 1
            ncols -= 1
        return kvect

    @njit
    def _rlcp_len_nb(a: np.ndarray, b: np.ndarray) -> int:
        L1 = a.shape[0]
        L2 = b.shape[0]
        m = L1 if L1 < L2 else L2
        k = 0
        for i in range(m):
            if a[L1 - 1 - i] == b[L2 - 1 - i]:
                k += 1
            else:
                break
        return k

    @njit(parallel=True)
    def _pairwise_rlcp_nb(seqmat: np.ndarray) -> np.ndarray:
        n = seqmat.shape[0]
        L = seqmat.shape[1]
        D = np.zeros((n, n), dtype=np.float64)
        for i in prange(n):
            for j in range(i + 1, n):
                k = _rlcp_len_nb(seqmat[i], seqmat[j])
                D[i, j] = D[j, i] = float(2 * L - 2 * k)
        return D

    @njit
    def _dhd_distance_nb(a: np.ndarray, b: np.ndarray, freq: np.ndarray, variant: int) -> float:
        L = a.shape[0]
        acc = 0.0
        for t in range(L):
            if a[t] == b[t]:
                continue
            p = freq[t]
            if variant == 0:  # mean
                acc += 1.0 - 0.5 * (p[a[t]] + p[b[t]])
            else:  # sum
                acc += 2.0 - (p[a[t]] + p[b[t]])
        return acc / L if L > 0 else 0.0

    @njit(parallel=True)
    def _pairwise_dhd_nb(seqmat: np.ndarray, freq: np.ndarray, variant: int) -> np.ndarray:
        n = seqmat.shape[0]
        D = np.zeros((n, n), dtype=np.float64)
        for i in prange(n):
            for j in range(i + 1, n):
                v = _dhd_distance_nb(seqmat[i], seqmat[j], freq, variant)
                D[i, j] = D[j, i] = v
        return D

    @njit
    def _om_distance_scalar_nb(a: np.ndarray, b: np.ndarray, indel: float, subm: float) -> float:
        """OM distance with scalar substitution cost - Numba accelerated."""
        L1 = len(a)
        L2 = len(b)
        dp = np.zeros((L1 + 1, L2 + 1), dtype=np.float64)
        for i in range(1, L1 + 1):
            dp[i, 0] = dp[i - 1, 0] + indel
        for j in range(1, L2 + 1):
            dp[0, j] = dp[0, j - 1] + indel
        for i in range(1, L1 + 1):
            for j in range(1, L2 + 1):
                sub_cost = 0.0 if a[i - 1] == b[j - 1] else subm
                dp[i, j] = min(
                    dp[i - 1, j] + indel, dp[i, j - 1] + indel, dp[i - 1, j - 1] + sub_cost
                )
        return dp[L1, L2]

    @njit(parallel=True)
    def _pairwise_om_scalar_nb(seqmat: np.ndarray, indel: float, subm: float) -> np.ndarray:
        """Pairwise OM with scalar costs - Numba accelerated."""
        n = seqmat.shape[0]
        D = np.zeros((n, n), dtype=np.float64)
        for i in prange(n):
            for j in range(i + 1, n):
                d = _om_distance_scalar_nb(seqmat[i], seqmat[j], indel, subm)
                D[i, j] = D[j, i] = d
        return D

    @njit
    def _om_distance_matrix_nb(
        a: np.ndarray, b: np.ndarray, indel: float, subm: np.ndarray
    ) -> float:
        """OM distance with matrix substitution costs - Numba accelerated."""
        L1 = len(a)
        L2 = len(b)
        dp = np.zeros((L1 + 1, L2 + 1), dtype=np.float64)
        for i in range(1, L1 + 1):
            dp[i, 0] = dp[i - 1, 0] + indel
        for j in range(1, L2 + 1):
            dp[0, j] = dp[0, j - 1] + indel
        for i in range(1, L1 + 1):
            for j in range(1, L2 + 1):
                sub_cost = subm[a[i - 1], b[j - 1]]
                dp[i, j] = min(
                    dp[i - 1, j] + indel, dp[i, j - 1] + indel, dp[i - 1, j - 1] + sub_cost
                )
        return dp[L1, L2]

    @njit(parallel=True)
    def _pairwise_om_matrix_nb(seqmat: np.ndarray, indel: float, subm: np.ndarray) -> np.ndarray:
        """Pairwise OM with matrix costs - Numba accelerated."""
        n = seqmat.shape[0]
        D = np.zeros((n, n), dtype=np.float64)
        for i in prange(n):
            for j in range(i + 1, n):
                d = _om_distance_matrix_nb(seqmat[i], seqmat[j], indel, subm)
                D[i, j] = D[j, i] = d
        return D

    @njit
    def _lcs_distance_nb(a: np.ndarray, b: np.ndarray) -> float:
        """LCS distance - Numba accelerated."""
        L1 = len(a)
        L2 = len(b)
        dp = np.zeros((L1 + 1, L2 + 1), dtype=np.int32)
        for i in range(L1):
            for j in range(L2):
                if a[i] == b[j]:
                    dp[i + 1, j + 1] = dp[i, j] + 1
                else:
                    dp[i + 1, j + 1] = max(dp[i, j + 1], dp[i + 1, j])
        lcs = dp[L1, L2]
        return float(L1 + L2 - 2 * lcs)

    @njit(parallel=True)
    def _pairwise_lcs_nb(seqmat: np.ndarray) -> np.ndarray:
        """Pairwise LCS distance - Numba accelerated."""
        n = seqmat.shape[0]
        D = np.zeros((n, n), dtype=np.float64)
        for i in prange(n):
            for j in range(i + 1, n):
                d = _lcs_distance_nb(seqmat[i], seqmat[j])
                D[i, j] = D[j, i] = d
        return D

    @njit
    def _lcp_distance_nb(a: np.ndarray, b: np.ndarray) -> float:
        """LCP distance - Numba accelerated."""
        L1 = min(len(a), len(b))
        k = 0
        while k < L1 and a[k] == b[k]:
            k += 1
        return float(len(a) + len(b) - 2 * k)

    @njit(parallel=True)
    def _pairwise_lcp_nb(seqmat: np.ndarray) -> np.ndarray:
        """Pairwise LCP distance - Numba accelerated."""
        n = seqmat.shape[0]
        D = np.zeros((n, n), dtype=np.float64)
        for i in prange(n):
            for j in range(i + 1, n):
                d = _lcp_distance_nb(seqmat[i], seqmat[j])
                D[i, j] = D[j, i] = d
        return D

    @njit
    def _chi2_distance_nb(a: np.ndarray, b: np.ndarray, K: int) -> float:
        """Chi-square distance on duration distributions - Numba accelerated."""
        ua = np.zeros(K, dtype=np.float64)
        ub = np.zeros(K, dtype=np.float64)
        for i in range(len(a)):
            ua[a[i]] += 1.0
        for i in range(len(b)):
            ub[b[i]] += 1.0
        suma = np.sum(ua)
        sumb = np.sum(ub)
        if suma > 0:
            ua = ua / suma
        if sumb > 0:
            ub = ub / sumb
        chi2 = 0.0
        for k in range(K):
            denom = ua[k] + ub[k] + 1e-12
            chi2 += (ua[k] - ub[k]) ** 2 / denom
        return chi2

    @njit(parallel=True)
    def _pairwise_chi2_nb(seqmat: np.ndarray, K: int) -> np.ndarray:
        """Pairwise chi-square distance - Numba accelerated."""
        n = seqmat.shape[0]
        D = np.zeros((n, n), dtype=np.float64)
        for i in prange(n):
            for j in range(i + 1, n):
                d = _chi2_distance_nb(seqmat[i], seqmat[j], K)
                D[i, j] = D[j, i] = d
        return D

    @njit
    def _jsd_transition_distance_nb(a: np.ndarray, b: np.ndarray, K: int, eps: float) -> float:
        """JSD transition distance - Numba accelerated."""
        if len(a) < 2 or len(b) < 2:
            return 0.0
        pa = np.zeros(K * K, dtype=np.float64)
        pb = np.zeros(K * K, dtype=np.float64)
        for i in range(len(a) - 1):
            idx = a[i] * K + a[i + 1]
            pa[idx] += 1.0
        for i in range(len(b) - 1):
            idx = b[i] * K + b[i + 1]
            pb[idx] += 1.0
        sa = np.sum(pa)
        sb = np.sum(pb)
        if sa > 0:
            pa = pa / sa
        if sb > 0:
            pb = pb / sb
        m = 0.5 * (pa + pb)
        kl_a = 0.0
        kl_b = 0.0
        for k in range(K * K):
            if pa[k] > 0:
                kl_a += pa[k] * (np.log(pa[k] + eps) - np.log(m[k] + eps))
            if pb[k] > 0:
                kl_b += pb[k] * (np.log(pb[k] + eps) - np.log(m[k] + eps))
        return 0.5 * kl_a + 0.5 * kl_b

    @njit(parallel=True)
    def _pairwise_jsd_transition_nb(seqmat: np.ndarray, K: int, eps: float) -> np.ndarray:
        """Pairwise JSD transition distance - Numba accelerated."""
        n = seqmat.shape[0]
        D = np.zeros((n, n), dtype=np.float64)
        for i in prange(n):
            for j in range(i + 1, n):
                d = _jsd_transition_distance_nb(seqmat[i], seqmat[j], K, eps)
                D[i, j] = D[j, i] = d
        return D


# Duration-vector helpers


def _dur_vector(seq: np.ndarray, K: int) -> np.ndarray:
    v = np.bincount(seq, minlength=K).astype(float)
    s = v.sum()
    return v / s if s > 0 else v


def dur_l1_distance(a: np.ndarray, b: np.ndarray, K: int) -> float:
    return float(np.abs(_dur_vector(a, K) - _dur_vector(b, K)).sum())


def dur_cosine_distance(a: np.ndarray, b: np.ndarray, K: int) -> float:
    u = _dur_vector(a, K)
    v = _dur_vector(b, K)
    nu = np.linalg.norm(u)
    nv = np.linalg.norm(v)
    if nu == 0.0 or nv == 0.0:
        return 0.0
    cos = float(np.dot(u, v) / (nu * nv))
    return float(1.0 - cos)


def dur_chi2_distance(a: np.ndarray, b: np.ndarray, K: int, eps: float = 1e-12) -> float:
    u = _dur_vector(a, K)
    v = _dur_vector(b, K)
    denom = u + v + eps
    return float(np.sum((u - v) ** 2 / denom))


def jsd_transition_distance(a: np.ndarray, b: np.ndarray, K: int, eps: float = 0.0) -> float:
    """Jensen-Shannon divergence between the transition distributions of ``a`` and ``b``.

    Each distribution is the share of each ordered pair of consecutive states among a
    sequence's transitions. A sequence with fewer than two positions has no transition;
    the distance is then 0.
    """
    # per-sequence transition distribution KxK flattened
    if len(a) < 2 or len(b) < 2:
        return 0.0
    Ta = np.zeros((K, K), dtype=float)
    Tb = np.zeros((K, K), dtype=float)
    ia = a[:-1]
    ja = a[1:]
    np.add.at(Ta, (ia, ja), 1.0)
    ib = b[:-1]
    jb = b[1:]
    np.add.at(Tb, (ib, jb), 1.0)
    pa = Ta.ravel()
    pb = Tb.ravel()
    sa = pa.sum()
    sb = pb.sum()
    if sa > 0:
        pa = pa / sa
    if sb > 0:
        pb = pb / sb
    m = 0.5 * (pa + pb)

    def _kl(p, q):
        positive = p > 0
        return np.sum(p[positive] * (np.log(p[positive] + eps) - np.log(q[positive] + eps)))

    jsd = 0.5 * _kl(pa, m) + 0.5 * _kl(pb, m)
    return float(jsd)


def dhd_distance(
    a: np.ndarray, b: np.ndarray, freq_by_pos: np.ndarray, *, variant: str = "mean"
) -> float:
    """Dynamic Hamming Distance variants.

    For a mismatch at position t between states x,y with per-position frequencies p_t,
    we define:
      - variant="mean" (default): cost_t(x,y) = 1 - 0.5 * (p_t[x] + p_t[y]); d = mean_t cost_t
      - variant="sum": cost_t(x,y) = 2 - (p_t[x] + p_t[y]); d = (1/L) * sum_t cost_t
    The "sum" variant equals 2 * "mean" variant at each t; both average over positions for scale.
    """
    L = len(a)
    acc = 0.0
    variant_l = normalize_dhd_variant(variant)
    for t in range(L):
        if a[t] == b[t]:
            continue
        p = freq_by_pos[t]
        if variant_l == "mean":
            acc += 1.0 - 0.5 * (p[a[t]] + p[b[t]])
        else:  # sum
            acc += 2.0 - (p[a[t]] + p[b[t]])
    return float(acc / L if L > 0 else 0.0)


def _freq_by_pos(seqmat: np.ndarray, *, smooth_window: int | None = None) -> np.ndarray:
    smooth_window = validate_smooth_window(smooth_window)
    L = seqmat.shape[1]
    K = int(seqmat.max()) + 1 if seqmat.size else 0
    freq = np.zeros((L, K), dtype=float)
    for t in range(L):
        vals, cnt = np.unique(seqmat[:, t], return_counts=True)
        s = cnt.sum()
        if s:
            freq[t, vals] = cnt / s
    if smooth_window and smooth_window > 1 and L > 0:
        # simple moving average along t for each state, then renormalize rows
        w = np.ones(smooth_window, dtype=float) / smooth_window
        left_pad = (smooth_window - 1) // 2
        right_pad = smooth_window // 2
        padded = np.pad(freq, ((left_pad, right_pad), (0, 0)), mode="edge")
        smoothed = np.zeros_like(freq)
        for k in range(K):
            smoothed[:, k] = np.convolve(padded[:, k], w, mode="valid")
        row_sums = smoothed.sum(axis=1, keepdims=True)
        freq = np.divide(smoothed, row_sums, out=np.zeros_like(smoothed), where=row_sums > 0)
    return freq


# ---------------------------------------------------------------------------
# TraMineR-faithful Dynamic Hamming Distance (Lesnard 2010).
#
# TraMineR's seqdist(method="DHD") builds time-varying substitution costs from
# transition rates: sm <- seqsubm(seqdata, "TRATE", cval=4, time.varying=TRUE)
# (seqdist.R line 451), then runs position-varying Hamming (sum of per-position
# substitution costs). The cost transform is seqcost.R::tratecostBoth:
#
#   interior t:  C_t(i,j) = cval - tr[i,j,t-1] - tr[j,i,t-1] - tr[i,j,t] - tr[j,i,t]
#   boundary t:  C_t(i,j) = cval + 2*c           (c uses only the available side)
#   then clamped at 0; diagonal = 0.
#
# tr[i,j,s] = P(X_{s+1}=j | X_s=i) is the time-varying transition rate
# (seqtrate(time.varying=TRUE); zero-count rows -> 0, no smoothing). cval=4 for
# DHD. This is distinct from the legacy marginal-frequency variant below.
# ---------------------------------------------------------------------------
def _timevarying_trate(
    seqmat: np.ndarray,
    *,
    weights: np.ndarray | None = None,
    ignore_code: int | None = None,
    n_states: int | None = None,
) -> np.ndarray:
    """Time-varying transition rates ``tr[K, K, L-1]`` with ``tr[i,j,s]=P(X_{s+1}=j|X_s=i)``.

    Matches TraMineR ``seqtrate(time.varying=TRUE)``: a zero-denominator row (state ``i``
    absent at position ``s``) yields an all-zero row; no smoothing is applied.

    ``ignore_code`` (optional) drops any transition pair whose source or target equals
    that code — used to exclude the structural ``void`` symbol so void padding does not
    contaminate the rates (R's ``seqtrate`` treats void as NA). Defaults to ``None``
    (count every pair), so existing callers are unchanged.
    """
    n, L = seqmat.shape
    observed = seqmat if ignore_code is None else seqmat[seqmat != ignore_code]
    observed_states = int(observed.max()) + 1 if observed.size else 0
    K = observed_states if n_states is None else int(n_states)
    if K < observed_states:
        raise ValueError("n_states cannot be smaller than an observed state code plus one")
    if K <= 0 or L < 2:
        return np.zeros((max(K, 0), max(K, 0), max(L - 1, 0)), dtype=float)
    w = np.ones(n, dtype=float) if weights is None else np.asarray(weights, dtype=float)
    tr = np.zeros((K, K, L - 1), dtype=float)
    for s in range(L - 1):
        a = seqmat[:, s]
        b = seqmat[:, s + 1]
        ws = w
        if ignore_code is not None:
            keep = (a != ignore_code) & (b != ignore_code)
            a, b, ws = a[keep], b[keep], w[keep]
        counts = np.zeros((K, K), dtype=float)
        np.add.at(counts, (a, b), ws)
        denom = counts.sum(axis=1)
        nz = denom > 0
        tr[nz, :, s] = counts[nz, :] / denom[nz, None]
    return tr


def _dhd_cost_string_message(subm: str) -> str:
    """The refusal of a cost string DHD cannot use; R's text for CONSTANT (seqdist.R:253),
    followed by the strings DHD does resolve."""
    if subm.upper() == "CONSTANT":
        return (
            "'sm = \"CONSTANT\"' is not relevant for DHD, consider HAM instead; DHD's cost "
            "strings are TRATE, INDELS and INDELSLOG"
        )
    return "DHD subm must be TRATE, INDELS, INDELSLOG, or a time-varying array"


def _dhd_costs_traminer(
    seqmat: np.ndarray,
    *,
    cval: float = 4.0,
    weights: np.ndarray | None = None,
    missing_code: int | None = None,
    n_states: int | None = None,
    ignore_code: int | None = None,
) -> np.ndarray:
    """TraMineR DHD per-position substitution costs ``costs[L, K, K]`` (transition='both')."""
    n, L = seqmat.shape
    observed = seqmat if ignore_code is None else seqmat[seqmat != ignore_code]
    observed_states = int(observed.max()) + 1 if observed.size else 0
    K = observed_states if n_states is None else int(n_states)
    if K < observed_states:
        raise ValueError("n_states cannot be smaller than an observed state code plus one")
    costs = np.zeros((L, max(K, 0), max(K, 0)), dtype=float)
    if K <= 0 or L == 0:
        return costs
    tr = _timevarying_trate(
        seqmat,
        weights=weights,
        ignore_code=ignore_code,
        n_states=K,
    )
    for t in range(L):
        debut = t == 0
        fin = t == L - 1
        for i in range(K):
            for j in range(i + 1, K):
                c = 0.0
                if not debut:
                    c -= tr[i, j, t - 1] + tr[j, i, t - 1]
                if not fin:
                    c -= tr[i, j, t] + tr[j, i, t]
                val = (c + cval) if (not debut and not fin) else (cval + 2.0 * c)
                if val < 0.0:
                    val = 0.0
                costs[t, i, j] = costs[t, j, i] = val
    # seqcost(..., method="TRATE", time.varying=TRUE) overwrites every
    # missing<->observed entry with the fixed miss.cost (defaulting to cval).
    if missing_code is not None and 0 <= missing_code < K:
        costs[:, missing_code, :] = cval
        costs[:, :, missing_code] = cval
        costs[:, missing_code, missing_code] = 0.0
    return costs


def _dhd_indels_costs(
    seqmat: np.ndarray,
    *,
    method: str,
    weights: np.ndarray | None = None,
    n_states: int | None = None,
    ignore_code: int | None = None,
) -> np.ndarray:
    """TraMineR time-varying INDELS/INDELSLOG costs for DHD."""
    n, length = seqmat.shape
    observed = seqmat if ignore_code is None else seqmat[seqmat != ignore_code]
    observed_states = int(observed.max()) + 1 if observed.size else 0
    state_count = observed_states if n_states is None else int(n_states)
    if state_count < observed_states:
        raise ValueError("n_states cannot be smaller than an observed state code plus one")
    row_weights = np.ones(n, dtype=float) if weights is None else np.asarray(weights, dtype=float)
    if row_weights.ndim != 1 or row_weights.shape[0] != n:
        raise ValueError("weights length must equal number of sequences")

    costs = np.zeros((length, state_count, state_count), dtype=float)
    for time in range(length):
        states = seqmat[:, time]
        keep = np.ones(n, dtype=bool) if ignore_code is None else states != ignore_code
        counts = np.bincount(
            states[keep],
            weights=row_weights[keep],
            minlength=state_count,
        ).astype(float, copy=False)
        total = float(counts.sum())
        frequencies = counts / total if total > 0 else np.full(state_count, np.nan)
        frequencies = np.where(np.isnan(frequencies), 1.0, frequencies)
        if method == "INDELSLOG":
            indels = np.log(2.0 / (1.0 + frequencies))
        else:
            with np.errstate(divide="ignore"):
                indels = 1.0 / frequencies
            indels = np.where(np.isinf(indels), np.finfo(float).max, indels)
        with np.errstate(over="ignore"):
            costs[time] = np.add.outer(indels, indels)
        np.fill_diagonal(costs[time], 0.0)
    return costs


#: seqcost's INDELS stand-in for a state of frequency 0 (``.Machine$double.xmax``).
_INDELS_STAND_IN = float(np.finfo(float).max)


def _neutral_costs(costs: np.ndarray, kept: list[int]) -> float:
    """The smallest cost between two distinct kept states, 0 with fewer than two."""
    if len(kept) < 2:
        return 0.0
    block = costs[np.ix_(kept, kept)]
    return float(block[~np.eye(len(kept), dtype=bool)].min())


def _without_unobserved_stand_ins(
    subm: np.ndarray, indel: np.ndarray, present: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Keep INDELS stand-in costs of states no compared sequence contains out of the scale.

    INDELS gives a state of frequency 0 the stand-in indel ``DBL_MAX`` (seqcost.R). When no
    compared sequence contains that state its costs enter no raw distance, yet R lets them
    set the normalization bound (OMdistance.cpp:23-33; seqdist.R:644-649) and OMloc's local
    costs, so every normalized distance becomes 0 (TX-AUD-20261019-002). Such a state's costs
    are replaced by the smallest cost among the present states, which leaves every maximum
    over the present states, and so every distance, as with that state left out of the
    alphabet. A state a compared sequence contains keeps its cost, stand-in or not.
    """
    costs = np.asarray(subm, dtype=float)
    indels = np.asarray(indel, dtype=float)
    size = costs.shape[0]
    present_set = {int(code) for code in np.asarray(present).ravel() if 0 <= int(code) < size}
    absent = [
        state
        for state in range(size)
        if state not in present_set and indels.ndim == 1 and indels[state] >= _INDELS_STAND_IN
    ]
    if not absent or not present_set:
        return subm, indel
    kept = sorted(present_set)
    neutral_cost = _neutral_costs(costs, kept)
    neutral_indel = float(indels[kept].min())
    new_costs = costs.copy()
    new_indels = indels.copy()
    for state in absent:
        new_costs[state, :] = neutral_cost
        new_costs[:, state] = neutral_cost
        new_costs[state, state] = 0.0
        new_indels[state] = neutral_indel
    return new_costs, new_indels


def _dhd_without_unobserved_stand_ins(costs: np.ndarray, seqmat: np.ndarray) -> np.ndarray:
    """Per position, :func:`_without_unobserved_stand_ins` for DHD's time-varying INDELS costs.

    A state absent from position ``t`` of every compared sequence has the stand-in there;
    its costs at ``t`` become the smallest cost among the states present at ``t``, so DHD's
    bound sums, per position, the largest cost between states present there.
    """
    out = np.asarray(costs, dtype=float).copy()
    for time in range(out.shape[0]):
        present = sorted({int(code) for code in np.asarray(seqmat)[:, time]})
        present = [code for code in present if 0 <= code < out.shape[1]]
        size = out.shape[1]
        absent = [
            state
            for state in range(size)
            if state not in present
            and all(
                out[time, state, other] >= _INDELS_STAND_IN
                for other in range(size)
                if other != state
            )
        ]
        if not absent or not present:
            continue
        neutral_cost = _neutral_costs(out[time], present)
        for state in absent:
            out[time, state, :] = neutral_cost
            out[time, :, state] = neutral_cost
            out[time, state, state] = 0.0
    return out


def _validate_dhd_cost_array(
    subm,
    *,
    time_width: int,
    n_states: int,
) -> np.ndarray:
    """Validate the R-facing DHD ``(time, state, state)`` cost contract.

    TraMineR accepts an explicit DHD cost only when the corresponding R array
    is exactly ``states x states x ncol(seqdata)``.  Traminex exposes the same
    data with time first, so neither a static matrix nor a longer array may be
    broadcast or truncated at the ``seqdist`` boundary.
    """
    try:
        costs = np.asarray(subm, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "DHD subm must be a time-varying array matching the active states "
            "and original stored width"
        ) from exc

    expected = (int(time_width), int(n_states), int(n_states))
    # TraMineR (and seqcost(time_varying=True)["sm"]) lays the same data out
    # states x states x time; name that layout and the conversion whenever the
    # array reads as it (audit S1a-03).
    r_layout = (int(n_states), int(n_states), int(time_width))
    layout_hint = (
        f"; this looks like TraMineR's (K, K, T) layout {r_layout} -- seqdist's DHD "
        "costs are time-first (T, K, K): pass np.moveaxis(sm, -1, 0)"
    )
    looks_like_r_layout = (
        costs.ndim == 3
        and costs.shape == r_layout
        and bool(np.all(np.diagonal(costs, axis1=0, axis2=1) == 0))
    )
    if costs.ndim != 3 or costs.shape != expected:
        raise ValueError(
            "DHD subm must be a time-varying array matching the active states "
            f"and original stored width; expected {expected}, got {costs.shape}"
            + (layout_hint if looks_like_r_layout else "")
        )
    if not np.all(np.isfinite(costs)):
        raise ValueError("DHD substitution costs must be finite")
    if np.any(costs < 0):
        raise ValueError("DHD substitution costs cannot be negative")
    if np.any(np.diagonal(costs, axis1=1, axis2=2) != 0):
        raise ValueError(
            "DHD substitution-cost diagonals must be zero"
            + (layout_hint if looks_like_r_layout else "")
        )

    _warn_asymmetric_costs(costs, stacklevel=3)
    return costs


def _pairwise_dhd_traminer(
    seqmat: np.ndarray,
    *,
    cval: float = 4.0,
    weights: np.ndarray | None = None,
    missing_code: int | None = None,
    n_states: int | None = None,
    void_code: int | None = None,
    norm=None,
    costs: np.ndarray | None = None,
) -> np.ndarray:
    """TraMineR DHD pairwise distances: ``sum_t costs[t, a_t, b_t]`` (position-varying HAM).

    ``costs`` are the :func:`_dhd_costs_traminer` costs for these arguments when the
    caller has built them already, so one call builds them once (TX-REP-20261019-007).
    """
    n, L = seqmat.shape
    if costs is None:
        costs = _dhd_costs_traminer(
            seqmat,
            cval=cval,
            weights=weights,
            missing_code=missing_code,
            n_states=n_states,
            ignore_code=void_code,
        )
    D = np.zeros((n, n), dtype=float)
    if n == 0 or L == 0:
        return D
    rows = [row[row != void_code] if void_code is not None else row for row in seqmat]
    logical_length = len(rows[0]) if rows else 0
    if len({len(row) for row in rows}) > 1:
        raise ValueError("DHD is not defined for sequences of different length")
    tpos = np.arange(logical_length)
    for a in range(n):
        for b in range(a + 1, n):
            with np.errstate(over="ignore"):
                s = float(costs[tpos, rows[a], rows[b]].sum())
            D[a, b] = D[b, a] = s
    maxdist = _dhd_max_cost_sum(costs[:logical_length])
    return _normalize_fixed_position_output(D, maxdist, "DHD", norm)


def _dhd_max_cost_sum(costs: np.ndarray) -> float:
    """Sum per-position maxima without NumPy's benign overflow warning."""
    if costs.shape[0] == 0:
        return 0.0
    return float(sum(float(value) for value in np.max(costs, axis=(1, 2))))


def _normalize_fixed_position_output(
    distances: np.ndarray, maxdist: float, method: str, norm
) -> np.ndarray:
    """Normalize HAM/DHD output using their shared fixed-position bound.

    Each cell is normalized on its own, as in ``seqdist``'s square normalizer:
    an asymmetric reference distance fills one cell of its pair
    (TX-REP-20261017-001; batch-6 review C-1).
    """
    code = resolve_norm_code(method, norm)
    if code == 0:
        return distances
    out = np.asarray(distances, dtype=float).copy()
    for i in range(out.shape[0]):
        for j in range(i + 1, out.shape[1]):
            upper = normalize_value(out[i, j], maxdist, maxdist, maxdist, code)
            lower = (
                upper
                if out[j, i] == out[i, j]
                else normalize_value(out[j, i], maxdist, maxdist, maxdist, code)
            )
            out[i, j], out[j, i] = upper, lower
    return out


def _validate_low_level_dhd_costs(subm, *, length: int, observed_states: int) -> np.ndarray:
    """Check caller-supplied DHD costs as R's checkcost does; return them time-first.

    The low-level ``pairwise`` takes a static ``(K, K)`` matrix (R's HAM-style
    use, broadcast to every position) or a time-varying ``(L, K, K)`` array.
    R's checkcost (TraMineR-checkcost.R:40-75) refuses a wrong size, a negative
    cost, a NaN and a non-zero diagonal; this layer has no declared alphabet, so
    the size check asks the states to cover every observed code
    (TX-AUD-20261010-012).
    """
    try:
        costs = np.asarray(subm, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("DHD subm must be a 2D (K,K) or 3D (L,K,K) array") from exc
    if costs.ndim == 2:
        costs = np.broadcast_to(costs, (length, *costs.shape))
    elif costs.ndim == 3:
        if costs.shape[0] != length:
            raise ValueError(
                f"time-varying DHD subm must have first axis L={length}, got {costs.shape[0]}"
            )
    else:
        raise ValueError("DHD subm must be a 2D (K,K) or 3D (L,K,K) array")
    if costs.shape[1] != costs.shape[2]:
        raise ValueError("DHD subm must be square in its two state axes")
    if costs.shape[1] < observed_states:
        raise ValueError(
            f"DHD subm must cover the observed states: {costs.shape[1]}x{costs.shape[1]} "
            f"is smaller than {observed_states}x{observed_states}"
        )
    if not np.all(np.isfinite(costs)):
        raise ValueError("DHD substitution costs must be finite")
    if np.any(costs < 0):
        raise ValueError("DHD substitution costs cannot be negative")
    if np.any(np.diagonal(costs, axis1=1, axis2=2) != 0):
        raise ValueError("DHD substitution-cost diagonals must be zero")
    return costs


def _pairwise_dhd_explicit_subm(seqmat: np.ndarray, subm: np.ndarray) -> np.ndarray:
    """Raw DHD/HAM sums with a caller-supplied substitution matrix (``seqdist(sm=...)``).

    ``subm`` is either a static ``(K, K)`` matrix (constant over time, i.e. R
    ``seqdist(method="HAM", sm=subm)``) or a time-varying ``(L, K, K)`` array
    (R ``seqdist(method="DHD", sm=subm)``). The distance is the position-wise
    substitution sum ``sum_t subm_t[a_t, b_t]`` (a diagonal of zero gives the
    usual "cost only on mismatch" behaviour). Callers normalise the result with
    R's bound, ``_normalize_fixed_position_output(D, _dhd_max_cost_sum(costs), ...)``.
    """
    n, L = seqmat.shape
    M = np.asarray(subm, dtype=float)
    if M.ndim == 2:
        sm_t = np.broadcast_to(M, (L, M.shape[0], M.shape[1]))
    elif M.ndim == 3:
        if M.shape[0] != L:
            raise ValueError(f"time-varying DHD subm must have first axis L={L}, got {M.shape[0]}")
        sm_t = M
    else:
        raise ValueError("DHD subm must be a 2D (K,K) or 3D (L,K,K) array")
    tpos = np.arange(L)

    def distance(a: int, b: int) -> float:
        with np.errstate(over="ignore"):
            return float(sm_t[tpos, seqmat[a], seqmat[b]].sum())

    if _cost_array_is_asymmetric(sm_t):
        return _pairwise_apply_stable_unique(seqmat, distance)
    return _pairwise_apply(n, distance)


def _is_legacy_dhd(kwargs: dict) -> bool:
    """True if DHD should use the legacy marginal-frequency metric instead of TraMineR parity.

    Legacy (opt-in) is selected by an explicit ``variant=``, a truthy ``smooth_window=``,
    or ``dhd_profile`` other than ``"default"``. Bare ``method="DHD"`` is TraMineR-faithful.
    """
    return (
        "variant" in kwargs
        or kwargs.get("smooth_window") is not None
        or str(kwargs.get("dhd_profile", "default")).lower() != "default"
    )


def _coerce_norm_flag(norm) -> bool:
    """CHI2/EUCLID use a boolean norm (TraMineR: 'auto' -> TRUE, 'none' -> FALSE)."""
    if norm is None:
        return False
    if isinstance(norm, bool):
        return norm
    # R (seqdist.R:313-314) accepts only "auto"/"none" for CHI2/EUCLID.
    s = str(norm).lower()
    if s == "auto":
        return True
    if s == "none":
        return False
    raise ValueError(
        "CHI2/EUCLID norm is not supported; only True/'auto' or False/None/'none' are valid"
    )


def _user_interval_breaks(breaks: object, L: int) -> list[tuple[int, int]]:
    """Convert R's ``breaks`` list to 0-indexed inclusive pairs.

    ``seqdist-CHI2.R:40`` requires a list of ordered pairs ``c(start, end)``,
    both **1-based and inclusive**, and the intervals may overlap. When it is
    supplied, R ignores ``step`` and ``overlap`` (``seqdist.R:256-258``).
    """
    if isinstance(breaks, (str, bytes)) or not isinstance(breaks, Iterable):
        raise ValueError(
            "breaks should be a list of ordered pairs: "
            "breaks=[(start1, end1), (start2, end2), ...]"
        )
    out: list[tuple[int, int]] = []
    for index, pair in enumerate(breaks):
        if isinstance(pair, (str, bytes)) or not isinstance(pair, Iterable):
            raise ValueError(
                "breaks should be a list of ordered pairs: "
                "breaks=[(start1, end1), (start2, end2), ...]"
            )
        raw_values = list(pair)
        if len(raw_values) != 2:
            raise ValueError(
                f"breaks[{index}] must hold exactly two positions, got {len(raw_values)}"
            )
        # R documents breaks as ordered pairs of integers (seqdist.Rd:209-214) but
        # does not check them; int() truncated 2.7 to 2 (TX-AUD-20261011-048).
        for value in raw_values:
            if (
                isinstance(value, (bool, np.bool_))
                or not isinstance(value, (Real, np.number))
                or not np.isfinite(float(value))
                or float(value) % 1 != 0
            ):
                raise ValueError(
                    f"breaks[{index}] = {tuple(raw_values)!r} must hold two integer positions"
                )
        values = [int(value) for value in raw_values]
        start, end = values
        if start > end:
            raise ValueError(f"breaks[{index}] is not ordered: start {start} > end {end}")
        if start < 1 or end > L:
            raise ValueError(
                f"breaks[{index}] = ({start}, {end}) falls outside the sequence "
                f"positions 1..{L}"
            )
        out.append((start - 1, end - 1))
    if not out:
        raise ValueError("breaks must name at least one interval")
    return out


def _resolve_global_pdotj(global_pdotj: object, K: int) -> np.ndarray | None:
    """Validate and renormalise a numeric ``global.pdotj`` (seqdist-CHI2.R:62-67).

    R divides a supplied vector by its own sum, so an unnormalised vector gives
    exactly the answer its normalised form would. The string ``"obs"`` is
    resolved earlier, by ``seqdist``, which has the sequence object ``seqmeant``
    needs; reaching here with a string means a lower-level caller supplied one.
    """
    if global_pdotj is None:
        return None
    if isinstance(global_pdotj, str):
        raise ValueError(
            "global_pdotj='obs' is resolved by seqdist, which has the sequence "
            "object seqmeant needs; pass a vector of proportions here"
        )
    values = np.asarray(global_pdotj, dtype=float).ravel()
    if values.ndim != 1 or values.shape[0] != K:
        raise ValueError("global_pdotj length must equal n_states")
    if np.any(values < 0) or not values.sum() > 0:
        raise ValueError(
            "When a vector, global_pdotj should be non negative with at least "
            "one strictly positive element."
        )
    return values / values.sum()


def _resolve_legacy_dhd_settings(seqmat: np.ndarray, kwargs: dict) -> tuple[np.ndarray, int, str]:
    """Resolve the legacy-DHD variant/window once for every dispatcher.

    Shared by this module's ``pairwise`` and the compiled dispatcher in
    ``distances/__init__``, so the calibrated-profile contract cannot drift
    between the two (cycle-8 audit F-41).
    """
    norm = kwargs.get("norm", None)
    if not (norm is None or (isinstance(norm, str) and norm.strip().lower() == "none")):
        # The legacy marginal-frequency metric has no R counterpart and no R
        # bound; ignoring the request silently was the defect (TX-AUD-20261010-023).
        raise ValueError(
            "the legacy DHD (variant/smooth_window/calibrated profile) has no "
            f"normalisation; norm must be None or 'none', got {norm!r}"
        )
    variant = kwargs.get("variant", None)
    smooth_window = kwargs.get("smooth_window", None)
    profile = str(kwargs.get("dhd_profile", "default")).lower()
    if profile not in {"default", "calibrated"}:
        raise ValueError("dhd_profile must be 'default' or 'calibrated'")
    if profile == "calibrated":
        # Fails loudly without a file, a "dhd" section or with an unknown key;
        # an explicit argument (not None) overrides the file (TX-AUD-20261010-032).
        from ..config import _calibrated_profile_values

        explicit = {
            key: value
            for key, value in (("variant", variant), ("smooth_window", smooth_window))
            if value is not None
        }
        resolved = _calibrated_profile_values(
            "dhd",
            explicit,
            kwargs.get("calibrated_params_paths", kwargs.get("calibrated_params_path")),
        )
        variant = resolved.get("variant", None)
        smooth_window = resolved.get("smooth_window", None)
    if variant is None:
        variant = "mean"
    variant_l = normalize_dhd_variant(variant)
    var_code = 0 if variant_l == "mean" else 1
    freq = _freq_by_pos(np.asarray(seqmat), smooth_window=smooth_window)
    return freq, var_code, variant_l


def _interval_breaks(L: int, step: int, overlap: bool) -> list[tuple[int, int]]:
    """Generate 0-indexed inclusive interval breaks matching seqdist-CHI2.R.

    ``step`` is a positive integer (checked by :func:`_chi2_euclid_pairwise`).
    """
    if overlap and step % 2 != 0:
        # R (seqdist.R:259-260); an odd step cannot be shifted by step/2.
        raise ValueError("'step' must be even when 'overlap = TRUE'")
    if step == 1:
        return [(i, i) for i in range(L)]
    # seqdist-CHI2.R takes the one-interval shortcut only when step == ncol; a
    # longer step with overlap still builds [1, 1+step/2], [1, ncol] and the
    # shifted interval (audit S1b-10). Without overlap that is one interval.
    if step == L or (step > L and not overlap):
        return [(0, L - 1)]
    half = step // 2
    breaks: list[tuple[int, int]] = []
    if overlap:
        breaks.append((0, min(half, L - 1)))
    bb = list(range(1, L + 1, step)) + [L + 1]
    for i in range(1, len(bb)):
        s1, e1 = bb[i - 1], bb[i] - 1
        breaks.append((s1 - 1, min(e1, L) - 1))
        if overlap:
            ps, pe = breaks[-1]
            breaks.append((min(ps + half, L - 1), min(pe + half, L - 1)))
    return breaks


def _chi2_euclid_pairwise(
    seqmat: np.ndarray,
    *,
    metric: str,
    step: float,
    overlap: bool,
    norm: bool,
    weights: np.ndarray | None,
    global_pdotj: np.ndarray | None,
    n_states: int | None = None,
    breaks: object = None,
) -> np.ndarray:
    """CHI2 / EUCLID pairwise distance on per-interval state-distribution vectors.

    Faithful port of seqdist-CHI2.R: per interval each sequence contributes its state
    time-proportion vector; pdotj is the (weighted) overall state proportion. The
    distance is sqrt(sum_col (P[i]-P[j])^2 / pdotj_col); EUCLID replaces pdotj with a
    constant (maxd=2 if norm else 1). When norm, distances are divided by sqrt(nbreaks).
    """
    # R refuses a step that is not a non-negative integer whatever the windows
    # (seqdist.R:189-190), and CHI2 stops on 0 inside seq() unless breaks replaces
    # the step. Every route reaches this check; the low-level layer used to
    # truncate 2.5 to 2 and read 0 as the whole sequence (batch-6 review I-3, M-1).
    if not is_r_positive_integer(step) or (breaks is None and float(step) < 1):
        raise ValueError("'step' must be a positive integer")
    n, L = seqmat.shape
    observed_states = int(seqmat.max()) + 1 if seqmat.size else 0
    K = observed_states if n_states is None else int(n_states)
    if K < observed_states:
        raise ValueError("n_states cannot be smaller than an observed state code plus one")
    if K < 0:
        raise ValueError("n_states must be non-negative")
    if n == 0 or L == 0 or K == 0:
        return np.zeros((n, n), dtype=float)
    w = np.ones(n, dtype=float)
    if weights is not None:
        # seqdist's case-weight contract (TX-AUD-20261011-008): one finite,
        # non-negative value per sequence. All-zero weights stay accepted: R
        # returns an all-zero CHI2 matrix for them.
        from ..sequences import _validate_sequence_weights

        if np.ndim(weights) != 1 or np.shape(weights)[0] != n:
            raise ValueError("weights length must equal number of sequences")
        validated = _validate_sequence_weights(weights, n)
        # R never weights EUCLID (seqdist-CHI2.R:7 sets weighted <- FALSE), so its
        # state support counts every sequence once (TX-AUD-20261010-067).
        if metric != "euclid":
            w = validated
    if breaks is None:
        intervals = _interval_breaks(L, int(step), bool(overlap))
    else:
        intervals = _user_interval_breaks(breaks, L)
    global_pdotj = _resolve_global_pdotj(global_pdotj, K)
    col_vectors: list[np.ndarray] = []
    pdotj_vals: list[float] = []
    for s, e in intervals:
        seg = seqmat[:, s : e + 1]
        counts = np.zeros((n, K), dtype=float)
        for k in range(K):
            counts[:, k] = (seg == k).sum(axis=1)
        ndot = (w[:, None] * counts).sum(axis=0)  # weighted column totals
        rs = counts.sum(axis=1, keepdims=True)
        P = np.divide(counts, rs, out=np.zeros_like(counts), where=rs > 0)
        if metric == "euclid":
            maxd = 2.0 if norm else 1.0
            last = np.full(K, maxd, dtype=float)
            last[ndot == 0] = 0.0
        else:  # chi2
            if global_pdotj is not None:
                last = np.asarray(global_pdotj, dtype=float).copy()
                if last.ndim != 1 or last.shape[0] != K:
                    raise ValueError("global_pdotj length must equal n_states")
                last[ndot == 0] = 0.0
            else:
                snd = ndot.sum()
                last = ndot / snd if snd > 0 else np.zeros(K, dtype=float)
            pos = last[ndot != 0]
            if norm and pos.size > 1:
                two_min = np.sort(pos)[:2]
                last = last * (1.0 / two_min[0] + 1.0 / two_min[1])
        for k in range(K):
            col_vectors.append(P[:, k])
            pdotj_vals.append(last[k])
    allmat = np.stack(col_vectors, axis=1)
    pdotj = np.asarray(pdotj_vals, dtype=float)
    keep = pdotj > 0
    allmat = allmat[:, keep]
    inv = 1.0 / pdotj[keep]
    D = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            diff = allmat[i] - allmat[j]
            D[i, j] = D[j, i] = float(np.sqrt(np.sum(diff * diff * inv)))
    if norm:
        D = D / np.sqrt(len(intervals))
    return D


def _chi2_euclid_rows(
    seqs: list[np.ndarray],
    *,
    metric: str,
    step: float,
    overlap: bool,
    norm: bool,
    weights: np.ndarray | None,
    global_pdotj: np.ndarray | None,
    n_states: int | None = None,
    time_width: int | None = None,
    breaks: object = None,
) -> np.ndarray:
    """CHI2/EUCLID on canonical rows, with absent ragged positions ignored."""
    n = len(seqs)
    canonical_maxlen = max((len(row) for row in seqs), default=0)
    maxlen = canonical_maxlen if time_width is None else int(time_width)
    if maxlen < canonical_maxlen:
        raise ValueError("time_width cannot be shorter than a canonical sequence row")
    padded = np.full((n, maxlen), -1, dtype=np.int64)
    for i, row in enumerate(seqs):
        padded[i, : len(row)] = row
    return _chi2_euclid_pairwise(
        padded,
        metric=metric,
        step=step,
        overlap=overlap,
        norm=norm,
        weights=weights,
        global_pdotj=global_pdotj,
        n_states=n_states,
        breaks=breaks,
    )


def _pad_kweights(kweights, L: int) -> np.ndarray:
    """Validate kweights and pad/truncate them to length L (seqdist.R:265-273, 655).

    R refuses a non-numeric or negative vector (a NaN fails the same check),
    recycles a scalar over the sequence length, and zero-pads a short vector.
    """
    if kweights is None:
        return np.ones(L, dtype=float)
    refusal = "'kweights' must be a vector of positive numbers"
    try:
        src = np.asarray(kweights, dtype=float)
    except (TypeError, ValueError):
        raise ValueError(refusal) from None
    if src.ndim == 0:
        src = np.full(L, float(src))
    if src.ndim != 1 or np.any(np.isnan(src)) or np.any(src < 0):
        raise ValueError(refusal)
    kw = np.zeros(L, dtype=float)
    nmin = min(L, src.shape[0])
    kw[:nmin] = src[:nmin]
    return kw


# A number beyond the double range as ``(mantissa, exponent)``: mantissa * 2**exponent,
# the mantissa in [0.5, 1) or 0. Additions align both to the larger exponent, which
# is exact, so each sum rounds as it would with an unlimited exponent range.
_Extended = tuple[float, int]


def _x_norm(mantissa: float, exponent: int) -> _Extended:
    if mantissa == 0.0 or not math.isfinite(mantissa):
        return mantissa, 0
    fraction, shift = math.frexp(mantissa)
    return fraction, exponent + shift


def _x_add(left: _Extended, right: _Extended) -> _Extended:
    if left[0] == 0.0:
        return right
    if right[0] == 0.0:
        return left
    top = max(left[1], right[1])
    return _x_norm(math.ldexp(left[0], left[1] - top) + math.ldexp(right[0], right[1] - top), top)


def _x_sqrt(value: _Extended) -> float:
    """The square root as a double: Inf beyond the double range, 0 for a non-positive value."""
    mantissa, exponent = value
    if math.isnan(mantissa):
        return math.nan
    if mantissa <= 0.0:
        return 0.0
    if exponent % 2:
        mantissa, exponent = mantissa * 2.0, exponent - 1
    try:
        return math.ldexp(math.sqrt(mantissa), exponent // 2)
    except OverflowError:
        return math.inf


def _x_spell_power(durations: np.ndarray, tpow: float) -> tuple[np.ndarray, np.ndarray]:
    """Spell lengths to the power ``tpow`` as mantissas and exponents.

    A finite power is :func:`spell_length_power`'s (R's ``^``); one beyond the double
    range is computed with 50 significant digits.
    """
    mantissas = np.empty(len(durations), dtype=float)
    exponents = np.empty(len(durations), dtype=np.int64)
    for index, duration in enumerate(np.asarray(durations, dtype=float)):
        value = float(spell_length_power(np.array([duration]), tpow)[0])
        if math.isfinite(value):
            mantissas[index], exponents[index] = _x_norm(value, 0)
            continue
        with decimal.localcontext() as context:
            context.prec = 50
            context.Emax, context.Emin = decimal.MAX_EMAX, decimal.MIN_EMIN
            power = decimal.Decimal(duration) ** decimal.Decimal(tpow)
            exponent = (
                int((power.ln() / decimal.Decimal(2).ln()).to_integral_value(decimal.ROUND_FLOOR))
                + 1
            )
            fraction = float(power / decimal.Decimal(2) ** exponent)
        mantissas[index], exponents[index] = _x_norm(fraction, exponent)
    return mantissas, exponents


def _subseq_rows(
    input_rows: list[np.ndarray], method: str, kwargs: dict, *, n_states: int | None = None
) -> np.ndarray:
    """Pairwise distance for the subsequence-vector family (NMS, NMSMST, SVRspell).

    NMS counts subsequences on the state sequence (hard match, or soft via ``prox``).
    SVRspell/NMSMST count subsequences on the spell (DSS) representation weighted by spell
    durations^tpow; SVRspell allows soft matching via ``prox``. Distance (distMethod=2,
    norm="none"): d(i,j) = sqrt( sum_k kw[k] * (self_i[k] + self_j[k] - 2*common_ij[k]) ).
    The final sqrt lives in TraMineR's R wrapper (seqdist.R), not in the C++ kernels.

    R computes distinct sequences only, so identical rows are at distance 0. Where R's
    double sums overflow (NaN or Inf for a finite distance), the pair is recomputed with
    counts rescaled by powers of two (TX-AUD-20261019-029, TX-AUD-20261011-047).
    """
    with np.errstate(over="ignore", invalid="ignore"):
        return _subseq_rows_impl(input_rows, method, kwargs, n_states=n_states)


def _subseq_rows_impl(
    input_rows: list[np.ndarray], method: str, kwargs: dict, *, n_states: int | None
) -> np.ndarray:
    rows = [np.asarray(row) for row in input_rows]
    n = len(rows)
    L = max((len(row) for row in rows), default=0)
    K = (
        int(n_states)
        if n_states is not None
        else max((int(row.max()) + 1 for row in rows if len(row)), default=0)
    )
    kw = _pad_kweights(kwargs.get("kweights", None), L)
    tpow = float(kwargs.get("tpow", 1.0))
    prox = kwargs.get("prox", None)

    if method == "NMS":
        P = (
            _validate_proximity_matrix(prox, n_states=K, method=method)
            if prox is not None
            else None
        )
        seqs = rows

        def counts(i: int, j: int, ks: int) -> np.ndarray:
            return subseq_counts(seqs[i], seqs[j], ks, prox=P)

    elif method == "NMSMST":
        spells = [spell_encode(row) for row in rows]
        seqs = [s for s, _ in spells]
        durs = [spell_length_power(d, tpow) for _, d in spells]

        def counts(i: int, j: int, ks: int) -> np.ndarray:
            return nmsmst_counts(seqs[i], seqs[j], durs[i], durs[j], ks)

        def counts_scaled(i, j, ks, di, dj, out):
            return _nmsmst_counts(seqs[i], seqs[j], di, dj, ks, out)

    else:  # SVRSPELL
        soft = (
            _validate_proximity_matrix(prox, n_states=K, method=method)
            if prox is not None
            else np.eye(K, dtype=float)
        )
        spells = [spell_encode(row) for row in rows]
        seqs = [s for s, _ in spells]
        durs = [spell_length_power(d, tpow) for _, d in spells]

        def counts(i: int, j: int, ks: int) -> np.ndarray:
            return svrspell_counts(seqs[i], seqs[j], durs[i], durs[j], ks, soft)

        def counts_scaled(i, j, ks, di, dj, out):
            return _svrspell_counts(seqs[i], seqs[j], di, dj, ks, soft, out)

    # R adds kweights[k] * count[k] one length after another, skipping zero
    # weights, up to the longest compared row (NMSdistance.cpp,
    # SUBSEQdistance::distance; the spell rows for NMSMST and SVRspell). np.sum's
    # pairwise summation rounded differently (TX-REP-20261022-001); np.cumsum
    # accumulates in order.
    count_length = _longest(seqs)
    weighted = np.flatnonzero(kw[:count_length] != 0.0)

    def kw_sum(vec: np.ndarray) -> float:
        if weighted.size == 0:
            return 0.0
        counts = np.zeros(count_length, dtype=float)
        ln = min(len(vec), count_length)
        counts[:ln] = vec[:ln]
        return float(np.cumsum(kw[weighted] * counts[weighted])[-1])

    batched = _USE_BATCHED and method == "NMS"
    if batched:
        # NMS counts for many pairs at once (TX-AUD-20261010-022).
        cells = (_longest(seqs) + 1) ** 2

        lengths = np.asarray([len(row) for row in seqs], dtype=np.int64)

        def count_batch(left: np.ndarray, right: np.ndarray):
            def raw(lo: np.ndarray, ro: np.ndarray):
                sizes = np.minimum(lengths[lo], lengths[ro])
                if P is not None:
                    return _batched.soft_counts_raw(seqs, lo, ro, sizes, P)
                return _batched.subseq_counts_raw(seqs, lo, ro, sizes)

            return _in_batches(raw, left, right, cells, depth=4)

        everyone = np.arange(n, dtype=np.int64)
        self_vals = [kw_sum(vec) for vec in count_batch(everyone, everyone)]
    else:
        self_vals = [kw_sum(counts(i, i, len(seqs[i]))) for i in range(n)]

    finite_weights = bool(np.all(np.isfinite(kw[weighted])))
    spell_powers: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    scaled_self: dict[int, tuple[np.ndarray, np.ndarray]] = {}

    def scaled_counts(i: int, j: int) -> tuple[np.ndarray, np.ndarray]:
        """Counts of pair (i, j) as mantissas and exponents (counts beyond the double range)."""
        size = min(len(seqs[i]), len(seqs[j]))
        if method == "NMS":
            left, right = np.array([i], dtype=np.int64), np.array([j], dtype=np.int64)
            sizes = np.array([size], dtype=np.int64)
            if P is not None:
                found, exps = _batched.soft_counts_scaled(seqs, left, right, sizes, P)
            else:
                found, exps = _batched.subseq_counts_scaled(seqs, left, right, sizes)
            return found[0], exps[0]
        for r in (i, j):
            if r not in spell_powers:
                spell_powers[r] = _x_spell_power(spells[r][1], tpow)
        (mi, ei), (mj, ej) = spell_powers[i], spell_powers[j]
        top = max([int(ei.max()) if ei.size else 0, int(ej.max()) if ej.size else 0])
        out: list[int] = []
        found_counts = counts_scaled(
            i, j, size, np.ldexp(mi, ei - top), np.ldexp(mj, ej - top), out
        )
        degree = 1 if method == "NMSMST" else 2  # durations per count: min(t_i, t_j), or t_i t_j
        return found_counts, np.asarray(out, dtype=np.int64) + degree * top

    def weighted_sum(found: tuple[np.ndarray, np.ndarray]) -> _Extended:
        mantissas, exponents = found
        total: _Extended = (0.0, 0)
        for k in weighted.tolist():  # R's order, zero weights skipped
            if k < len(mantissas):
                # Normalize the count first, so weight * mantissa stays below the
                # double range; only the exponent moves, so the product rounds as before.
                fraction, exponent = _x_norm(float(mantissas[k]), int(exponents[k]))
                total = _x_add(total, _x_norm(float(kw[k]) * fraction, exponent))
        return total

    def extended(i: int, j: int) -> float:
        for r in (i, j):
            if r not in scaled_self:
                scaled_self[r] = scaled_counts(r, r)
        ival, jval = weighted_sum(scaled_self[i]), weighted_sum(scaled_self[j])
        aval = weighted_sum(scaled_counts(i, j))
        return _x_sqrt(_x_add(_x_add(ival, jval), (-2.0 * aval[0], aval[1])))

    def identical(i: int, j: int) -> float | None:
        return 0.0 if np.array_equal(rows[i], rows[j]) else None

    def finish(i: int, j: int, vec) -> float:
        sq = self_vals[i] + self_vals[j] - 2.0 * kw_sum(vec)
        # Counts or duration powers beyond the double range made the sums NaN or
        # Inf, as R's do; with finite weights the distance is finite or truly Inf.
        if not np.isfinite(sq) and finite_weights:
            return extended(i, j)
        # TraMineR applies sqrt in the R wrapper (seqdist.R) for NMS/NMSMST/SVRspell;
        # clamp tiny negative float noise of the squared form before the root, but
        # let a NaN (e.g. from an infinite kweight) through, as R does.
        if np.isnan(sq):
            return float("nan")
        return float(np.sqrt(sq)) if sq > 0.0 else 0.0

    if batched:
        return _pairwise_apply_batched(seqs, False, count_batch, finish, identical)

    def dist(i: int, j: int) -> float:
        answer = identical(i, j)
        if answer is not None:
            return answer
        ks = min(len(seqs[i]), len(seqs[j]))
        return finish(i, j, counts(i, j, ks))

    return _pairwise_apply(n, dist)


def _subseq_pairwise(seqmat: np.ndarray, method: str, kwargs: dict) -> np.ndarray:
    return _subseq_rows([seqmat[i] for i in range(seqmat.shape[0])], method, kwargs)


def _pairwise_rows_classic(
    seqs: list[np.ndarray],
    method: str,
    kwargs: dict,
    *,
    n_states: int | None = None,
) -> np.ndarray:
    """Classic TraMineR distances on already canonical variable-length rows."""
    n = len(seqs)
    mu = method.upper()
    code = resolve_norm_code(mu, kwargs.get("norm", None))

    if mu == "OM":
        subm = kwargs.get("subm", 2.0)
        if isinstance(subm, str):
            raise TypeError("OM subm must be a float or matrix")
        indel_arg = kwargs.get("indel", 1.0)
        indel_vec = None
        if np.isscalar(indel_arg):
            indel = _scalar_float(indel_arg)
            # The dense OM path always refused this; the row engine returned
            # negative distances, even between identical sequences (audit S1b-01).
            if (
                np.isnan(indel)
                or indel < 0
                or (indel == 0 and not kwargs.get("_allow_zero_indel", False))
            ):
                raise ValueError("OM indel cost must be strictly positive")
        else:
            indel_vec = np.asarray(indel_arg, dtype=float)
            if indel_vec.ndim != 1 or np.any(np.isnan(indel_vec)) or np.any(indel_vec <= 0):
                raise ValueError("OM indel vector must be one-dimensional and strictly positive")
            active_count = (
                int(n_states)
                if n_states is not None
                else max((int(row.max()) + 1 for row in seqs if len(row)), default=0)
            )
            if len(indel_vec) != active_count:
                raise ValueError("OM indel vector length must equal the active state count")
            indel = float(indel_vec.max()) if len(indel_vec) else 0.0
        if kwargs.get("indel_by_pos") is not None or kwargs.get("subm_by_pos") is not None:
            raise NotImplementedError(
                "position-dependent OM costs are not supported with ragged sequences"
            )
        maxscost = _maxscost(_subm_max_bound(subm), indel, code)
        asymmetric_costs = isinstance(subm, np.ndarray) and _cost_array_is_asymmetric(subm)
        if asymmetric_costs:
            _warn_asymmetric_costs(subm, stacklevel=4)

        def om_finish(i: int, j: int, raw) -> float:
            m, nn = len(seqs[i]), len(seqs[j])
            maxdist = abs(nn - m) * indel + maxscost * min(m, nn)
            return normalize_value(float(raw), maxdist, m * indel, nn * indel, code)

        if _USE_BATCHED and _worth_batching(n, 5):

            def om_raw_batch(left: np.ndarray, right: np.ndarray):
                return _in_batches(
                    lambda lo, ro: _batched.om_raw(
                        seqs, lo, ro, indel=indel, subm=subm, indel_by_state=indel_vec
                    ),
                    left,
                    right,
                    _longest(seqs),
                )

            return _pairwise_apply_batched(seqs, asymmetric_costs, om_raw_batch, om_finish)

        def om_dist(i: int, j: int) -> float:
            raw = om_distance(seqs[i], seqs[j], indel=indel, subm=subm, indel_by_state=indel_vec)
            return om_finish(i, j, raw)

        return (
            _pairwise_apply_stable_unique(seqs, om_dist)
            if asymmetric_costs
            else _pairwise_apply(n, om_dist)
        )

    kernels = {
        "LCS": lcs_distance,
        "LCP": lcp_distance,
        "RLCP": rlcp_distance,
    }
    fn = kernels.get(mu)
    if fn is None:
        raise ValueError(f"unsupported classic row distance: {method!r}")

    def classic_finish(i: int, j: int, raw) -> float:
        m, nn = len(seqs[i]), len(seqs[j])
        return normalize_value(float(raw), m + nn, m, nn, code)

    if mu == "LCS" and _USE_BATCHED and _worth_batching(n, 8):

        def lcs_raw_batch(left: np.ndarray, right: np.ndarray):
            return _in_batches(
                lambda lo, ro: _batched.lcs_raw(seqs, lo, ro), left, right, _longest(seqs)
            )

        return _pairwise_apply_batched(seqs, False, lcs_raw_batch, classic_finish)

    def classic_dist(i: int, j: int) -> float:
        return classic_finish(i, j, fn(seqs[i], seqs[j]))

    return _pairwise_apply(n, classic_dist)


def _refuse_nan_parameter(method: str, name: str, value: float) -> float:
    """Refuse a NaN method parameter before its sign check.

    ``value < 0`` is False for NaN, so the sign checks alone let NaN through and
    every distance became NaN; R stops on these parameters ("missing value
    where TRUE/FALSE needed", TX-AUD-20261010-013).
    """
    if math.isnan(value):
        raise ValueError(f"{method} '{name}' must not be NaN")
    return value


def _pairwise_rows_omslen(
    seqs: list[np.ndarray], kwargs: dict, *, n_states: int | None = None
) -> np.ndarray:
    """TraMineR OMslen on canonical variable-length rows."""
    n = len(seqs)
    state_count = (
        max((int(row.max()) + 1 for row in seqs if len(row)), default=0)
        if n_states is None
        else int(n_states)
    )
    indel, indel_max = _state_indel_spec(
        kwargs.get("indel", 1.0), n_states=state_count, method="OMslen"
    )
    subm = kwargs.get("subm", 2.0)
    link = str(kwargs.get("link", "mean")).lower()
    # R checks the link before computing, whatever the number of sequences
    # (seqdist.R:185-186, 229-233); the pair loop alone missed one sequence or
    # identical rows (TX-AUD-20261011-051).
    if link == "none":
        raise ValueError("'link' must be specified for OMslen: 'mean' or 'gmean'")
    if link not in {"mean", "gmean"}:
        raise ValueError("OMslen link must be 'mean' or 'gmean'")
    h = _refuse_nan_parameter("OMslen", "h", float(kwargs.get("h", 0.5)))
    if h < 0:
        raise ValueError("OMslen 'h' must be greater than or equal to 0")
    durs = [spell_length_per_position(row) for row in seqs]
    code = resolve_norm_code("OMslen", kwargs.get("norm", None))
    max_off = _subm_max_bound(subm)
    scost_off = max_off / 2.0 if link == "mean" else max_off
    maxscost = _maxscost(scost_off, indel_max, code)
    factor = 2.0 if link == "mean" else 1.0

    has_infinite_indel = (
        bool(np.isposinf(_scalar_float(indel)))
        if np.isscalar(indel)
        else bool(np.any(np.isposinf(np.asarray(indel, dtype=float))))
    )

    def special(i: int, j: int) -> float | None:
        if np.array_equal(seqs[i], seqs[j]):
            return 0.0
        if has_infinite_indel:
            return float("nan")
        return None

    def finish(i: int, j: int, raw) -> float:
        m, nn = len(seqs[i]), len(seqs[j])
        maxdist = abs(nn - m) * indel_max + factor * maxscost * min(m, nn)
        return normalize_value(float(raw), maxdist, m * indel_max, nn * indel_max, code)

    if _USE_BATCHED and _worth_batching(len(seqs), 5):
        weights = [spell_length_power(dur, -float(h)) for dur in durs]

        def raw_batch(left: np.ndarray, right: np.ndarray):
            return _in_batches(
                lambda lo, ro: _batched.omslen_raw(
                    seqs, lo, ro, indel=indel, subm=subm, link=link, weights=weights
                ),
                left,
                right,
                _longest(seqs),
            )

        return _pairwise_apply_for_costs_batched(seqs, subm, raw_batch, finish, special)

    def dist(i: int, j: int) -> float:
        answer = special(i, j)
        if answer is not None:
            return answer
        raw = omslen_distance(
            seqs[i],
            seqs[j],
            indel=indel,
            subm=subm,
            link=link,
            dur_a=durs[i],
            dur_b=durs[j],
            h=h,
        )
        return finish(i, j, raw)

    return _pairwise_apply_for_costs(seqs, subm, dist)


def _refuse_negative_spell_indels(
    seqs: list[np.ndarray], indel: float | np.ndarray, timecost: float, tpow: float
) -> None:
    """Refuse data in which a spell's indel cost ``indel[s] + timecost * (t**tpow - 1)`` is negative.

    With ``tpow < 0`` the cost falls as the spell lengthens (OMPerdistance.h), so each
    state's longest spell has its lowest cost. R computes on and returns negative
    distances, which no normalization bounds (TX-REP-20261022-003); a zero cost stays.
    """
    longest: dict[int, int] = {}
    for row in seqs:
        states, durations = spell_encode(row)
        for state, length in zip(states.tolist(), durations.tolist()):
            longest[int(state)] = max(longest.get(int(state), 0), int(length))
    for state, length in sorted(longest.items()):
        base = _scalar_float(indel) if np.isscalar(indel) else float(np.asarray(indel)[state])
        power = float(spell_length_power(np.array([float(length)]), tpow)[0])
        cost = base + timecost * (power - 1.0)
        if cost < 0:
            raise ValueError(
                "OMspell spell indel cost must be non-negative: with "
                f"tpow={tpow!r} a spell of length {length} (state code {state}) costs "
                f"indel + expcost * (t**tpow - 1) = {cost!r}; R returns negative "
                "distances there"
            )


def _pairwise_rows_omspell(
    seqs: list[np.ndarray], kwargs: dict, *, n_states: int | None = None
) -> np.ndarray:
    """TraMineR OMspell on canonical variable-length rows."""
    n = len(seqs)
    state_count = (
        max((int(row.max()) + 1 for row in seqs if len(row)), default=0)
        if n_states is None
        else int(n_states)
    )
    indel, indel_max = _state_indel_spec(
        kwargs.get("indel", 1.0), n_states=state_count, method="OMspell"
    )
    subm = kwargs.get("subm", 2.0)
    timecost = _refuse_nan_parameter("OMspell", "expcost", float(kwargs.get("expcost", 0.5)))
    if timecost < 0:
        # R (seqdist.R): "'expcost' must be positive" for OMloc and OMspell;
        # a negative expcost yields negative distances.
        raise ValueError("OMspell expcost must be non-negative")
    tpow = float(kwargs.get("tpow", 1.0))
    code = resolve_norm_code("OMspell", kwargs.get("norm", None))
    if tpow < 0 and timecost > 0:
        _refuse_negative_spell_indels(seqs, indel, timecost, tpow)
    bound_cache: dict[tuple[int, int], float] = {}

    has_infinite_indel = (
        bool(np.isposinf(_scalar_float(indel)))
        if np.isscalar(indel)
        else bool(np.any(np.isposinf(np.asarray(indel, dtype=float))))
    )

    def special(i: int, j: int) -> float | None:
        if np.array_equal(seqs[i], seqs[j]):
            return 0.0
        if has_infinite_indel:
            return float("nan")
        return None

    # gmean keeps R's bound, R's maxscost read from the upper triangle (owner decision of
    # 2026-10-22, TX-REP-20261022-002).
    r_maxscost = _maxscost(_subm_max_offdiag(subm), indel_max, code)

    def finish(i: int, j: int, raw) -> float:
        m, nn = len(seqs[i]), len(seqs[j])
        if code == 2:
            r_bound = abs(nn - m) * indel_max + r_maxscost * min(m, nn)
            return normalize_value(float(raw), r_bound, m * indel_max, nn * indel_max, code)
        if (m, nn) not in bound_cache:
            # A true bound, with the spell-duration costs (TX-REP-20261021-004).
            bound_cache[(m, nn)] = omspell_bound(
                m,
                nn,
                indel=indel_max,
                substitution=_subm_max_bound(subm),
                timecost=timecost,
                tpow=tpow,
                yujianbo=code == 4,
            )
        return normalize_value(
            float(raw), bound_cache[(m, nn)], m * indel_max, nn * indel_max, code
        )

    if _USE_BATCHED and _worth_batching(len(seqs), 6):
        spell_states, spell_durations = _batched.spell_rows(seqs, tpow)

        def raw_batch(left: np.ndarray, right: np.ndarray):
            return _in_batches(
                lambda lo, ro: _batched.omspell_raw(
                    spell_states, spell_durations, lo, ro, indel=indel, subm=subm, timecost=timecost
                ),
                left,
                right,
                _longest(spell_states),
            )

        return _pairwise_apply_for_costs_batched(seqs, subm, raw_batch, finish, special)

    def dist(i: int, j: int) -> float:
        answer = special(i, j)
        if answer is not None:
            return answer
        raw = omspell_distance(
            seqs[i],
            seqs[j],
            indel=indel,
            subm=subm,
            timecost=timecost,
            tpow=tpow,
        )
        return finish(i, j, raw)

    return _pairwise_apply_for_costs(seqs, subm, dist)


def _pairwise_rows_twed(
    seqs: list[np.ndarray], kwargs: dict, *, n_states: int | None = None
) -> np.ndarray:
    """TraMineR TWED on canonical variable-length rows."""
    n = len(seqs)
    if "nu" not in kwargs:
        raise ValueError("TWED requires 'nu' (a strictly positive stiffness parameter)")
    nu = _refuse_nan_parameter("TWED", "nu", float(kwargs["nu"]))
    if nu <= 0:
        raise ValueError("TWED 'nu' must be strictly greater than 0")
    h = _refuse_nan_parameter("TWED", "h", float(kwargs.get("h", 0.5)))
    if h < 0:
        raise ValueError("TWED 'h' must be greater than or equal to 0")
    subm = kwargs.get("subm", 2.0)
    code = resolve_norm_code("TWED", kwargs.get("norm", None))
    max_off = _subm_max_bound(subm)
    state_count = (
        max((int(row.max()) + 1 for row in seqs if len(row)), default=0)
        if n_states is None
        else int(n_states)
    )
    indel = _twed_indel_value(
        kwargs.get("indel", None), subm=subm, nu=nu, h=h, n_states=state_count
    )
    # R's bound caps the largest cost at two indels (OMdistance.cpp:23-34), but TWED's
    # matches and warps pay sm directly, so with a small indel normalized values
    # exceeded 1; and against an empty sequence R sets d = n * indel
    # (TWEDdistance.cpp:144-149), above its bound (TX-REP-20261021-005). The bound keeps
    # R's form (TWEDdistance.cpp:153) with the uncapped largest cost (at least 2 * indel
    # under YujianBo, R's convention), and is R's fixed distance with an empty sequence.
    bound_sub = max(max_off, 2.0 * indel) if code == 4 else max_off
    if code == 2:
        # gmean keeps R's bound, maxscost from the upper triangle capped at 2 * indel (owner
        # decision of 2026-10-22, TX-REP-20261022-002); against an empty sequence R's
        # normalizer gives 1 whatever the bound.
        bound_sub = _maxscost(_subm_max_offdiag(subm), indel, code)

    def finish(i: int, j: int, raw) -> float:
        m, nn = len(seqs[i]), len(seqs[j])
        if min(m, nn) == 0:
            maxdist = abs(nn - m) * indel
        else:
            maxdist = abs(nn - m) * (nu + h + bound_sub)
            maxdist += 2.0 * (bound_sub + nu) * min(m, nn)
        return normalize_value(float(raw), maxdist, m * indel, nn * indel, code)

    if _USE_BATCHED and _worth_batching(len(seqs), 9):
        # twed_distance's border value for an infinite indel (TWEDdistance.cpp).
        big = 1000.0 * _subm_max_offdiag(subm) if np.isposinf(indel) else float("inf")

        def raw_batch(left: np.ndarray, right: np.ndarray):
            return _in_batches(
                lambda lo, ro: _batched.twed_raw(
                    seqs, lo, ro, nu=nu, lam=h, subm=subm, indel=indel, big=big
                ),
                left,
                right,
                _longest(seqs),
            )

        return _pairwise_apply_for_costs_batched(seqs, subm, raw_batch, finish)

    def dist(i: int, j: int) -> float:
        raw = twed_distance(seqs[i], seqs[j], nu=nu, lam=h, subm=subm, indel=indel)
        return finish(i, j, raw)

    return _pairwise_apply_for_costs(seqs, subm, dist)


def _pairwise_rows_omstran(seqs: list[np.ndarray], kwargs: dict, *, n_states: int) -> np.ndarray:
    """TraMineR OMstran on canonical variable-length rows."""
    from ._omstran import omstran_transform_rows

    n = len(seqs)
    if "otto" not in kwargs:
        raise ValueError("OMstran requires 'otto' (origin weight in (0, 1])")
    otto = float(kwargs["otto"])
    if not (0.0 <= otto <= 1.0):
        raise ValueError("OMstran 'otto' must be in [0, 1]")
    transindel = str(kwargs.get("transindel", "constant"))
    previous = bool(kwargs.get("previous", False))
    add_column = bool(kwargs.get("add_column", True))
    rows, indels, new_sm = omstran_transform_rows(
        seqs,
        otto=otto,
        transindel=transindel,
        previous=previous,
        add_column=add_column,
        indel=kwargs.get("indel", 1.0),
        subm=kwargs.get("subm", 2.0),
        weights=kwargs.get("weights", None),
        n_states=n_states,
    )
    code = resolve_norm_code("OMstran", kwargs.get("norm", None))
    indel_norm = float(np.max(indels)) if len(indels) else 0.0
    maxscost = _maxscost(_subm_max_offdiag(new_sm), indel_norm, code)

    def finish(i: int, j: int, raw) -> float:
        m, nn = len(rows[i]), len(rows[j])
        maxdist = abs(nn - m) * indel_norm + maxscost * min(m, nn)
        return normalize_value(float(raw), maxdist, m * indel_norm, nn * indel_norm, code)

    if _USE_BATCHED and _worth_batching(len(seqs), 6):

        def raw_batch(left: np.ndarray, right: np.ndarray):
            return _in_batches(
                lambda lo, ro: _batched.om_raw(
                    rows, lo, ro, indel=1.0, subm=new_sm, indel_by_state=indels
                ),
                left,
                right,
                _longest(rows),
            )

        return _pairwise_apply_for_costs_batched(seqs, kwargs.get("subm", 2.0), raw_batch, finish)

    def dist(i: int, j: int) -> float:
        raw = om_distance(rows[i], rows[j], subm=new_sm, indel_by_state=indels)
        return finish(i, j, raw)

    return _pairwise_apply_for_costs(seqs, kwargs.get("subm", 2.0), dist)


def _omloc_distance(
    a: np.ndarray,
    b: np.ndarray,
    *,
    subm: float | np.ndarray,
    expcost: float,
    context: float,
    maxscost: float,
) -> float:
    """Raw OMloc recurrence used by TraMineR's ``OMVIdistance``."""
    m, n = len(a), len(b)
    if m == 0 or n == 0:
        raise ValueError("OMloc is not defined for empty sequences")
    if np.isscalar(subm):
        scalar = _scalar_float(subm)

        def sc(x: int, y: int) -> float:
            return 0.0 if x == y else scalar

    else:
        matrix = np.asarray(subm, dtype=float)

        def sc(x: int, y: int) -> float:
            return float(matrix[x, y])

    base = expcost * maxscost

    def gap(state: int, prev: int, current: int) -> float:
        return base + context * (sc(prev, state) + sc(current, state)) / 2.0

    dp = np.zeros((m + 1, n + 1), dtype=float)
    for i in range(1, m + 1):
        dp[i, 0] = dp[i - 1, 0] + gap(int(a[i - 1]), int(b[0]), int(b[0]))
    for j in range(1, n + 1):
        dp[0, j] = dp[0, j - 1] + gap(int(b[j - 1]), int(a[0]), int(a[0]))
    for j in range(1, n + 1):
        prev_b = int(b[max(j - 2, 0)])
        current_b = int(b[j - 1])
        for i in range(1, m + 1):
            prev_a = int(a[max(i - 2, 0)])
            current_a = int(a[i - 1])
            delete = dp[i - 1, j] + gap(current_a, prev_b, current_b)
            insert = dp[i, j - 1] + gap(current_b, prev_a, current_a)
            substitute = dp[i - 1, j - 1] + sc(current_a, current_b)
            dp[i, j] = min(delete, insert, substitute)
    return float(dp[m, n])


def _pairwise_rows_omloc(seqs: list[np.ndarray], kwargs: dict) -> np.ndarray:
    """TraMineR OMloc on canonical rows; separate from exported exploratory ``omloc``."""
    subm = kwargs.get("subm", 2.0)
    expcost = _refuse_nan_parameter("OMloc", "expcost", float(kwargs.get("expcost", 0.5)))
    context_arg = kwargs.get("context", None)
    context = (
        1.0 - 2.0 * expcost
        if context_arg is None
        else _refuse_nan_parameter("OMloc", "context", float(context_arg))
    )
    if expcost < 0 or context < 0:
        raise ValueError("OMloc expcost and context must be non-negative")
    code = resolve_norm_code("OMloc", kwargs.get("norm", None))
    max_off = _subm_max_offdiag(subm)
    # R's normalization indel (seqdist.R:723) and maxscost (OMdistance.cpp:23-34); the
    # latter also sets the local indels, so both stay R's for the raw distance.
    indel = _subm_max_all(subm) * expcost + context
    maxscost = _maxscost(max_off, indel, code)
    # R's bound prices indels at ``indel``, but one costs up to
    # timecost + context * max(sm) (OMVIdistance.h), so normalized values exceeded 1
    # (TX-REP-20261021-003). The bound is the alignment R's OM bound describes, with
    # the largest indel and substitution costs this recurrence can charge.
    local_indel = expcost * maxscost + context * _subm_max_all(subm)
    bound_sub = _maxscost(_subm_max_bound(subm), local_indel, code)

    def finish(i: int, j: int, raw) -> float:
        m, nn = len(seqs[i]), len(seqs[j])
        if code == 2:
            # gmean keeps R's bound (owner decision of 2026-10-22, TX-REP-20261022-002): its
            # lengths keep it non-negative only for a bound close to the true maximum.
            maxdist = abs(nn - m) * indel + maxscost * min(m, nn)
        else:
            maxdist = abs(nn - m) * local_indel + bound_sub * min(m, nn)
        return normalize_value(float(raw), maxdist, m * indel, nn * indel, code)

    if _USE_BATCHED and _worth_batching(len(seqs), 8):

        def raw_batch(left: np.ndarray, right: np.ndarray):
            return _in_batches(
                lambda lo, ro: _batched.omloc_raw(
                    seqs, lo, ro, subm=subm, expcost=expcost, context=context, maxscost=maxscost
                ),
                left,
                right,
                _longest(seqs),
            )

        return _pairwise_apply_for_costs_batched(seqs, kwargs.get("subm", 2.0), raw_batch, finish)

    def dist(i: int, j: int) -> float:
        raw = _omloc_distance(
            seqs[i],
            seqs[j],
            subm=subm,
            expcost=expcost,
            context=context,
            maxscost=maxscost,
        )
        return finish(i, j, raw)

    return _pairwise_apply_for_costs(seqs, kwargs.get("subm", 2.0), dist)


def _pairwise_rows(
    input_rows: list[np.ndarray],
    method: str,
    kwargs: dict | None = None,
    *,
    n_states: int | None = None,
    missing_code: int | None = None,
    weights: np.ndarray | None = None,
    reference_index: int | None = None,
    time_width: int | None = None,
) -> np.ndarray:
    """Compute an R-facing distance matrix from canonical variable-length rows.

    ``input_rows`` must already have structural void padding removed.  Unlike the
    historical dense :func:`pairwise` dispatcher, this adapter uses each pair's
    actual lengths for normalization and owns dataset-level preprocessing needed
    by CHI2/EUCLID, DHD, and OMstran.  A caller that appends external reference
    rows can therefore append zero reference weights and slice the returned cross
    block without changing the fitted distribution of the source dataset.

    ``n_states`` is the active alphabet width (including the retained missing
    state, when enabled). ``missing_code`` identifies that state for DHD's fixed
    missing-transition cost. ``reference_index`` identifies an appended
    external reference: R checks fixed-position source lengths before appending
    it and then compares only through the shorter row. These are internal
    metadata rather than public distance kwargs. ``time_width`` preserves the
    original rectangular column count after structural void deletion.
    """
    options = {} if kwargs is None else dict(kwargs)

    # Accept the existing internal transport spelling while keeping it out of
    # method-kwarg validation and the low-level kernels.
    transported_n_states = options.pop("_n_states", None)
    if n_states is None and transported_n_states is not None:
        n_states = int(transported_n_states)
    elif transported_n_states is not None:
        assert n_states is not None
        if int(transported_n_states) != int(n_states):
            raise ValueError("conflicting n_states values")
    transported_missing = options.pop("_missing_code", None)
    if missing_code is None and transported_missing is not None:
        missing_code = int(transported_missing)
    elif transported_missing is not None:
        assert missing_code is not None
        if int(transported_missing) != int(missing_code):
            raise ValueError("conflicting missing_code values")

    rows: list[np.ndarray] = []
    observed_states = 0
    for row in input_rows:
        raw = np.asarray(row)
        if raw.ndim != 1:
            raise ValueError("canonical sequence rows must be one-dimensional")
        if raw.dtype.kind not in {"i", "u"}:
            raise TypeError("canonical sequence rows must contain integer state codes")
        _validate_dense_sequence_matrix(raw[None, :])
        arr = raw.astype(np.int32, copy=False)
        if arr.size:
            if int(arr.min()) < 0:
                raise ValueError("canonical sequence rows cannot contain negative state codes")
            observed_states = max(observed_states, int(arr.max()) + 1)
        rows.append(arr)

    state_count = observed_states if n_states is None else int(n_states)
    if state_count < observed_states:
        raise ValueError("n_states cannot be smaller than an observed state code plus one")
    if state_count < 0:
        raise ValueError("n_states must be non-negative")
    if missing_code is not None and not (0 <= int(missing_code) < state_count):
        raise ValueError("missing_code must identify an active state")
    canonical_maxlen = max((len(row) for row in rows), default=0)
    if time_width is not None and int(time_width) < canonical_maxlen:
        raise ValueError("time_width cannot be shorter than a canonical sequence row")

    n = len(rows)
    if reference_index is not None and not (0 <= int(reference_index) < n):
        raise ValueError("reference_index must identify an appended sequence row")
    row_weights = None if weights is None else np.asarray(weights, dtype=float)
    if row_weights is not None and (row_weights.ndim != 1 or row_weights.shape[0] != n):
        raise ValueError("weights length must equal number of sequences")

    mu = method.upper()
    if mu in {"NMS", "NMSMST", "SVRSPELL"}:
        # R's public wrapper exposes norm="none" as the absence of a norm
        # request, while refusing every actual normalization for this family.
        norm = options.get("norm", None)
        if isinstance(norm, str) and norm.strip().lower() == "none":
            options.pop("norm")

    # Explicit row weights are dataset metadata.  They override an ordinary
    # method kwarg only for methods whose preprocessing consumes weights.
    if mu in {"CHI2", "OMSTRAN"} and row_weights is not None:
        options["weights"] = row_weights
    allow_zero_indel = bool(options.pop("_allow_zero_indel", False))
    validate_pairwise_kwargs(mu, options, dispatcher=False)
    if mu in _CHECKED_COST_METHODS:
        _refuse_invalid_fixed_costs(mu, options.get("subm", None), rows)
    if allow_zero_indel:
        options["_allow_zero_indel"] = True

    if mu in {"OM", "LCS", "LCP", "RLCP"}:
        return _pairwise_rows_classic(rows, mu, options, n_states=state_count)
    if mu == "OMLOC":
        return _pairwise_rows_omloc(rows, options)
    if mu == "OMSLEN":
        return _pairwise_rows_omslen(rows, options, n_states=state_count)
    if mu == "OMSPELL":
        return _pairwise_rows_omspell(rows, options, n_states=state_count)
    if mu == "OMSTRAN":
        return _pairwise_rows_omstran(rows, options, n_states=state_count)
    if mu == "TWED":
        return _pairwise_rows_twed(rows, options, n_states=state_count)
    if mu in {"NMS", "NMSMST", "SVRSPELL"}:
        return _subseq_rows(rows, mu, options, n_states=state_count)
    if mu in {"CHI2", "EUCLID"}:
        return _chi2_euclid_rows(
            rows,
            metric="euclid" if mu == "EUCLID" else "chi2",
            step=options.get("step", 1),
            overlap=bool(options.get("overlap", False)),
            norm=_coerce_norm_flag(options.get("norm", None)),
            # TraMineR forces ``weighted=FALSE`` for EUCLID, so an appended
            # external reference participates in interval support even though
            # its stored sequence weight is zero.
            weights=None if mu == "EUCLID" else options.get("weights", None),
            global_pdotj=options.get("global_pdotj", None),
            n_states=state_count,
            time_width=time_width,
            breaks=options.get("breaks", None),
        )

    if mu not in {"HAM", "DHD"}:
        raise ValueError(
            f"Unknown distance method: {method!r}. See documentation for supported methods."
        )

    source_indices = [
        index for index in range(n) if reference_index is None or index != reference_index
    ]
    source_lengths = {len(rows[index]) for index in source_indices}
    if len(source_lengths) > 1:
        raise ValueError(f"{mu} is not defined for sequences of different length")
    length = next(iter(source_lengths), 0)
    if reference_index is not None and len(rows[reference_index]) > length:
        raise ValueError("external refseq cannot be longer than the longest sequence")
    dense = (
        np.vstack([rows[index] for index in source_indices]).astype(np.int32, copy=False)
        if source_indices
        else np.empty((0, length), dtype=np.int32)
    )

    if mu == "HAM":
        subm = options.get("subm", None)
        if subm is None:
            distances = _pairwise_apply(
                n,
                lambda i, j: float(
                    np.count_nonzero(
                        rows[i][: min(len(rows[i]), len(rows[j]))]
                        != rows[j][: min(len(rows[i]), len(rows[j]))]
                    )
                ),
            )
            max_cost = 1.0
        else:
            matrix = np.asarray(subm, dtype=float)
            if matrix.ndim != 2:
                raise ValueError("HAM subm must be a 2D (K,K) array")
            asymmetric_costs = _cost_array_is_asymmetric(matrix)
            if asymmetric_costs:
                _warn_asymmetric_costs(matrix, stacklevel=4)

            def ham_distance_for_rows(i: int, j: int) -> float:
                return float(
                    matrix[
                        rows[i][: min(len(rows[i]), len(rows[j]))],
                        rows[j][: min(len(rows[i]), len(rows[j]))],
                    ].sum()
                )

            distances = (
                _pairwise_apply_stable_unique(rows, ham_distance_for_rows)
                if asymmetric_costs
                else _pairwise_apply(n, ham_distance_for_rows)
            )
            max_cost = _subm_max_all(matrix)
        return _normalize_fixed_position_output(
            distances,
            float(length) * max_cost,
            "HAM",
            options.get("norm", None),
        )

    if _is_legacy_dhd(options):
        # The opt-in exploratory DHD variants are intentionally delegated to
        # their existing dense implementation unchanged.
        if reference_index is not None:
            raise NotImplementedError(
                "short external references are supported only for R-faithful DHD"
            )
        return pairwise(dense, method="DHD", **options)

    subm = options.get("subm", None)
    norm = options.get("norm", None)
    source_weights = (
        None if row_weights is None else np.asarray(row_weights, dtype=float)[source_indices]
    )
    fit_width = length if time_width is None else int(time_width)
    if fit_width > length:
        fit_dense = np.full((dense.shape[0], fit_width), -1, dtype=np.int32)
        fit_dense[:, :length] = dense
        ignored_state = -1
    else:
        fit_dense = dense
        ignored_state = None

    cost_method = subm.upper() if isinstance(subm, str) else None
    if subm is None or cost_method == "TRATE":
        costs = _dhd_costs_traminer(
            fit_dense,
            cval=4.0,
            weights=source_weights,
            missing_code=missing_code,
            n_states=state_count,
            ignore_code=ignored_state,
        )
    elif cost_method in {"INDELS", "INDELSLOG"}:
        costs = _dhd_indels_costs(
            fit_dense,
            method=cost_method,
            weights=source_weights,
            n_states=state_count,
            ignore_code=ignored_state,
        )
    elif isinstance(subm, str):
        raise ValueError(_dhd_cost_string_message(subm))
    else:
        costs = _validate_dhd_cost_array(
            subm,
            time_width=fit_width,
            n_states=state_count,
        )

    logical_costs = costs[:length]
    if cost_method == "INDELS":
        # TX-AUD-20261019-002 (batch-7 review I-2): a state absent from a position of every
        # compared row, an external reference included, keeps its stand-in out of the bound.
        compared = np.full((n, length), -1, dtype=np.int32)
        for index, row in enumerate(rows):
            compared[index, : min(len(row), length)] = row[:length]
        logical_costs = _dhd_without_unobserved_stand_ins(logical_costs, compared)

    def dhd_pair_distance(i: int, j: int) -> float:
        pair_length = min(len(rows[i]), len(rows[j]))
        with np.errstate(over="ignore"):
            return float(
                logical_costs[
                    np.arange(pair_length),
                    rows[i][:pair_length],
                    rows[j][:pair_length],
                ].sum()
            )

    distances = (
        _pairwise_apply_stable_unique(rows, dhd_pair_distance)
        if _cost_array_is_asymmetric(logical_costs)
        else _pairwise_apply(n, dhd_pair_distance)
    )
    maxdist = _dhd_max_cost_sum(logical_costs)
    return _normalize_fixed_position_output(distances, maxdist, "DHD", norm)


def _as_cost_matrix(subm: object) -> object:
    """A substitution matrix given as nested lists or tuples, as a float array.

    Every route reads a matrix through ``np.ndarray`` checks, so a nested-list
    spelling took the symmetric path of the ragged route: it was applied in the
    reverse direction and without R's asymmetry warning (TX-AUD-20261019-003).
    A number -- Python, a NumPy scalar or a 0-d array -- becomes a Python float,
    since the routes test ``isinstance(subm, (int, float))`` and sent NumPy scalars
    down the matrix path (TX-AUD-20261019-028); a bool or a complex number is refused,
    as R refuses a logical or complex ``sm`` (batch-9 review M-5). Strings and arrays are
    returned unchanged.
    """
    if isinstance(subm, (list, tuple)):
        return np.asarray(subm, dtype=float)
    if is_refused_cost_scalar(subm):
        raise ValueError("'sm' invalid type or unknown value")
    if isinstance(subm, (Real, np.number)) or (isinstance(subm, np.ndarray) and subm.ndim == 0):
        return float(subm)
    return subm


def pairwise(seqmat: np.ndarray, method: str = "OM", **kwargs) -> np.ndarray:
    allow_zero_indel = bool(kwargs.pop("_allow_zero_indel", False))
    if "subm" in kwargs:
        kwargs["subm"] = _as_cost_matrix(kwargs["subm"])
    validate_pairwise_kwargs(method, kwargs, dispatcher=False)
    seqmat = _validate_dense_sequence_matrix(seqmat)
    n = seqmat.shape[0]
    mu = method.upper()
    if mu in _CHECKED_COST_METHODS - {"OM"}:
        # Dense OM validates its costs below, with its scalar/matrix dispatch.
        _refuse_invalid_fixed_costs(mu, kwargs.get("subm", None), list(seqmat))
    if mu == "HAM":
        # R-faithful Hamming (TraMineR seqdist method="HAM"): the RAW per-position
        # substitution cost sum (default norm="none"), NOT the length-normalized
        # proportion. The default substitution cost is unit (1 per mismatch); an
        # explicit ``subm`` matrix gives sum_t subm[a_t, b_t]. ``norm="maxlength"``
        # (== "auto" for fixed length) divides by L*maxcost.
        subm = kwargs.get("subm", None)
        norm = kwargs.get("norm", None)
        L = seqmat.shape[1]
        if subm is None:
            D = _pairwise_apply(n, lambda i, j: float(np.count_nonzero(seqmat[i] != seqmat[j])))
            maxcost = 1.0
        else:
            M = np.asarray(
                _validate_substitution_costs(
                    subm,
                    seqmat=seqmat,
                    method="HAM",
                    allow_scalar=False,
                ),
                dtype=float,
            )

            def ham_fn(i: int, j: int) -> float:
                return float(M[seqmat[i], seqmat[j]].sum())

            D = (
                _pairwise_apply_stable_unique(seqmat, ham_fn)
                if _cost_array_is_asymmetric(M)
                else _pairwise_apply(n, ham_fn)
            )
            maxcost = _subm_max_all(M)
        if norm in (None, "none"):
            return D
        if norm in ("maxlength", "auto", "length", "max"):
            denom = float(L) * maxcost
            return D / denom if denom > 0 else D
        raise ValueError(
            "HAM supports norm in {None,'none','maxlength','auto'}; "
            "'gmean'/'maxdist' are not implemented (use OM or norm='maxlength')"
        )
    if mu == "WHAM":
        # weighted Hamming with optional per-position weights; normalized by sum of weights
        w = kwargs.get("weights", None)
        if w is not None:
            w = np.asarray(w, dtype=float)
            if w.ndim != 1 or w.shape[0] != seqmat.shape[1]:
                raise ValueError(
                    "weights must be one-dimensional with length equal to sequence length"
                )
            if not np.all(np.isfinite(w)) or np.any(w < 0) or not np.any(w > 0):
                raise ValueError("weights must be finite, non-negative, and have a positive sum")
            with np.errstate(over="ignore"):
                s = float(w.sum())
            if not np.isfinite(s):
                w = w / float(w.max())
                s = float(w.sum())
            return _pairwise_apply(
                n, lambda i, j: float(np.dot(w, (seqmat[i] != seqmat[j]).astype(float))) / s
            )
        else:
            return _pairwise_apply(n, lambda i, j: ham_distance(seqmat[i], seqmat[j]))
    if mu == "OM":
        if np.ndim(kwargs.get("indel", 1.0)) == 1:
            # A per-state indel vector: the row engine computes it, as the
            # dispatcher's void route does (TX-AUD-20261011-050).
            from . import _pairwise_void_delete

            return _pairwise_void_delete(seqmat, -1, "OM", **kwargs)
        indel = float(kwargs.get("indel", 1.0))
        if np.isnan(indel) or indel < 0 or (indel == 0 and not allow_zero_indel):
            raise ValueError("OM indel cost must be strictly positive")
        subm = kwargs.get("subm", 2.0)
        if subm is None:
            subm = 2.0
        if isinstance(subm, str):
            raise TypeError("OM subm accepts float or ndarray in python kernels")
        subm = _validate_substitution_costs(
            subm,
            seqmat=seqmat,
            method="OM",
            allow_scalar=True,
        )
        indel_by_pos, subm_by_pos = _validate_om_position_costs(
            seqmat,
            indel_by_pos=kwargs.get("indel_by_pos", None),
            subm_by_pos=kwargs.get("subm_by_pos", None),
        )
        norm = kwargs.get("norm", None)
        if norm not in (None, "length", "max"):
            raise ValueError("OM supports norm=None|'length'|'max'")
        L = seqmat.shape[1]
        asymmetric_costs = (isinstance(subm, np.ndarray) and _cost_array_is_asymmetric(subm)) or (
            subm_by_pos is not None and _cost_array_is_asymmetric(subm_by_pos)
        )
        # Use Numba acceleration if available and no position-dependent costs
        plain_costs = indel_by_pos is None and subm_by_pos is None
        if (
            _USE_BATCHED
            and _worth_batching(n, 7)
            and plain_costs
            and not np.isposinf(indel)
            and (asymmetric_costs or not _HAVE_NUMBA)
        ):
            dense_rows = [seqmat[i] for i in range(n)]

            def _om_batch(left: np.ndarray, right: np.ndarray):
                return _in_batches(
                    lambda lo, ro: _batched.om_raw(dense_rows, lo, ro, indel=indel, subm=subm),
                    left,
                    right,
                    L,
                )

            D = _pairwise_apply_batched(
                dense_rows, asymmetric_costs, _om_batch, lambda i, j, raw: float(raw)
            )
        elif asymmetric_costs:
            D = _pairwise_apply_stable_unique(
                seqmat,
                lambda i, j: om_distance(
                    seqmat[i],
                    seqmat[j],
                    indel,
                    subm,
                    indel_by_pos=indel_by_pos,
                    subm_by_pos=subm_by_pos,
                ),
            )
        elif np.isposinf(indel) and indel_by_pos is None:
            D = _pairwise_apply(
                n,
                lambda i, j: om_distance(
                    seqmat[i],
                    seqmat[j],
                    indel,
                    subm,
                    subm_by_pos=subm_by_pos,
                ),
            )
        elif _HAVE_NUMBA and indel_by_pos is None and subm_by_pos is None:
            if isinstance(subm, (int, float)):
                D = _pairwise_om_scalar_nb(seqmat, indel, float(subm))
            elif isinstance(subm, np.ndarray):
                D = _pairwise_om_matrix_nb(seqmat, indel, np.asarray(subm, dtype=np.float64))
            else:
                # Fallback to Python for other subm types
                D = _pairwise_apply(
                    n,
                    lambda i, j: om_distance(
                        seqmat[i],
                        seqmat[j],
                        indel,
                        subm,
                        indel_by_pos=indel_by_pos,
                        subm_by_pos=subm_by_pos,
                    ),
                )
        else:
            D = _pairwise_apply(
                n,
                lambda i, j: om_distance(
                    seqmat[i],
                    seqmat[j],
                    indel,
                    subm,
                    indel_by_pos=indel_by_pos,
                    subm_by_pos=subm_by_pos,
                ),
            )
        if norm == "length" and L > 0:
            D = D / float(L)
        elif norm == "max":
            if subm_by_pos is not None:
                # Position-varying costs are a Python extension without an R
                # counterpart; the bound stays the largest cost anywhere.
                max_sub = float(np.max(subm_by_pos)) if subm_by_pos.size else 0.0
            else:
                # The whole off-diagonal maximum, not R's upper triangle
                # (TX-AUD-20261019-008).
                max_sub = _subm_max_bound(subm)
            if indel_by_pos is not None:
                max_insdel = float(np.max(indel_by_pos)) * 2.0 if len(indel_by_pos) > 0 else 0.0
            else:
                max_insdel = 2.0 * indel
            # Max possible per-position OM cost is min(max_sub, 2*indel): the DP never pays
            # more than an indel pair when a substitution is dearer (R OMdistance.cpp:
            # maxscost = fmin2(maxscost, 2*indel); maxpossiblecost = maxscost * min(m,n)).
            denom = float(L) * min(max_sub, max_insdel)
            if denom > 0:
                D = D / denom
        return D
    if mu == "OMLOC":
        return _pairwise_rows_omloc([seqmat[i] for i in range(n)], kwargs)
    if mu == "LCS":
        norm = kwargs.get("norm", None)
        if norm not in (None, "length"):
            raise ValueError("LCS supports norm=None|'length'")
        if _HAVE_NUMBA:
            D = _pairwise_lcs_nb(seqmat)
        elif _USE_BATCHED and _worth_batching(n, 8):
            # Without Numba (the core install), the row engine ragged LCS uses: the
            # per-pair reference was about 30 times slower (TX-AUD-20261019-022).
            dense_rows = [seqmat[i] for i in range(n)]

            def _lcs_batch(left: np.ndarray, right: np.ndarray):
                return _in_batches(
                    lambda lo, ro: _batched.lcs_raw(dense_rows, lo, ro),
                    left,
                    right,
                    seqmat.shape[1],
                )

            D = _pairwise_apply_batched(dense_rows, False, _lcs_batch, lambda i, j, raw: float(raw))
        else:
            D = _pairwise_apply(n, lambda i, j: lcs_distance(seqmat[i], seqmat[j]))
        if norm == "length":
            L = seqmat.shape[1]
            if L > 0:
                D = D / float(2 * L)
        return D
    if mu == "LCP":
        norm = kwargs.get("norm", None)
        if norm not in (None, "length"):
            raise ValueError("LCP supports norm=None|'length'")
        if _HAVE_NUMBA:
            D = _pairwise_lcp_nb(seqmat)
        else:
            D = _pairwise_apply(n, lambda i, j: lcp_distance(seqmat[i], seqmat[j]))
        if norm == "length":
            L = seqmat.shape[1]
            if L > 0:
                D = D / float(2 * L)
        return D
    if mu == "RLCP":
        norm = kwargs.get("norm", None)
        if norm not in (None, "length"):
            raise ValueError("RLCP supports norm=None|'length'")
        if _HAVE_NUMBA:
            D = _pairwise_rlcp_nb(seqmat)
        else:
            D = _pairwise_apply(n, lambda i, j: rlcp_distance(seqmat[i], seqmat[j]))
        if norm == "length":
            L = seqmat.shape[1]
            if L > 0:
                D = D / float(2 * L)
        return D
    if mu == "DUR":
        K = int(seqmat.max()) + 1 if seqmat.size else 0
        metric = str(kwargs.get("metric", "l1")).lower()
        if metric == "l1":
            return _pairwise_apply(n, lambda i, j: dur_l1_distance(seqmat[i], seqmat[j], K))
        if metric == "cos":
            return _pairwise_apply(n, lambda i, j: dur_cosine_distance(seqmat[i], seqmat[j], K))
        if metric == "chi2":
            if _HAVE_NUMBA:
                return _pairwise_chi2_nb(seqmat, K)
            return _pairwise_apply(n, lambda i, j: dur_chi2_distance(seqmat[i], seqmat[j], K))
        raise ValueError("Unsupported DUR metric: %r" % metric)
    if mu == "JSD_TRANS":
        K = int(seqmat.max()) + 1 if seqmat.size else 0
        # Only positive shares enter the logarithm (a > 0, m >= a / 2), so no
        # constant is needed: with the default 0 a pair of disjoint transition
        # distributions is ln 2 apart, the divergence's definition (batch-15 review M-2).
        eps = float(kwargs.get("eps", 0.0))
        if not np.isfinite(eps) or eps < 0:
            raise ValueError("eps must be finite and non-negative")
        if _HAVE_NUMBA:
            return _pairwise_jsd_transition_nb(seqmat, K, eps)
        return _pairwise_apply(
            n, lambda i, j: jsd_transition_distance(seqmat[i], seqmat[j], K, eps=eps)
        )
    if mu == "DHD":
        subm = kwargs.get("subm", None)
        if isinstance(subm, str):
            # The cost strings are resolved by seqdist against the corpus; this
            # layer takes an array. Without this guard the string fell through to
            # np.asarray(subm, dtype=float) and produced a bare "could not
            # convert string to float", naming neither the parameter nor the
            # layer -- the same gap cycle-6 F-009 closed for `indel`. The OM
            # branch already states its own accepted set; DHD did not.
            raise TypeError(
                f"DHD subm={subm!r} is not an array; pairwise accepts a (K, K) or "
                "(L, K, K) array of substitution costs. The cost strings TRATE, "
                "INDELS and INDELSLOG are resolved by seqdist, which has the "
                "corpus to derive them from."
            )
        if subm is not None:
            # Caller-supplied substitution costs: R seqdist(method="DHD"/"HAM", sm=...),
            # checked as checkcost does and normalised by R's bound, the sum over
            # positions of the largest cost (seqdist.R:645-649; TX-AUD-20261010-012/-023).
            costs = _validate_low_level_dhd_costs(
                subm,
                length=seqmat.shape[1],
                observed_states=int(seqmat.max()) + 1 if seqmat.size else 0,
            )
            _warn_asymmetric_costs(costs, stacklevel=3)
            raw = _pairwise_dhd_explicit_subm(seqmat, costs)
            return _normalize_fixed_position_output(
                raw, _dhd_max_cost_sum(costs), "DHD", kwargs.get("norm", None)
            )
        # Default DHD is TraMineR-faithful (transition-rate time-varying costs, cval=4).
        # The legacy marginal-frequency metric is opt-in (variant/smooth_window/calibrated).
        if not _is_legacy_dhd(kwargs):
            if str(kwargs.get("dhd_profile", "default")).lower() not in {"default", "calibrated"}:
                raise ValueError("dhd_profile must be 'default' or 'calibrated'")
            # Normalised as seqdist does, by R's bound (TX-AUD-20261010-023).
            return _pairwise_dhd_traminer(seqmat, cval=4.0, norm=kwargs.get("norm", None))
        freq, var_code, variant_l = _resolve_legacy_dhd_settings(seqmat, kwargs)
        if _HAVE_NUMBA:
            return _pairwise_dhd_nb(seqmat, freq, var_code)
        return _pairwise_apply(
            n, lambda i, j: dhd_distance(seqmat[i], seqmat[j], freq, variant=variant_l)
        )
    if mu in ("NMS", "NMSMST", "SVRSPELL"):
        return _subseq_pairwise(seqmat, mu, kwargs)
    if mu == "OMSTRAN":
        from ._omstran import omstran_transform

        if "otto" not in kwargs:
            raise ValueError("OMstran requires 'otto' (origin weight in (0, 1])")
        otto = float(kwargs["otto"])
        if not (0.0 <= otto <= 1.0):
            raise ValueError("OMstran 'otto' must be in [0, 1]")
        new_seqmat, indels, new_sm = omstran_transform(
            seqmat,
            otto=otto,
            transindel=str(kwargs.get("transindel", "constant")),
            previous=bool(kwargs.get("previous", False)),
            add_column=bool(kwargs.get("add_column", True)),
            indel=kwargs.get("indel", 1.0),
            subm=kwargs.get("subm", 2.0),
            weights=kwargs.get("weights", None),
        )
        code = resolve_norm_code("OMstran", kwargs.get("norm", None))

        def _raw_stran(i: int, j: int) -> float:
            return om_distance(new_seqmat[i], new_seqmat[j], subm=new_sm, indel_by_state=indels)

        # After the transform OMstran runs OM with a per-token indel vector (OMVIdistance.cpp):
        # the normalization indel is max(indels) (seqdist.R:696) and maxscost comes from the
        # transformed substitution matrix. maxpossiblecost = |n-m|*indel + maxscost*min(m,n);
        # lengths l = (augmented length)*indel.
        indel_norm = float(np.asarray(indels, dtype=float).max()) if len(indels) else 0.0
        maxscost = _maxscost(_subm_max_offdiag(new_sm), indel_norm, code)

        def _stran_finish(i: int, j: int, raw) -> float:
            if code == 0:
                return float(raw)
            m_i = len(new_seqmat[i])
            n_j = len(new_seqmat[j])
            maxdist = abs(n_j - m_i) * indel_norm + maxscost * min(m_i, n_j)
            return normalize_value(float(raw), maxdist, m_i * indel_norm, n_j * indel_norm, code)

        if _USE_BATCHED and _worth_batching(new_seqmat.shape[0], 6):
            augmented = [new_seqmat[i] for i in range(new_seqmat.shape[0])]

            def _stran_batch(left: np.ndarray, right: np.ndarray):
                return _in_batches(
                    lambda lo, ro: _batched.om_raw(
                        augmented, lo, ro, indel=1.0, subm=new_sm, indel_by_state=indels
                    ),
                    left,
                    right,
                    _longest(augmented),
                )

            return _pairwise_apply_for_costs_batched(
                list(seqmat), kwargs.get("subm", 2.0), _stran_batch, _stran_finish
            )

        def _stran_dist(i: int, j: int) -> float:
            return _stran_finish(i, j, _raw_stran(i, j))

        return _pairwise_apply_for_costs(list(seqmat), kwargs.get("subm", 2.0), _stran_dist)
    if mu == "OMSLEN":
        state_count = int(seqmat.max()) + 1 if seqmat.size else 0
        return _pairwise_rows_omslen([seqmat[i] for i in range(n)], kwargs, n_states=state_count)
    if mu == "OMSPELL":
        state_count = int(seqmat.max()) + 1 if seqmat.size else 0
        return _pairwise_rows_omspell([seqmat[i] for i in range(n)], kwargs, n_states=state_count)
    if mu == "TWED":
        state_count = int(seqmat.max()) + 1 if seqmat.size else 0
        return _pairwise_rows_twed([seqmat[i] for i in range(n)], kwargs, n_states=state_count)
    if mu in ("CHI2", "EUCLID"):
        # R parity (seqdist-CHI2.R): default step=1 (one interval per position).
        # EUCLID never reads the weights (its bound replaces pdotj), matching
        # R's forced weighted=FALSE; only CHI2 consumes them (audit F-49).
        return _chi2_euclid_pairwise(
            seqmat,
            metric="euclid" if mu == "EUCLID" else "chi2",
            step=kwargs.get("step", 1),
            overlap=bool(kwargs.get("overlap", False)),
            norm=_coerce_norm_flag(kwargs.get("norm", None)),
            weights=kwargs.get("weights", None),
            global_pdotj=kwargs.get("global_pdotj", None),
            breaks=kwargs.get("breaks", None),
        )
    raise ValueError(
        f"Unknown distance method: {method!r}. See documentation for supported methods."
    )
