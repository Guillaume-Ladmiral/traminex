"""Raw sequence distances for many pairs at once (TX-AUD-20261010-022).

The row engine used to run one pure-Python dynamic program per pair, 80-200
times slower than TraMineR for OM with per-state indels, OMslen, OMloc, OMstran,
TWED and NMS. Each function here evaluates the recurrence of its per-pair
reference in :mod:`traminex.distances.python_kernels` (named in its docstring)
for a batch of pairs: the table is swept cell by cell -- for the NMS counts,
anti-diagonal by anti-diagonal -- and every cell is computed by the same
floating-point expression, operand for operand, as the reference, only on
arrays holding one value per pair. Per-row quantities (spell weights and
durations) are computed by the reference's own calls, once per row, and
Python's ``min`` is reproduced exactly (:func:`_pymin`), so the results are
bit-identical to the reference (``tests/test_distance_batched_engine.py``).

Rows are padded to the longest row of the batch; a pair's cells beyond its own
lengths are computed from padding but never read by its own cells, and each
pair's result is taken at its own lengths.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence

import numpy as np

from ._spells import spell_encode, spell_length_power

#: Upper bound on the cells one batch holds in an (pairs x width) table.
_BATCH_CELLS = 4_000_000


def batch_size(width: int, depth: int = 1) -> int:
    """Pairs per batch for tables of ``width`` cells per pair (``depth`` tables)."""
    return max(1, _BATCH_CELLS // max(1, (width + 1) * max(1, depth)))


def _pymin(first: np.ndarray, *others: np.ndarray) -> np.ndarray:
    """Python's ``min(first, *others)`` per element.

    ``min`` keeps the running value unless a later one compares smaller, which
    differs from ``np.minimum`` when a NaN is involved.
    """
    result = first
    for other in others:
        result = np.where(other < result, other, result)
    return result


def padded(
    rows: Sequence[np.ndarray], index: np.ndarray, *, fill: float = 0, dtype=np.int64
) -> tuple[np.ndarray, np.ndarray]:
    """Rows ``index`` stacked into a ``(pairs, width)`` array, and their lengths."""
    lengths = np.fromiter((len(rows[k]) for k in index), dtype=np.int64, count=len(index))
    width = int(lengths.max()) if lengths.size else 0
    out = np.full((len(index), width), fill, dtype=dtype)
    for position, k in enumerate(index):
        out[position, : lengths[position]] = rows[k]
    return out, lengths


def _substitution(subm, *, scalar_test: Callable[[object], bool], scalar_value) -> Callable:
    """The reference's ``sc(x, y)``: 0 on equal states or the scalar, else ``M[x, y]``."""
    if scalar_test(subm):
        value = scalar_value(subm)

        def sc(x: np.ndarray, y: np.ndarray) -> np.ndarray:
            return np.where(x == y, 0.0, value)

        return sc
    matrix = np.asarray(subm, dtype=float)

    def sc_matrix(x: np.ndarray, y: np.ndarray) -> np.ndarray:
        return matrix[x, y]

    return sc_matrix


def _is_python_number(value: object) -> bool:
    return isinstance(value, (int, float))


def _sweep(
    a_len: np.ndarray,
    b_len: np.ndarray,
    width_a: int,
    width_b: int,
    first_row: np.ndarray,
    first_column: np.ndarray,
    cell: Callable[[int, int, np.ndarray, np.ndarray], np.ndarray],
) -> np.ndarray:
    """Sweep a ``(len_a + 1) x (len_b + 1)`` table for every pair; return its corner.

    ``first_row`` is ``dp[0, :]`` and ``first_column`` ``dp[:, 0]`` (``(pairs,
    width + 1)`` each); ``cell(i, j, previous_row, current_row)`` returns
    ``dp[i, j]`` from ``dp[i - 1, :]`` and the cells of row ``i`` left of ``j``.
    """
    pairs = a_len.shape[0]
    rows = np.arange(pairs)
    result = np.empty(pairs, dtype=float)
    done = a_len == 0
    result[done] = first_row[rows[done], b_len[done]]
    previous = first_row
    for i in range(1, width_a + 1):
        current = np.empty((pairs, width_b + 1), dtype=float)
        current[:, 0] = first_column[:, i]
        for j in range(1, width_b + 1):
            current[:, j] = cell(i, j, previous, current)
        done = a_len == i
        result[done] = current[rows[done], b_len[done]]
        previous = current
    return result


def _cumulative(costs: np.ndarray) -> np.ndarray:
    """``out[:, 0] = 0`` and ``out[:, k] = out[:, k - 1] + costs[:, k - 1]``, in order."""
    out = np.zeros((costs.shape[0], costs.shape[1] + 1), dtype=float)
    for k in range(1, costs.shape[1] + 1):
        out[:, k] = out[:, k - 1] + costs[:, k - 1]
    return out


def om_raw(
    rows: Sequence[np.ndarray],
    left: np.ndarray,
    right: np.ndarray,
    *,
    indel: float,
    subm,
    indel_by_state: np.ndarray | None = None,
) -> np.ndarray:
    """Batched ``python_kernels.om_distance`` (no position-dependent costs)."""
    if indel_by_state is None and np.isposinf(indel):
        return np.asarray(
            [0.0 if np.array_equal(rows[i], rows[j]) else float("inf") for i, j in zip(left, right)]
        )
    a, a_len = padded(rows, left)
    b, b_len = padded(rows, right)
    if indel_by_state is not None:
        states = np.asarray(indel_by_state, dtype=float)
        indel_a, indel_b = states[a], states[b]
    else:
        indel_a = np.full(a.shape, indel, dtype=float)
        indel_b = np.full(b.shape, indel, dtype=float)
    sc = _substitution(subm, scalar_test=_is_python_number, scalar_value=float)

    def cell(i: int, j: int, previous: np.ndarray, current: np.ndarray) -> np.ndarray:
        return _pymin(
            previous[:, j] + indel_a[:, i - 1],
            current[:, j - 1] + indel_b[:, j - 1],
            previous[:, j - 1] + sc(a[:, i - 1], b[:, j - 1]),
        )

    return _sweep(
        a_len, b_len, a.shape[1], b.shape[1], _cumulative(indel_b), _cumulative(indel_a), cell
    )


def lcs_raw(rows: Sequence[np.ndarray], left: np.ndarray, right: np.ndarray) -> np.ndarray:
    """Batched ``python_kernels.lcs_distance``: ``|a| + |b| - 2 * LCS(a, b)``."""
    a, a_len = padded(rows, left)
    b, b_len = padded(rows, right)
    pairs = a.shape[0]
    index = np.arange(pairs)
    lcs = np.zeros(pairs, dtype=np.int64)
    previous = np.zeros((pairs, b.shape[1] + 1), dtype=np.int64)
    for i in range(a.shape[1]):
        current = np.zeros_like(previous)
        for j in range(b.shape[1]):
            current[:, j + 1] = np.where(
                a[:, i] == b[:, j],
                previous[:, j] + 1,
                np.maximum(previous[:, j + 1], current[:, j]),
            )
        done = a_len == i + 1
        lcs[done] = current[index[done], b_len[done]]
        previous = current
    return (a_len + b_len - 2 * lcs).astype(float)


def twed_raw(
    rows: Sequence[np.ndarray],
    left: np.ndarray,
    right: np.ndarray,
    *,
    nu: float,
    lam: float,
    subm,
    indel: float,
    big: float,
) -> np.ndarray:
    """Batched ``python_kernels.twed_distance`` for a resolved ``indel``.

    ``big`` is the reference's border value: ``1000 * max off-diagonal cost``
    when ``indel`` is infinite, else ``inf``.
    """
    a, a_len = padded(rows, left)
    b, b_len = padded(rows, right)
    pairs = a.shape[0]
    result = np.empty(pairs, dtype=float)
    empty = (a_len == 0) | (b_len == 0)
    result[empty] = [float(int(m + n) * indel) for m, n in zip(a_len[empty], b_len[empty])]
    keep = ~empty
    if not keep.any():
        return result
    a, a_len, b, b_len = a[keep], a_len[keep], b[keep], b_len[keep]
    sc = _substitution(subm, scalar_test=_is_python_number, scalar_value=float)
    zero = np.zeros(a.shape[0], dtype=np.int64)
    width_b = b.shape[1]
    first_row = np.zeros((a.shape[0], width_b + 1), dtype=float)
    first_column = np.zeros((a.shape[0], a.shape[1] + 1), dtype=float)
    if np.isposinf(indel):
        first_row[:, 0] = np.nan
        first_column[:, 0] = np.nan
    border = np.full(a.shape[0], big, dtype=float)

    def cell(i: int, j: int, previous: np.ndarray, current: np.ndarray) -> np.ndarray:
        i_state, i_m1 = a[:, i - 1], (zero if i == 1 else a[:, i - 2])
        j_state, j_m1 = b[:, j - 1], (zero if j == 1 else b[:, j - 2])
        same = (i_state == j_state) & (i_m1 == j_m1)
        cost = np.where(same, 0.0, sc(i_state, j_state) + sc(i_m1, j_m1))
        sub = previous[:, j - 1] + cost + 2.0 * nu * abs(i - j)
        if j - 1 == 0:
            i_warp = border
            if i - 1 > 0:
                sub = border
        else:
            i_warp = current[:, j - 1] + sc(j_state, j_m1) + nu + lam
        if i - 1 == 0:
            j_warp = border
            if j - 1 > 0:
                sub = border
        else:
            j_warp = previous[:, j] + sc(i_state, i_m1) + nu + lam
        return _pymin(i_warp, j_warp, sub)

    result[keep] = _sweep(a_len, b_len, a.shape[1], width_b, first_row, first_column, cell)
    return result


def omslen_raw(
    rows: Sequence[np.ndarray],
    left: np.ndarray,
    right: np.ndarray,
    *,
    indel,
    subm,
    link: str,
    weights: Sequence[np.ndarray],
) -> np.ndarray:
    """Batched ``python_kernels.omslen_distance``.

    ``weights[k]`` is row ``k``'s ``spell_length_power(durations, -h)``, computed by the
    caller with the reference's own call.
    """
    a, a_len = padded(rows, left)
    b, b_len = padded(rows, right)
    w_a, _ = padded(weights, left, fill=1.0, dtype=float)
    w_b, _ = padded(weights, right, fill=1.0, dtype=float)
    sc = _substitution(subm, scalar_test=_is_python_number, scalar_value=float)
    if np.isscalar(indel):
        scalar = float(np.asarray(indel).item())  # the reference's _scalar_float
        cost_a, cost_b = scalar * w_a, scalar * w_b
    else:
        values = np.asarray(indel, dtype=float)
        cost_a, cost_b = values[a] * w_a, values[b] * w_b
    if link == "mean":

        def sub_cost(x, y, wi, wj):
            return (sc(x, y) / 2.0) * (wi + wj)

    elif link == "gmean":

        def sub_cost(x, y, wi, wj):
            return sc(x, y) * np.sqrt(wi * wj)

    else:
        raise ValueError("OMslen link must be 'mean' or 'gmean'")

    def cell(i: int, j: int, previous: np.ndarray, current: np.ndarray) -> np.ndarray:
        return _pymin(
            previous[:, j] + cost_a[:, i - 1],
            current[:, j - 1] + cost_b[:, j - 1],
            previous[:, j - 1] + sub_cost(a[:, i - 1], b[:, j - 1], w_a[:, i - 1], w_b[:, j - 1]),
        )

    return _sweep(
        a_len, b_len, a.shape[1], b.shape[1], _cumulative(cost_b), _cumulative(cost_a), cell
    )


def spell_rows(
    rows: Sequence[np.ndarray], tpow: float
) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """Each row's DSS states and ``durations ** tpow - 1``, as ``omspell_distance`` derives them."""
    states, durations = [], []
    for row in rows:
        spell_states, spell_durations = spell_encode(row)
        states.append(spell_states)
        durations.append(spell_length_power(spell_durations, float(tpow)) - 1.0)
    return states, durations


def omspell_raw(
    states: Sequence[np.ndarray],
    durations: Sequence[np.ndarray],
    left: np.ndarray,
    right: np.ndarray,
    *,
    indel,
    subm,
    timecost: float,
) -> np.ndarray:
    """Batched ``python_kernels.omspell_distance`` on rows from :func:`spell_rows`."""
    a, a_len = padded(states, left)
    b, b_len = padded(states, right)
    d_a, _ = padded(durations, left, fill=0.0, dtype=float)
    d_b, _ = padded(durations, right, fill=0.0, dtype=float)
    sc = _substitution(subm, scalar_test=_is_python_number, scalar_value=float)
    base_a: np.ndarray
    base_b: np.ndarray
    if np.isscalar(indel):
        scalar = float(np.asarray(indel).item())  # the reference's _scalar_float
        base_a = np.full(a.shape, scalar, dtype=float)
        base_b = np.full(b.shape, scalar, dtype=float)
    else:
        values = np.asarray(indel, dtype=float)
        base_a, base_b = values[a], values[b]
    gap_a = base_a + timecost * d_a
    gap_b = base_b + timecost * d_b

    def cell(i: int, j: int, previous: np.ndarray, current: np.ndarray) -> np.ndarray:
        x, y = a[:, i - 1], b[:, j - 1]
        di, dj = d_a[:, i - 1], d_b[:, j - 1]
        sub = np.where(x == y, timecost * np.abs(di - dj), sc(x, y) + (di + dj) * timecost)
        return _pymin(
            previous[:, j] + gap_a[:, i - 1],
            current[:, j - 1] + gap_b[:, j - 1],
            previous[:, j - 1] + sub,
        )

    return _sweep(
        a_len, b_len, a.shape[1], b.shape[1], _cumulative(gap_b), _cumulative(gap_a), cell
    )


def omloc_raw(
    rows: Sequence[np.ndarray],
    left: np.ndarray,
    right: np.ndarray,
    *,
    subm,
    expcost: float,
    context: float,
    maxscost: float,
) -> np.ndarray:
    """Batched ``python_kernels._omloc_distance``."""
    a, a_len = padded(rows, left)
    b, b_len = padded(rows, right)
    if (a_len == 0).any() or (b_len == 0).any():
        raise ValueError("OMloc is not defined for empty sequences")
    sc = _substitution(
        subm, scalar_test=np.isscalar, scalar_value=lambda value: float(np.asarray(value).item())
    )
    base = expcost * maxscost

    def gap(state, previous_state, current_state):
        return base + context * (sc(previous_state, state) + sc(current_state, state)) / 2.0

    first_a, first_b = a[:, 0], b[:, 0]
    border_a = (
        np.stack([gap(a[:, i], first_b, first_b) for i in range(a.shape[1])], axis=1)
        if a.shape[1]
        else np.zeros((a.shape[0], 0))
    )
    border_b = (
        np.stack([gap(b[:, j], first_a, first_a) for j in range(b.shape[1])], axis=1)
        if b.shape[1]
        else np.zeros((b.shape[0], 0))
    )

    def cell(i: int, j: int, previous: np.ndarray, current: np.ndarray) -> np.ndarray:
        previous_b, current_b = b[:, max(j - 2, 0)], b[:, j - 1]
        previous_a, current_a = a[:, max(i - 2, 0)], a[:, i - 1]
        delete = previous[:, j] + gap(current_a, previous_b, current_b)
        insert = current[:, j - 1] + gap(current_b, previous_a, current_a)
        substitute = previous[:, j - 1] + sc(current_a, current_b)
        return _pymin(delete, insert, substitute)

    return _sweep(
        a_len, b_len, a.shape[1], b.shape[1], _cumulative(border_b), _cumulative(border_a), cell
    )


def subseq_counts_raw(
    rows: Sequence[np.ndarray], left: np.ndarray, right: np.ndarray, sizes: np.ndarray
) -> list[np.ndarray]:
    """Batched hard-matching ``python_kernels.subseq_counts`` (``prox=None``)."""
    return _subseq_counts(rows, left, right, sizes, None)


def subseq_counts_scaled(
    rows: Sequence[np.ndarray], left: np.ndarray, right: np.ndarray, sizes: np.ndarray
) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """:func:`subseq_counts_raw` for counts beyond the double range.

    After each subsequence length every pair's table is divided by a power of two,
    so count ``k`` of a pair is ``counts[k] * 2**exponents[k]``. Dividing by a power
    of two is exact, so the mantissas are those of the same arithmetic with an
    unlimited exponent range (TX-AUD-20261011-047).
    """
    exponents: list[np.ndarray] = []
    counts = _subseq_counts(rows, left, right, sizes, exponents)
    return counts, exponents


def _subseq_counts(
    rows: Sequence[np.ndarray],
    left: np.ndarray,
    right: np.ndarray,
    sizes: np.ndarray,
    exponents_out: list[np.ndarray] | None,
) -> list[np.ndarray]:
    """Hard-matching subsequence counts; ``exponents_out`` given, rescaled per length.

    Returns one count vector of length ``sizes[p]`` per pair. Each table is
    stored by anti-diagonal (cell ``(i, j)`` at ``[i + j, i]``), so a diagonal
    of the recurrence is one slice; the totals add the matched cells in the
    reference's row-major order, and pairs that have stopped are dropped.
    """
    a, a_len = padded(rows, left, fill=-1)
    b, b_len = padded(rows, right, fill=-2)
    pairs, height, width = a.shape[0], a.shape[1], b.shape[1]
    sizes = np.asarray(sizes, dtype=np.int64)
    counts = np.zeros((pairs, int(sizes.max()) if pairs else 0), dtype=float)
    exponents = np.zeros(counts.shape, dtype=np.int64)
    shift = np.zeros(pairs, dtype=np.int64)
    if exponents_out is not None:
        exponents_out.extend(exponents[p, : sizes[p]] for p in range(pairs))
    started = (a_len > 0) & (b_len > 0) & (sizes > 0)
    # A slice in which every pair has an empty side has no count to fill, and its
    # table may have no column at all (every size 0), as in _soft_counts.
    if pairs == 0 or height == 0 or width == 0 or not started.any():
        return [counts[p, : sizes[p]] for p in range(pairs)]
    match = a[:, :, None] == b[:, None, :]  # padding codes -1 and -2 never match
    flat_match = match.reshape(pairs, height * width)
    matched = flat_match.sum(axis=1)
    total = matched.astype(float)  # whole ones: exact in any order
    counts[started, 0] = total[started]
    running = started & (a_len > 1) & (b_len > 1)
    if not running.any():
        return [counts[p, : sizes[p]] for p in range(pairs)]

    diagonals = height + width - 1
    rows_i, columns_j = np.indices((height, width))
    skewed = ((rows_i + columns_j) * height + rows_i).ravel()
    # Each pair's matches in row-major order (the order in which the reference
    # adds them), as positions in the skewed table; padding points at cell
    # (diagonals - 1, 0), which lies outside the table and stays 0.
    pair_of, cell_of = np.nonzero(flat_match)
    rank = np.arange(pair_of.size) - np.repeat(np.cumsum(matched) - matched, matched)
    gather = np.full((pairs, int(matched.max())), (diagonals - 1) * height, dtype=np.int64)
    gather[pair_of, rank] = skewed[cell_of]
    # Tables are (diagonal, row, pair): a diagonal of every pair is one block.
    is_match = np.zeros((diagonals * height, pairs), dtype=bool)
    is_match[skewed, :] = flat_match.T
    table = np.zeros((diagonals, height, pairs), dtype=float)
    # The next k's weights: the table at the matches, 0 elsewhere. The masked
    # copy below never writes an unmatched cell, so those stay 0.
    weights = is_match.reshape(table.shape).astype(float)
    sweep = [
        (d, max(0, d - (width - 2)), min(height - 2, d) + 1)
        for d in range(height + width - 4, -1, -1)
    ]
    ids = np.arange(pairs)
    flat_index = gather * pairs + ids[:, None]
    k = 1
    while True:
        active = running & (k < sizes[ids]) & (total > 0)
        if not active.any():
            break
        if 2 * int(active.sum()) <= ids.size:
            keep = np.flatnonzero(active)
            ids, active, gather = ids[keep], active[keep], gather[keep]
            # Contiguous copies keep the pair axis innermost and make the
            # reshapes below views, which the masked copy writes through.
            table = np.ascontiguousarray(table[:, :, keep])
            weights = np.ascontiguousarray(weights[:, :, keep])
            is_match = np.ascontiguousarray(is_match[:, keep])
            flat_index = gather * ids.size + np.arange(ids.size)[:, None]
        for d, lo, hi in sweep:  # the reference's (((v[i+1,j] + v[i,j+1]) - v[i+1,j+1]) + h)
            cell = table[d, lo:hi]
            np.add(table[d + 1, lo + 1 : hi + 1], table[d + 1, lo:hi], out=cell)
            np.subtract(cell, table[d + 2, lo + 1 : hi + 1], out=cell)
            np.add(cell, weights[d + 2, lo + 1 : hi + 1], out=cell)
        proceed = active & (table[0, 0] != 0)
        flat_weights = weights.reshape(diagonals * height, ids.size)
        np.copyto(flat_weights, table.reshape(flat_weights.shape), where=is_match)
        total = np.cumsum(table.take(flat_index), axis=1)[:, -1]  # sequential, row-major
        counts[ids[proceed], k] = total[proceed]
        if exponents_out is not None:
            exponents[ids[proceed], k] = shift[ids[proceed]]
            scale = np.frexp(np.max(np.abs(weights), axis=(0, 1)))[1].astype(np.int64)
            np.ldexp(weights, -scale, out=weights)
            shift[ids] += scale
        running = proceed
        k += 1
    return [counts[p, : sizes[p]] for p in range(pairs)]


def soft_counts_raw(
    rows: Sequence[np.ndarray],
    left: np.ndarray,
    right: np.ndarray,
    sizes: np.ndarray,
    prox: np.ndarray,
) -> list[np.ndarray]:
    """Batched soft-matching ``python_kernels.subseq_counts`` (``prox`` given)."""
    return _soft_counts(rows, left, right, sizes, prox, None)


def soft_counts_scaled(
    rows: Sequence[np.ndarray],
    left: np.ndarray,
    right: np.ndarray,
    sizes: np.ndarray,
    prox: np.ndarray,
) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """:func:`soft_counts_raw` for counts beyond the double range, as :func:`subseq_counts_scaled`."""
    exponents: list[np.ndarray] = []
    counts = _soft_counts(rows, left, right, sizes, prox, exponents)
    return counts, exponents


def _soft_counts(
    rows: Sequence[np.ndarray],
    left: np.ndarray,
    right: np.ndarray,
    sizes: np.ndarray,
    prox: np.ndarray,
    exponents_out: list[np.ndarray] | None,
) -> list[np.ndarray]:
    """Soft-matching subsequence counts; ``exponents_out`` given, rescaled per order.

    Mirrors ``_soft_subseq_counts``: per order, an exclusive suffix sum along
    each row, then one down each column times the proximity weights. The
    reference shrinks its region by the last row and column after every order,
    but those cells are then already 0 (their suffix sums cover only zero
    cells), so the passes run over the whole padded table: the extra cells add
    signed zeros, which leave every running sum unchanged.
    """
    a, a_len = padded(rows, left, fill=-1)
    b, b_len = padded(rows, right, fill=-1)
    pairs, height, width = a.shape[0], a.shape[1], b.shape[1]
    sizes = np.asarray(sizes, dtype=np.int64)
    counts = np.zeros((pairs, int(sizes.max()) if pairs else 0), dtype=float)
    exponents = np.zeros(counts.shape, dtype=np.int64)
    shift = np.zeros(pairs, dtype=np.int64)
    if exponents_out is not None:
        exponents_out.extend(exponents[p, : sizes[p]] for p in range(pairs))
    started = (a_len > 0) & (b_len > 0) & (sizes > 0)
    if not started.any():
        return [counts[p, : sizes[p]] for p in range(pairs)]
    valid = (a >= 0)[:, :, None] & (b >= 0)[:, None, :]
    weights = np.zeros((pairs, height + 1, width + 1), dtype=float)
    weights[:, :height, :width] = np.where(
        valid,
        np.asarray(prox, dtype=float)[np.maximum(a, 0)[:, :, None], np.maximum(b, 0)[:, None, :]],
        0.0,
    )
    # The reference adds the weights in row-major order.
    total = np.cumsum(weights.reshape(pairs, -1), axis=1)[:, -1]
    counts[started, 0] = total[started] + 1.0
    running = started & (total != 0.0)
    recurrence = weights.copy()
    ids = np.arange(pairs)
    order = 0
    while True:
        order += 1
        # The reference stops at order ksize, or once its region has no row or column.
        shortest = np.minimum(a_len[ids], b_len[ids])
        active = running & (order < sizes[ids]) & (order <= shortest + 1)
        if not active.any():
            break
        if 2 * int(active.sum()) <= ids.size:
            keep = np.flatnonzero(active)
            ids, active = ids[keep], active[keep]
            weights, recurrence = weights[keep], recurrence[keep]
        # The largest region any pair still has; cells past a pair's own region are 0.
        n_rows = min(height + 1, int(a_len[ids].max()) + 2 - order)
        n_columns = min(width + 1, int(b_len[ids].max()) + 2 - order)
        region = recurrence[:, :n_rows, :n_columns]
        region_weights = weights[:, :n_rows, :n_columns]
        # Exclusive suffix sums: cumsum adds in sequence, from the last cell back
        # as the reference does; it starts from x rather than 0.0 + x, which
        # differs only in the sign of a zero, and that sign never reaches a count.
        along_rows = np.cumsum(region[:, :, ::-1], axis=2)
        region[:, :, : n_columns - 1] = along_rows[:, :, : n_columns - 1][:, :, ::-1]
        region[:, :, n_columns - 1] = 0.0
        down_columns = np.cumsum(region[:, ::-1, :], axis=1)
        region[:, : n_rows - 1, :] = (
            down_columns[:, : n_rows - 1, :][:, ::-1, :] * region_weights[:, : n_rows - 1, :]
        )
        region[:, n_rows - 1, :] = 0.0 * region_weights[:, n_rows - 1, :]
        # Columns in turn, each from its last row up.
        column_order = region[:, ::-1, :].transpose(0, 2, 1).reshape(ids.size, -1)
        total = np.cumsum(column_order, axis=1)[:, -1]
        proceed = active & (total != 0.0)
        counts[ids[proceed], order] = total[proceed]
        if exponents_out is not None:
            exponents[ids[proceed], order] = shift[ids[proceed]]
            scale = np.frexp(np.max(np.abs(recurrence), axis=(1, 2)))[1].astype(np.int64)
            np.ldexp(recurrence, -scale[:, None, None], out=recurrence)
            shift[ids] += scale
        running = proceed
    return [counts[p, : sizes[p]] for p in range(pairs)]
