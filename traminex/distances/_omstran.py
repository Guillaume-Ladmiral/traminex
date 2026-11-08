"""OMstran transition-augmented OM transform (TraMineR seqdist-OMstran.R).

Each position is replaced by a transition token ``(from, to)`` (or ``(prev, from, to)`` when
``previous=True``); a new substitution matrix and per-token indel vector are derived from the
origin-state costs weighted by ``otto`` plus a transition penalty weighted by ``1 - otto``.
The transformed sequences are then compared with ordinary OM (per-state indels).
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def _transition_rates_rows(
    rows: list[np.ndarray], K: int, weights: np.ndarray | None
) -> np.ndarray:
    counts = np.zeros((K, K), dtype=float)
    if weights is None:
        row_weights = np.ones(len(rows), dtype=float)
    else:
        row_weights = np.asarray(weights, dtype=float)
        if row_weights.ndim != 1 or row_weights.shape[0] != len(rows):
            raise ValueError("weights length must equal number of sequences")
    for row, weight in zip(rows, row_weights):
        if len(row) > 1:
            np.add.at(counts, (row[:-1], row[1:]), float(weight))
    rs = counts.sum(axis=1, keepdims=True)
    return np.divide(counts, rs, out=np.zeros_like(counts), where=rs > 0)


def _transition_rates(seqmat: np.ndarray, K: int, weights: np.ndarray | None) -> np.ndarray:
    """Dense compatibility wrapper around the canonical-row implementation."""
    return _transition_rates_rows(
        [np.asarray(row, dtype=np.int64) for row in np.asarray(seqmat)], K, weights
    )


def omstran_transform_rows(
    input_rows: list[np.ndarray],
    *,
    otto: float,
    transindel: str = "constant",
    previous: bool = False,
    add_column: bool = True,
    indel: float | np.ndarray = 1.0,
    subm: float | np.ndarray = 2.0,
    weights: np.ndarray | None = None,
    n_states: int | None = None,
) -> tuple[list[np.ndarray], np.ndarray, np.ndarray]:
    """Row-native OMstran transform for canonical, void-deleted sequences."""
    if transindel not in ("constant", "prob", "subcost"):
        raise ValueError("transindel must be one of: 'constant', 'prob', 'subcost'")

    rows: list[np.ndarray] = []
    observed_states = 0
    for raw in input_rows:
        row = np.asarray(raw)
        if row.ndim != 1 or row.dtype.kind not in {"i", "u"}:
            raise TypeError("OMstran rows must be one-dimensional integer arrays")
        row = row.astype(np.int64, copy=False)
        if row.size:
            if int(row.min()) < 0:
                raise ValueError("OMstran state codes cannot be negative")
            observed_states = max(observed_states, int(row.max()) + 1)
        rows.append(row)

    declared_widths = [observed_states]
    if n_states is not None:
        declared_widths.append(int(n_states))
    if not np.isscalar(indel):
        declared_widths.append(len(np.asarray(indel)))
    if not np.isscalar(subm):
        matrix = np.asarray(subm)
        if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
            raise ValueError("OMstran subm must be a square matrix")
        declared_widths.append(int(matrix.shape[0]))
    K = max(declared_widths, default=0)
    if n_states is not None and int(n_states) < observed_states:
        raise ValueError("n_states cannot be smaller than the observed state width")

    indel_vec: NDArray[np.float64] = (
        np.full(K, float(np.asarray(indel).item()), dtype=float)
        if np.isscalar(indel)
        else np.asarray(indel, dtype=float)
    )
    if indel_vec.ndim != 1 or len(indel_vec) != K:
        raise ValueError("OMstran indel vector length must equal the active state count")
    if not np.all(np.isfinite(indel_vec)) or np.any(indel_vec <= 0):
        raise ValueError("OMstran indel cost must be finite and strictly positive")
    if np.isscalar(subm):
        sm: NDArray[np.float64] = np.full((K, K), float(np.asarray(subm).item()), dtype=float)
        np.fill_diagonal(sm, 0.0)
    else:
        sm = np.asarray(subm, dtype=float)
        if sm.shape != (K, K):
            raise ValueError("OMstran subm shape must match the active state count")
    if transindel == "prob" and not any(len(row) > 1 for row in rows):
        # Transition rates need a transition: seqtrate refuses such data, and traminex
        # refuses TRATE costs on them (TX-AUD-20261019-080; R stops here too).
        raise ValueError(
            "OMstran transindel='prob' needs transition rates: absolute lag must be smaller "
            "than the sequence width"
        )
    tr = _transition_rates_rows(rows, K, weights) if transindel == "prob" else None

    token_rows: list[list[tuple[int, ...]]] = []
    for row in rows:
        length = len(row)
        tokens: list[tuple[int, ...]] = []
        if add_column:
            for position in range(length):
                current = int(row[position])
                following = int(row[min(position + 1, length - 1)])
                if previous:
                    preceding = int(row[max(position - 1, 0)])
                    tokens.append((preceding, current, following))
                else:
                    tokens.append((current, following))
        elif previous:
            for position in range(1, length - 1):
                tokens.append(
                    (
                        int(row[position - 1]),
                        int(row[position]),
                        int(row[position + 1]),
                    )
                )
        else:
            for position in range(max(0, length - 1)):
                tokens.append((int(row[position]), int(row[position + 1])))
        token_rows.append(tokens)

    alphabet = sorted({token for tokens in token_rows for token in tokens})
    code = {token: index for index, token in enumerate(alphabet)}
    transformed = [
        np.asarray([code[token] for token in tokens], dtype=np.int64) for tokens in token_rows
    ]

    max_indel = float(indel_vec.max()) if indel_vec.size else 1.0
    max_sm = float(sm.max()) if sm.size else 1.0
    indelrate = max_indel / max_sm if max_sm > 0 else 0.0
    transweight = 1.0 - otto
    if previous:
        transweight /= 2.0
    transweight *= indelrate
    indel_norm = indelrate * indel_vec / max_indel if max_indel > 0 else indel_vec
    sm_norm = sm / max_sm if max_sm > 0 else sm
    stateindel = indel_norm * otto
    compare = 1 if previous else 0

    def raw_transindel(token: tuple[int, ...]) -> float:
        if previous:
            preceding, current, following = token
            if transindel == "constant":
                return float((preceding != current) + (current != following))
            if transindel == "prob":
                assert tr is not None
                return float(2.0 - tr[preceding, current] - tr[current, following])
            return float(sm_norm[preceding, current] + sm_norm[current, following])
        current, following = token
        if current == following:
            return 0.0
        if transindel == "constant":
            return 1.0
        if transindel == "prob":
            assert tr is not None
            return float(1.0 - tr[current, following])
        return float(sm_norm[current, following])

    indels = np.zeros(len(alphabet), dtype=float)
    for token, index in code.items():
        current = token[compare]
        indels[index] = stateindel[current] + transweight * raw_transindel(token)

    new_sm = np.zeros((len(alphabet), len(alphabet)), dtype=float)
    for left in range(len(alphabet)):
        left_state = alphabet[left][compare]
        for right in range(left + 1, len(alphabet)):
            right_state = alphabet[right][compare]
            cost = otto * sm_norm[left_state, right_state] + (
                indels[left] + indels[right] - stateindel[left_state] - stateindel[right_state]
            )
            new_sm[left, right] = new_sm[right, left] = cost
    return transformed, indels, new_sm


def omstran_transform(
    seqmat: np.ndarray,
    *,
    otto: float,
    transindel: str = "constant",
    previous: bool = False,
    add_column: bool = True,
    indel: float | np.ndarray = 1.0,
    subm: float | np.ndarray = 2.0,
    weights: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``(new_seqmat, indels_vec, new_sm)`` for OMstran (see seqdist-OMstran.R)."""
    seqmat = np.asarray(seqmat)
    if seqmat.ndim != 2:
        raise ValueError("seqmat must be a two-dimensional sequence matrix")
    transformed, indels, new_sm = omstran_transform_rows(
        [row for row in seqmat],
        otto=otto,
        transindel=transindel,
        previous=previous,
        add_column=add_column,
        indel=indel,
        subm=subm,
        weights=weights,
    )
    width = len(transformed[0]) if transformed else 0
    if any(len(row) != width for row in transformed):
        raise AssertionError("dense OMstran transform produced unequal row lengths")
    dense = np.vstack(transformed) if transformed else np.empty((0, width), dtype=np.int64)
    return dense, indels, new_sm
