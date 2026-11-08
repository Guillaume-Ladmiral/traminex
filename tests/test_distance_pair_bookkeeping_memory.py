"""The row engine's pair bookkeeping stays bounded (batch-2 review I-2, TX-AUD-20261010-022).

The batched row engine and the per-pair fallback listed every pair as Python
tuples (the batched engine also a pending list and two index lists, the
asymmetric route a member list per pair), so their memory grew as n squared
Python objects: at n = 3000 about 980 MiB of bookkeeping beside a 69 MiB result.
The pairs are now generated in bounded blocks of index arrays. Each test gives
the helper a trivial pair function, so the traced peak is the result matrix
plus the bookkeeping alone.
"""

from __future__ import annotations

import tracemalloc

import numpy as np
import pytest

from traminex.distances import python_kernels as pk

N = 1200  # 719,400 pairs
RESULT_BYTES = N * N * 8
# Blocks of 2**17 pairs hold a few MiB of index and value arrays; listing every
# pair as Python objects takes about 100 MiB at this size.
ALLOWANCE = 16 * 2**20


def _peak_bytes(call) -> tuple[int, np.ndarray]:
    tracemalloc.start()
    try:
        result = call()
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    return peak, result


def _distinct_rows() -> list[np.ndarray]:
    return [np.array([k], dtype=np.int64) for k in range(N)]


def test_per_pair_route_bookkeeping_is_bounded() -> None:
    peak, result = _peak_bytes(lambda: pk._pairwise_apply(N, lambda i, j: float(i + j)))
    assert result[3, 5] == result[5, 3] == 8.0
    assert peak < RESULT_BYTES + ALLOWANCE, peak


def test_stable_unique_route_bookkeeping_is_bounded() -> None:
    rows = _distinct_rows()
    peak, result = _peak_bytes(
        lambda: pk._pairwise_apply_stable_unique(rows, lambda i, j: float(10 * i + j))
    )
    assert result[3, 5] == result[5, 3] == 35.0
    assert peak < RESULT_BYTES + ALLOWANCE, peak


@pytest.mark.parametrize("asymmetric", [False, True], ids=["symmetric", "asymmetric"])
def test_batched_route_bookkeeping_is_bounded(asymmetric: bool) -> None:
    rows = _distinct_rows()

    def raw_batch(left: np.ndarray, right: np.ndarray) -> np.ndarray:
        return (10 * left + right).astype(float)

    def finish(i: int, j: int, raw) -> float:
        return float(raw)

    peak, result = _peak_bytes(
        lambda: pk._pairwise_apply_batched(rows, asymmetric, raw_batch, finish)
    )
    assert result[3, 5] == result[5, 3] == 35.0
    assert peak < RESULT_BYTES + ALLOWANCE, peak
