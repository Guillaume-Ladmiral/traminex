"""The void route's docstring names the engine it calls (TX-AUD-20261019-032).

``_pairwise_void_delete`` said it routed ragged input "to the per-pair Python kernels"; it calls the row engine,
``python_kernels._pairwise_rows``, which evaluates pairs in vectorized batches where the method has a batched kernel.
"""

from __future__ import annotations

import numpy as np

import traminex as tx
from traminex import distances
from traminex.distances import python_kernels


def test_ragged_input_reaches_the_row_engine(monkeypatch) -> None:
    calls: list[str] = []
    original = python_kernels._pairwise_rows

    def spy(rows, method, *args, **kwargs):
        calls.append(method)
        return original(rows, method, *args, **kwargs)

    monkeypatch.setattr(python_kernels, "_pairwise_rows", spy)
    rows = np.array([["A", "B", None], ["A", "A", "B"], ["B", None, None]], dtype=object)
    D = np.asarray(tx.seqdist(tx.seqdef(rows, alphabet=["A", "B"]), method="LCS"), dtype=float)
    assert calls == ["LCS"]
    assert D[0, 1] == 1.0  # LCS 'AB' vs 'AAB' = 2 + 3 - 2*2


def test_the_void_route_docstring_names_the_row_engine() -> None:
    doc = distances._pairwise_void_delete.__doc__
    assert "``python_kernels._pairwise_rows``" in doc
    assert "per-pair Python kernels" not in doc
