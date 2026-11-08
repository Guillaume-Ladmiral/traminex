"""NMS and SVRspell with a proximity matrix holding negative entries match R (TX-REP-20261019-001).

The batch-2 review (M-6) found no test passing a ``prox`` with a negative entry. The fixture
``scripts/r_goldens/make_goldens_negative_prox.R`` records R's full matrix and the distances to reference sequence
2; traminex equals both.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx

GOLDEN = pd.read_csv(Path(__file__).resolve().parent / "r_goldens" / "negative_prox.csv", dtype=str)
X = tx.seqdef(np.array([list("AABCC"), list("ABBCA"), list("CCBAA"), list("BBCAB")], dtype=object))
PROX = np.array([[1, -0.4, 0.2], [-0.4, 1, 0.3], [0.2, 0.3, 1]])


def _r(method: str, refseq: str) -> np.ndarray:
    row = GOLDEN[(GOLDEN["method"] == method) & (GOLDEN["refseq"] == refseq)]
    assert len(row) == 1
    return np.array([float(v) for v in row["value"].iloc[0].split()])


@pytest.mark.parametrize("method", ["NMS", "SVRspell"])
def test_negative_proximities_match_r(method: str) -> None:
    assert (PROX < 0).any()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        full = np.asarray(tx.seqdist(X, method=method, prox=PROX), dtype=float)
        ref = np.asarray(tx.seqdist(X, method=method, prox=PROX, refseq=2), dtype=float).ravel()
    np.testing.assert_allclose(full.ravel(order="F"), _r(method, "none"), rtol=0, atol=1e-10)
    np.testing.assert_allclose(ref, _r(method, "2"), rtol=0, atol=1e-10)
