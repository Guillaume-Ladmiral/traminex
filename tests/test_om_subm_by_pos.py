import numpy as np
import pandas as pd
import pytest

import traminex as tx
import traminex.distances as distances
from traminex.distances import pairwise as pw


def test_om_with_subm_by_pos_matches_constant_matrix_when_slices_identical():
    df = pd.DataFrame(
        {
            "t1": ["A", "B", "A"],
            "t2": ["B", "B", "A"],
            "t3": ["A", "A", "B"],
        }
    )
    X = tx.seqdef(df)
    K = len(X.alphabet)
    L = X.seqlen
    # constant subm
    M = np.full((K, K), 2.0, dtype=float)
    np.fill_diagonal(M, 0.0)
    # subm_by_pos: identical slices
    S = np.repeat(M[None, :, :], L, axis=0)
    D1 = pw(X.data, method="OM", indel=1.0, subm=M)
    D2 = pw(X.data, method="OM", indel=1.0, subm_by_pos=S)
    assert np.allclose(D1, D2)


def test_compiled_dispatch_does_not_ignore_subm_by_pos(monkeypatch):
    sequences = np.array([[0, 0, 1], [0, 1, 1]], dtype=int)
    costs = np.array(
        [
            [[0.0, 1.0], [1.0, 0.0]],
            [[0.0, 4.0], [4.0, 0.0]],
            [[0.0, 1.0], [1.0, 0.0]],
        ]
    )
    expected = distances._py.pairwise(
        sequences,
        method="OM",
        indel=3.0,
        subm=2.0,
        subm_by_pos=costs,
    )

    class RefusingCompiledStub:
        def pairwise_om(self, *args):
            pytest.fail("position-dependent OM must not use the scalar compiled kernel")

        def pairwise_om_mat(self, *args):
            pytest.fail("position-dependent OM must not use the matrix compiled kernel")

    monkeypatch.setattr(distances, "_pb", RefusingCompiledStub())
    monkeypatch.setattr(distances, "_impl_name", lambda: "pybind11")

    actual = pw(
        sequences,
        method="OM",
        indel=3.0,
        subm=2.0,
        subm_by_pos=costs,
    )

    np.testing.assert_array_equal(actual, expected)
