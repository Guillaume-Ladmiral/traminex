"""Structure / seeded test for seqgen (T17 VERIFY).

seqgen draws each cell iid from the alphabet (R: sample(alphabet, length,
replace=TRUE, p)). It is *stochastic* and R's RNG stream is not reproducible
from Python, so there is no numeric R golden -- we pin the structural contract
(shape, alphabet membership, marginal proportions, seeded reproducibility),
matching the deterministic guarantees R itself makes.
"""

from __future__ import annotations

import numpy as np
import pytest

from traminex import seqgen


def test_seqgen_shape_and_alphabet_membership():
    out = seqgen(10, 5, ["A", "B", "C"], seed=42)
    assert out.shape == (10, 5)
    assert set(np.unique(out)).issubset({"A", "B", "C"})


def test_seqgen_seeded_reproducibility():
    a = seqgen(8, 6, ["A", "B", "C", "D"], seed=123)
    b = seqgen(8, 6, ["A", "B", "C", "D"], seed=123)
    np.testing.assert_array_equal(a, b)


def test_seqgen_uniform_marginals_large_n():
    out = seqgen(2000, 10, ["A", "B", "C", "D"], seed=7)
    _, counts = np.unique(out, return_counts=True)
    props = counts / out.size
    np.testing.assert_allclose(props, 0.25, atol=0.02)


def test_seqgen_weighted_marginals_large_n():
    p = [0.7, 0.2, 0.1]
    out = seqgen(3000, 10, ["A", "B", "C"], p=p, seed=11)
    props = {s: (out == s).mean() for s in ["A", "B", "C"]}
    np.testing.assert_allclose([props["A"], props["B"], props["C"]], p, atol=0.02)


def test_seqgen_rejects_bad_probabilities():
    with pytest.raises(ValueError):
        seqgen(5, 3, ["A", "B", "C"], p=[0.5, 0.5])  # wrong length
    # Weights that do not sum to 1 are normalized, as R's sample() reads them
    # (TX-AUD-20261011-072).
    weights = np.array([0.5, 0.4, 0.4])
    np.testing.assert_array_equal(
        seqgen(5, 3, ["A", "B", "C"], p=weights, seed=3),
        seqgen(5, 3, ["A", "B", "C"], p=weights / weights.sum(), seed=3),
    )
