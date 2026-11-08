import numpy as np

from traminex.distances.costs import trate_substitution_matrix
from traminex.distances.python_kernels import pairwise as py_pairwise


def test_trate_modes_and_normalize():
    """`normalize="mean_offdiag"` scales the whole off-diagonal by 1/mean(off-diag).

    Load-bearing: asserts the scaling holds for *every* off-diagonal entry (not
    just one cell) and unconditionally (the fixture is a real TRATE matrix whose
    off-diagonal mean is strictly positive), and that the normalized off-diagonal
    mean is exactly 1. (Explicit smoothing keeps this independent of the
    R-faithful default; R parity of the raw matrix is covered by
    test_r_goldens_costs_trate.py.)
    """
    mat = np.array(
        [
            [0, 1, 1, 2],
            [1, 0, 1, 2],
            [1, 1, 2, 2],
        ],
        dtype=int,
    )
    C = trate_substitution_matrix(mat, smoothing=1e-6, mode="symmetric", normalize="none")
    Cn = trate_substitution_matrix(mat, smoothing=1e-6, mode="symmetric", normalize="mean_offdiag")
    assert C.shape == (3, 3)
    off_mask = ~np.eye(3, dtype=bool)
    off = C[off_mask]
    assert off.mean() > 0  # real TRATE costs -> the scaling is always exercised
    # Every off-diagonal entry is scaled by exactly 1 / mean(off-diagonal).
    np.testing.assert_allclose(Cn[off_mask], C[off_mask] / off.mean(), rtol=0, atol=1e-12)
    # By construction the normalized off-diagonal mean is 1 (to float precision).
    np.testing.assert_allclose(Cn[off_mask].mean(), 1.0, rtol=0, atol=1e-12)
    # Diagonal stays zero under normalization.
    np.testing.assert_array_equal(np.diag(Cn), np.zeros(3))


def test_dhd_variants_relation():
    """DHD `sum` variant is *exactly* twice the `mean` variant for all pairs.

    Per-position cost is 2 - (p_x + p_y) for "sum" vs 1 - 0.5(p_x + p_y) for
    "mean" (see DHDConfig), i.e. exactly 2x at every position; DHD sums those
    costs, so the whole distance matrix is exactly doubled. Asserted across the
    entire matrix at machine precision (was a loose 1.8<=ratio<=2.2 on one cell).
    """
    mat = np.array(
        [
            [0, 1, 2, 2, 1],
            [0, 1, 1, 2, 2],
            [0, 2, 2, 2, 2],
        ],
        dtype=int,
    )
    D_mean = py_pairwise(mat, method="DHD", variant="mean")
    D_sum = py_pairwise(mat, method="DHD", variant="sum")
    assert D_mean.max() > 0  # the relation is non-trivial (distances are not all zero)
    np.testing.assert_allclose(D_sum, 2.0 * D_mean, rtol=0, atol=1e-12)
