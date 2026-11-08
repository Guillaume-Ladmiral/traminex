import numpy as np
import pytest

from traminex import seqdef
from traminex.seqtree import seqtree_cv_fit


def test_seqtree_cv_smoke():
    try:
        import sklearn  # noqa: F401
    except Exception:
        pytest.skip("sklearn not available")
    mat = np.array(
        [
            [0, 1, 2, 2],
            [0, 1, 1, 2],
            [2, 2, 1, 0],
            [2, 2, 2, 0],
            [0, 1, 2, 1],
            [2, 1, 0, 0],
        ],
        dtype=int,
    )
    x = seqdef(mat)
    y = np.array([0, 0, 1, 1, 0, 1])
    best, names, cvres = seqtree_cv_fit(x, y, cv=2)
    assert best is not None
