import numpy as np
import pytest

from traminex import seqdef
from traminex.distances import pairwise as auto_pairwise
from traminex.distances.python_kernels import pairwise as py_pairwise


def test_om_normalize_and_indel_by_pos():
    # Two sequences differing in two positions
    mat = np.array(
        [
            [0, 1, 2, 2],
            [0, 2, 1, 2],
        ],
        dtype=int,
    )
    x = seqdef(mat)
    # Base OM with indel=1, subm=2
    # Two mismatches are optimally handled by two indels (1+1) rather than two substitutions (2+2)
    D = py_pairwise(x.data, method="OM", indel=1.0, subm=2.0)
    assert np.isclose(D[0, 1], 2.0)
    # Normalized should divide by L=4
    Dn = py_pairwise(x.data, method="OM", indel=1.0, subm=2.0, norm="length")
    assert np.isclose(Dn[0, 1], D[0, 1] / 4.0)
    # Position-dependent indel: make later positions cheaper to insert/delete; still same here due to symmetry
    indel_by_pos = np.array([2.0, 2.0, 0.1, 0.1], dtype=float)
    Dpos = py_pairwise(x.data, method="OM", indel=1.0, subm=2.0, indel_by_pos=indel_by_pos)
    assert Dpos[0, 1] >= 0.0
    # Dispatcher normalization (compiles fallback to python here) should match python kernel
    Dn_auto = auto_pairwise(x.data, method="OM", indel=1.0, subm=2.0, norm="length")
    assert np.isclose(Dn_auto[0, 1], Dn[0, 1])


def test_om_max_normalization_uses_position_dependent_substitution_maximum():
    mat = np.array([[0, 0, 0], [1, 1, 1]], dtype=np.int32)
    subm_by_pos = np.array(
        [
            [[0.0, 1.0], [1.0, 0.0]],
            [[0.0, 4.0], [4.0, 0.0]],
            [[0.0, 1.0], [1.0, 0.0]],
        ]
    )
    indel_by_pos = np.full(3, 3.0)

    raw = auto_pairwise(
        mat,
        method="OM",
        indel=3.0,
        indel_by_pos=indel_by_pos,
        subm_by_pos=subm_by_pos,
    )
    normalized = auto_pairwise(
        mat,
        method="OM",
        indel=3.0,
        indel_by_pos=indel_by_pos,
        subm_by_pos=subm_by_pos,
        norm="max",
    )

    # The regular-OM max contract is L * min(max substitution, 2 * max indel).
    np.testing.assert_allclose(normalized, raw / (3.0 * 4.0), rtol=0, atol=0)


def test_dispatch_rejects_legacy_normalize_kwarg():
    mat = np.array(
        [
            [0, 1, 2, 2],
            [0, 2, 1, 2],
        ],
        dtype=int,
    )
    x = seqdef(mat)
    with pytest.raises(TypeError, match="normalize is no longer supported"):
        auto_pairwise(x.data, method="OM", indel=1.0, subm=2.0, normalize=True)
