import numpy as np
import pytest

from traminex import seqdef


def test_seqtree_smoke():
    pytest.importorskip("sklearn")
    # small labeled example: first two vs last two
    mat = np.array(
        [
            [0, 1, 2, 2],
            [0, 1, 1, 2],
            [2, 2, 1, 0],
            [2, 2, 2, 0],
        ],
        dtype=int,
    )
    x = seqdef(mat)
    y = np.array([0, 0, 1, 1])
    from traminex import seqtree_fit, seqtree_leaf_summary

    res = seqtree_fit(x, y, max_depth=2)
    summary = seqtree_leaf_summary(res)
    assert summary.shape[0] >= 1

    with pytest.warns(DeprecationWarning, match="pass only result"):
        legacy_summary = seqtree_leaf_summary(x, res)
    assert legacy_summary.equals(summary)


@pytest.mark.parametrize("dtype", [np.int8, np.uint8, np.int32, np.uint64])
def test_exploratory_transition_features_keep_narrow_state_pairs(dtype):
    from traminex import STSList
    from traminex.seqtree_explore import _features_from_sequences

    data = np.array([[19, 19], [19, 0]], dtype=dtype)
    data.setflags(write=False)
    source = STSList(data, [f"S{i}" for i in range(20)])
    values, names = _features_from_sequences(source, kind="transitions")
    expected = np.zeros((2, 400))
    expected[0, 399] = 1
    expected[1, 380] = 1
    np.testing.assert_array_equal(values, expected)
    assert names[399] == "tr_S19_S19" and names[380] == "tr_S19_S0"
    np.testing.assert_array_equal(source.data, data)
