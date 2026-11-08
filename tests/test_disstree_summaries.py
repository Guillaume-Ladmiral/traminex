import numpy as np

from traminex.disstree import disstree_feature_importance, disstree_fit, disstree_leaf_summary


def test_disstree_summaries():
    D = np.array(
        [
            [0, 1, 1, 5],
            [1, 0, 1, 5],
            [1, 1, 0, 5],
            [5, 5, 5, 0],
        ],
        dtype=float,
    )
    y = np.array([0, 0, 0, 1])
    res = disstree_fit(D, y, n_components=2, max_depth=2, random_state=0)
    df = disstree_leaf_summary(res)
    assert set(df.columns) == {"leaf", "size"}
    imp = disstree_feature_importance(res)
    assert set(imp.columns) == {"feature", "importance"}
    assert imp.shape[0] == res.coords.shape[1]
