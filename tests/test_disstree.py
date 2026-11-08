import numpy as np

from traminex.disstree import disstree_fit, disstree_predict


def test_disstree_fit_and_predict():
    # Two clusters in distance space
    D = np.array(
        [
            [0, 1, 1, 5, 5, 5],
            [1, 0, 1, 5, 5, 5],
            [1, 1, 0, 5, 5, 5],
            [5, 5, 5, 0, 1, 1],
            [5, 5, 5, 1, 0, 1],
            [5, 5, 5, 1, 1, 0],
        ],
        dtype=float,
    )
    y = np.array([0, 0, 0, 1, 1, 1])
    res = disstree_fit(D, y, n_components=2, max_depth=2, random_state=0)
    # Predict for two new points closer to cluster 0 and 1 respectively
    D_new = np.array(
        [
            [0.5, 0.5, 0.5, 6, 6, 6],
            [6, 6, 6, 0.5, 0.5, 0.5],
        ]
    )
    pred = disstree_predict(res, D_new)
    assert pred.shape == (2,)
    # Points are clearly closer to each cluster; projection should preserve this separation
    assert np.array_equal(pred, np.array([0, 1]))


def test_disstree_predict_matches_training_rows_when_distances_match():
    D = np.array(
        [
            [0, 1, 1, 5, 5, 5],
            [1, 0, 1, 5, 5, 5],
            [1, 1, 0, 5, 5, 5],
            [5, 5, 5, 0, 1, 1],
            [5, 5, 5, 1, 0, 1],
            [5, 5, 5, 1, 1, 0],
        ],
        dtype=float,
    )
    y = np.array([0, 0, 0, 1, 1, 1])
    res = disstree_fit(D, y, n_components=2, max_depth=2, random_state=0)

    train_pred = res.model.predict(res.coords)
    idx = np.array([1, 4], dtype=int)
    pred_from_dnew = disstree_predict(res, D[idx, :])

    assert np.array_equal(pred_from_dnew, train_pred[idx])
