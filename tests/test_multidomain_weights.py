import warnings

import numpy as np
import pytest

import traminex as tx


def create_test_multidomain():
    """Create simple 2-channel test data with clear group structure."""
    # Channel 1: Groups differ primarily here
    ch1 = tx.seqdef(
        np.array(
            [
                [0, 0, 1],  # Group 0
                [0, 0, 1],  # Group 0
                [1, 1, 0],  # Group 1
                [1, 1, 0],  # Group 1
            ]
        )
    )
    # Channel 2: Groups differ less here
    ch2 = tx.seqdef(
        np.array(
            [
                [0, 1, 1],  # Group 0
                [0, 1, 1],  # Group 0
                [1, 0, 0],  # Group 1
                [1, 0, 1],  # Group 1 (slightly different)
            ]
        )
    )
    md = tx.SeqMD([ch1, ch2], names=["primary", "secondary"])
    groups = [0, 0, 1, 1]
    return md, groups


def test_gridsearch_basic():
    """Test basic gridsearch functionality."""
    md, groups = create_test_multidomain()

    # Coarse grid for speed
    weights = tx.seqMD_gridsearch_weights(md, groups, method="HAM", weight_step=0.5)

    assert isinstance(weights, list)
    assert len(weights) == 2
    assert all(isinstance(w, (int, float)) for w in weights)
    assert all(w >= 0 for w in weights)


def test_gridsearch_return_all():
    """Test gridsearch with return_all=True."""
    md, groups = create_test_multidomain()

    result = tx.seqMD_gridsearch_weights(md, groups, method="HAM", weight_step=0.5, return_all=True)

    assert isinstance(result, dict)
    assert "best_weights" in result
    assert "best_score" in result
    assert "all_weights" in result
    assert "all_scores" in result

    assert len(result["best_weights"]) == 2
    assert isinstance(result["best_score"], (int, float))
    assert len(result["all_weights"]) > 0
    assert len(result["all_scores"]) == len(result["all_weights"])


def test_gridsearch_silhouette_range():
    """Test that silhouette scores are in valid range."""
    md, groups = create_test_multidomain()

    result = tx.seqMD_gridsearch_weights(md, groups, method="HAM", weight_step=0.5, return_all=True)

    # Silhouette should be in [-1, 1]
    assert -1 <= result["best_score"] <= 1
    assert all(-1 <= s <= 1 for s in result["all_scores"])


def test_gridsearch_custom_grid():
    """Test gridsearch with custom weight grid."""
    md, groups = create_test_multidomain()

    # Provide explicit grid
    custom_grid = [
        [1.0, 0.0],
        [0.5, 0.5],
        [0.0, 1.0],
    ]

    weights = tx.seqMD_gridsearch_weights(md, groups, method="HAM", weight_grid=custom_grid)

    assert isinstance(weights, list)
    assert len(weights) == 2


def test_gridsearch_empty_weight_grid_raises():
    """Empty explicit weight grids must be rejected."""
    md, groups = create_test_multidomain()

    with pytest.raises(ValueError, match="contain at least one"):
        tx.seqMD_gridsearch_weights(md, groups, method="HAM", weight_grid=[])


@pytest.mark.parametrize(
    "consumer",
    [
        lambda md: tx.seqMD_gridsearch_weights(md, [], method="HAM"),
        lambda md: tx.seqMD_optimize_weights(md, [], method="HAM"),
    ],
)
def test_weight_search_consumers_revalidate_empty_channels_after_mutation(consumer):
    md, _ = create_test_multidomain()
    md.channels.clear()

    with pytest.raises(ValueError, match="channels must contain at least one channel"):
        consumer(md)


@pytest.mark.parametrize(
    "consumer",
    [
        lambda md: tx.seqMD_gridsearch_weights(
            md,
            [0, 1, 0, 1],
            method="HAM",
            weight_grid=[[1.0, 1.0]],
        ),
        lambda md: tx.seqMD_optimize_weights(
            md,
            [0, 1, 0, 1],
            method="HAM",
            max_iter=1,
        ),
    ],
)
def test_weight_search_consumers_revalidate_unaligned_channels_after_mutation(consumer):
    md, _ = create_test_multidomain()
    md.channels[1] = tx.seqdef([["A", "B", "C"]])

    with pytest.raises(ValueError, match="channels must contain the same number of sequences"):
        consumer(md)


def test_combiner_refuses_broadcastable_distance_shapes():
    from traminex.multidomain import _combine_channel_distances

    with pytest.raises(ValueError, match="same square shape"):
        _combine_channel_distances(
            [np.zeros((2, 2)), np.zeros((1, 1))],
            [1.0, 1.0],
        )


def test_gridsearch_invalid_groups():
    """Test that invalid groups raise errors."""
    md, groups = create_test_multidomain()

    # Wrong length
    with pytest.raises(ValueError, match="groups length"):
        tx.seqMD_gridsearch_weights(md, [0, 0], method="HAM")

    # Single group
    with pytest.raises(ValueError, match="at least 2 distinct"):
        tx.seqMD_gridsearch_weights(md, [0, 0, 0, 0], method="HAM")


def test_gridsearch_unsupported_metric():
    """Test that unsupported metrics raise errors."""
    md, groups = create_test_multidomain()

    with pytest.raises(ValueError, match="Unsupported metric"):
        tx.seqMD_gridsearch_weights(md, groups, method="HAM", metric="unknown")


def test_gridsearch_channel_preference():
    """Test that gridsearch can identify more discriminative channel."""
    # Channel 1 perfectly separates groups
    ch1 = tx.seqdef(
        np.array(
            [
                [0, 0, 0],  # Group 0
                [0, 0, 0],  # Group 0
                [1, 1, 1],  # Group 1
                [1, 1, 1],  # Group 1
            ]
        )
    )
    # Channel 2 doesn't separate at all
    ch2 = tx.seqdef(
        np.array(
            [
                [0, 1, 0],  # Group 0
                [1, 0, 1],  # Group 0
                [0, 1, 0],  # Group 1
                [1, 0, 1],  # Group 1
            ]
        )
    )
    md = tx.SeqMD([ch1, ch2])
    groups = [0, 0, 1, 1]

    result = tx.seqMD_gridsearch_weights(md, groups, method="HAM", weight_step=0.5, return_all=True)

    # Channel 1 should get more weight (first component larger)
    assert result["best_weights"][0] > result["best_weights"][1]
    # Silhouette should be high
    assert result["best_score"] > 0.5


def test_optimize_basic():
    """Test basic optimize functionality."""
    md, groups = create_test_multidomain()

    result = tx.seqMD_optimize_weights(md, groups, method="HAM")

    assert isinstance(result, dict)
    assert "weights" in result
    assert "score" in result
    assert "iterations" in result
    assert "converged" in result

    assert len(result["weights"]) == 2
    assert isinstance(result["score"], (int, float))
    assert isinstance(result["iterations"], int)
    assert isinstance(result["converged"], bool)


def test_optimize_convergence():
    """Test that optimize converges for simple problem."""
    md, groups = create_test_multidomain()

    result = tx.seqMD_optimize_weights(md, groups, method="HAM", max_iter=20)

    # Should converge in reasonable iterations
    assert result["iterations"] <= 20
    # If converged, iterations should be less than max
    if result["converged"]:
        assert result["iterations"] < 20


def test_optimize_init_weights():
    """Test optimize with custom initial weights."""
    md, groups = create_test_multidomain()

    init = [1.5, 0.5]
    result = tx.seqMD_optimize_weights(md, groups, method="HAM", init_weights=init, max_iter=10)

    assert isinstance(result["weights"], list)
    assert len(result["weights"]) == 2


def test_optimize_invalid_init_weights():
    """Test that invalid init_weights raises error."""
    md, groups = create_test_multidomain()

    with pytest.raises(ValueError, match="must have"):
        tx.seqMD_optimize_weights(md, groups, init_weights=[1.0])

    with pytest.raises(ValueError, match="positive sum"):
        tx.seqMD_optimize_weights(md, groups, init_weights=[0.0, 0.0])

    with pytest.raises(ValueError, match="non-negative"):
        tx.seqMD_optimize_weights(md, groups, init_weights=[1.0, -1.0])


@pytest.mark.parametrize("max_iter", [0, -1])
def test_optimize_rejects_nonpositive_max_iter(max_iter):
    md, groups = create_test_multidomain()

    with pytest.raises(ValueError, match="max_iter must be a positive integer"):
        tx.seqMD_optimize_weights(md, groups, method="HAM", max_iter=max_iter)


@pytest.mark.parametrize("max_iter", [True, 1.5, "2"])
def test_optimize_rejects_noninteger_max_iter(max_iter):
    md, groups = create_test_multidomain()

    with pytest.raises(TypeError, match="max_iter must be a positive integer"):
        tx.seqMD_optimize_weights(md, groups, method="HAM", max_iter=max_iter)


@pytest.mark.parametrize("tol", [-0.1, np.nan, np.inf, -np.inf])
def test_optimize_rejects_invalid_numeric_tolerance(tol):
    md, groups = create_test_multidomain()

    with pytest.raises(ValueError, match="tol must be finite and non-negative"):
        tx.seqMD_optimize_weights(md, groups, method="HAM", tol=tol)


def test_optimize_rejects_nonnumeric_tolerance():
    md, groups = create_test_multidomain()

    with pytest.raises(TypeError, match="tol must be a real number"):
        tx.seqMD_optimize_weights(md, groups, method="HAM", tol="small")


@pytest.mark.parametrize(
    "init_weights,match",
    [
        ([1.0, np.nan], "finite"),
        ([1.0, np.inf], "finite"),
        ([1.0, "heavy"], "numeric"),
        (np.array([[1.0, 1.0]]), "one-dimensional"),
    ],
)
def test_optimize_rejects_invalid_initial_weight_values(init_weights, match):
    md, groups = create_test_multidomain()

    with pytest.raises((TypeError, ValueError), match=match):
        tx.seqMD_optimize_weights(md, groups, method="HAM", init_weights=init_weights)


def test_optimize_vs_gridsearch():
    """Test that optimize finds similar or better solution than coarse grid."""
    md, groups = create_test_multidomain()

    # Coarse grid search
    grid_result = tx.seqMD_gridsearch_weights(
        md, groups, method="HAM", weight_step=0.5, return_all=True
    )

    # Optimization
    opt_result = tx.seqMD_optimize_weights(md, groups, method="HAM")

    # Optimize should find at least as good a solution
    # (allowing small tolerance for numerical differences)
    assert opt_result["score"] >= grid_result["best_score"] - 0.01


def test_weight_normalization():
    """Test that weights are properly normalized."""
    md, groups = create_test_multidomain()

    result = tx.seqMD_gridsearch_weights(md, groups, method="HAM", weight_step=0.5, return_all=True)

    # Weights should sum to n_channels (2 in this case)
    weight_sum = sum(result["best_weights"])
    assert abs(weight_sum - 2.0) < 0.01


def test_silhouette_computation():
    """Test internal silhouette computation with known case."""
    from traminex.multidomain import _silhouette_from_distance_matrix

    # Simple case: two well-separated groups
    D = np.array(
        [
            [0.0, 0.1, 1.0, 1.0],
            [0.1, 0.0, 1.0, 1.0],
            [1.0, 1.0, 0.0, 0.1],
            [1.0, 1.0, 0.1, 0.0],
        ]
    )
    groups = np.array([0, 0, 1, 1])

    score = _silhouette_from_distance_matrix(D, groups)

    # Should be high (close to 1) for well-separated groups
    assert score > 0.8
    from sklearn.metrics import silhouette_score

    assert score == pytest.approx(silhouette_score(D, groups, metric="precomputed"))


def test_silhouette_poor_separation():
    """Test silhouette with poorly separated groups."""
    from traminex.multidomain import _silhouette_from_distance_matrix

    # Poor separation: mixed distances
    D = np.array(
        [
            [0.0, 0.5, 0.5, 0.5],
            [0.5, 0.0, 0.5, 0.5],
            [0.5, 0.5, 0.0, 0.5],
            [0.5, 0.5, 0.5, 0.0],
        ]
    )
    groups = np.array([0, 0, 1, 1])

    score = _silhouette_from_distance_matrix(D, groups)

    # Should be low or negative for poor separation
    assert score < 0.5


def test_three_channel_gridsearch():
    """Test gridsearch with 3 channels."""
    ch1 = tx.seqdef(np.array([[0, 1], [0, 1], [1, 0], [1, 0]]))
    ch2 = tx.seqdef(np.array([[0, 0], [0, 0], [1, 1], [1, 1]]))
    ch3 = tx.seqdef(np.array([[1, 1], [1, 1], [0, 0], [0, 0]]))
    md = tx.SeqMD([ch1, ch2, ch3])
    groups = [0, 0, 1, 1]

    # Very coarse grid for 3 channels
    weights = tx.seqMD_gridsearch_weights(md, groups, method="HAM", weight_step=0.5)

    assert len(weights) == 3
    assert all(w >= 0 for w in weights)


def test_three_channel_optimize():
    """Test optimization with 3 channels."""
    ch1 = tx.seqdef(np.array([[0, 1], [0, 1], [1, 0], [1, 0]]))
    ch2 = tx.seqdef(np.array([[0, 0], [0, 0], [1, 1], [1, 1]]))
    ch3 = tx.seqdef(np.array([[1, 1], [1, 1], [0, 0], [0, 0]]))
    md = tx.SeqMD([ch1, ch2, ch3])
    groups = [0, 0, 1, 1]

    result = tx.seqMD_optimize_weights(md, groups, method="HAM", max_iter=10)

    assert len(result["weights"]) == 3
    assert result["score"] > -1  # At least some separation


def test_per_channel_kwargs():
    """Test that per_channel_kwargs are passed through."""
    md, groups = create_test_multidomain()

    # Different indel costs per channel
    per_ch_kwargs = [
        {"indel": 1.0, "subm": 2.0},
        {"indel": 2.0, "subm": 2.0},
    ]

    # Should not raise error
    weights = tx.seqMD_gridsearch_weights(
        md, groups, method="OM", weight_step=0.5, per_channel_kwargs=per_ch_kwargs
    )

    assert len(weights) == 2


def test_singleton_cluster_silhouette():
    """Test silhouette computation with singleton clusters."""
    from traminex.multidomain import _silhouette_from_distance_matrix

    D = np.array(
        [
            [0.0, 1.0, 1.0],
            [1.0, 0.0, 0.1],
            [1.0, 0.1, 0.0],
        ]
    )
    groups = np.array([0, 1, 1])  # Group 0 is singleton

    score = _silhouette_from_distance_matrix(D, groups)

    from sklearn.metrics import silhouette_samples

    independent_samples = silhouette_samples(D, groups, metric="precomputed")
    np.testing.assert_allclose(independent_samples, [0.0, 0.9, 0.9])
    assert score == pytest.approx(0.6)
    assert score == pytest.approx(float(independent_samples.mean()))


def test_weight_grid_generation():
    """Test internal weight grid generation."""
    from traminex.multidomain import _generate_weight_grid

    grid = _generate_weight_grid(2, 0.5)

    assert len(grid) > 0
    assert all(len(w) == 2 for w in grid)
    # All should sum to approximately 1 (normalized)
    assert all(abs(sum(w) - 1.0) < 0.01 for w in grid)


def test_gridsearch_with_string_groups():
    """Test that string group labels work."""
    md, _ = create_test_multidomain()
    groups = ["A", "A", "B", "B"]

    weights = tx.seqMD_gridsearch_weights(md, groups, method="HAM", weight_step=0.5)

    assert len(weights) == 2


@pytest.mark.parametrize("weight_step", [0.0, -0.1, np.nan, np.inf, -np.inf])
def test_gridsearch_rejects_invalid_numeric_weight_step(weight_step):
    md, groups = create_test_multidomain()

    with pytest.raises(ValueError, match="weight_step must be finite and strictly positive"):
        tx.seqMD_gridsearch_weights(md, groups, method="HAM", weight_step=weight_step)


@pytest.mark.parametrize("weight_step", [True, "coarse"])
def test_gridsearch_rejects_nonnumeric_weight_step(weight_step):
    md, groups = create_test_multidomain()

    with pytest.raises(TypeError, match="weight_step must be a real number"):
        tx.seqMD_gridsearch_weights(md, groups, method="HAM", weight_step=weight_step)


@pytest.mark.parametrize(
    "weight_grid,match",
    [
        ([[1.0, np.nan]], "finite"),
        ([[1.0, np.inf]], "finite"),
        ([[1.0, "heavy"]], "numeric"),
        ([np.array([[1.0, 1.0]])], "one-dimensional"),
        ([[0.0, 0.0]], "positive sum"),
        ([[-1.0, 2.0]], "non-negative"),
    ],
)
def test_gridsearch_rejects_invalid_explicit_weight_vectors(weight_grid, match):
    md, groups = create_test_multidomain()

    with pytest.raises((TypeError, ValueError), match=match):
        tx.seqMD_gridsearch_weights(md, groups, method="HAM", weight_grid=weight_grid)


def test_additive_optimizers_preserve_inputs_and_exact_result_schemas():
    md, groups = create_test_multidomain()
    grid = [[1.0, 0.0], [0.5, 0.5], [0.0, 1.0]]
    initial = [1.5, 0.5]
    grid_before = [weights.copy() for weights in grid]
    initial_before = initial.copy()

    grid_result = tx.seqMD_gridsearch_weights(
        md, groups, method="HAM", weight_grid=grid, return_all=True
    )
    optimize_result = tx.seqMD_optimize_weights(
        md,
        groups,
        method="HAM",
        init_weights=initial,
        max_iter=1,
        tol=0.0,
    )

    assert set(grid_result) == {"best_weights", "best_score", "all_weights", "all_scores"}
    assert set(optimize_result) == {"weights", "score", "iterations", "converged"}
    assert sum(grid_result["best_weights"]) == pytest.approx(2.0)
    assert sum(optimize_result["weights"]) > 0.0
    assert grid == grid_before
    assert initial == initial_before


@pytest.mark.parametrize(
    "grid_kwargs",
    [{"weight_step": 0.5}, {"weight_grid": [[1.0, 0.0], [0.5, 0.5], [0.0, 1.0]]}],
    ids=["generated", "explicit"],
)
def test_gridsearch_reports_all_weights_on_the_scale_it_tested(grid_kwargs):
    """Audit 2026-09-22 S2-16: ``best_weights`` was not among ``all_weights``.

    Every grid vector is rescaled to sum to the number of channels before it is
    tested, and ``best_weights`` is reported on that scale, but ``all_weights``
    listed the grid as given (summing to 1). The documented contract is "all
    tested weight combinations", so they are now reported as tested.
    """
    md, groups = create_test_multidomain()

    result = tx.seqMD_gridsearch_weights(md, groups, method="HAM", return_all=True, **grid_kwargs)

    assert all(sum(weights) == pytest.approx(2.0) for weights in result["all_weights"])
    best = result["all_scores"].index(result["best_score"])
    assert result["all_weights"][best] == result["best_weights"]


@pytest.mark.parametrize("scale", [1.0, 1e308], ids=["ordinary", "large-finite"])
def test_equal_finite_grid_weights_preserve_normalization_and_silhouette(scale):
    md, groups = create_test_multidomain()
    with warnings.catch_warnings():
        # Check the numerical result independently of the validator's warning.
        warnings.simplefilter("ignore", RuntimeWarning)
        result = tx.seqMD_gridsearch_weights(
            md, groups, method="HAM", weight_grid=[[scale, scale]], return_all=True
        )

    np.testing.assert_allclose(result["best_weights"], [1.0, 1.0], rtol=0, atol=1e-15)
    np.testing.assert_allclose(result["all_weights"], [[1.0, 1.0]], rtol=0, atol=1e-15)
    # Equal-weight HAM: between-group distances are 6 and 5, within-group
    # distances are 0 and 1. Silhouettes are [1, 1, 5/6, 4/5].
    assert result["best_score"] == pytest.approx(109 / 120)
    assert result["all_scores"] == pytest.approx([109 / 120])


@pytest.mark.parametrize("scale", [1.0, 1e308], ids=["ordinary", "large-finite"])
def test_coordinate_search_has_scale_equivalent_starting_weights(scale):
    md, groups = create_test_multidomain()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        result = tx.seqMD_optimize_weights(
            md, groups, method="HAM", init_weights=[scale, scale], max_iter=1, tol=0
        )
    # The first coordinate grows to 4, then the weaker channel falls to 0.
    # Only starting weights are normalized; coordinate updates retain their
    # documented scale, so the optimized vector need not sum to two.
    np.testing.assert_allclose(result["weights"], [4.0, 0.0], rtol=0, atol=1e-15)
    assert result["score"] == pytest.approx(1.0)
    assert result["iterations"] == 1
    assert result["converged"] is False


@pytest.mark.parametrize("consumer", ["grid", "optimize"])
def test_finite_extreme_search_weights_work_with_runtime_warnings_as_errors(consumer):
    md, groups = create_test_multidomain()
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        if consumer == "grid":
            result = tx.seqMD_gridsearch_weights(
                md, groups, method="HAM", weight_grid=[[1e308, 1e308]], return_all=True
            )
            assert result["best_weights"] == pytest.approx([1.0, 1.0])
        else:
            result = tx.seqMD_optimize_weights(
                md, groups, method="HAM", init_weights=[1e308, 1e308], max_iter=1
            )
            assert result["weights"] == pytest.approx([4.0, 0.0])


def test_raw_additive_distances_allow_finite_weights_with_unrepresentable_total():
    # Each weighted channel distance is exactly zero, so every output remains
    # representable even though the sum of the unused weight magnitudes does not.
    channel = tx.seqdef([["A", "B"], ["A", "B"]])
    md = tx.SeqMD([channel, channel.copy()])
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        distances = tx.seqMD_distances(md, method="HAM", weights=[1e308, 1e308])
    np.testing.assert_array_equal(distances, [[0.0, 0.0], [0.0, 0.0]])


def test_raw_additive_distances_keep_supplied_weight_scale():
    md, _ = create_test_multidomain()
    distances = tx.seqMD_distances(md, method="HAM", weights=[2.0, 3.0])
    np.testing.assert_array_equal(
        distances, [[0, 0, 15, 12], [0, 0, 15, 12], [15, 15, 0, 3], [12, 12, 3, 0]]
    )
