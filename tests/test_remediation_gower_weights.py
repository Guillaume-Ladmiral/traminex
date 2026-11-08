"""OBS-002: validate raw Gower weights against its bounded-average contract.

Finite nonnegative feature weights define weighted dissimilarities in [0, 1].
The supplied vector must contain one weight per feature and positive total
weight. Individual zero weights remain useful. Raw Gower's existing zero
fallback for pairs with no comparable positive weight, and its dropping of
constant numeric columns, are distinct from the public FEATURES/daisy path.
"""

import numpy as np
import pytest

import traminex as tx
from traminex.distances.costs import gower_distance

_FEATURES = np.array([[0.0, 0.0], [1.0, 1.0], [0.0, 2.0]])


@pytest.mark.parametrize(
    "weights",
    [
        pytest.param([1.0, -0.9], id="retained-negative-out-of-range-case"),
        pytest.param([-1.0, 2.0], id="negative-first-weight"),
        pytest.param([1.0, np.nan], id="nan"),
        pytest.param([np.inf, 1.0], id="positive-infinity"),
        pytest.param([1.0, -np.inf], id="negative-infinity"),
        pytest.param([0.0, 0.0], id="zero-total"),
    ],
)
def test_raw_gower_rejects_weights_that_do_not_define_bounded_averages(weights):
    with pytest.raises(ValueError):
        gower_distance(_FEATURES, weights=weights)


@pytest.mark.parametrize(
    "weights",
    [
        pytest.param(1.0, id="scalar"),
        pytest.param([[1.0], [2.0]], id="column-matrix"),
        pytest.param([[1.0, 2.0], [3.0, 4.0]], id="square-matrix"),
        pytest.param([1.0], id="too-short"),
        pytest.param([1.0, 2.0, 3.0], id="too-long"),
    ],
)
def test_raw_gower_requires_one_dimensional_feature_weight_vector(weights):
    with pytest.raises(ValueError):
        gower_distance(_FEATURES, weights=weights)


def test_raw_gower_nonnegative_fractional_weights_form_a_bounded_mean():
    # Numeric distances have ranges 1 and 2. The pairwise weighted means are
    # (.25*1+.75*.5), (.25*0+.75*1), and (.25*1+.75*.5).
    actual = gower_distance(_FEATURES, weights=[0.25, 0.75])
    expected = [[0.0, 0.625, 0.75], [0.625, 0.0, 0.625], [0.75, 0.625, 0.0]]
    np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-14)
    assert np.isfinite(actual).all()
    assert np.all((actual >= 0) & (actual <= 1))


def test_raw_gower_individual_zero_weight_excludes_feature_from_denominator():
    actual = gower_distance(_FEATURES, weights=[1.0, 0.0])
    np.testing.assert_array_equal(actual, [[0, 1, 0], [1, 0, 1], [0, 1, 0]])


def test_raw_gower_missing_values_use_pairwise_comparable_weight():
    features = [[0.0, np.nan], [1.0, 0.0], [np.nan, 2.0]]
    actual = gower_distance(features, weights=[1.0, 3.0])
    # Pair 0/1 compares only the first feature; 1/2 only the second. Pair 0/2
    # has no comparable features and preserves the raw extension's zero.
    np.testing.assert_array_equal(actual, [[0, 1, 0], [1, 0, 1], [0, 1, 0]])


def test_raw_gower_zero_weight_and_missing_pair_preserve_zero_fallback():
    features = [[np.nan, 0.0], [1.0, 1.0], [0.0, 2.0]]
    actual = gower_distance(features, weights=[1.0, 0.0])
    # Missing first-feature comparisons cannot borrow the second feature's
    # zero weight. The complete pair still differs maximally on feature one.
    np.testing.assert_array_equal(actual, [[0, 0, 0], [0, 0, 1], [0, 1, 0]])


def test_raw_gower_asymmetric_absence_uses_pair_specific_weight_denominator():
    features = [[0, 0], [0, 1], [1, 0], [1, 1]]
    actual = gower_distance(features, weights=[2.0, 1.0], types=["asymm", "asymm"])
    expected = [
        [0, 1, 1, 1],
        [1, 0, 1, 2 / 3],
        [1, 1, 0, 1 / 3],
        [1, 2 / 3, 1 / 3, 0],
    ]
    np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-14)


def test_raw_gower_constant_feature_distinction_from_public_features_is_preserved():
    features = np.array([[5.0, 0.0], [5.0, 10.0]])
    raw = gower_distance(features, weights=[1.0, 1.0])
    sequences = tx.seqdef([["a", "b"], ["b", "a"]], alphabet=["a", "b"])
    public = tx.seqcost(
        sequences, method="FEATURES", state_features=features, feature_weights=[1.0, 1.0]
    )["sm"]
    np.testing.assert_array_equal(raw, [[0, 1], [1, 0]])
    np.testing.assert_array_equal(np.asarray(public), [[0, 0.5], [0.5, 0]])


@pytest.mark.parametrize("weights", [[1e308, 1e308], [1e308, 5e307]])
def test_raw_gower_large_finite_weights_preserve_scale_invariance(weights):
    with np.errstate(all="raise"):
        actual = gower_distance(_FEATURES, weights=weights)
    expected = gower_distance(_FEATURES, weights=np.asarray(weights) / 1e308)
    np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-14)


def test_raw_gower_empty_rows_return_empty_matrix():
    assert gower_distance(np.empty((0, 2))).shape == (0, 0)


def test_raw_gower_all_missing_numeric_column_is_not_comparable():
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        actual = gower_distance([[np.nan, 0], [np.nan, 1]])
    np.testing.assert_array_equal(actual, [[0, 1], [1, 0]])


def test_raw_gower_finite_extreme_interval_preserves_scaled_distances():
    with np.errstate(all="raise"):
        actual = gower_distance([[-1e308], [0.0], [1e308]])
    np.testing.assert_array_equal(actual, [[0, 0.5, 1], [0.5, 0, 0.5], [1, 0.5, 0]])


@pytest.mark.parametrize("weights", [[1e308, 1e308], [1e308, 5e307]])
@pytest.mark.parametrize("dispatcher", [True, False])
def test_weighted_hamming_large_finite_weights_preserve_bounded_scale(weights, dispatcher):
    from traminex.distances import pairwise, python_kernels

    rows = np.array([[0, 0], [1, 1], [0, 1]], dtype=np.int32)
    function = pairwise if dispatcher else python_kernels.pairwise
    with np.errstate(all="raise"):
        actual = function(rows, method="WHAM", weights=weights)
    expected = function(rows, method="WHAM", weights=np.asarray(weights) / 1e308)
    np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-14)
    assert np.isfinite(actual).all()
    assert np.all((actual >= 0) & (actual <= 1))


@pytest.mark.parametrize("ignored", [5.0, np.nan])
@pytest.mark.parametrize("reverse_features", [False, True])
def test_raw_gower_only_comparable_positive_weight_cancels(ignored, reverse_features):
    features = np.array([[ignored, 0.0], [ignored, 1.0]])
    weights = np.array([1e308, 1e-300])
    if reverse_features:
        features = features[:, ::-1]
        weights = weights[::-1]
    original_weights = weights.copy()
    # The first column is excluded. The remaining distance is w*1/w = 1,
    # regardless of the weight of the excluded column or column order.
    with np.errstate(all="raise"):
        actual = gower_distance(features, weights=weights)
    np.testing.assert_array_equal(actual, [[0, 1], [1, 0]])
    np.testing.assert_array_equal(weights, original_weights)
    assert actual.dtype == np.dtype(float)


def test_raw_gower_large_unequal_weights_respect_pair_specific_missingness():
    # Each complete pair has a full-range difference in its only comparable
    # feature. No cross-pair has a comparable feature, so those entries are 0.
    features = [[np.nan, 0.0], [np.nan, 1.0], [0.0, np.nan], [1.0, np.nan]]
    with np.errstate(all="raise"):
        actual = gower_distance(features, weights=[1e308, 1e-300])
    np.testing.assert_array_equal(actual, [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])


def test_raw_gower_large_unequal_weights_exclude_joint_binary_absence():
    with np.errstate(all="raise"):
        actual = gower_distance([[0, 0], [0, 1]], weights=[1e308, 1e-300], types=["asymm", "asymm"])
    # The first binary feature is jointly absent and excluded for this pair.
    np.testing.assert_array_equal(actual, [[0, 1], [1, 0]])


@pytest.mark.parametrize("ignored", [5.0, np.nan])
def test_raw_array_features_substitution_costs_retain_small_effective_weight(ignored):
    with np.errstate(all="raise"):
        actual = tx.seqsubm(
            np.array([[0, 1], [1, 0]], dtype=np.int32),
            method="FEATURES",
            state_features=[[ignored, 0.0], [ignored, 1.0]],
            feature_weights=[1e308, 1e-300],
        )
    np.testing.assert_array_equal(actual, [[0, 1], [1, 0]])


def test_raw_gower_subnormal_positive_weights_preserve_fractional_mean():
    smallest = np.nextafter(0.0, 1.0)
    # Equal positive weights must give (1 + .5)/2, not round the .5 weighted
    # contribution to zero before division. Other pairs give .5 and .75.
    with np.errstate(all="raise"):
        actual = gower_distance(_FEATURES, weights=[smallest, smallest])
    np.testing.assert_array_equal(actual, [[0, 0.75, 0.5], [0.75, 0, 0.75], [0.5, 0.75, 0]])


def test_raw_gower_nearby_extreme_numeric_values_preserve_equal_spacing():
    values = [1e308]
    for _ in range(5):
        values.append(np.nextafter(values[-1], np.inf))
    # These six float64 values are equally spaced. Their normalized pairwise
    # distances are therefore the index differences divided by five.
    expected = np.abs(np.arange(6)[:, None] - np.arange(6)[None, :]) / 5
    with np.errstate(all="raise"):
        actual = gower_distance(np.array(values)[:, None])
    np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize(
    ("features", "weights", "expected"),
    [
        pytest.param([[0.0, 0.0], [1.0, 1.0]], [1e308, 1e308], 1.0, id="weights-near-the-maximum"),
        pytest.param([[-1e308], [1e308]], [1.0], 1.0, id="range-of-the-whole-double-line"),
        pytest.param(
            [[0.0, 0.0], [1.0, 0.5], [0.5, 1.0]], [5e-324, 5e-324], 0.75, id="subnormal-weights"
        ),
    ],
)
def test_raw_gower_hand_values_need_no_extended_precision(monkeypatch, features, weights, expected):
    # TX-AUD-20261011-035: with ``np.longdouble`` replaced by float64 inside the module, as on a
    # platform without an extended long double, the hand values still hold. The extension used
    # to rely on an 80-bit long double and gave nan, nan and 0.5 there.
    import types

    import traminex.distances.costs as costs

    proxy = types.ModuleType("numpy_with_a_float64_longdouble")
    proxy.__getattr__ = lambda name: getattr(np, name)  # type: ignore[attr-defined]
    proxy.longdouble = np.float64  # type: ignore[attr-defined]
    monkeypatch.setattr(costs, "np", proxy)
    with np.errstate(all="raise"):
        actual = costs.gower_distance(np.array(features, dtype=object), weights=np.array(weights))
    assert actual[0, 1] == expected


def test_new_distance_methods_md_describes_the_double_precision_scaling():
    from pathlib import Path

    doc = Path(__file__).resolve().parents[1] / "docs" / "NEW_DISTANCE_METHODS.md"
    text = " ".join(doc.read_text(encoding="utf-8").split())
    assert "extended precision" not in text
    assert "exact powers of two" in text
