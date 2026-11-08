"""Tests for Gower mixed-type dissimilarity and the FEATURES substitution-cost method.

This tests the retained raw-array extension, which drops constant numeric
columns. The STSList FEATURES path follows cluster::daisy and tests its distinct
constant-column/missing-pair semantics separately. Expected values here are
computed by hand from the raw extension's Gower formula:

    d(i,j) = sum_f w_f * delta_ijf * d_ijf / sum_f w_f * delta_ijf

with numeric d = |xi-xj|/range_f, nominal d = 0/1, ordinal = rank-scaled numeric,
asymmetric-binary delta = 0 when both absent.
"""

import numpy as np
import pandas as pd
import pytest

from traminex.distances.costs import gower_distance


def test_gower_single_numeric_feature():
    # range = 10; pairwise |delta|/10
    feats = np.array([[0.0], [10.0], [5.0]])
    D = gower_distance(feats, types=["numeric"])
    assert D.shape == (3, 3)
    assert np.allclose(np.diag(D), 0.0)
    assert np.allclose(D, D.T)
    assert D[0, 1] == pytest.approx(1.0)
    assert D[0, 2] == pytest.approx(0.5)
    assert D[1, 2] == pytest.approx(0.5)


def test_gower_single_nominal_feature():
    feats = np.array([["a"], ["b"], ["a"]], dtype=object)
    D = gower_distance(feats, types=["nominal"])
    assert D[0, 1] == pytest.approx(1.0)
    assert D[0, 2] == pytest.approx(0.0)
    assert D[1, 2] == pytest.approx(1.0)


def test_gower_mixed_numeric_and_nominal_equal_weights():
    # col0 numeric (range 10), col1 nominal
    feats = np.array([[0.0, "a"], [10.0, "b"], [0.0, "a"]], dtype=object)
    D = gower_distance(feats, types=["numeric", "nominal"])
    assert D[0, 1] == pytest.approx((1.0 + 1.0) / 2)  # both features maximally differ
    assert D[0, 2] == pytest.approx(0.0)
    assert D[1, 2] == pytest.approx((1.0 + 1.0) / 2)


def test_gower_feature_weights():
    # numeric range 10 differs fully; nominal identical -> weighted mean
    feats = np.array([[0.0, "a"], [10.0, "a"]], dtype=object)
    D = gower_distance(feats, types=["numeric", "nominal"], weights=[2.0, 1.0])
    assert D[0, 1] == pytest.approx((2.0 * 1.0 + 1.0 * 0.0) / 3.0)  # = 2/3


def test_gower_ordinal_rank_scaled():
    feats = np.array([["low"], ["med"], ["high"]], dtype=object)
    D = gower_distance(feats, types=["ordinal"], ordinal_levels=[["low", "med", "high"]])
    assert D[0, 2] == pytest.approx(1.0)
    assert D[0, 1] == pytest.approx(0.5)
    assert D[1, 2] == pytest.approx(0.5)


def test_gower_asymmetric_binary_ignores_joint_absence():
    # two asymmetric-binary features; state2 vs state0 both-zero on feat1 -> ignored
    feats = np.array([[1, 0], [0, 1], [0, 0]], dtype=float)
    D = gower_distance(feats, types=["asymm", "asymm"])
    # s0=(1,0) s1=(0,1): feat0 differ(1), feat1 differ(1) -> both counted -> (1+1)/2=1
    assert D[0, 1] == pytest.approx(1.0)
    # s0=(1,0) s2=(0,0): feat0 differ(1) counted; feat1 both-zero -> ignored -> 1/1 = 1
    assert D[0, 2] == pytest.approx(1.0)
    # s1=(0,1) s2=(0,0): feat0 both-zero ignored; feat1 differ(1) -> 1/1 = 1
    assert D[1, 2] == pytest.approx(1.0)


def test_gower_constant_numeric_feature_contributes_zero():
    # zero-range numeric feature: contributes nothing (excluded from denominator)
    feats = np.array([[5.0, 0.0], [5.0, 10.0]], dtype=float)
    D = gower_distance(feats, types=["numeric", "numeric"])
    # col0 range 0 -> skipped; col1 range 10 -> |0-10|/10 = 1 -> mean over comparable = 1
    assert D[0, 1] == pytest.approx(1.0)


@pytest.mark.parametrize("missing", [np.nan, None, pd.NA])
def test_gower_missing_nominal_feature_is_excluded_for_affected_pairs(missing):
    actual = gower_distance([[0, "a"], [1, missing], [1, "b"]], types=["numeric", "nominal"])
    # Pair 1/2 has equal numeric values and no comparable nominal value.
    np.testing.assert_array_equal(actual, [[0, 1, 1], [1, 0, 0], [1, 0, 0]])


@pytest.mark.parametrize("missing", [np.nan, None, pd.NA])
@pytest.mark.parametrize("levels", [None, [["low", "high"]]])
def test_gower_missing_ordinal_feature_is_excluded_from_ranks_and_pairs(missing, levels):
    actual = gower_distance(
        [[0, "low"], [1, missing], [1, "high"]], types=["numeric", "ordinal"], ordinal_levels=levels
    )
    np.testing.assert_array_equal(actual, [[0, 1, 1], [1, 0, 0], [1, 0, 0]])


@pytest.mark.parametrize("missing", [np.nan, None, pd.NA])
def test_gower_missing_asymmetric_binary_feature_is_not_a_presence(missing):
    actual = gower_distance([[0, 1], [0, missing], [1, 0]], types=["asymm", "asymm"])
    # Pair 0/1 has joint absence in feature one and missing feature two.
    np.testing.assert_array_equal(actual, [[0, 0, 1], [0, 0, 1], [1, 1, 0]])


@pytest.mark.parametrize("value", [np.inf, -np.inf])
def test_gower_rejects_infinite_numeric_features(value):
    with pytest.raises(ValueError):
        gower_distance([[0.0], [value]])


@pytest.mark.parametrize("value", [-1.0, 2.0, np.inf, -np.inf])
def test_gower_asymmetric_binary_requires_observed_values_in_zero_one(value):
    with pytest.raises(ValueError):
        gower_distance([[0.0], [value]], types=["asymm"])
