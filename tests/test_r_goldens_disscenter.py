"""R-golden tests for ``disscenter`` (Phase 2.3 / A7) and the name-clash fix.

Anchors traminex's ``disscenter`` (distance of each case to its group centre)
to TraMineR 2.2-14 ``disscenter(diss, group=, medoids.index=, allcenter=,
weights=, squared=)`` on the shared Phase 2 fixture
(``scripts/r_goldens/make_goldens_disscenter.R``), and pins the resolution of
the historical name clash: matrix double-centering is now
``double_center_matrix`` (with deprecated alias ``disscenter_double_center``).
"""

from __future__ import annotations

import numpy as np
import pytest
from r_goldens_utils import R_GOLDEN_ATOL, load_golden

from traminex.diss import disscenter, disscenter_double_center, double_center_matrix

GROUPS = np.repeat([1, 2, 3], 5)


def _inputs():
    D = load_golden("phase2_input_D")
    wreal = load_golden("phase2_input_weights_real").ravel()
    return D, wreal


def test_disscenter_no_group_unweighted_matches_r():
    D, _ = _inputs()
    np.testing.assert_allclose(
        disscenter(D), load_golden("disscenter_nogrp_unw").ravel(), rtol=0, atol=R_GOLDEN_ATOL
    )


def test_disscenter_grouped_unweighted_matches_r():
    D, _ = _inputs()
    np.testing.assert_allclose(
        disscenter(D, group=GROUPS),
        load_golden("disscenter_grp_unw").ravel(),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


def test_disscenter_grouped_weighted_matches_r():
    D, wreal = _inputs()
    np.testing.assert_allclose(
        disscenter(D, group=GROUPS, weights=wreal),
        load_golden("disscenter_grp_wreal").ravel(),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


def test_disscenter_grouped_weighted_squared_matches_r():
    D, wreal = _inputs()
    np.testing.assert_allclose(
        disscenter(D, group=GROUPS, weights=wreal, squared=True),
        load_golden("disscenter_grp_wreal_sq").ravel(),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


def test_disscenter_allcenter_weighted_matches_r():
    D, wreal = _inputs()
    res = disscenter(D, group=GROUPS, weights=wreal, allcenter=True)
    assert list(res.columns) == ["1", "2", "3"]
    np.testing.assert_allclose(
        res.to_numpy(),
        load_golden("disscenter_grp_wreal_allcenter"),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


def test_disscenter_medoids_first_matches_r():
    D, wreal = _inputs()
    med = disscenter(D, group=GROUPS, weights=wreal, medoids_index="first")
    golden_1based = load_golden("disscenter_grp_wreal_medoids_first").ravel().astype(int)
    actual = np.array([med[lv] for lv in [1, 2, 3]], dtype=int)
    np.testing.assert_array_equal(actual, golden_1based - 1)  # R is 1-based
    # medoids_index="all" returns arrays containing the same minimisers.
    med_all = disscenter(D, group=GROUPS, weights=wreal, medoids_index="all")
    for lv in [1, 2, 3]:
        assert med[lv] in np.atleast_1d(med_all[lv])


def test_name_clash_resolved_double_centering_renamed():
    # The matrix double-centering helper formerly named disscenter.
    D, _ = _inputs()
    r = D.mean(axis=1, keepdims=True)
    c = D.mean(axis=0, keepdims=True)
    expected = D - r - c + D.mean()
    np.testing.assert_allclose(double_center_matrix(D), expected, rtol=0, atol=1e-12)
    with pytest.warns(DeprecationWarning, match="double_center_matrix"):
        np.testing.assert_allclose(disscenter_double_center(D), expected, rtol=0, atol=1e-12)
    # And disscenter itself now returns a per-case vector (R semantics), not
    # an (n, n) matrix.
    assert disscenter(D).shape == (D.shape[0],)
