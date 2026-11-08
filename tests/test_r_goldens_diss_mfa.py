"""R-golden tests for dissmfacw (weighted multi-factor analysis) -- B-1 / E-1.

Anchors traminex's weighted distance/covariate association to TraMineR 2.2-14
``dissmfacw`` (see scripts/r_goldens/make_goldens_diss_mfa.R). The deterministic
decomposition (PseudoF, PseudoR2, and explained SS per term + Total) is goldened; the
permutation p-values are stochastic and checked separately with a fixed seed
(plan decision 0.2.2).
"""

from __future__ import annotations

import numpy as np
import pytest
from r_goldens_utils import R_GOLDEN_ATOL, load_golden

from traminex import dissmfacw


def _inputs():
    D = load_golden("phase2_input_D")
    w = load_golden("phase2_input_weights_real").ravel()
    grp = load_golden("phase2_input_grp").ravel().astype(int)
    i = np.arange(1, len(grp) + 1)
    h = np.where(i % 2 == 0, "even", "odd")
    return D, w, grp, h


def _deterministic_table(result):
    return np.column_stack([result.pseudo_F, result.pseudo_R2, result.SS])


def test_dissmfacw_deterministic_decomposition_matches_r():
    D, w, grp, h = _inputs()
    # Numeric arrays are continuous.  Strings request R-style factor coding.
    res = dissmfacw(
        D,
        factors={"g": grp.astype(str), "h": h},
        weights=w,
        nperm=1,
    )
    # Rows in factor order then Total: g, h, Total.
    assert list(res.variables) == ["g", "h", "Total"]
    np.testing.assert_allclose(
        _deterministic_table(res),
        load_golden("dissmfacw_det"),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


def test_dissmfacw_numeric_predictor_matches_r_as_one_continuous_column():
    D, _, grp, _ = _inputs()
    x = np.linspace(-1, 1, len(grp))
    res = dissmfacw(D, factors={"x": x}, nperm=1)

    assert list(res.variables) == ["x", "Total"]
    np.testing.assert_allclose(
        _deterministic_table(res),
        load_golden("dissmfacw_numeric_unweighted"),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


@pytest.mark.parametrize(
    ("weights_tag", "golden_name"),
    [
        ("unweighted", "dissmfacw_mixed_unweighted"),
        ("weighted", "dissmfacw_mixed_weighted"),
    ],
)
def test_dissmfacw_mixed_predictors_match_r_in_insertion_order(weights_tag, golden_name):
    D, w, grp, _ = _inputs()
    x = np.linspace(-1, 1, len(grp))
    weights = None if weights_tag == "unweighted" else w
    res = dissmfacw(
        D,
        factors={"g": grp.astype(str), "x": x},
        weights=weights,
        nperm=1,
    )

    assert list(res.variables) == ["g", "x", "Total"]
    np.testing.assert_allclose(
        _deterministic_table(res),
        load_golden(golden_name),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


def test_dissmfacw_exact_clustassoc_mixed_case_matches_r():
    D = np.ones((6, 6), dtype=float)
    np.fill_diagonal(D, 0)
    D[:3, :3] = 0.1
    D[3:, 3:] = 0.1
    np.fill_diagonal(D, 0)
    clusters = np.array(["1", "1", "1", "2", "2", "2"])
    covariate = np.array([0, 0.1, 0.2, 2, 2.1, 2.2])

    res = dissmfacw(
        D,
        factors={"clustf": clusters, "covar": covariate},
        nperm=1,
    )

    assert list(res.variables) == ["clustf", "covar", "Total"]
    np.testing.assert_allclose(
        _deterministic_table(res),
        load_golden("dissmfacw_clustassoc_mixed"),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


def test_dissmfacw_predictor_and_weight_validation_matches_r_contract():
    D = np.ones((6, 6), dtype=float)
    np.fill_diagonal(D, 0)

    # R accepts continuous constants but rejects a factor with one level. A
    # constant adds no rank to the intercept: R's values for it are an
    # artefact of its column count (TX-AUD-20261010-040), here it has no test.
    with pytest.warns(UserWarning, match="rank 1 for 2 columns"):
        constant = dissmfacw(D, factors={"x": np.ones(6)}, nperm=1)
    assert np.isnan(constant.mfac["PseudoF"]).all()
    with pytest.raises(ValueError, match="categorical.*two levels"):
        dissmfacw(D, factors={"g": np.repeat("one", 6)}, nperm=1)

    with pytest.raises(ValueError, match="length 5, expected 6"):
        dissmfacw(D, factors={"x": np.arange(5)}, nperm=1)
    with pytest.raises(ValueError, match="weights must have length 6"):
        dissmfacw(D, factors={"x": np.arange(6)}, weights=np.ones(5), nperm=1)
    with pytest.raises(ValueError, match="weights must be finite"):
        dissmfacw(
            D,
            factors={"x": np.arange(6)},
            weights=np.array([1, 1, np.nan, 1, 1, 1]),
            nperm=1,
        )
    with pytest.raises(ValueError, match="weights must be non-negative"):
        dissmfacw(
            D,
            factors={"x": np.arange(6)},
            weights=np.array([1, 1, -1, 1, 1, 1]),
            nperm=1,
        )
    with pytest.raises(ValueError, match="numeric predictor 'x' must be finite"):
        dissmfacw(D, factors={"x": np.array([0, 1, np.inf, 3, 4, 5])}, nperm=1)

    # TraMineR drops predictor-NA rows and subsets D and weights identically.
    with_missing = np.array([0, 1, np.nan, 3, 4, 5])
    weights = np.array([1, 2, np.nan, 4, 5, 6], dtype=float)
    kept = np.array([0, 1, 3, 4, 5])
    dropped = dissmfacw(D, factors={"x": with_missing}, weights=weights, nperm=1)
    explicit = dissmfacw(
        D[np.ix_(kept, kept)],
        factors={"x": with_missing[kept]},
        weights=weights[kept],
        nperm=1,
    )
    np.testing.assert_allclose(
        _deterministic_table(dropped),
        _deterministic_table(explicit),
        rtol=0,
        atol=1e-12,
    )

    # Zero case weights are accepted by TraMineR.
    zero_weight = dissmfacw(
        D,
        factors={"x": np.arange(6)},
        weights=np.array([1, 1, 0, 1, 1, 1]),
        nperm=1,
    )
    assert np.isfinite(_deterministic_table(zero_weight)).all()


def test_dissmfacw_permutation_pvalue_fixed_seed():
    # p-values are stochastic: pin with a fixed seed + loose checks (decision 0.2.2).
    D, w, grp, h = _inputs()
    factors = {"g": grp.astype(str), "h": h}
    res1 = dissmfacw(D, factors=factors, weights=w, nperm=200, seed=1)
    res2 = dissmfacw(D, factors=factors, weights=w, nperm=200, seed=1)
    assert np.all((res1.p_value >= 0) & (res1.p_value <= 1))
    # Same seed -> identical p-values (deterministic given the seed).
    np.testing.assert_array_equal(res1.p_value, res2.p_value)
    # The strongly-associated factor g should be far more significant than h.
    g_p, h_p = res1.p_value[0], res1.p_value[1]
    assert g_p < h_p


def test_dissmfacw_deterministic_F_is_seed_independent():
    # The deterministic PseudoF/PseudoR2 must not depend on the permutation seed.
    D, w, grp, h = _inputs()
    factors = {"g": grp.astype(str), "h": h}
    a = dissmfacw(D, factors=factors, weights=w, nperm=50, seed=1)
    b = dissmfacw(D, factors=factors, weights=w, nperm=50, seed=99)
    np.testing.assert_allclose(a.pseudo_F, b.pseudo_F, rtol=0, atol=R_GOLDEN_ATOL)
    np.testing.assert_allclose(a.pseudo_R2, b.pseudo_R2, rtol=0, atol=R_GOLDEN_ATOL)


def test_dissvar_grp_matches_live_r():
    # RP4: per-group discrepancy (TraMineRextras dissvar.grp), on phase2 inputs.
    from traminex import dissvar_grp

    D = load_golden("phase2_input_D")
    grp = load_golden("phase2_input_grp").ravel().astype(int)
    w = load_golden("phase2_input_weights_real").ravel()

    uw = dissvar_grp(D, grp)
    np.testing.assert_allclose(
        list(uw.values()), load_golden("dissvar_grp_uw").ravel(), rtol=0, atol=R_GOLDEN_ATOL
    )
    sq = dissvar_grp(D, grp, squared=True)
    np.testing.assert_allclose(
        list(sq.values()), load_golden("dissvar_grp_uw_sq").ravel(), rtol=0, atol=R_GOLDEN_ATOL
    )
    wt = dissvar_grp(D, grp, weights=w)
    np.testing.assert_allclose(
        list(wt.values()), load_golden("dissvar_grp_w").ravel(), rtol=0, atol=R_GOLDEN_ATOL
    )


def test_dissmergegroups_matches_live_r():
    # RP4: greedy ASW group merging (TraMineR dissmergegroups + WeightedCluster ASW).
    from traminex import dissmergegroups
    from traminex.diss_merge import _wc_asw

    D = load_golden("phase2_input_D")
    g0 = load_golden("dissmerge_input_group").ravel().astype(int)
    w = load_golden("phase2_input_weights_real").ravel()

    # Direct ASW pin: traminex _wc_asw == WeightedCluster wcClusterQuality ASW.
    assert _wc_asw(D, g0, np.ones(len(g0))) == pytest.approx(
        float(load_golden("dissmerge_asw_uw").ravel()[0]), abs=R_GOLDEN_ATOL
    )
    assert _wc_asw(D, g0, w) == pytest.approx(
        float(load_golden("dissmerge_asw_w").ravel()[0]), abs=R_GOLDEN_ATOL
    )

    # End-to-end merged partition + final quality.
    for tag, weights in [("uw", None), ("w", w)]:
        mg = dissmergegroups(D, g0, weights=weights, min_group=2, small=0.05)
        np.testing.assert_array_equal(
            mg["group"], load_golden(f"dissmerge_{tag}_group").ravel().astype(int)
        )
        assert mg["quality"] == pytest.approx(
            float(load_golden(f"dissmerge_{tag}_qual").ravel()[0]), abs=R_GOLDEN_ATOL
        )


def test_dissindic_matches_live_r():
    # RP4: individual marginality / gain (TraMineRextras dissindic), phase2 inputs.
    from traminex import dissindic

    D = load_golden("phase2_input_D")
    grp = load_golden("phase2_input_grp").ravel().astype(int)
    w = load_golden("phase2_input_weights_real").ravel()

    for tag, weights, squared in [("uw", None, False), ("w", w, False), ("w_sq", w, True)]:
        r = dissindic(D, grp, weights=weights, squared=squared)
        np.testing.assert_allclose(
            r["marginality"].to_numpy(dtype=float),
            load_golden(f"dissindic_{tag}_marg").ravel(),
            rtol=0,
            atol=R_GOLDEN_ATOL,
        )
        np.testing.assert_allclose(
            r["gain"].to_numpy(dtype=float),
            load_golden(f"dissindic_{tag}_gain").ravel(),
            rtol=0,
            atol=R_GOLDEN_ATOL,
        )
