"""R-golden tests for dissCompare / seqLRT / seqBIC (RP4, s=0).

Anchors the LRT / BIC group comparison to live TraMineRextras 0.6.9 (see
scripts/r_goldens/make_goldens_group_compare.R). Only the deterministic s=0
path is pinned. The seqCompare/seqLRT/seqBIC cases are added in their commit.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDEN_ATOL, R_GOLDENS_DIR, load_golden

from traminex import dissCompare, seqBIC, seqCompare, seqdef, seqLRT
from traminex import group_compare as group_compare_module


def _inputs():
    D = load_golden("phase2_input_D")
    grp = load_golden("phase2_input_grp").ravel().astype(int)
    w = load_golden("phase2_input_weights_real").ravel()
    g2 = np.where(grp <= 1, 1, 2)
    return D, g2, w


@pytest.mark.parametrize("tag,weighted", [("uw", False), ("w", True)])
def test_disscompare_matches_live_r(tag, weighted):
    D, g2, w = _inputs()
    weights = w if weighted else None
    out = dissCompare(D, g2, s=0, stat="all", weighted=weighted, weights=weights)
    assert list(out.columns) == ["LRT", "p-value", "Delta BIC", "Bayes Factor"]
    np.testing.assert_allclose(
        out.to_numpy(dtype=float).ravel(),
        load_golden(f"dissCompare_{tag}").ravel(),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


def test_disscompare_rejects_three_groups():
    D, _, _ = _inputs()
    g3 = np.array([0, 1, 2] * 5)
    with pytest.raises(ValueError):
        dissCompare(D, g3)


def _constant_inputs():
    D = np.ones((12, 12), dtype=float)
    np.fill_diagonal(D, 0)
    group = np.array(["a"] * 7 + ["b"] * 5, dtype=object)
    weights = np.array([1, 2, 1.5, 0.5, 3, 1, 2.5, 1, 4, 0.75, 2, 1.25])
    return D, group, weights


def _r3_expected(case_name: str) -> pd.DataFrame:
    records = pd.read_csv(
        R_GOLDENS_DIR / "group_compare_r3_cases.csv",
        dtype=str,
        keep_default_na=False,
    )
    selected = records.loc[(records["case"] == case_name) & (records["record_type"] == "result")]
    row_order = list(dict.fromkeys(selected["row"]))
    column_order = list(dict.fromkeys(selected["column"]))
    result = pd.DataFrame(index=row_order, columns=column_order, dtype=float)
    for record in selected.itertuples(index=False):
        result.loc[record.row, record.column] = float(record.value)
    return result


@pytest.mark.parametrize(
    "case_name,kwargs",
    [
        ("s0_unweighted_default", {"s": 0, "weighted": False}),
        ("s0_weighted_default", {"s": 0, "weighted": True, "use_weights": True}),
        ("s0_by_group", {"s": 0, "weighted": "by.group", "use_weights": True}),
        ("s0_squared_false", {"s": 0, "squared": False}),
        ("s0_squared_true", {"s": 0, "squared": True}),
        ("s0_lrt_inc_df", {"s": 0, "stat": "LRT", "inc_df": True}),
        ("s0_bic_inc_df", {"s": 0, "stat": "BIC", "inc_df": True}),
        ("sample_below_both_bf", {"s": 3, "seed": 1401}),
        ("sample_below_avg_bf", {"s": 3, "seed": 1401, "BFopt": 1}),
        ("sample_below_from_bic", {"s": 3, "seed": 1401, "BFopt": 2}),
        ("sample_at_min", {"s": 5, "seed": 1402, "inc_df": True}),
        ("sample_at_max", {"s": 7, "seed": 1403, "inc_df": True}),
        ("sample_above_max", {"s": 9, "seed": 1404, "inc_df": True}),
        ("sample_weighted_unit", {"s": 3, "seed": 1405, "weighted": True}),
        ("sample_by_group_unit", {"s": 3, "seed": 1405, "weighted": "by.group"}),
    ],
)
def test_disscompare_r3_routes_match_exchangeable_r_goldens(case_name, kwargs):
    D, group, weights = _constant_inputs()
    call_kwargs = dict(kwargs)
    use_weights = call_kwargs.pop("use_weights", False)
    if use_weights:
        call_kwargs["weights"] = weights
    elif case_name in {"sample_weighted_unit", "sample_by_group_unit"}:
        call_kwargs["weights"] = np.ones(12)
    result = dissCompare(D, group, **call_kwargs)
    expected = _r3_expected(case_name)
    assert list(result.index) == list(expected.index)
    assert list(result.columns) == list(expected.columns)
    if case_name.startswith("sample_below"):
        # R and NumPy intentionally use different random bytes. The exact R
        # values remain frozen above; compare deterministic aggregation here.
        from scipy import stats as spstats

        assert result.loc["1", "p-value"] == pytest.approx(
            spstats.chi2.sf(result.loc["1", "LRT"], 3), rel=0, abs=R_GOLDEN_ATOL
        )
        if "Bayes Factor (From Avg BIC)" in result:
            assert result.loc["1", "Bayes Factor (From Avg BIC)"] == pytest.approx(
                np.exp(result.loc["1", "Delta BIC"] / 2), rel=0, abs=R_GOLDEN_ATOL
            )
        return
    np.testing.assert_allclose(result, expected, rtol=0, atol=2e-12)


def _multiple_set_inputs():
    D, _, _ = _constant_inputs()
    group = np.array(
        ["a", "a", "a", "a", "b", "b", "b", "a", "a", "a", "b", "b"],
        dtype=object,
    )
    sets = np.array(["early"] * 7 + ["late"] * 5, dtype=object)
    return D, group, sets


@pytest.mark.parametrize(
    "case_name,kwargs",
    [
        ("multiple_sets_s0", {"s": 0}),
        ("multiple_sets_sample", {"s": 2, "seed": 1406, "inc_df": True}),
    ],
)
def test_disscompare_multiple_sets_row_order_and_aggregation(case_name, kwargs):
    D, group, sets = _multiple_set_inputs()
    result = dissCompare(D, group, set=sets, **kwargs)
    expected = _r3_expected(case_name)
    assert list(result.index) == ["early", "late"]
    assert list(result.columns) == list(expected.columns)
    if case_name == "multiple_sets_sample":
        from scipy import stats as spstats

        np.testing.assert_array_equal(result["df"], [3, 2])
        np.testing.assert_allclose(
            result["p-value"],
            [
                spstats.chi2.sf(result.loc["early", "LRT"], 3),
                spstats.chi2.sf(result.loc["late", "LRT"], 2),
            ],
            rtol=0,
            atol=2e-15,
        )
        np.testing.assert_allclose(
            result["Bayes Factor (From Avg BIC)"],
            np.exp(result["Delta BIC"] / 2),
            rtol=0,
            atol=2e-15,
        )
        return
    np.testing.assert_allclose(result, expected, rtol=0, atol=2e-12)


def test_disscompare_mixed_case_set_order_matches_r() -> None:
    D = np.ones((8, 8), dtype=float)
    np.fill_diagonal(D, 0)
    group = np.asarray(["g1", "g2"] * 4, dtype=object)
    sets = np.asarray(["a"] * 4 + ["A"] * 4, dtype=object)

    result = dissCompare(D, group, set=sets, s=0, weighted=False)
    expected = _r3_expected("mixed_case_sets")

    assert list(result.index) == list(expected.index) == ["a", "A"]
    np.testing.assert_allclose(result, expected, rtol=0, atol=2e-12)


def test_disscompare_locale_set_order_matches_r() -> None:
    """Set levels take R's collation (TX-AUD-20261011-015); they used a portable key."""
    r_order = ["_", "-", "!", "?", ".", "*", "%", "a", "A", "á", "ä", "b"]
    D = np.ones((4 * len(r_order), 4 * len(r_order)), dtype=float)
    np.fill_diagonal(D, 0)
    group = np.asarray(["g1", "g2"] * (2 * len(r_order)), dtype=object)
    sets = np.repeat(r_order, 4)

    result = dissCompare(D, group, set=sets, s=0, weighted=False)
    expected = _r3_expected("locale_case_sets")

    assert list(expected.index) == r_order
    assert list(result.index) == r_order
    np.testing.assert_allclose(result, expected, rtol=0, atol=2e-12)


def test_disscompare_na_removal_message_and_result_match_r():
    D, group, sets = _multiple_set_inputs()
    group[1] = None
    sets[8] = None
    transcript = (R_GOLDENS_DIR / "group_compare_r3_errors.txt").read_text(encoding="utf-8")
    assert "na_removed|message|[!!] 2 cases removed" in transcript
    with pytest.warns(UserWarning, match="2 cases removed"):
        result = dissCompare(D, group, set=sets, s=0)
    np.testing.assert_allclose(result, _r3_expected("na_removed"), rtol=0, atol=2e-12)


def test_sampling_plan_repeats_balanced_draws_and_records_r_index_quirk():
    # R's plan draws the first level's sample from 1:max(n) even when that level
    # is the smaller group (the facts below); traminex draws the larger group
    # first, as seqCompare does (tests/test_r_goldens_disscompare_order.py,
    # TX-AUD-20261011-022). The draw mechanics themselves are shared.
    records = pd.read_csv(
        R_GOLDENS_DIR / "group_compare_r3_cases.csv",
        dtype=str,
        keep_default_na=False,
    )
    plan = records.loc[records["case"] == "reversed_unbalanced_plan"]
    facts = dict(zip(plan["column"], plan["value"], strict=True))
    assert facts == {
        "first_group_size": "3",
        "second_group_size": "5",
        "sample_size": "3",
        "draw_count": "2",
        "first_relative_max": "5",
        "second_relative_max": "3",
    }

    first, second = group_compare_module._build_sampling_plan(5, 3, 3, np.random.default_rng(1410))
    assert first.shape == second.shape == (2, 3)
    assert first.max() == 4
    assert second.max() == 2
    assert sorted(np.bincount(first.ravel(), minlength=5).tolist()) == [1, 1, 1, 1, 2]
    assert sorted(np.bincount(second.ravel(), minlength=3).tolist()) == [2, 2, 2]


def test_all_sampling_plans_are_built_before_first_statistic(monkeypatch):
    D, group, sets = _multiple_set_inputs()
    built = []
    original_build = group_compare_module._build_sampling_plan
    original_stat = group_compare_module._seqxcomp

    def tracking_build(*args, **kwargs):
        built.append(args[:3])
        return original_build(*args, **kwargs)

    def checking_stat(*args, **kwargs):
        assert len(built) == 2
        return original_stat(*args, **kwargs)

    monkeypatch.setattr(group_compare_module, "_build_sampling_plan", tracking_build)
    monkeypatch.setattr(group_compare_module, "_seqxcomp", checking_stat)
    dissCompare(D, group, set=sets, s=2, seed=22)


def test_sampling_seed_generator_and_global_rng_isolation():
    points = np.array([[i, (i % 3) * 0.7] for i in range(12)], dtype=float)
    D = np.linalg.norm(points[:, None, :] - points[None, :, :], axis=2)
    _, group, _ = _constant_inputs()
    kwargs = {"s": 3, "stat": "all", "weighted": False}
    first = dissCompare(D, group, seed=31, **kwargs)
    second = dissCompare(D, group, seed=31, **kwargs)
    np.testing.assert_array_equal(first, second)

    rng_a = np.random.default_rng(32)
    rng_b = np.random.default_rng(32)
    a1 = dissCompare(D, group, seed=rng_a, **kwargs)
    a2 = dissCompare(D, group, seed=rng_a, **kwargs)
    b1 = dissCompare(D, group, seed=rng_b, **kwargs)
    b2 = dissCompare(D, group, seed=rng_b, **kwargs)
    np.testing.assert_array_equal(a1, b1)
    np.testing.assert_array_equal(a2, b2)

    np.random.seed(713)
    expected = np.random.random(4)
    np.random.seed(713)
    dissCompare(D, group, seed=33, **kwargs)
    np.testing.assert_array_equal(np.random.random(4), expected)


@pytest.mark.parametrize(
    "kwargs,match",
    [
        ({"s": -1}, "s"),
        ({"s": 1.5}, "s"),
        ({"weighted": "bad"}, "weighted"),
        ({"squared": "bad"}, "squared"),
        ({"stat": "bad"}, "stat"),
        ({"BFopt": 3}, "BFopt"),
        ({"inc_df": 1}, "inc_df"),
    ],
)
def test_disscompare_option_validation(kwargs, match):
    D, group, _ = _constant_inputs()
    with pytest.raises((TypeError, ValueError), match=match):
        dissCompare(D, group, **kwargs)


def test_disscompare_shape_and_length_validation():
    D, group, weights = _constant_inputs()
    with pytest.raises(ValueError, match="square"):
        dissCompare(D[:, :-1], group)
    with pytest.raises(ValueError, match="group length"):
        dissCompare(D, group[:-1])
    with pytest.raises(ValueError, match="set length"):
        dissCompare(D, group, set=np.ones(11))
    with pytest.raises(ValueError, match="weights"):
        dissCompare(D, group, weights=weights[:-1])


def test_bfopt_three_is_an_explicit_python_refusal_of_r_fallthrough():
    transcript = (R_GOLDENS_DIR / "group_compare_r3_errors.txt").read_text(encoding="utf-8")
    assert "bfopt_three|success" in transcript
    D, group, _ = _constant_inputs()
    with pytest.raises(ValueError, match="BFopt"):
        dissCompare(D, group, s=3, BFopt=3)


# seqCompare / seqLRT / seqBIC (sequence inputs, s=0).
_SALPH = ["1", "2", "3", "4"]


def _seq_fixture():
    codes = load_golden("seqcompare_input_codes").astype(int)
    w = load_golden("seqcompare_input_weights").ravel()
    chars = np.array([[_SALPH[c] for c in row] for row in codes], dtype=object)
    x = seqdef(chars, alphabet=_SALPH, weights=w)
    g = load_golden("seqcompare_input_group").ravel().astype(int)
    return x, g


@pytest.mark.parametrize("tag,weighted", [("uw", False), ("w", True)])
def test_seqlrt_seqbic_seqcompare_match_live_r(tag, weighted):
    x, g = _seq_fixture()
    lrt = seqLRT(x, group=g, s=0, method="LCS", weighted=weighted)
    bic = seqBIC(x, group=g, s=0, method="LCS", weighted=weighted)
    cmp = seqCompare(x, group=g, s=0, stat="all", method="LCS", weighted=weighted)

    np.testing.assert_allclose(
        lrt.to_numpy(dtype=float).ravel(),
        load_golden(f"seqLRT_{tag}").ravel(),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )
    np.testing.assert_allclose(
        bic.to_numpy(dtype=float).ravel(),
        load_golden(f"seqBIC_{tag}").ravel(),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )
    np.testing.assert_allclose(
        cmp.to_numpy(dtype=float).ravel(),
        load_golden(f"seqCompare_{tag}").ravel(),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


def _r3_sequence_fixture():
    states = [f"S{index:02d}" for index in range(1, 13)]
    chars = np.asarray([[state, state] for state in states], dtype=object)
    _, _, weights = _constant_inputs()
    x = seqdef(chars, alphabet=states, weights=weights)
    group = np.array(["a"] * 7 + ["b"] * 5, dtype=object)
    multi_group = np.array(
        ["a", "a", "a", "a", "b", "b", "b", "a", "a", "a", "b", "b"],
        dtype=object,
    )
    sets = np.array(["early"] * 7 + ["late"] * 5, dtype=object)
    return x, group, multi_group, sets


@pytest.mark.parametrize(
    "case_name,kwargs",
    [
        ("seq_group_s0", {"s": 0}),
        ("seq_group_s0_unweighted", {"s": 0, "weighted": False}),
        ("seq_group_s0_by_group", {"s": 0, "weighted": "by.group"}),
        ("seq_group_sample_opt1", {"s": 3, "seed": 1501, "opt": 1, "weighted": False}),
        ("seq_group_sample_opt2", {"s": 3, "seed": 1501, "opt": 2, "weighted": False}),
        ("seq_group_sample_auto1", {"s": 3, "seed": 1501, "weighted": False}),
        ("seq_group_sample_auto2", {"s": 7, "seed": 1502, "weighted": False}),
    ],
)
def test_seqcompare_group_sampling_options_match_exchangeable_r(case_name, kwargs):
    x, group, _, _ = _r3_sequence_fixture()

    result = seqCompare(x, group=group, method="LCS", **kwargs)

    expected = _r3_expected(case_name)
    assert list(result.index) == list(expected.index)
    assert list(result.columns) == list(expected.columns)
    if kwargs["s"] == 0:
        np.testing.assert_allclose(result, expected, rtol=0, atol=3e-12)
    else:
        # R and NumPy intentionally use different random bytes. The complete R
        # values remain frozen; compare exact deterministic aggregation facts.
        draw_count = 3 if kwargs["s"] == 3 else 1
        from scipy import stats as spstats

        assert result.loc["1", "p-value"] == pytest.approx(
            spstats.chi2.sf(result.loc["1", "LRT"], draw_count), rel=0, abs=R_GOLDEN_ATOL
        )
        bf_column = "Bayes Factor (From Avg BIC)" if draw_count > 1 else "Bayes Factor"
        assert result.loc["1", bf_column] == pytest.approx(
            np.exp(result.loc["1", "Delta BIC"] / 2), rel=0, abs=R_GOLDEN_ATOL
        )


def test_seqcompare_sampling_options_are_numerically_equivalent_in_python():
    x, group, _, _ = _r3_sequence_fixture()
    common = {"group": group, "s": 3, "seed": 1501, "weighted": False, "method": "LCS"}

    option_one = seqCompare(x, opt=1, **common)
    option_two = seqCompare(x, opt=2, **common)
    automatic = seqCompare(x, **common)

    pd.testing.assert_frame_equal(
        option_one, option_two, check_exact=False, rtol=0, atol=R_GOLDEN_ATOL
    )
    pd.testing.assert_frame_equal(
        option_one, automatic, check_exact=False, rtol=0, atol=R_GOLDEN_ATOL
    )


def test_seqcompare_group_set_route_matches_r():
    x, _, group, sets = _r3_sequence_fixture()

    result = seqCompare(
        x,
        group=group,
        set=sets,
        s=2,
        seed=1503,
        inc_df=True,
        weighted=False,
        method="LCS",
    )

    expected = _r3_expected("seq_group_sets_sample")
    assert list(result.index) == ["early", "late"]
    assert list(result.columns) == list(expected.columns)
    np.testing.assert_array_equal(result["df"], expected["df"])
    from scipy import stats as spstats

    np.testing.assert_allclose(
        result["p-value"],
        [
            spstats.chi2.sf(result.loc["early", "LRT"], 3),
            spstats.chi2.sf(result.loc["late", "LRT"], 2),
        ],
        rtol=0,
        atol=3e-15,
    )
    np.testing.assert_allclose(
        result["Bayes Factor (From Avg BIC)"],
        np.exp(result["Delta BIC"] / 2),
        rtol=0,
        atol=3e-15,
    )


def test_seqcompare_two_object_and_paired_list_routes_match_r():
    x, _, _, _ = _r3_sequence_fixture()
    chars = np.asarray(x.alphabet, dtype=object)[x.data]
    x1 = seqdef(chars[:7], alphabet=x.alphabet, weights=x.weights[:7])
    x2 = seqdef(chars[7:], alphabet=x.alphabet, weights=x.weights[7:])

    two_s0 = seqCompare(x1, x2, s=0, method="LCS")
    two_sample = seqCompare(x1, x2, s=3, seed=1504, opt=2, BFopt=2, weighted=False, method="LCS")
    list1 = [
        seqdef(chars[[0, 1, 2, 3]], alphabet=x.alphabet, weights=x.weights[[0, 1, 2, 3]]),
        seqdef(chars[[4, 5, 6]], alphabet=x.alphabet, weights=x.weights[[4, 5, 6]]),
    ]
    list2 = [
        seqdef(chars[[7, 8, 9]], alphabet=x.alphabet, weights=x.weights[[7, 8, 9]]),
        seqdef(chars[[10, 11]], alphabet=x.alphabet, weights=x.weights[[10, 11]]),
    ]
    paired_s0 = seqCompare(list1, list2, s=0, method="LCS")
    paired_sample = seqCompare(
        list1,
        list2,
        s=2,
        seed=1505,
        BFopt=1,
        inc_df=True,
        weighted=False,
        method="LCS",
    )

    for actual, case_name in [
        (two_s0, "seq_two_object_s0"),
        (two_sample, "seq_two_object_sample"),
        (paired_s0, "seq_paired_lists_s0"),
        (paired_sample, "seq_paired_lists_sample"),
    ]:
        expected = _r3_expected(case_name)
        assert list(actual.index) == list(expected.index)
        assert list(actual.columns) == list(expected.columns)
        if case_name.endswith("s0"):
            np.testing.assert_allclose(actual, expected, rtol=0, atol=3e-12)
        else:
            if "df" in actual:
                np.testing.assert_array_equal(actual["df"], expected["df"])
            bf_column = next(column for column in actual if column.startswith("Bayes Factor"))
            if "From Avg BIC" in bf_column:
                np.testing.assert_allclose(
                    actual[bf_column], np.exp(actual["Delta BIC"] / 2), rtol=0
                )


def test_seqlrt_and_seqbic_are_exact_seqcompare_routes_and_match_r():
    x, group, _, _ = _r3_sequence_fixture()
    common = {
        "group": group,
        "s": 3,
        "seed": 1506,
        "opt": 2,
        "inc_df": True,
        "weighted": False,
        "method": "LCS",
    }

    lrt = seqLRT(x, **common)
    bic = seqBIC(x, BFopt=2, **common)

    pd.testing.assert_frame_equal(
        lrt, seqCompare(x, stat="LRT", **common), check_exact=False, rtol=0, atol=R_GOLDEN_ATOL
    )
    pd.testing.assert_frame_equal(
        bic,
        seqCompare(x, stat="BIC", BFopt=2, **common),
        check_exact=False,
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )
    assert list(lrt.columns) == list(_r3_expected("seq_lrt_wrapper").columns)
    assert list(bic.columns) == list(_r3_expected("seq_bic_wrapper").columns)
    np.testing.assert_array_equal(lrt["df"], [3])
    np.testing.assert_array_equal(bic["df"], [3])
    np.testing.assert_allclose(
        bic["Bayes Factor (From Avg BIC)"], np.exp(bic["Delta BIC"] / 2), rtol=0
    )


@pytest.mark.parametrize(
    "s,opt,expected_calls",
    [(3, 1, 3), (3, 2, 1), (3, None, 3), (7, None, 1), (0, 1, 1)],
)
def test_seqcompare_option_controls_public_seqdist_construction(
    monkeypatch, s, opt, expected_calls
):
    x, group, _, _ = _r3_sequence_fixture()
    calls = []
    original = group_compare_module.seqdist

    def tracking_seqdist(sequence, method="OM", **kwargs):
        calls.append((sequence.nseq, method, dict(kwargs)))
        return original(sequence, method=method, **kwargs)

    monkeypatch.setattr(group_compare_module, "seqdist", tracking_seqdist)
    seqCompare(x, group=group, s=s, seed=1501, opt=opt, method="LCS")

    assert len(calls) == expected_calls
    assert all(call[1] == "LCS" for call in calls)


def test_seqcompare_forwards_distance_kwargs_and_default_lcs(monkeypatch):
    x, group, _, _ = _r3_sequence_fixture()
    calls = []
    original = group_compare_module.seqdist

    def tracking_seqdist(sequence, method="OM", **kwargs):
        calls.append((method, dict(kwargs)))
        return original(sequence, method=method, **kwargs)

    monkeypatch.setattr(group_compare_module, "seqdist", tracking_seqdist)
    implicit = seqCompare(x, group=group, s=0, weighted=False)
    explicit = seqCompare(x, group=group, s=0, weighted=False, method="LCS")
    pd.testing.assert_frame_equal(implicit, explicit, check_exact=False, rtol=0, atol=R_GOLDEN_ATOL)
    assert calls[0][0] == "LCS"

    subm = np.ones((12, 12), dtype=float) - np.eye(12)
    seqCompare(
        x,
        group=group,
        s=0,
        weighted=False,
        method="HAM",
        subm=subm,
        norm="maxlength",
    )
    assert calls[-1] == (
        "HAM",
        {"weighted": False, "subm": subm, "norm": "maxlength"},
    )

    default_sampling = seqCompare(x, group=group, weighted=False, seed=1513)
    explicit_sampling = seqCompare(x, group=group, s=100, seed=1513, weighted=False, method="LCS")
    pd.testing.assert_frame_equal(
        default_sampling, explicit_sampling, check_exact=False, rtol=0, atol=R_GOLDEN_ATOL
    )


def test_seqcompare_builds_all_set_plans_before_any_distance(monkeypatch):
    x, _, group, sets = _r3_sequence_fixture()
    built = []
    original_build = group_compare_module._build_sampling_plan
    original_distance = group_compare_module.seqdist

    def tracking_build(*args, **kwargs):
        built.append(args[:3])
        return original_build(*args, **kwargs)

    def checking_distance(*args, **kwargs):
        assert len(built) == 2
        return original_distance(*args, **kwargs)

    monkeypatch.setattr(group_compare_module, "_build_sampling_plan", tracking_build)
    monkeypatch.setattr(group_compare_module, "seqdist", checking_distance)

    seqCompare(x, group=group, set=sets, s=2, opt=2, method="LCS")


def test_seqcompare_seed_reproducibility_and_global_rng_isolation():
    x, group, _, _ = _r3_sequence_fixture()
    data_snapshot = x.data.copy()
    weights_snapshot = x.weights.copy()
    group_snapshot = group.copy()
    kwargs = {"group": group, "s": 3, "opt": 1, "method": "LCS"}
    first = seqCompare(x, seed=1510, **kwargs)
    second = seqCompare(x, seed=1510, **kwargs)
    np.testing.assert_array_equal(first, second)

    rng_a = np.random.default_rng(1511)
    rng_b = np.random.default_rng(1511)
    np.testing.assert_array_equal(
        seqCompare(x, seed=rng_a, **kwargs), seqCompare(x, seed=rng_b, **kwargs)
    )

    np.random.seed(714)
    expected = np.random.random(4)
    np.random.seed(714)
    seqCompare(x, seed=1512, **kwargs)
    np.testing.assert_array_equal(np.random.random(4), expected)
    np.testing.assert_array_equal(x.data, data_snapshot)
    np.testing.assert_array_equal(x.weights, weights_snapshot)
    np.testing.assert_array_equal(group, group_snapshot)


def test_seqcompare_na_removal_and_input_validation():
    x, group, _, sets = _r3_sequence_fixture()
    group[1] = None
    sets[8] = None
    with pytest.warns(UserWarning, match="2 sequences removed"):
        result = seqCompare(x, group=group, set=sets, s=0, method="LCS")
    assert list(result.index) == ["early", "late"]

    transcript = (R_GOLDENS_DIR / "group_compare_r3_errors.txt").read_text(encoding="utf-8")
    assert "seq_list_length|simpleError" in transcript
    assert "seq_list_element|simpleError" in transcript
    assert "seq_bad_opt|success" in transcript
    with pytest.raises(ValueError, match="cannot both be None"):
        seqCompare(x, s=0)
    with pytest.raises(ValueError, match="set.*group"):
        seqCompare(x, set=sets, s=0)
    with pytest.raises(ValueError, match="same length"):
        seqCompare([x, x], [x], s=0)
    with pytest.raises(TypeError, match="element"):
        seqCompare([x, object()], [x, x], s=0)
    with pytest.raises(TypeError, match="seqdata2"):
        seqCompare(x, object(), s=0)
    with pytest.raises(ValueError, match="opt"):
        seqCompare(x, group=group, s=3, opt=9)
    with pytest.raises(ValueError, match="group.*seqdata2"):
        seqCompare(x, x, group=group, s=0)
    with pytest.raises(ValueError, match="group length"):
        seqCompare(x, group=group[:-1], s=0)
    with pytest.raises(ValueError, match="set length"):
        seqCompare(x, group=group, set=sets[:-1], s=0)
