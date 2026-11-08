"""Batch-12 dissimilarity-analysis, tree and representative items against R
(``scripts/r_goldens/make_goldens_diss_trees_reps.R``).

TX-AUD-20261019-046 and TX-AUD-20261011-004: R's group functions call ``factor()``, whose levels are the group's sorted
values written as R's ``as.character`` writes them, cases matched by that text: ``1``, not ``1.0``; ``1e+05``; and
0.3 and 0.1 + 0.2, equal to 15 significant digits, one level.
"""

from __future__ import annotations

import csv
import importlib
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx

GOLDEN_PATH = Path(__file__).resolve().parent / "r_goldens" / "diss_trees_reps.csv"
with open(GOLDEN_PATH, encoding="utf-8", newline="") as handle:
    GOLDEN = {(row["case"], row["field"]): row["value"] for row in csv.DictReader(handle)}


def _numbers(case: str, field: str) -> np.ndarray:
    return np.array([np.nan if v == "NA" else float(v) for v in GOLDEN[(case, field)].split()])


def _texts(case: str, field: str) -> list[str]:
    return GOLDEN[(case, field)].split(" ;; ")


def _square(case: str) -> np.ndarray:
    values = _numbers(case, "diss")
    size = int(round(values.size**0.5))
    return values.reshape(size, size, order="F")  # R's as.vector is column-major


D8 = _square("d8")
D16 = _square("d16")
GROUP_CASES = ["g_float", "g_decimal", "g_scientific", "g_merged", "g_factor_numeric"]


def _group(case: str):
    values = _numbers(case, "group")
    if (case, "levels") in GOLDEN:
        return pd.Categorical(values, categories=[float(v) for v in _texts(case, "levels")])
    return values


@pytest.mark.parametrize("case", GROUP_CASES)
def test_dissassoc_names_and_groups_levels_as_r(case: str) -> None:
    groups = tx.dissassoc(D8, _group(case), R=1).groups
    assert list(groups.index) == _texts(case, "dissassoc_rows")
    np.testing.assert_allclose(groups["n"].to_numpy(), _numbers(case, "dissassoc_n"))


@pytest.mark.parametrize("case", GROUP_CASES)
def test_disscenter_names_levels_as_r(case: str) -> None:
    assert list(tx.disscenter(D8, _group(case), allcenter=True).columns) == _texts(
        case, "allcenter_cols"
    )
    medoids = tx.disscenter(D8, _group(case), medoids_index="first")
    # The dict stays keyed by group value; its values in R's level order are R's medoids (1-based in R).
    assert [int(index) + 1 for index in medoids.values()] == [
        int(v) for v in _numbers(case, "medoids")
    ]


@pytest.mark.parametrize("case", GROUP_CASES)
def test_dissvar_grp_names_levels_as_r(case: str) -> None:
    result = tx.dissvar_grp(D8, _group(case))
    assert list(result) == _texts(case, "dissvar_grp_names")
    np.testing.assert_allclose(list(result.values()), _numbers(case, "dissvar_grp"), rtol=1e-12)


def test_an_integer_beyond_r_s_integer_range_is_a_double_level() -> None:
    # Batch-12 review M-5: R holds such a number as a double, so 123456789012345678 and 123456789012345679 are one level.
    from traminex._categorical_order import r_factor

    group = np.array([123456789012345678] * 4 + [123456789012345679] * 2 + [5] * 2, dtype=np.int64)
    assert r_factor(group)[0] == _texts("g_int64", "levels")
    groups = tx.dissassoc(D8, group, R=1).groups
    assert list(groups.index) == _texts("g_int64", "dissassoc_rows")
    np.testing.assert_allclose(groups["n"].to_numpy(), _numbers("g_int64", "dissassoc_n"))


def test_dissvar_grp_without_a_group_is_one_group_named_1() -> None:
    result = tx.dissvar_grp(D8)
    assert list(result) == _texts("dissvar_grp_default", "names")
    np.testing.assert_allclose(
        list(result.values()), _numbers("dissvar_grp_default", "value"), rtol=1e-12
    )


# Batch-12 review M-3: a rule's levels are read with R's string escapes, and an escape R does not know stops it.
ESCAPE_CASES = sorted(
    {case for case, field in GOLDEN if case.startswith("rule_escape_")},
    key=lambda case: int(case.rsplit("_", 1)[1]),
)


@pytest.mark.parametrize("case", ESCAPE_CASES)
def test_rule_levels_read_r_s_string_escapes(case: str) -> None:
    literal = GOLDEN[(case, "literal")]
    rules = [f"g %in% c('{literal}')", "g %in% c('w')"]
    if GOLDEN[(case, "codes")] == "ERROR":
        assert GOLDEN[(case, "assign")] == "ERROR"
        with pytest.raises(ValueError, match="cannot parse rule term"):
            tx.disstree_assign(rules, pd.DataFrame({"g": ["x", "w"]}))
    else:
        level = "".join(chr(int(code)) for code in _numbers(case, "codes"))
        result = tx.disstree_assign(rules, pd.DataFrame({"g": [level, "w"]}))
        np.testing.assert_array_equal(result, _numbers(case, "assign"))


# TX-AUD-20261010-079: R's dissassoc computes every group whose weight total is not zero, negative totals included (a
# zero total gives discrepancy 0), and its weighted residuals use weighted.mean, sum((x * w)[w != 0]) / sum(w), unguarded.
D9 = _square("d9")


@pytest.mark.parametrize("case", ["assoc_neg_total", "assoc_zero_total"])
def test_dissassoc_with_negative_or_zero_group_weight_totals_as_r(case: str) -> None:
    result = tx.dissassoc(
        D9,
        _numbers("d9", "group"),
        weights=_numbers(case, "weights"),
        R=1,
        weight_permutation="diss",
    )
    for actual, field in [
        (result.groups["n"], "groups_n"),
        (result.groups["discrepancy"], "groups_disc"),
        (result.anova_table["SS"], "anova_ss"),
        (result.anova_table["MSE"], "anova_mse"),
        (result.stat["t0"], "stat"),
    ]:
        np.testing.assert_allclose(
            actual.to_numpy(), _numbers(case, field), rtol=1e-12, err_msg=field
        )


# TX-AUD-20261010-083: R's tmrWeightedInertiaDist (weightedinertia.cpp:17-34) adds w[j] * w[i] * D[j, i], i outer and
# j > i inner (positions in the index list), to one running double and sums the weights the same way. On this input a
# pairwise sum differs from that order in the last bits.
D64 = _square("d64")
W64 = _numbers("d64", "weights")


def _c_loop_sum(D: np.ndarray, indices: list[int], w: np.ndarray) -> tuple[float, float]:
    """weightedinertia.cpp:17-31 written out: the running sum and the running weight total."""
    result = totweights = 0.0
    for position, i in enumerate(indices):
        iweight = w[i]
        totweights += iweight
        for j in indices[position + 1 :]:
            result += w[j] * iweight * D[j, i]
    return result, totweights


@pytest.mark.parametrize("block", [None, 5])
def test_weighted_inertia_and_dissvar_sum_in_r_s_order(block, monkeypatch) -> None:
    if block is not None:  # carry the running sum across many small product blocks
        diss_module = importlib.import_module("traminex.diss")
        monkeypatch.setattr(diss_module, "_PRODUCT_BLOCK", block)
        # The blocks are those of the NumPy path; the compiled one has none.
        monkeypatch.setattr(diss_module, "_compiled_pair_sum", lambda: None)
    every = list(range(64))
    subset = [int(v) - 1 for v in _numbers("inertia_order", "subset")]
    for indices, field in [(every, "all"), (subset, "subset_value")]:
        result, totweights = _c_loop_sum(D64, indices, W64)
        assert result / totweights == _numbers("inertia_order", field)[0]  # the loop is R's order
        assert tx.weighted_inertia(D64, indices, W64) == result / totweights
        # The undivided sum over i < j of w_i * w_j * D[i, j]: the same order on the transposed matrix.
        assert tx.weighted_pairwise_sum(D64.T, indices, W64) == result
    result, totweights = _c_loop_sum(D64, every, W64)
    assert tx.dissvar(D64, weights=W64) == result / totweights / totweights
    assert tx.dissvar(D64, weights=W64) == _numbers("inertia_order", "dissvar")[0]


# Batch-12 review M-4: disscenter subtracts R's weighted.mean, whose sums are R's long-double sum(); with three groups on
# the 64-case input a double sum differs from it in the last bits of 21 of the 64 values.
@pytest.mark.skipif(
    np.finfo(np.longdouble).nmant != 63, reason="R's sum() is an x87 80-bit long-double sum"
)
def test_disscenter_subtracts_r_s_weighted_mean() -> None:
    from traminex.diss import _disscenter_values

    group = _numbers("disscenter_64_3", "group")
    expected = _numbers("disscenter_64_3", "values")
    np.testing.assert_array_equal(np.asarray(tx.disscenter(D64, group, weights=W64)), expected)
    # dissassoc's own route to the same values.
    np.testing.assert_array_equal(_disscenter_values(D64, (group - 1).astype(int), W64), expected)


# TX-AUD-20261011-066: TraMineRextras' dissindic computes with a zero weight (that case's residual is 0/0, NaN), and stops
# on a character group, which traminex computes as R computes the same group as a factor (R-defect correction).
INDIC_GROUP = np.array(_texts("indic", "group"), dtype=object)


def test_dissindic_accepts_a_zero_weight_as_r() -> None:
    result = tx.dissindic(D8, INDIC_GROUP, weights=_numbers("indic", "weights"))
    for column in ("marginality", "gain"):
        np.testing.assert_allclose(
            result[column].to_numpy(),
            _numbers("indic_zero_weight", column),
            rtol=1e-10,
            atol=1e-12,
            err_msg=column,
        )


def test_dissindic_computes_a_character_group_where_r_stops() -> None:
    assert GOLDEN[("indic_character", "result")].startswith("ERROR")
    result = tx.dissindic(D8, INDIC_GROUP)
    for column in ("marginality", "gain"):
        np.testing.assert_allclose(
            result[column].to_numpy(),
            _numbers("indic_factor", column),
            rtol=1e-10,
            atol=1e-12,
            err_msg=column,
        )


# Batch-12 review I-1: with a group whose weights are all zero, the zero-weight cases take no part (their indicators
# are NaN) and the others' are dissindic's on the subset without them (the definition). R agrees when the group is not
# the reference level; when it is, R's qr.Q adds a direction outside the column space (R defect, recorded).
@pytest.mark.parametrize("case", ["indic_zero_group_ref", "indic_zero_group_other"])
def test_dissindic_with_a_zero_weight_group_follows_the_definition(case: str) -> None:
    weights = _numbers(case, "weights")
    positive = weights > 0
    result = tx.dissindic(D8, INDIC_GROUP, weights=weights)
    for column in ("marginality", "gain"):
        values = result[column].to_numpy()
        definition = _numbers(case, "subset_" + column)
        assert np.isnan(values[~positive]).all()
        np.testing.assert_allclose(values[positive], definition, rtol=1e-10, atol=1e-12)
        r_values = _numbers(case, column)
        if case == "indic_zero_group_other":
            np.testing.assert_allclose(values, r_values, rtol=1e-10, atol=1e-12, equal_nan=True)
        else:
            assert not np.allclose(r_values[positive], definition, rtol=1e-6)


@pytest.mark.parametrize("case", ["set_float", "set_merged"])
def test_dissCompare_names_sets_as_r(case: str) -> None:
    groups = np.array(["a", "b"] * 8, dtype=object)
    result = tx.dissCompare(D16, groups, set=_numbers(case, "set"), s=0)
    assert list(result.index) == _texts(case, "dissCompare_rows")


@pytest.mark.parametrize("case", ["rep_float", "rep_text"])
def test_seqrep_grp_names_levels_without_r_s_leading_space(case: str) -> None:
    # R's seqrep.grp names its levels through TraMineR:::group, which pastes "" and the level: " 1", " x"
    # (seqrep-grp.R:16). traminex names them as the levels themselves (R-defect correction, seqrep.grp row).
    r_names = _texts(case, "seqrep_grp_names")
    assert all(name.startswith(" ") for name in r_names)
    cells = np.array([list(row) for row in _texts("seq16", "cells")], dtype=object)
    sequences = tx.seqdef(cells, alphabet=["A", "B", "C"])
    group = np.repeat([1.0, 2.0], 8) if case == "rep_float" else np.repeat(["x", "y"], 8)
    assert list(tx.seqrep_grp(sequences, group=group, diss=D16)) == [name[1:] for name in r_names]


def test_dissmergegroups_numbers_groups_as_r() -> None:
    # dissmergegroups numbers groups with as.integer(factor(group)) (dissmergegroups.R:105): 0.3 and 0.1 + 0.2 are one.
    result = tx.dissmergegroups(D16, _numbers("merge_merged", "group"), silent=True)
    assert [int(v) for v in np.asarray(result["group"])] == [
        int(v) for v in _numbers("merge_merged", "membership")
    ]


# TX-AUD-20261011-055: R's disstree factorises a predictor (DTNGroupFactorBinary, disstree.R:404), so numbers equal to
# 15 significant digits are one level, and a numeric split's break is as.numeric of the level text (disstree.R:544).
D20 = _square("d20")


def test_a_tree_predictor_merges_levels_as_r() -> None:
    tree = tx.disstree(D20, pd.DataFrame({"x": _numbers("tree_merged", "x")}), R=1, min_size=2)
    assert (tree.root.split is not None) == (GOLDEN[("tree_merged", "split")] == "TRUE")


def test_a_numeric_break_is_read_back_from_its_level_text() -> None:
    tree = tx.disstree(D20, pd.DataFrame({"x": _numbers("tree_third", "x")}), R=1, min_size=2)
    assert tree.root.split.breaks == _numbers("tree_third", "breaks")[0]
    assert tx.disstree_get_rules(tree) == _texts("tree_third", "rules")


# TX-AUD-20261019-047: a root-only tree (max_depth 1) has R's NA leaf labels and one NA rule, which assigns no profile.
def _root_only_tree():
    x = _numbers("tree_third", "x")
    return tx.disstree(D20, pd.DataFrame({"x": x}), R=1, min_size=2, max_depth=1), x


def test_a_root_only_tree_has_r_s_na_labels_and_rule() -> None:
    tree, _ = _root_only_tree()
    labels = tx.disstreeleaf(tree, label=True)
    assert [label is None for label in labels] == [
        flag == "TRUE" for flag in GOLDEN[("tree_root", "leaf_labels_na")].split()
    ]
    rules = tx.disstree_get_rules(tree)
    assert [rule is None for rule in rules] == [
        flag == "TRUE" for flag in GOLDEN[("tree_root", "rules_na")].split()
    ]


def test_a_root_only_tree_s_rule_assigns_no_profile() -> None:
    tree, x = _root_only_tree()
    # R's rules carry their covariates (here x) as an attribute; a Python list names them through covar.
    profile = pd.DataFrame({"x": x[[0, 10, 19]]})
    expected = _numbers("tree_root", "assign")
    covariates = _texts("tree_root", "rules_covariates")
    for covar in (None, covariates):
        np.testing.assert_array_equal(
            tx.disstree_assign(tx.disstree_get_rules(tree), profile, covar=covar), expected
        )


# TX-AUD-20261010-065: R parses each rule (disstree.R:166-167), so ' & ' inside a quoted level is part of the level.
D6 = _square("d6")


@pytest.mark.parametrize("case", ["rule_amp", "rule_quote"])
def test_a_tree_s_rules_on_a_level_with_an_ampersand_assign_as_r(case: str) -> None:
    group = np.array(_texts(case, "g"), dtype=object)
    tree = tx.disstree(D6, pd.DataFrame({"g": group}), R=1, min_size=2)
    rules = tx.disstree_get_rules(tree)
    assert rules == _texts(case, "rules")
    profile = pd.DataFrame({"g": group})
    if GOLDEN[(case, "assign")] == "ERROR":  # R cannot parse its own rule for a level with a quote
        with pytest.raises(ValueError, match="cannot parse rule term"):
            tx.disstree_assign(rules, profile)
    else:
        np.testing.assert_array_equal(tx.disstree_assign(rules, profile), _numbers(case, "assign"))


@pytest.mark.parametrize("case", ["rule_escaped", "rule_double", "rule_doubled"])
def test_hand_written_rules_read_r_string_literals(case: str) -> None:
    profile = pd.DataFrame({"g": np.array(_texts("rule_quote", "g"), dtype=object)})
    rules = _texts(case, "rules")
    if GOLDEN[(case, "assign")] == "ERROR":
        with pytest.raises(ValueError, match="cannot parse rule term"):
            tx.disstree_assign(rules, profile)
    else:
        np.testing.assert_array_equal(tx.disstree_assign(rules, profile), _numbers(case, "assign"))


# TX-REP-20261022-006: the grouped plots' panels and legend follow factor(group) (seqedplot.R:53,
# seqplot-tentrop.R:10): 0.3 and 0.1 + 0.2 are one level, 1 is named 1.
def test_plot_panels_and_legend_follow_r_factor() -> None:
    from traminex.visualization.series import _factor_panels, grouped_transversal_entropy

    group = _numbers("plot_levels", "group")
    rows, titles = _factor_panels(group)
    assert titles == _texts("plot_levels", "levels")
    assert rows == [[2, 3], [0, 1]]
    sequences = tx.seqdef(
        np.array([list("AAB"), list("ABB"), list("BBA"), list("BAA")], dtype=object)
    )
    labels, _ = grouped_transversal_entropy(sequences, group)
    assert list(labels) == _texts("plot_levels", "levels")


# TX-AUD-20261011-023: R's tmrWeightedInertiaContrib (weightedinertia.cpp:68-101) adds each case's terms in index order,
# so identical sequences tie exactly and the first is the medoid. The weighted mvad tree (OM, constant costs) and dissrf
# on an index-order sortv.
@pytest.fixture(scope="module")
def mvad():
    alphabet = _texts("mvad", "alphabet")
    cells = np.array(
        [[alphabet[int(c) - 1] for c in row] for row in GOLDEN[("mvad", "codes")].split()],
        dtype=object,
    )
    weights = _numbers("mvad", "weights")
    sequences = tx.seqdef(cells, alphabet=alphabet, weights=weights)
    diss = np.asarray(tx.seqdist(sequences, method="OM", subm="CONSTANT", indel=1), dtype=float)
    position = np.arange(1, diss.shape[0] + 1, dtype=float)
    sums = [diss.sum(), (diss**2).sum(), (diss * np.outer(position, position)).sum()]
    assert sums == list(_numbers("mvad", "diss_sums"))  # integer distances: R's matrix exactly
    predictors = pd.DataFrame(
        {
            name: pd.Categorical(_texts("mvad", name), categories=_texts("mvad", name + "_levels"))
            for name in ["male", "funemp", "gcse5eq", "Grammar", "fmpr", "livboth"]
        }
    )
    return diss, weights, predictors


def test_tree_medoids_break_ties_as_r(mvad) -> None:
    diss, weights, predictors = mvad
    tree = tx.disstree(
        diss,
        predictors,
        weights=weights,
        R=1,
        weight_permutation="diss",
        min_size=0.05,
        max_depth=5,
    )
    names, medoids = [], []

    def walk(node, name):
        names.append(name)
        medoids.append(int(node.info["medoid"]) + 1)
        if node.split is not None:
            walk(node.kids[0], name + "L")
            walk(node.kids[1], name + "R")

    walk(tree.root, "N")
    assert names == _texts("tree_medoids", "names")
    assert medoids == [int(v) for v in _numbers("tree_medoids", "medoids")]


@pytest.mark.parametrize("case", ["rf_weighted", "rf_unweighted"])
def test_dissrf_medoids_break_ties_as_r(mvad, case: str) -> None:
    diss, weights, _ = mvad
    result = tx.dissrf(
        diss,
        k=10,
        weights=weights if case == "rf_weighted" else None,
        sortv=np.arange(1, diss.shape[0] + 1),
    )
    assert [int(v) + 1 for v in result["medoids"]] == [int(v) for v in _numbers(case, "medoids")]
    np.testing.assert_allclose(result["R2"], _numbers(case, "R2")[0], rtol=1e-12)


# The controls of dissrf's former near-tie tolerance (np.isclose, rtol and atol 1e-12): R takes the exact first minimum,
# the middle of [0, s, 2s] at every scale and the third case where the second's contribution is larger by 4e-15.
def test_rf_medoid_is_r_s_exact_first_minimum_at_every_scale() -> None:
    base = np.abs(np.arange(3)[:, None] - np.arange(3)[None, :]).astype(float)
    for scale, medoid in zip(_numbers("rf_near_tie", "scales"), _numbers("rf_near_tie", "medoids")):
        result = tx.dissrf(base * scale, k=1, sortv=None, weights=np.ones(3))
        assert int(result["medoids"][0]) + 1 == int(medoid), scale


def test_a_contribution_larger_by_4e_15_is_not_a_tie() -> None:
    from traminex.repseq_rf import _medoid_index

    diss = _square("medoid_last_ulp")
    medoid = int(_numbers("medoid_last_ulp", "medoid")[0])
    assert int(tx.disscenter(diss, medoids_index="first", weights=np.ones(4))) + 1 == medoid
    assert _medoid_index(diss, weights=np.ones(4), squared=False) + 1 == medoid


# TX-AUD-20261011-052: R's cmdscale has no tie test, so dissrf's default axis does not depend on the scale of the
# distances; traminex's repeated-root test was absolute (atol 1e-12) and fell back to index order at scale 1e-10.
def test_dissrf_default_axis_is_scale_free() -> None:
    alphabet = _texts("mvad", "alphabet")
    cells = np.array(
        [[alphabet[int(c) - 1] for c in row] for row in GOLDEN[("mvad", "codes")].split()[:80]],
        dtype=object,
    )
    diss = np.asarray(tx.seqdist(tx.seqdef(cells, alphabet=alphabet), method="LCS"), dtype=float)
    assert [diss.sum(), (diss**2).sum()] == list(_numbers("rf_scale", "diss_sums"))
    for scale, r2 in zip(_numbers("rf_scale", "scales"), _numbers("rf_scale", "R2")):
        np.testing.assert_allclose(
            tx.dissrf(diss * scale, k=8)["R2"], r2, rtol=1e-10, err_msg=str(scale)
        )


# TX-REP-20261021-001: dissrf orders a text or categorical sortv as R's order(sortv) does (dissrf.R:94): text in R's
# collation, a factor by its level codes, missing values last.
@pytest.fixture(scope="module")
def mvad80_lcs():
    alphabet = _texts("mvad", "alphabet")
    cells = np.array(
        [[alphabet[int(c) - 1] for c in row] for row in GOLDEN[("mvad", "codes")].split()[:80]],
        dtype=object,
    )
    return np.asarray(tx.seqdist(tx.seqdef(cells, alphabet=alphabet), method="LCS"), dtype=float)


@pytest.mark.parametrize(
    "case", ["rf_sortv_text", "rf_sortv_factor", "rf_sortv_text_na", "rf_sortv_factor_na"]
)
def test_dissrf_orders_text_and_categorical_sortv_as_r(mvad80_lcs, case: str) -> None:
    values = np.array(_texts("rf_sortv", "text"), dtype=object)
    if case.endswith("_na"):
        values[[int(v) - 1 for v in _numbers("rf_sortv", "na")]] = None
    sortv = (
        pd.Categorical(values, categories=_texts("rf_sortv", "levels"))
        if "factor" in case
        else values
    )
    result = tx.dissrf(mvad80_lcs, k=4, sortv=sortv)
    assert [int(v) + 1 for v in result["medoids"]] == [int(v) for v in _numbers(case, "medoids")]
    np.testing.assert_allclose(result["R2"], _numbers(case, "R2")[0], rtol=1e-12)


# TX-AUD-20261010-080: R's seqtree passes seqdist.args as given (seqtree.R:15-16): list() calls seqdist without a method
# and R stops; only the default is list(method = "LCS", norm = "auto").
def test_seqtree_passes_an_empty_seqdist_args_through() -> None:
    cells = np.array([list(row) for row in _texts("seqtree_args", "cells")], dtype=object)
    sequences = tx.seqdef(cells, alphabet=["A", "B", "C"])
    predictors = pd.DataFrame({"g": _texts("seqtree_args", "g")})
    assert GOLDEN[("seqtree_args", "empty")].startswith("ERROR")
    with pytest.raises(ValueError, match="must be specified"):
        tx.seqtree(sequences, predictors, R=10, seqdist_args={})
    for kwargs in ({}, {"seqdist_args": None}):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            tree = tx.seqtree(sequences, predictors, R=10, **kwargs)
        assert tree.root.info["n"] == _numbers("seqtree_args", "default_root_n")[0]


def test_dissrf_reads_a_category_dtype_series_sortv_by_its_codes(mvad80_lcs) -> None:
    sortv = pd.Series(
        pd.Categorical(_texts("rf_sortv", "text"), categories=_texts("rf_sortv", "levels"))
    )
    result = tx.dissrf(mvad80_lcs, k=4, sortv=sortv)
    assert [int(v) + 1 for v in result["medoids"]] == [
        int(v) for v in _numbers("rf_sortv_factor", "medoids")
    ]
