"""R-golden tests for dissrf / seqrf (RP4 T13b + the star seqrf).

Anchors traminex's relative-frequency grouping (grp_meth="prop") to live
TraMineR 2.2-14 ``dissrf`` (see scripts/r_goldens/make_goldens_diss_rf.R), using
R's cmdscale sort key (dumped) to avoid MDS sign drift. The seqrf cases are added
in their commit.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest
from r_goldens_utils import R_GOLDEN_ATOL, load_golden, load_golden_lines

from traminex import seqdef, seqrf
from traminex.repseq_rf import dissrf

R2_CASES = Path(__file__).with_name("r_goldens") / "diss_rf_r2_cases.csv"


def _inputs():
    D = load_golden("phase2_input_D")
    sortv = load_golden("dissrf_input_sortv").ravel()
    w = load_golden("phase2_input_weights_real").ravel()
    return D, sortv, w


@pytest.mark.parametrize("tag,weighted", [("uw", False), ("w", True)])
def test_dissrf_prop_matches_live_r(tag, weighted):
    D, sortv, w = _inputs()
    weights = w if weighted else None
    rf = dissrf(D, k=3, sortv=sortv, weights=weights, grp_meth="prop")

    np.testing.assert_array_equal(
        rf["medoids"], load_golden(f"dissrf_{tag}_medoids").ravel().astype(int)
    )
    stats = np.array([rf["R2"], rf["Fstat"], rf["pvalue"]], dtype=float)
    np.testing.assert_allclose(
        stats, load_golden(f"dissrf_{tag}_stats").ravel(), rtol=0, atol=R_GOLDEN_ATOL
    )
    np.testing.assert_allclose(
        rf["heights"], load_golden(f"dissrf_{tag}_heights").ravel(), rtol=0, atol=R_GOLDEN_ATOL
    )


def test_dissrf_k1_returns_overall_medoid_without_partition_statistics():
    D, sortv, w = _inputs()
    rf = dissrf(D, k=1, sortv=sortv, weights=w, grp_meth="prop")

    np.testing.assert_array_equal(
        rf["medoids"], load_golden("dissrf_w_k1_medoid").ravel().astype(int)
    )
    assert rf["R2"] == 0.0
    assert rf["Fstat"] == 0.0
    assert np.isnan(rf["pvalue"])
    sizes = load_golden("dissrf_w_k1_sizes").ravel()
    np.testing.assert_allclose(
        [rf["sizes"]["wsum"], rf["sizes"]["gsize"]], sizes, rtol=0, atol=R_GOLDEN_ATOL
    )


# seqrf (the star item): wraps dissrf on a sequence dissimilarity matrix.
_RF_ALPH = ["1", "2", "3", "4"]


def _seqrf_fixture():
    codes = load_golden("seqrf_input_codes").astype(int)
    w = load_golden("seqrf_input_weights").ravel()
    chars = np.array([[_RF_ALPH[c] for c in row] for row in codes], dtype=object)
    return seqdef(chars, alphabet=_RF_ALPH, weights=w)


@pytest.mark.parametrize("tag,weighted", [("uw", False), ("w", True)])
def test_seqrf_matches_live_r(tag, weighted):
    x = _seqrf_fixture()
    D = load_golden("seqrf_input_diss")
    sortv = load_golden("seqrf_input_sortv").ravel()
    srf = seqrf(x, D, k=3, sortv=sortv, weighted=weighted)

    np.testing.assert_array_equal(
        srf["medoids"], load_golden(f"seqrf_{tag}_medoids").ravel().astype(int)
    )
    np.testing.assert_allclose(
        srf["heights"], load_golden(f"seqrf_{tag}_heights").ravel(), rtol=0, atol=R_GOLDEN_ATOL
    )
    np.testing.assert_array_equal(
        srf["seqtoplot"], load_golden(f"seqrf_{tag}_repcodes").astype(int)
    )


# String sortv resolution (deterministic, no MDS): "from.start" / "from.end" order
# the sequences lexicographically from the start / end (R order(as.data.frame(...))).
@pytest.mark.parametrize("sortv,tag", [("from.start", "fromstart"), ("from.end", "fromend")])
def test_seqrf_string_sortv_matches_live_r(sortv, tag):
    x = _seqrf_fixture()
    D = load_golden("seqrf_input_diss")
    srf = seqrf(x, D, k=3, sortv=sortv, weighted=True)

    np.testing.assert_array_equal(
        srf["medoids"], load_golden(f"seqrf_{tag}_medoids").ravel().astype(int)
    )
    np.testing.assert_allclose(
        srf["heights"], load_golden(f"seqrf_{tag}_heights").ravel(), rtol=0, atol=R_GOLDEN_ATOL
    )
    np.testing.assert_array_equal(
        srf["seqtoplot"], load_golden(f"seqrf_{tag}_repcodes").astype(int)
    )


# R2 complete grouping contract ------------------------------------------------
_D_CASES = {
    "d_prop_unweighted_numeric": dict(grp_meth="prop", sort="numeric", k=3),
    "d_prop_integer_numeric": dict(grp_meth="prop", sort="numeric", k=3, weights="integer"),
    "d_prop_fractional_numeric": dict(grp_meth="prop", sort="numeric", k=3, weights="fractional"),
    "d_prop_fractional_default_mds": dict(
        grp_meth="prop",
        sort="mds",
        k=3,
        weights="mds_counterexample",
        exact_mds=True,
        default_mds=True,
    ),
    "d_first_numeric_sq0_default": dict(grp_meth="first", sort="numeric"),
    "d_first_ties_sq1_default": dict(grp_meth="first", sort="ties", squared=True),
    "d_first_null_sq0_pow2": dict(grp_meth="first", sort="null", pow=2),
    "d_first_mds_sq1_pow1": dict(grp_meth="first", sort="mds", squared=True, pow=1),
    "d_random_numeric_sq0_default": dict(grp_meth="random", sort="numeric", seed=6101),
    "d_random_ties_sq1_default": dict(grp_meth="random", sort="ties", squared=True, seed=6102),
    "d_random_null_sq0_pow2": dict(grp_meth="random", sort="null", pow=2, seed=6103),
    "d_random_mds_sq1_pow1": dict(grp_meth="random", sort="mds", squared=True, pow=1, seed=6104),
    "d_k1_forces_prop": dict(grp_meth="first", sort="numeric", k=1),
}

_S_CASES = {
    "s_prop_unweighted_numeric": dict(grp_meth="prop", sort="numeric", k=3),
    "s_prop_integer_numeric": dict(grp_meth="prop", sort="numeric", k=3, weights="integer"),
    "s_prop_fractional_numeric": dict(grp_meth="prop", sort="numeric", k=3, weighted=True),
    "s_prop_fractional_default_mds": dict(
        grp_meth="prop",
        sort="mds",
        k=3,
        weights="mds_counterexample",
        weighted=True,
        exact_mds=True,
        default_mds=True,
    ),
    "s_first_numeric_sq0_default": dict(grp_meth="first", sort="numeric"),
    "s_first_ties_sq1_default": dict(grp_meth="first", sort="ties", squared=True),
    "s_first_null_sq0_pow2": dict(grp_meth="first", sort="null", pow=2),
    "s_first_mds_sq1_pow1": dict(grp_meth="first", sort="mds", squared=True, pow=1),
    "s_random_numeric_sq0_default": dict(grp_meth="random", sort="numeric", seed=7101),
    "s_random_ties_sq1_default": dict(grp_meth="random", sort="ties", squared=True, seed=7102),
    "s_random_null_sq0_pow2": dict(grp_meth="random", sort="null", pow=2, seed=7103),
    "s_random_mds_sq1_pow1": dict(grp_meth="random", sort="mds", squared=True, pow=1, seed=7104),
    "s_k1_forces_prop": dict(grp_meth="random", sort="numeric", k=1),
}


def _r2_rows():
    with R2_CASES.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _case_rows(case):
    return [row for row in _r2_rows() if row["case"] == case]


def _section(rows, name):
    return sorted(
        (row for row in rows if row["section"] == name),
        key=lambda row: (int(row["row"]), int(row["col"])),
    )


def _vector(rows, name):
    return np.asarray([float(row["value"]) for row in _section(rows, name)])


def _group_lists(rows, name):
    selected = _section(rows, name)
    groups = sorted({int(row["row"]) for row in selected})
    return [
        np.asarray([float(row["value"]) for row in selected if int(row["row"]) == group])
        for group in groups
    ]


def _sort_value(kind, numeric, n):
    if kind == "numeric":
        return numeric
    if kind == "ties":
        return np.resize(np.asarray([0, 0, 1, 1, 2], dtype=float), n)
    if kind == "null":
        return None
    return "mds"


def _run_d_case(config):
    D, numeric, fractional = _inputs()
    if config.get("exact_mds"):
        D = load_golden("dissrf_weighted_mds_input_D")
    weights = None
    if config.get("weights") == "integer":
        weights = np.resize(np.asarray([1, 2, 3], dtype=float), D.shape[0])
    elif config.get("weights") == "fractional":
        weights = fractional
    elif config.get("weights") == "mds_counterexample":
        weights = load_golden("dissrf_weighted_mds_input_weights").ravel()
    kwargs = {
        "k": config.get("k", 4),
        "sortv": _sort_value(config["sort"], numeric, D.shape[0]),
        "weights": weights,
        "grp_meth": config["grp_meth"],
        "squared": config.get("squared", False),
        "pow": config.get("pow"),
    }
    if "seed" in config:
        kwargs["random_state"] = config["seed"]
    if config.get("default_mds"):
        kwargs.pop("sortv")
    return dissrf(D, **kwargs)


def _run_s_case(config):
    x = _seqrf_fixture()
    D = load_golden("seqrf_input_diss")
    numeric = load_golden("seqrf_input_sortv").ravel()
    weights = None
    if config.get("exact_mds"):
        codes = load_golden("seqrf_input_codes").astype(int)[:10]
        chars = np.array([[_RF_ALPH[c] for c in row] for row in codes], dtype=object)
        x = seqdef(chars, alphabet=_RF_ALPH)
        D = load_golden("dissrf_weighted_mds_input_D")
    if config.get("weights") == "integer":
        weights = np.resize(np.asarray([3, 2, 1], dtype=float), D.shape[0])
    elif config.get("weights") == "mds_counterexample":
        weights = load_golden("dissrf_weighted_mds_input_weights").ravel()
    kwargs = {
        "k": config.get("k", 4),
        "sortv": _sort_value(config["sort"], numeric, D.shape[0]),
        "weights": weights,
        "weighted": config.get("weighted", False),
        "grp_meth": config["grp_meth"],
        "squared": config.get("squared", False),
        "pow": config.get("pow"),
    }
    if "seed" in config:
        kwargs["random_state"] = config["seed"]
    if config.get("default_mds"):
        kwargs.pop("sortv")
    return seqrf(x, D, **kwargs)


def _assert_exact_rf(out, rows):
    np.testing.assert_array_equal(out["medoids"], _vector(rows, "medoid").astype(int))
    np.testing.assert_allclose(out["heights"], _vector(rows, "height"), atol=R_GOLDEN_ATOL, rtol=0)
    np.testing.assert_allclose(out["at"], _vector(rows, "at"), atol=R_GOLDEN_ATOL, rtol=0)
    for key, section in (
        ("index_list", "index"),
        ("weights_list", "weight"),
        ("dist_list", "distance"),
    ):
        expected = _group_lists(rows, section)
        assert len(out[key]) == len(expected)
        for actual, wanted in zip(out[key], expected):
            np.testing.assert_allclose(actual, wanted, rtol=0, atol=R_GOLDEN_ATOL)
    np.testing.assert_allclose(
        [out["R2"], out["Fstat"], out["pvalue"]],
        _vector(rows, "statistic"),
        rtol=0,
        atol=R_GOLDEN_ATOL,
        equal_nan=True,
    )
    assert list(out["sizes"]) == [row["text"] for row in _section(rows, "size")]
    np.testing.assert_allclose(
        list(out["sizes"].values()), _vector(rows, "size"), rtol=0, atol=R_GOLDEN_ATOL
    )
    expected_method = _section(rows, "attribute")[0]["text"].split("=", 1)[1]
    assert out["grp_meth"] == expected_method
    if _section(rows, "membership"):
        np.testing.assert_array_equal(out["mdsk"], _vector(rows, "membership").astype(int))
        np.testing.assert_array_equal(
            out["kmedoid_index"], _vector(rows, "case_medoid").astype(int)
        )
        np.testing.assert_allclose(
            out["kmedoid_dist"], _vector(rows, "case_distance"), atol=R_GOLDEN_ATOL, rtol=0
        )


@pytest.mark.parametrize(
    "case",
    [
        case
        for case, config in _D_CASES.items()
        if config["grp_meth"] != "random" and (config["sort"] != "mds" or config.get("exact_mds"))
    ],
)
def test_dissrf_r2_deterministic_groups_match_live_r(case):
    rows = _case_rows(case)
    assert {row["case"] for row in _r2_rows()} == set(_D_CASES) | set(_S_CASES)
    _assert_exact_rf(_run_d_case(_D_CASES[case]), rows)


@pytest.mark.parametrize(
    "case",
    [
        case
        for case, config in _S_CASES.items()
        if (config["grp_meth"] != "random" or config.get("k") == 1)
        and (config["sort"] != "mds" or config.get("exact_mds"))
    ],
)
def test_seqrf_r2_deterministic_groups_match_live_r(case):
    rows = _case_rows(case)
    out = _run_s_case(_S_CASES[case])
    _assert_exact_rf(out["rf"], rows)
    expected = _group_lists(rows, "representative")
    np.testing.assert_array_equal(out["seqtoplot"], np.vstack(expected).astype(int))


@pytest.mark.parametrize(
    "case,config,runner",
    [
        (case, config, _run_d_case)
        for case, config in _D_CASES.items()
        if config["grp_meth"] == "random"
    ]
    + [
        (case, config, _run_s_case)
        for case, config in _S_CASES.items()
        if config["grp_meth"] == "random" and config.get("k") != 1
    ],
)
def test_random_groups_match_r_invariants_and_cover_each_case_once(case, config, runner):
    rows = _case_rows(case)
    result = runner(config)
    out = result["rf"] if case.startswith("s_") else result
    expected_sizes = sorted(len(group) for group in _group_lists(rows, "index"))
    assert sorted(len(group) for group in out["index_list"]) == expected_sizes
    np.testing.assert_array_equal(
        np.sort(np.concatenate(out["index_list"])), np.arange(len(out["mdsk"]))
    )
    assert all(np.all(group == 1) for group in out["weights_list"])
    assert np.isfinite([out["R2"], out["Fstat"], out["pvalue"]]).all()


@pytest.mark.parametrize(
    "runner,config",
    [
        (_run_d_case, _D_CASES["d_random_numeric_sq0_default"]),
        (_run_s_case, _S_CASES["s_random_numeric_sq0_default"]),
    ],
)
def test_random_state_reproducibility_and_generator_passthrough(runner, config):
    first = runner(config)
    second = runner(config)
    first_rf = first.get("rf", first)
    second_rf = second.get("rf", second)
    np.testing.assert_array_equal(first_rf["mdsk"], second_rf["mdsk"])

    different_result = runner(dict(config, seed=config["seed"] + 1))
    different_rf = different_result.get("rf", different_result)
    assert not np.array_equal(first_rf["mdsk"], different_rf["mdsk"])

    seeded = dict(config)
    seeded.pop("seed")
    seeded["seed"] = np.random.default_rng(9127)
    generated = runner(seeded)
    generated_rf = generated.get("rf", generated)
    assert sorted(np.bincount(generated_rf["mdsk"])) == sorted(np.bincount(first_rf["mdsk"]))


def test_seqrf_passes_generator_object_to_dissrf_unchanged(monkeypatch):
    import traminex.repseq_rf as module

    captured = []
    original = module.dissrf

    def recording_dissrf(*args, random_state=None, **kwargs):
        captured.append(random_state)
        return original(*args, random_state=random_state, **kwargs)

    monkeypatch.setattr(module, "dissrf", recording_dissrf)
    generator = np.random.default_rng(4402)
    config = dict(_S_CASES["s_random_numeric_sq0_default"], seed=generator)
    _run_s_case(config)
    assert captured == [generator]


def test_none_random_state_uses_a_fresh_generator(monkeypatch):
    calls = []
    real_default_rng = np.random.default_rng

    def recording_default_rng(seed=None):
        calls.append(seed)
        return real_default_rng(123 + len(calls))

    monkeypatch.setattr(np.random, "default_rng", recording_default_rng)
    config = dict(_D_CASES["d_random_numeric_sq0_default"])
    config.pop("seed")
    _run_d_case(config)
    _run_d_case(config)
    assert calls == [None, None]


def test_sortv_none_is_original_order_while_default_is_mds():
    D, _, _ = _inputs()
    null = dissrf(D, k=4, sortv=None, grp_meth="first")
    np.testing.assert_array_equal(np.concatenate(null["index_list"]), np.arange(D.shape[0]))
    implicit = dissrf(D, k=4, grp_meth="first")
    explicit = dissrf(D, k=4, sortv="mds", grp_meth="first")
    np.testing.assert_array_equal(implicit["mdsk"], explicit["mdsk"])


def test_weighted_nonprop_refusals_and_method_validation_match_r_contract():
    assert load_golden_lines("diss_rf_r2_errors") == [
        "dissrf_weighted_first=ERROR:Selected grp.meth 'first' does not apply to weighted data! Use grp.meth='prop'.",
        "dissrf_weighted_random=ERROR:Selected grp.meth 'random' does not apply to weighted data! Use grp.meth='prop'.",
        "seqrf_weighted_first=ERROR:Selected grp.meth 'first' does not apply to weighted data! Use grp.meth='prop'.",
        "seqrf_explicit_weighted_random=ERROR:Selected grp.meth 'random' does not apply to weighted data! Use grp.meth='prop'.",
        "dissrf_invalid_method=ERROR: grp.meth must be one of 'first', 'random', 'prop'",
    ]
    D, sortv, weights = _inputs()
    for method in ("first", "random"):
        with pytest.raises(ValueError, match="does not apply to weighted data"):
            dissrf(D, k=3, sortv=sortv, weights=weights, grp_meth=method)
    with pytest.raises(ValueError, match="grp_meth"):
        dissrf(D, k=3, sortv=sortv, grp_meth="bogus")
