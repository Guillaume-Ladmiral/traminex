"""R-golden tests for dissdomassoc (RP4 T14).

Anchors traminex's weighted inter-domain Pearson/Spearman correlation + R2 +
Cronbach to live TraMineR 2.2-14 ``dissdomassoc`` (see
scripts/r_goldens/make_goldens_dissdomassoc.R), on the shared 3-channel
multidomain fixture (scripts/r_goldens/make_goldens_md_inputs.R).
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest
from r_goldens_utils import R_GOLDEN_ATOL, load_golden, load_golden_lines

from traminex import dissdomassoc
from traminex.diss import _weighted_midrank

R2_RANKS = Path(__file__).with_name("r_goldens") / "dissdomassoc_r2_rank.csv"


def _domain_dists():
    return [load_golden(f"md_input_d{k}") for k in (1, 2, 3)]


@pytest.mark.parametrize("tag,weighted", [("uw", False), ("w", True)])
def test_dissdomassoc_matches_live_r(tag, weighted):
    dlist = _domain_dists()
    weights = load_golden("md_input_weights").ravel() if weighted else None
    res = dissdomassoc(dlist, what=("pearson", "spearman", "R2", "cronbach"), weights=weights)

    np.testing.assert_allclose(
        res["Pearson"], load_golden(f"dissdomassoc_{tag}_pearson"), rtol=0, atol=R_GOLDEN_ATOL
    )
    np.testing.assert_allclose(
        res["p.pearson"], load_golden(f"dissdomassoc_{tag}_ppearson"), rtol=0, atol=R_GOLDEN_ATOL
    )
    np.testing.assert_allclose(
        res["Pearson.Rsquare"],
        load_golden(f"dissdomassoc_{tag}_r2").ravel(),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )
    np.testing.assert_allclose(
        res["Spearman"], load_golden(f"dissdomassoc_{tag}_spearman"), rtol=0, atol=R_GOLDEN_ATOL
    )
    assert res["Cronbach"] == pytest.approx(
        float(load_golden(f"dissdomassoc_{tag}_cronbach").ravel()[0]), abs=R_GOLDEN_ATOL
    )


def _r2_rows():
    with R2_RANKS.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _rows(case, section, domain=None):
    return [
        row
        for row in _r2_rows()
        if row["case"] == case
        and row["section"] == section
        and (domain is None or row["domain"] == domain)
    ]


def _matrix(case, section):
    rows = _rows(case, section)
    result = np.full(
        (max(int(row["row"]) for row in rows), max(int(row["col"]) for row in rows)),
        np.nan,
    )
    for row in rows:
        result[int(row["row"]) - 1, int(row["col"]) - 1] = float(row["value"])
    return result


def _matrix_from_lower(values, n=5):
    result = np.zeros((n, n), dtype=float)
    upper_row, upper_col = np.triu_indices(n, 1)
    result[upper_col, upper_row] = values
    return result + result.T


def _synthetic(kind):
    if kind == "unique":
        values = [
            np.arange(1, 11),
            np.asarray([10, 2, 8, 4, 6, 1, 9, 3, 7, 5]),
            np.asarray([3, 9, 1, 8, 2, 10, 4, 7, 5, 6]),
        ]
    elif kind == "ties":
        values = [
            np.asarray([1, 1, 2, 2, 3, 3, 4, 5, 5, 5]),
            np.asarray([2, 1, 2, 1, 3, 4, 3, 4, 4, 4]),
            np.asarray([1, 2, 1, 3, 2, 3, 5, 5, 5, 5]),
        ]
    else:
        values = [np.full(10, value) for value in (1, 2, 3)]
    return [_matrix_from_lower(value) for value in values]


_CASES = {
    "base_unweighted_all": dict(kind="base", what=("all",), dnames=["Home", "Work", "Care"]),
    "base_fractional_weighted_rank_joint": dict(
        kind="base",
        weights="base_fractional",
        what=("pearson", "spearman", "R2"),
        w_rank=True,
        joint=True,
        dnames=["Home", "Work", "Care"],
    ),
    "base_integer_weighted_rank_spearman": dict(
        kind="base", weights="base_integer", what=("spearman",), w_rank=True
    ),
    "unique_fractional_weighted_rank": dict(
        kind="unique", weights="small_fractional", what=("spearman", "R2"), w_rank=True
    ),
    "ties_fractional_weighted_rank": dict(
        kind="ties", weights="small_fractional", what=("pearson", "spearman"), w_rank=True
    ),
    "ties_fractional_plain_rank": dict(
        kind="ties", weights="small_fractional", what=("pearson", "spearman"), w_rank=False
    ),
    "ties_zero_weights": dict(kind="ties", weights="zero", what=("spearman",), w_rank=True),
    "all_equal_rank": dict(
        kind="equal", weights="small_fractional", what=("spearman",), w_rank=True
    ),
    "route_pearson": dict(kind="unique", what=("pearson",)),
    "route_r2": dict(kind="unique", what=("R2",)),
    "route_cronbach": dict(kind="unique", what=("cronbach",)),
    "route_cron_subsets": dict(kind="unique", what=("cron.subsets",)),
}


def _case_inputs(config):
    config = dict(config)
    kind = config.pop("kind")
    distances = _domain_dists() if kind == "base" else _synthetic(kind)
    weight_kind = config.pop("weights", None)
    weights = {
        None: None,
        "base_fractional": load_golden("md_input_weights").ravel(),
        "base_integer": np.resize(np.asarray([1, 2, 3], dtype=float), 12),
        "small_fractional": np.asarray([0.5, 1, 1.5, 2, 3]),
        "zero": np.asarray([0, 1, 0, 2, 3], dtype=float),
    }[weight_kind]
    if config.pop("joint", False):
        config["jointdiss"] = sum(distances)
    config.setdefault("dnames", ["D1", "D2", "D3"])
    config["weights"] = weights
    return distances, config


def _expected_keys(case):
    excluded = {"distance", "ordinary_rank", "weighted_rank", "condition", "result_error"}
    ordered = []
    for row in _r2_rows():
        if row["case"] == case and row["section"] not in excluded:
            if row["section"] not in ordered:
                ordered.append(row["section"])
    return ordered


@pytest.mark.parametrize("case", [case for case in _CASES if case != "all_equal_rank"])
def test_dissdomassoc_r2_all_routes_match_r(case):
    assert {row["case"] for row in _r2_rows()} == set(_CASES)
    distances, kwargs = _case_inputs(_CASES[case])
    out = dissdomassoc(distances, **kwargs)
    assert list(out) == _expected_keys(case)
    for key in out:
        if key == "dnames":
            assert out[key] == [row["text"] for row in _rows(case, "dnames")]
        elif isinstance(out[key], dict):
            expected_rows = _rows(case, key)
            assert list(out[key]) == [row["text"].replace(";", ",") for row in expected_rows]
            np.testing.assert_allclose(
                list(out[key].values()),
                [float(row["value"]) for row in expected_rows],
                atol=R_GOLDEN_ATOL,
                rtol=0,
            )
        else:
            expected = _matrix(case, key).squeeze()
            np.testing.assert_allclose(
                out[key], expected, rtol=0, atol=R_GOLDEN_ATOL, equal_nan=True
            )


@pytest.mark.parametrize(
    "case",
    [
        "base_fractional_weighted_rank_joint",
        "base_integer_weighted_rank_spearman",
        "unique_fractional_weighted_rank",
        "ties_fractional_weighted_rank",
        "ties_zero_weights",
        "all_equal_rank",
    ],
)
def test_weighted_midrank_vectors_match_r_formula(case):
    distances, kwargs = _case_inputs(_CASES[case])
    weights = kwargs["weights"]
    pair_indices = np.triu_indices(len(weights), 1)
    pair_weights = weights[pair_indices[0]] * weights[pair_indices[1]]
    all_distances = distances + ([kwargs["jointdiss"]] if "jointdiss" in kwargs else [])
    domains = list(dict.fromkeys(row["domain"] for row in _rows(case, "weighted_rank")))
    for domain, distance in zip(domains, all_distances):
        values = distance[np.triu_indices(distance.shape[0], 1)]
        expected = np.asarray([float(row["value"]) for row in _rows(case, "weighted_rank", domain)])
        np.testing.assert_allclose(
            _weighted_midrank(values, pair_weights), expected, rtol=0, atol=R_GOLDEN_ATOL
        )


def test_w_rank_changes_only_spearman_construction_and_preserves_key_order():
    distances, kwargs = _case_inputs(_CASES["base_fractional_weighted_rank_joint"])
    weighted = dissdomassoc(distances, **kwargs)
    plain = dissdomassoc(distances, **dict(kwargs, w_rank=False))
    assert list(weighted) == list(plain)
    for key in ("Pearson", "p.pearson", "Pearson.Rsquare", "dnames"):
        if key == "dnames":
            assert weighted[key] == plain[key]
        else:
            np.testing.assert_array_equal(weighted[key], plain[key])
    assert not np.array_equal(weighted["Spearman"], plain["Spearman"])


def test_all_equal_weighted_rank_boundary_matches_r_refusal():
    assert (
        "missing value where TRUE/FALSE needed"
        in _rows("all_equal_rank", "result_error")[0]["text"]
    )
    distances, kwargs = _case_inputs(_CASES["all_equal_rank"])
    with pytest.raises(ValueError, match="finite|constant"):
        dissdomassoc(distances, **kwargs)


def test_dissdomassoc_r2_refusals_include_invalid_ppearson_spelling():
    assert load_golden_lines("dissdomassoc_r2_errors") == [
        "invalid_Ppearson=ERROR:bad what values, allowed are pearson,spearman,R2,cronbach,cron.subsets,all",
        "not_list=ERROR:domdiss should be a list of at least two distance matrices or objects",
        "one_domain=ERROR:domdiss should be a list of at least two distance matrices or objects",
        "bad_what=ERROR:bad what values, allowed are pearson,spearman,R2,cronbach,cron.subsets,all",
        "weights_short=ERROR:length of weights not equal to number of cases!",
        "joint_nonnumeric=ERROR:when not NULL, jointdiss must be a distance matrix or object",
        "nonnumeric_domain=ERROR:non numeric values in distance matrices",
    ]
    with pytest.raises(ValueError, match="what"):
        dissdomassoc(_domain_dists(), what=("Ppearson",))
