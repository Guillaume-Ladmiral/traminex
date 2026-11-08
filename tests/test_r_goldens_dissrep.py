"""R-golden tests for dissrep / seqrep (RP4 T13).

Anchors traminex's representative-set selection AND quality statistics to live
TraMineR 2.2-14 ``dissrep`` (see scripts/r_goldens/make_goldens_dissrep.R), on
the shared deterministic phase2 15x15 distance matrix.
"""

from __future__ import annotations

import csv

import numpy as np
import pytest
from r_goldens_utils import R_GOLDEN_ATOL, R_GOLDENS_DIR, load_golden, load_golden_lines

from traminex import dissrep, seqdef, seqdist
from traminex.cluster import seqrep


def _inputs():
    D = load_golden("phase2_input_D")
    w = load_golden("phase2_input_weights_real").ravel()
    return D, w


# tag -> (criterion, weighted)
_CASES = {
    "density_uw": ("density", False),
    "density_w": ("density", True),
    "dist_uw": ("dist", False),
    "dist_w": ("dist", True),
}

_STAT_COLS = ["na", "na(%)", "nb", "nb(%)", "SD", "MD", "DC", "V", "Q"]


@pytest.mark.parametrize("tag", list(_CASES))
def test_dissrep_matches_live_r(tag):
    crit, weighted = _CASES[tag]
    D, w = _inputs()
    weights = w if weighted else None
    out = dissrep(D, criterion=crit, coverage=0.5, weights=weights)

    idx_golden = load_golden(f"dissrep_{tag}_idx").ravel().astype(int)
    np.testing.assert_array_equal(out["representatives"], idx_golden)

    stats_golden = load_golden(f"dissrep_{tag}_stats")
    stats = out["statistics"][_STAT_COLS].to_numpy(dtype=float)
    np.testing.assert_allclose(stats, stats_golden, rtol=0, atol=R_GOLDEN_ATOL)

    q_golden = float(load_golden(f"dissrep_{tag}_quality").ravel()[0])
    assert out["quality"] == pytest.approx(q_golden, abs=R_GOLDEN_ATOL)


# seqrep (T13c): wraps dissrep on a sequence dissimilarity matrix.
_ALPH = ["1", "2", "3", "4"]
_SEQREP_CASES = {
    "density_uw": ("density", False),
    "density_w": ("density", True),
    "dist_w": ("dist", True),
}


def _seqrep_fixture():
    codes = load_golden("seqrep_input_codes").astype(int)
    w = load_golden("seqrep_input_weights").ravel()
    chars = np.array([[_ALPH[c] for c in row] for row in codes], dtype=object)
    return seqdef(chars, alphabet=_ALPH, weights=w)


def _ragged_seqrep_fixture():
    rows = [
        ["C", "B", "B", "B"],
        ["A", "C"],
        ["C", "B"],
        ["B"],
        ["C", "C", "C", "B"],
        ["A", "A", "B", "B", "C", "A"],
        ["A", "A", "C"],
        ["B", "A", "A"],
        ["B"],
        ["B", "B", "A", "B", "A", "C"],
    ]
    weights = np.asarray([1.0, 2.5, 1.5, 3.0, 1.0, 2.0, 1.0, 4.0, 2.0, 1.0])
    return seqdef(rows, alphabet=["A", "B", "C"], weights=weights)


@pytest.mark.parametrize("tag", list(_SEQREP_CASES))
def test_seqrep_matches_live_r(tag):
    crit, weighted = _SEQREP_CASES[tag]
    x = _seqrep_fixture()
    D = load_golden("seqrep_input_diss")
    out = seqrep(x, diss=D, criterion=crit, coverage=0.25, weighted=weighted)

    np.testing.assert_array_equal(
        out["rep_indices"], load_golden(f"seqrep_{tag}_idx").ravel().astype(int)
    )
    stats = out["statistics"][_STAT_COLS].to_numpy(dtype=float)
    np.testing.assert_allclose(
        stats, load_golden(f"seqrep_{tag}_stats"), rtol=0, atol=R_GOLDEN_ATOL
    )
    q_golden = float(load_golden(f"seqrep_{tag}_quality").ravel()[0])
    assert out["quality"] == pytest.approx(q_golden, abs=R_GOLDEN_ATOL)


_EXPLICIT_TIES = np.asarray([4, 4, 3, 3, 2, 2, 1, 1, 0, 0], dtype=float)
_R2_CASES = {
    "prob_uw": {"criterion": "prob", "weighted": False},
    "prob_w": {"criterion": "prob", "weighted": True},
    "mscore_uw": {"criterion": "mscore", "weighted": False},
    "mscore_w": {"criterion": "mscore", "weighted": True},
    "explicit_mscore_desc_ties": {
        "criterion": "mscore",
        "weighted": True,
        "score": _EXPLICIT_TIES,
        "decreasing": True,
        "coverage": 0.4,
    },
    "explicit_prob_asc_ties_nrep2": {
        "criterion": "prob",
        "weighted": False,
        "score": _EXPLICIT_TIES,
        "decreasing": False,
        "nrep": 2,
    },
    "density_uw": {"criterion": "density", "weighted": False},
    "density_w": {"criterion": "density", "weighted": True},
    "dist_uw": {"criterion": "dist", "weighted": False},
    "dist_w": {"criterion": "dist", "weighted": True},
    "density_nrep2_radius0_dmax12": {
        "criterion": "density",
        "weighted": True,
        "nrep": 2,
        "pradius": 0,
        "dmax": 12,
    },
    "dist_nrep3_radius20_dmax12": {
        "criterion": "dist",
        "weighted": False,
        "nrep": 3,
        "pradius": 0.2,
        "dmax": 12,
    },
}


def _r2_rows():
    with (R_GOLDENS_DIR / "seqrep_r2_cases.csv").open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _section(rows, name):
    return [row for row in rows if row["section"] == name]


def _matrix(rows, section):
    selected = _section(rows, section)
    nrow = max(int(row["row"]) for row in selected)
    ncol = max(int(row["col"]) for row in selected)
    result = np.full((nrow, ncol), np.nan, dtype=float)
    for row in selected:
        result[int(row["row"]) - 1, int(row["col"]) - 1] = float(row["value"])
    return result


@pytest.mark.parametrize("case", list(_R2_CASES))
def test_seqrep_r2_criteria_options_and_attributes_match_r(case):
    all_rows = _r2_rows()
    assert {row["case"] for row in all_rows} == set(_R2_CASES)
    rows = [row for row in all_rows if row["case"] == case]
    x = _seqrep_fixture()
    D = load_golden("seqrep_input_diss")
    out = seqrep(x, diss=D, **_R2_CASES[case])

    expected_indices = _matrix(rows, "rep_index").ravel().astype(int)
    np.testing.assert_array_equal(out["rep_indices"], expected_indices)
    np.testing.assert_array_equal(
        out["representatives"], _matrix(rows, "representative").astype(int)
    )
    np.testing.assert_allclose(
        out["scores"],
        _matrix(rows, "score").ravel(),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )
    np.testing.assert_allclose(
        out["statistics"][_STAT_COLS].to_numpy(dtype=float),
        _matrix(rows, "statistic"),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )
    assert out["quality"] == pytest.approx(float(_matrix(rows, "quality")[0, 0]), abs=R_GOLDEN_ATOL)
    np.testing.assert_array_equal(out["rep_group"], _matrix(rows, "rep_group").ravel().astype(int))
    np.testing.assert_allclose(
        out["distance_matrix"],
        _matrix(rows, "distance"),
        rtol=0,
        atol=R_GOLDEN_ATOL,
        equal_nan=True,
    )

    expected_idx_rep = []
    idx_rows = _section(rows, "idx_rep")
    for group in range(1, len(expected_indices) + 1):
        expected_idx_rep.append(
            np.asarray(
                [int(float(row["value"])) for row in idx_rows if int(row["row"]) == group],
                dtype=int,
            )
        )
    assert len(out["idx_rep"]) == len(expected_idx_rep)
    for actual, expected in zip(out["idx_rep"], expected_idx_rep):
        np.testing.assert_array_equal(actual, expected)

    attributes = _section(rows, "attribute")
    assert out["metadata"]["weights_sum"] == pytest.approx(
        float(attributes[0]["value"]), rel=0, abs=R_GOLDEN_ATOL
    )
    assert out["metadata"]["dmax"] == pytest.approx(
        float(attributes[1]["value"]), rel=0, abs=R_GOLDEN_ATOL
    )
    assert attributes[2]["text"] == f"criterion={_R2_CASES[case]['criterion']}"
    assert out["metadata"]["weighted"] == (attributes[3]["text"] == "weighted=TRUE")


def test_seqrep_r2_oracle_messages_and_errors_are_frozen():
    lines = load_golden_lines("seqrep_r2_errors")
    assert len(lines) == 28
    assert any(
        line.startswith("mscore_uw:messages=") and "criterion still under development" in line
        for line in lines
    )
    assert lines[-4:] == [
        "invalid_criterion=ERROR:Unknown criterion / no score provided",
        "invalid_score_length=ERROR:Score must be a vector of length equal to10",
        "invalid_pradius=ERROR:pradius must be between 0 and 1",
        "invalid_coverage_zero_nrep_null=ERROR:missing value where TRUE/FALSE needed",
    ]

    x = _seqrep_fixture()
    D = load_golden("seqrep_input_diss")
    with pytest.raises(ValueError, match="criterion"):
        seqrep(x, diss=D, criterion="bogus")
    with pytest.raises(ValueError, match="score.*length 10"):
        seqrep(x, diss=D, criterion="prob", score=np.asarray([1, 2]))
    with pytest.raises(ValueError, match="pradius"):
        seqrep(x, diss=D, pradius=1.1)
    with pytest.raises(ValueError):
        seqrep(x, diss=D, coverage=0, nrep=None)


@pytest.mark.parametrize("weighted", [False, True])
def test_seqrep_internal_ragged_distance_preserves_sts_metadata(weighted):
    x = _ragged_seqrep_fixture()
    data_before = x.data.copy()
    weights_before = x.weights.copy()
    aware_distance = seqdist(x, method="LCS")

    internal = seqrep(
        x,
        method="LCS",
        criterion="dist",
        coverage=0.25,
        weighted=weighted,
    )
    precomputed = seqrep(
        x,
        diss=aware_distance,
        criterion="dist",
        coverage=0.25,
        weighted=weighted,
    )

    expected_indices = np.asarray([2, 3]) if not weighted else precomputed["rep_indices"]
    np.testing.assert_array_equal(internal["rep_indices"], expected_indices)
    np.testing.assert_array_equal(internal["rep_indices"], precomputed["rep_indices"])
    np.testing.assert_array_equal(internal["representatives"], precomputed["representatives"])
    np.testing.assert_allclose(internal["scores"], precomputed["scores"], rtol=0, atol=0)
    np.testing.assert_allclose(
        internal["statistics"][_STAT_COLS],
        precomputed["statistics"][_STAT_COLS],
        rtol=0,
        atol=0,
    )
    assert internal["quality"] == precomputed["quality"]
    np.testing.assert_array_equal(x.data, data_before)
    np.testing.assert_array_equal(x.weights, weights_before)
