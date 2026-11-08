"""RP6b: live-R golden (deterministic part) + structure test for seqpolyads.

Anchors the deterministic ``observed.dist`` / ``mean.dist['Obs']`` of traminex
``seqpolyads`` to TraMineRextras 0.6.9 (golden by
``scripts/r_goldens/make_goldens_seqpolyads.R``). The Monte-Carlo permutation part
(``random.dist`` / ``U`` / ``V``) is stochastic -- R's ``set.seed`` + ``sample.int`` are
not reproducible from numpy -- so it is only structure-tested.
"""

from __future__ import annotations

import csv
import re
from functools import lru_cache
from pathlib import Path

import numpy as np
import pytest
from r_goldens_utils import R_DIFFERENCE_ATOL, R_GOLDEN_ATOL, load_golden

from traminex.format import seqdef
from traminex.stats_p2 import seqpolyads

R2_CASES = Path(__file__).with_name("r_goldens") / "seqpolyads_r2_cases.csv"
R2_ERRORS = Path(__file__).with_name("r_goldens") / "seqpolyads_r2_errors.txt"
R2_MESSAGES = Path(__file__).with_name("r_goldens") / "seqpolyads_r2_messages.txt"

ALPH = ["A", "B", "C", "D"]
S1 = [
    ["A", "A", "B", "C", "C", "C"],
    ["A", "B", "B", "B", "C", "A"],
    ["D", "D", "D", "D", "D", "D"],
    ["A", "B", "A", "B", "A", "B"],
]
S2 = [
    ["B", "B", "A", "A", "D", "D"],
    ["C", "C", "C", "A", "A", "A"],
    ["A", "A", "B", "B", "C", "C"],
    ["D", "C", "B", "A", "D", "C"],
]
S3 = [
    ["C", "A", "C", "A", "C", "A"],
    ["D", "D", "A", "A", "B", "B"],
    ["B", "B", "B", "C", "C", "D"],
    ["A", "B", "C", "D", "A", "B"],
]


def _seqlist():
    return [
        seqdef(np.array(S1, dtype=object), alphabet=ALPH),
        seqdef(np.array(S2, dtype=object), alphabet=ALPH),
        seqdef(np.array(S3, dtype=object), alphabet=ALPH),
    ]


def _weighted_seqlist(count=3, width=6):
    weights = [1, 2, 3, 4]
    matrices = [S1, S2, S3][:count]
    return [
        seqdef(
            np.asarray(matrix, dtype=object)[:, :width],
            alphabet=ALPH,
            weights=weights,
        )
        for matrix in matrices
    ]


def _missing_seqlist():
    first = np.asarray(S1, dtype=object).copy()
    second = np.asarray(S2, dtype=object).copy()
    first[0, 2] = None
    first[3, 5] = None
    second[1, 1] = None
    second[3, 5] = None
    weights = [1, 2, 3, 4]
    return [seqdef(matrix, alphabet=ALPH, weights=weights) for matrix in (first, second)]


def test_seqpolyads_observed_dist():
    # method="OM": traminex seqdist OM == R OM (raw), so observed.dist matches R exactly.
    res = seqpolyads(_seqlist(), a=1, method="OM", subm=2.0, indel=1.0, T=500)
    np.testing.assert_allclose(
        res["observed.dist"], load_golden("seqpolyads_obs").ravel(), rtol=0, atol=R_GOLDEN_ATOL
    )
    np.testing.assert_allclose(
        res["mean.dist"][0],
        load_golden("seqpolyads_meanobs").ravel()[0],
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


def test_seqpolyads_observed_dist_ham():
    # method="HAM" is R's default. Now that traminex's seqdist HAM is R-faithful (raw
    # count, audit P1), seqpolyads(method="HAM") matches R exactly (it diverged by 1/L
    # before, which is why the OM golden above existed as a workaround).
    res = seqpolyads(_seqlist(), a=1, method="HAM", T=500)
    np.testing.assert_allclose(
        res["observed.dist"], load_golden("seqpolyads_obs_ham").ravel(), rtol=0, atol=R_GOLDEN_ATOL
    )
    np.testing.assert_allclose(
        res["mean.dist"][0],
        load_golden("seqpolyads_meanobs_ham").ravel()[0],
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


def test_seqpolyads_monte_carlo_structure():
    T = 300
    res = seqpolyads(_seqlist(), a=1, method="OM", subm=2.0, indel=1.0, T=T, s=123)
    n = 4
    assert res["random.dist"].shape == (T,)
    assert np.all(res["random.dist"] >= 0)
    assert res["U"].shape == (n,)
    assert res["V"].shape == (n,)
    assert np.all((res["V"] >= 0) & (res["V"] <= 1))
    assert set(np.unique(res["V.95"])).issubset({0.0, 1.0})
    assert np.isfinite(res["mean.dist"][1])
    # seeded RNG -> reproducible within Python (not byte-identical to R)
    res2 = seqpolyads(_seqlist(), a=1, method="OM", subm=2.0, indel=1.0, T=T, s=123)
    np.testing.assert_array_equal(res["random.dist"], res2["random.dist"])


# R2 complete controls and a=2 ------------------------------------------------
@lru_cache(maxsize=1)
def _r2_rows():
    with R2_CASES.open(encoding="utf-8", newline="") as stream:
        return tuple(csv.DictReader(stream))


def _values(case, section):
    rows = sorted(
        (row for row in _r2_rows() if row["case"] == case and row["section"] == section),
        key=lambda row: int(row["index"]),
    )
    return np.asarray(
        [np.nan if row["value"] in {"NA", "NaN"} else float(row["value"]) for row in rows],
        dtype=float,
    )


def _run_r2_case(case, *, core=1, show_time=False):
    if case == "a1_ham_two":
        return seqpolyads(_weighted_seqlist(2), a=1, method="HAM", T=12, s=11, core=core)
    if case == "a1_om_three_pair_coreauto":
        return seqpolyads(
            _weighted_seqlist(),
            a=1,
            method="OM",
            subm=2.0,
            pair_weights=[1, 2, 0.5],
            indel=1,
            weighted=False,
            T=12,
            s=22,
            core=core,
        )
    if case in {"a1_weighted_type1", "a1_weighted_type2"}:
        return seqpolyads(
            _weighted_seqlist(),
            a=1,
            method="HAM",
            pair_weights=[2, 1, 3],
            role_weights=[0.2, 0.3, 0.5],
            weighted=True,
            rand_weight_type=1 if case.endswith("type1") else 2,
            T=12,
            s=33,
            core=core,
        )
    if case == "a2_ham_replace_true":
        return seqpolyads(
            _weighted_seqlist(2), a=2, method="HAM", replace=True, T=12, s=44, core=core
        )
    if case == "a2_ham_replace_false":
        return seqpolyads(
            _weighted_seqlist(),
            a=2,
            method="HAM",
            pair_weights=[1, 2, 1],
            replace=False,
            weighted=False,
            T=12,
            s=55,
            core=core,
        )
    if case == "a2_om_weighted":
        return seqpolyads(
            _weighted_seqlist(),
            a=2,
            method="OM",
            subm=2.0,
            indel=1,
            role_weights=[0.5, 0.3, 0.2],
            rand_weight_type=2,
            replace=True,
            weighted=True,
            T=12,
            s=66,
            core=core,
        )
    if case == "a2_missing_void":
        return seqpolyads(
            _missing_seqlist(),
            a=2,
            method="OM",
            subm=2.0,
            indel=1,
            replace=True,
            with_missing=True,
            T=8,
            s=77,
            core=core,
        )
    if case in {"seed_defect_s1", "seed_defect_s999"}:
        return seqpolyads(
            _weighted_seqlist(2),
            a=1,
            method="HAM",
            T=10,
            s=1 if case.endswith("s1") else 999,
            core=core,
        )
    if case == "seed_width5":
        return seqpolyads(_weighted_seqlist(2, width=5), a=1, method="HAM", T=10, s=1, core=core)
    if case == "show_time":
        return seqpolyads(
            _weighted_seqlist(2),
            a=1,
            method="HAM",
            T=4,
            s=88,
            core=core,
            show_time=show_time,
        )
    raise AssertionError(case)


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*Found\ voids\ in\ input\ sequences.\ dropping\ all\ right\ missings\ in\ the\ resampled\ sequences!:UserWarning"
)
@pytest.mark.parametrize(
    "case",
    [
        "a1_ham_two",
        "a1_om_three_pair_coreauto",
        "a1_weighted_type1",
        "a1_weighted_type2",
        "a2_ham_replace_true",
        "a2_ham_replace_false",
        "a2_om_weighted",
        "a2_missing_void",
        "seed_defect_s1",
        "seed_defect_s999",
        "seed_width5",
        "show_time",
    ],
)
def test_seqpolyads_r2_observed_components_match_r(case):
    result = _run_r2_case(case, core="auto" if "coreauto" in case else 1)
    np.testing.assert_allclose(
        result["observed.dist"], _values(case, "observed"), rtol=0, atol=R_GOLDEN_ATOL
    )
    np.testing.assert_allclose(
        result["mean.dist"][0], _values(case, "mean_dist")[0], rtol=0, atol=R_GOLDEN_ATOL
    )
    assert list(result) == ["mean.dist", "U", "U.tp", "V", "V.95", "observed.dist", "random.dist"]
    assert result["random.dist"].shape == _values(case, "random").shape


@pytest.mark.parametrize("case", ["a1_om_three_pair_coreauto", "a2_ham_replace_true"])
def test_seqpolyads_core_is_execution_only_and_rng_plan_is_reproducible(case):
    serial = _run_r2_case(case, core=1)
    automatic = _run_r2_case(case, core="auto")
    repeated = _run_r2_case(case, core=1)
    for key in serial:
        np.testing.assert_array_equal(serial[key], automatic[key])
        np.testing.assert_array_equal(serial[key], repeated[key])


def test_seqpolyads_python_seed_contract_is_intentional_r_divergence():
    # The R oracle freezes its overwritten-argument defect: both supplied seeds
    # produce the same stream because sequence width 6 becomes set.seed(6).
    np.testing.assert_array_equal(
        _values("seed_defect_s1", "random"), _values("seed_defect_s999", "random")
    )
    assert not np.array_equal(_values("seed_defect_s1", "random"), _values("seed_width5", "random"))

    first = _run_r2_case("seed_defect_s1")
    other = _run_r2_case("seed_defect_s999")
    assert not np.array_equal(first["random.dist"], other["random.dist"])
    np.testing.assert_array_equal(
        first["random.dist"], _run_r2_case("seed_defect_s1")["random.dist"]
    )


def test_seqpolyads_a2_replace_and_random_weight_controls_are_effective():
    replacement = _run_r2_case("a2_ham_replace_true")
    permutation = _run_r2_case("a2_ham_replace_false")
    assert np.all(replacement["random.dist"] >= 0)
    assert np.all(permutation["random.dist"] >= 0)

    type1 = _run_r2_case("a1_weighted_type1")
    type2 = _run_r2_case("a1_weighted_type2")
    np.testing.assert_array_equal(type1["random.dist"], type2["random.dist"])
    assert type1["mean.dist"][1] != type2["mean.dist"][1]


def test_seqpolyads_show_time_is_one_line_and_does_not_change_result(capsys):
    quiet = _run_r2_case("show_time", show_time=False)
    assert capsys.readouterr().out == ""
    timed = _run_r2_case("show_time", show_time=True)
    output = capsys.readouterr().out
    assert re.fullmatch(r"seqpolyads elapsed: \d+(?:\.\d+)? seconds\n", output)
    assert float(output.split()[2]) >= 0
    for key in quiet:
        np.testing.assert_array_equal(quiet[key], timed[key])

    messages = R2_MESSAGES.read_text(encoding="utf-8")
    assert "timing_header:user  system elapsed" in messages
    assert "timing_values:<elapsed>" in messages


@lru_cache(maxsize=1)
def _r2_errors():
    lines = R2_ERRORS.read_text(encoding="utf-8").splitlines()[1:]
    return dict(line.split("=", 1) for line in lines)


def test_seqpolyads_r2_control_errors_match_r_categories():
    with pytest.raises(ValueError, match="Bad 'a' value"):
        seqpolyads(_weighted_seqlist(2), a=3, T=2)
    with pytest.raises(ValueError, match='core must be numeric or "auto"'):
        seqpolyads(_weighted_seqlist(2), core="many", T=2)
    with pytest.raises(ValueError, match="role_weights"):
        seqpolyads(_weighted_seqlist(), role_weights=[1, 2], T=2)
    with pytest.raises(ValueError, match="at least two"):
        seqpolyads([_weighted_seqlist(2)[0]], T=2)
    with pytest.raises(TypeError, match="stslist"):
        seqpolyads([_weighted_seqlist(2)[0], np.asarray(S2)], T=2)
    with pytest.raises(ValueError, match="same size and alphabet"):
        seqpolyads([_weighted_seqlist(2)[0], _weighted_seqlist(2, width=5)[1]], T=2)
    assert set(_r2_errors()) == {
        "bad_a",
        "bad_core",
        "bad_role_weights",
        "not_list",
        "not_sts",
        "mismatched_width",
    }


def _four_roles():
    def role(rows):
        return seqdef([list(row) for row in rows], alphabet=["A", "B"])

    return [role(["AAA", "AAA"]), role(["AAB", "AAA"]), role(["ABB", "BBB"]), role(["BBB", "ABA"])]


def test_seqpolyads_pair_weights_follow_the_rd_row_wise_order():
    """Owner decision D3 (TX-AUD-20261019-011): the pair weights follow the Rd's row-wise order.

    Worked by hand. Polyad 1 holds AAA, AAB, ABB, BBB: HAM d12 = 1, d13 = 2, d14 = 3, d23 = 1,
    d24 = 2, d34 = 1; polyad 2 holds AAA, AAA, BBB, ABA: d12 = 0, d13 = 3, d14 = 1, d23 = 3,
    d24 = 1, d34 = 2. The pairs (1,2), (1,3), (1,4), (2,3), (2,4), (3,4) take the weights 1..6 and
    the sum is divided by the 6 pairs: (1 + 4 + 9 + 4 + 10 + 6) / 6 = 34 / 6 and
    (0 + 6 + 3 + 12 + 5 + 12) / 6 = 38 / 6. R's code pairs the weights column-wise ((1,2), (1,3),
    (2,3), (1,4), (2,4), (3,4)) and gives 36 / 6 for both, recorded by the R golden.
    """
    result = seqpolyads(
        _four_roles(), a=1, method="HAM", pair_weights=[1, 2, 3, 4, 5, 6], weighted=False, T=4, s=99
    )
    np.testing.assert_allclose(result["observed.dist"], [34 / 6, 38 / 6], rtol=0, atol=1e-12)
    r_values = _values("four_roles_pair_weights", "observed")
    np.testing.assert_allclose(r_values, [6.0, 6.0], rtol=0, atol=R_GOLDEN_ATOL)
    assert not np.allclose(result["observed.dist"], r_values, rtol=0, atol=R_DIFFERENCE_ATOL)
