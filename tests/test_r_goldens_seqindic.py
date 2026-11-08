"""R-golden parity for seqindic (T11) vs TraMineR 2.2-14.

seqindic is the composite per-sequence indicator table (re-targeted from the
previous one-hot indicator matrix). Pins the default indicator set
{visited, trans, entr, cplx, turb2n} across the with.missing / ragged paths, plus
a broad basic+diversity+complexity set. Goldens:
scripts/r_goldens/make_goldens_seqindic.R.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from corestats_fixtures import fixture_a, fixture_b, fixture_c
from r_goldens_utils import R_GOLDEN_ATOL, load_golden

from traminex.format import seqdef
from traminex.sequences import STSList
from traminex.stats import seqindic

_DEFAULT = ["visited", "trans", "entr", "cplx", "turb2n"]
_BROAD = [
    "lgth",
    "nonm",
    "dlgth",
    "visited",
    "visitp",
    "recu",
    "trans",
    "transp",
    "entr",
    "cplx",
    "turb",
    "turbn",
    "turb2",
    "turb2n",
]
_DEFAULT_LABELS = ["Visited", "Trans", "Entr", "Cplx", "Turb2n"]
GOLDEN_DIR = Path(__file__).parent / "r_goldens"

_DURATION_CASES = {
    "together_scrambled_w05": {
        "indic": ["volat", "dustd2", "meand2", "nsubs", "dustd", "meand"],
        "w": 0.5,
    },
    "meand": {"indic": "meand"},
    "meand2": {"indic": "meand2"},
    "dustd": {"indic": "dustd"},
    "dustd2": {"indic": "dustd2"},
    "nsubs": {"indic": "nsubs"},
    "volat_w0": {"indic": "volat", "w": 0.0},
    "volat_w05": {"indic": "volat", "w": 0.5},
    "volat_w1": {"indic": "volat", "w": 1.0},
}

_BINARY_ALL = ["ppos", "nvolat", "vpos", "integr"]
_BINARY_CASES = {
    "all_dss_false": (
        fixture_a,
        False,
        _BINARY_ALL,
        {
            "pos_states": ["A", "B"],
            "neg_states": ["C", "D"],
            "dss": False,
            "pow": 2,
            "w": 0.25,
        },
    ),
    "all_dss_true": (
        fixture_a,
        False,
        _BINARY_ALL,
        {
            "pos_states": ["A", "B"],
            "neg_states": ["C", "D"],
            "dss": True,
            "pow": 2,
            "w": 0.25,
        },
    ),
    "ppos_forces_nondss": (
        fixture_a,
        False,
        ["ppos"],
        {"pos_states": ["A", "B"], "neg_states": ["C", "D"], "dss": True},
    ),
    "ppos_unclassified_nan": (
        fixture_a,
        False,
        ["ppos"],
        {"pos_states": ["A"], "neg_states": ["B"]},
    ),
    "ppos_phantom_missing_state": (
        fixture_a,
        True,
        ["ppos"],
        {"pos_states": ["*"], "neg_states": ["A"]},
    ),
    "nvolat_forces_dss": (
        fixture_a,
        False,
        ["nvolat"],
        {"pos_states": ["A", "B"], "neg_states": ["C", "D"], "dss": False},
    ),
    "vpos_dss_false_w0": (
        fixture_a,
        False,
        ["vpos"],
        {"pos_states": ["A", "B"], "neg_states": ["C", "D"], "dss": False, "w": 0},
    ),
    "vpos_dss_true_w1": (
        fixture_a,
        False,
        ["vpos"],
        {"pos_states": ["A", "B"], "neg_states": ["C", "D"], "dss": True, "w": 1},
    ),
    "integr_dss_false_pow2": (
        fixture_a,
        False,
        ["integr"],
        {
            "pos_states": ["A", "B"],
            "neg_states": ["C", "D"],
            "dss": False,
            "pow": 2,
        },
    ),
    "missing_all": (
        fixture_b,
        True,
        _BINARY_ALL,
        {
            "pos_states": ["A", "*"],
            "neg_states": ["B", "C"],
            "dss": False,
            "pow": 2,
            "w": 0.75,
        },
    ),
}

_RANKED_ALL = ["degrad", "bad", "prec", "insec"]
_RANKED_REVERSE_ARGS = {
    "state_order": ["D", "C", "B", "A"],
    "method": "RANK",
    "pow": 2,
    "spow": 1,
    "otto": 0.35,
    "a": 0.8,
    "b": 1.1,
    "unused_argument": "filtered",
}
_RANKED_CASES = {
    "all_default": (fixture_a, False, _RANKED_ALL, None),
    "all_explicit_reverse": (
        fixture_a,
        False,
        _RANKED_ALL,
        _RANKED_REVERSE_ARGS,
    ),
    "all_partial_order": (
        fixture_a,
        False,
        _RANKED_ALL,
        {"state_order": ["B", "C"]},
    ),
    "all_state_equiv": (
        fixture_a,
        False,
        _RANKED_ALL,
        {"state_equiv": [["B", "C"]]},
    ),
    "all_partial_equiv": (
        fixture_a,
        False,
        _RANKED_ALL,
        {"state_order": ["A", "D"], "state_equiv": [["D", "B"]]},
    ),
    "all_explicit_stprec": (
        fixture_a,
        False,
        _RANKED_ALL,
        {"stprec": [0.0, 0.5, 0.5, -1.0]},
    ),
    "all_stprec_state_equiv": (
        fixture_a,
        False,
        _RANKED_ALL,
        {"stprec": [0.0, 0.2, 0.8, 1.0], "state_equiv": [["B", "C"]]},
    ),
    "phantom_missing_all": (
        fixture_a,
        True,
        _RANKED_ALL,
        {"unused_argument": 42},
    ),
    "phantom_missing_equiv": (
        fixture_a,
        True,
        _RANKED_ALL,
        {"state_equiv": [["*", "A"]]},
    ),
    "degrad_default": (fixture_a, False, ["degrad"], None),
    "degrad_single_state_order": (
        fixture_a,
        False,
        ["degrad"],
        {"state_order": ["A"]},
    ),
    "degrad_constant_stprec": (
        fixture_a,
        False,
        ["degrad"],
        {"stprec": [1.0, 1.0, 1.0, 1.0]},
    ),
    "bad_default": (fixture_a, False, ["bad"], None),
    "bad_single_state_order": (
        fixture_a,
        False,
        ["bad"],
        {"state_order": ["A"]},
    ),
    "bad_constant_stprec": (
        fixture_a,
        False,
        ["bad"],
        {"stprec": [1.0, 1.0, 1.0, 1.0]},
    ),
    "prec_default": (fixture_a, False, ["prec"], None),
    "insec_default_spow0": (fixture_a, False, ["insec"], None),
    "insec_explicit_spow1": (
        fixture_a,
        False,
        ["insec"],
        {"pow": 2, "spow": 1},
    ),
    "prec_explicit_correction": (
        fixture_a,
        False,
        ["prec"],
        {"correction": [-0.1, 0.2, 0.0, 0.5]},
    ),
    "insec_explicit_correction": (
        fixture_a,
        False,
        ["insec"],
        {"correction": [-0.1, 0.2, 0.0, 0.5]},
    ),
    "private_degrad_options": (
        fixture_a,
        False,
        ["prec", "insec"],
        {
            "penalized": "POS",
            "method": "TRATEDSS",
            "weight_type": "INV",
            "border_effect": 20,
            "pow": 2,
        },
    ),
    "missing_all": (
        fixture_b,
        True,
        _RANKED_ALL,
        {"state_order": ["*", "C", "B", "A"], "method": "RANK", "pow": 2},
    ),
    "filtered_unknown": (
        fixture_a,
        False,
        _RANKED_ALL,
        {"unused_argument": 42},
    ),
    "prec_type_silently_overwritten": (
        fixture_a,
        False,
        ["prec"],
        {"type": 2},
    ),
}

_GROUP_CASES = {
    "basic": (False, {"indic": "basic"}),
    "diversity": (False, {"indic": "diversity"}),
    "diversity_plus_visited": (False, {"indic": ["diversity", "visited"]}),
    "visited_visitp_phantom_missing": (
        True,
        {"indic": ["visited", "visitp"]},
    ),
    "complexity": (False, {"indic": "complexity"}),
    "binary": (
        False,
        {
            "indic": "binary",
            "ipos_args": {
                "pos_states": ["A", "B"],
                "neg_states": ["C", "D"],
                "dss": False,
                "pow": 2,
                "w": 0.25,
            },
        },
    ),
    "ranked": (False, {"indic": "ranked", "prec_args": {"method": "RANK"}}),
    "all_only": (False, {"indic": "all"}),
    "aliases": (
        False,
        {
            "indic": ["prec2", "inpos", "visit"],
            "ipos_args": {
                "pos_states": ["A", "B"],
                "neg_states": ["C", "D"],
                "dss": False,
                "pow": 2,
                "w": 0.25,
            },
            "prec_args": {"method": "RANK"},
        },
    ),
    "mixed_scrambled_duplicates": (
        False,
        {"indic": ["volat", "meand", "trans", "meand", "cplx", "dustd2"]},
    ),
    "direct_visitp": (False, {"indic": "visitp"}),
}


def _check(x, tag: str, with_missing: bool, indic):
    df = seqindic(x, indic=indic, with_missing=with_missing)
    np.testing.assert_allclose(
        df.to_numpy(), load_golden(f"seqindic_{tag}"), rtol=0, atol=R_GOLDEN_ATOL
    )


def _check_case_table(
    filename: str,
    case: str,
    x,
    *,
    with_missing: bool,
    kwargs: dict,
) -> None:
    table = pd.read_csv(GOLDEN_DIR / filename, keep_default_na=False)
    rows = table.loc[table["case"] == case]
    assert not rows.empty
    row_ids = rows["row"].drop_duplicates().tolist()
    columns = rows["column"].drop_duplicates().tolist()
    assert row_ids == list(range(1, x.nseq + 1))
    assert len(rows) == len(row_ids) * len(columns)
    assert not rows[["row", "column"]].duplicated().any()
    actual = seqindic(x, with_missing=with_missing, **kwargs)
    if columns == [""]:
        assert actual.shape == (x.nseq, 0)
        assert actual.columns.tolist() == []
        return
    expected = (
        rows.pivot(index="row", columns="column", values="value")
        .reindex(index=row_ids, columns=columns)
        .apply(pd.to_numeric, errors="coerce")
        .to_numpy(dtype=float)
    )
    assert actual.columns.tolist() == columns
    np.testing.assert_allclose(actual.to_numpy(), expected, rtol=0, atol=R_GOLDEN_ATOL)


@pytest.mark.parametrize(
    ("filename", "expected_cases"),
    [
        ("seqindic_A_duration.csv", set(_DURATION_CASES)),
        ("seqindic_C_duration_wm1.csv", set(_DURATION_CASES)),
        ("seqindic_binary.csv", set(_BINARY_CASES) | {"overwritten_binary_args"}),
        (
            "seqindic_ranked.csv",
            set(_RANKED_CASES) | {"overwritten_ranked_args", "ipos_type_source_bug"},
        ),
        ("seqindic_groups.csv", set(_GROUP_CASES)),
    ],
)
def test_seqindic_golden_case_manifests_are_exact(filename: str, expected_cases: set[str]) -> None:
    table = pd.read_csv(GOLDEN_DIR / filename, keep_default_na=False)
    assert set(table["case"]) == expected_cases


def test_seqindic_default_a():
    df = seqindic(fixture_a(), indic=_DEFAULT)
    assert list(df.columns) == _DEFAULT_LABELS
    _check(fixture_a(), "A", False, _DEFAULT)


def test_seqindic_default_b_with_and_without_missing():
    _check(fixture_b(), "B_wm0", False, _DEFAULT)
    _check(fixture_b(), "B_wm1", True, _DEFAULT)


def test_seqindic_default_c_ragged_with_and_without_missing():
    _check(fixture_c(), "C_wm0", False, _DEFAULT)
    _check(fixture_c(), "C_wm1", True, _DEFAULT)


def test_seqindic_broad_a():
    _check(fixture_a(), "A_broad", False, _BROAD)


@pytest.mark.parametrize("case", list(_DURATION_CASES))
@pytest.mark.parametrize(
    ("filename", "builder", "with_missing"),
    [
        ("seqindic_A_duration.csv", fixture_a, False),
        ("seqindic_C_duration_wm1.csv", fixture_c, True),
    ],
)
def test_seqindic_duration_subsequence_and_volatility_match_r(
    case: str,
    filename: str,
    builder,
    with_missing: bool,
):
    _check_case_table(
        filename,
        case,
        builder(),
        with_missing=with_missing,
        kwargs=_DURATION_CASES[case],
    )


@pytest.mark.parametrize("case", list(_BINARY_CASES))
def test_seqindic_binary_indicators_match_r(case: str):
    builder, with_missing, indic, ipos_args = _BINARY_CASES[case]
    _check_case_table(
        "seqindic_binary.csv",
        case,
        builder(),
        with_missing=with_missing,
        kwargs={"indic": indic, "ipos_args": ipos_args},
    )


def test_seqindic_ipos_args_accept_r_spellings_deterministically():
    _check_case_table(
        "seqindic_binary.csv",
        "all_dss_false",
        fixture_a(),
        with_missing=False,
        kwargs={
            "indic": _BINARY_ALL,
            "ipos_args": {
                "pos.states": ["A", "B"],
                "neg.states": ["C", "D"],
                "dss": False,
                "pow": 2,
                "w": 0.25,
            },
        },
    )


def test_seqindic_ipos_overrides_warn_and_do_not_mutate_input():
    replacement = fixture_b()
    ipos_args = {
        "x": replacement,
        "include_missing": True,
        "index": "integr",
        "pos_states": ["A", "B"],
        "neg_states": ["C", "D"],
        "dss": False,
    }
    original_keys = list(ipos_args)
    original_pos = list(ipos_args["pos_states"])
    original_neg = list(ipos_args["neg_states"])
    with pytest.warns(UserWarning) as caught:
        _check_case_table(
            "seqindic_binary.csv",
            "overwritten_binary_args",
            fixture_a(),
            with_missing=False,
            kwargs={"indic": _BINARY_ALL, "ipos_args": ipos_args},
        )
    messages = [str(item.message) for item in caught]
    expected = pd.read_csv(GOLDEN_DIR / "seqindic_binary_errors.txt", sep="|")
    expected = expected.loc[expected["case"] == "overwritten_binary_args", "message"]
    assert len(messages) == len(expected) == 3
    for token in ("seqdata", "with.missing", "index"):
        assert any(token in message for message in messages)
    assert list(ipos_args) == original_keys
    assert ipos_args["x"] is replacement
    assert ipos_args["include_missing"] is True
    assert ipos_args["index"] == "integr"
    assert ipos_args["pos_states"] == original_pos
    assert ipos_args["neg_states"] == original_neg


def test_seqindic_binary_oracle_refusals():
    conditions = pd.read_csv(GOLDEN_DIR / "seqindic_binary_errors.txt", sep="|")
    errors = conditions.loc[conditions["kind"] == "error"]
    assert set(errors["case"]) == {
        "empty_ipos_args",
        "invalid_indicator",
        "invalid_pos_state",
        "duplicate_pos_state",
    }
    with pytest.raises(ValueError, match="non.empty ipos_args|non empty ipos_args"):
        seqindic(fixture_a(), indic=["ppos"])
    with pytest.raises(ValueError, match="invalid.*indic"):
        seqindic(fixture_a(), indic=["not-an-indicator"])
    with pytest.raises(ValueError, match="invalid.*pos_states|invalid.*pos.states"):
        seqindic(fixture_a(), indic="ppos", ipos_args={"pos_states": ["Z"]})
    with pytest.raises(ValueError, match="Multiple occurrences.*pos_states|duplicate"):
        seqindic(fixture_a(), indic="ppos", ipos_args={"pos_states": ["A", "A"]})


def test_seqindic_binary_accepts_scalar_multicharacter_state_names():
    x = seqdef(
        np.array(
            [["Good", "Good", "Bad"], ["Bad", "Other", "Good"]],
            dtype=object,
        ),
        alphabet=["Good", "Bad", "Other"],
    )
    scalar = seqindic(
        x,
        indic="ppos",
        ipos_args={"pos.states": "Good", "neg.states": "Bad"},
    )
    listed = seqindic(
        x,
        indic="ppos",
        ipos_args={"pos_states": ["Good"], "neg_states": ["Bad"]},
    )
    np.testing.assert_allclose(scalar, listed, rtol=0, atol=R_GOLDEN_ATOL)


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*No\ prec\.args\ list\.\ Alphabet\ order\ used\ as\ state\ order\.:UserWarning"
)
@pytest.mark.parametrize("case", list(_RANKED_CASES))
def test_seqindic_ranked_indicators_match_r(case: str):
    builder, with_missing, indic, prec_args = _RANKED_CASES[case]
    _check_case_table(
        "seqindic_ranked.csv",
        case,
        builder(),
        with_missing=with_missing,
        kwargs={"indic": indic, "prec_args": prec_args},
    )


def test_seqindic_prec_args_accept_r_spellings_deterministically():
    _check_case_table(
        "seqindic_ranked.csv",
        "private_degrad_options",
        fixture_a(),
        with_missing=False,
        kwargs={
            "indic": ["prec", "insec"],
            "prec_args": {
                "penalized": "POS",
                "method": "TRATEDSS",
                "weight.type": "INV",
                "border.effect": 20,
                "pow": 2,
            },
        },
    )


@pytest.mark.parametrize("case", list(_GROUP_CASES))
def test_seqindic_groups_aliases_order_and_quirks_match_r(case: str):
    with_missing, kwargs = _GROUP_CASES[case]
    _check_case_table(
        "seqindic_groups.csv",
        case,
        fixture_a(),
        with_missing=with_missing,
        kwargs=kwargs,
    )


def test_seqindic_ranked_overrides_warn_and_do_not_mutate_input():
    replacement = fixture_b()
    prec_args = {
        "x": replacement,
        "with_missing": True,
        "method": "RANK",
    }
    original = dict(prec_args)
    with pytest.warns(UserWarning) as caught:
        _check_case_table(
            "seqindic_ranked.csv",
            "overwritten_ranked_args",
            fixture_a(),
            with_missing=False,
            kwargs={"indic": _RANKED_ALL, "prec_args": prec_args},
        )
    messages = [str(item.message) for item in caught]
    expected = pd.read_csv(GOLDEN_DIR / "seqindic_ranked_warnings.txt", sep="|")
    expected = expected.loc[expected["case"] == "overwritten_ranked_args", "message"]
    assert len(messages) == len(expected) == 2
    assert any("seqdata" in message for message in messages)
    assert any("with.missing" in message for message in messages)
    assert prec_args == original
    assert prec_args["x"] is replacement


def test_seqindic_empty_prec_args_and_ipos_type_source_bug_match_r():
    conditions = pd.read_csv(GOLDEN_DIR / "seqindic_ranked_warnings.txt", sep="|")
    with pytest.warns(UserWarning, match="No prec.args"):
        _check_case_table(
            "seqindic_ranked.csv",
            "all_default",
            fixture_a(),
            with_missing=False,
            kwargs={"indic": _RANKED_ALL},
        )

    with pytest.warns(UserWarning, match="type argument.*overwritten"):
        _check_case_table(
            "seqindic_ranked.csv",
            "ipos_type_source_bug",
            fixture_a(),
            with_missing=False,
            kwargs={
                "indic": "prec",
                "ipos_args": {"type": 2},
                "prec_args": {"method": "TRATEDSS"},
            },
        )
    warnings = conditions.loc[conditions["kind"] == "warning"]
    errors = conditions.loc[conditions["kind"] == "error"]
    assert set(warnings["case"]) == {
        "empty_prec_args",
        "overwritten_ranked_args",
        "ipos_type_source_bug",
    }
    assert set(errors["case"]) == {
        "prec_spell_integr_duplicate",
        "insec_spell_integr_duplicate",
        "prec_single_state_order",
        "insec_single_state_order",
        "prec_constant_stprec",
        "insec_constant_stprec",
    }
    with pytest.raises(TypeError, match="spell.integr.*multiple"):
        seqindic(
            fixture_a(),
            indic="prec",
            prec_args={"spell.integr": False},
        )
    with pytest.raises(TypeError, match="spell.integr.*multiple"):
        seqindic(
            fixture_a(),
            indic="insec",
            prec_args={"spell_integr": True},
        )
    for indic in ("prec", "insec"):
        with pytest.raises(ValueError, match="NA/NaN argument"):
            seqindic(
                fixture_a(),
                indic=indic,
                prec_args={"state_order": ["A"]},
            )
        with pytest.raises(ValueError, match="NA/NaN argument"):
            seqindic(
                fixture_a(),
                indic=indic,
                prec_args={"stprec": [1.0, 1.0, 1.0, 1.0]},
            )


def test_seqindic_visitp_reuses_loose_object_alphabet_boundary():
    x = STSList(
        np.array([[0, 1], [1, 1]], dtype=int),
        alphabet=["A", "*"],
        nr="*",
    )
    actual = seqindic(x, indic=["visited", "visitp"], with_missing=True)
    np.testing.assert_allclose(actual["Visited"], [2.0, 1.0], rtol=0, atol=R_GOLDEN_ATOL)
    np.testing.assert_allclose(actual["Visitp"], [1.0, 0.5], rtol=0, atol=R_GOLDEN_ATOL)


_ROWNAMES = pd.read_csv(GOLDEN_DIR / "seqindic_rownames.csv", keep_default_na=False, dtype=str)
_IDS_DATA = [["A", "A", "B"], ["A", "B", "B"], ["B", "B", "A"]]
_ROWNAME_OBJECTS = {
    "character_ids": lambda: seqdef(_IDS_DATA, alphabet=["A", "B"], id=["p1", "p2", "p3"]),
    "default_ids": lambda: seqdef(_IDS_DATA, alphabet=["A", "B"]),
    "numeric_ids": lambda: seqdef(_IDS_DATA, alphabet=["A", "B"], id=[10, 20, 30]),
    # R's x[c(3, 1), ]: the rows keep their names.
    "reordered_subset": lambda: seqdef(
        [_IDS_DATA[2], _IDS_DATA[0]], alphabet=["A", "B"], id=["p3", "p1"]
    ),
}


@pytest.mark.parametrize("case", sorted(_ROWNAME_OBJECTS))
def test_seqindic_rows_are_named_by_the_sequence_ids(case: str) -> None:
    # TX-AUD-20261010-093: R names the rows by rownames(seqdata), "1".."n" without ids.
    expected = _ROWNAMES.loc[_ROWNAMES["case"] == case].sort_values("row")["rowname"].tolist()
    assert expected
    actual = seqindic(_ROWNAME_OBJECTS[case](), indic=["visited"])
    assert [str(value) for value in actual.index] == expected
