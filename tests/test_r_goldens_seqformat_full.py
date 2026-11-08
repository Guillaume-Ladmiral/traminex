"""TraMineR 2.2-14 goldens for the complete seqformat conversion contract.

R returns a raw matrix/data frame for ``to="STS"``. Traminex deliberately
re-expresses that table as an ``STSList``: the CSV pins the raw R cells, row
identities, and axis names, while its metadata columns document the Python
integer-encoding policy for R's missing mask.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal
from r_goldens_utils import R_GOLDEN_ATOL

from traminex.seqformat import seqformat
from traminex.sequences import STSList

GOLDEN_DIR = Path(__file__).resolve().parent / "r_goldens"
CASES = pd.read_csv(
    GOLDEN_DIR / "seqformat_input_cases.csv",
    dtype=str,
    keep_default_na=False,
)
OUTPUT_CASES = pd.read_csv(
    GOLDEN_DIR / "seqformat_output_cases.csv",
    dtype=str,
    keep_default_na=False,
)


def _output_schema() -> dict[str, dict[str, str]]:
    fields = (
        "case",
        "source",
        "target",
        "r_class",
        "nrow",
        "ncol",
        "columns",
        "row_names",
        "column_types",
        "factor_levels",
        "compressed",
    )
    out: dict[str, dict[str, str]] = {}
    for line in (GOLDEN_DIR / "seqformat_output_schema.txt").read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        values = line.split("|")
        assert len(values) == len(fields), line
        row = dict(zip(fields, values, strict=True))
        out[row["case"]] = row
    return out


OUTPUT_SCHEMA = _output_schema()


def _conditions(kind: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in (GOLDEN_DIR / "seqformat_input_errors.txt").read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        row_kind, name, message = line.split("|", 2)
        if row_kind == kind:
            out[name] = message.replace(r"\n", "\n")
    return out


ERRORS = _conditions("ERROR")
WARNINGS = _conditions("WARN")


def _output_conditions(kind: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in (GOLDEN_DIR / "seqformat_output_errors.txt").read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        row_kind, name, message = line.split("|", 2)
        if row_kind == kind:
            out[name] = message.replace(r"\n", "\n")
    return out


OUTPUT_ERRORS = _output_conditions("ERROR")
OUTPUT_WARNINGS = _output_conditions("WARN")


def _case(name: str) -> pd.DataFrame:
    case = CASES.loc[CASES["case"] == name].copy()
    assert not case.empty, f"missing golden case {name!r}"
    case["row"] = case["row"].astype(int)
    case["position"] = case["position"].astype(int)
    return case.sort_values(["row", "position"], kind="stable")


def _decoded_missing_mask(x: STSList) -> list[list[str]]:
    special = {code for code in (x.nr_code, x.void_code) if code is not None}
    return [
        ["<NA>" if int(code) in special else x.alphabet[int(code)] for code in row]
        for row in x.data
    ]


def _optional_code(value: str) -> int | None:
    return None if value == "<NA>" else int(value)


def _issue_indexes(value: str) -> list[int]:
    return [] if value == "" else [int(item) for item in value.split(";")]


def assert_normalized_case(
    actual: STSList,
    name: str,
    *,
    check_ids: bool = True,
    check_issues: bool = True,
) -> None:
    expected = _case(name)
    rows = [
        group.sort_values("position", kind="stable")["value"].tolist()
        for _, group in expected.groupby("row", sort=True)
    ]
    assert _decoded_missing_mask(actual) == rows

    first = expected.iloc[0]
    assert actual.alphabet == first["alphabet"].split(";")
    assert actual.nr == first["nr"]
    assert actual.void == first["void"]
    assert actual.nr_code == _optional_code(first["nr_code"])
    assert actual.void_code == _optional_code(first["void_code"])
    assert actual.start == int(first["start"])

    if check_issues and first["issues_present"] == "TRUE":
        assert actual.issues == {
            "truncated": _issue_indexes(first["issues_truncated"]),
            "empty.seq": _issue_indexes(first["issues_empty_seq"]),
            "start.before.birth": _issue_indexes(first["issues_start_before_birth"]),
            "spell.before.birth": _issue_indexes(first["issues_spell_before_birth"]),
        }
    elif check_issues:
        assert actual.issues is None

    columns = (
        expected.loc[expected["row"] == expected["row"].min()]
        .sort_values("position", kind="stable")["column"]
        .tolist()
    )
    assert actual.xtlab == columns

    if check_ids:
        ids = expected.drop_duplicates("row", keep="first")["id"].tolist()
        assert actual.ids is not None
        assert [str(value) for value in actual.ids] == ids


def _sts_extended() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "ignored": ["x", "y"],
            "t1": ["A", "B"],
            "t2": ["A", "@"],
            "t3": ["B", "C"],
            "t4": ["@", "C"],
        },
        index=["case-z", "case-a"],
    )


def _compressed_sts() -> pd.DataFrame:
    return pd.DataFrame(
        {"sequence": ["A;B;@", "C;@"]},
        index=["compressed-1", "compressed-2"],
    )


def _custom_stslist() -> STSList:
    return STSList(
        data=np.array([[0, 2, 1], [1, 1, 2]], dtype=int),
        alphabet=["A", "B", "M", "V"],
        labels=["Alpha", "Beta", "Missing", "Void"],
        cpal=["red", "blue", "grey", "white"],
        xtlab=["T1", "T2", "T3"],
        nr="M",
        void="V",
        nr_code=2,
        void_code=3,
        weights=np.array([0.25, 1.75]),
        ids=["object-b", "object-a"],
        missing_color="grey",
        issues={
            "truncated": [2],
            "empty.seq": [],
            "start.before.birth": [],
            "spell.before.birth": [],
        },
    )


def _spell_process() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "who": ["b", "b", "b", "a", "a"],
            "lo": [1, 3, 6, 2, 4],
            "hi": [3, 4, 8, 3, 7],
            "state": ["A", "B", "C", "C", "@"],
        }
    )


def _spell_calendar() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "state": ["A", "B", "C"],
            "who": [10, 10, 20],
            "hi": [2002, 2004, 2003],
            "lo": [2001, 2004, 2002],
        }
    )


def _spell_years() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "who": ["p2", "p2", "p1"],
            "begin": [2002, 2004, 1999],
            "end": [2003, 2005, 2001],
            "state": ["A", "B", "C"],
        }
    )


def _births() -> pd.DataFrame:
    return pd.DataFrame({"person": ["p1", "p2"], "birth": [1998, 2001]})


def _spell_ages() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "who": ["p2", "p2", "p1"],
            "begin": [1, 3, 1],
            "end": [2, 3, 2],
            "state": ["A", "B", "C"],
        }
    )


def _spell_birth_edges() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "who": ["crossing", "before", "normal"],
            "begin": [1999, 1998, 2001],
            "end": [2001, 1999, 2002],
            "state": ["X", "Y", "Z"],
        }
    )


def _birth_edges() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "person": ["crossing", "before", "normal"],
            "birth": [2000, 2000, 2000],
        }
    )


def test_current_wide_long_routes_remain_exact():
    wide = pd.DataFrame({"id": ["p2", "p1"], "left": ["A", "B"], "right": ["C", "D"]})
    expected_long = pd.DataFrame(
        {
            "id": ["p1", "p1", "p2", "p2"],
            "time": [1, 2, 1, 2],
            "state": ["B", "D", "A", "C"],
        }
    )
    long = seqformat(wide, from_="wide", to="long", id_col="id")
    assert_frame_equal(long, expected_long, check_exact=False, rtol=0, atol=R_GOLDEN_ATOL)

    expected_wide = pd.DataFrame({"id": ["p1", "p2"], 1: ["B", "A"], 2: ["D", "C"]})
    back = seqformat(
        long,
        from_="long",
        to="wide",
        id_col="id",
        time_col="time",
        state_col="state",
    )
    assert_frame_equal(back, expected_wide, check_exact=False, rtol=0, atol=R_GOLDEN_ATOL)


def test_wide_long_validation_and_copy_boundaries_are_explicit():
    frame = pd.DataFrame({"id": [1], "state": ["A"]})
    copied = seqformat(frame, from_="wide", to="wide")
    assert_frame_equal(copied, frame, check_exact=False, rtol=0, atol=R_GOLDEN_ATOL)
    assert copied is not frame

    with pytest.raises(ValueError, match="id_col, time_col, and state_col"):
        seqformat(frame, from_="long", to="wide")
    with pytest.raises(ValueError, match="Unsupported format conversion"):
        seqformat(frame, from_="long", to="longer")


def test_sts_extended_and_compressed_inputs_match_r():
    extended = seqformat(
        _sts_extended(),
        var=["t1", "t2", "t3", "t4"],
        from_="STS",
        to="STS",
        missing="@",
    )
    assert_normalized_case(extended, "sts_extended")

    right_na = seqformat(
        _compressed_sts(),
        from_="STS",
        to="STS",
        stsep=";",
        missing="@",
        right=None,
    )
    right_del = seqformat(
        _compressed_sts(),
        from_="STS",
        to="STS",
        stsep=";",
        missing="@",
        right="DEL",
    )
    assert_normalized_case(right_na, "sts_compressed_custom_right_na")
    assert_normalized_case(right_del, "sts_compressed_custom_right_del")
    np.testing.assert_array_equal(right_na.data, right_del.data)
    assert right_na.alphabet == right_del.alphabet

    auto = seqformat(
        pd.DataFrame(
            {"sequence": ["A-B-B", "C-@-A"]},
            index=["auto-1", "auto-2"],
        ),
        from_="STS",
        to="STS",
        missing="@",
    )
    assert_normalized_case(auto, "sts_compressed_auto")

    empty = seqformat(
        pd.DataFrame(
            {"sequence": ["A-B", ""]},
            index=["nonempty", "empty"],
        ),
        from_="STS",
        to="STS",
        missing="@",
    )
    assert_normalized_case(empty, "sts_compressed_empty_row")

    trailing = seqformat(
        pd.DataFrame(
            {"sequence": ["A-", "-", "--", "A--B"]},
            index=["a-trailing", "one-separator", "two-separators", "interior"],
        ),
        from_="STS",
        to="STS",
        missing="@",
    )
    assert_normalized_case(trailing, "sts_compressed_trailing_separator")

    partial_missing_with_separator = seqformat(
        pd.DataFrame(
            {"sequence": ["A-B", None]},
            index=["observed", "missing"],
        ),
        from_="STS",
        to="STS",
    )
    assert partial_missing_with_separator.data.shape == (2, 2)
    assert _decoded_missing_mask(partial_missing_with_separator) == [
        ["A", "B"],
        ["<NA>", "<NA>"],
    ]

    marker_states = seqformat(
        pd.DataFrame(
            [["A", "*", "%"]],
            index=["marker-states"],
            columns=["t1", "t2", "t3"],
        ),
        from_="STS",
        to="STS",
        missing="@",
    )
    assert_normalized_case(marker_states, "sts_marker_states")
    assert marker_states.nr_code is None
    assert marker_states.void_code is None

    numeric = seqformat(
        pd.DataFrame(
            [[1.0, 2.0], [2.0, 1.0]],
            index=["numeric-1", "numeric-2"],
            columns=["t1", "t2"],
        ),
        from_="STS",
        to="STS",
        missing="1",
    )
    assert_normalized_case(numeric, "sts_numeric_missing_token")

    rendered = seqformat(
        pd.DataFrame(
            [[1e5, 100001.0, 123456.0, np.inf, -np.inf, 1.234567890123456]],
            index=["numeric-rendering"],
            columns=["t1", "t2", "t3", "t4", "t5", "t6"],
        ),
        from_="STS",
        to="STS",
    )
    assert_normalized_case(rendered, "sts_numeric_rendering")

    promoted = seqformat(
        pd.DataFrame(
            {
                "integer_state": [100000, 100001],
                "double_state": [1.5, 2.5],
            },
            index=["promoted-1", "promoted-2"],
        ),
        from_="STS",
        to="STS",
    )
    assert_normalized_case(promoted, "sts_numeric_frame_promotion")

    logical = seqformat(
        pd.DataFrame(
            {
                "logical_state": [True, False],
                "numeric_state": [1.0, 2.5],
            },
            index=["logical-1", "logical-2"],
        ),
        from_="STS",
        to="STS",
    )
    assert_normalized_case(logical, "sts_logical_numeric_promotion")

    collision = seqformat(
        pd.DataFrame(
            [["A", "*", "@", "B"], ["%", "A", "@", "@"]],
            index=["collision-1", "collision-2"],
            columns=["t1", "t2", "t3", "t4"],
        ),
        from_="STS",
        to="STS",
        missing="@",
    )
    assert_normalized_case(collision, "sts_marker_collision")
    assert collision.nr == "<NR>"
    assert collision.void == "<VOID>"


def test_stslist_infers_from_and_preserves_source_metadata():
    source = _custom_stslist()
    with pytest.warns(UserWarning, match="nr.*void"):
        result = seqformat(source, to="STS", right=None)

    assert_normalized_case(result, "sts_object_inferred", check_issues=False)
    assert result.labels == source.labels
    assert result.cpal == source.cpal
    np.testing.assert_allclose(result.weights, source.weights, rtol=0, atol=R_GOLDEN_ATOL)
    assert result.missing_color == source.missing_color
    assert result.ids is not source.ids
    assert result.issues == source.issues
    assert result.issues is not source.issues
    assert result.data[0, 1] == result.nr_code
    assert result.data[1, 2] == result.void_code


def test_sps_default_and_custom_inputs_match_r_and_drop_source_ids():
    default = pd.DataFrame(
        {
            "spell1": ["(A,2)", "(C,1)"],
            "spell2": ["(B,1)", "(@,2)"],
        },
        index=["must-not-survive-1", "must-not-survive-2"],
    )
    result = seqformat(default, from_="SPS", to="STS", missing="@")
    assert_normalized_case(result, "sps_default", check_ids=False)
    assert result.ids is None

    custom = pd.DataFrame(
        {"sequence": ["{A;2}~{B;1}", "{C;1}~{@;2}"]},
        index=["also-dropped-1", "also-dropped-2"],
    )
    result = seqformat(
        custom,
        from_="SPS",
        to="STS",
        stsep="~",
        sps_in={"xfix": "{}", "sdsep": ";"},
        missing="@",
    )
    assert_normalized_case(result, "sps_custom_compressed", check_ids=False)
    assert result.ids is None

    edges = pd.DataFrame(
        {
            "spell1": ["(A-B,1)", "(C,0)"],
            "spell2": ["(NA,1)", "(D,2.9)"],
            "spell3": ["(F-NA,1)", "(E,-2)"],
        },
        index=["discarded-1", "discarded-2"],
    )
    result = seqformat(edges, from_="SPS", to="STS")
    assert_normalized_case(result, "sps_coercion_edges", check_ids=False)

    trailing = seqformat(
        pd.DataFrame(
            [["(A-,1)", None]],
            index=["sps-trailing"],
            columns=["episode-1", "episode-2"],
        ),
        from_="SPS",
        to="STS",
    )
    assert_normalized_case(
        trailing,
        "sps_trailing_state_separator",
        check_ids=False,
    )
    assert trailing.ids is None


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*'missing'\ set\ as\ .c\('\*','%'\).,\ the\ 'nr'\ and\ 'void'\ code\ from\ the\ 'data'\ state\ sequence\ object:UserWarning"
)
@pytest.mark.parametrize(
    "source",
    [
        STSList(
            data=np.array([[0, 1, 1]], dtype=int),
            alphabet=["A", "B"],
        ),
        STSList(
            data=np.array(
                [
                    [0, 0, 2, 1],
                    [1, 2, 0, 3],
                ],
                dtype=int,
            ),
            alphabet=["A", "B", "*", "%"],
            nr="*",
            void="%",
            nr_code=2,
            void_code=3,
        ),
    ],
)
def test_canonical_sps_output_is_accepted_directly_on_inverse_route(source):
    emitted = seqformat(source, from_="STS", to="SPS", right=None)
    assert isinstance(emitted, list)

    direct = seqformat(emitted, from_="SPS", to="STS")
    matrix = seqformat(np.asarray(emitted, dtype=object)[:, None], from_="SPS", to="STS")

    assert _decoded_missing_mask(direct) == _decoded_missing_mask(matrix)
    assert direct.alphabet == matrix.alphabet
    assert direct.nr == matrix.nr and direct.void == matrix.void
    assert direct.nr_code == matrix.nr_code and direct.void_code == matrix.void_code
    assert direct.xtlab == matrix.xtlab
    assert direct.ids is matrix.ids is None


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*'missing'\ set\ as\ .c\('\*','%'\).,\ the\ 'nr'\ and\ 'void'\ code\ from\ the\ 'data'\ state\ sequence\ object:UserWarning"
)
def test_canonical_sps_output_roundtrips_with_custom_delimiters():
    source = STSList(
        data=np.array([[0, 0, 1], [1, 0, 0]], dtype=int),
        alphabet=["A", "B"],
    )
    spec = {"xfix": "{}", "sdsep": ";"}
    emitted = seqformat(source, from_="STS", to="SPS", sps_out=spec)

    direct = seqformat(emitted, from_="SPS", to="STS", sps_in=spec)
    matrix = seqformat(
        np.asarray(emitted, dtype=object)[:, None],
        from_="SPS",
        to="STS",
        sps_in=spec,
    )

    assert (
        _decoded_missing_mask(direct)
        == _decoded_missing_mask(matrix)
        == [
            ["A", "A", "B"],
            ["B", "A", "A"],
        ]
    )


def test_sps_input_accepts_an_unambiguous_tuple_of_spell_strings():
    actual = seqformat(("(A,2)-(B,1)", "(B,1)-(A,2)"), from_="SPS", to="STS")

    assert _decoded_missing_mask(actual) == [["A", "A", "B"], ["B", "A", "A"]]


@pytest.mark.parametrize(
    "data",
    [
        [["(A,1)"], ["(B,1)"]],
        ["(A,1)", ["(B,1)"]],
        ["(A,1)", 2],
        [np.array(["(A,1)"]), "(B,1)"],
    ],
)
def test_sps_input_refuses_nested_or_mixed_python_sequences(data):
    with pytest.raises(ValueError, match="one-dimensional.*strings"):
        seqformat(data, from_="SPS", to="STS")


def test_sps_input_does_not_reinterpret_a_bare_string_as_a_sequence_collection():
    with pytest.raises(ValueError, match="data frame or a matrix"):
        seqformat("(A,1)", from_="SPS", to="STS")


def test_spell_process_overlap_gap_missing_and_limit_match_r():
    common = dict(
        from_="SPELL",
        to="STS",
        id="who",
        begin="lo",
        end="hi",
        status="state",
        process=True,
        limit=6,
        fillblanks="_",
        missing="@",
    )
    with pytest.warns(UserWarning, match="truncated at limit 6"):
        later_wins = seqformat(_spell_process(), overwrite=True, right=None, **common)
    assert_normalized_case(later_wins, "spell_process_overwrite")

    with pytest.warns(UserWarning, match="truncated at limit 6"):
        earlier_wins = seqformat(_spell_process(), overwrite=False, **common)
    assert_normalized_case(earlier_wins, "spell_process_keep_earlier")


def test_spell_factor_status_retains_declared_levels_and_appends_fill_state():
    data = pd.DataFrame(
        {
            "who": ["p", "p", "q"],
            "begin": [1, 3, 1],
            "end": [1, 3, 1],
            "state": pd.Categorical(
                ["B", "A", "A"],
                categories=["A", "B", "UNUSED"],
            ),
        }
    )
    result = seqformat(
        data,
        from_="SPELL",
        to="STS",
        id="who",
        begin="begin",
        end="end",
        status="state",
        process=True,
        limit=3,
        fillblanks="G",
    )
    assert_normalized_case(result, "spell_factor_status")


def test_spell_numeric_status_uses_r_character_missing_coercion():
    data = pd.DataFrame(
        {
            "who": ["p", "q"],
            "begin": [1, 1],
            "end": [1, 1],
            "state": [1.0, 2.0],
        }
    )
    result = seqformat(
        data,
        from_="SPELL",
        to="STS",
        id="who",
        begin="begin",
        end="end",
        status="state",
        process=True,
        limit=2,
        missing="1",
    )
    assert_normalized_case(result, "spell_numeric_missing_token")


def test_spell_calendar_axes_selectors_and_bound_recomputation_match_r():
    explicit = seqformat(
        _spell_calendar(),
        from_="SPELL",
        to="STS",
        id=2,
        begin=4,
        end=3,
        status=1,
        process=False,
        tmin=2000,
        tmax=2004,
        fillblanks="_",
    )
    assert_normalized_case(explicit, "spell_calendar")

    relative = seqformat(
        _spell_calendar(),
        var=["who", "lo", "hi", "state"],
        from_="SPELL",
        to="STS",
        id=1,
        begin=2,
        end=3,
        status=4,
        process=False,
        tmin=2000,
        tmax=2004,
        fillblanks="_",
    )
    assert_normalized_case(relative, "spell_calendar_var_relative")

    inferred = seqformat(
        _spell_calendar(),
        from_="SPELL",
        to="STS",
        id=2,
        begin=4,
        end=3,
        status=1,
        process=False,
        tmin=1900,
        tmax=None,
        fillblanks="_",
    )
    assert_normalized_case(inferred, "spell_calendar_inferred_one_bound")


def test_spell_pdata_modes_and_source_fill_quirk_match_r():
    years_common = dict(
        from_="SPELL",
        to="STS",
        id="who",
        begin="begin",
        end="end",
        status="state",
        process=True,
        limit=5,
    )
    year_to_age = seqformat(
        _spell_years(),
        pdata=_births(),
        pvar=["person", "birth"],
        **years_common,
    )
    assert_normalized_case(year_to_age, "spell_year_to_age")

    fractional = seqformat(
        pd.DataFrame(
            {
                "who": ["fractional"],
                "begin": [2001],
                "end": [2002],
                "state": ["A"],
            }
        ),
        pdata=pd.DataFrame({"person": ["fractional"], "birth": [2000.5]}),
        pvar=["person", "birth"],
        **{**years_common, "limit": 4},
    )
    assert_normalized_case(fractional, "spell_fractional_birth")

    births_missing = _births()
    births_missing.loc[births_missing["person"] == "p2", "birth"] = np.nan
    with pytest.warns(UserWarning, match="empty sequences because of missing start time"):
        missing_birth = seqformat(
            _spell_years(),
            pdata=births_missing,
            pvar=["person", "birth"],
            **years_common,
        )
    assert_normalized_case(missing_birth, "spell_year_to_age_missing_birth")

    auto = seqformat(_spell_years(), pdata="auto", **years_common)
    assert_normalized_case(auto, "spell_year_to_age_auto")

    year_gap = pd.DataFrame(
        {
            "who": ["gap", "gap"],
            "begin": [2000, 2002],
            "end": [2000, 2002],
            "state": ["A", "B"],
        }
    )
    gap_birth = pd.DataFrame({"person": ["gap"], "birth": [2000]})
    fill_quirk = seqformat(
        year_gap,
        from_="SPELL",
        to="STS",
        id="who",
        begin="begin",
        end="end",
        status="state",
        process=True,
        pdata=gap_birth,
        pvar=["person", "birth"],
        fillblanks="G",
        limit=3,
    )
    assert_normalized_case(fill_quirk, "spell_year_to_age_fill_quirk")

    age_to_year = seqformat(
        _spell_ages(),
        from_="SPELL",
        to="STS",
        id="who",
        begin="begin",
        end="end",
        status="state",
        process=False,
        pdata=_births(),
        pvar=["person", "birth"],
        tmin=1998,
        tmax=2005,
    )
    assert_normalized_case(age_to_year, "spell_age_to_year")


def test_spell_before_birth_diagnostics_and_issue_indexes_match_r():
    with pytest.warns(UserWarning) as caught:
        result = seqformat(
            _spell_birth_edges(),
            from_="SPELL",
            to="STS",
            id="who",
            begin="begin",
            end="end",
            status="state",
            process=True,
            pdata=_birth_edges(),
            pvar=["person", "birth"],
            limit=3,
        )
    assert [str(item.message) for item in caught] == [
        WARNINGS["spell_birth_edges_1"].removeprefix(" [!!] "),
        WARNINGS["spell_birth_edges_2"].removeprefix(" [!!] "),
    ]
    assert_normalized_case(result, "spell_birth_edges")


def test_spell_no_positive_end_warning_and_issue_indexes_match_r():
    data = pd.DataFrame(
        {
            "who": ["empty"],
            "begin": [np.nan],
            "end": [np.nan],
            "state": ["A"],
        }
    )
    with pytest.warns(UserWarning) as caught:
        result = seqformat(
            data,
            from_="SPELL",
            to="STS",
            id="who",
            begin="begin",
            end="end",
            status="state",
            process=True,
            limit=2,
        )
    assert [str(item.message) for item in caught] == [
        WARNINGS["spell_no_positive_end_1"].removeprefix(" [!!] "),
        WARNINGS["spell_no_positive_end_2"].removeprefix(" [!!] "),
    ]
    assert_normalized_case(result, "spell_no_positive_end")


def _error_calls() -> dict[str, object]:
    base = _spell_process()
    bad_character = base.copy()
    bad_character["lo"] = bad_character["lo"].astype(str)
    bad_float = base.copy()
    bad_float["lo"] = bad_float["lo"].astype(float)
    bad_float.loc[0, "lo"] = 1.5
    bad_start = base.copy()
    bad_start.loc[0, "lo"] = 0
    bad_end = base.copy()
    bad_end.loc[1, "hi"] = 2
    bad_births = _births().copy()
    bad_births["birth"] = bad_births["birth"].astype(str)
    missing_id = base.copy()
    missing_id.loc[1, "who"] = None

    spell_kwargs = dict(
        from_="SPELL",
        to="STS",
        id="who",
        begin="lo",
        end="hi",
        status="state",
    )
    years_kwargs = dict(
        from_="SPELL",
        to="STS",
        id="who",
        begin="begin",
        end="end",
        status="state",
    )
    return {
        "missing_from_non_stslist": lambda: seqformat(
            np.array([["A", "B"]], dtype=object), to="STS"
        ),
        "missing_to": lambda: seqformat(np.array([["A", "B"]], dtype=object), from_="STS"),
        "unknown_from": lambda: seqformat(
            np.array([["A", "B"]], dtype=object), from_="BAD", to="STS"
        ),
        "unknown_to": lambda: seqformat(
            np.array([["A", "B"]], dtype=object), from_="STS", to="BAD"
        ),
        "compressed_separator_required": lambda: seqformat(
            np.array([["ABC"]], dtype=object), from_="STS", to="STS"
        ),
        "compressed_numeric_noncharacter": lambda: seqformat(
            np.array([[1], [2]]), from_="STS", to="STS"
        ),
        "compressed_numeric_negative_separator": lambda: seqformat(
            np.array([[-1]]), from_="STS", to="STS"
        ),
        "compressed_all_missing": lambda: seqformat(
            np.array([[None]], dtype=object), from_="STS", to="STS"
        ),
        "compressed_partial_missing": lambda: seqformat(
            np.array([["A"], [None]], dtype=object), from_="STS", to="STS"
        ),
        "compressed_partial_missing_numeric": lambda: seqformat(
            np.array([[1.0], [np.nan]]), from_="STS", to="STS"
        ),
        "compressed_empty_only": lambda: seqformat(
            np.array([[""]], dtype=object), from_="STS", to="STS"
        ),
        "sts_bad_data": lambda: seqformat(1, from_="STS", to="STS"),
        "sps_rejects_stslist": lambda: seqformat(_custom_stslist(), from_="SPS", to="STS"),
        "spell_requires_dataframe": lambda: seqformat(
            np.array([[1, 1, 2, "A"]], dtype=object), from_="SPELL", to="STS"
        ),
        "bad_missing": lambda: seqformat(
            np.array([["A", "B"]], dtype=object),
            from_="STS",
            to="STS",
            missing=1,
        ),
        "bad_stsep": lambda: seqformat(
            np.array([["A-B"]], dtype=object), from_="STS", to="STS", stsep=1
        ),
        "multichar_stsep": lambda: seqformat(
            np.array([["A--B"]], dtype=object), from_="STS", to="STS", stsep="--"
        ),
        "spell_bad_selector": lambda: seqformat(
            base,
            **{**spell_kwargs, "id": ["who", "lo"]},
        ),
        "spell_missing_id": lambda: seqformat(missing_id, **spell_kwargs),
        "spell_bad_limit": lambda: seqformat(base, limit="six", **spell_kwargs),
        "spell_bad_fillblanks": lambda: seqformat(base, fillblanks=0, **spell_kwargs),
        "spell_empty_fillblanks": lambda: seqformat(base, fillblanks="", **spell_kwargs),
        "spell_multichar_fillblanks": lambda: seqformat(base, fillblanks="GG", **spell_kwargs),
        "spell_character_time": lambda: seqformat(bad_character, **spell_kwargs),
        "spell_noninteger_time": lambda: seqformat(bad_float, **spell_kwargs),
        "spell_start_before_one": lambda: seqformat(bad_start, **spell_kwargs),
        "spell_end_before_start": lambda: seqformat(bad_end, **spell_kwargs),
        "pdata_requires_pvar": lambda: seqformat(_spell_years(), pdata=_births(), **years_kwargs),
        "pdata_birth_must_be_integer": lambda: seqformat(
            _spell_years(), pdata=bad_births, pvar=["person", "birth"], **years_kwargs
        ),
        "pdata_missing_id": lambda: seqformat(
            _spell_years(),
            pdata=_births().loc[lambda frame: frame["person"] == "p1"],
            pvar=["person", "birth"],
            **years_kwargs,
        ),
        "pdata_invalid_string": lambda: seqformat(
            _spell_years(),
            pdata="bogus",
            **years_kwargs,
        ),
        "spell_before_birth_missing_end": lambda: seqformat(
            pd.DataFrame(
                {
                    "who": ["missing-end"],
                    "begin": [1999],
                    "end": [np.nan],
                    "state": ["A"],
                }
            ),
            pdata=pd.DataFrame({"person": ["missing-end"], "birth": [2000]}),
            pvar=["person", "birth"],
            limit=3,
            **years_kwargs,
        ),
        "spell_fill_after_missing_end": lambda: seqformat(
            pd.DataFrame(
                {
                    "who": ["gap", "gap"],
                    "begin": [1, 3],
                    "end": [np.nan, 3],
                    "state": ["A", "B"],
                }
            ),
            fillblanks="G",
            limit=3,
            **years_kwargs,
        ),
        "pdata_duplicate_id": lambda: seqformat(
            _spell_years(),
            pdata=pd.concat(
                [
                    _births(),
                    pd.DataFrame({"person": ["p2"], "birth": [2002]}),
                ],
                ignore_index=True,
            ),
            pvar=["person", "birth"],
            **years_kwargs,
        ),
        "pdata_missing_id_age2year": lambda: seqformat(
            _spell_ages(),
            from_="SPELL",
            to="STS",
            id="who",
            begin="begin",
            end="end",
            status="state",
            process=False,
            pdata=_births().loc[lambda frame: frame["person"] == "p1"],
            pvar=["person", "birth"],
            tmin=1998,
            tmax=2005,
        ),
        "sps_bad_duration": lambda: seqformat(
            np.array([["(A,x)"]], dtype=object),
            from_="SPS",
            to="STS",
            stsep="-",
        ),
        "sps_in_missing_sdsep": lambda: seqformat(
            np.array([["(A,1)", "(B,1)"]], dtype=object),
            from_="SPS",
            to="STS",
            sps_in={"xfix": "()"},
        ),
        "sps_in_empty": lambda: seqformat(
            np.array([["(A,1)", "(B,1)"]], dtype=object),
            from_="SPS",
            to="STS",
            sps_in={},
        ),
        "sps_all_missing": lambda: seqformat(
            np.array([[None, None]], dtype=object), from_="SPS", to="STS"
        ),
        "deprecated_compressed_conflict": lambda: seqformat(
            np.array([["A", "B"]], dtype=object),
            from_="STS",
            to="STS",
            compress=True,
            compressed=False,
        ),
        "deprecated_nr_conflict": lambda: seqformat(
            np.array([["A", "@"]], dtype=object),
            from_="STS",
            to="STS",
            missing="@",
            nr="*",
        ),
    }


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(r"ignore:\s*No\ positive\ value\ in\ 'end'\ column!:UserWarning")
@pytest.mark.filterwarnings(
    r"ignore:\s*skipping\ episode\ \d+\ for\ case\ \d+\ due\ to\ missing\ start/end\ time:UserWarning"
)
@pytest.mark.parametrize("name", sorted(ERRORS))
def test_recorded_r_input_errors(name):
    calls = _error_calls()
    assert name in calls, f"no Python refusal probe for recorded R condition {name!r}"
    with pytest.raises((TypeError, ValueError)) as caught:
        calls[name]()
    assert str(caught.value) == ERRORS[name]


def test_deprecated_aliases_warn_and_retain_input_semantics():
    with pytest.warns(DeprecationWarning, match="compressed is deprecated") as caught:
        compressed = seqformat(
            np.array([["A", "B"]], dtype=object),
            from_="STS",
            to="STS",
            compressed=False,
        )
    assert "compressed is deprecated" in str(caught[0].message)
    assert isinstance(compressed, STSList)

    with pytest.warns(DeprecationWarning, match="nr is deprecated") as caught:
        nr_alias = seqformat(
            np.array([["A", "@"]], dtype=object),
            from_="STS",
            to="STS",
            nr="@",
        )
    assert "nr is deprecated" in str(caught[0].message)
    assert _decoded_missing_mask(nr_alias) == [["A", "<NA>"]]
    assert nr_alias.nr == "*"


def test_pvar_without_dataframe_emits_recorded_warning_and_uses_auto_mode():
    expected = WARNINGS["pvar_ignored_without_dataframe"]
    with pytest.warns(UserWarning) as caught:
        result = seqformat(
            _spell_years(),
            from_="SPELL",
            to="STS",
            id="who",
            begin="begin",
            end="end",
            status="state",
            process=True,
            pdata="auto",
            pvar=["person", "birth"],
            limit=5,
        )
    assert expected.removeprefix(" [!!] ") in str(caught[0].message)
    assert_normalized_case(result, "spell_year_to_age_auto")


def test_condition_artifact_has_expected_alias_and_metadata_messages():
    assert "compressed is deprecated" in WARNINGS["deprecated_compressed"]
    assert "nr is deprecated" in WARNINGS["deprecated_nr"]
    assert "c('M','V')" in WARNINGS["stslist_custom_missing_codes"]
    assert "truncated at limit 6" in WARNINGS["spell_truncated"]
    assert "empty sequences because of missing start time" in WARNINGS["spell_missing_birth"]


# ---------------------------------------------------------------------------
# Task 14: complete three-input by six-output conversion hub.
#
# Python representation contract:
# - STS and DSS targets are STSList objects, including compressed one-column
#   results;
# - SPS targets are one canonical, hyphen-joined spell string per sequence;
# - SRS, TSE, and SPELL targets are DataFrames with a RangeIndex.
# ---------------------------------------------------------------------------


def _output_sts() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "group": ["g2", "g1"],
            "id": ["subject-2", "subject-1"],
            "t1": ["A", "B"],
            "t2": ["A", "C"],
            "t3": ["@", "C"],
            "t4": ["B", "@"],
        },
        index=["p2", "p1"],
    )


def _output_sps() -> pd.DataFrame:
    return pd.DataFrame(
        [
            ["(A,2)", "(@,1)", "(B,1)"],
            ["(B,1)", "(C,2)", "(@,1)"],
        ],
        index=["p2", "p1"],
        columns=["spell1", "spell2", "spell3"],
    )


def _output_spell() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "who": ["p2", "p2", "p2", "p1", "p1", "p1"],
            "begin": [1, 3, 4, 1, 2, 4],
            "end": [2, 3, 4, 1, 3, 4],
            "state": ["A", "@", "B", "B", "C", "@"],
        }
    )


def _output_tevent() -> pd.DataFrame:
    states = ["A", "B", "C", "@"]
    value = pd.DataFrame(None, index=states, columns=states, dtype=object)
    # Base R's strsplit() discards the terminal empty component. These trailing
    # separators therefore still produce exactly one event apiece.
    value.loc["A", "A"] = "startA,"
    value.loc["B", "B"] = "startB"
    value.loc["C", "C"] = "startC"
    value.loc["A", "B"] = "AtoB,"
    value.loc["B", "C"] = "leaveB,enterC"
    value.loc["C", "A"] = "CtoA"
    value.loc["A", "@"] = "rawEnterMissing"
    value.loc["@", "B"] = "rawLeaveMissing"
    value.loc["C", "@"] = "rawEnterMissingFromC"
    return value


def _output_pdata() -> pd.DataFrame:
    return pd.DataFrame({"person": ["birth-p2", "birth-p1"], "birth": [1990, 1980]})


def _output_route_args(source: str, target: str) -> dict[str, object]:
    if source == "STS":
        args: dict[str, object] = {
            "data": _output_sts(),
            "var": ["t1", "t2", "t3", "t4"],
            "from_": "STS",
            "to": target,
            "missing": "@",
        }
    elif source == "SPS":
        args = {
            "data": _output_sps(),
            "from_": "SPS",
            "to": target,
            "missing": "@",
        }
    else:
        args = {
            "data": _output_spell(),
            "from_": "SPELL",
            "to": target,
            "id": "who",
            "begin": "begin",
            "end": "end",
            "status": "state",
            "limit": 4,
            "missing": "@",
        }

    if source == "STS" and target == "STS":
        args["compress"] = True
    if source == "SPS" and target == "DSS":
        args["compress"] = True
    if source == "SPELL" and target in {"DSS", "SPELL"}:
        args["right"] = None
    if source == "SPELL" and target == "SPS":
        args["right"] = "NEUTRAL"

    if target == "SPS":
        args["sps_out"] = {"xfix": "{}", "sdsep": ";"}
    if source == "SPS" and target == "SPS":
        args["sps_out"] = {"xfix": "[", "sdsep": ":"}
    if source == "SPELL" and target == "SPS":
        args["compress"] = True

    if target == "SRS":
        args["nrep"] = 2 if source == "SPS" else 1
    if source == "STS" and target == "SRS":
        args["covar"] = ["group", "id"]

    if target == "TSE":
        args["tevent"] = _output_tevent()
        if source == "STS":
            args["id"] = ["alpha", "beta"]
        if source == "SPS":
            args["id"] = ["sps-z", "sps-a"]

    if target == "SPELL":
        args["with_missing"] = source != "STS"
        if source == "STS":
            args["pdata"] = _output_pdata()
            args["pvar"] = ["person", "birth"]
    return args


OUTPUT_ROUTES = [
    (source, target)
    for source in ("STS", "SPS", "SPELL")
    for target in ("STS", "DSS", "SPS", "SRS", "TSE", "SPELL")
]


def _output_case(name: str) -> pd.DataFrame:
    case = OUTPUT_CASES.loc[OUTPUT_CASES["case"] == name].copy()
    assert not case.empty, f"missing output golden case {name!r}"
    case["row"] = case["row"].astype(int)
    case["position"] = case["position"].astype(int)
    return case.sort_values(["row", "position"], kind="stable")


def _expected_output_rows(name: str) -> list[list[str]]:
    expected = _output_case(name)
    return [
        group.sort_values("position", kind="stable")["value"].tolist()
        for _, group in expected.groupby("row", sort=True)
    ]


def _render_output_cell(value: object) -> str:
    try:
        missing = pd.isna(value)
    except (TypeError, ValueError):
        missing = False
    if isinstance(missing, (bool, np.bool_)) and bool(missing):
        return "<NA>"
    if isinstance(value, (int, np.integer)) and not isinstance(value, (bool, np.bool_)):
        return str(int(value))
    if isinstance(value, (float, np.floating)) and float(value).is_integer():
        return str(int(value))
    return str(value)


def _schema_vector(value: str) -> list[str]:
    return [] if value == "<NONE>" else value.split(";")


def _decoded_output_sts(x: STSList, target: str) -> list[list[str]]:
    rows: list[list[str]] = []
    for row in x.data:
        decoded: list[str] = []
        for value in row:
            code = int(value)
            if x.void_code is not None and code == x.void_code:
                decoded.append("<NA>")
            elif target == "STS" and x.nr_code is not None and code == x.nr_code:
                decoded.append("<NA>")
            else:
                decoded.append(x.alphabet[code])
        rows.append(decoded)
    return rows


def _assert_output_route(actual: object, source: str, target: str) -> None:
    name = f"route_{source.lower()}_{target.lower()}"
    schema = OUTPUT_SCHEMA[name]
    expected_rows = _expected_output_rows(name)
    expected_row_names = _schema_vector(schema["row_names"])
    expected_columns = _schema_vector(schema["columns"])
    python_columns = (
        expected_columns
        if expected_columns
        else [f"[{position}]" for position in range(1, int(schema["ncol"]) + 1)]
    )

    if target in {"STS", "DSS"}:
        assert isinstance(actual, STSList)
        assert actual.nseq == int(schema["nrow"])
        assert actual.seqlen == int(schema["ncol"])
        assert _decoded_output_sts(actual, target) == expected_rows
        assert actual.xtlab == python_columns
        assert actual.start == 1
        if "<NA>" in {value for row in expected_rows for value in row}:
            assert actual.void_code is not None
            assert actual.alphabet[actual.void_code] == "%"
        if schema["compressed"] == "FALSE" and (
            target == "STS" or any("*" in row for row in expected_rows)
        ):
            assert actual.nr_code is not None
            assert actual.alphabet[actual.nr_code] == "*"
        if target == "STS" and source == "SPS" and schema["compressed"] == "FALSE":
            # TraMineR's SPS_to_STS discards source row names.
            assert actual.ids is None
        else:
            assert actual.ids is not None
            assert [str(value) for value in actual.ids] == expected_row_names
        return

    if target == "SPS":
        assert isinstance(actual, list)
        expected_strings = [
            "-".join(value for value in row if value != "<NA>") for row in expected_rows
        ]
        assert actual == expected_strings
        return

    assert isinstance(actual, pd.DataFrame)
    assert actual.shape == (int(schema["nrow"]), int(schema["ncol"]))
    assert [str(column) for column in actual.columns] == python_columns
    actual_rows = [
        [_render_output_cell(value) for value in row]
        for row in actual.itertuples(index=False, name=None)
    ]
    assert actual_rows == expected_rows

    # R's automatic data.frame row.names are represented idiomatically.
    assert isinstance(actual.index, pd.RangeIndex)
    if schema["factor_levels"] != "<NONE>":
        state_levels = schema["factor_levels"].removeprefix("states=").split(",")
        assert actual.attrs.get("statl") == state_levels


def test_output_artifacts_cover_each_route_once():
    expected = {f"route_{source.lower()}_{target.lower()}" for source, target in OUTPUT_ROUTES}
    auxiliary = {"matrix_covar", "stslist_covar"}
    assert set(OUTPUT_SCHEMA) == expected | auxiliary
    assert set(OUTPUT_CASES["case"]) == expected | auxiliary
    for name, schema in OUTPUT_SCHEMA.items():
        rows = OUTPUT_CASES.loc[OUTPUT_CASES["case"] == name]
        assert len(rows) == int(schema["nrow"]) * int(schema["ncol"])


@pytest.mark.parametrize(
    ("source", "target"),
    OUTPUT_ROUTES,
    ids=lambda value: value.lower(),
)
def test_all_eighteen_output_routes_match_r(source: str, target: str):
    actual = seqformat(**_output_route_args(source, target))
    _assert_output_route(actual, source, target)


def _assert_recorded_dataframe_case(actual: pd.DataFrame, name: str) -> None:
    schema = OUTPUT_SCHEMA[name]
    assert actual.shape == (int(schema["nrow"]), int(schema["ncol"]))
    assert [str(column) for column in actual.columns] == _schema_vector(schema["columns"])
    actual_rows = [
        [_render_output_cell(value) for value in row]
        for row in actual.itertuples(index=False, name=None)
    ]
    assert actual_rows == _expected_output_rows(name)


def test_srs_covar_supports_matrix_and_stslist_sources():
    matrix_source = np.array(
        [["g2", "A", "B"], ["g1", "B", "C"]],
        dtype=object,
    )
    matrix_actual = seqformat(
        matrix_source,
        var=[2, 3],
        from_="STS",
        to="SRS",
        covar=[1],
    )
    assert isinstance(matrix_actual, pd.DataFrame)
    _assert_recorded_dataframe_case(matrix_actual, "matrix_covar")

    stslist_source = STSList(
        data=np.array([[0, 1], [1, 2]], dtype=int),
        alphabet=["A", "B", "C"],
        xtlab=["T1", "T2"],
    )
    with pytest.warns(UserWarning, match="nr.*void"):
        stslist_actual = seqformat(
            stslist_source,
            from_="STS",
            to="SRS",
            covar=[1],
        )
    assert isinstance(stslist_actual, pd.DataFrame)
    _assert_recorded_dataframe_case(stslist_actual, "stslist_covar")


def _output_error_calls() -> dict[str, object]:
    base = {
        "data": _output_sts(),
        "var": ["t1", "t2", "t3", "t4"],
        "from_": "STS",
        "missing": "@",
    }
    incomplete_tevent = _output_tevent().loc[["A", "B", "C"], ["A", "B", "C"]]
    return {
        "missing_tevent": lambda: seqformat(**base, to="TSE", id=["a", "b"]),
        "duplicate_tse_ids": lambda: seqformat(
            **base,
            to="TSE",
            id=["a", "a"],
            tevent=_output_tevent(),
        ),
        "empty_tse_ids": lambda: seqformat(
            **base,
            to="TSE",
            id=[],
            tevent=_output_tevent(),
        ),
        "incomplete_tevent": lambda: seqformat(
            **base,
            to="TSE",
            id=["a", "b"],
            tevent=incomplete_tevent,
        ),
        "compress_srs": lambda: seqformat(**base, to="SRS", compress=True),
        "compress_tse": lambda: seqformat(
            **base,
            to="TSE",
            id=["a", "b"],
            tevent=_output_tevent(),
            compress=True,
        ),
        "compress_spell": lambda: seqformat(**base, to="SPELL", compress=True),
        "spell_bad_pdata_length": lambda: seqformat(
            **base,
            to="SPELL",
            pdata=_output_pdata().iloc[[0]],
            pvar=["person", "birth"],
        ),
        "bad_covar": lambda: seqformat(**base, to="SRS", covar=["notthere"]),
        "srs_short_nrep": lambda: seqformat(
            np.array([["A", "B"], ["B", "C"]], dtype=object),
            from_="STS",
            to="SRS",
            nrep=1,
        ),
        "srs_negative_nrep": lambda: seqformat(
            np.array([["A", "B", "C"], ["B", "C", "A"]], dtype=object),
            from_="STS",
            to="SRS",
            nrep=-1,
        ),
    }


OUTPUT_ERROR_PATTERNS = {
    "missing_tevent": r"tevent.*(specified|required)",
    "duplicate_tse_ids": r"id.*unique",
    "empty_tse_ids": r"id.*(invalid|non-empty)",
    "incomplete_tevent": r"(tevent|transition|argument is of length zero)",
    "compress_srs": r"SRS is not a compressible format",
    "compress_tse": r"TSE is not a compressible format",
    "compress_spell": r"SPELL is not a compressible format",
    "spell_bad_pdata_length": r"pdata.*one entry per sequence",
    "bad_covar": r"undefined columns selected",
    "srs_short_nrep": r"invalid 'times' argument",
    "srs_negative_nrep": r"numbers of columns of arguments do not match",
}


@pytest.mark.parametrize("name", sorted(OUTPUT_ERRORS))
def test_recorded_r_output_refusals(name: str):
    calls = _output_error_calls()
    assert set(calls) == set(OUTPUT_ERRORS) == set(OUTPUT_ERROR_PATTERNS)
    with pytest.raises(
        (TypeError, ValueError, IndexError, KeyError), match=OUTPUT_ERROR_PATTERNS[name]
    ):
        calls[name]()


def test_output_validation_rejects_malformed_non_numeric_controls():
    base = {
        "data": _output_sts(),
        "var": ["t1", "t2", "t3", "t4"],
        "from_": "STS",
        "missing": "@",
    }

    with pytest.raises(ValueError, match="'pvar' required"):
        seqformat(**base, to="SPELL", pdata=_output_pdata())
    with pytest.raises(ValueError, match="'sps_out' must be a mapping"):
        seqformat(**base, to="SPS", sps_out="not-a-mapping")
    with pytest.raises(ValueError, match="'tevent' must be a data frame"):
        seqformat(**base, to="TSE", id=["a", "b"], tevent={})


def _warning_messages(caught: pytest.WarningsRecorder) -> list[str]:
    return [str(item.message) for item in caught]


def test_output_covar_warning_is_pinned_and_covar_is_ignored():
    expected = OUTPUT_WARNINGS["covar_ignored_non_srs"].removeprefix(" [!!] ")
    with pytest.warns(UserWarning) as caught:
        actual = seqformat(
            _output_sts(),
            var=["t1", "t2", "t3", "t4"],
            from_="STS",
            to="DSS",
            missing="@",
            covar=["group"],
        )
    assert expected in _warning_messages(caught)
    _assert_output_route(actual, "STS", "DSS")


def test_output_deprecated_aliases_warn_and_keep_route_semantics():
    with pytest.warns(DeprecationWarning) as compressed_caught:
        compressed = seqformat(
            _output_sts(),
            var=["t1", "t2", "t3", "t4"],
            from_="STS",
            to="DSS",
            missing="@",
            compressed=True,
        )
    expected_compressed = OUTPUT_WARNINGS["deprecated_compressed_output"].removeprefix(" [!!] ")
    assert expected_compressed in _warning_messages(compressed_caught)
    assert isinstance(compressed, STSList)
    assert compressed.xtlab == ["Sequence"]
    assert _decoded_output_sts(compressed, "DSS") == [["A-*-B"], ["B-C"]]

    with pytest.warns(DeprecationWarning) as nr_caught:
        nr_result = seqformat(
            _output_sts(),
            var=["t1", "t2", "t3", "t4"],
            from_="STS",
            to="SPS",
            nr="@",
        )
    expected_nr = OUTPUT_WARNINGS["deprecated_nr_output"].removeprefix(" [!!] ")
    assert expected_nr in _warning_messages(nr_caught)
    assert isinstance(nr_result, list)
    assert nr_result == ["(A,2)-(*,1)-(B,1)", "(B,1)-(C,2)"]


def test_tse_omitted_id_warns_twice_and_uses_sequence_indexes():
    expected = [
        OUTPUT_WARNINGS["tse_default_ids_1"].removeprefix(" [!!] "),
        OUTPUT_WARNINGS["tse_default_ids_2"].removeprefix(" [!!] "),
    ]
    with pytest.warns(UserWarning) as caught:
        actual = seqformat(
            _output_sts(),
            var=["t1", "t2", "t3", "t4"],
            from_="STS",
            to="TSE",
            missing="@",
            tevent=_output_tevent(),
        )
    assert _warning_messages(caught) == expected
    assert isinstance(actual, pd.DataFrame)
    assert actual["id"].tolist() == [1, 1, 1, 2, 2, 2, 2]


def test_spell_output_warns_when_pvar_has_no_dataframe_pdata():
    expected = OUTPUT_WARNINGS["spell_pvar_ignored_without_dataframe"].removeprefix(" [!!] ")
    with pytest.warns(UserWarning) as caught:
        actual = seqformat(
            _output_sts(),
            var=["t1", "t2", "t3", "t4"],
            from_="STS",
            to="SPELL",
            missing="@",
            pdata="auto",
            pvar=[1, 2],
        )
    assert expected in _warning_messages(caught)
    assert isinstance(actual, pd.DataFrame)
