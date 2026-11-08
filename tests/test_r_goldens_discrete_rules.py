"""R-golden coverage for TraMineRextras 0.6.9 ``createdatadiscrete``."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDEN_ATOL, R_GOLDENS_DIR

import traminex
from traminex.events import SubSeqElist


def _inputs():
    ids = [101, 102, 103, 104]
    data = pd.DataFrame(
        {
            "IDPERS": ids,
            "A": [20, 22, 19, 30],
            "AST": [1, 0, 1, 0],
            "B": [23, 21, 24, 20],
            "BST": [1, 1, 0, 1],
            "C": [25, 24, 22, 21],
            "CST": [0, 1, 1, 0],
            "sex": [1, 2, 1, 2],
            "cohort": [7, 7, 8, 8],
        }
    )
    agemin = pd.DataFrame({"IDPERS": ids, "AGE": [18, 19, 18, 20]})
    agemax = pd.DataFrame({"IDPERS": ids, "AGE": [25, 24, 24, 26]})
    return ids, data, ["A", "B", "C"], agemin, agemax


def _golden_frame(case_name: str, event: str) -> pd.DataFrame:
    records = pd.read_csv(
        R_GOLDENS_DIR / "discrete_rules_r3_cases.csv",
        dtype=str,
        keep_default_na=False,
    )
    selected = records.loc[
        (records["case"] == case_name)
        & (records["event"] == event)
        & (records["record_type"] == "cell")
    ]
    row_order = list(dict.fromkeys(selected["row"]))
    column_order = list(dict.fromkeys(selected["column"]))
    result = pd.DataFrame(index=row_order, columns=column_order, dtype=float)
    for record in selected.itertuples(index=False):
        result.loc[record.row, record.column] = float(record.value)
    return result.reset_index(drop=True)


@pytest.mark.parametrize(
    "case_name,supvar",
    [("no_supvar", None), ("one_supvar", "sex")],
)
def test_createdatadiscrete_frames_mapping_and_messages_match_r(case_name, supvar, capsys):
    ids, data, variables, agemin, agemax = _inputs()
    snapshots = [frame.copy(deep=True) for frame in (data, agemin, agemax)]

    result = traminex.createdatadiscrete(ids, data, variables, agemin, agemax, supvar=supvar)

    assert list(result) == ["data", "vars", "ids"]
    assert list(result["data"]) == variables
    assert result["vars"] == []
    assert result["ids"] == ids
    for event in variables:
        frame = result["data"][event]
        pd.testing.assert_frame_equal(
            frame, _golden_frame(case_name, event), check_exact=False, rtol=0, atol=R_GOLDEN_ATOL
        )
        assert all(dtype == np.dtype("float64") for dtype in frame.dtypes)
    assert capsys.readouterr().err.splitlines() == [
        "Creating data set for variable A",
        "Creating data set for variable B",
        "Creating data set for variable C",
    ]
    for frame, snapshot in zip((data, agemin, agemax), snapshots, strict=True):
        pd.testing.assert_frame_equal(
            frame, snapshot, check_exact=False, rtol=0, atol=R_GOLDEN_ATOL
        )


def test_multiple_supvars_are_corrected_after_freezing_r_dimension_defect(capsys):
    transcript = (R_GOLDENS_DIR / "discrete_rules_r3_errors.txt").read_text(encoding="utf-8")
    assert "multiple_supvar|simpleError|length of 'dimnames'" in transcript
    ids, data, variables, agemin, agemax = _inputs()

    result = traminex.createdatadiscrete(
        ids,
        data,
        variables,
        agemin,
        agemax,
        supvar=["sex", "cohort"],
    )

    assert list(result["data"]["A"]) == [
        "IDPERS",
        "AGE",
        "A",
        "B",
        "C",
        "sex",
        "cohort",
    ]
    assert set(result["data"]["A"]["sex"]) == {1.0, 2.0}
    assert set(result["data"]["A"]["cohort"]) == {7.0, 8.0}
    assert "Creating data set for variable A" in capsys.readouterr().err


def test_single_variable_is_corrected_and_empty_variables_are_refused(capsys):
    transcript = (R_GOLDENS_DIR / "discrete_rules_r3_errors.txt").read_text(encoding="utf-8")
    assert "single_variable|simpleError|missing value where TRUE/FALSE needed" in transcript
    assert "empty_variables|simpleError|missing value where TRUE/FALSE needed" in transcript
    ids, data, _, agemin, agemax = _inputs()

    result = traminex.createdatadiscrete(ids, data, ["A"], agemin, agemax)
    assert list(result["data"]) == ["A"]
    assert list(result["data"]["A"]) == ["IDPERS", "AGE", "A"]
    assert result["data"]["A"]["A"].sum() == 2
    capsys.readouterr()

    with pytest.raises(ValueError, match="vars must contain at least one"):
        traminex.createdatadiscrete(ids, data, [], agemin, agemax)


def test_ids_without_valid_windows_are_skipped_without_r_rest_reuse(capsys):
    transcript = (R_GOLDENS_DIR / "discrete_rules_r3_errors.txt").read_text(encoding="utf-8")
    assert "missing_window_id|simpleError|missing value where TRUE/FALSE needed" in transcript
    assert "na_window|success" in transcript
    ids, data, _, agemin, agemax = _inputs()
    missing_min = agemin.loc[agemin["IDPERS"] != 102].copy()
    na_max = agemax.copy()
    na_max.loc[na_max["IDPERS"] == 103, "AGE"] = np.nan

    result = traminex.createdatadiscrete(ids, data, ["A", "B"], missing_min, na_max)

    for frame in result["data"].values():
        assert 102.0 not in set(frame["IDPERS"])
        assert 103.0 not in set(frame["IDPERS"])
        assert set(frame["IDPERS"]) <= {101.0, 104.0}
    capsys.readouterr()


def test_ids_absent_from_data_are_skipped(capsys):
    ids, data, _, agemin, agemax = _inputs()
    result = traminex.createdatadiscrete(
        ids + [999],
        data,
        ["A", "B"],
        pd.concat([agemin, pd.DataFrame({"IDPERS": [999], "AGE": [18]})]),
        pd.concat([agemax, pd.DataFrame({"IDPERS": [999], "AGE": [20]})]),
    )
    assert all(999.0 not in set(frame["IDPERS"]) for frame in result["data"].values())
    capsys.readouterr()


@pytest.mark.parametrize(
    "target,column,match",
    [
        ("data", "BST", "BST"),
        ("agemin", "AGE", "agemin"),
        ("agemax", "IDPERS", "agemax"),
    ],
)
def test_malformed_schemas_are_refused(target, column, match):
    ids, data, _, agemin, agemax = _inputs()
    call_data, call_min, call_max = data, agemin, agemax
    if target == "data":
        call_data = data.drop(columns=column)
    elif target == "agemin":
        call_min = agemin.drop(columns=column)
    else:
        call_max = agemax.drop(columns=column)
    with pytest.raises(ValueError, match=match):
        traminex.createdatadiscrete(ids, call_data, ["A", "B"], call_min, call_max)


def test_duplicate_ids_status_and_window_validation():
    ids, data, _, agemin, agemax = _inputs()
    with pytest.raises(ValueError, match="ids must be unique"):
        traminex.createdatadiscrete(ids + [101], data, ["A", "B"], agemin, agemax)
    with pytest.raises(ValueError, match="data IDPERS must be unique"):
        traminex.createdatadiscrete(
            ids, pd.concat([data, data.iloc[[0]]]), ["A", "B"], agemin, agemax
        )
    invalid_status = data.copy()
    invalid_status.loc[0, "AST"] = 2
    with pytest.raises(ValueError, match="AST"):
        traminex.createdatadiscrete(ids, invalid_status, ["A", "B"], agemin, agemax)
    with pytest.raises(ValueError, match="agemin IDPERS must be unique"):
        traminex.createdatadiscrete(
            ids,
            data,
            ["A", "B"],
            pd.concat([agemin, agemin.iloc[[0]]]),
            agemax,
        )


def test_progress_messages_match_r_boundary(capsys):
    transcript = (R_GOLDENS_DIR / "discrete_rules_r3_errors.txt").read_text(encoding="utf-8")
    assert "progress|message|100% completed" in transcript
    ids = list(range(1, 1001))
    data = pd.DataFrame({"IDPERS": ids, "A": 2, "AST": 0, "B": 2, "BST": 0})
    agemin = pd.DataFrame({"IDPERS": ids, "AGE": 1})
    agemax = pd.DataFrame({"IDPERS": ids, "AGE": 1})

    traminex.createdatadiscrete(ids, data, ["A", "B"], agemin, agemax)

    assert capsys.readouterr().err.splitlines() == [
        "Creating data set for variable A",
        "100% completed",
        "Creating data set for variable B",
        "100% completed",
    ]


def test_createdatadiscrete_is_the_public_r_derived_name():
    assert traminex.createdatadiscrete.__module__ == "traminex.discrete_rules"


_RULE_LABELS = [
    "(A)",
    "(B)",
    "(C)",
    "(D)",
    "(A)-(B)",
    "(B)-(C)",
    "(A)-(B)-(C)",
    "(A,B)-(C)",
    "(A)-(D)",
]


def _rule_inputs(capsys):
    ids = list(range(1, 49))
    values = np.asarray(ids)
    data = pd.DataFrame(
        {
            "IDPERS": ids,
            "A": 19 + values % 4,
            "AST": (values % 4 != 0).astype(int),
            "B": 21 + values % 4,
            "BST": (values % 5 != 0).astype(int),
            "C": 23 + values % 3,
            "CST": (values % 3 != 0).astype(int),
            "D": 25,
            "DST": (values == 1).astype(int),
            "sex": 1 + values % 2,
            "cohort": 1 + values % 3,
        }
    )
    agemin = pd.DataFrame({"IDPERS": ids, "AGE": 18})
    agemax = pd.DataFrame({"IDPERS": ids, "AGE": 26})
    discrete = traminex.createdatadiscrete(
        ids, data, ["A", "B", "C", "D"], agemin, agemax, supvar="sex"
    )
    capsys.readouterr()

    event_rows = []
    for row in data.itertuples(index=False):
        for event in ("A", "B", "C", "D"):
            if getattr(row, f"{event}ST") == 1:
                event_rows.append(
                    {"id": row.IDPERS, "timestamp": getattr(row, event), "event": event}
                )
    eseq = traminex.seqecreate(pd.DataFrame(event_rows))
    fsub = traminex.seqefsub(eseq, str_subseq=_RULE_LABELS)
    return fsub, discrete, data


def _subset_fsub(fsub: SubSeqElist, positions: list[int]) -> SubSeqElist:
    return SubSeqElist(
        fsub.eseq,
        fsub.constraint,
        [fsub._subseq_objs[position] for position in positions],
        fsub.data.iloc[positions].reset_index(drop=True),
        type=fsub.type,
    )


def _rule_golden() -> pd.DataFrame:
    return pd.read_csv(
        R_GOLDENS_DIR / "discrete_seq_rules_r3_cases.csv",
        dtype=str,
        keep_default_na=False,
    )


def _expected_rule_case(case_name: str):
    records = _rule_golden()
    selected = records.loc[records["case"] == case_name]
    hazards = selected.loc[selected["record_type"] == "frame", "value"].tolist()
    graph = {}
    for key in ("events", "hr", "hrindep", "pvalues"):
        rows = selected.loc[selected["record_type"] == f"graph_{key}", "value"]
        values = []
        for value in rows:
            parts = value.split("\x1e") if value else []
            values.append(parts if key == "events" else [float(part) for part in parts])
        graph[key] = values
    columns_record = selected.loc[
        (selected["record_type"] == "metadata") & (selected["key"] == "columns"),
        "value",
    ].iloc[0]
    columns = columns_record.split("\x1e") if columns_record else []
    return hazards, graph, columns


def _assert_rule_case(result: pd.DataFrame, case_name: str) -> None:
    hazards, expected_graph, columns = _expected_rule_case(case_name)
    assert list(result.columns) == columns
    if columns:
        assert result["Hazards"].tolist() == hazards
    else:
        assert result.empty
    graph = result.attrs["graph"]
    assert graph["events"] == expected_graph["events"]
    assert graph["hrindep"] == expected_graph["hrindep"]
    for key in ("hr", "pvalues"):
        assert len(graph[key]) == len(expected_graph[key])
        for actual, expected in zip(graph[key], expected_graph[key], strict=True):
            np.testing.assert_allclose(actual, expected, rtol=0, atol=R_GOLDEN_ATOL)


@pytest.mark.parametrize(
    "case_name,kwargs",
    [
        ("cloglog_all", {"link": "cloglog"}),
        ("logit_all", {"link": "logit"}),
        ("topt", {"topt": True}),
        ("dep_c", {"dep": "C"}),
        ("dep_unknown", {"dep": "Z"}),
    ],
)
def test_seqerulesdisc_main_routes_match_r(case_name, kwargs, capsys):
    fsub, discrete, tsef = _rule_inputs(capsys)

    result = traminex.seqerulesdisc(fsub, discrete, tsef, pvalue=1, adjust=False, **kwargs)

    _assert_rule_case(result, case_name)
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""


def test_seqerulesdisc_adjustment_and_r_display_boundary(capsys):
    transcript = (R_GOLDENS_DIR / "discrete_seq_rules_r3_errors.txt").read_text(encoding="utf-8")
    assert 'cloglog_all|output|[1] "END"' in transcript
    fsub, discrete, tsef = _rule_inputs(capsys)

    result = traminex.seqerulesdisc(fsub, discrete, tsef, pvalue=1)

    _assert_rule_case(result, "adjusted")
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == "Adjusted p-value threshold : 0.111111111111111\n"


@pytest.mark.parametrize(
    "case_name,positions",
    [("duplicates", [4, 4, 5]), ("low_event", [8]), ("dep_unknown", [])],
)
def test_seqerulesdisc_subsequence_boundaries_match_r(case_name, positions, capsys):
    fsub, discrete, tsef = _rule_inputs(capsys)

    result = traminex.seqerulesdisc(
        _subset_fsub(fsub, positions), discrete, tsef, pvalue=1, adjust=False
    )

    _assert_rule_case(result, case_name)


@pytest.mark.parametrize("link", ["cloglog", "logit"])
def test_binomial_irls_matches_corrected_r_supvars_reference(link, capsys):
    from traminex.discrete_rules import _binomial_irls

    _, discrete, tsef = _rule_inputs(capsys)
    frame = discrete["data"]["B"].copy()
    frame["cohort"] = frame["IDPERS"].map(tsef.set_index("IDPERS")["cohort"])
    terms = ["(Intercept)", "A", "AGE", "sex", "cohort"]
    design = np.column_stack([np.ones(len(frame)), frame[["A", "AGE", "sex", "cohort"]].to_numpy()])

    fit = _binomial_irls(frame["B"].to_numpy(), design, link=link)

    records = _rule_golden()
    expected = records.loc[records["case"] == f"corrected_supvars_{link}"].pivot(
        index="row", columns="key", values="value"
    )
    expected = expected.loc[terms]
    np.testing.assert_allclose(
        fit.coefficients, expected["estimate"].astype(float), rtol=0, atol=R_GOLDEN_ATOL
    )
    np.testing.assert_allclose(
        fit.pvalues, expected["pvalue"].astype(float), rtol=0, atol=R_GOLDEN_ATOL
    )


def test_seqerulesdisc_corrects_broken_r_supvars_branch(capsys):
    transcript = (R_GOLDENS_DIR / "discrete_seq_rules_r3_errors.txt").read_text(encoding="utf-8")
    assert "supvars|simpleError|object 'model' not found" in transcript
    fsub, discrete, tsef = _rule_inputs(capsys)

    result = traminex.seqerulesdisc(
        fsub,
        discrete,
        tsef,
        pvalue=1,
        adjust=False,
        supvars=["sex", "cohort"],
    )

    assert not result.empty
    graph = result.attrs["graph"]
    assert len(graph["hrindep"]) == len(graph["hr"])
    assert all(
        len(rule_controls) == len(rule_edges)
        for rule_controls, rule_edges in zip(graph["hrindep"], graph["hr"], strict=True)
    )
    assert all(len(edge_controls) == 2 for rule in graph["hrindep"] for edge_controls in rule)


def test_seqerulesdisc_safely_returns_empty_for_r_broken_empty_input(capsys):
    transcript = (R_GOLDENS_DIR / "discrete_seq_rules_r3_errors.txt").read_text(encoding="utf-8")
    assert "empty_subseqs|subscriptOutOfBoundsError|subscript out of bounds" in transcript
    fsub, discrete, tsef = _rule_inputs(capsys)

    result = traminex.seqerulesdisc(_subset_fsub(fsub, []), discrete, tsef, pvalue=1, adjust=False)

    assert result.empty
    assert result.attrs["graph"] == {
        "events": [],
        "hr": [],
        "hrindep": [],
        "pvalues": [],
    }


def test_seqerulesdisc_validates_inputs(capsys):
    fsub, discrete, tsef = _rule_inputs(capsys)
    with pytest.raises(ValueError, match="link"):
        traminex.seqerulesdisc(fsub, discrete, tsef, link="bad")
    with pytest.raises(ValueError, match="pvalue"):
        traminex.seqerulesdisc(fsub, discrete, tsef, pvalue=1.1)
    with pytest.raises(ValueError, match="supvars"):
        traminex.seqerulesdisc(fsub, discrete, tsef, supvars=["unknown"])
    with pytest.raises(TypeError, match="SubSeqElist"):
        traminex.seqerulesdisc(object(), discrete, tsef)
    with pytest.raises(ValueError, match="datadiscr"):
        traminex.seqerulesdisc(fsub, {}, tsef)


def test_seqerulesdisc_is_public_with_the_full_signature():
    import inspect

    assert traminex.seqerulesdisc.__module__ == "traminex.discrete_rules"
    assert str(inspect.signature(traminex.seqerulesdisc)) == (
        "(fsubseq: 'SubSeqElist', datadiscr: 'Mapping[str, object]', "
        "tsef: 'pd.DataFrame', *, pvalue: 'float' = 0.1, "
        "supvars: 'Sequence[str] | None' = None, adjust: 'bool' = True, "
        "topt: 'bool' = False, link: \"Literal['cloglog', 'logit']\" = 'cloglog', "
        "dep: 'str | None' = None) -> 'pd.DataFrame'"
    )


def _fractional_inputs():
    ids = [201, 202, 203]
    data = pd.DataFrame(
        {
            "IDPERS": ids,
            "A": [20.5, 22, 19.75],
            "AST": [1, 0, 1],
            "B": [23.25, 21, 30],
            "BST": [1, 1, 0],
        }
    )
    agemin = pd.DataFrame({"IDPERS": ids, "AGE": [18, 19.5, 18]})
    agemax = pd.DataFrame({"IDPERS": ids, "AGE": [25, 24, 24.6]})
    return ids, data, ["A", "B"], agemin, agemax


def test_createdatadiscrete_fractional_ages_match_r(capsys):
    # TX-AUD-20261019-015: R accepts fractional ages; rows as.integer((end - start) + 1),
    # ages start, start + 1, ..., an independent event from row as.integer((age - start) + 1).
    ids, data, variables, agemin, agemax = _fractional_inputs()
    result = traminex.createdatadiscrete(ids, data, variables, agemin, agemax)
    for event in variables:
        pd.testing.assert_frame_equal(
            result["data"][event],
            _golden_frame("fractional_ages", event),
            check_exact=False,
            rtol=0,
            atol=R_GOLDEN_ATOL,
        )


@pytest.mark.parametrize("variables", [["A", "B"], ["B", "A"]])
def test_an_event_before_agemin_is_refused_as_r_stops(variables, capsys):
    conditions = (R_GOLDENS_DIR / "discrete_rules_r3_errors.txt").read_text(encoding="utf-8")
    assert "event_before_agemin|simpleError|invalid 'nrow' value (< 0)" in conditions
    assert (
        "independent_event_before_agemin|simpleError|only 0's may be mixed with negative subscripts"
        in conditions
    )
    data = pd.DataFrame({"IDPERS": [301], "A": [15], "AST": [1], "B": [21], "BST": [1]})
    with pytest.raises(ValueError, match="precedes"):
        traminex.createdatadiscrete(
            [301],
            data,
            variables,
            pd.DataFrame({"IDPERS": [301], "AGE": [18]}),
            pd.DataFrame({"IDPERS": [301], "AGE": [25]}),
        )
