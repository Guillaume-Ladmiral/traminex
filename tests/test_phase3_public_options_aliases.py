from __future__ import annotations

import inspect
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex.distances as distances
from traminex import (
    dissrep,
    disstree,
    seqdef,
    seqici,
    seqindic_dyn,
    seqrep,
    seqrep_grp,
    seqST,
    seqtree,
)
from traminex.events import (
    is_eseq,
    is_seqelist,
    seqeconstraint,
    seqecontain,
    seqecreate,
    seqefsub,
    seqeid,
    seqelength,
    seqetm,
    seqeweight,
)


def _distance_fixture() -> np.ndarray:
    return np.array(
        [
            [0.0, 0.2, 0.8, 1.0],
            [0.2, 0.0, 0.7, 0.9],
            [0.8, 0.7, 0.0, 0.3],
            [1.0, 0.9, 0.3, 0.0],
        ]
    )


def _sequence_fixture():
    return seqdef(
        np.array(
            [
                ["A", "A", "B"],
                ["A", "B", "B"],
                ["B", "B", "A"],
                ["B", "A", "A"],
            ]
        ),
        weights=[1.0, 2.0, 1.0, 1.0],
    )


def _event_fixture():
    return seqecreate(
        id=[1, 1, 2, 2],
        timestamp=[0, 2, 0, 1],
        event=["A", "B", "B", "A"],
    )


def _assert_rep_equal(actual: dict, expected: dict) -> None:
    np.testing.assert_array_equal(actual["representatives"], expected["representatives"])
    np.testing.assert_allclose(actual["distances"], expected["distances"])
    np.testing.assert_array_equal(actual["rep_group"], expected["rep_group"])
    pd.testing.assert_frame_equal(actual["statistics"], expected["statistics"])
    assert actual["quality"] == pytest.approx(expected["quality"])


@pytest.mark.parametrize(
    ("alias", "canonical", "value"),
    [("trep", "coverage", 0.75), ("tsim", "pradius", 0.55)],
)
def test_dissrep_deprecated_aliases_match_canonical_and_do_not_mutate(
    alias: str, canonical: str, value: float
) -> None:
    diss = _distance_fixture()
    original = diss.copy()
    expected = dissrep(diss, **{canonical: value})

    with pytest.warns(DeprecationWarning, match=rf"{alias}.*{canonical}") as caught:
        actual = dissrep(diss, **{alias: value})

    _assert_rep_equal(actual, expected)
    np.testing.assert_array_equal(diss, original)
    # Resolved: through a symbolic link the code's file name and __file__ can name
    # the same file by two paths (TX-AUD-20261010-090 (a)).
    assert Path(caught[0].filename).resolve() == Path(__file__).resolve()


def test_dissrep_alias_default_precedence_conflict_and_unknown() -> None:
    diss = _distance_fixture()
    with pytest.warns(DeprecationWarning):
        actual = dissrep(diss, coverage=0.25, trep=0.75)
    expected = dissrep(diss, coverage=0.75)
    _assert_rep_equal(actual, expected)

    with pytest.raises(ValueError, match="coverage.*trep"):
        dissrep(diss, coverage=0.5, trep=0.75)
    with pytest.raises(TypeError, match="unexpected keyword"):
        dissrep(diss, not_an_option=True)


@pytest.mark.parametrize(
    ("alias", "canonical", "value"),
    [("minSize", "min_size", 0.2), ("maxdepth", "max_depth", 1)],
)
def test_disstree_deprecated_aliases_match_canonical(
    alias: str, canonical: str, value: float | int
) -> None:
    diss = _distance_fixture()
    predictors = pd.DataFrame({"group": ["a", "a", "b", "b"]})
    expected = disstree(diss, predictors, R=1, **{canonical: value})

    with pytest.warns(DeprecationWarning, match=rf"{alias}.*{canonical}"):
        actual = disstree(diss, predictors, R=1, **{alias: value})

    np.testing.assert_array_equal(actual.fitted, expected.fitted)
    np.testing.assert_array_equal(actual.root.info["ind"], expected.root.info["ind"])
    for key in ("depth", "vardis", "n", "medoid", "splitschedule"):
        assert actual.root.info[key] == expected.root.info[key]


def test_disstree_alias_default_precedence_and_conflict() -> None:
    diss = _distance_fixture()
    predictors = pd.DataFrame({"group": ["a", "a", "b", "b"]})
    with pytest.warns(DeprecationWarning):
        actual = disstree(diss, predictors, min_size=0.05, minSize=0.2, R=1)
    expected = disstree(diss, predictors, min_size=0.2, R=1)
    np.testing.assert_array_equal(actual.fitted, expected.fitted)

    with pytest.raises(ValueError, match="min_size.*minSize"):
        disstree(diss, predictors, min_size=0.1, minSize=0.2, R=1)


@pytest.mark.parametrize(
    ("alias", "canonical", "value"),
    [
        ("trep", "coverage", 0.75),
        ("tsim", "pradius", 0.55),
        ("dist_matrix", "diss", _distance_fixture()),
    ],
)
def test_seqrep_deprecated_aliases_match_canonical(
    alias: str, canonical: str, value: float | np.ndarray
) -> None:
    seqdata = _sequence_fixture()
    # A supplied matrix leaves seqdist's options nothing to act on (they are refused).
    distance = {} if canonical == "diss" else {"method": "LCS"}
    expected = seqrep(seqdata, criterion="density", **distance, **{canonical: value})

    with pytest.warns(DeprecationWarning, match=rf"{alias}.*{canonical}"):
        actual = seqrep(seqdata, criterion="density", **distance, **{alias: value})

    np.testing.assert_array_equal(actual["rep_indices"], expected["rep_indices"])
    np.testing.assert_allclose(actual["distance_matrix"], expected["distance_matrix"])
    pd.testing.assert_frame_equal(actual["statistics"], expected["statistics"])


def test_seqrep_alias_default_precedence_and_conflict() -> None:
    seqdata = _sequence_fixture()
    diss = _distance_fixture()
    with pytest.warns(DeprecationWarning):
        actual = seqrep(seqdata, diss=None, dist_matrix=diss)
    expected = seqrep(seqdata, diss=diss)
    np.testing.assert_array_equal(actual["rep_indices"], expected["rep_indices"])

    with pytest.raises(ValueError, match="diss.*dist_matrix"):
        seqrep(seqdata, diss=diss, dist_matrix=diss.copy())
    with pytest.raises(ValueError, match="coverage.*trep"):
        seqrep(seqdata, diss=diss, coverage=0.5, trep=0.75)


@pytest.mark.parametrize(
    ("alias", "canonical", "value"),
    [
        ("minSize", "min_size", 0.2),
        ("maxdepth", "max_depth", 1),
        ("seqdist_arg", "seqdist_args", {"method": "HAM"}),
    ],
)
def test_seqtree_deprecated_aliases_match_canonical(
    alias: str, canonical: str, value: object
) -> None:
    seqdata = _sequence_fixture()
    predictors = pd.DataFrame({"group": ["a", "a", "b", "b"]})
    kwargs = {"diss": _distance_fixture()} if canonical != "seqdist_args" else {}
    expected = seqtree(seqdata, predictors, R=1, **kwargs, **{canonical: value})

    with pytest.warns(DeprecationWarning, match=rf"{alias}.*{canonical}"):
        actual = seqtree(seqdata, predictors, R=1, **kwargs, **{alias: value})

    np.testing.assert_array_equal(actual.fitted, expected.fitted)


def test_seqtree_alias_default_precedence_and_conflict() -> None:
    seqdata = _sequence_fixture()
    predictors = pd.DataFrame({"group": ["a", "a", "b", "b"]})
    diss = _distance_fixture()
    with pytest.warns(DeprecationWarning):
        actual = seqtree(seqdata, predictors, diss=diss, min_size=0.05, minSize=0.2, R=1)
    expected = seqtree(seqdata, predictors, diss=diss, min_size=0.2, R=1)
    np.testing.assert_array_equal(actual.fitted, expected.fitted)

    with pytest.raises(ValueError, match="seqdist_args.*seqdist_arg"):
        seqtree(
            seqdata,
            predictors,
            seqdist_args={"method": "OM"},
            seqdist_arg={"method": "HAM"},
            R=1,
        )


def test_seqrep_grp_mdis_alias_matches_canonical_and_does_not_mutate() -> None:
    seqdata = _sequence_fixture()
    group = np.array(["a", "a", "b", "b"], dtype=object)
    diss = _distance_fixture()
    original_group = group.copy()
    original_diss = diss.copy()
    expected = seqrep_grp(seqdata, group, diss=diss)

    with pytest.warns(DeprecationWarning, match=r"mdis.*diss"):
        actual = seqrep_grp(seqdata, group, mdis=diss)

    assert actual.keys() == expected.keys()
    for key in actual:
        pd.testing.assert_frame_equal(actual[key], expected[key])
    np.testing.assert_array_equal(group, original_group)
    np.testing.assert_array_equal(diss, original_diss)


def test_seqrep_grp_mdis_default_precedence_and_conflict() -> None:
    seqdata = _sequence_fixture()
    group = ["a", "a", "b", "b"]
    diss = _distance_fixture()
    with pytest.warns(DeprecationWarning):
        actual = seqrep_grp(seqdata, group, diss=None, mdis=diss)
    expected = seqrep_grp(seqdata, group, diss=diss)
    for key in actual:
        pd.testing.assert_frame_equal(actual[key], expected[key])

    with pytest.raises(ValueError, match="diss.*mdis"):
        seqrep_grp(seqdata, group, diss=diss, mdis=diss.copy())


def test_seqrep_grp_with_missing_is_forwarded_and_typed(monkeypatch) -> None:
    seqdata = _sequence_fixture()
    captured: dict[str, object] = {}

    def fake_seqdist(arg, **kwargs):
        captured["arg"] = arg
        captured["kwargs"] = kwargs
        return _distance_fixture()

    monkeypatch.setattr(distances, "seqdist", fake_seqdist)
    seqrep_grp(
        seqdata,
        ["a", "a", "b", "b"],
        with_missing=True,
        method="HAM",
    )
    assert captured["arg"] is seqdata
    assert captured["kwargs"] == {"with_missing": True, "method": "HAM"}

    with pytest.raises(TypeError, match="with_missing.*boolean"):
        seqrep_grp(seqdata, ["a", "a", "b", "b"], with_missing=1)


@pytest.mark.parametrize(
    ("function", "names"),
    [
        (dissrep, {"trep", "tsim"}),
        (disstree, {"minSize", "maxdepth"}),
        (seqrep, {"trep", "tsim", "dist_matrix"}),
        (seqtree, {"minSize", "maxdepth", "seqdist_arg"}),
        (seqrep_grp, {"with_missing", "mdis"}),
    ],
)
def test_phase3_representative_tree_runtime_signatures(function, names: set[str]) -> None:
    assert names <= set(inspect.signature(function).parameters)


@pytest.mark.parametrize(
    ("function", "canonical", "value"),
    [
        (is_eseq, "x", lambda eseq: eseq[0]),
        (is_seqelist, "x", lambda eseq: eseq),
        (seqeid, "eseq", lambda eseq: eseq),
        (seqelength, "eseq", lambda eseq: eseq),
        (seqeweight, "eseq", lambda eseq: eseq),
    ],
)
def test_event_s_alias_matches_required_canonical(function, canonical, value) -> None:
    eseq = _event_fixture()
    argument = value(eseq)
    expected = function(**{canonical: argument})

    with pytest.warns(DeprecationWarning, match=rf"s.*{canonical}") as caught:
        actual = function(s=argument)

    if isinstance(expected, np.ndarray):
        np.testing.assert_array_equal(actual, expected)
    else:
        assert actual == expected
    # Resolved: through a symbolic link the code's file name and __file__ can name
    # the same file by two paths (TX-AUD-20261010-090 (a)).
    assert Path(caught[0].filename).resolve() == Path(__file__).resolve()

    with pytest.raises(ValueError, match=rf"{canonical}.*s"):
        function(**{canonical: argument, "s": argument})


def test_seqetm_seq_alias_matches_canonical_and_conflicts() -> None:
    seqdata = _sequence_fixture()
    expected = seqetm(seqdata, method="state")
    original = seqdata.data.copy()

    with pytest.warns(DeprecationWarning, match=r"seq.*x"):
        actual = seqetm(seq=seqdata, method="state")

    pd.testing.assert_frame_equal(actual, expected)
    np.testing.assert_array_equal(seqdata.data, original)
    with pytest.raises(ValueError, match=r"x.*seq"):
        seqetm(seqdata, seq=seqdata)


@pytest.mark.parametrize(
    ("alias", "canonical", "value"),
    [
        ("maxGap", "max_gap", 2),
        ("windowSize", "window_size", 4),
        ("ageMin", "age_min", 1),
        ("ageMax", "age_max", 6),
        ("ageMaxEnd", "age_max_end", 7),
        ("countMethod", "count_method", "CDIST"),
    ],
)
def test_seqeconstraint_camel_aliases_match_canonical(alias, canonical, value) -> None:
    expected = seqeconstraint(**{canonical: value})
    with pytest.warns(DeprecationWarning, match=rf"{alias}.*{canonical}"):
        actual = seqeconstraint(**{alias: value})
    assert actual == expected


def test_seqeconstraint_alias_default_precedence_and_conflict() -> None:
    with pytest.warns(DeprecationWarning):
        actual = seqeconstraint(max_gap=-1, maxGap=2)
    assert actual == seqeconstraint(max_gap=2)
    with pytest.raises(ValueError, match=r"max_gap.*maxGap"):
        seqeconstraint(max_gap=1, maxGap=2)


@pytest.mark.parametrize(
    ("alias", "canonical", "value"),
    [
        ("seq", "eseq", lambda eseq: eseq),
        ("eventList", "event_list", lambda eseq: ["A"]),
        ("exclude", "unknown_exclude", lambda eseq: True),
    ],
)
def test_seqecontain_aliases_match_canonical(alias, canonical, value) -> None:
    eseq = _event_fixture()
    base = {"eseq": eseq, "event_list": ["A"]}
    if canonical in base:
        base.pop(canonical)
    resolved = value(eseq)
    expected = seqecontain(**base, **{canonical: resolved})

    with pytest.warns(DeprecationWarning, match=rf"{alias}.*{canonical}"):
        actual = seqecontain(**base, **{alias: resolved})

    np.testing.assert_array_equal(actual, expected)


def test_seqecontain_alias_default_precedence_conflicts_and_unknown() -> None:
    eseq = _event_fixture()
    with pytest.warns(DeprecationWarning):
        actual = seqecontain(eseq, ["A"], unknown_exclude=False, exclude=True)
    np.testing.assert_array_equal(actual, seqecontain(eseq, ["A"], unknown_exclude=True))
    with pytest.raises(ValueError, match=r"unknown_exclude.*exclude"):
        seqecontain(eseq, ["A"], unknown_exclude=True, exclude=False)
    with pytest.raises(TypeError, match="unexpected keyword"):
        seqecontain(eseq, ["A"], discarded=True)


def test_seqecreate_endEvent_alias_matches_canonical_without_mutation() -> None:
    frame = pd.DataFrame(
        {
            "id": [1, 1, 2, 2],
            "timestamp": [0, 3, 0, 4],
            "event": ["A", "END", "B", "END"],
        }
    )
    original = frame.copy(deep=True)
    expected = seqecreate(frame, end_event="END")

    with pytest.warns(DeprecationWarning, match=r"endEvent.*end_event"):
        actual = seqecreate(frame, endEvent="END")

    assert actual.as_character() == expected.as_character()
    np.testing.assert_array_equal(actual.lengths, expected.lengths)
    pd.testing.assert_frame_equal(frame, original)


def test_seqecreate_endEvent_default_precedence_and_conflict() -> None:
    kwargs = {
        "id": [1, 1],
        "timestamp": [0, 2],
        "event": ["A", "END"],
    }
    with pytest.warns(DeprecationWarning):
        actual = seqecreate(**kwargs, end_event=None, endEvent="END")
    expected = seqecreate(**kwargs, end_event="END")
    assert actual.as_character() == expected.as_character()
    with pytest.raises(ValueError, match=r"end_event.*endEvent"):
        seqecreate(**kwargs, end_event="STOP", endEvent="END")


@pytest.mark.parametrize(
    ("alias", "canonical", "value", "base"),
    [
        ("seq", "eseq", lambda eseq: eseq, {"str_subseq": ["(A)"]}),
        (
            "strsubseq",
            "str_subseq",
            lambda eseq: ["(A)", "(B)"],
            {},
        ),
        ("minSupport", "min_support", lambda eseq: 1.0, {"max_k": 1}),
        ("pMinSupport", "pmin_support", lambda eseq: 0.5, {"max_k": 1}),
        ("maxK", "max_k", lambda eseq: 1, {"min_support": 1.0}),
    ],
)
def test_seqefsub_aliases_match_canonical(alias, canonical, value, base) -> None:
    eseq = _event_fixture()
    call = {"eseq": eseq, **base}
    if canonical in call:
        call.pop(canonical)
    resolved = value(eseq)
    expected = seqefsub(**call, **{canonical: resolved})

    with pytest.warns(DeprecationWarning, match=rf"{alias}.*{canonical}"):
        actual = seqefsub(**call, **{alias: resolved})

    assert actual.subseq == expected.subseq
    pd.testing.assert_frame_equal(actual.data, expected.data)


def test_seqefsub_alias_default_precedence_and_conflict() -> None:
    eseq = _event_fixture()
    with pytest.warns(DeprecationWarning):
        actual = seqefsub(
            eseq,
            str_subseq=None,
            strsubseq=["(A)"],
        )
    expected = seqefsub(eseq, str_subseq=["(A)"])
    pd.testing.assert_frame_equal(actual.data, expected.data)
    with pytest.raises(ValueError, match=r"min_support.*minSupport"):
        seqefsub(eseq, min_support=2, minSupport=1)


@pytest.mark.parametrize(
    ("function", "names"),
    [
        (is_eseq, {"x", "s"}),
        (is_seqelist, {"x", "s"}),
        (seqeid, {"eseq", "s"}),
        (seqelength, {"eseq", "s"}),
        (seqeweight, {"eseq", "s"}),
        (seqetm, {"x", "seq"}),
        (
            seqeconstraint,
            {
                "maxGap",
                "windowSize",
                "ageMin",
                "ageMax",
                "ageMaxEnd",
                "countMethod",
            },
        ),
        (seqecontain, {"seq", "eventList", "exclude"}),
        (seqecreate, {"endEvent"}),
        (
            seqefsub,
            {"seq", "strsubseq", "minSupport", "pMinSupport", "maxK"},
        ),
    ],
)
def test_phase3_event_runtime_signatures(function, names: set[str]) -> None:
    assert names <= set(inspect.signature(function).parameters)


def test_seqici_silent_controls_exact_r_message_and_preserves_output(capsys) -> None:
    seqdata = _sequence_fixture()
    original = seqdata.data.copy()
    quiet = seqici(seqdata, silent=True)
    assert capsys.readouterr().err == ""

    verbose = seqici(seqdata, silent=False)
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == " [>] computing complexity index for 4 sequences ...\n"
    np.testing.assert_allclose(verbose, quiet)
    np.testing.assert_array_equal(seqdata.data, original)

    with pytest.raises(TypeError, match="silent.*boolean"):
        seqici(seqdata, silent=1)


def test_seqST_silent_controls_exact_r_messages_and_preserves_output(capsys) -> None:
    seqdata = _sequence_fixture()
    original = seqdata.data.copy()
    quiet = seqST(seqdata, silent=True)
    assert capsys.readouterr().err == ""

    verbose = seqST(seqdata, silent=False)
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == (
        " [>] extracting symbols and durations ...\n"
        " [>] computing turbulence type 1 for 4 sequence(s) ...\n"
    )
    np.testing.assert_allclose(verbose, quiet)
    np.testing.assert_array_equal(seqdata.data, original)

    with pytest.raises(TypeError, match="silent.*boolean"):
        seqST(seqdata, silent="yes")


def test_seqindic_dyn_silent_indic_suppresses_only_indicator_messages(monkeypatch, capsys) -> None:
    import traminex.stats_utils as stats_utils

    seqdata = _sequence_fixture()
    original = seqdata.data.copy()

    def noisy_seqindic(x, **kwargs):
        print("indicator message", file=sys.stderr)
        return pd.DataFrame({"C": np.arange(x.nseq, dtype=float)})

    monkeypatch.setattr(stats_utils, "seqindic", noisy_seqindic)
    quiet = seqindic_dyn(
        seqdata,
        window_size=2,
        wstep=1,
        silent_indic=True,
    )
    assert capsys.readouterr().err == ""

    verbose = seqindic_dyn(
        seqdata,
        window_size=2,
        wstep=1,
        silent_indic=False,
    )
    captured = capsys.readouterr()
    assert captured.err == "indicator message\nindicator message\n"
    np.testing.assert_array_equal(verbose, quiet)
    np.testing.assert_array_equal(seqdata.data, original)

    with pytest.raises(TypeError, match="silent_indic.*boolean"):
        seqindic_dyn(seqdata, window_size=2, silent_indic=1)


@pytest.mark.parametrize(
    ("function", "name"),
    [(seqici, "silent"), (seqST, "silent"), (seqindic_dyn, "silent_indic")],
)
def test_phase3_silent_option_runtime_signatures(function, name: str) -> None:
    assert name in inspect.signature(function).parameters
