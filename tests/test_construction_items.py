"""Batch-14 sequence-construction, format and multidomain items against R (``make_goldens_construction_items.R``).

The golden records what TraMineR and TraMineRextras return for the inputs these tests use, R's errors included where
traminex corrects R; each test states whether it expects R's value or a correction derived from a stated rule.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx

GOLDEN_PATH = Path(__file__).resolve().parent / "r_goldens" / "construction_items.csv"
with open(GOLDEN_PATH, encoding="utf-8", newline="") as handle:
    GOLDEN = {(row["case"], row["field"]): row["value"] for row in csv.DictReader(handle)}


def _numbers(case: str, field: str) -> np.ndarray:
    return np.array([np.nan if v == "NA" else float(v) for v in GOLDEN[(case, field)].split()])


def _texts(case: str, field: str) -> list[str]:
    return GOLDEN[(case, field)].split(" ;; ")


def _cells(case: str) -> np.ndarray:
    return np.array(
        [[None if v == "NA" else v for v in row.split()] for row in _texts(case, "cells")],
        dtype=object,
    )


# TX-AUD-20261019-033 (owner decision 2026-10-23): R's otherwise recodes every unlisted cell, the missing and void
# markers included (seqrecode.R:69-77), so the void padding becomes a state and the shorter sequence lengthens (the
# golden). traminex applies otherwise to the real states only: the markers stay markers unless recodes names them.
def test_seqrecode_otherwise_leaves_the_markers_alone() -> None:
    x = tx.seqdef(_cells("recode"))
    assert list(tx.seqlength(x)) == list(_numbers("recode", "input_lengths"))
    assert _texts("recode_otherwise", "sequences") == ["A-O-O-O", "A-O-O-O"]  # R
    y = tx.seqrecode(x, {"A": ["a"]}, otherwise="O")
    assert y.to_strings() == ["A-O-*-O", "A-O-%-%"]
    assert list(tx.seqlength(y)) == list(_numbers("recode", "input_lengths"))
    assert tx.alphabet(y) == ["A", "O"]


def test_seqrecode_recodes_a_marker_named_in_recodes() -> None:
    x = tx.seqdef(_cells("recode"))
    assert _texts("recode_missing_named", "sequences") == ["A-O-M-O", "A-O-O-O"]  # R
    y = tx.seqrecode(x, {"A": ["a"], "M": ["*"]}, otherwise="O")
    assert y.to_strings() == ["A-O-M-O", "A-O-%-%"]
    assert list(tx.seqlength(y)) == list(_numbers("recode", "input_lengths"))
    assert tx.alphabet(y) == ["A", "M", "O"]


# Batch-14 review I-3: otherwise changes only the real states recodes does not list, so with it a recode gives the
# result it gives without it, the unlisted real states replaced. A recode named like a used marker (supported, the
# seqrecode row) raised a duplicated-alphabet error with otherwise; and an otherwise spelled like an unused marker but
# naming a real state turned the void padding into that state.
@pytest.mark.parametrize("marker", ["*", "%"])
def test_seqrecode_otherwise_with_a_recode_named_like_a_used_marker(marker: str) -> None:
    x = tx.seqdef(_cells("recode"))
    plain = tx.seqrecode(x, {marker: ["a"]})
    y = tx.seqrecode(x, {marker: ["a"]}, otherwise="O")
    assert y.to_strings() == [s.replace("b", "O").replace("c", "O") for s in plain.to_strings()]
    assert list(tx.seqlength(y)) == list(tx.seqlength(plain))


def test_seqrecode_otherwise_naming_a_real_state_leaves_the_void_padding() -> None:
    from traminex.sequences import STSList

    # A real state "*" spelled like the nr marker, which no cell uses; void pads the second row.
    x = STSList(data=np.array([[1, 0, 2], [1, 3, 3]]), alphabet=["*", "a", "b", "%"], void_code=3)
    assert x.nr == "*" and x.nr_code is None
    y = tx.seqrecode(x, {"A": ["a"]}, otherwise="*")
    assert y.to_strings() == ["A-*-*", "A-%-%"]
    assert list(tx.seqlength(y)) == list(tx.seqlength(x)) == [3, 1]


# TX-AUD-20261019-035 and TX-REP-20261019-013: R's route for long data is reshape to wide, then seqdef, which pairs
# labels and cpal with the real states and gives the markers their own (missing.color; void white in traminex's
# seqdef). The long builders and seqformat's STS output carry the same labels and colours as seqdef on the wide table.
def _long(case: str) -> pd.DataFrame:
    states = [None if v == "NA" else v for v in _texts(case, "state")]
    return pd.DataFrame(
        {
            "id": _numbers(case, "id").astype(int),
            "time": _numbers(case, "time").astype(int),
            "state": states,
        }
    )


def _wide(long: pd.DataFrame) -> pd.DataFrame:
    rows = {}
    for person, time, state in long.itertuples(index=False):
        rows.setdefault(person, {})[time] = state
    times = sorted(long["time"].unique())
    return pd.DataFrame(
        [[rows[p].get(t) for t in times] for p in rows], columns=[f"t{t}" for t in times]
    )


def _metadata(x: tx.STSList) -> dict:
    from traminex.sequences import _real_state_codes

    real = _real_state_codes(x)
    return {
        "labels": [x.labels[code] for code in real],
        "cpal": [x.cpal[code].upper() for code in real],
        "nr": None if x.nr_code is None else (x.labels[x.nr_code], x.cpal[x.nr_code]),
        "void": None if x.void_code is None else (x.labels[x.void_code], x.cpal[x.void_code]),
    }


def _long_builders(long: pd.DataFrame, **kwargs) -> dict:
    builders = {
        "seqdef_from_long": tx.seqdef_from_long(long, **kwargs),
        "TSE_to_STS": tx.TSE_to_STS(long, **kwargs),
    }
    if not kwargs:  # R's seqformat takes no labels; its result goes to seqdef
        builders["seqformat"] = tx.seqformat(_wide(long), from_="STS", to="STS")
    return builders


def test_long_builders_carry_seqdef_s_labels_and_colours() -> None:
    long = _long("long_defaults")
    reference = _metadata(tx.seqdef(_wide(long)))
    assert reference["labels"] == _texts("long_defaults", "stlab")  # R
    assert reference["cpal"] == _texts("long_defaults", "cpal")  # R
    for name, x in _long_builders(long).items():
        assert _metadata(x) == reference, name
    assert reference["void"] == ("%", "white")


def test_tse_to_sts_event_mode_carries_seqdef_s_colours() -> None:
    events = pd.DataFrame({"id": [1, 1, 2], "time": [1, 3, 2], "event": ["a", "b", "a"]})
    x = tx.TSE_to_STS(events, tmax=4)
    cells = [[None if s in (x.nr, x.void) else s for s in row.split("-")] for row in x.to_strings()]
    reference = _metadata(tx.seqdef(pd.DataFrame(cells)))
    assert _metadata(x)["labels"] == reference["labels"]
    assert _metadata(x)["cpal"] == reference["cpal"]


def test_long_builders_take_labels_of_the_real_states() -> None:
    long = _long("long_real_labels")
    labels = _texts("long_real_labels", "labels_given")
    assert _texts("long_real_labels", "stlab") == labels  # R accepts them
    for name, x in _long_builders(long, labels=labels).items():
        meta = _metadata(x)
        assert meta["labels"] == labels, name
        assert meta["nr"] == ("*", _texts("long_real_labels", "missing_color")[0]), name


# TX-AUD-20261019-037: HSPELL_to_STS's per-row tmin_col/tmax_col (traminex extensions) mark each person's observation
# window. As R's endObs, the periods outside it are unobserved and fillblanks does not fill them: leading ones become
# missing and trailing ones void by seqdef's defaults. Inside the window fillblanks fills the gaps between spells.
_HSPELL = pd.DataFrame(
    {
        "id": ["p1", "p2"],
        "b1": [1, 2],
        "e1": [3, 3],
        "s1": ["A", "B"],
        "b2": [6, 5],
        "e2": [8, 8],
        "s2": ["C", "C"],
        "lo": [1, 4],
        "hi": [7, 8],
        "obs": [8, 7],
    }
)
_HSPELL_ARGS = dict(
    id_col="id", begin=["b1", "b2"], end=["e1", "e2"], status=["s1", "s2"], tmin=1, tmax=8
)


def _hspell_rows(x: tx.STSList) -> list[list[str]]:
    return [row.split("-") for row in x.to_strings()]


def test_hspell_bounds_are_not_filled() -> None:
    x = tx.HSPELL_to_STS(_HSPELL, **_HSPELL_ARGS, tmin_col="lo", tmax_col="hi", fillblanks="G")
    assert _hspell_rows(x) == [
        ["A", "A", "A", "G", "G", "C", "C", "%"],
        ["*", "*", "*", "G", "C", "C", "C", "C"],
    ]


def test_hspell_bounds_and_end_of_observation_together() -> None:
    x = tx.HSPELL_to_STS(
        _HSPELL, **_HSPELL_ARGS, tmin_col="lo", tmax_col="hi", endObs="obs", fillblanks="G"
    )
    assert _hspell_rows(x) == [
        ["A", "A", "A", "G", "G", "C", "C", "%"],
        ["*", "*", "*", "G", "C", "C", "%", "%"],
    ]


# TX-AUD-20261019-039: R's route for long data (reshape to wide, then seqdef) and TSE_to_STS's event mode keep the
# persons in order of first appearance, for text and numeric ids (the golden). traminex's state mode sorted them.
@pytest.mark.parametrize("case", ["order_text", "order_numeric"])
def test_long_builders_keep_first_appearance_order(case: str) -> None:
    ids = _numbers(case, "id").astype(int) if case == "order_numeric" else _texts(case, "id")
    long = pd.DataFrame(
        {"id": ids, "time": _numbers(case, "time").astype(int), "state": _texts(case, "state")}
    )
    route = _texts(case, "route_ids")
    assert [str(i) for i in tx.seqdef_from_long(long).ids] == route
    assert [str(i) for i in tx.TSE_to_STS(long).ids] == route
    events = pd.DataFrame(
        {"id": list(dict.fromkeys(ids)), "time": [0, 0, 0], "event": ["A", "B", "A"]}
    )
    assert [str(i) for i in tx.TSE_to_STS(events, tmax=2).ids] == _texts(case, "event_ids")


# Batch-14 review I-2: R's route reshapes a person with no id like the others, then seqdef with the ids as row names
# stops (the golden); the long builders dropped that person without a message. They refuse it with R's message.
@pytest.mark.parametrize("builder", ["seqdef_from_long", "TSE_to_STS"])
def test_long_builders_refuse_a_missing_id_as_r_does(builder: str) -> None:
    long = pd.DataFrame(
        {
            "id": _numbers("long_missing_id", "id"),
            "time": _numbers("long_missing_id", "time").astype(int),
            "state": _texts("long_missing_id", "state"),
        }
    )
    message = GOLDEN[("long_missing_id", "result")].removeprefix("ERROR: ")
    with pytest.raises(ValueError, match=re.escape(message)):
        getattr(tx, builder)(long)


# Batch-14 review I-1: the persons' order is found in one pass; a set of the ids was rebuilt for every person, so
# 8,000 persons took about 5 seconds and 20,000 about 30. Linear, 20,000 persons take well under a second.
def test_long_builders_order_many_persons_in_linear_time() -> None:
    import time

    n = 20_000
    long = pd.DataFrame(
        {
            "id": np.repeat(np.arange(n, 0, -1), 2),
            "time": np.tile([1, 2], n),
            "state": np.tile(["A", "B"], n),
        }
    )
    start = time.perf_counter()
    x = tx.seqdef_from_long(long)
    elapsed = time.perf_counter() - start
    assert [int(i) for i in x.ids[:3]] == [n, n - 1, n - 2]
    assert elapsed < 10.0, elapsed


# TX-AUD-20261019-040:R's seqdef takes left, gaps or right equal to the nr marker "*" and lists "*" as a state of the
# alphabet (twice when two of them name it; the golden). traminex refuses a replacement state spelled like the nr
# marker, and the seqdef ledger row records the correction.
def _policy_frame() -> pd.DataFrame:
    rows = [[None if v == "NA" else v for v in row.split()] for row in _texts("policy", "cells")]
    return pd.DataFrame(rows, columns=["a", "b", "c"])


@pytest.mark.parametrize("case", ["policy_left", "policy_gaps_left", "policy_right"])
def test_seqdef_refuses_the_nr_marker_as_a_replacement_state(case: str) -> None:
    options = dict(item.split(" = ") for item in _texts(case, "options"))
    assert "*" in _texts(case, "alphabet")  # R lists the marker as a state
    with pytest.raises(
        ValueError, match="missing replacement states must differ from the 'nr' marker"
    ):
        tx.seqdef(_policy_frame(), **options)


def test_seqdef_ledger_row_records_the_replacement_refusal() -> None:
    assert _texts("policy_gaps_left", "alphabet") == ["A", "B", "*", "*"]  # R's duplicated entry
    ledger = Path(__file__).resolve().parents[1] / "docs" / "API_PARITY_DECISIONS.csv"
    with open(ledger, encoding="utf-8", newline="") as handle:
        (row,) = [r for r in csv.DictReader(handle) if r["r_name"] == "seqdef"]
    assert "TX-AUD-20261019-040" in row["notes"]
    assert "construction_items" in row["evidence"]


# TX-AUD-20261019-041 and TX-AUD-20261010-075: the multidomain weight helpers are traminex extensions. Their inputs are
# checked up front with named errors, their docstring examples run, and docs/INDEX.md names every Config class.
def _two_channels() -> tx.SeqMD:
    ch1 = tx.seqdef(np.array([[0, 1, 2], [0, 1, 2], [2, 1, 0], [2, 1, 0]]))
    ch2 = tx.seqdef(np.array([[0, 0, 1], [0, 0, 1], [1, 1, 0], [1, 1, 0]]))
    return tx.SeqMD([ch1, ch2], names=["work", "family"])


@pytest.mark.parametrize("step", [1.5, 2.0, 5.0])
def test_gridsearch_refuses_a_weight_step_above_one(step: float) -> None:
    # The grid runs from 0 to 1 in steps of weight_step, so a step above 1 leaves no usable vector.
    with pytest.raises(ValueError, match="weight_step must be at most 1"):
        tx.seqMD_gridsearch_weights(_two_channels(), [0, 0, 1, 1], method="HAM", weight_step=step)


@pytest.mark.parametrize("helper", ["seqMD_gridsearch_weights", "seqMD_optimize_weights"])
def test_weight_helpers_refuse_missing_group_labels(helper: str) -> None:
    with pytest.raises(ValueError, match="groups must not contain missing values"):
        getattr(tx, helper)(_two_channels(), [0, 0, np.nan, np.nan], method="HAM")


def test_seqmd_refuses_names_of_the_wrong_length() -> None:
    md = _two_channels()
    with pytest.raises(ValueError, match="names must give one name per channel"):
        tx.SeqMD(md.channels, names=["only-one"])


def test_multidomain_examples_run() -> None:
    import doctest

    import traminex.multidomain as module

    result = doctest.testmod(module, verbose=False)
    assert (result.failed, result.attempted) == (0, result.attempted) and result.attempted > 0


def test_index_names_every_config_class() -> None:
    import traminex.config_dataclasses as configs

    index = (Path(__file__).resolve().parents[1] / "docs" / "INDEX.md").read_text(encoding="utf-8")
    classes = sorted(name for name in dir(configs) if name.endswith("Config"))
    assert classes == ["DHDConfig", "DURConfig", "LCSLCPConfig", "OMConfig", "TRATEConfig"]
    assert [name for name in classes if name not in index] == []


# TX-AUD-20261019-042 and TX-AUD-20261010-094: both marker representations of the same data give the same results
# (tests/test_external_marker_equivalence.py). An object declaring the nr marker's text in its alphabet beside an
# external nr_code (the 10-094 record's object) is read through its canonical form by every per-state accessor and
# statistic; and an STSList refuses nr text equal to void text, as seqdef does.
def _record_twins():
    from traminex.sequences import STSList

    x = STSList(
        data=np.array([[0, 1, 3], [1, 1, 0], [0, 3, 1]]), alphabet=["A", "B", "*"], nr_code=3
    )
    return x, x._canonicalized()


def _same(a, b) -> bool:
    from traminex.sequences import STSList

    if a is None or b is None:
        return a is b
    if isinstance(a, STSList):
        return isinstance(b, STSList) and all(
            _same(getattr(a, name), getattr(b, name))
            for name in (
                "alphabet",
                "data",
                "nr",
                "void",
                "nr_code",
                "void_code",
                "labels",
                "cpal",
                "ids",
            )
        )
    if isinstance(a, (pd.DataFrame, pd.Series)):
        return a.equals(b)
    if hasattr(a, "keys"):
        return all(_same(a[key], b[key]) for key in a.keys())
    if isinstance(a, (list, str)):
        return a == b
    return bool(np.array_equal(np.asarray(a), np.asarray(b), equal_nan=True))


_TWIN_CALLS = {
    **{
        f"{name}(include_missing={im})": (
            lambda o, name=name, im=im: getattr(tx, name)(o, include_missing=im)
        )
        for name in ("seqstatf", "seqtrate", "seqstatd", "seqmeant", "seqistatd")
        for im in (False, True)
    },
    "alphabet()": lambda o: tx.alphabet(o),
    "alphabet(with_missing=True)": lambda o: tx.alphabet(o, with_missing=True),
    "stlab": tx.stlab,
    "cpal": tx.cpal,
    "seqlength": tx.seqlength,
    # TX-REP-20261023-003 (batch-14 review M-1): the transition matrix, seqformat's outputs from an STSList, and
    # the spells of STS_to_SPELL with right, which re-run seqdef on the decoded states.
    "seqetm": tx.seqetm,
    **{
        f"seqformat(to={to!r})": (lambda o, to=to: tx.seqformat(o, to=to))
        for to in ("STS", "DSS", "SPS", "SRS", "SPELL")
    },
    "STS_to_SPELL(right='DEL')": lambda o: tx.STS_to_SPELL(o, right="DEL"),
}


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*'missing'\ set\ as\ .c\('\*','%'\).,\ the\ 'nr'\ and\ 'void'\ code\ from\ the\ 'data'\ state\ sequence\ object:UserWarning"
)
@pytest.mark.parametrize("call", sorted(_TWIN_CALLS))
def test_both_marker_representations_agree(call: str) -> None:
    x, canonical = _record_twins()
    assert _same(_TWIN_CALLS[call](x), _TWIN_CALLS[call](canonical))


def test_stslist_refuses_nr_text_equal_to_void_text() -> None:
    from traminex.sequences import STSList

    with pytest.raises(ValueError, match="'nr' and 'void' must be distinct character markers"):
        STSList(
            data=np.array([[0, 1, 2]]), alphabet=["A"], nr="m", void="m", nr_code=1, void_code=2
        )


# TX-AUD-20261019-062: R's seqdef runs as.matrix on a data frame, which pads a numeric column's numbers (" 1" beside
# "1") or makes NA cells and a duplicated alphabet (the golden). traminex reads each cell with R's as.character, so a
# frame mixing numeric and text columns gives the same object as its all-text version (a recorded correction of R).
@pytest.mark.parametrize("case", ["mixed_int_chr", "mixed_int_chr2", "mixed_dbl_chr"])
def test_seqdef_reads_mixed_columns_cell_by_cell(case: str) -> None:
    from traminex._rformat import r_as_character

    numbers = _numbers(case, "column_a")
    column_a = numbers if case == "mixed_dbl_chr" else numbers.astype(int)  # R's integer columns
    mixed = pd.DataFrame({"a": column_a, "b": _texts(case, "column_b")})
    as_text = pd.DataFrame(
        {"a": [r_as_character(v) for v in column_a], "b": _texts(case, "column_b")}
    )
    x, y = tx.seqdef(mixed), tx.seqdef(as_text)
    assert (x.alphabet, x.to_strings()) == (y.alphabet, y.to_strings())
    assert _texts(case, "alphabet") != tx.alphabet(x)  # R pads, duplicates or loses cells


def test_seqdef_ledger_row_records_the_mixed_column_reading() -> None:
    ledger = Path(__file__).resolve().parents[1] / "docs" / "API_PARITY_DECISIONS.csv"
    with open(ledger, encoding="utf-8", newline="") as handle:
        (row,) = [r for r in csv.DictReader(handle) if r["r_name"] == "seqdef"]
    assert "TX-AUD-20261019-062" in row["notes"]


# TX-AUD-20261019-070: seqfpos's per-state table (a traminex extension; R's seqfpos takes one state) lists the real
# states and, with include_missing, the missing state; never void, which is not part of a sequence (ruling R3-q).
# docs/API_CONSISTENCY.md states the three per-state width rules the ledger records.
def test_seqfpos_table_has_no_void_column() -> None:
    x = tx.seqdef([["b", "a", None, "c"], ["a", "a", "b", None], ["c", "b", "b", "b"]], right="DEL")
    assert x.void_code is not None and x.nr_code is not None
    table = tx.seqfpos(x, state=None, include_missing=True)
    assert list(table.columns) == ["a", "b", "c", "*"]
    assert table["*"].tolist() == [3, 0, 0]
    assert list(tx.seqfpos(x, state=None).columns) == ["a", "b", "c"]


# TX-REP-20261023-004 (batch-14 review M-2): seqfposend's table mode (a traminex extension, like seqfpos's) follows
# the same rule; with include_missing it had a void column, giving the position where a sequence ends.
def test_seqfposend_table_has_no_void_column() -> None:
    x = tx.seqdef([["b", "a", None, "c"], ["a", "a", "b", None], ["c", "b", "b", "b"]], right="DEL")
    table = tx.seqfposend(x, include_missing=True)
    assert list(table.columns) == ["a", "b", "c", "*"]
    assert table["*"].tolist() == [3, 0, 0]
    assert list(tx.seqfposend(x).columns) == ["a", "b", "c"]


def test_api_consistency_states_the_per_state_width_rules() -> None:
    doc = (Path(__file__).resolve().parents[1] / "docs" / "API_CONSISTENCY.md").read_text(
        encoding="utf-8"
    )
    text = " ".join(doc.split())
    assert "Per-**state** results are a different matter and are marker-padded" not in text
    for name in (
        "seqstatd",
        "seqmeant",
        "seqistatd",
        "seqstatf",
        "seqtrate",
        "seqfpos",
        "seqfposend",
    ):
        assert f"`{name}`" in text, name
    # The table modes are tested here, not in tests/test_per_state_shape_and_selection.py (review M-2).
    assert "`tests/test_construction_items.py`" in text


# TX-AUD-20261010-073: R's route for long data prints numeric states with as.character (1, 2 and 2.5; the golden), as
# traminex's seqdef does; the long builders gave '1.0' and '2.0'.
@pytest.mark.parametrize("case", ["float_whole", "float_fraction"])
def test_long_builders_read_numeric_states_as_r_prints_them(case: str) -> None:
    long = pd.DataFrame(
        {"id": [1, 1, 2, 2], "time": [1, 2, 1, 2], "state": _numbers(case, "state")}
    )
    expected = _texts(case, "alphabet")
    assert tx.alphabet(tx.seqdef_from_long(long)) == expected
    assert tx.alphabet(tx.TSE_to_STS(long)) == expected


# TX-AUD-20261011-060: with sm given as cost-method strings R's seqMD never copies a list indel into indel_list and
# stops ("subscript out of bounds"; seqMD.R:119-124, 226-264); with the same costs as matrices it computes. traminex
# computes the string route; its distances equal R's matrix route (the golden), and the seqMD row records it.
def test_seqmd_list_indels_with_cost_method_strings_equal_r_s_matrix_route() -> None:
    channels = [
        tx.seqdef([row.split("-") for row in _texts("md_indel", f"channel_{k}")], alphabet=alphabet)
        for k, alphabet in ((1, ["A", "B"]), (2, ["X", "Y"]))
    ]
    indel = [list(_numbers("md_indel", "indel_1")), list(_numbers("md_indel", "indel_2"))]
    assert GOLDEN[("md_indel", "string_route")] == "ERROR: subscript out of bounds"  # R
    # The per-state costs change R's distances (one indel per channel gives others), so the comparison below
    # shows that traminex uses them (batch-14 review M-3).
    assert not np.allclose(
        _numbers("md_indel", "matrix_route_diss"), _numbers("md_indel", "matrix_route_scalar_diss")
    )
    diss = tx.seqMD(channels, method="OM", sm=["CONSTANT", "CONSTANT"], indel=indel, what="diss")
    np.testing.assert_allclose(
        np.asarray(diss, dtype=float).ravel(order="F"),
        _numbers("md_indel", "matrix_route_diss"),
        rtol=0,
        atol=1e-12,
    )


def test_seqmd_ledger_row_records_list_indels_with_cost_method_strings() -> None:
    ledger = Path(__file__).resolve().parents[1] / "docs" / "API_PARITY_DECISIONS.csv"
    with open(ledger, encoding="utf-8", newline="") as handle:
        (row,) = [r for r in csv.DictReader(handle) if r["r_name"] == "seqMD"]
    assert "TX-AUD-20261011-060" in row["notes"]


# TX-AUD-20261011-062: with no rows R's min and max of the empty index are Inf and -Inf (seqstart.R:18-25), so a default
# tmax stops with "no data" and a default tmin leaves a frame matrix() refuses; a frame of negative length is refused
# too (seqstart.R:31), and one of length 0 is returned (the golden). traminex raises R's messages, not numpy's.
def _seqstart_rows(rows: int):
    from traminex.group_compare import _subset_sequence

    source = tx.seqdef([row.split() for row in _texts("seqstart", "cells")])
    return _subset_sequence(source, np.arange(rows))


def _bound(case: str, field: str) -> int | None:
    value = GOLDEN[(case, field)]
    return None if value == "NULL" else int(value)


@pytest.mark.parametrize("case", ["seqstart_empty", "seqstart_empty_tmax", "seqstart_inverted"])
def test_seqstart_stops_with_r_s_messages_on_an_empty_or_inverted_window(case: str) -> None:
    x = _seqstart_rows(int(GOLDEN[(case, "rows")]))
    message = GOLDEN[(case, "result")].removeprefix("ERROR: ").removeprefix("[!] ").rstrip(".")
    with pytest.raises(ValueError, match=re.escape(message)):
        tx.seqstart(x, 1, 1, tmin=_bound(case, "tmin"), tmax=_bound(case, "tmax"))


def test_seqstart_returns_a_zero_width_frame_as_r_does() -> None:
    case = "seqstart_zero_width"
    x = _seqstart_rows(int(GOLDEN[(case, "rows")]))
    result = tx.seqstart(x, 1, 1, tmin=_bound(case, "tmin"), tmax=_bound(case, "tmax"))
    assert GOLDEN[(case, "result")] == f"dim {result.nseq} {result.seqlen}"


def test_seqstart_ledger_row_records_its_stslist_input() -> None:
    ledger = Path(__file__).resolve().parents[1] / "docs" / "API_PARITY_DECISIONS.csv"
    with open(ledger, encoding="utf-8", newline="") as handle:
        (row,) = [r for r in csv.DictReader(handle) if r["r_name"] == "seqstart"]
    assert "TX-AUD-20261011-062" in row["notes"]
