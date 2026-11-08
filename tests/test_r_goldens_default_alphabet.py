"""The default alphabet is R's seqstatl order on every construction route (owner decision D2).

R goldens: ``scripts/r_goldens/make_goldens_default_alphabet.R`` (``default_alphabet_inputs.csv``,
``default_alphabet.csv``). Without ``alphabet=``, R's seqdef orders the states as seqstatl returns
them (seqdef.R:59-66; seqstatl.R:21-30): numbers by value when every value, a ``missing=`` code
included, is an R number, R's collation otherwise. Audit findings TX-AUD-20261019-071 (the
first-appearance default and its inconsistent application), -038 (seqMD's combined alphabet sorted
as text) and -013 (seqmodst and seqgranularity tie rules over the stored order).
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDEN_ATOL, R_GOLDENS_DIR

import traminex as tx

INPUTS = pd.read_csv(
    R_GOLDENS_DIR / "default_alphabet_inputs.csv",
    dtype=str,
    keep_default_na=False,
    encoding="utf-8",
)
RESULTS = pd.read_csv(
    R_GOLDENS_DIR / "default_alphabet.csv", dtype=str, keep_default_na=False, encoding="utf-8"
)


def _table(case: str, table: str) -> pd.DataFrame:
    rows = INPUTS[(INPUTS["case"] == case) & (INPUTS["table"] == table)]
    assert not rows.empty, (case, table)
    columns = list(dict.fromkeys(rows["column"]))
    frame = pd.DataFrame(
        {
            column: [
                None if value == "<NA>" else value
                for value in rows[rows["column"] == column].sort_values(
                    "row", key=lambda s: s.astype(int)
                )["value"]
            ]
            for column in columns
        }
    )
    return frame


def _r(case: str, route: str, kind: str = "alphabet") -> list[str]:
    rows = RESULTS[
        (RESULTS["case"] == case) & (RESULTS["route"] == route) & (RESULTS["kind"] == kind)
    ]
    assert not rows.empty, (case, route, kind)
    return rows.sort_values("position", key=lambda s: s.astype(int))["value"].tolist()


def _quiet(function, *args, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return function(*args, **kwargs)


STS_CASES = {
    "numeric": {},
    "text_case": {},
    "text_accents": {},
    "text_punct": {},
    "mixed_number_text": {},
    "missing_text_code": {"missing": "x"},
    "missing_numeric_code": {"missing": "99"},
    "left_policy": {"left": "LL"},
    "right_na": {"right": None},
    "census_wide": {},
    "events_grid": {},
}


@pytest.mark.parametrize("case", sorted(STS_CASES))
def test_seqdef_sts_default_alphabet_is_r_order(case: str) -> None:
    x = _quiet(tx.seqdef, _table(case, "sts"), **STS_CASES[case])
    assert tx.alphabet(x) == _r(case, "seqdef")


def test_logical_column_is_read_as_a_number() -> None:
    frame = _table("logical_promotion", "frame")
    frame = pd.DataFrame(
        {
            "logical": frame["logical"].map({"TRUE": True, "FALSE": False}),
            "integer": frame["integer"].astype(int),
        }
    )
    assert tx.alphabet(_quiet(tx.seqdef, frame)) == _r("logical_promotion", "seqdef")


def test_seqdef_and_seqformat_sps_routes() -> None:
    sps = _table("sps", "sps")
    assert tx.alphabet(_quiet(tx.seqdef, sps, informat="SPS")) == _r("sps", "seqdef_SPS")
    formatted = _quiet(tx.seqformat, sps, from_="SPS", to="STS")
    assert tx.alphabet(formatted) == _r("sps", "seqformat_SPS_STS")


def _spell(case: str) -> pd.DataFrame:
    frame = _table(case, "spell")
    for column in ("id", "begin", "end"):
        frame[column] = frame[column].astype(float)
    return frame


def test_seqdef_spell_route_text_and_factor_status() -> None:
    text = _spell("spell_text")
    assert tx.alphabet(_quiet(tx.seqdef, text, informat="SPELL", process=False)) == _r(
        "spell_text", "seqdef_SPELL"
    )
    factor = _spell("spell_factor")
    levels = _table("spell_factor", "levels")["level"].tolist()
    factor["status"] = pd.Categorical(factor["status"], categories=levels)
    # R's seqdef sorts the cells as seqstatl does and ignores the factor levels.
    assert tx.alphabet(_quiet(tx.seqdef, factor, informat="SPELL", process=False)) == _r(
        "spell_factor", "seqdef_SPELL"
    )


def test_seqformat_spell_factor_keeps_levels() -> None:
    """seqformat's own output keeps a factor status's levels, as R's SPELL_to_STS columns do
    (seqformat-SPELL_to_STS.R:347-353); plain text status takes R's seqdef order."""
    factor = _spell("spell_factor")
    factor["status"] = pd.Categorical(
        factor["status"], categories=_table("spell_factor", "levels")["level"].tolist()
    )
    formatted = _quiet(tx.seqformat, factor, from_="SPELL", to="STS", process=False)
    assert tx.alphabet(formatted) == _r("spell_factor", "seqformat_SPELL_STS", "levels")
    text = _quiet(tx.seqformat, _spell("spell_text"), from_="SPELL", to="STS", process=False)
    assert tx.alphabet(text) == _r("spell_text", "seqdef_SPELL")


def test_long_builders_default_alphabet() -> None:
    long = _table("long", "long")
    long["id"] = long["id"].astype(int)
    long["time"] = long["time"].astype(int)
    expected = _r("long", "seqdef_reshape")
    assert tx.alphabet(_quiet(tx.seqdef_from_long, long)) == expected
    assert tx.alphabet(_quiet(tx.TSE_to_STS, long)) == expected


@pytest.mark.parametrize("case, lengths", [("long_star", [3, 3]), ("long_star_gap", [3, 2])])
def test_long_builders_keep_a_star_state_unlike_r_reshape(case: str, lengths: list[int]) -> None:
    # Batch-7 review M-3, ruling R7-e: R's reshape-then-seqdef lists * as a state yet counts
    # its cells as missing; the long builders keep * a real state. Lengths counted by hand:
    # person 1 is A, *, B; person 2 is A, A, B (long_star) or A, missing, B (long_star_gap).
    long = _table(case, "long")
    long["id"] = long["id"].astype(int)
    long["time"] = long["time"].astype(int)
    r_lengths = [int(value) for value in _r(case, "seqdef_reshape", "seqlength_without_missing")]
    for build in (tx.seqdef_from_long, tx.TSE_to_STS):
        x = _quiet(build, long)
        assert tx.alphabet(x) == _r(case, "seqdef_reshape")
        actual = np.asarray(tx.seqlength(x, with_missing=False)).ravel().tolist()
        assert actual == lengths
        assert actual != r_lengths


def test_tse_to_sts_event_mode_default_alphabet() -> None:
    events = _table("tse_event", "events")
    events["id"] = events["id"].astype(int)
    events["time"] = events["time"].astype(float)
    assert tx.alphabet(_quiet(tx.TSE_to_STS, events, tmax=3)) == _r(
        "tse_event", "TSE_to_STS_seqdef"
    )


def test_events_to_state_sequences_default_order() -> None:
    """The grid of p1 Zeta@0, alpha@1, Zeta@2 and p2 alpha@1 with fill "_" is the input table of
    case events_grid; its states take seqdef's order (code points would give Zeta, _, alpha)."""
    events = tx.EventSequence(
        id=["p1", "p1", "p1", "p2"],
        time=np.array([0.0, 1.0, 2.0, 1.0]),
        label=["Zeta", "alpha", "Zeta", "alpha"],
    )
    x = _quiet(tx.events_to_state_sequences, events, fill="_")
    assert x.to_strings() == ["Zeta-alpha-Zeta", "_-alpha-alpha"]
    assert tx.alphabet(x) == _r("events_grid", "seqdef")


def test_seqrecode_and_seqdss_follow_the_input_order() -> None:
    census = _quiet(tx.seqdef, _table("census_wide", "sts"))
    assert tx.alphabet(_quiet(tx.seqrecode, census, {"Z": ["c"]})) == _r(
        "census_wide", "seqrecode_Z_c"
    )
    assert tx.alphabet(_quiet(tx.seqdss, census)) == _r("census_wide", "seqdss")


def _channel(case: str, table: str) -> tx.STSList:
    return _quiet(tx.seqdef, _table(case, table).to_numpy(dtype=object))


def test_seqmd_numeric_alphabet_038() -> None:
    channels = [_channel("seqmd_numeric", "channel1"), _channel("seqmd_numeric", "channel2")]
    combined = _quiet(tx.seqMD, channels, what="MDseq", ch_sep=".")
    assert tx.alphabet(combined) == _r("seqmd_numeric", "seqMD_MDseq")
    cost = _quiet(tx.seqMD, channels, what="cost", sm="CONSTANT", ch_sep=".")
    assert list(cost["alphabet"]) == _r("seqmd_numeric", "seqMD_cost", "rownames")


def test_seqmd_text_alphabet_uses_r_collation() -> None:
    channels = [_channel("seqmd_text", "channel1"), _channel("seqmd_text", "channel2")]
    combined = _quiet(tx.seqMD, channels, what="MDseq")
    assert tx.alphabet(combined) == _r("seqmd_text", "seqMD_MDseq")


def test_modal_ties_013() -> None:
    """seqmodst and seqgranularity(mostfreq) give a tie to the first tied state of the alphabet
    (seqmodst.R:28-32); with R's default order they match R."""
    ties = _quiet(tx.seqdef, _table("ties", "sts"))
    assert tx.alphabet(ties) == _r("ties", "seqdef")
    assert _quiet(tx.seqmodst, ties).to_strings() == _r("ties", "seqmodst", "string")
    granular = _quiet(tx.seqgranularity, ties, tspan=2, method="mostfreq")
    assert granular.to_strings() == _r("ties", "seqgranularity_mostfreq_2", "string")
    reversed_order = _quiet(tx.seqdef, _table("ties", "sts"), alphabet=["b", "a"])
    assert _quiet(tx.seqmodst, reversed_order).to_strings() == ["b-b-b"]


def test_order_dependent_defaults_match_r() -> None:
    x = _quiet(tx.seqdef, _table("indicators", "sts"))
    assert tx.alphabet(x) == _r("indicators", "seqdef")
    for name in ("seqibad", "seqprecstart", "seqidegrad", "seqprecarity", "seqinsecurity"):
        actual = np.asarray(_quiet(getattr(tx, name), x), dtype=float).ravel()
        expected = np.asarray(_r("indicators", name, "value"), dtype=float)
        np.testing.assert_allclose(actual, expected, rtol=0, atol=R_GOLDEN_ATOL, err_msg=name)
    assert _quiet(tx.seqmodst, x).to_strings() == _r("indicators", "seqmodst", "string")
    granular = _quiet(tx.seqgranularity, x, tspan=2, method="mostfreq")
    assert granular.to_strings() == _r("indicators", "seqgranularity_mostfreq_2", "string")
    r_rows = [
        label.removeprefix("[").removesuffix(" ->]")
        for label in _r("indicators", "seqtrate", "rownames")
    ]
    assert list(_quiet(tx.seqtrate, x).index) == r_rows


@pytest.mark.parametrize(
    "build",
    [
        lambda rows, alphabet: tx.seqdef(rows, alphabet=alphabet),
        lambda rows, alphabet: tx.seqdef(
            pd.DataFrame({"x": ["(c,1)-(a,1)-(b,1)", "(b,2)-(a,1)"]}),
            informat="SPS",
            alphabet=alphabet,
        ),
        lambda rows, alphabet: tx.seqdef(
            pd.DataFrame(
                {
                    "id": [1, 1, 1, 2, 2],
                    "begin": [1, 2, 3, 1, 3],
                    "end": [1, 2, 3, 2, 3],
                    "status": ["c", "a", "b", "b", "a"],
                }
            ),
            informat="SPELL",
            process=False,
            alphabet=alphabet,
        ),
    ],
    ids=["STS", "SPS", "SPELL"],
)
def test_explicit_alphabet_keeps_caller_order(build) -> None:
    rows = [["c", "a", "b"], ["b", "b", "a"]]
    x = _quiet(build, rows, ["c", "a", "b"])
    assert tx.alphabet(x) == ["c", "a", "b"]
