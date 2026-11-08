"""Batch-15 sequence-statistics items against R (``make_goldens_statistics_items.R``).

The golden records what TraMineR and TraMineRextras return for the inputs these tests use, R's errors included where
traminex stops as R does or corrects R; each test states whether it expects R's value or a value derived from a stated
rule.
"""

from __future__ import annotations

import csv
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx

GOLDEN_PATH = Path(__file__).resolve().parent / "r_goldens" / "statistics_items.csv"
with open(GOLDEN_PATH, encoding="utf-8", newline="") as handle:
    GOLDEN = {(row["case"], row["field"]): row["value"] for row in csv.DictReader(handle)}


def _numbers(case: str, field: str) -> np.ndarray:
    return np.array([np.nan if v == "NA" else float(v) for v in GOLDEN[(case, field)].split()])


def _texts(case: str, field: str) -> list[str]:
    return GOLDEN[(case, field)].split(" ;; ")


# TX-AUD-20261010-062: seqdomassoc.Rd takes assoc as a character string or a vector with both measures; a string was
# read character by character and refused. R's values for assoc = "LRT" (the golden).
def _domains() -> list:
    return [
        tx.seqdef([row.split("-") for row in _texts("domassoc", f"domain_{k}")]) for k in (1, 2)
    ]


def test_seqdomassoc_takes_assoc_as_one_string() -> None:
    result = tx.seqdomassoc(_domains(), assoc="LRT")
    assert list(result.columns) == _texts("domassoc", "columns")
    np.testing.assert_allclose(
        result.to_numpy(dtype=float).ravel(), _numbers("domassoc", "values"), rtol=0, atol=1e-12
    )
    pd.testing.assert_frame_equal(result, tx.seqdomassoc(_domains(), assoc=["LRT"]))


# TX-AUD-20261010-095: Shannon entropy by its definition, a zero frequency contributing 0 (R's own entropy.R:5-7 drops
# them). entropy took log(p + 1e-12), so uniform data gave 0.9999999999971146 normalized; the same bias sat in
# seqtree_explore's summary entropy and local_indel_by_pos(method="entropy"). The missing state is one of the K states.
@pytest.mark.parametrize("k", [2, 4])
def test_entropy_of_uniform_data_is_exactly_its_maximum(k: int) -> None:
    from traminex.stats_basic import entropy

    states = list("ABCD"[:k])
    x = tx.seqdef([states, states[::-1]])
    assert abs(entropy(x) - 1.0) <= 1e-15
    assert abs(entropy(x, normalized=False) - math.log(k)) <= 1e-15


def test_entropy_counts_the_missing_state_as_documented() -> None:
    from traminex.stats_basic import entropy

    # Two A, two B and two missing cells: uniform over K = 3 states.
    x = tx.seqdef([["A", None, "B"], ["B", "A", None]], right=None)
    assert abs(entropy(x) - 1.0) <= 1e-15
    assert "missing state" in " ".join((entropy.__doc__ or "").split())


def test_seqtree_explore_summary_entropy_of_uniform_rows_is_one() -> None:
    from traminex.seqtree_explore import _features_from_sequences

    summary, names = _features_from_sequences(tx.seqdef([["A", "B"], ["B", "A"]]), kind="summary")
    assert names[0] == "entropy_norm"
    assert np.all(np.abs(summary[:, 0] - 1.0) <= 1e-15)


def test_local_indel_entropy_of_uniform_positions() -> None:
    from traminex.distances.costs import local_indel_by_pos

    # Each position holds one A and one B: normalized entropy 1, indel 1 - 0.5 * 1.
    indel = local_indel_by_pos(
        tx.seqdef([["A", "B"], ["B", "A"]]).data, method="entropy", scale=0.5
    )
    assert np.all(np.abs(indel - 0.5) <= 1e-15)


# TX-AUD-20261010-096: a group mixing text and numbers has no R counterpart (R's c() makes it text before seqdiff sees
# it); seqdiff refuses it as traminex's other group functions do (r_factor_levels, TX-AUD-20261010-039) instead of
# failing while sorting. On no sequences R's seqindic.dyn stops (the golden) and seqindic_dyn stops with R's message.
def test_seqdiff_refuses_a_group_mixing_text_and_numbers() -> None:
    x = tx.seqdef([["A", "B", "A"], ["B", "B", "A"], ["A", "A", "B"], ["B", "A", "A"]])
    with pytest.raises(TypeError, match="group levels mix text and numbers"):
        tx.seqdiff(x, np.array(["a", 1, "a", 1], dtype=object))


def test_seqindic_dyn_stops_on_no_sequences_as_r_does() -> None:
    from traminex.group_compare import _subset_sequence

    x = tx.seqdef([row.split("-") for row in _texts("indic_dyn_empty", "cells")])
    message = GOLDEN[("indic_dyn_empty", "result")].removeprefix("ERROR: ")
    with pytest.raises(ValueError, match=re.escape(message)):
        tx.seqindic_dyn(_subset_sequence(x, np.array([], dtype=int)))


# TX-AUD-20261011-001: R's seqgen returns the alphabet's own type (sample() keeps it; the golden's typeof); seqgen gave
# text. R reads one number x >= 1 as 1:x, so alphabet = 5 draws 1..5 (the golden), against seqgen.Rd ("the alphabet
# from which the sequences are generated"): traminex draws 5, a recorded correction.
_R_TYPES = {"integer": "i", "double": "f", "character": "O"}  # text: an object array of str


@pytest.mark.parametrize(
    ("field", "alphabet"),
    [
        ("typeof_integer", [1, 2, 3]),
        ("typeof_double", [1.0, 2.0, 3.0]),
        ("typeof_character", ["A", "B"]),
    ],
)
def test_seqgen_keeps_the_alphabet_s_type(field: str, alphabet: list) -> None:
    drawn = tx.seqgen(2, 3, alphabet, seed=1)
    assert drawn.dtype.kind == _R_TYPES[GOLDEN[("seqgen_type", field)]]
    assert all(isinstance(value, type(alphabet[0])) for value in drawn.ravel().tolist())
    assert set(drawn.ravel().tolist()) <= set(alphabet)


@pytest.mark.parametrize("alphabet", [[5], 5])
def test_seqgen_draws_a_one_value_alphabet_s_value(alphabet: object) -> None:
    drawn = tx.seqgen(50, 20, alphabet, seed=1)
    assert drawn.dtype.kind == "i" and set(drawn.ravel().tolist()) == {5}
    assert _numbers("seqgen_type", "one_value_5").tolist() == [1, 2, 3, 4, 5]  # R


# TX-AUD-20261011-072: seqgen passes p to R's sample(), which reads probability weights: under one seed c(5, 3, 2)
# draws what c(.5, .3, .2) draws (the golden), and a weight that is negative, missing or infinite, or no positive weight,
# stops with sample()'s message (the golden). seqgen refused weights that did not sum to 1.
def test_seqgen_reads_p_as_probability_weights() -> None:
    assert GOLDEN[("seqgen_weights", "same_draws")] == "TRUE"  # R
    weights = tx.seqgen(50, 10, ["A", "B", "C"], p=[5, 3, 2], seed=7)
    probabilities = tx.seqgen(50, 10, ["A", "B", "C"], p=[0.5, 0.3, 0.2], seed=7)
    np.testing.assert_array_equal(weights, probabilities)
    assert "B" not in set(tx.seqgen(50, 10, ["A", "B", "C"], p=[1, 0, 1], seed=7).ravel().tolist())


@pytest.mark.parametrize(
    "case", ["negative", "missing", "infinite", "zero", "negative_then_missing"]
)
def test_seqgen_refuses_weights_as_r_does(case: str) -> None:
    message = GOLDEN[("seqgen_weights", f"result_{case}")].removeprefix("ERROR: ")
    with pytest.raises(ValueError, match=re.escape(message)):
        tx.seqgen(2, 3, ["A", "B", "C"], p=_numbers("seqgen_weights", f"p_{case}").tolist(), seed=1)


# TX-AUD-20261011-002: R's seqpcfilter matches the method with match.arg: an exact name or a unique prefix,
# case-sensitive (the golden); seqpcfilter took exact names only. R checks level against the unmatched name
# (seqpcplot.R:1101), so ("cum", NULL) escapes the check "cumfreq" gets (the golden); traminex checks the matched name,
# a recorded correction.
def _pcfilter_call(case: str) -> tuple[str, float | None]:
    level = GOLDEN[(f"pcfilter_{case}", "level")]
    return GOLDEN[(f"pcfilter_{case}", "method")], None if level == "NULL" else float(level)


@pytest.mark.parametrize("case", ["cum", "m", "lin", "cumfreq"])
def test_seqpcfilter_matches_the_method_as_match_arg(case: str) -> None:
    result = tx.seqpcfilter(*_pcfilter_call(case))
    assert [result.type, result.value] == GOLDEN[(f"pcfilter_{case}", "result")].split()[:2]


@pytest.mark.parametrize("case", ["empty", "x", "upper"])
def test_seqpcfilter_refuses_what_match_arg_refuses(case: str) -> None:
    assert GOLDEN[(f"pcfilter_{case}", "result")].startswith("ERROR: 'arg' should be one of")  # R
    with pytest.raises(ValueError, match="method must be one of minfreq, cumfreq, linear"):
        tx.seqpcfilter(*_pcfilter_call(case))


def test_seqpcfilter_checks_the_level_of_the_matched_method() -> None:
    assert GOLDEN[("pcfilter_cum_null", "result")] == "function cumfreq NULL"  # R
    with pytest.raises(ValueError, match="requires a valid 'level' value"):
        tx.seqpcfilter(*_pcfilter_call("cum_null"))


# TX-AUD-20261011-068: seqtransn_count's R analogue is seqtransn: missing cells are left out of the spells and the
# count is divided by each sequence's length (missing cells counted, void not) minus 1 (the golden, three objects).
# seqtransn_count counted changes into and out of the missing marker and divided by the padded width.
@pytest.mark.parametrize("case", ["transn_ragged", "transn_missing_row", "transn_length_one"])
def test_seqtransn_count_counts_and_normalizes_as_seqtransn(case: str) -> None:
    from traminex.stats_complexity import seqtransn_count

    rows = [[None if v == "NA" else v for v in row.split()] for row in _texts(case, "cells")]
    right = GOLDEN[(case, "right")]
    states = sorted({v for row in rows for v in row if v is not None})
    x = tx.seqdef(rows, alphabet=states, right=None if right == "NA" else right)
    assert [s.replace("-%", "") for s in x.to_strings()] == _texts(case, "seqconc")  # R's object
    np.testing.assert_allclose(
        seqtransn_count(x, normalized=False), _numbers(case, "transn"), rtol=0, atol=1e-15
    )
    np.testing.assert_allclose(
        seqtransn_count(x), _numbers(case, "transn_norm"), rtol=0, atol=1e-15
    )


# TX-AUD-20261011-069: R's seqsha stops on a missing time or event with its message and on a time of 0 or below (the
# golden); seqsha.Rd says event is logical, and R keeps text events as a text column (the golden). seqsha computed with
# time 0 (person 1's event written at index -1), raised numpy's errors, and read "no" and NaN events as True.
def _sha_call(case: str) -> tuple[list, list]:
    events = _texts(case, "event")
    if case == "sha_event_text":
        parsed: list = events
    elif case == "sha_event_numeric":
        parsed = [float(v) for v in events]
    else:
        parsed = [{"TRUE": True, "FALSE": False, "NA": None}[v] for v in events]
    return _numbers(case, "time").tolist(), parsed


def _sha(case: str):
    x = tx.seqdef([row.split("-") for row in _texts("sha", "cells")])
    time, event = _sha_call(case)
    return tx.seqsha(x, time=time, event=event)


@pytest.mark.parametrize("case", ["sha_control", "sha_event_numeric"])
def test_seqsha_reads_logical_events_as_r_does(case: str) -> None:
    rows, _, values = GOLDEN[(case, "result")].split(" | ")
    result = _sha(case)
    assert f"rows {len(result)}" == rows
    expected = [{"TRUE": 1, "FALSE": 0, "1": 1, "0": 0}[v] for v in values.split(" ;; ")]
    assert result["event"].astype(int).tolist() == expected


@pytest.mark.parametrize("case", ["sha_time_missing", "sha_event_missing"])
def test_seqsha_stops_on_missing_values_with_r_s_message(case: str) -> None:
    message = GOLDEN[(case, "result")].removeprefix("ERROR:").strip()
    with pytest.raises(ValueError, match=re.escape(message)):
        _sha(case)


@pytest.mark.parametrize("case", ["sha_time_zero", "sha_time_negative"])
def test_seqsha_stops_on_a_time_below_one_as_r_does(case: str) -> None:
    assert GOLDEN[(case, "result")].startswith("ERROR:")  # R
    with pytest.raises(ValueError, match="time must be at least 1"):
        _sha(case)


def test_seqsha_refuses_text_events() -> None:
    assert "| event character |" in GOLDEN[("sha_event_text", "result")]  # R keeps the text
    with pytest.raises(ValueError, match="event must be logical"):
        _sha("sha_event_text")


# TX-AUD-20261011-070: R stops on valid input where its own definitions give a value (the golden records R's four
# errors): a sequence without transitions keeps seqidegrad's penalty 0 (seqidegrad.R:310-332); the share of positive
# states of an all-A sequence with pos A is 1 (seqipos.Rd); with with.missing = TRUE the one-state precarity is
# otto * 0 + (1 - otto) * 0 = 0 (seqprecstart gives A 0, the complexity is 0; seqprecarity.R:137); and without a
# missing cell, with.missing adds a state of frequency 0, so seqlogp equals its with.missing = FALSE value (the
# golden). traminex returns these values; the four ledger rows record the correction.
def _degenerate_objects():
    one = tx.seqdef([row.split("-") for row in _texts("degenerate", "one_state")], alphabet=["A"])
    two = tx.seqdef(
        [row.split("-") for row in _texts("degenerate", "two_state")], alphabet=["A", "B"]
    )
    return one, two


def test_degenerate_inputs_take_the_values_r_s_definitions_give() -> None:
    for field in ("seqidegrad", "seqipos", "seqprecarity_with_missing", "seqlogp_with_missing"):
        assert GOLDEN[("degenerate", field)].startswith("ERROR:"), field  # R
    one, two = _degenerate_objects()
    np.testing.assert_array_equal(np.asarray(tx.seqidegrad(one), dtype=float).ravel(), [0.0, 0.0])
    np.testing.assert_array_equal(
        np.asarray(tx.seqipos(one, pos_states=["A"]), dtype=float).ravel(), [1.0, 1.0]
    )
    np.testing.assert_array_equal(
        np.asarray(tx.seqprecarity(one, include_missing=True), dtype=float).ravel(), [0.0, 0.0]
    )
    np.testing.assert_allclose(
        np.asarray(
            tx.seqlogp(two, prob="freq", begin="global.freq", include_missing=True), dtype=float
        ).ravel(),
        _numbers("degenerate", "seqlogp_without_missing"),
        rtol=0,
        atol=1e-12,
    )


@pytest.mark.parametrize("r_name", ["seqidegrad", "seqipos", "seqprecarity", "seqlogp"])
def test_degenerate_ledger_rows_record_where_r_stops(r_name: str) -> None:
    ledger = Path(__file__).resolve().parents[1] / "docs" / "API_PARITY_DECISIONS.csv"
    with open(ledger, encoding="utf-8", newline="") as handle:
        (row,) = [r for r in csv.DictReader(handle) if r["r_name"] == r_name]
    assert "TX-AUD-20261011-070" in row["notes"]


# TX-AUD-20261011-073: seqimplic.Rd gives group as the membership of each sequence. R stops on a longer group, and on a
# shorter one recycles its row mask and returns indices without one finite value (the golden). seqimplic raised a raw
# IndexError on a shorter group and computed silently on a longer one; both are refused.
def _implic_object():
    return tx.seqdef([row.split("-") for row in _texts("implic", "cells")])


def test_seqimplic_with_one_group_entry_per_sequence_matches_r() -> None:
    result = tx.seqimplic(_implic_object(), _texts("implic", "group"))
    np.testing.assert_allclose(
        np.asarray(result.indices, dtype=float).ravel(order="F"),
        _numbers("implic_full", "indices"),
        rtol=0,
        atol=1e-12,
    )


@pytest.mark.parametrize("case", ["implic_short", "implic_long"])
def test_seqimplic_refuses_a_group_of_the_wrong_length(case: str) -> None:
    r_result = GOLDEN[(case, "result")]
    assert r_result.endswith("finite 0") or r_result.startswith(
        "ERROR:"
    )  # R: nothing finite, or an error
    group = _texts("implic", "group")
    group = group[:5] if case == "implic_short" else [*group, "g"]
    with pytest.raises(ValueError, match="group must have one entry per sequence"):
        tx.seqimplic(_implic_object(), group)


# TX-AUD-20261011-074: seqalign's indices are two zero-based positions (its docstring); R's are one-based numbers that
# R coerces (TRUE is 1, 1.7 is 1). Booleans, fractions and a wrong number of positions raised raw TypeError,
# IndexError or unpacking errors; they are refused with a TypeError naming indices. NumPy integers stay accepted.
def _three_sequences():
    return tx.seqdef([["A", "B", "A"], ["B", "B", "A"], ["A", "A", "B"]], alphabet=["A", "B"])


@pytest.mark.parametrize(
    "indices", [(True, 2), (0.5, 1), (0,), (0, 1, 2)], ids=["bool", "float", "one", "three"]
)
def test_seqalign_refuses_indices_that_are_not_two_integer_positions(indices: tuple) -> None:
    with pytest.raises(TypeError, match="indices must be two integer positions"):
        tx.seqalign(_three_sequences(), indices)


def test_seqalign_accepts_numpy_integer_positions() -> None:
    x = _three_sequences()
    assert (
        tx.seqalign(x, (np.int64(0), np.int32(1))).total_cost == tx.seqalign(x, (0, 1)).total_cost
    )


# TX-REP-20261019-011 (the 2026-10-19 audit's oracle map): the ledger rows of STS_to_SRS, seqibad, seqpcfilter,
# dissvar.grp and seqsha cited no R generator that calls the R function. Each row cites one (STS_to_SRS's through
# seqformat(to = "SRS"), which calls it at seqformat.R:246), and cites the generator writing the fixtures its tests read.
_R_CALLS = {
    "STS_to_SRS": r'seqformat\([^)]*to\s*=\s*"SRS"',
    "seqibad": r"seqibad\s*\(",
    "seqpcfilter": r"seqpcfilter\s*\(",
    "dissvar.grp": r"dissvar\.grp\s*\(",
    "seqsha": r"seqsha\s*\(",
}


def _ledger_rows() -> dict[str, dict[str, str]]:
    ledger = Path(__file__).resolve().parents[1] / "docs" / "API_PARITY_DECISIONS.csv"
    with open(ledger, encoding="utf-8", newline="") as handle:
        return {row["r_name"]: row for row in csv.DictReader(handle)}


@pytest.mark.parametrize("r_name", sorted(_R_CALLS))
def test_ledger_row_cites_a_generator_that_calls_the_r_function(r_name: str) -> None:
    root = Path(__file__).resolve().parents[1]
    cited = [
        e
        for e in _ledger_rows()[r_name]["evidence"].split(";")
        if e.startswith("scripts/r_goldens/")
    ]

    def calls(path: str) -> bool:
        lines = (root / path).read_text(encoding="utf-8").splitlines()
        code = "\n".join(line for line in lines if not line.lstrip().startswith("#"))
        return re.search(_R_CALLS[r_name], code) is not None

    assert any(calls(path) for path in cited), cited


def test_ledger_rows_cite_the_generators_of_their_fixtures() -> None:
    rows = _ledger_rows()
    assert "scripts/r_goldens/make_goldens_seqsha.R" in rows["seqsha"]["evidence"].split(";")
    assert "scripts/r_goldens/make_goldens_seqpcfilter.R" in rows["seqpcfilter"]["evidence"].split(
        ";"
    )
    assert 'seqformat(to = "SRS")' in rows["STS_to_SRS"]["notes"]


# Batch-15 review I-1 (TX-AUD-20261010-096's ZeroDivisionError sub-claim, first recorded as not reproduced because the
# attempts did not call seqdiff): a window without a complete case stops R's seqdiff in seqdist's matrix() (the
# golden), and with all-zero weights R's dissassoc divides 0 by 0 and returns NaN (the golden). seqdiff and dissassoc
# raised Python's ZeroDivisionError.
def _no_complete_case(case: str):
    from traminex.group_compare import _subset_sequence

    if case == "no_rows":
        base = tx.seqdef([["A", "B", "A"], ["B", "B", "A"]])
        return _subset_sequence(base, np.array([], dtype=int)), np.array([], dtype=object)
    if case == "all_void_row":
        return tx.seqdef([["A", None], [None, None]]), np.array(["g", "h"], dtype=object)
    return tx.seqdef([[None, None], [None, None]], right=None, alphabet=["A"]), np.array(
        ["g", "h"], dtype=object
    )


@pytest.mark.parametrize("case", ["no_rows", "all_void_row", "all_missing"])
def test_seqdiff_stops_on_a_window_without_a_complete_case_as_r_does(case: str) -> None:
    message = GOLDEN[(f"seqdiff_{case}", "result")].removeprefix("ERROR: ")
    x, group = _no_complete_case(case)
    with pytest.raises(ValueError, match=re.escape(message)):
        tx.seqdiff(x, group)


def test_dissassoc_with_zero_weights_gives_r_s_values() -> None:
    diss = _numbers("diss_zero_weights", "diss").reshape(2, 2)
    result = tx.dissassoc(
        diss,
        np.array(_texts("diss_zero_weights", "group")),
        weights=[0.0, 0.0],
        R=0,
        weight_permutation="diss",
    )
    for column, field in (("n", "groups_n"), ("discrepancy", "groups_discrepancy")):
        np.testing.assert_array_equal(
            result.groups[column].to_numpy(dtype=float), _numbers("diss_zero_weights", field)
        )
    np.testing.assert_array_equal(
        result.stat["t0"].to_numpy(dtype=float), _numbers("diss_zero_weights", "stat")
    )
    for column, field in (("SS", "anova_ss"), ("df", "anova_df"), ("MSE", "anova_mse")):
        np.testing.assert_array_equal(
            result.anova_table[column].to_numpy(dtype=float), _numbers("diss_zero_weights", field)
        )


def test_seqdiff_with_zero_weights_gives_r_s_values() -> None:
    x = tx.seqdef(
        [row.split("-") for row in _texts("diss_zero_weights", "cells")], weights=[0.0, 0.0]
    )
    result = tx.seqdiff(x, np.array(_texts("diss_zero_weights", "group"), dtype=object))
    np.testing.assert_array_equal(
        result["stat"].to_numpy(dtype=float).ravel(order="F"),
        _numbers("diss_zero_weights", "seqdiff_stat"),
    )
    np.testing.assert_array_equal(
        result["discrepancy"].to_numpy(dtype=float).ravel(order="F"),
        _numbers("diss_zero_weights", "seqdiff_discrepancy"),
    )


# TX-REP-20261023-005 (batch-15 review I-1): with no sequence, or only missing cells, no situation is observed and R's
# seqemlt stops (the golden); seqemlt raised Python's ZeroDivisionError. It stops with R's message.
@pytest.mark.parametrize("case", ["no_rows", "all_missing"])
def test_seqemlt_stops_without_an_observed_state_as_r_does(case: str) -> None:
    message = GOLDEN[(f"seqemlt_{case}", "result")].removeprefix("ERROR: ")
    x, _ = _no_complete_case(case)
    with pytest.raises(ValueError, match=re.escape(message)):
        tx.seqemlt(x)


# Batch-15 review M-1 (TX-AUD-20261011-069): R's seqsha also stops on a fractional or an infinite time (the golden);
# seqsha truncated 2.5 to 2 without a message and gave numpy's error for an infinite time.
@pytest.mark.parametrize("case", ["sha_time_fraction", "sha_time_infinite"])
def test_seqsha_stops_on_a_time_that_is_not_a_whole_number_as_r_does(case: str) -> None:
    assert GOLDEN[(case, "result")].startswith("ERROR:")  # R
    with pytest.raises(ValueError, match="time must be a whole number"):
        _sha(case)


# Batch-15 review M-2 (the root cause of TX-AUD-20261010-095): JSD_TRANS took ln(a + 1e-12) - ln(b + 1e-12) by default,
# though only positive shares enter the logarithm (a > 0, b = (a + c)/2 > 0); two sequences with disjoint transitions
# gave ln 2 - 2.0e-12. The Jensen-Shannon divergence of disjoint distributions is ln 2 (its definition).
def test_jsd_trans_of_disjoint_transitions_is_ln_2() -> None:
    x = tx.seqdef([["A", "A", "A"], ["B", "B", "B"]], alphabet=["A", "B"])
    diss = np.asarray(tx.seqdist(x, method="JSD_TRANS"), dtype=float)
    assert abs(diss[0, 1] - math.log(2)) <= 1e-15


# Batch-15 review M-3 (TX-AUD-20261011-074): the index check of T11 refused two-element containers the earlier
# unpacking accepted (range, a pandas Series); any non-text sequence of two integer positions is accepted again.
@pytest.mark.parametrize(
    "indices", [range(2), pd.Series([0, 1]), np.array([0, 1]), [0, 1]], ids=str
)
def test_seqalign_accepts_any_sequence_of_two_integer_positions(indices: object) -> None:
    x = _three_sequences()
    assert tx.seqalign(x, indices).total_cost == tx.seqalign(x, (0, 1)).total_cost


def test_seqalign_refuses_text_indices() -> None:
    with pytest.raises(TypeError, match="indices must be two integer positions"):
        tx.seqalign(_three_sequences(), "01")


# Batch-15 review M-4 (TX-AUD-20261010-096): a group of booleans and numbers was refused with "mix text and numbers",
# though it holds no text; the refusal names the mix it finds.
def test_group_refusal_names_the_mix_it_finds() -> None:
    from traminex._categorical_order import r_factor_levels

    x = tx.seqdef([["A", "B", "A"], ["B", "B", "A"], ["A", "A", "B"], ["B", "A", "A"]])
    with pytest.raises(TypeError, match="group levels mix numbers and booleans"):
        tx.seqdiff(x, np.array([True, 2, True, 2], dtype=object))
    with pytest.raises(TypeError, match="group levels mix text, numbers and booleans"):
        r_factor_levels(["a", 1, True])
    with pytest.raises(TypeError, match="group levels mix text and numbers"):
        r_factor_levels(["a", 1])


# Batch-15 review M-5 (TX-AUD-20261011-070): with a second, unobserved state in the alphabet, R computes seqidegrad 0
# and seqipos 1 (the golden), so R stops for a one-state alphabet only; its seqprecarity stops there for both
# with.missing values (the golden), where the definition gives otto * 0 + (1 - otto) * 0 = 0 (A ranks first; no
# transition, complexity 0). The ledger notes say so.
def test_one_state_data_with_a_two_state_alphabet() -> None:
    x = tx.seqdef(
        [row.split("-") for row in _texts("degenerate", "one_state")], alphabet=["A", "B"]
    )
    np.testing.assert_array_equal(
        np.asarray(tx.seqidegrad(x), dtype=float).ravel(),
        _numbers("degenerate_two_alphabet", "seqidegrad"),
    )
    np.testing.assert_array_equal(
        np.asarray(tx.seqipos(x, pos_states=["A"]), dtype=float).ravel(),
        _numbers("degenerate_two_alphabet", "seqipos"),
    )
    for field, include_missing in (("seqprecarity", False), ("seqprecarity_with_missing", True)):
        assert GOLDEN[("degenerate_two_alphabet", field)].startswith("ERROR:")  # R
        np.testing.assert_array_equal(
            np.asarray(tx.seqprecarity(x, include_missing=include_missing), dtype=float).ravel(),
            [0.0, 0.0],
        )


def test_degenerate_ledger_notes_name_the_one_state_alphabet() -> None:
    rows = _ledger_rows()
    for r_name in ("seqidegrad", "seqipos", "seqprecarity"):
        assert "one-state alphabet" in rows[r_name]["notes"], r_name
    assert "a two-state alphabet of which one state is observed" in rows["seqprecarity"]["notes"]


# TX-REP-20261023-007 (batch-15 review M-7): plot.stslist reads sample.meth with match.arg (plot.stslist.R:16), an
# exact name or a unique prefix, case-sensitive (the golden); select_index_rows (and seqIplot through it) took exact
# names only.
@pytest.mark.parametrize("value", ["systematic", "sys", "r"])
def test_select_index_rows_matches_sample_meth_as_match_arg(value: str) -> None:
    from traminex.visualization.series import select_index_rows

    matched = GOLDEN[(f"sample_meth_{value}", "result")]
    x = tx.seqdef(
        [["A", "B", "A"], ["B", "B", "A"], ["A", "A", "B"], ["B", "A", "A"], ["A", "B", "B"]]
    )
    got = select_index_rows(x, sampling=3, sample_meth=value, random_state=5)
    want = select_index_rows(x, sampling=3, sample_meth=matched, random_state=5)
    np.testing.assert_array_equal(got.display_rows, want.display_rows)


@pytest.mark.parametrize("value", ["", "SYS", "x"])
def test_select_index_rows_refuses_what_match_arg_refuses(value: str) -> None:
    from traminex.visualization.series import select_index_rows

    assert GOLDEN[(f"sample_meth_{value or 'empty'}", "result")].startswith("ERROR:")  # R
    x = tx.seqdef([["A", "B"], ["B", "A"]])
    with pytest.raises(ValueError, match="sample_meth must be 'systematic' or 'random'"):
        select_index_rows(x, sampling=1, sample_meth=value)
