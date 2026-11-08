"""Batch-11 event and discrete-rule items against R (``scripts/r_goldens/make_goldens_event_arguments.R``).

TX-AUD-20261019-043: seqerulesdisc builds its rule labels with R's ``paste`` of ``round(HR, 2)`` and ``round(p, 4)``
(TraMineRextras seqerulesdisc.R:89), so a p-value rounding to 0.0002 reads ``2e-04``, as R prints it.
"""

from __future__ import annotations

import contextlib
import csv
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx
from traminex.events import SeqElist

GOLDEN_PATH = Path(__file__).resolve().parent / "r_goldens" / "event_arguments.csv"
with open(GOLDEN_PATH, encoding="utf-8", newline="") as handle:
    GOLDEN = {(row["case"], row["field"]): row["value"] for row in csv.DictReader(handle)}


def _tsef(case: str) -> pd.DataFrame:
    rows = [list(map(int, row.split(","))) for row in GOLDEN[(case, "tsef")].split("|")]
    return pd.DataFrame(rows, columns=["IDPERS", "A", "AST", "B", "BST"])


def _rules(case: str, link: str = "cloglog") -> tuple[list[str], list[str]]:
    """seqerulesdisc's labels on a golden case's data, and the warnings seqerulesdisc emits."""
    tsef = _tsef(case)
    ids = list(tsef.IDPERS)
    agemin = pd.DataFrame({"IDPERS": ids, "AGE": 15})
    agemax = pd.DataFrame({"IDPERS": ids, "AGE": 30})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        data = tx.createdatadiscrete(ids, tsef, ["A", "B"], agemin, agemax)
        events = sorted(
            (person, time, name)
            for person, a, ast, b, bst in tsef.itertuples(index=False, name=None)
            for time, name, seen in ((a, "A", ast), (b, "B", bst))
            if seen == 1
        )
        sequences = tx.seqecreate(
            id=[e[0] for e in events],
            timestamp=[e[1] for e in events],
            event=[e[2] for e in events],
        )
        subsequences = tx.seqefsub(sequences, pmin_support=0.01)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        rules = tx.seqerulesdisc(subsequences, data, tsef, pvalue=1, adjust=False, link=link)
    return list(rules["Hazards"]), [
        str(w.message) for w in caught if issubclass(w.category, UserWarning)
    ]


@pytest.mark.parametrize(
    "case", ["erd_seed11_n80", "erd_seed3_n150", "erd_seed7_n60", "erd_large_hr"]
)
def test_seqerulesdisc_labels_print_r_s_numbers(case: str) -> None:
    assert _rules(case)[0] == GOLDEN[(case, "labels")].split(" ;; ")


@pytest.mark.parametrize("link", ["cloglog", "logit"])
def test_seqerulesdisc_labels_print_a_large_hazard_ratio_as_r_does(link: str) -> None:
    # Batch-11 review M-5: a well-conditioned fit with HR about 1.5e3 (cloglog) and 2.0e3 (logit), both links.
    suffix = "" if link == "cloglog" else "_logit"
    labels, messages = _rules("erd_large_hr", link)
    assert labels == GOLDEN[("erd_large_hr", f"labels{suffix}")].split(" ;; ")
    assert messages == [m for m in GOLDEN[("erd_large_hr", f"warnings{suffix}")].split(" | ") if m]


LABEL = re.compile(r"^(?P<rule>\w+)\(HR=(?P<hr>[^,]+),p=(?P<p>[^)]+)\) => (?P<to>\w+)$")


@pytest.mark.parametrize("link", ["cloglog", "logit"])
def test_seqerulesdisc_labels_on_separated_data_agree_to_the_fit_s_precision(link: str) -> None:
    # On perfectly separated data the fit's last digits depend on floating-point order (R's dqrls and OpenBLAS are not
    # reproduced bit for bit; TX-REP-20261022-004's record): R prints HR about 2.05e10 with 13 significant digits.
    # The coefficient (about 24) agrees to about 5e-11 relative, so HR = exp(coefficient) to about 1.2e-9: the
    # hazard ratio is compared at rtol 1e-8 and the rule, the p-value text and the warnings exactly.
    suffix = "" if link == "cloglog" else "_logit"
    labels, messages = _rules("erd_separated", link)
    (expected,) = GOLDEN[("erd_separated", f"labels{suffix}")].split(" ;; ")
    ours, r = LABEL.match(labels[0]), LABEL.match(expected)
    assert len(labels) == 1 and ours and r
    assert (ours["rule"], ours["to"], ours["p"]) == (r["rule"], r["to"], r["p"])
    assert float(ours["hr"]) == pytest.approx(float(r["hr"]), rel=1e-8)
    assert messages == [m for m in GOLDEN[("erd_separated", f"warnings{suffix}")].split(" | ") if m]


# TX-AUD-20261019-044: seqecontain matches factor(event.list, levels = dict) (seqecontain.R:10), so the requested
# values are read as R's as.character reads them: 1.0 and "1" both name the event 1. Scalars as R passes them.
CONTAIN = {
    "contain_double_1": (1.0, False),
    "contain_integer_1": (1, False),
    "contain_text_1": ("1", False),
    "contain_double_2_5": (2.5, False),
    "contain_double_1_and_2": ([1.0, 2.0], True),
    "contain_double_2_exclude": (2.0, True),
    "contain_double_1_and_unknown_exclude": ([1.0, 3.0], True),
}


@pytest.mark.parametrize("case", sorted(CONTAIN))
def test_seqecontain_reads_event_values_as_r_does(case: str) -> None:
    value, exclude = CONTAIN[case]
    sequences = tx.seqecreate(id=[1, 1, 2], timestamp=[0, 1, 0], event=[1, 2, 2])
    actual = [bool(v) for v in tx.seqecontain(sequences, value, unknown_exclude=exclude)]
    assert actual == [v == "TRUE" for v in GOLDEN[(case, "contain")].split()]


# TX-AUD-20261019-083: seqerules sorts with R's order() (seqerules.R:77-78), so rule text follows R's collation
# (a before B), not code points; numeric columns sort as before. The data of the generator's rules_ cases.
RULE_ROWS = (
    (1, 0, "a"),
    (1, 1, "B"),
    (1, 2, "c"),
    (2, 0, "a"),
    (2, 2, "c"),
    (3, 0, "B"),
    (3, 1, "c"),
    (3, 3, "a"),
    (4, 0, "a"),
    (4, 1, "B"),
    (5, 1, "c"),
    (5, 2, "a"),
    (5, 3, "B"),
    (6, 0, "B"),
    (6, 2, "a"),
    (6, 3, "c"),
)
RULE_CASES = {
    "rules_sort_rule": dict(sortv="rule"),
    "rules_sort_rule_decreasing": dict(sortv="rule", decreasing=True),
    "rules_sort_support_decreasing": dict(sortv="support", decreasing=True),
    "rules_unsorted": dict(),
}


@pytest.mark.parametrize("case", sorted(RULE_CASES))
def test_seqerules_sorts_rule_text_in_r_s_order(case: str) -> None:
    sequences = tx.seqecreate(
        id=[r[0] for r in RULE_ROWS],
        timestamp=[r[1] for r in RULE_ROWS],
        event=[r[2] for r in RULE_ROWS],
    )
    rules = tx.seqerules(tx.seqefsub(sequences, pmin_support=0.3), **RULE_CASES[case])
    assert list(rules["Rules"]) == GOLDEN[(case, "rules")].split(" ;; ")


# TX-AUD-20261010-077: seqeconstraint stops unless count.method is in 1:6 or a known name (seqeconstraint.R:18-22), so
# a fractional number is refused with R's message rather than truncated; TRUE counts as 1. The values as R passes them.
COUNT_METHODS = {
    "cm_double_2_5": 2.5,
    "cm_double_0": 0.0,
    "cm_double_7": 7.0,
    "cm_double_6": 6.0,
    "cm_double_2": 2.0,
    "cm_logical_true": True,
    "cm_name_cwin": "CWIN",
}


@pytest.mark.parametrize("case", sorted(COUNT_METHODS))
def test_seqeconstraint_takes_count_methods_as_r_does(case: str) -> None:
    expected = GOLDEN[(case, "count_method")]
    if expected.startswith(" [!]"):
        with pytest.raises(ValueError, match=re.escape(expected)):
            tx.seqeconstraint(count_method=COUNT_METHODS[case])
    else:
        assert tx.seqeconstraint(count_method=COUNT_METHODS[case]).count_method == float(expected)


def test_seqeconstraint_refuses_a_count_method_number_given_as_text() -> None:
    # R accepts "3" but stores no count method, and seqefsub then fails; traminex refuses it at once.
    assert GOLDEN[("cm_text_3", "count_method")] == "NULL"
    assert GOLDEN[("seqefsub_text_3", "seqefsub_error")] == "argument is of length zero"
    with pytest.raises(ValueError, match=re.escape(" [!] unknown count.method input")):
        tx.seqeconstraint(count_method="3")


def test_seqefsub_names_cwin_when_its_window_is_missing() -> None:
    # R's message names 'CMINWIN' for count method 3, which R's own seqeconstraint maps from "CWIN"; traminex names CWIN.
    r_message = GOLDEN[("seqefsub_cwin_without_window", "seqefsub_error")]
    assert "'CMINWIN'" in r_message
    sequences = tx.seqecreate(id=[1, 1, 2], timestamp=[0, 1, 0], event=[1, 2, 2])
    with pytest.raises(ValueError, match=re.escape(r_message.replace("'CMINWIN'", "'CWIN'"))):
        tx.seqefsub(sequences, pmin_support=0.1, constraint=tx.seqeconstraint(count_method="CWIN"))


# TX-AUD-20261011-057: after the names, seqeapplysub uses a method only when method %in% 1:5 (seqeapplysub.R:28), a
# match by value: 2.5, 0 and 6 leave the constraint's count method, and the text "2" (not "2.0") counts as 2. The c2
# cases count under count_method 2. The values as R passes them.
APPLY_METHODS = {
    "apply_method_null": (None, None),
    "apply_method_double_2_5": (2.5, None),
    "apply_method_double_2": (2.0, None),
    "apply_method_double_0": (0.0, None),
    "apply_method_double_6": (6.0, None),
    "apply_method_logical_true": (True, None),
    "apply_method_text_2": ("2", None),
    "apply_method_text_2_0": ("2.0", None),
    "apply_method_name_count": ("count", None),
    "apply_c2_method_null": (None, 2),
    "apply_c2_method_double_1_5": (1.5, 2),
    "apply_c2_method_double_1": (1.0, 2),
}


@pytest.mark.parametrize("case", sorted(APPLY_METHODS))
def test_seqeapplysub_uses_a_method_as_r_does(case: str) -> None:
    method, count_method = APPLY_METHODS[case]
    sequences = tx.seqecreate(
        id=[1, 1, 1, 2, 2, 3], timestamp=[0, 1, 2, 0, 3, 1], event=list("ABAABB")
    )
    subsequences = tx.seqefsub(sequences, str_subseq=["(A)", "(A)-(B)"])
    constraint = None if count_method is None else tx.seqeconstraint(count_method=count_method)
    counts = tx.seqeapplysub(subsequences, method=method, constraint=constraint).to_numpy()
    expected = [
        [float(cell) for cell in row.split(",")] for row in GOLDEN[(case, "applysub")].split("|")
    ]
    assert counts.tolist() == expected


# TX-AUD-20261011-058: R's "[.seqelist" takes any index vector (seqe.R:22-32); a SeqElist takes lists, integer arrays
# and boolean masks with NumPy's 0-based rules (a traminex contract) and returns a SeqElist with the sequences' weights
# and events. Sequence 1 (0-based) holds only its end event. R's positions are 1-based; these are the same selections.
SUBSETS = {
    "subset_positions_1_3": ([0, 2], np.array([0, 2])),
    "subset_positions_3_1": ([2, 0], np.array([-1, 0])),
    "subset_mask_true_false_true": ([True, False, True], np.array([True, False, True])),
    "subset_positions_none": ([], np.array([], dtype=int)),
}


def _with_end() -> SeqElist:
    sequences = tx.seqecreate(
        id=[1, 1, 2, 3, 3],
        timestamp=[0, 1, 5, 0, 2],
        event=["a", "b", "END", "a", "c"],
        end_event="END",
    )
    sequences.weights = [1, 2.5, 4]
    return sequences


@pytest.mark.parametrize("index_form", [0, 1])
@pytest.mark.parametrize("case", sorted(SUBSETS))
def test_seqelist_takes_index_lists_arrays_and_masks(case: str, index_form: int) -> None:
    source = _with_end()
    subset = source[SUBSETS[case][index_form]]
    assert isinstance(subset, SeqElist)
    assert len(subset) == int(GOLDEN[(case, "length")])
    assert [str(i) for i in tx.seqeid(subset)] == GOLDEN[(case, "ids")].split()
    assert subset.weights.tolist() == [float(w) for w in GOLDEN[(case, "weights")].split()]
    events = [f"{int(float(i))},{t:g},{e}" for i, t, e in zip(subset.id, subset.time, subset.label)]
    assert (" ".join(events) or "none") == GOLDEN[(case, "events")]
    assert all(s is not t for s in subset.sequences for t in source.sequences)


def test_seqelist_indexing_follows_numpy_rules() -> None:
    source = _with_end()
    assert source[-1].as_character() == source.sequences[2].as_character()
    assert source[[-1]].as_character() == [source.sequences[2].as_character()]
    for bad in ([True, False], [0, 3], np.array([0.0, 1.0])):
        with pytest.raises(IndexError):
            source[bad]


# TX-AUD-20261010-076: traminex corrects two R defects of seqefsub and copies a slice's sequences, which R shares; the
# goldens record R's behaviour and the ledger rows (seqefsub, [.seqelist) record the differences. traminex's values
# follow the definitions: a support no subsequence reaches finds none, and weighted=False changes only the counting.
def _weighted_four() -> SeqElist:
    sequences = tx.seqecreate(
        id=[1, 1, 2, 2, 2, 3, 4, 4], timestamp=[0, 3, 1, 2, 5, 4, 0, 2], event=list("ABACBABC")
    )
    sequences.weights = [1, 2.5, 1, 1]
    return sequences


def test_seqefsub_finds_nothing_where_r_stops() -> None:
    assert (
        GOLDEN[("fsub_min_support_100", "seqefsub_error")]
        == "attempt access index 0/0 in VECTOR_ELT"
    )
    result = tx.seqefsub(_weighted_four(), min_support=100)
    assert len(result) == 0 and result.data.shape == (0, 2)


@pytest.mark.parametrize("case", ["fsub_str_unweighted", "fsub_min_support_unweighted"])
def test_unweighted_seqefsub_leaves_the_caller_s_weights(case: str) -> None:
    r_weights = {"fsub_str_unweighted": "1 1 1 1", "fsub_min_support_unweighted": "1 2.5 1 1"}[case]
    assert GOLDEN[(case, "caller_weights")] == r_weights
    sequences = _weighted_four()
    route = (
        dict(str_subseq=["(A)-(B)", "(B)"])
        if case == "fsub_str_unweighted"
        else dict(min_support=1)
    )
    result = tx.seqefsub(sequences, weighted=False, **route)
    assert sequences.weights.tolist() == [1, 2.5, 1, 1]
    if case == "fsub_str_unweighted":
        expected = [[float(v) for v in row.split(",")] for row in GOLDEN[(case, "data")].split("|")]
        assert result.data[["Support", "Count"]].values.tolist() == expected


def test_a_slice_holds_copies_of_its_sequences() -> None:
    assert GOLDEN[("slice_weights_written", "source_weights")] == "9 9 1 1"
    sequences = _weighted_four()
    head = sequences[0:2]
    head.weights = [9, 9]
    assert sequences.weights.tolist() == [1, 2.5, 1, 1]


def test_the_ledger_records_seqefsub_s_corrections_and_the_slice_copy() -> None:
    ledger = Path(__file__).resolve().parents[1] / "docs" / "API_PARITY_DECISIONS.csv"
    with open(ledger, encoding="utf-8", newline="") as handle:
        notes = {row["r_name"]: row["notes"] for row in csv.DictReader(handle)}
    assert "TX-AUD-20261010-076" in notes["seqefsub"]
    assert "TX-AUD-20261010-076" in notes["[.seqelist"]


@pytest.mark.parametrize(
    "events",
    [
        {1.0},
        frozenset({1.0}),
        (value for value in [1.0]),
        {1.0: "x"}.keys(),
        (1.0,),
        np.array([1.0]),
        pd.Series([1.0]),
        np.array(1.0),
    ],
    ids=[
        "set",
        "frozenset",
        "generator",
        "dict_keys",
        "tuple",
        "array",
        "series",
        "zero_dim_array",
    ],
)
def test_seqecontain_reads_any_collection_of_events(events) -> None:
    # Only a single value is wrapped; NumPy counts a set, a generator or a dict view as zero-dimensional, which made
    # them one value that matched nothing (batch-11 review M-2). The list [1.0] gives R's TRUE FALSE.
    sequences = tx.seqecreate(id=[1, 1, 2], timestamp=[0, 1, 0], event=[1, 2, 2])
    expected = [v == "TRUE" for v in GOLDEN[("contain_double_1", "contain")].split()]
    assert [bool(v) for v in tx.seqecontain(sequences, events)] == expected


@pytest.mark.parametrize("case", ["rules_sort_rule", "rules_sort_rule_decreasing"])
def test_seqerules_sorts_rule_text_in_r_s_order_as_pandas_strings(case: str) -> None:
    # With pandas' string inference (pandas 3's default) the rule column is a string dtype, not object; it is
    # still sorted in R's collation (batch-11 review M-1). pandas 3 infers it by default and deprecates the
    # option; pandas 2 infers it only through the option, and pandas 2.1 and 2.2 build such strings with
    # pyarrow only (pandas 2.3 needs none), so the case is skipped only where pandas cannot build them
    # (TX-REP-20261025-003, batch-19 review I-1).
    if int(pd.__version__.split(".")[0]) >= 3:
        inference = contextlib.nullcontext()
    else:
        try:
            with pd.option_context("future.infer_string", True):
                pd.Series(["a"])
        except ImportError as error:
            pytest.skip(
                f"pandas {pd.__version__} infers a string dtype only with pyarrow ({error})"
            )
        inference = pd.option_context("future.infer_string", True)
    with inference:
        sequences = tx.seqecreate(
            id=[r[0] for r in RULE_ROWS],
            timestamp=[r[1] for r in RULE_ROWS],
            event=[r[2] for r in RULE_ROWS],
        )
        rules = tx.seqerules(tx.seqefsub(sequences, pmin_support=0.3), **RULE_CASES[case])
        assert not pd.api.types.is_object_dtype(rules["Rules"])
    assert list(rules["Rules"]) == GOLDEN[(case, "rules")].split(" ;; ")


# Batch-11 review M-4: the labels' numbers are R's paste(round(x, d)); R's round() (nmath/fround.c) picks the nearer of
# floor(x * 10^d) / 10^d and ceiling(x * 10^d) / 10^d, the even one on a tie, which is not Python's round at decimal
# ties (6.165 -> 6.16; 5e-05 -> 0 at 4 digits).
@pytest.mark.parametrize(
    "case,digits",
    [("round_ties_d2", 2), ("round_ties_d4", 4), ("round_edges_d2", 2), ("round_edges_d4", 4)],
)
def test_label_numbers_round_as_r_does(case: str, digits: int) -> None:
    from traminex.discrete_rules import _r_number

    values = [float(text) for text in GOLDEN[(case, "x")].split()]
    expected = GOLDEN[(case, "rounded")].split()
    mismatches = [
        (v, ours, r) for v, r in zip(values, expected) if (ours := _r_number(v, digits)) != r
    ]
    assert not mismatches, f"{len(mismatches)} of {len(values)} differ, e.g. {mismatches[:3]}"
