"""The event and group-comparison surfaces on a WIDE corpus (W-04).

Cycle 6 adjudicated these against 5 sequences / 13 events and two 8x6 matrices,
labelled SAMPLE. That was enough to find three defects; it is not enough to
conclude the rest is clean. This module re-runs the same differentials against
live R on a corpus an order of magnitude larger and deliberately more awkward:

* 38 event sequences of uneven length (1 to 14 events), over 8 event types,
  with non-uniform weights and simultaneous events in 9 of them -- two of the
  40 drawn sequences produce no events and never reach the TSE table, which is
  the shape a real extract produces;
* an unbalanced two-group split (26 / 12) and an unbalanced three-group split
  (12 / 13 / 13), so ``seqecmpgroup`` meets a table that is not square;
* a separate 30 x 10 state-sequence corpus with a ragged right edge, weights
  and a 7 / 23 group split for the group-comparison side.

Generator: ``scripts/r_goldens/make_goldens_events_wide_corpus.R``.
"""

from __future__ import annotations

import csv

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDEN_ATOL, R_GOLDENS_DIR, load_golden_lines

import traminex as tx

# --------------------------------------------------------------------------- #
# Corpora, rebuilt from the fixtures.
# --------------------------------------------------------------------------- #
LEVELS = load_golden_lines("events_wide_levels")
WEIGHTS = np.array([float(v) for v in load_golden_lines("events_wide_weights")])
GROUP2 = load_golden_lines("events_wide_group2")
GROUP3 = load_golden_lines("events_wide_group3")
IDCOST = np.array([float(v) for v in load_golden_lines("events_wide_idcost")])


def _values(name: str) -> dict[str, str]:
    with (R_GOLDENS_DIR / name).open() as handle:
        return {row["key"]: row["value"] for row in csv.DictReader(handle)}


EVENTS = _values("events_wide_values.csv")
GROUP_COMPARE = _values("events_wide_gc_values.csv")


def _eseq():
    frame = pd.read_csv(R_GOLDENS_DIR / "events_wide_tse.csv")
    sequences = tx.seqecreate(frame)
    sequences.weights = WEIGHTS
    return sequences


def _floats(blob: str) -> np.ndarray:
    return np.array([float(token) for token in blob.split(";")], dtype=float)


def _assert_floats(key: str, values, atol: float = R_GOLDEN_ATOL) -> None:
    np.testing.assert_allclose(
        np.asarray(values, dtype=float).ravel(order="F"),
        _floats(EVENTS[key]),
        rtol=0,
        atol=atol,
        equal_nan=True,
    )


def _assert_labels(key: str, labels) -> None:
    assert list(labels) == EVENTS[key].split(";")


# --------------------------------------------------------------------------- #
# 0. The corpus really is wide, and really is awkward.
# --------------------------------------------------------------------------- #
def test_the_corpus_is_as_wide_as_the_module_claims():
    sequences = _eseq()
    frame = pd.read_csv(R_GOLDENS_DIR / "events_wide_tse.csv")

    assert len(sequences) == 38
    assert len(frame) == 234
    assert list(sequences.dictionary) == LEVELS == [f"E{k}" for k in range(1, 9)]

    lengths = [len(item.codes) for item in sequences]
    assert (min(lengths), max(lengths)) == (1, 14)

    simultaneous = sum(1 for item in sequences if len(set(item.times)) < len(item.times))
    assert simultaneous == 9

    assert sorted(pd.Series(GROUP2).value_counts().tolist(), reverse=True) == [26, 12]
    assert sorted(pd.Series(GROUP3).value_counts().tolist(), reverse=True) == [13, 13, 12]


def test_the_rendered_sequences_match_r():
    assert _eseq().as_character() == load_golden_lines("events_wide_aschar")


# --------------------------------------------------------------------------- #
# 1. Mining.
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("min_support", [4, 8, 15])
def test_mining_at_several_support_thresholds_matches_r(min_support):
    mined = tx.seqefsub(_eseq(), min_support=min_support)

    _assert_labels(f"fsub|ms={min_support:g}|subseq", mined.subseq)
    _assert_floats(f"fsub|ms={min_support:g}|support", mined.data["Support"])
    _assert_floats(f"fsub|ms={min_support:g}|count", mined.data["Count"])
    assert len(mined) == int(EVENTS[f"fsub|ms={min_support:g}|n"])


@pytest.mark.parametrize("min_support", [4, 8])
def test_unweighted_mining_matches_r(min_support):
    mined = tx.seqefsub(_eseq(), min_support=min_support, weighted=False)

    _assert_labels(f"fsub_unweighted|ms={min_support:g}|subseq", mined.subseq)
    _assert_floats(f"fsub_unweighted|ms={min_support:g}|count", mined.data["Count"])


@pytest.mark.parametrize("max_k", [1, 2, 3])
def test_depth_limited_mining_matches_r(max_k):
    mined = tx.seqefsub(_eseq(), min_support=6, max_k=max_k)

    _assert_labels(f"fsub|maxk={max_k}|subseq", mined.subseq)
    _assert_floats(f"fsub|maxk={max_k}|count", mined.data["Count"])


@pytest.mark.parametrize("count_method", ["COBJ", "CDIST_O", "CMINWIN", "CDIST"])
def test_every_count_method_matches_r_on_the_wide_corpus(count_method):
    constraint = tx.seqeconstraint(count_method=count_method)
    mined = tx.seqefsub(_eseq(), min_support=6, constraint=constraint)

    _assert_labels(f"fsub|cm={count_method}|subseq", mined.subseq)
    _assert_floats(f"fsub|cm={count_method}|count", mined.data["Count"])


def test_the_windowed_count_method_matches_r():
    constraint = tx.seqeconstraint(count_method="CWIN", window_size=6)
    mined = tx.seqefsub(_eseq(), min_support=6, constraint=constraint)

    _assert_labels("fsub|cm=CWIN|subseq", mined.subseq)
    _assert_floats("fsub|cm=CWIN|count", mined.data["Count"])


@pytest.mark.parametrize(
    "tag,constraint_kwargs",
    [
        ("maxgap3", {"max_gap": 3}),
        ("window8", {"window_size": 8}),
        ("agemax10", {"age_max": 10}),
        ("agemin5", {"age_min": 5}),
        ("gap2_win5", {"max_gap": 2, "window_size": 5}),
    ],
)
def test_constrained_mining_matches_r(tag, constraint_kwargs):
    constraint = tx.seqeconstraint(**constraint_kwargs)
    mined = tx.seqefsub(_eseq(), min_support=5, constraint=constraint)

    _assert_labels(f"fsub|{tag}|subseq", mined.subseq)
    _assert_floats(f"fsub|{tag}|count", mined.data["Count"])


def test_the_count_methods_do_not_all_agree():
    """Otherwise the cases above would prove nothing about the method."""
    counts = {
        method: EVENTS[f"fsub|cm={method}|count"]
        for method in ("COBJ", "CDIST_O", "CMINWIN", "CDIST")
    }
    assert len(set(counts.values())) > 1


# --------------------------------------------------------------------------- #
# 2. The count matrix.
# --------------------------------------------------------------------------- #
def _base_mined():
    return tx.seqefsub(_eseq(), min_support=8)


def test_the_base_mining_matches_r():
    _assert_labels("base|subseq", _base_mined().subseq)


@pytest.mark.parametrize(
    "method", ["count", "presence", "age", "COBJ", "CDIST_O", "CMINWIN", "CDIST"]
)
def test_the_count_matrix_matches_r_under_every_method(method):
    matrix = tx.seqeapplysub(_base_mined(), method=method)

    rows, columns = (int(v) for v in EVENTS[f"applysub|{method}|dim"].split("x"))
    assert matrix.shape == (rows, columns)
    _assert_floats(f"applysub|{method}", matrix.to_numpy())


def test_the_count_matrix_matches_r_under_a_windowed_constraint():
    matrix = tx.seqeapplysub(
        _base_mined(), constraint=tx.seqeconstraint(count_method="CWIN", window_size=6)
    )

    _assert_floats("applysub|CWIN", matrix.to_numpy())


# --------------------------------------------------------------------------- #
# 3. Association rules.
# --------------------------------------------------------------------------- #
def test_the_rule_table_matches_r_column_by_column():
    rules = tx.seqerules(_base_mined())

    rows, columns = (int(v) for v in EVENTS["rules|dim"].split("x"))
    assert rules.shape == (rows, columns)
    assert list(rules["Rules"]) == EVENTS["rules|Rules"].split(";")
    for column in (
        "Support",
        "Conf",
        "Lift",
        "Standardlift",
        "JMeasure",
        "ImplicStat",
        "p.value",
        "p.valueB1",
        "p.valueB2",
    ):
        _assert_floats(f"rules|{column}", rules[column])


@pytest.mark.parametrize("sortv", ["support", "lift", "jmeasure"])
def test_rule_sorting_matches_r(sortv):
    rules = tx.seqerules(_base_mined(), sortv=sortv, decreasing=True)

    assert list(rules["Rules"]) == EVENTS[f"rules_sorted|{sortv}"].split(";")


# --------------------------------------------------------------------------- #
# 4. Discriminating subsequences, including a three-group table.
# --------------------------------------------------------------------------- #
CMPGROUP_CASES = [
    ("g2_w", "GROUP2", True, "chisq"),
    ("g2_uw", "GROUP2", False, "chisq"),
    ("g3_w", "GROUP3", True, "chisq"),
    ("g2_bonf", "GROUP2", True, "bonferroni"),
]


def _r_cmpgroup(tag: str) -> pd.DataFrame:
    labels = EVENTS[f"cmpgroup|{tag}|subseq"].split(";")
    columns = EVENTS[f"cmpgroup|{tag}|columns"].split(";")
    values = _floats(EVENTS[f"cmpgroup|{tag}|values"]).reshape(
        (len(labels), len(columns)), order="F"
    )
    return pd.DataFrame(values, columns=columns, index=labels)


@pytest.mark.parametrize("tag,group_name,weighted,method", CMPGROUP_CASES)
def test_group_comparison_of_subsequences_matches_r(tag, group_name, weighted, method):
    """Compared row for row by ``index``, which identifies the subsequence.

    The row ORDER is asserted separately below (the fixture keeps 15
    significant digits, so values are compared to 1e-10 here).
    """
    group = {"GROUP2": GROUP2, "GROUP3": GROUP3}[group_name]
    expected = _r_cmpgroup(tag)

    compared = tx.seqecmpgroup(_base_mined(), group=group, weighted=weighted, method=method)

    assert sorted(compared.subseq) == sorted(expected.index)
    assert list(compared.data.columns) == list(expected.columns)
    np.testing.assert_allclose(
        compared.data.sort_values("index").to_numpy(dtype=float),
        expected.sort_values("index").to_numpy(),
        rtol=0,
        atol=1e-10,
    )


@pytest.mark.parametrize("tag,group_name,weighted,method", CMPGROUP_CASES)
def test_row_order_matches_r(tag, group_name, weighted, method):
    """Both sides order by the chi-square statistic, descending, ties stable.

    Where the statistic is zero in exact arithmetic, R's Yates-corrected
    arithmetic leaves residue of order 1e-29 and sorts on it. traminex used to
    cancel it to exact 0.0 and fall back to enumeration order there; since
    audit S5-16 it computes the statistic as R's chisq.test does, residue
    included, so the whole order is R's -- the tied blocks too.
    """
    group = {"GROUP2": GROUP2, "GROUP3": GROUP3}[group_name]
    expected = _r_cmpgroup(tag)

    compared = tx.seqecmpgroup(_base_mined(), group=group, weighted=weighted, method=method)

    assert list(compared.subseq) == list(expected.index)


def test_only_the_two_group_cases_carry_numerically_zero_statistics():
    """Pins that the order test above covers tied zero-statistic blocks."""
    zeros = {
        tag: int((_r_cmpgroup(tag)["statistic"].to_numpy() < 1e-12).sum())
        for tag, *_ in CMPGROUP_CASES
    }
    assert zeros == {"g2_w": 11, "g2_uw": 22, "g3_w": 0, "g2_bonf": 11}


def test_the_bonferroni_correction_changes_the_p_values():
    plain = _r_cmpgroup("g2_w").sort_values("index")["p.value"].to_numpy()
    corrected = _r_cmpgroup("g2_bonf").sort_values("index")["p.value"].to_numpy()
    assert not np.allclose(plain, corrected, rtol=0, atol=1e-9)


# --------------------------------------------------------------------------- #
# 5. seqentrans, seqecontain and the accessors.
# --------------------------------------------------------------------------- #
def test_transition_and_event_counts_match_r():
    augmented = tx.seqentrans(_base_mined())

    _assert_floats("entrans|ntrans", augmented.data["ntrans"])
    _assert_floats("entrans|nevent", augmented.data["nevent"])


@pytest.mark.parametrize(
    "key,unknown_exclude", [("contain|any_E1E2", False), ("contain|all_E1E2", True)]
)
def test_event_containment_matches_r(key, unknown_exclude):
    contained = tx.seqecontain(_eseq(), ["E1", "E2"], unknown_exclude=unknown_exclude)

    _assert_floats(key, contained.astype(float))


def test_the_accessors_match_r():
    sequences = _eseq()
    _assert_floats("accessor|seqeid", tx.seqeid(sequences))
    _assert_floats("accessor|seqeweight", tx.seqeweight(sequences))
    _assert_floats("accessor|seqelength", tx.seqelength(sequences))


# --------------------------------------------------------------------------- #
# 6. Event-sequence edit distance across the whole option grid.
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("interval", ["No", "previous", "next"])
@pytest.mark.parametrize("norm", ["none", "max", "YujianBo"])
def test_event_edit_distance_matches_r_on_the_wide_corpus(interval, norm):
    distance = tx.seqedist(_eseq(), IDCOST, 0.5, interval=interval, norm=norm)

    _assert_floats(f"edist|{interval}|{norm}", distance)


# --------------------------------------------------------------------------- #
# 7. The wider group-comparison corpus, with missing data.
# --------------------------------------------------------------------------- #
GC_ALPHABET = load_golden_lines("events_wide_gc_alphabet")
GC_WEIGHTS = np.array([float(v) for v in load_golden_lines("events_wide_gc_weights")])
GC_GROUP = load_golden_lines("events_wide_gc_group")
GC_SET = load_golden_lines("events_wide_gc_set")


def _gc_diss() -> np.ndarray:
    return pd.read_csv(R_GOLDENS_DIR / "events_wide_gc_diss.csv").to_numpy(dtype=float)


def _gc_sequences():
    frame = pd.read_csv(R_GOLDENS_DIR / "events_wide_gc_states.csv", dtype=str)
    frame = frame.replace("<NA>", None)
    return tx.seqdef(frame, alphabet=list(GC_ALPHABET), weights=GC_WEIGHTS)


def _gc_floats(key: str) -> np.ndarray:
    return np.array([float(token) for token in GROUP_COMPARE[key].split(";")], dtype=float)


def test_the_group_comparison_corpus_is_wide_and_ragged():
    sequences = _gc_sequences()
    assert (sequences.nseq, sequences.seqlen) == (30, 10)
    assert sorted(pd.Series(GC_GROUP).value_counts().tolist(), reverse=True) == [23, 7]
    # a ragged right edge really is present
    assert sequences.void_code is not None
    assert bool(np.any(sequences.data == sequences.void_code))


@pytest.mark.parametrize("stat", ["all", "LRT", "BIC"])
@pytest.mark.parametrize("weighted_name", ["FALSE", "TRUE", "by.group"])
def test_diss_compare_matches_r_on_the_wider_corpus(stat, weighted_name):
    weighted = {"FALSE": False, "TRUE": True, "by.group": "by.group"}[weighted_name]

    result = tx.dissCompare(
        _gc_diss(),
        GC_GROUP,
        s=0,
        stat=stat,
        weighted=weighted,
        weights=GC_WEIGHTS,
        inc_df=True,
    )

    key = f"diss|stat={stat}|weighted={weighted_name}"
    if stat == "all":
        assert list(result.columns) == GROUP_COMPARE[f"{key}|columns"].split(";")
    else:
        # R drops to a bare vector for a single stat -- no dim, no colnames --
        # while traminex keeps a labelled frame throughout (the "result
        # orientation" re-expression in docs/API_CONSISTENCY.md).
        assert GROUP_COMPARE[f"{key}|columns"] == ""
    np.testing.assert_allclose(
        result.to_numpy(dtype=float).ravel(order="F"),
        _gc_floats(f"{key}|values"),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


@pytest.mark.parametrize("level", ["first", "second"])
def test_each_set_matches_r_run_on_that_set_alone(level):
    """Derived without R's set= route, whose index defect W-02 recorded."""
    rows = np.flatnonzero(np.asarray(GC_SET) == level)
    matrix = _gc_diss()

    result = tx.dissCompare(
        matrix[np.ix_(rows, rows)],
        [GC_GROUP[i] for i in rows],
        s=0,
        stat="all",
        weights=GC_WEIGHTS[rows],
        inc_df=True,
    )

    np.testing.assert_allclose(
        result.to_numpy(dtype=float).ravel(order="F"),
        _gc_floats(f"diss_set_alone|{level}|values"),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ seqdist.\ 'with\.missing'\ set\ as\ FALSE\ as\ 'seqdata'\ has\ no\ non\-void\ missing\ values:UserWarning"
)
@pytest.mark.parametrize("method", ["LCS", "OM"])
@pytest.mark.parametrize("weighted_name", ["FALSE", "TRUE"])
def test_seq_compare_matches_r_on_the_ragged_corpus(method, weighted_name):
    weighted = weighted_name == "TRUE"
    kwargs: dict[str, object] = {"with_missing": True}
    if method == "OM":
        kwargs.update(subm="CONSTANT", indel=1.0)

    result = tx.seqCompare(
        _gc_sequences(),
        group=GC_GROUP,
        s=0,
        stat="all",
        method=method,
        weighted=weighted,
        **kwargs,
    )

    np.testing.assert_allclose(
        result.to_numpy(dtype=float).ravel(order="F"),
        _gc_floats(f"seqCompare|{method}|weighted={weighted_name}"),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


def test_how_many_quantities_this_module_compares():
    """The denominator, so the report's number is checkable."""
    event_values = sum(len(blob.split(";")) for blob in EVENTS.values())
    group_values = sum(len(blob.split(";")) for blob in GROUP_COMPARE.values())

    assert len(EVENTS) == 100
    assert len(GROUP_COMPARE) == 24
    assert event_values == 56_939
    assert group_values == 80
    assert event_values + group_values == 57_019
