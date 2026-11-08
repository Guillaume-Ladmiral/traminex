"""Deep and degenerate inputs to the ported event-mining engine (cycle-7 W-01).

``traminex/events/_engine.py`` is a hand-written port of TraMineR's vendored C++
miner. Ordinary input was golden-pinned long ago; these are the two shapes it
had never been fed:

**Deep.** Long runs of *simultaneous* events (zero inter-event gaps). That run
length is the recursion depth of ``SequenceEventNode::getNextWithGap``
(``SOURCES/TraMineR/src/eventseq.h:100``) -- the one genuinely self-recursive
routine the port reproduces -- reached only from ``tmrseqedist`` with
``interval="next"`` (``seqedist.cpp:171,186``). C recurses on the machine stack
and answers; the port must answer too.

**Degenerate.** One event, everything simultaneous, unsorted timestamps,
duplicate events, negative timestamps, boundary constraints, and support and
``max_k`` limits that keep nothing.

Where R itself aborts, traminex deliberately does *not* copy the abort; each
such divergence is asserted against the recorded R error text so it stays a
decision rather than an accident.

Generator: ``scripts/r_goldens/make_goldens_event_engine_edges.R``.
"""

from __future__ import annotations

import csv
import warnings

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDEN_ATOL, R_GOLDENS_DIR, load_golden_lines

import traminex as tx

# --------------------------------------------------------------------------- #
# Fixtures loaded from the R goldens.
# --------------------------------------------------------------------------- #
RUN_LEN = int(load_golden_lines("event_edges_deep_runlen")[0])
DEEP_LEVELS = load_golden_lines("event_edges_deep_levels")
DEEP_IDCOST = [float(v) for v in load_golden_lines("event_edges_deep_idcost")]


def _r_floats(blob: str) -> np.ndarray:
    """Parse R ``format(digits = 17)`` output; it right-pads with spaces."""
    return np.array([float(token.strip()) for token in blob.split(";")], dtype=float)


def _deep_seqedist_goldens() -> dict[tuple[str, str], np.ndarray]:
    path = R_GOLDENS_DIR / "event_edges_deep_seqedist.csv"
    out: dict[tuple[str, str], np.ndarray] = {}
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        interval, norm, blob = line.split(",", 2)
        out[(interval, norm)] = _r_floats(blob)
    return out


DEEP_SEQEDIST = _deep_seqedist_goldens()


def _degenerate_goldens() -> dict[str, str]:
    path = R_GOLDENS_DIR / "event_edges_degenerate.csv"
    with path.open() as handle:
        return {row["key"]: row["value"] for row in csv.DictReader(handle)}


DEGENERATE = _degenerate_goldens()


# --------------------------------------------------------------------------- #
# Corpora -- byte-for-byte the ones the generator builds.
# --------------------------------------------------------------------------- #
def _deep_eseq():
    frames = [
        pd.DataFrame(
            {
                "id": "a",
                "timestamp": 1.0,
                "event": ["E%d" % (k % 3) for k in range(1, RUN_LEN + 1)],
            }
        ),
        pd.DataFrame({"id": ["a"], "timestamp": [5.0], "event": ["F"]}),
        pd.DataFrame(
            {
                "id": "b",
                "timestamp": 1.0,
                "event": ["E%d" % (k % 2) for k in range(1, RUN_LEN + 1)],
            }
        ),
        pd.DataFrame({"id": ["b"], "timestamp": [7.0], "event": ["F"]}),
        pd.DataFrame({"id": "c", "timestamp": [1.0, 2.0, 3.0], "event": ["E0", "E1", "F"]}),
    ]
    return tx.seqecreate(pd.concat(frames, ignore_index=True))


def _mk(ids, times, events):
    return tx.seqecreate(pd.DataFrame({"id": ids, "timestamp": times, "event": events}))


def _base():
    return _mk(
        ["a", "a", "a", "b", "b", "c"],
        [1, 2, 3, 1, 2, 1],
        ["A", "B", "C", "A", "C", "B"],
    )


# --------------------------------------------------------------------------- #
# 1. Deep simultaneous runs.
# --------------------------------------------------------------------------- #
def test_deep_corpus_matches_the_r_event_alphabet():
    eseq = _deep_eseq()
    assert list(eseq.dictionary) == DEEP_LEVELS
    assert RUN_LEN > 1000, "the run must exceed CPython's default recursion limit"


@pytest.mark.parametrize("interval", ["No", "previous", "next"])
@pytest.mark.parametrize("norm", ["none", "max", "YujianBo"])
def test_seqedist_over_a_1200_event_simultaneous_run_matches_r(interval, norm):
    """``interval="next"`` is the only caller of the recursive getNextWithGap.

    Before the de-recursion this raised ``RecursionError`` for that interval at
    this run length while R returned a value.
    """
    eseq = _deep_eseq()

    got = tx.seqedist(eseq, DEEP_IDCOST, 0.25, interval=interval, norm=norm)

    expected = DEEP_SEQEDIST[(interval, norm)].reshape(got.shape, order="F")
    np.testing.assert_allclose(got, expected, rtol=0, atol=R_GOLDEN_ATOL)


def test_seqedist_next_differs_from_the_default_interval_on_the_deep_corpus():
    """Guards the test above against passing for the wrong reason."""
    assert not np.allclose(
        DEEP_SEQEDIST[("next", "none")], DEEP_SEQEDIST[("No", "none")], rtol=0, atol=1e-9
    )


def test_seqedist_next_is_independent_of_ambient_python_stack_depth():
    """The threshold used to move with how deep the caller already was."""
    eseq = _deep_eseq()
    expected = DEEP_SEQEDIST[("next", "none")]

    def nest(depth):
        if depth:
            return nest(depth - 1)
        return tx.seqedist(eseq, DEEP_IDCOST, 0.25, interval="next", norm="none")

    got = nest(300)

    np.testing.assert_allclose(
        got, expected.reshape(got.shape, order="F"), rtol=0, atol=R_GOLDEN_ATOL
    )


# --------------------------------------------------------------------------- #
# 2. Degenerate corpora that R can also run.
# --------------------------------------------------------------------------- #
def _assert_subseq(key, subseq):
    assert list(subseq) == DEGENERATE[key].split(";")


def _assert_floats(key, values):
    np.testing.assert_allclose(
        np.asarray(values, dtype=float),
        _r_floats(DEGENERATE[key]),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


def test_base_corpus_rendering_and_mining_match_r():
    eseq = _base()
    assert eseq.as_character() == DEGENERATE["base|as.character"].split(";")
    assert list(eseq.dictionary) == DEGENERATE["base|levels"].split(";")
    mined = tx.seqefsub(eseq, min_support=1)
    _assert_subseq("base|fsub_subseq", mined.subseq)
    _assert_floats("base|fsub_count", mined.data["Count"])
    _assert_floats("base|fsub_support", mined.data["Support"])


def test_single_event_sequence_matches_r():
    eseq = _mk(["a"], [1], ["A"])
    assert eseq.as_character() == DEGENERATE["one_event|as.character"].split(";")
    mined = tx.seqefsub(eseq, min_support=1)
    _assert_subseq("one_event|fsub_subseq", mined.subseq)
    _assert_floats("one_event|fsub_count", mined.data["Count"])
    _assert_floats("one_event|seqedist", tx.seqedist(eseq, [1.0], 1.0).ravel(order="F"))


def test_all_events_at_one_timestamp_match_r():
    eseq = _mk(["a"] * 3 + ["b"] * 3, [1] * 6, ["A", "B", "C"] * 2)
    assert eseq.as_character() == DEGENERATE["same_ts|as.character"].split(";")
    mined = tx.seqefsub(eseq, min_support=1)
    _assert_subseq("same_ts|fsub_subseq", mined.subseq)
    _assert_floats("same_ts|fsub_count", mined.data["Count"])
    _assert_floats("same_ts|fsub_support", mined.data["Support"])
    _assert_floats(
        "same_ts|seqedist_next",
        tx.seqedist(eseq, [1.0] * 3, 1.0, interval="next").ravel(order="F"),
    )


def test_unsorted_timestamps_are_sorted_like_r():
    eseq = _mk(["a"] * 3, [3, 1, 2], ["C", "A", "B"])
    assert eseq.as_character() == DEGENERATE["unsorted|as.character"].split(";")
    _assert_subseq("unsorted|fsub_subseq", tx.seqefsub(eseq, min_support=1).subseq)


def test_duplicate_events_at_one_timestamp_match_r():
    eseq = _mk(["a"] * 3, [1, 1, 2], ["A", "A", "B"])
    assert eseq.as_character() == DEGENERATE["dup_events|as.character"].split(";")
    mined = tx.seqefsub(eseq, min_support=1)
    _assert_subseq("dup_events|fsub_subseq", mined.subseq)
    _assert_floats("dup_events|fsub_count", mined.data["Count"])


def test_negative_timestamps_match_r():
    eseq = _mk(["a"] * 3, [-2, -1, 3], ["A", "B", "C"])
    assert eseq.as_character() == DEGENERATE["neg_ts|as.character"].split(";")
    mined = tx.seqefsub(eseq, min_support=1)
    _assert_subseq("neg_ts|fsub_subseq", mined.subseq)
    _assert_floats("neg_ts|fsub_count", mined.data["Count"])
    _assert_floats("neg_ts|seqedist", tx.seqedist(eseq, [1.0] * 3, 1.0).ravel(order="F"))


@pytest.mark.parametrize(
    "key,constraint_kwargs",
    [
        ("cst|maxgap_zero", {"max_gap": 0}),
        ("cst|windowsize_zero", {"window_size": 0}),
        ("cst|cwin_windowsize_zero", {"count_method": "CWIN", "window_size": 0}),
    ],
)
def test_zero_width_constraints_still_mine_single_events_like_r(key, constraint_kwargs):
    mined = tx.seqefsub(_base(), min_support=1, constraint=tx.seqeconstraint(**constraint_kwargs))
    _assert_subseq(key, mined.subseq)


def test_cwin_with_zero_window_counts_like_r():
    mined = tx.seqefsub(
        _base(),
        min_support=1,
        constraint=tx.seqeconstraint(count_method="CWIN", window_size=0),
    )
    _assert_floats("cst|cwin_windowsize_zero_count", mined.data["Count"])


@pytest.mark.parametrize(
    "key,kwargs",
    [
        ("ms|zero", {"min_support": 0}),
        ("ms|negative", {"min_support": -1}),
        ("pms|zero", {"pmin_support": 0}),
        ("maxk|one", {"min_support": 1, "max_k": 1}),
    ],
)
def test_support_and_depth_boundaries_match_r(key, kwargs):
    _assert_subseq(key, tx.seqefsub(_base(), **kwargs).subseq)


@pytest.mark.parametrize(
    "key,idcost,vparam",
    [
        ("edist|neg_idcost", [-1.0, -1.0, -1.0], 1.0),
        ("edist|neg_vparam", [1.0, 1.0, 1.0], -1.0),
    ],
)
def test_degenerate_seqedist_arguments_match_r(key, idcost, vparam):
    got = tx.seqedist(_base(), idcost, vparam)
    expected = _r_floats(DEGENERATE[key]).reshape(got.shape, order="F")
    np.testing.assert_allclose(got, expected, rtol=0, atol=R_GOLDEN_ATOL, equal_nan=True)


def test_zero_indel_cost_gives_zero_distances_where_r_gives_nan():
    """With every indel cost 0, every raw distance and Yujian-Bo's bound are 0.

    R divides 0 by 0 and returns NaN off the diagonal; TraMineR's own normaliser
    maps a raw distance of 0 to 0 (distancecalculator.h), so traminex returns 0
    (R defect corrected, TX-AUD-20261010-037; ledger seqedist row). The golden
    keeps R's NaN so the divergence stays visible. No runtime warning either.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        got = tx.seqedist(_base(), [0.0, 0.0, 0.0], 1.0)

    np.testing.assert_array_equal(got, np.zeros_like(got))
    r_values = _r_floats(DEGENERATE["edist|zero_idcost"]).reshape(got.shape, order="F")
    off_diagonal = ~np.eye(got.shape[0], dtype=bool)
    assert np.all(np.isnan(r_values[off_diagonal]))


# --------------------------------------------------------------------------- #
# 3. Refusals -- traminex must name the argument, and must never copy R's abort.
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "key,kwargs",
    [
        ("cst|maxgap_gt_window", {"max_gap": 5, "window_size": 2}),
        ("cst|agemin_gt_agemax", {"age_min": 5, "age_max": 2}),
    ],
)
def test_contradictory_constraints_refuse_with_r_wording(key, kwargs):
    with pytest.raises(ValueError) as excinfo:
        tx.seqeconstraint(**kwargs)
    assert str(excinfo.value).strip() == DEGENERATE[key].removeprefix("ERROR:").strip()


def test_cwin_without_a_window_size_refuses_naming_the_right_method():
    """R's own message misnames method 3 as 'CMINWIN'; traminex names CWIN.

    Deliberate divergence: the wording is a diagnostic, not a computed value,
    and R's is factually wrong (method 3 is CWIN, method 4 is CMINWIN).
    """
    assert "CMINWIN" in DEGENERATE["cst|cwin_no_windowsize"]

    with pytest.raises(ValueError) as excinfo:
        tx.seqefsub(_base(), min_support=1, constraint=tx.seqeconstraint(count_method="CWIN"))
    assert "CWIN" in str(excinfo.value)
    assert "CMINWIN" not in str(excinfo.value)


def test_missing_minimum_support_refuses_like_r():
    with pytest.raises(ValueError) as excinfo:
        tx.seqefsub(_base())
    assert "minimum support" in str(excinfo.value)
    assert "ERROR" in DEGENERATE["ms|none"]


def test_wrong_indel_cost_length_refuses_like_r():
    with pytest.raises(ValueError) as excinfo:
        tx.seqedist(_base(), [1.0, 1.0], 1.0)
    assert "one indel cost per event" in str(excinfo.value)


# --------------------------------------------------------------------------- #
# 4. Deliberate divergences: R aborts in C, traminex returns an empty result.
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "key,kwargs",
    [
        ("ms|huge", {"min_support": 1000}),
        ("ms|equal_n", {"min_support": 3}),
        ("pms|one", {"pmin_support": 1}),
        ("pms|gt_one", {"pmin_support": 1.5}),
        ("maxk|zero", {"min_support": 1, "max_k": 0}),
    ],
)
def test_mining_that_keeps_nothing_returns_an_empty_list_where_r_aborts(key, kwargs):
    """R fails in C with "attempt access index 0/0 in VECTOR_ELT".

    Copying that would turn "no subsequence met the threshold" -- a legitimate
    analytical outcome -- into a crash. traminex returns the empty list.
    """
    assert DEGENERATE[key].startswith("ERROR: attempt access index 0/0")

    mined = tx.seqefsub(_base(), **kwargs)

    assert list(mined.subseq) == []
    assert len(mined) == 0
    assert mined.data.shape[0] == 0


@pytest.mark.parametrize(
    "key,constraint_kwargs",
    [
        ("cst|agemax_zero", {"age_max": 0}),
        ("cst|excludes_everything", {"age_min": 100, "age_max": 200}),
    ],
)
def test_constraints_that_exclude_everything_return_empty_where_r_aborts(key, constraint_kwargs):
    assert DEGENERATE[key].startswith("ERROR: attempt access index 0/0")

    mined = tx.seqefsub(_base(), min_support=1, constraint=tx.seqeconstraint(**constraint_kwargs))

    assert list(mined.subseq) == []


def test_seqeapplysub_on_an_empty_mining_returns_a_zero_column_matrix():
    matrix = tx.seqeapplysub(tx.seqefsub(_base(), min_support=1000))
    assert matrix.shape == (3, 0)


def test_seqerules_without_a_multi_transition_subsequence_returns_no_rules():
    rules = tx.seqerules(tx.seqefsub(_base(), min_support=1, max_k=1))
    assert rules.shape[0] == 0
    assert list(rules.columns)[:2] == ["Rules", "Support"]


def test_empty_event_corpus_is_handled_throughout():
    empty = tx.seqecreate(pd.DataFrame({"id": [], "timestamp": [], "event": []}))
    assert len(empty) == 0
    assert empty.as_character() == []
    assert list(tx.seqefsub(empty, min_support=1).subseq) == []
    assert tx.seqedist(empty, [], 1.0).shape == (0, 0)


# --------------------------------------------------------------------------- #
# 5. seqecmpgroup degeneracies -- R aborts; traminex must refuse, not empty out.
# --------------------------------------------------------------------------- #
def test_seqecmpgroup_refuses_a_single_group_instead_of_silently_dropping_every_row():
    """A one-level grouping gives df = 0, so every p-value is NaN.

    ``NaN <= pvalue_limit`` is False, so the result silently came back empty --
    a wrong count from a public function. R aborts with "incorrect number of
    dimensions"; traminex names the argument instead.
    """
    assert DEGENERATE["cmpgroup|single_group"].startswith("ERROR:")
    mined = tx.seqefsub(_base(), min_support=1)

    with pytest.raises(ValueError) as excinfo:
        tx.seqecmpgroup(mined, group=["g", "g", "g"])

    message = str(excinfo.value)
    assert "group" in message
    assert "two" in message or "2" in message


def test_seqecmpgroup_refuses_a_group_of_the_wrong_length():
    """R: "variable lengths differ (found for 'group')"; a numpy broadcast
    error from inside a private helper is not an acceptable substitute."""
    assert "variable lengths differ" in DEGENERATE["cmpgroup|wrong_length_group"]
    mined = tx.seqefsub(_base(), min_support=1)

    with pytest.raises(ValueError) as excinfo:
        tx.seqecmpgroup(mined, group=["g", "h"])

    message = str(excinfo.value)
    assert "group" in message
    assert "3" in message and "2" in message


def test_seqecmpgroup_refuses_a_group_level_carrying_no_weight():
    """A zero-weight level makes the expected counts 0, so the statistic is NaN
    and every row was silently dropped. R aborts with "NA/NaN argument"."""
    assert DEGENERATE["cmpgroup|zero_weight_level"].startswith("ERROR:")
    eseq = _base()
    eseq.weights = [0.0, 1.0, 1.0]
    mined = tx.seqefsub(eseq, min_support=1)

    with pytest.raises(ValueError) as excinfo:
        tx.seqecmpgroup(mined, group=["g", "h", "h"])

    assert "weight" in str(excinfo.value).lower()


def test_seqecmpgroup_with_two_usable_groups_still_matches_the_r_shape():
    """The refusals above must not narrow the ordinary path."""
    mined = tx.seqefsub(_base(), min_support=1)

    result = tx.seqecmpgroup(mined, group=["g", "g", "h"])

    expected_shape = tuple(int(v) for v in DEGENERATE["cmpgroup|two_groups_shape"].split(";"))
    assert result.data.shape == expected_shape
