"""The findings of the 2026-10-25 study, fixed on 2026-11-07, pinned against what R does.

The supervision repository's study of 2026-10-25 (its ``docs/history/2026-10-25/``) found:

- F-01: ``seqdef(states=...)`` renamed one state at a time, as R does (seqdef.R:78-81), so a cell
  renamed early was renamed again and distinct states merged, against seqdef.Rd:42. traminex now
  maps all states at once; a ``DIVERGENCE`` row of docs/API_PARITY_DECISIONS.csv records it.
- F-04: ``seqdist(method="NMS")`` raised ``IndexError`` when every pair of a slice had an empty
  side; R returns a matrix.
- F-05: ``seqrep`` and ``dissrep`` with weights of zero total raised ``ZeroDivisionError``; R stops
  too, so traminex now names the cause.
- One sequence: R's ``plot.stslist`` takes any string other than a method name as a sort vector of
  length 1 (plot.stslist.R:62-76). A Python string is a method name, so traminex refuses it and
  takes a one-element list as the sort vector (a ``DIVERGENCE`` row of the ledger).
- O-7: ``STSList.copy()`` was said to share ``issues`` with the original; it does not.

R's values and behaviour come from ``scripts/r_goldens/make_goldens_study_findings.R``.
"""

from __future__ import annotations

import csv

import numpy as np
import pytest
from r_goldens_utils import R_GOLDEN_ATOL, R_GOLDENS_DIR

import traminex as tx
from traminex import dissrep, seqdef, seqdist
from traminex.cluster import seqrep
from traminex.visualization.series import resolve_sort_order


def _read(name: str) -> dict[tuple[str, str], str]:
    with (R_GOLDENS_DIR / name).open(encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(line for line in handle if not line.startswith("#")))
    return {(case, key): value for case, key, value in rows[1:]}


ITEMS = _read("study_findings.csv")


def _matrix(case: str, n: int) -> np.ndarray:
    return np.array(ITEMS[(case, "distances")].split(" "), dtype=float).reshape(n, n, order="F")


def _cells(seq) -> list[list[str]]:
    return [[seq.alphabet[code] for code in row] for row in np.asarray(seq.data)]


THREE = [[None], ["B"], [None]]
FOUR = [[None, None], ["B", None], [None, None], ["B", "C"]]


@pytest.mark.parametrize("prox", [None, np.array([[1.0]])])
def test_nms_with_only_empty_pairs_in_a_slice_matches_r(prox):
    seq = seqdef(THREE, right="DEL")
    case = "nms_empty_three" if prox is None else "nms_soft_empty_three"
    kwargs = {} if prox is None else {"prox": prox}
    # R warns too: " [!!] Empty sequence(s)  1, 3  may cause inconsistent results!"
    with pytest.warns(UserWarning, match=r"Empty sequence\(s\)  1, 3"):
        found = seqdist(seq, method="NMS", **kwargs)
    np.testing.assert_allclose(found, _matrix(case, 3), rtol=0, atol=R_GOLDEN_ATOL)


@pytest.mark.parametrize("method", ["NMS", "NMSMST", "SVRspell"])
def test_subsequence_distances_with_empty_sequences_match_r(method):
    seq = seqdef(FOUR, right="DEL")
    with pytest.warns(UserWarning, match=r"Empty sequence\(s\)  1, 3"):
        found = seqdist(seq, method=method)
    np.testing.assert_allclose(
        found,
        _matrix(f"{method.lower()}_empty_four", 4),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


def _r_states(case: str) -> tuple[list[str], list[list[str]]]:
    alphabet = ITEMS[(case, "alphabet")].split(" ")
    rows = [ITEMS[(case, f"row{i}")].split(" ") for i in (1, 2)]
    return alphabet, rows


def test_states_with_names_disjoint_from_the_alphabet_match_r():
    seq = seqdef([["1", "2", "3"], ["3", "2", "1"]], states=["x", "y", "z"])
    alphabet, rows = _r_states("states_disjoint")
    assert seq.alphabet == alphabet
    assert _cells(seq) == rows


@pytest.mark.parametrize(
    ("case", "data", "kwargs", "expected"),
    [
        (
            "states_shift",
            [["1", "2", "3"], ["3", "2", "1"]],
            {"states": ["2", "3", "4"]},
            [["2", "3", "4"], ["4", "3", "2"]],
        ),
        (
            "states_swap",
            [["A", "B"], ["B", "A"]],
            {"alphabet": ["A", "B"], "states": ["B", "A"]},
            [["B", "A"], ["A", "B"]],
        ),
    ],
)
def test_states_rename_every_state_once_where_r_merges_them(case, data, kwargs, expected):
    seq = seqdef(data, **kwargs)
    alphabet, r_rows = _r_states(case)
    assert seq.alphabet == alphabet
    # R renames one state at a time, so its rows hold fewer distinct states than the input.
    assert len({cell for row in r_rows for cell in row}) < len(
        {cell for row in data for cell in row}
    )
    assert _cells(seq) == expected
    assert _cells(seq) != r_rows


@pytest.mark.parametrize("criterion", ["density", "freq", "dist"])
def test_seqrep_refuses_weights_of_zero_total_as_r_does(criterion):
    assert ITEMS[(f"seqrep_zero_weights_{criterion}", "refusal")] == (
        "missing value where TRUE/FALSE needed"
    )
    seq = seqdef([["A", "B"], ["B", "A"], ["A", "A"]], weights=[0.0, 0.0, 0.0])
    with pytest.raises(ValueError, match="positive total"):
        seqrep(seq, diss=seqdist(seq, method="HAM"), criterion=criterion)


@pytest.mark.parametrize("weights", [[0.0, 0.0, 0.0], [1.0, -1.0, 0.0], [np.nan, 1.0, 1.0]])
def test_dissrep_refuses_weights_without_a_positive_total(weights):
    assert ITEMS[("dissrep_zero_weights", "refusal")] == "missing value where TRUE/FALSE needed"
    with pytest.raises(ValueError, match="positive total"):
        dissrep(np.zeros((3, 3)), weights=np.array(weights))


def test_one_sequence_takes_a_one_element_list_as_its_sort_vector_and_refuses_a_stray_string():
    assert ITEMS[("sortv_one_sequence_unknown_string", "outcome")] == "drawn"
    one = seqdef([["A", "B", "B"]], alphabet=["A", "B"])
    assert resolve_sort_order(one, ["not-a-method"]).tolist() == [0]
    with pytest.raises(ValueError, match="'from.start' or 'from.end'"):
        resolve_sort_order(one, "not-a-method")


def test_a_copy_does_not_share_its_issues_with_the_original():
    seq = seqdef([["A", "B"], ["B", "A"]])
    original = tx.STSList(data=seq.data, alphabet=seq.alphabet, issues={"x": [1]})
    copy = original.copy()
    copy.issues["x"].append(2)
    copy.issues["y"] = [1]
    assert original.issues == {"x": [1]}
