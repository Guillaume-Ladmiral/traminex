"""Batch-18 ledger items against R (``make_goldens_ledger_items.R``; TX-AUD-20261011-005).

Five ledger rows said EXACT although traminex departs from R. ``gower_matrix`` refused a
weight total of zero or below where R's formula (dissmfacw.R:11-22) gives NaN or finite
values: traminex now returns R's values. The other four departures are recorded in the
ledger as DIVERGENCE and tested here against values derived from their definitions, with
R's values recorded beside them.
"""

from __future__ import annotations

import csv
import warnings
from pathlib import Path

import numpy as np
from r_goldens_utils import R_GOLDEN_ATOL

import traminex as tx

GOLDEN_PATH = Path(__file__).resolve().parent / "r_goldens" / "ledger_items.csv"
with open(GOLDEN_PATH, encoding="utf-8", newline="") as handle:
    GOLDEN = {(row["case"], row["field"]): row["value"] for row in csv.DictReader(handle)}


def _numbers(case: str, field: str) -> np.ndarray:
    return np.array([float(v) if v != "NA" else np.nan for v in GOLDEN[(case, field)].split()])


def _texts(case: str, field: str) -> list[str]:
    return GOLDEN[(case, field)].split(" ;; ")


def _gower(weights: list[float]) -> np.ndarray:
    D = _numbers("gower_totals", "D").reshape(3, 3)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        return tx.gower_matrix(D, weights=np.array(weights, dtype=float))


def test_gower_matrix_gives_r_values_for_a_zero_weight_total() -> None:
    result = _gower([1, 1, -2])

    expected = _numbers("gower_totals", "zero_total").reshape(3, 3)
    assert np.isnan(expected).all()
    assert np.isnan(result).all()


def test_gower_matrix_gives_r_values_for_a_negative_weight_total() -> None:
    result = _gower([-1, -1, -2])

    expected = _numbers("gower_totals", "negative_total").reshape(3, 3)
    np.testing.assert_allclose(result, expected, rtol=0, atol=R_GOLDEN_ATOL)


# The four departures below are recorded as DIVERGENCE (rulings R18-j, owner decision of 2026-10-24
# for seqecontain). Each test asserts traminex's value from the definition and that R's differs.

ALPHABET = ["A", "B", "C"]


def _cells(case: str, field: str = "cells") -> list[list[str | None]]:
    return [[None if c == "NA" else c for c in row.split("-")] for row in _texts(case, field)]


def test_is_stslist_accepts_a_seqdef_object_with_ids_and_weights() -> None:
    # A seqdef result is a state-sequence object; R's is.stslist rejects it because seqdef names the
    # weights 1..n and the rows by id, and is.stslist requires the names to match.
    s = tx.seqdef(
        _cells("is_stslist"),
        alphabet=ALPHABET,
        id=_texts("is_stslist", "row_names"),
        weights=[1, 2, 3, 4],
    )

    assert tx.is_stslist(s) is True
    assert GOLDEN[("is_stslist", "id_and_weights")] == "FALSE"
    assert _texts("is_stslist", "weight_names") != _texts("is_stslist", "row_names")
    assert (GOLDEN[("is_stslist", "weights_only")], GOLDEN[("is_stslist", "id_only")]) == (
        "TRUE",
        "TRUE",
    )


def _recoded_rows(cells: list[list[str | None]], case: str) -> list:
    s = tx.seqrecode(tx.seqdef(cells, alphabet=ALPHABET), recodes={"%": "B"})
    shown = ["-".join(s.alphabet[c] if c != s.void_code else "%" for c in row) for row in s.data]
    assert shown == _texts(case, "cells")
    rows = []
    for i in range(len(s.data)):
        row = s.copy()
        row.data = s.data[i : i + 1]
        rows.append(row)
    return rows


# Without the voids the rows are A A C, A C C, C C A and C C A (by hand from the recorded cells).
_LLCS_BY_DEFINITION = np.array([[3, 2, 1, 1], [2, 3, 2, 2], [1, 2, 3, 3], [1, 2, 3, 3]])
_LLCP_BY_DEFINITION = np.array([[3, 1, 0, 0], [1, 3, 0, 0], [0, 0, 3, 3], [0, 0, 3, 3]])


def test_seqllcs_and_seqllcp_compare_the_non_void_states() -> None:
    rows = _recoded_rows(_cells("is_stslist"), "llc_interior_void")
    llcs = np.array([[tx.seqLLCS(a, b) for b in rows] for a in rows])
    llcp = np.array([[tx.seqLLCP(a, b) for b in rows] for a in rows])

    np.testing.assert_array_equal(llcs, _LLCS_BY_DEFINITION)
    np.testing.assert_array_equal(llcp, _LLCP_BY_DEFINITION)
    # R compares the first seqlength() cells, so the void counts as a state and a state drops out.
    r_llcs = _numbers("llc_interior_void", "seqLLCS").reshape(4, 4)
    r_llcp = _numbers("llc_interior_void", "seqLLCP").reshape(4, 4)
    assert (r_llcs[1, 3], r_llcs[2, 3], r_llcp[2, 3]) == (1, 2, 0)
    assert not np.array_equal(r_llcs, _LLCS_BY_DEFINITION)
    assert not np.array_equal(r_llcp, _LLCP_BY_DEFINITION)


# A missing cell beside each interior void (batch-18 review M-3). Without the voids the rows are
# A * C, A * C, * C A and C * A, the missing state "*" kept as a state (by hand from the recorded
# cells).
_LLCS_MISSING_BY_DEFINITION = np.array([[3, 3, 2, 1], [3, 3, 2, 1], [2, 2, 3, 2], [1, 1, 2, 3]])
_LLCP_MISSING_BY_DEFINITION = np.array([[3, 3, 0, 0], [3, 3, 0, 0], [0, 0, 3, 0], [0, 0, 0, 3]])


def test_seqllcs_and_seqllcp_keep_a_missing_cell_beside_an_interior_void() -> None:
    rows = _recoded_rows(_cells("llc_void_beside_missing", "input"), "llc_void_beside_missing")
    llcs = np.array([[tx.seqLLCS(a, b) for b in rows] for a in rows])
    llcp = np.array([[tx.seqLLCP(a, b) for b in rows] for a in rows])

    np.testing.assert_array_equal(llcs, _LLCS_MISSING_BY_DEFINITION)
    np.testing.assert_array_equal(llcp, _LLCP_MISSING_BY_DEFINITION)
    # R's first seqlength() cells take the void in and leave a state out here too.
    r_llcs = _numbers("llc_void_beside_missing", "seqLLCS").reshape(4, 4)
    r_llcp = _numbers("llc_void_beside_missing", "seqLLCP").reshape(4, 4)
    assert (r_llcs[0, 1], r_llcp[0, 1]) == (2, 0)
    assert not np.array_equal(r_llcs, _LLCS_MISSING_BY_DEFINITION)
    assert not np.array_equal(r_llcp, _LLCP_MISSING_BY_DEFINITION)


def test_seqecontain_takes_a_single_event_sequence() -> None:
    es = tx.seqecreate(id=[1, 1, 2], timestamp=[0, 1, 0], event=["a", "b", "b"])

    assert tx.seqecontain(es[0], "a") is True
    assert GOLDEN[("seqecontain_single", "single")].startswith("ERROR: eseq should be a seqelist")
    listed = [str(bool(v)).upper() for v in tx.seqecontain(es, "a")]
    assert listed == _texts("seqecontain_single", "list")


def test_the_four_departures_are_recorded_in_the_ledger() -> None:
    with open(
        Path(__file__).resolve().parents[1] / "docs" / "API_PARITY_DECISIONS.csv", encoding="utf-8"
    ) as handle:
        rows = {row["r_name"]: row for row in csv.DictReader(handle)}

    for name in ("is.stslist", "seqLLCS", "seqLLCP", "seqecontain"):
        assert rows[name]["status"] == "DIVERGENCE", name
        assert "scripts/r_goldens/make_goldens_ledger_items.R" in rows[name]["evidence"], name
        assert "tests/test_ledger_items.py" in rows[name]["evidence"], name


def test_the_four_docstrings_do_not_call_the_port_faithful() -> None:
    # Each is recorded as DIVERGENCE; its docstring names the departure (batch-18 review M-10).
    for function in (tx.is_stslist, tx.seqLLCS, tx.seqLLCP, tx.seqecontain):
        doc = " ".join((function.__doc__ or "").split())

        assert "faithful" not in doc.lower(), function.__name__
        assert "TX-AUD-20261011-005" in doc, function.__name__
