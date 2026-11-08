"""Void positions stay out of the integration weights (TX-AUD-20261019-012).

seqintegr.Rd defines the index for a sequence of length L over positions 1..L, and
``with.missing=TRUE`` makes only non-void missing values a state. R's seqintegr divides by the
weights of every column when ``with.missing=TRUE`` (seqintegr.R:30-34), so the void padding of a
short sequence enters the denominator and its value depends on the other sequences' lengths.
traminex corrects this (DIVERGENCE rows of seqintegr, seqibad, seqipos). R's values are recorded
by ``scripts/r_goldens/make_goldens_audit20260922_s3b.R`` (cases ``s8_012_*``); the expected
values are worked by hand.
"""

from __future__ import annotations

import csv
import warnings
from pathlib import Path

import numpy as np
import pytest

import traminex as tx

GOLDEN = Path(__file__).parent / "r_goldens" / "audit20260922_s3b.txt"
ROWS = [["a", "b", "b", "a"], ["a", "b", None, None], ["b", None, "a", "a"]]


def _r(case: str) -> np.ndarray:
    for line in GOLDEN.read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        if fields[0] == case:
            assert fields[1] == "value", line
            return np.array([float(v) for v in fields[2].split()])
    raise AssertionError(case)


def _quiet(function, *args, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return function(*args, **kwargs)


def _x(rows=ROWS):
    return _quiet(tx.seqdef, rows)


def test_seqintegr_counts_only_the_sequence_positions() -> None:
    # Row 2 is a-b and two void cells: positions 1 and 2 weigh 1 and 2 of 3. Rows 1 and 3 have no
    # void (row 3's nonresponse at position 2 is a state with include_missing=True).
    expected = np.array([[0.5, 0.5, 0.0], [1 / 3, 2 / 3, 0.0], [0.7, 0.1, 0.2]])
    actual = np.asarray(_quiet(tx.seqintegr, _x(), include_missing=True), dtype=float)
    np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-12)
    r_values = _r("s8_012_integr_with_missing").reshape(3, 3, order="F")
    np.testing.assert_allclose(r_values[[0, 2]], expected[[0, 2]], rtol=0, atol=1e-12)
    np.testing.assert_allclose(r_values[1], [0.1, 0.2, 0.0], rtol=0, atol=1e-12)


def test_the_value_of_a_sequence_does_not_depend_on_the_data_width() -> None:
    narrow = np.asarray(
        _quiet(tx.seqintegr, _x([["a", "b", None], ["a", "b", "b"]]), include_missing=True), float
    )
    np.testing.assert_allclose(narrow[0], [1 / 3, 2 / 3, 0.0], rtol=0, atol=1e-12)
    r_values = _r("s8_012_integr_width3_with_missing").reshape(2, 3, order="F")
    np.testing.assert_allclose(r_values[0], [1 / 6, 1 / 3, 0.0], rtol=0, atol=1e-12)


def test_without_missing_values_are_unchanged_and_r_s() -> None:
    actual = np.asarray(_quiet(tx.seqintegr, _x(), include_missing=False), dtype=float)
    np.testing.assert_allclose(
        actual.ravel(order="F"), _r("s8_012_integr_without_missing"), rtol=0, atol=1e-12
    )


def test_seqibad_and_seqipos_follow() -> None:
    # stprec of a, b, * is 0, 0.5, 1 (R's equally spaced default); row 2: 2/3 * 0.5 = 1/3.
    ibad = np.asarray(_quiet(tx.seqibad, _x(), include_missing=True), dtype=float)
    np.testing.assert_allclose(ibad, [0.25, 1 / 3, 0.25], rtol=0, atol=1e-12)
    np.testing.assert_allclose(
        _r("s8_012_ibad_with_missing"), [0.25, 0.1, 0.25], rtol=0, atol=1e-12
    )
    ipos = np.asarray(
        _quiet(tx.seqipos, _x(), pos_states=["a"], index="integr", include_missing=True), float
    ).ravel()
    np.testing.assert_allclose(ipos, [0.5, 1 / 3, 0.7], rtol=0, atol=1e-12)
    np.testing.assert_allclose(
        _r("s8_012_ipos_integr_with_missing"), [0.5, 0.1, 0.7], rtol=0, atol=1e-12
    )


@pytest.mark.parametrize("pow_", [0.5, 2.0])
def test_recency_exponent_uses_the_same_positions(pow_: float) -> None:
    actual = np.asarray(_quiet(tx.seqintegr, _x(), include_missing=True, pow=pow_), dtype=float)
    weights = np.array([1.0, 2.0]) ** pow_
    np.testing.assert_allclose(actual[1, :2], weights / weights.sum(), rtol=0, atol=1e-12)


EMPTY = [["a", "b", None, None], [None, None, None, None], ["b", "b", "a", "b"]]


def test_an_empty_sequence_has_no_positions() -> None:
    # Row 2 is all missing, so seqdef voids it: no positions, the index is 0/0. R gives NaN without
    # missing values and 0 with them, from the every-column denominator corrected above
    # (TX-REP-20261025-006). Rows 1 and 3: a-b over positions 1-2, b-b-a-b over 1-4; stprec a 0,
    # b 1 without missing values, a 0, b 0.5, * 1 with them.
    x = _quiet(tx.seqdef, EMPTY, alphabet=["a", "b"])
    nan = np.nan
    for include_missing, ibad in ((False, [2 / 3, nan, 0.7]), (True, [1 / 3, nan, 0.35])):
        integr = _quiet(tx.seqintegr, x, state="a", include_missing=include_missing)
        np.testing.assert_allclose(np.ravel(integr), [1 / 3, nan, 0.3], rtol=0, atol=1e-12)
        ipos = _quiet(
            tx.seqipos, x, pos_states=["a"], index="integr", include_missing=include_missing
        )
        np.testing.assert_allclose(np.ravel(ipos), [1 / 3, nan, 0.3], rtol=0, atol=1e-12)
        actual = _quiet(tx.seqibad, x, include_missing=include_missing)
        np.testing.assert_allclose(np.ravel(actual), ibad, rtol=0, atol=1e-12)
    for case, expected in (
        ("s8_rep006_integr_empty_without_missing", [1 / 3, nan, 0.3]),
        ("s8_rep006_integr_empty_with_missing", [0.1, 0.0, 0.3]),
        ("s8_rep006_ibad_empty_without_missing", [2 / 3, nan, 0.7]),
        ("s8_rep006_ibad_empty_with_missing", [0.1, 0.0, 0.35]),
        ("s8_rep006_ipos_integr_empty_with_missing", [0.1, 0.0, 0.3]),
    ):
        np.testing.assert_allclose(_r(case), expected, rtol=0, atol=1e-12, err_msg=case)


def test_the_ledger_notes_state_the_empty_sequence() -> None:
    ledger = Path(__file__).resolve().parents[1] / "docs" / "API_PARITY_DECISIONS.csv"
    with open(ledger, encoding="utf-8") as handle:
        notes = {row["r_name"]: row["notes"] for row in csv.DictReader(handle)}

    for name in ("seqintegr", "seqibad", "seqipos"):
        assert "An empty sequence (all void) gives NaN" in notes[name], name
        assert "TX-REP-20261025-006" in notes[name], name
    assert "all-void rows without missing values" in notes["seqintegr"]


def test_seqindic_integration_columns_follow_seqintegr_and_seqibad() -> None:
    # seqindic's Integr and Bad columns are seqintegr and seqibad, so they carry the same correction
    # and the same empty-sequence value (batch-19 review M-1). R: cases s8_rep006_seqindic_*.
    x = _quiet(tx.seqdef, EMPTY, alphabet=["a", "b"])
    nan = np.nan
    for with_missing, bad in ((False, [2 / 3, nan, 0.7]), (True, [1 / 3, nan, 0.35])):
        table = _quiet(
            tx.seqindic, x, indic=["integr", "bad"], ipos_args={"pos_states": ["a"]},
            with_missing=with_missing,
        )  # fmt: skip
        np.testing.assert_allclose(table["Integr"], [1 / 3, nan, 0.3], rtol=0, atol=1e-12)
        np.testing.assert_allclose(table["Bad"], bad, rtol=0, atol=1e-12)
    # R's matrices, column by column (Integr, then Bad).
    np.testing.assert_allclose(
        _r("s8_rep006_seqindic_empty_without_missing"),
        [1 / 3, nan, 0.3, 2 / 3, nan, 0.7],
        rtol=0,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        _r("s8_rep006_seqindic_empty_with_missing"),
        [0.1, 0.0, 0.3, 0.1, 0.0, 0.35],
        rtol=0,
        atol=1e-12,
    )


def test_the_seqindic_note_and_changelog_state_the_correction() -> None:
    root = Path(__file__).resolve().parents[1]
    with open(root / "docs" / "API_PARITY_DECISIONS.csv", encoding="utf-8") as handle:
        notes = {row["r_name"]: row["notes"] for row in csv.DictReader(handle)}
    changelog = " ".join((root / "CHANGELOG.md").read_text(encoding="utf-8").split())

    assert "TX-AUD-20261019-012" in notes["seqindic"]
    assert "An empty sequence (all void) gives NaN" in notes["seqindic"]
    assert "`seqindic`'s `Integr` and `Bad` columns with `with_missing=True`" in changelog
