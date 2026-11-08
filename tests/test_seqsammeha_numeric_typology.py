"""``seqsammeha`` names a numeric typology's columns and levels as R prints them.

``seqsammeha`` (TraMineRextras ``seqsamm.R``) stores ``as.character(typology)``
and names the indicator columns from the typology values, so a typology given
as the numbers 1, 2 gives columns ``SAMM1``/``SAMM2`` and levels ``None 1 2``.
traminex used Python's ``str``, so the same typology passed as floats (as a
NumPy or pandas column usually is) gave ``SAMM1.0``/``SAMM2.0`` and levels
``1.0 2.0`` (2026-09-22 audit S5-12).

R fills the indicator columns by position (``ret[, tt]`` with a numeric
``tt``), so a numeric typology whose values are not 1, 2, ... in order of
first appearance gets mislabelled indicators (2, 1 swaps them) or fails
(10, 20: "subscript out of bounds"). traminex matches by value; that case is
a recorded divergence and is pinned here against R's output.

Expected values: ``scripts/r_goldens/make_goldens_audit20260922_s5.R``
(TraMineR 2.2-14, TraMineRextras 0.6.9).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import traminex as tx

GOLDEN = Path(__file__).parent / "r_goldens" / "audit20260922_s5_samm_typology.txt"
LINES = [line.split("\t") for line in GOLDEN.read_text().splitlines() if not line.startswith("#")]
EHA = {fields[1]: fields[2:] for fields in LINES if fields[0] == "eha"}
INDICATOR = {(fields[1], fields[2]): fields[3] for fields in LINES if fields[0] == "indicator"}
TYPOLOGIES = {
    "whole": np.array([1.0, 2.0, 1.0, 2.0]),
    "fractional": np.array([1.5, 2.5, 1.5, 2.5]),
    "reversed": np.array([2.0, 1.0, 2.0, 1.0]),
}


def _samm():
    rows = ["AABBCCA", "BBAAACC", "CABBBCC", "ABCABCA", "BCCCAAA"]
    sequences = tx.seqdef([list(row) for row in rows], alphabet=["A", "B", "C"])
    return tx.seqsamm(sequences, sublength=3)


def _indicators(frame):
    return [column for column in frame.columns if column.startswith("SAMM")][1:]


@pytest.mark.parametrize("tag", ["whole", "fractional"])
def test_numeric_typology_names_follow_r(tag):
    result = tx.seqsammeha(_samm(), spell="A", typology=TYPOLOGIES[tag])
    columns, levels, values = EHA[tag]
    assert _indicators(result) == columns.split(" ")
    assert list(result["SAMMtypology"].cat.categories) == levels.split(" ")
    assert [str(value) for value in result["SAMMtypology"]] == values.split(" ")
    for column in _indicators(result):
        np.testing.assert_array_equal(
            result[column].to_numpy(), np.array(INDICATOR[(tag, column)].split(" "), dtype=float)
        )


def test_indicators_are_matched_by_value_where_r_matches_by_position():
    result = tx.seqsammeha(_samm(), spell="A", typology=TYPOLOGIES["reversed"])
    columns, levels, values = EHA["reversed"]
    assert _indicators(result) == columns.split(" ") == ["SAMM2", "SAMM1"]
    assert list(result["SAMMtypology"].cat.categories) == levels.split(" ")
    typology = [str(value) for value in result["SAMMtypology"]]
    assert typology == values.split(" ")
    for column in ("SAMM1", "SAMM2"):
        expected = np.array([value == column[4:] for value in typology], dtype=float)
        np.testing.assert_array_equal(result[column].to_numpy(), expected)
    # R's SAMM2 column holds the indicator of typology 1 (and vice versa).
    np.testing.assert_array_equal(
        np.array(INDICATOR[("reversed", "SAMM2")].split(" "), dtype=float),
        result["SAMM1"].to_numpy(),
    )
