"""seqdef pairs states, labels and colours with R's default order (TX-AUD-20261010-001).

Without ``alphabet=``, TraMineR pairs ``states=``, ``labels=``, ``cpal=`` and its
default palette with the order ``seqstatl`` returns (numeric when every state
is an R number, R's collation otherwise). Since owner decision D2 (2026-10-19)
traminex stores that same order, so the alphabet equals R's by position and each
*state* receives R's name, label and colour, against the fixtures written by
``scripts/r_goldens/make_goldens_seqdef_r_pairing.R``.

The two ``defect_*`` cases record an R defect: R re-renders numeric levels with
``as.character`` and loses cells spelled ``01`` or ``0x10``. There the expected
value comes from the definition, and the R cells are checked to stay recorded.
"""

from __future__ import annotations

import csv
from collections.abc import Callable
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx
from traminex import seqdef

GOLDENS = Path(__file__).resolve().parent / "r_goldens"

BIOFAM_LABELS = [
    "Parent",
    "Left",
    "Married",
    "Left+Marr",
    "Child",
    "Left+Child",
    "Left+Marr+Child",
    "Divorced",
]
BIOFAM_STATES = ["P", "L", "M", "LM", "C", "LC", "LMC", "D"]

# The arguments of each TraMineR::seqdef call in the generator (facts of the
# call, not expected values) and how its input cells are typed: biofam is an
# integer matrix in R, every other input a character matrix.
CALLS: dict[str, tuple[Callable[[str], object], dict[str, object]]] = {
    "biofam_labels_states": (int, {"labels": BIOFAM_LABELS, "states": BIOFAM_STATES}),
    "char_states_labels_cpal": (
        str,
        {"states": ["x", "y", "z"], "labels": ["Lx", "Ly", "Lz"], "cpal": ["red", "green", "blue"]},
    ),
    "char_default": (str, {}),
    "numeric_text_states": (str, {"states": ["two", "nine", "ten"]}),
    "underscore_text_labels": (str, {"labels": ["L1", "L2", "L3"]}),
    "left_policy_labels_cpal": (
        str,
        {"left": "LL", "labels": ["La", "Lb", "Lc", "Lleft"], "cpal": ["c1", "c2", "c3", "c4"]},
    ),
    # seqstatl sees the missing code too: a text code gives character order.
    "missing_text_code_states": (
        str,
        {"missing": "x", "states": ["s1", "s2", "s3"], "labels": ["L1", "L2", "L3"]},
    ),
    "missing_numeric_code_states": (
        str,
        {"missing": "99", "states": ["s1", "s2", "s3"], "labels": ["L1", "L2", "L3"]},
    ),
}
# Spellings of seqdef's default nonresponse and void markers in R's output.
R_MARKERS = {"*", "%"}


def _read(kind: str) -> list[dict[str, str]]:
    with (GOLDENS / f"seqdef_pairing_{kind}.csv").open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


INPUTS = _read("inputs")
CELLS = _read("cells")
STATES = _read("states")


def _matrix(
    rows: list[dict[str, str]], case: str, field: str, convert: Callable[[str], object] = str
) -> np.ndarray:
    cells = {(int(row["row"]), int(row["col"])): row[field] for row in rows if row["case"] == case}
    assert cells, f"no golden rows for {case}"
    shape = (max(i for i, _ in cells), max(j for _, j in cells))
    out = np.empty(shape, dtype=object)
    for (i, j), value in cells.items():
        out[i - 1, j - 1] = None if value == "<NA>" else convert(value)
    return out


def _r_states(case: str) -> list[dict[str, str]]:
    return sorted(
        (row for row in STATES if row["case"] == case), key=lambda row: int(row["position"])
    )


def _build(case: str):
    convert, kwargs = CALLS[case]
    return seqdef(_matrix(INPUTS, case, "value", convert), **kwargs)


def _decoded(seq) -> list[list[str]]:
    return [[str(seq.alphabet[int(code)]) for code in row] for row in np.asarray(seq.data)]


def _real_alphabet(seq) -> list[str]:
    markers = {seq.nr_code, seq.void_code}
    return [str(state) for code, state in enumerate(seq.alphabet) if code not in markers]


@pytest.mark.parametrize("case", sorted(CALLS))
def test_decoded_cells_match_r(case: str):
    assert _decoded(_build(case)) == _matrix(CELLS, case, "state").tolist()


@pytest.mark.parametrize("case", sorted(CALLS))
def test_every_state_gets_r_label_and_colour(case: str):
    seq = _build(case)
    alphabet = [str(state) for state in seq.alphabet]
    for row in _r_states(case):
        position = alphabet.index(row["state"])
        assert (row["state"], seq.labels[position], seq.cpal[position]) == (
            row["state"],
            row["label"],
            row["cpal"],
        )


@pytest.mark.parametrize("case", sorted(CALLS))
def test_alphabet_is_r_order(case: str):
    # Owner decision D2 (2026-10-19): the stored alphabet is R's, position by
    # position, including the policy states R appends (seqdef.R:59-66, 116-132).
    # Until then this test pinned first-appearance order (TX-AUD-20261019-071).
    assert _real_alphabet(_build(case)) == [row["state"] for row in _r_states(case)]


def test_noncanonical_numeric_spelling_keeps_the_cell_state():
    # seqstatl orders 01 (the number 1) before 2, so states= gives 01 -> one.
    seq = seqdef(np.array([["2", "01", "2"]], dtype=object), states=["one", "two"])
    assert _decoded(seq) == [["two", "one", "two"]]
    # R loses the cell (its level is re-rendered as "1"); the divergence stays recorded.
    assert _matrix(CELLS, "defect_noncanonical_states", "state").tolist() == [["two", None, "two"]]


def test_hexadecimal_spelling_is_ordered_by_its_numeric_value():
    # R's as.numeric reads 0x10 as 16, so the order is 2, 9, 0x10.
    seq = seqdef(np.array([["0x10", "9", "2"]], dtype=object), labels=["L1", "L2", "L3"])
    assert dict(zip(seq.alphabet, seq.labels)) == {"2": "L1", "9": "L2", "0x10": "L3"}
    assert _matrix(CELLS, "defect_hex_text_labels", "state").tolist() == [[None, "9", "2"]]


def test_states_of_another_length_than_the_alphabet_are_refused():
    # seqdef.Rd: states "Must have a length equal to the size of the alphabet".
    # R does not check: it pairs the first entries with the observed states and
    # keeps the rest as unobserved states (recorded, divergence kept visible).
    assert _matrix(CELLS, "defect_states_longer", "state").tolist() == [["P", "L", "M"]]
    assert [row["state"] for row in _r_states("defect_states_longer")] == ["P", "L", "M", "LM"]
    with pytest.raises(
        ValueError, match=r"'states' must have the same length as the alphabet \(3 states\), got 4"
    ):
        seqdef(np.array([["0", "1", "3"]], dtype=object), states=["P", "L", "M", "LM"])


def test_marker_labels_after_the_real_states_keep_their_place():
    # traminex also accepts labels for the markers after the real states.
    seq = seqdef(np.array([[None, "C", "A"]], dtype=object), labels=["La", "Lc", "Lnr"])
    by_state = dict(zip(seq.alphabet, seq.labels))
    assert (by_state["A"], by_state["C"]) == ("La", "Lc")
    assert seq.labels[seq.nr_code] == "Lnr"


@pytest.mark.parametrize("constructor", ["seqdef_from_long", "TSE_to_STS"])
def test_long_constructors_pair_labels_and_colours_as_seqdef(constructor: str):
    # TX-REP-20261012-001: both constructors document seqdef's order, so the
    # same data must give each state the same label and colour as seqdef;
    # for these character states R's order is alphabetical (A, B, C).
    rows = [["C", "A", "B", "A"], ["B", "C", "C", "A"]]
    long = pd.DataFrame(
        [
            {"id": person + 1, "time": time + 1, "state": state}
            for person, row in enumerate(rows)
            for time, state in enumerate(row)
        ]
    )
    options = {"labels": ["Lx", "Ly", "Lz"], "cpal": ["red", "green", "blue"]}
    build = getattr(tx, constructor)
    seq = build(long, id_col="id", time_col="time", state_col="state", **options)
    pairs = {
        state: (label, colour) for state, label, colour in zip(seq.alphabet, seq.labels, seq.cpal)
    }
    expected = {"A": ("Lx", "red"), "B": ("Ly", "green"), "C": ("Lz", "blue")}
    assert {state: pairs[state] for state in expected} == expected
    wide = seqdef(np.array(rows, dtype=object), **options)
    assert dict(zip(wide.alphabet, zip(wide.labels, wide.cpal))) == expected
