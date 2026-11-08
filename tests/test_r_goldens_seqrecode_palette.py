"""seqrecode keeps R's colours and takes vector labels (TX-AUD-20261011-024).

R's seqrecode gives a recoded state the colour of its first source state, the
otherwise state that of the first original state no recode uses, and keeps a
retained state's colour (seqrecode.R:88-104); a source that is the missing code
gives missing.color. traminex coloured retained and otherwise states grey
(#808080), and refused R's plain-vector labels and colours with AttributeError. Fixture:
``scripts/r_goldens/make_goldens_seqrecode_palette.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "seqrecode_palette.csv",
    dtype=str,
    keep_default_na=False,
)
SEQS = tx.seqdef(
    np.array(
        [
            ["A", "A", "B", "B", "C", "C"],
            ["A", "B", "B", None, "C", "D"],
            ["B", "B", "C", None, None, None],
            ["C", "C", "C", "D", "D", "D"],
        ],
        dtype=object,
    ),
    alphabet=["A", "B", "C", "D"],
)
CASES = {
    "merge_ab": dict(recodes={"X": ["A", "B"]}),
    "merge_ab_otherwise": dict(recodes={"X": ["A", "B"]}, otherwise="Z"),
    "two_recodes": dict(recodes={"X": ["C"], "Y": ["D"]}),
    "missing_source": dict(recodes={"M": ["*"], "X": ["A", "B"]}),
    "vector_labels": dict(recodes={"X": ["A", "B"]}, labels=["x", "c", "d"]),
    "vector_cpal": dict(recodes={"X": ["A", "B"]}, cpal=["red", "green", "blue"]),
    "nr_named_last": dict(recodes={"X": ["A", "B"], "*": ["C"]}),
    "void_source": dict(recodes={"X": ["A", "B"], "Y": ["%"]}),
}
WARNINGS = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "seqrecode_palette_warnings.csv",
    dtype=str,
    keep_default_na=False,
)


def _real(x: tx.STSList) -> list[tuple[str, str, str]]:
    markers = {x.nr_code, x.void_code} - {None}
    labels = x.labels if x.labels is not None else x.alphabet
    return [
        (x.alphabet[code], x.cpal[code], labels[code])
        for code in range(len(x.alphabet))
        if code not in markers
    ]


def _r(case: str) -> list[tuple[str, str, str]]:
    rows = GOLDEN[GOLDEN["case"] == case].sort_values("position", key=lambda c: c.astype(int))
    return list(zip(rows["state"], rows["colour"], rows["label"]))


@pytest.mark.parametrize("case", sorted(CASES))
def test_seqrecode_states_colours_and_labels_match_r(case: str) -> None:
    if case == "void_source":
        message = WARNINGS.set_index("case").loc[case, "warning"]
        with pytest.warns(
            UserWarning, match=message.replace("[!]", r"\[!\]").replace("(%)", r"\(%\)")
        ):
            recoded = tx.seqrecode(SEQS, **CASES[case])
    else:
        recoded = tx.seqrecode(SEQS, **CASES[case])
    assert _real(recoded) == _r(case)


@pytest.mark.parametrize("recodes", [{"X": ["A", "B"], "*": ["C"]}, {"*": ["C"], "X": ["A", "B"]}])
def test_the_marker_slots_keep_their_colours(recodes) -> None:
    recoded = tx.seqrecode(SEQS, recodes)
    assert recoded.cpal[recoded.nr_code] == SEQS.missing_color
    assert recoded.cpal[recoded.void_code] == "white"


def test_a_recode_named_like_nr_does_not_shift_the_others_colours() -> None:
    # By name (the definition): X takes A's colour, D keeps its own. R picks by
    # position, so a first recode named "*" gives X the colour of C.
    recoded = tx.seqrecode(SEQS, {"*": ["C"], "X": ["A", "B"]})
    by_name = [("X", SEQS.cpal[0], "X"), ("D", SEQS.cpal[3], "D")]
    assert _real(recoded) == by_name
    assert _r("nr_named_first") != by_name
    assert _r("nr_named_first")[0][1] == SEQS.cpal[2]  # R: X gets C's colour


def test_a_vector_cpal_skips_the_marker_slots() -> None:
    recoded = tx.seqrecode(SEQS, {"X": ["A", "B"], "*": ["C"]}, cpal=["red", "blue"])
    assert _real(recoded) == [("X", "red", "X"), ("D", "blue", "D")]
