"""Event labels and group levels follow R's collation (TX-AUD-20261011-015, -017).

The package's portable label key ordered text by code point after folding case,
so seqecreate listed 'A>%' before 'A>*' and 'é' after 'f', which changed
seqefsub's tie order; group_compare, samm and group_p used the same key, and
with dissCompare's draw rule (the first level is drawn from) that changed a
sampled statistic. R orders with levels(factor()), i.e. ICU collation.
Fixtures: ``scripts/r_goldens/make_goldens_seqecreate_label_order.R`` and
``make_goldens_factor_level_order.R`` (R's levels of twelve labels).
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import traminex as tx
from traminex import group_compare
from traminex._categorical_order import stable_sorted_unique

GOLDENS = Path(__file__).resolve().parent / "r_goldens"


def _read(name: str) -> dict[tuple[str, str], list[str]]:
    out: dict[tuple[str, str], list[str]] = defaultdict(list)
    with (GOLDENS / name).open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            out[(row["case"], row["kind"])].append(row["value"])
    return out


R = _read("seqecreate_label_order.csv")
R_LEVELS = _read("factor_level_order.csv")[("labels", "levels")]
CELLS = R[("transitions", "cells")]


def _transitions() -> tx.STSList:
    rows = [[None if c in ("*", "%") else c for c in CELLS[4 * i : 4 * i + 4]] for i in range(4)]
    return tx.seqdef(rows, alphabet=["A", "B", "C"], right="DEL")


def test_transition_dictionary_matches_r() -> None:
    seqs = _transitions()
    assert [s for row in seqs.data for s in (seqs._state_for_code(int(c)) for c in row)] == CELLS
    assert list(tx.seqecreate(seqs).dictionary) == R[("transitions", "dictionary")]


def test_seqefsub_order_matches_r() -> None:
    mined = tx.seqefsub(tx.seqecreate(_transitions()), pmin_support=0.25)
    assert list(mined.subseq) == R[("transitions", "seqefsub")]


def test_label_dictionary_matches_r() -> None:
    events = tx.seqecreate(
        id=list(range(1, 9)),
        timestamp=[0] * 8,
        event=["e", "é", "f", "_x", "B", "a", "A>%", "A>*"],
    )
    assert list(events.dictionary) == R[("labels", "dictionary")]


def test_text_levels_take_r_order() -> None:
    labels = ["b", "a", "B", "_c", "é", "e", "A>%", "A>*", "10", "9", "faible", "élevé"]
    assert stable_sorted_unique(labels) == R_LEVELS


def test_dissCompare_draws_from_r_first_group() -> None:
    # R's factor(c("faible", "élevé")) puts "élevé" first; dissCompare draws its
    # "first group" sample from that level (TX-AUD-20261011-022's draw rule).
    numpy = __import__("numpy")
    levels, _ = group_compare._compare_factor(
        numpy.array(["faible", "élevé"], dtype=object), numpy.array([True, True])
    )
    assert levels == ["élevé", "faible"]
