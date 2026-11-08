"""Event sequences keep TraMineR's gap chain when event times are fractional.

TraMineR builds each event sequence by inserting its events one at a time, in
their input order, into a chain of gaps (``tmrsequence.cpp:113-127``,
``Sequence::addEvent`` in ``eventseq.cpp:21-48``); every count method, the
event ages and the printed gaps are then read back from that chain. traminex
rebuilt the chain from the differences of the sorted absolute times. With
whole-number times the two agree exactly, but with fractional times the
repeated subtractions round differently: R's ``(c)`` at 2.0 after events at
0.1 and 0.2 has age 1.9999999999999998, so method ``"age"`` reports 1 where
traminex reported 2, and CDIST/CWIN counts differed as well (2026-09-22 audit
S5-01).

Expected values: ``scripts/r_goldens/make_goldens_audit20260922_s5.R``
(TraMineR 2.2-14, TraMineRextras 0.6.9).
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import numpy as np
import pytest

import traminex as tx
from traminex.events import seqeage

GOLDEN = Path(__file__).parent / "r_goldens" / "audit20260922_s5_fractional_event_times.txt"


def _load() -> dict[str, list[list[str]]]:
    rows: dict[str, list[list[str]]] = defaultdict(list)
    for line in GOLDEN.read_text().splitlines():
        if line.startswith("#"):
            continue
        kind, *fields = line.split("\t")
        rows[kind].append(fields)
    return rows


R = _load()


def _numbers(text: str) -> np.ndarray:
    return np.array([] if text == "none" else text.split(" "), dtype=float)


def _corpus(name: str) -> tx.SeqElist:
    rows = [fields for fields in R["corpus"] if fields[0] == name]
    return tx.seqecreate(
        id=[int(float(fields[1])) for fields in rows],
        timestamp=[float(fields[2]) for fields in rows],
        event=[fields[3] for fields in rows],
    )


CORPORA = ["minimal_age", "minimal_overlap", "shuffled"]


@pytest.mark.parametrize("name", CORPORA)
def test_printed_gaps_follow_the_insertion_chain(name):
    eseq = _corpus(name)
    expected = [fields[2] for fields in R["aschar"] if fields[0] == name]
    (levels,) = [fields[1] for fields in R["levels"] if fields[0] == name]
    assert eseq.dictionary == levels.split(" ")
    assert eseq.as_character() == expected


@pytest.mark.parametrize("name", CORPORA)
def test_event_ages_are_the_sums_of_the_chain_gaps(name):
    eseq = _corpus(name)
    for fields in R["seqeage"]:
        if fields[0] != name:
            continue
        code = int(fields[1])
        flat = _numbers(fields[2])
        expected = flat.reshape((len(eseq), -1), order="F")
        np.testing.assert_array_equal(seqeage(eseq, code), expected, err_msg=f"event {code}")


def _applied_constraints():
    cwin_small = tx.seqeconstraint(window_size=1, count_method="CWIN")
    cwin = tx.seqeconstraint(window_size=1.1, count_method="CWIN")
    maxgap = tx.seqeconstraint(max_gap=0.3, count_method="CDIST_O")
    # Same order as the ``applied`` list in the generator: (method, constraint).
    return [
        ("age", None),
        ("CDIST", None),
        (None, cwin_small),
        ("age", None),
        ("CDIST", None),
        (None, cwin),
        (None, maxgap),
    ]


@pytest.mark.parametrize("index", range(7))
def test_counts_of_named_subsequences_match_r(index):
    fields = {int(row[0]): row[1:] for row in R["apply"]}[index + 1]
    name, subsequences, flat = fields
    method, constraint = _applied_constraints()[index]
    eseq = _corpus(name)
    fsub = tx.seqefsub(eseq, str_subseq=subsequences.split(" "))
    if constraint is None:
        counts = tx.seqeapplysub(fsub, method=method)
    else:
        counts = tx.seqeapplysub(fsub, constraint=constraint)
    expected = _numbers(flat).reshape(counts.shape, order="F")
    np.testing.assert_array_equal(counts.to_numpy(), expected)


MINED = {
    "cdist": dict(count_method="CDIST"),
    "cwin": dict(window_size=1.1, count_method="CWIN"),
    "maxgap": dict(max_gap=0.7, count_method="CDIST_O"),
    "agemin": dict(age_min=0.3, age_max=4.7, count_method="CDIST_O"),
    "cminwin": dict(window_size=2.3, count_method="CMINWIN"),
}


@pytest.mark.parametrize("tag", sorted(MINED))
def test_mined_subsequences_match_r(tag):
    eseq = _corpus("shuffled")
    fsub = tx.seqefsub(eseq, min_support=3, constraint=tx.seqeconstraint(**MINED[tag]), max_k=3)
    (count,) = [int(fields[1]) for fields in R["mined"] if fields[0] == tag]
    rows = [fields[1:] for fields in R["mined_row"] if fields[0] == tag]
    assert len(fsub.subseq) == count == len(rows)
    assert fsub.subseq == [row[1] for row in rows]
    np.testing.assert_array_equal(
        fsub.data["Support"].to_numpy(), np.array([row[2] for row in rows], dtype=float)
    )
    np.testing.assert_array_equal(
        fsub.data["Count"].to_numpy(), np.array([row[3] for row in rows], dtype=float)
    )


@pytest.mark.parametrize("interval", ["No", "previous", "next"])
def test_event_edit_distance_walks_the_same_chain(interval):
    eseq = _corpus("shuffled")
    (flat,) = [fields[1] for fields in R["seqedist"] if fields[0] == interval]
    result = tx.seqedist(
        eseq, idcost=np.ones(len(eseq.dictionary)), vparam=0.3, interval=interval, norm="none"
    )
    np.testing.assert_array_equal(result, _numbers(flat).reshape(result.shape, order="F"))


def test_whole_number_times_are_unaffected():
    """Integer times subtract exactly, so both constructions give one chain."""
    eseq = tx.seqecreate(
        id=[1, 1, 1, 2, 2], timestamp=[5, 1, 3, 2, 2], event=["c", "a", "b", "b", "a"]
    )
    assert eseq.as_character() == ["1-(a)-2-(b)-2-(c)", "2-(a,b)"]
    np.testing.assert_array_equal(seqeage(eseq, 1), np.array([[1.0], [2.0]]))
