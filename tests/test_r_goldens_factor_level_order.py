"""Text group levels follow R's factor() order (TX-AUD-20261010-039).

R orders a character group with ``factor()``, i.e. by collation (ICU), so
dissassoc's group rows, disstree's split levels, seqrep.grp's groups and
dissmergegroups' group numbers all list ``_c, a, B, é``. traminex sorted them
by code point (``B, _c, a, é``). Numbers, booleans and declared pandas
categories keep their own order. Fixture:
``scripts/r_goldens/make_goldens_factor_level_order.R`` (with R's distance matrix).
"""

from __future__ import annotations

import csv
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx
from traminex import _categorical_order


def r_factor_levels(values, declared=None):
    return _categorical_order.r_factor_levels(values, declared)


GOLDENS = Path(__file__).resolve().parent / "r_goldens"


def _golden() -> dict[tuple[str, str], list[str]]:
    out: dict[tuple[str, str], list[str]] = defaultdict(list)
    with (GOLDENS / "factor_level_order.csv").open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            out[(row["case"], row["kind"])].append(row["value"])
    return out


R = _golden()
DISS = pd.read_csv(GOLDENS / "factor_level_order_diss.csv").to_numpy(float)
GROUP = np.array(["a", "B", "_c", "é"] * 3, dtype=object)
STRINGS = ["ABCAB", "BBCCA", "CABBA", "AACCB", "BCABC", "CCABA"]
STRINGS += ["ABBBC", "BACCA", "CBAAB", "AABBC", "BCCAA", "CABCB"]
SEQS = tx.seqdef(np.array([list(s) for s in STRINGS], dtype=object), alphabet=["A", "B", "C"])
LABELS = ["b", "a", "B", "_c", "é", "e", "A>%", "A>*", "10", "9", "faible", "élevé"]


def test_r_factor_levels_matches_r_levels() -> None:
    assert r_factor_levels(LABELS) == R[("labels", "levels")]


def test_dissassoc_group_rows_follow_r() -> None:
    rows = list(tx.dissassoc(DISS, GROUP, R=1).groups.index)
    assert rows == R[("dissassoc", "group_rows")]


def test_disstree_rule_levels_follow_r() -> None:
    tree = tx.disstree(DISS, pd.DataFrame({"g": GROUP}), R=1, min_size=0.1, max_depth=2)
    assert list(tx.disstree_get_rules(tree)) == R[("disstree", "rules")]


def test_seqrep_grp_order_follows_r() -> None:
    # R prefixes each name with a space (" _c"); the order is the contract here.
    expected = [name.strip() for name in R[("seqrep_grp", "groups")]]
    assert list(tx.seqrep_grp(SEQS, group=GROUP, diss=DISS, ret="stat")) == expected


@pytest.mark.parametrize(
    ("case", "min_group"), [("dissmergegroups_no_merge", 4), ("dissmergegroups_merged", 2)]
)
def test_dissmergegroups_numbering_matches_r(case: str, min_group: int) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # the numbering is the contract here
        result = tx.dissmergegroups(DISS, GROUP, min_group=min_group)
    assert [str(v) for v in np.asarray(result["group"]).tolist()] == R[(case, "group")]


def test_numeric_and_boolean_levels_do_not_need_icu(monkeypatch) -> None:
    def no_icu():
        raise RuntimeError("ICU unavailable")

    monkeypatch.setattr(_categorical_order, "_icu_root_collation_api", no_icu)
    assert r_factor_levels([10, 9, 2, 9]) == [2, 9, 10]
    assert r_factor_levels(np.array([2.5, 1.0, np.nan])) == [1.0, 2.5]
    assert r_factor_levels([True, False, True]) == [False, True]
    assert r_factor_levels(np.array([True, False])) == [False, True]
    with pytest.raises(RuntimeError, match="ICU unavailable"):
        r_factor_levels(["b", "a"])


def test_declared_categories_keep_their_order_and_drop_unused() -> None:
    assert r_factor_levels(["b", "a", "b"], declared=["z", "b", "a"]) == ["b", "a"]


def test_mixed_text_and_numbers_are_refused() -> None:
    with pytest.raises(TypeError, match="mix text and numbers"):
        r_factor_levels(["a", 1])
