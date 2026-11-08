"""Tied events and person ids follow R's collation (audit S2-07).

``TSE_to_STS`` and ``FCE_to_TSE`` ordered events tied in time, and
``toPersonPeriod`` ordered person ids, by Python code points, so ``C`` came
before ``b`` and ``_`` after the capitals. R orders them with ``order()`` on a
character key (``TSE_TO_STS.R:35``, ``FCE_to_TSE.R:40``,
``toPersonPeriod.R:15``), which sorts by the session collation: ICU under a
UTF-8 locale (``_x, a, b, B, C``), code points only under ``LC_COLLATE=C``.
With an explicit transition matrix the order of tied events changes the final
state, and a person-period table came back in a different row order.

Source of truth: live TraMineRextras 0.6.9 under ``LC_COLLATE=C.UTF-8``, the
collation ``traminex._categorical_order.r_character_order`` reproduces; values
in ``tests/r_goldens/audit20260922_s2_cases.txt``
(``scripts/r_goldens/make_goldens_audit20260922_s2.R``). The ICU guard mirrors
``tests/test_icu_collation_parity.py``.
"""

from __future__ import annotations

import ctypes.util
from pathlib import Path

import pandas as pd
import pytest

from traminex import FCE_to_TSE, TSE_to_STS, seqdef, toPersonPeriod

pytestmark = pytest.mark.skipif(
    ctypes.util.find_library("icui18n") is None,
    reason="R-compatible collation needs the system ICU library",
)

FIXTURE = Path(__file__).resolve().parent / "r_goldens" / "audit20260922_s2_cases.txt"


def _r(case: str, field: str) -> list[str]:
    for line in FIXTURE.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        name, key, value = line.split("\t")
        if (name, key) == (case, field):
            return value.split("|")
    raise AssertionError(f"no R record {case}/{field}")


def test_fce_to_tse_orders_tied_events_by_r_collation():
    fce = pd.DataFrame({"pid": [1, 2], "Cev": [3, 1], "bev": [3, 2], "aev": [5, 1], "_ev": [3, 1]})

    actual = FCE_to_TSE(fce, id="pid", cols=["Cev", "bev", "aev", "_ev"])

    assert [str(value) for value in actual["id"]] == _r("fce_tied_events", "id")
    assert [str(value) for value in actual["time"]] == _r("fce_tied_events", "time")
    assert actual["event"].tolist() == _r("fce_tied_events", "event")


def test_fce_to_tse_orders_character_times_by_r_collation():
    fce = pd.DataFrame({"pid": [1], "e1": ["b"], "e2": ["B"], "e3": ["a"]})

    actual = FCE_to_TSE(fce, id="pid", cols=["e1", "e2", "e3"])

    assert actual["time"].tolist() == _r("fce_character_times", "time")
    assert actual["event"].tolist() == _r("fce_character_times", "event")


def test_tse_to_sts_applies_tied_events_in_r_collation_order():
    stm = pd.DataFrame(
        [
            ["S_b", "S_bC"],
            ["S_b", "S_b_then_C"],
            ["S_bC", "S_C"],
            ["S_bC", "S_bC"],
            ["S_b_then_C", "S_b_then_C"],
        ],
        index=["None", "S_b", "S_C", "S_bC", "S_b_then_C"],
        columns=["b", "C"],
    )
    tse = pd.DataFrame({"id": [1, 1], "time": [2, 2], "event": ["C", "b"]})

    actual = TSE_to_STS(tse, stm=stm, tmax=4)

    decoded = [actual._state_for_code(int(code)) for code in actual.data[0]]
    assert decoded == _r("tse_tied_events", "cells")


def test_to_person_period_orders_ids_by_r_collation():
    x = seqdef(
        pd.DataFrame({"a": ["A", "B", "A", "B"], "b": ["B", "B", "A", "A"]}),
        id=["b", "C", "a", "_z"],
    )

    actual = toPersonPeriod(x)

    assert actual["id"].tolist() == _r("toPersonPeriod_id_order", "id")
    by_id = {person: group["state"].tolist() for person, group in actual.groupby("id")}
    assert by_id == {"b": ["A", "B"], "C": ["B", "B"], "a": ["A", "A"], "_z": ["B", "A"]}
