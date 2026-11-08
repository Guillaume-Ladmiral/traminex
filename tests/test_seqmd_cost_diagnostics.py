"""seqMD emits TraMineR 2.2-14's cost diagnostics for every channel, in R's order.

``seqMD.R:266-267`` runs ``checkcost`` on each channel's costs before weighting
them, for ``what="cost"`` and ``"diss"``: on an explicit matrix or array, and
equally on the costs seqMD generates from a method string, which it builds with
``seqcost(..., check.costs = FALSE)`` "because costs are checked below"
(``seqMD.R:241-243``). The call is seqMD's own and is in 2.2-13 too; it is not
the check 2.2-14 added inside ``seqcost`` (docs/API_PARITY.md). With
``what="diss"`` the combined costs then go through ``seqdist``, which runs its
own check (none for HAM). traminex's seqMD emitted none of the channel
diagnostics, and ``tests/test_r_goldens_seqMD.py`` never asserted the ones
``seqMD_r2_cases.csv`` records (2026-09-22 audit S10-11).

Two R records are compared. ``tests/r_goldens/seqMD_checkcost_diagnostics.txt``
(``scripts/r_goldens/make_goldens_audit20260922_s10b.R``) keeps R's exact text
and channel: ``message()`` diagnostics arrive as ``CostDiagnosticMessage``
warnings, ``warning()`` as ``UserWarning``, and the newline ``message()``
appends is removed. ``tests/r_goldens/seqMD_r2_cases.csv``
(``scripts/r_goldens/make_goldens_seqMD.R``) keeps its cases' diagnostics with
commas and newlines turned into spaces (that generator's ``r2_add``), so
traminex's text is put through the same transformation before comparing.
"""

from __future__ import annotations

import csv
import re
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pytest

from traminex import seqcost, seqdef, seqMD
from traminex.distances.python_kernels import CostDiagnosticMessage

GOLDENS = Path(__file__).parent / "r_goldens"


def _matrix(*values: float) -> np.ndarray:
    return np.array(values, dtype=float).reshape(3, 3)


SYM = _matrix(0, 2, 2, 2, 0, 2, 2, 2, 0)
GT = _matrix(0, 5, 5, 5, 0, 5, 5, 5, 0)
TRI = _matrix(0, 10, 2, 10, 0, 2, 2, 2, 0)
ASY = _matrix(0, 3, 2, 1, 0, 2, 2, 2, 0)
ZAB = _matrix(0, 0, 2, 0, 0, 2, 2, 2, 0)
ZALL = np.zeros((3, 3))


def _a():
    return seqdef([list("ABAB"), list("BABA"), list("AABB")], alphabet=["A", "B", "U"])


def _x():
    return seqdef([list("XYXY"), list("YXYX"), list("XXYY")], alphabet=["X", "Y", "V"])


def _c():
    return seqdef([list("PQQP"), list("QPPQ"), list("PPQQ")], alphabet=["P", "Q", "R"])


def _m():
    return seqdef([["A", None, "B", "A"], list("BABA"), list("ABAB")], alphabet=["A", "B", "U"])


def _t():
    return seqdef(
        [list("ABAB"), list("BABA"), list("CCAB"), list("AACC")], alphabet=["A", "B", "C"]
    )


def _s():
    return seqdef(
        [list("XYXY"), list("YXYX"), list("XXZZ"), list("ZYXX")], alphabet=["X", "Y", "Z"]
    )


def _missing_costs() -> np.ndarray:
    costs = np.full((4, 4), 2.0)
    np.fill_diagonal(costs, 0.0)
    costs[2, 3] = costs[3, 2] = 0.0  # (U, *)
    return costs


def _arrays() -> tuple[np.ndarray, np.ndarray]:
    """R's A_zero and A_asym: states x states x time, as seqMD takes them."""
    base = np.full((3, 3, 4), 2.0)
    for time in range(4):
        np.fill_diagonal(base[:, :, time], 0.0)
    zero = base.copy()
    zero[0, 1, 1] = zero[1, 0, 1] = 0.0
    asymmetric = base.copy()
    asymmetric[0, 1, 0] = asymmetric[0, 1, 2] = 3.0
    return zero, asymmetric


def _dhd_sm(k: int, width: int, scale: float) -> np.ndarray:
    out = np.zeros((k, k, width))
    for time in range(width):
        for left in range(k):
            for right in range(k):
                out[left, right, time] = scale * abs(left - right) * (time + 1)
    return out


def _time_varying(channel, method):
    has_missing = bool(channel.nr_code is not None and np.any(channel.data == channel.nr_code))
    # The generator builds these costs with suppressMessages(seqcost(...)) before the
    # recorded call (make_goldens_seqMD.R:241-246); seqcost now gives R's messages for
    # the matrices it generates (TX-AUD-20261011-080), so silence them the same way.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", CostDiagnosticMessage)
        costs = seqcost(channel, method, with_missing=has_missing, time_varying=True)["sm"]
    return np.asarray(costs, dtype=float)


EXACT_CASES = {
    "cost_explicit_zero_then_triangle": lambda: seqMD(
        [_a(), _x()], what="cost", sm=[ZAB, TRI], indel=[1, 6]
    ),
    "cost_explicit_twice_indel_then_asymmetric": lambda: seqMD(
        [_a(), _x()], what="cost", sm=[GT, ASY], indel=[1, 2]
    ),
    "cost_explicit_indel_vectors": lambda: seqMD(
        [_a(), _x()], what="cost", sm=[SYM, SYM], indel=[[0.1, 5, 5], [1, 1, 1]]
    ),
    "cost_explicit_all_zero_channel": lambda: seqMD(
        [_a(), _x()], what="cost", sm=[ZALL, SYM], indel=[1, 1]
    ),
    "cost_explicit_auto_indel": lambda: seqMD([_a(), _x()], what="cost", sm=[GT, TRI]),
    "cost_explicit_missing_state_label": lambda: seqMD(
        [_m(), _x()],
        what="cost",
        sm=[_missing_costs(), SYM],
        indel=[1, 1],
        with_missing=[True, False],
    ),
    "cost_generated_trate": lambda: seqMD([_t(), _s()], what="cost", sm="TRATE"),
    "cost_method_dhd_generated_trate": lambda: seqMD(
        [_m(), _x()], method="DHD", what="cost", sm="TRATE", with_missing=[True, False]
    ),
    "cost_explicit_arrays": lambda: seqMD(
        [_a(), _x()], method="DHD", what="cost", sm=list(_arrays()), indel=[1, 1]
    ),
    "diss_om_explicit_zero": lambda: seqMD(
        [_a(), _x()], method="OM", what="diss", sm=[ZAB, SYM], indel=[1, 1]
    ),
    "diss_ham_explicit_triangle": lambda: seqMD(
        [_a(), _x()], method="HAM", what="diss", sm=[TRI, SYM]
    ),
    "diss_dhd_generated_trate": lambda: seqMD(
        [_m(), _x()], method="DHD", what="diss", sm="TRATE", with_missing=[True, False]
    ),
    "mdseq_checks_no_costs": lambda: seqMD([_a(), _x()], what="MDseq", sm=[ZAB, TRI]),
}

#: Cost and distance cases of make_goldens_seqMD.R, rebuilt as its tests do.
R2_CASES = {
    "cost_trate_declared_2": lambda: seqMD(
        [_a(), _x()], what="cost", sm="TRATE", indel=[1, 2], cweight=[2, 1]
    ),
    "cost_indels_vector_mean": lambda: seqMD(
        [_a(), _x()],
        what="cost",
        sm=[np.asarray(seqcost(ch, "INDELS")["sm"]) for ch in (_a(), _x())],
        indel=[[1, 2, 3], [4, 5, 6]],
        cweight=[2, 1],
        link="mean",
    ),
    "cost_indelslog_auto": lambda: seqMD([_a(), _x()], what="cost", sm="INDELSLOG"),
    "cost_trate_declared_3": lambda: seqMD(
        [_a(), _x(), _c()], what="cost", sm="TRATE", indel=1.5, cweight=[1, 2, 3]
    ),
    "dhd_trate_cost_missing": lambda: seqMD(
        [_m(), _x()],
        method="DHD",
        what="cost",
        sm=[_time_varying(_m(), "TRATE"), _time_varying(_x(), "TRATE")],
        indel=[1, 2],
        with_missing=[True, False],
        cweight=[1, 2],
    ),
    "dhd_indels_cost_3": lambda: seqMD(
        [_a(), _x(), _c()],
        method="DHD",
        what="cost",
        sm=[_time_varying(ch, "INDELS") for ch in (_a(), _x(), _c())],
        indel=[1, 1, 1],
        cweight=[1, 0.5, 2],
    ),
    "dhd_indelslog_cost_mean": lambda: seqMD(
        [_a(), _x()],
        method="DHD",
        what="cost",
        sm=[_time_varying(ch, "INDELSLOG") for ch in (_a(), _x())],
        indel=[1, 1],
        cweight=[2, 1],
        link="mean",
    ),
    "dhd_explicit_cost_missing": lambda: seqMD(
        [_m(), _x()],
        method="DHD",
        what="cost",
        sm=[_dhd_sm(4, 4, 0.5), _dhd_sm(3, 4, 0.25)],
        indel=[1, 2],
        with_missing=[True, False],
        cweight=[2, 1],
    ),
    "dhd_trate_diss_full": lambda: seqMD(
        [_m(), _x()], method="DHD", what="diss", sm="TRATE", with_missing=[True, False]
    ),
    "dhd_indels_diss_full_3": lambda: seqMD(
        [_a(), _x(), _c()], method="DHD", what="diss", sm="INDELS", indel=1, cweight=[1, 0.5, 2]
    ),
    "dhd_indelslog_diss_condensed": lambda: seqMD(
        [_a(), _x()],
        method="DHD",
        what="diss",
        sm="INDELSLOG",
        indel=1,
        cweight=[2, 1],
        link="mean",
        full_matrix=False,
    ),
    "dhd_explicit_diss_condensed_3": lambda: seqMD(
        [_a(), _x(), _c()],
        method="DHD",
        what="diss",
        sm=[_dhd_sm(3, 4, 0.5), _dhd_sm(3, 4, 0.25), _dhd_sm(3, 4, 0.75)],
        indel=1,
        cweight=[1, 2, 1],
        full_matrix=False,
    ),
}


def _exact_record() -> dict[str, list[tuple[str, str]]]:
    recorded: dict[str, list[tuple[str, str]]] = defaultdict(list)
    path = GOLDENS / "seqMD_checkcost_diagnostics.txt"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = csv.reader((line for line in handle if not line.startswith("#")), delimiter="\t")
        for case, _order, channel, text in rows:
            if channel == "none":
                recorded[case] = []
                continue
            decoded = text.replace("\\n", "\n")
            if channel == "message":
                assert decoded.endswith("\n"), (case, decoded)
                decoded = decoded[:-1]
            recorded[case].append((channel, decoded))
    return dict(recorded)


def _r2_record() -> dict[str, list[str]]:
    recorded: dict[str, list[str]] = defaultdict(list)
    with (GOLDENS / "seqMD_r2_cases.csv").open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["case"] in R2_CASES and row["section"] == "warning":
                recorded[row["case"]].append((int(row["col"]), row["text"]))
    return {case: [text for _, text in sorted(recorded[case])] for case in R2_CASES}


EXACT = _exact_record()
R2 = _r2_record()


def _diagnostics(call) -> list[tuple[str, str]]:
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        call()
    out = []
    for item in caught:
        if issubclass(item.category, CostDiagnosticMessage):
            out.append(("message", str(item.message)))
        elif item.category is UserWarning:
            out.append(("warning", str(item.message)))
    return out


def test_every_exactly_recorded_case_is_exercised():
    assert set(EXACT_CASES) == set(EXACT)


@pytest.mark.parametrize("case", sorted(EXACT_CASES))
def test_seqmd_diagnostics_match_traminer_2_2_14(case):
    """``diss_ham_explicit_triangle`` included: R's seqdist never checks HAM,
    and with two combined states traminex's documented HAM triangle
    diagnostic (docs/API_PARITY.md) cannot fire either."""
    assert _diagnostics(EXACT_CASES[case]) == EXACT[case]


def _as_r2(channel: str, text: str) -> str:
    """make_goldens_seqMD.R's r2_add: message()'s newline kept, [\\r\\n,] -> space."""
    return re.sub(r"[\r\n,]", " ", text + ("\n" if channel == "message" else ""))


@pytest.mark.parametrize("case", sorted(R2_CASES))
def test_seqmd_diagnostics_match_the_r2_golden(case):
    actual = [_as_r2(channel, text) for channel, text in _diagnostics(R2_CASES[case])]
    assert actual == R2[case]


def test_the_r2_golden_records_diagnostics_for_three_of_its_cases():
    """Guard the comparison above against an empty record."""
    assert sorted(case for case, texts in R2.items() if texts) == [
        "cost_indels_vector_mean",
        "dhd_trate_cost_missing",
        "dhd_trate_diss_full",
    ]
