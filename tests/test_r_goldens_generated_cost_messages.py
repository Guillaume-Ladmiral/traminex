"""checkcost's messages for generated cost matrices, as R prints them (TX-AUD-20261011-080).

TraMineR 2.2-14 runs checkcost on every matrix seqcost generates (seqcost.R:286),
and seqdist resolves a cost-method string through seqcost once (seqdist.R:426-431).
traminex printed these messages only for a matrix passed explicitly. R's
``message()`` diagnostics are CostDiagnosticMessage warnings and R's ``warning()``
diagnostics UserWarning, with R's text (the newline ``message()`` appends removed).
A reference is computed by a second seqdist call, which printed the messages
again (batch-6 review I-1). Fixture
``scripts/r_goldens/make_goldens_generated_cost_messages.R``.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx
from traminex.distances.python_kernels import CostDiagnosticMessage
from traminex.sequences import STSList

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "generated_cost_messages.csv",
    dtype=str,
    keep_default_na=False,
)
TRI = tx.seqdef(
    np.array([list("ABABABAB"), list("BCBCBCBC")], dtype=object), alphabet=["A", "B", "C"]
)
TRI3 = tx.seqdef(
    np.array([list("ABABABAB"), list("BCBCBCBC"), list("ABBCABBC")], dtype=object),
    alphabet=["A", "B", "C"],
)
# Breaks the triangle inequality: d(A, C) = 3 > d(A, B) + d(B, C) = 2.
TRI_SM = np.array([[0, 1, 3], [1, 0, 1], [3, 1, 0]], dtype=float)
ALT = tx.seqdef(
    np.array([list("ABAB"), list("BABA"), list("ABAB")], dtype=object), alphabet=["A", "B"]
)
ONE_FULL = tx.seqdef(np.array([["A", "A"], ["A", "A"]], dtype=object), alphabet=["A"])
ONE = tx.seqdef(
    np.array([["A"] * 4, ["A", "A", None, None], ["A", "A", "A", None]], dtype=object),
    alphabet=["A"],
)
CASES = {
    "seqcost_trate": lambda: tx.seqcost(TRI, "TRATE"),
    "seqsubm_trate": lambda: tx.seqsubm(TRI, method="TRATE"),
    "seqdist_om_trate_indel1": lambda: tx.seqdist(TRI, method="OM", subm="TRATE", indel=1),
    "seqdist_om_trate_auto": lambda: tx.seqdist(TRI, method="OM", subm="TRATE"),
    "seqcost_constant_cval0": lambda: tx.seqcost(TRI, "CONSTANT", cval=0),
    "seqcost_constant": lambda: tx.seqcost(TRI, "CONSTANT"),
    "seqcost_constant_cval5": lambda: tx.seqcost(TRI, "CONSTANT", cval=5),
    "seqdist_om_constant": lambda: tx.seqdist(TRI, method="OM", subm="CONSTANT"),
    "seqcost_indels": lambda: tx.seqcost(TRI, "INDELS"),
    "seqcost_indelslog": lambda: tx.seqcost(TRI, "INDELSLOG"),
    "seqdist_om_indels": lambda: tx.seqdist(TRI, method="OM", subm="INDELS"),
    "seqdist_om_indelslog": lambda: tx.seqdist(TRI, method="OM", subm="INDELSLOG"),
    "seqcost_indels_one_state": lambda: tx.seqcost(ONE, "INDELS"),
    "seqcost_indelslog_one_state": lambda: tx.seqcost(ONE, "INDELSLOG"),
    "seqdist_om_indels_one_state": lambda: tx.seqdist(ONE, method="OM", subm="INDELS"),
    "seqdist_om_indelslog_one_state": lambda: tx.seqdist(ONE, method="OM", subm="INDELSLOG"),
    "seqcost_future": lambda: tx.seqcost(TRI, "FUTURE"),
    "seqcost_trate_time_varying": lambda: tx.seqcost(TRI, "TRATE", time_varying=True),
    "seqdist_dhd": lambda: tx.seqdist(TRI, method="DHD"),
    "seqdist_ham_trate": lambda: tx.seqdist(TRI, method="HAM", subm="TRATE"),
    "seqdist_lcs_one_state": lambda: tx.seqdist(ONE, method="LCS"),
    "seqdist_ham_one_state": lambda: tx.seqdist(ONE_FULL, method="HAM"),
    "seqdist_ham_default": lambda: tx.seqdist(TRI, method="HAM"),
    "seqmd_trate": lambda: tx.seqMD([TRI, TRI], method="OM", sm=["TRATE", "TRATE"], what="diss"),
    "seqcost_trate_time_varying_alternating": lambda: tx.seqcost(ALT, "TRATE", time_varying=True),
    "seqdist_dhd_alternating": lambda: tx.seqdist(ALT, method="DHD"),
    "seqdist_dhd_trate_alternating": lambda: tx.seqdist(ALT, method="DHD", subm="TRATE"),
    "seqdist_dhd_indels_alternating": lambda: tx.seqdist(ALT, method="DHD", subm="INDELS"),
    "seqdist_om_trate_refseq_1": lambda: tx.seqdist(TRI3, method="OM", subm="TRATE", refseq=1),
    "seqdist_om_trate_refseq_0": lambda: tx.seqdist(TRI3, method="OM", subm="TRATE", refseq=0),
    "seqdist_om_trate_sets": lambda: tx.seqdist(
        TRI3, method="OM", subm="TRATE", refseq=([1, 2], [3])
    ),
    "seqdist_om_trate_external": lambda: tx.seqdist(
        TRI3, method="OM", subm="TRATE", refseq=_row(TRI3, 1)
    ),
    "seqdist_om_trate_norm_refseq_1": lambda: tx.seqdist(
        TRI3, method="OM", subm="TRATE", norm="auto", refseq=1
    ),
    "seqdist_om_constant_refseq_1": lambda: tx.seqdist(
        TRI3, method="OM", subm="CONSTANT", refseq=1
    ),
    "seqdist_ham_trate_refseq_1": lambda: tx.seqdist(TRI3, method="HAM", subm="TRATE", refseq=1),
    "seqdist_omspell_trate_refseq_1": lambda: tx.seqdist(
        TRI3, method="OMspell", subm="TRATE", refseq=1
    ),
    "seqdist_dhd_refseq_1": lambda: tx.seqdist(TRI3, method="DHD", refseq=1),
    "seqdist_om_matrix_refseq_1": lambda: tx.seqdist(
        TRI3, method="OM", subm=TRI_SM, indel=2, refseq=1
    ),
    "seqdist_om_matrix_sets": lambda: tx.seqdist(
        TRI3, method="OM", subm=TRI_SM, indel=2, refseq=([1, 2], [3])
    ),
    "seqdist_lcs_one_state_refseq_1": lambda: tx.seqdist(ONE, method="LCS", refseq=1),
    "seqdist_ham_one_state_refseq_1": lambda: tx.seqdist(ONE_FULL, method="HAM", refseq=1),
}


def _row(seqs: STSList, index: int) -> STSList:
    return STSList(
        data=seqs.data[index : index + 1].copy(),
        alphabet=list(seqs.alphabet),
        nr_code=seqs.nr_code,
        void_code=seqs.void_code,
    )


def _r_messages(case: str) -> list[tuple[str, str]]:
    rows = GOLDEN[GOLDEN["case"] == case].sort_values("index", key=lambda c: c.astype(int))
    if rows["n"].iloc[0] == "0":
        return []
    return [(row["kind"], row["text"]) for _, row in rows.iterrows()]


def _traminex_messages(call) -> list[tuple[str, str]]:
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        call()
    found = []
    for warning in caught:
        text = str(warning.message)
        if "[!" not in text:
            continue
        kind = "message" if issubclass(warning.category, CostDiagnosticMessage) else "warning"
        found.append((kind, text))
    return found


def test_every_r_case_is_rebuilt() -> None:
    assert set(CASES) == set(GOLDEN["case"])


@pytest.mark.parametrize("case", sorted(CASES))
def test_messages_match_r(case: str) -> None:
    assert _traminex_messages(CASES[case]) == _r_messages(case)
