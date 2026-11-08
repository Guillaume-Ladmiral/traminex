"""seqdist's empty-sequence and with.missing messages, as R gives them (TX-AUD-20261011-042).

R warns "Empty sequence(s) <indices> may cause inconsistent results!" for a
sequence with no element to compare (seqdist.R:101-109), stops for OMloc, and
says so when it resets with.missing to FALSE because nothing is missing
(:114-117). Both are R ``message()`` calls with R's " [!!]" prefix
(TraMineR-logging_helpers.R:20-25, 89-91); traminex printed neither. They are
UserWarning here, with R's text. A reference is computed by a second seqdist
call, which printed them again (batch-6 review I-1). Fixture
``scripts/r_goldens/make_goldens_seqdist_messages.R``.
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
    Path(__file__).resolve().parent / "r_goldens" / "seqdist_messages.csv",
    dtype=str,
    keep_default_na=False,
)
EMP = tx.seqdef(np.array([["A", "B"], ["B", "B"], [None, None]], dtype=object), alphabet=["A", "B"])
EMP2 = tx.seqdef(
    np.array([[None, None], ["B", "B"], [None, None], ["A", "A"]], dtype=object),
    alphabet=["A", "B"],
)
S3 = tx.seqdef(
    np.array([["A", "B", "A"], ["B", "B", "A"], ["A", "A", "B"]], dtype=object), alphabet=["A", "B"]
)


def _row(seqs: STSList, index: int) -> STSList:
    return STSList(
        data=seqs.data[index : index + 1].copy(),
        alphabet=list(seqs.alphabet),
        nr_code=seqs.nr_code,
        void_code=seqs.void_code,
    )


CASES = {
    "empty_lcs": lambda: tx.seqdist(EMP, method="LCS"),
    "empty_om": lambda: tx.seqdist(EMP, method="OM", subm="CONSTANT"),
    "empty_ham": lambda: tx.seqdist(EMP, method="HAM"),
    "empty_omloc": lambda: tx.seqdist(EMP, method="OMloc", subm="CONSTANT"),
    "two_empty_lcs": lambda: tx.seqdist(EMP2, method="LCS"),
    "with_missing_reset_lcs": lambda: tx.seqdist(S3, method="LCS", with_missing=True),
    "with_missing_reset_om": lambda: tx.seqdist(
        S3, method="OM", subm="CONSTANT", with_missing=True
    ),
    "empty_lcs_refseq_1": lambda: tx.seqdist(EMP, method="LCS", refseq=1),
    "empty_lcs_refseq_0": lambda: tx.seqdist(EMP, method="LCS", refseq=0),
    "empty_lcs_sets": lambda: tx.seqdist(EMP, method="LCS", refseq=([1, 2], [3])),
    "empty_lcs_external": lambda: tx.seqdist(EMP, method="LCS", refseq=_row(EMP, 0)),
    "with_missing_reset_lcs_refseq_1": lambda: tx.seqdist(
        S3, method="LCS", with_missing=True, refseq=1
    ),
    "with_missing_reset_lcs_sets": lambda: tx.seqdist(
        S3, method="LCS", with_missing=True, refseq=([1, 2], [3])
    ),
}


def _r(case: str) -> tuple[bool, list[str]]:
    rows = GOLDEN[GOLDEN["case"] == case]
    texts = [text for text in rows["text"] if text]
    return rows["ok"].iloc[0] == "TRUE", texts


def test_every_r_case_is_rebuilt() -> None:
    assert set(CASES) == set(GOLDEN["case"])


@pytest.mark.parametrize("case", sorted(CASES))
def test_messages_and_outcome_match_r(case: str) -> None:
    ok, expected = _r(case)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        if ok:
            CASES[case]()
        else:
            with pytest.raises(ValueError):
                CASES[case]()
    found = [
        str(w.message)
        for w in caught
        if "[!" in str(w.message) and not issubclass(w.category, CostDiagnosticMessage)
    ]
    assert found == expected
    assert all(issubclass(w.category, UserWarning) for w in caught if "[!" in str(w.message))


def test_omloc_names_the_empty_sequences() -> None:
    with pytest.raises(ValueError, match="OMloc does not support empty sequences 3"):
        tx.seqdist(EMP, method="OMloc", subm="CONSTANT")


@pytest.mark.parametrize("refseq", [None, 1, 0, ([1, 2], [3])], ids=["full", "index", "0", "sets"])
def test_the_ignored_step_is_reported_once_per_call(refseq) -> None:
    # R reports it once (seqdist.R:256-258); the reference routes' second call
    # reported it again (batch-6 review I-1).
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        tx.seqdist(S3, method="CHI2", breaks=[(1, 2)], step=2, refseq=refseq)
    assert sum("'step' has been ignored" in str(w.message) for w in caught) == 1
