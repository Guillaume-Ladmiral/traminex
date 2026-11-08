"""Spellings traminex accepts beside R's, recorded as a traminex convenience (TX-AUD-20261019-061, -10-074).

Owner decision of 2026-10-22 (option B of triage R9-g). ``seqdist`` reads ``method``, ``subm`` and ``norm`` in any
letter case and takes the ``norm`` aliases ``None``, ``'length'`` and ``'max'``; ``seqMD`` reads ``method`` in any case.
R refuses every one of these spellings; ``scripts/r_goldens/make_goldens_spelling_convention.R`` records R's error and
R's value for the canonical spelling, which traminex must give. The seqdist and seqMD ledger rows and
docs/DISTANCES_NORMALIZATION.md list the spellings.
"""

from __future__ import annotations

import csv
import warnings
from pathlib import Path

import numpy as np
import pytest

import traminex as tx

REPO = Path(__file__).resolve().parents[1]
with open(
    REPO / "tests" / "r_goldens" / "spelling_convention.csv", encoding="utf-8", newline=""
) as handle:
    GOLDEN = {row["case"]: row for row in csv.DictReader(handle)}
X = tx.seqdef(
    np.array([list("AABC"), list("ABBC"), list("CCBA"), list("BBCA")], dtype=object),
    alphabet=["A", "B", "C"],
)
CASES = {
    "method_om": dict(method="om", subm="CONSTANT"),
    "method_ham": dict(method="ham"),
    "method_lcs": dict(method="lcs"),
    "sm_trate": dict(method="OM", subm="trate"),
    "sm_constant": dict(method="OM", subm="constant"),
    "norm_upper_case": dict(method="OM", subm="CONSTANT", norm="MAXLENGTH"),
    "norm_length": dict(method="OM", subm="CONSTANT", norm="length"),
    "norm_max": dict(method="OM", subm="CONSTANT", norm="max"),
    "norm_null": dict(method="OM", subm="CONSTANT", norm=None),
    "norm_yujianbo": dict(method="OM", subm="CONSTANT", norm="yujianbo"),
    # Follow-up review m-7.
    "method_omspell": dict(method="omspell", subm="CONSTANT"),
    "method_twed": dict(method="twed", subm="CONSTANT", nu=0.5),
    "method_mixed_case": dict(method="Om", subm="CONSTANT"),
    "sm_indels": dict(method="OM", subm="indels"),
    "sm_indelslog": dict(method="OM", subm="indelslog"),
    "norm_gmean_upper": dict(method="OM", subm="CONSTANT", norm="GMEAN"),
    "norm_auto_upper": dict(method="OM", subm="CONSTANT", norm="AUTO"),
    "norm_none_upper": dict(method="OM", subm="CONSTANT", norm="NONE"),
}
# seqMD passes norm to seqdist, so it reads the same spellings (follow-up review I-1).
SEQMD_NORM_CASES = {
    "seqmd_norm_upper_case": "GMEAN",
    "seqmd_norm_length": "length",
    "seqmd_norm_max": "max",
    "seqmd_norm_null": None,
}


def _canonical(case: str) -> np.ndarray:
    values = np.array([float(v) for v in GOLDEN[case]["canonical"].split()])
    size = int(round(values.size**0.5))
    return values.reshape(size, size, order="F")


def test_every_golden_case_is_checked() -> None:
    assert sorted(GOLDEN) == sorted([*CASES, *SEQMD_NORM_CASES, "seqmd_method_om"])


@pytest.mark.parametrize("case", sorted(CASES))
def test_seqdist_gives_r_s_canonical_value_for_a_spelling_r_refuses(case: str) -> None:
    assert GOLDEN[case]["ok"] == "FALSE", GOLDEN[case]["error"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        actual = np.asarray(tx.seqdist(X, **CASES[case]), dtype=float)
    np.testing.assert_allclose(actual, _canonical(case), rtol=1e-12, atol=1e-14)


def test_seqmd_gives_r_s_canonical_value_for_a_lowercase_method() -> None:
    assert GOLDEN["seqmd_method_om"]["ok"] == "FALSE", GOLDEN["seqmd_method_om"]["error"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        actual = np.asarray(
            tx.seqMD([X, X], method="om", sm=["CONSTANT", "CONSTANT"], what="diss"), dtype=float
        )
    np.testing.assert_allclose(actual, _canonical("seqmd_method_om"), rtol=1e-12, atol=1e-14)


@pytest.mark.parametrize("case", sorted(SEQMD_NORM_CASES))
def test_seqmd_gives_r_s_canonical_value_for_a_norm_spelling_r_refuses(case: str) -> None:
    assert GOLDEN[case]["ok"] == "FALSE", GOLDEN[case]["error"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        actual = np.asarray(
            tx.seqMD(
                [X, X],
                method="OM",
                sm=["CONSTANT", "CONSTANT"],
                norm=SEQMD_NORM_CASES[case],
                what="diss",
            ),
            dtype=float,
        )
    np.testing.assert_allclose(actual, _canonical(case), rtol=1e-12, atol=1e-14)


def test_the_ledger_and_the_documents_record_the_spellings() -> None:
    with open(REPO / "docs" / "API_PARITY_DECISIONS.csv", encoding="utf-8", newline="") as handle:
        rows = {row[1]: " ".join(row) for row in csv.reader(handle) if len(row) > 1}
    for name in ("seqdist", "seqMD"):
        assert (
            "traminex convenience (owner decision 2026-10-22, TX-AUD-20261019-061" in rows[name]
        ), name
    text = " ".join(
        (REPO / "docs" / "DISTANCES_NORMALIZATION.md").read_text(encoding="utf-8").split()
    )
    assert "`method` and `subm` are read in any letter case" in text
    assert "seqMD reads method and norm in any letter case" in rows["seqMD"]
    assert "`seqMD`'s `method` and `norm`" in text
