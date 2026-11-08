"""Indicators skip an internal void, where R reads the first seqlength cells (TX-AUD-20261010-055).

seqrecode(otherwise=void) can leave void cells inside a sequence (here
p-n-%-%-p-n). A void cell is not part of a sequence, so the sequence is
p-n-p-n. R's indicators read the first seqlength cells of each row (seqdss.R:12,37;
seqdur.R:11,29-31): they read p-n-%-% and lose the final p-n. traminex skips
void wherever it is, a recorded correction (ledger rows of the functions below).
The expected values are R's own on the void-free sequences, where R is right.
The second row's void is at the end, where R and traminex agree. traminex
already skipped internal void, so the tests characterize it. R's seqformat
writes an internal void as a missing spell (it sets nr and void alike to
missing); traminex's seqformat does too, while the string writers to_SPS and
to_DSS leave void out, as R's forms of the twin do. Fixture:
``scripts/r_goldens/make_goldens_internal_void_indicators.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx
from traminex.seqformat import to_DSS, to_SPS

# seqindic repeats R's message that no prec.args are given; it is not the subject here.
pytestmark = pytest.mark.filterwarnings("ignore:No prec.args list:UserWarning")

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "internal_void_indicators.csv",
    dtype=str,
    keep_default_na=False,
)
BASE = tx.seqdef(np.array([list("ABCCAB"), list("ABAACC")], dtype=object), alphabet=["A", "B", "C"])
RECODED = tx.seqrecode(BASE, {"p": ["A"], "n": ["B"]}, otherwise="%")
NUMERIC = [
    "seqtransn",
    "seqST",
    "seqici",
    "seqivolatility",
    "seqidegrad",
    "seqprecarity",
    "seqinsecurity",
    "seqsubsn",
]
INDIC = [
    "dlgth",
    "recu",
    "trans",
    "transp",
    "nsubs",
    "volat",
    "cplx",
    "turb",
    "turbn",
    "prec",
    "insec",
]


def _r(obj: str, function: str) -> list[str]:
    rows = GOLDEN[(GOLDEN["object"] == obj) & (GOLDEN["function_name"] == function)]
    return rows.sort_values("position", key=lambda c: c.astype(int))["value"].tolist()


def _r_numbers(obj: str, function: str) -> np.ndarray:
    return np.array([np.nan if v == "NA" else float(v) for v in _r(obj, function)])


def _ours(function: str) -> np.ndarray:
    return np.asarray(getattr(tx, function)(RECODED), dtype=float).ravel()


def test_the_recoded_rows_hold_an_internal_and_a_trailing_void() -> None:
    assert RECODED.to_strings() == ["p-n-%-%-p-n", "p-n-p-p-%-%"]
    assert _r("recoded", "seqconc") == _r("twin", "seqconc") == ["p-n-p-n", "p-n-p-p"]


def test_seqdss_is_r_seqdss_of_the_void_free_sequences() -> None:
    assert [text.rstrip("-%") for text in tx.seqdss(RECODED).to_strings()] == _r("twin", "seqdss")
    assert _r("recoded", "seqdss") == ["p-n", "p-n-p"]  # R drops the final p-n of row 1


def test_seqdur_is_r_seqdur_of_the_void_free_sequences() -> None:
    twin = _r_numbers("twin", "seqdur").reshape(2, -1, order="F")
    ours = np.asarray(tx.seqdur(RECODED), dtype=float)
    np.testing.assert_array_equal(ours[:, : twin.shape[1]], twin)
    recoded = _r_numbers("recoded", "seqdur").reshape(2, -1, order="F")
    assert not np.array_equal(recoded[0], twin[0, : recoded.shape[1]], equal_nan=True)


@pytest.mark.parametrize("function", NUMERIC)
def test_indicator_is_r_value_of_the_void_free_sequences(function: str) -> None:
    np.testing.assert_allclose(_ours(function), _r_numbers("twin", function), rtol=0, atol=1e-12)


@pytest.mark.parametrize("function", NUMERIC)
def test_r_differs_only_on_the_row_with_the_internal_void(function: str) -> None:
    recoded, twin = _r_numbers("recoded", function), _r_numbers("twin", function)
    assert not np.isclose(recoded[0], twin[0], rtol=0, atol=1e-12)
    assert np.isclose(recoded[1], twin[1], rtol=0, atol=1e-12)


@pytest.mark.parametrize(
    "column",
    [
        "Dlgth",
        "Recu",
        "Trans",
        "Transp",
        "Nsubs",
        "Volat",
        "Cplx",
        "Turb",
        "Turbn",
        "Prec",
        "Insec",
    ],
)
def test_seqindic_column_is_r_value_of_the_void_free_sequences(column: str) -> None:
    ours = np.asarray(tx.seqindic(RECODED, indic=INDIC)[column], dtype=float)
    twin = _r_numbers("twin", f"seqindic:{column}")
    np.testing.assert_allclose(ours, twin, rtol=0, atol=1e-12)
    assert not np.isclose(
        _r_numbers("recoded", f"seqindic:{column}")[0], twin[0], rtol=0, atol=1e-12
    )


def test_seqlogp_agrees_with_r_here_and_is_not_claimed() -> None:
    np.testing.assert_allclose(
        _ours("seqlogp"), _r_numbers("recoded", "seqlogp"), rtol=0, atol=1e-12
    )


def test_seqformat_writes_an_internal_void_as_r_does() -> None:
    # R's seqformat says so: "'missing' set as c('*','%')".
    with pytest.warns(UserWarning, match="'missing' set as"):
        sps = tx.seqformat(RECODED, from_="STS", to="SPS", compress=True)
    assert sps == _r("recoded", "seqformat_SPS")
    with pytest.warns(UserWarning, match="'missing' set as"):
        dss = tx.seqformat(RECODED, from_="STS", to="DSS", compress=True)
    assert dss.to_strings() == _r("recoded", "seqformat_DSS")


def test_the_string_writers_leave_void_out_as_r_forms_of_the_twin() -> None:
    sps = [
        "-".join(f"({token.replace('/', ',')})" for token in text.split("-"))
        for text in to_SPS(RECODED)
    ]
    assert sps == _r("twin", "seqformat_SPS")
    assert to_DSS(RECODED) == _r("twin", "seqformat_DSS")
