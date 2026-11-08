"""Degenerate OMstran and CHI2/EUCLID inputs: R's failures recorded, traminex's values asserted
(TX-AUD-20261019-080; ``scripts/r_goldens/make_goldens_degenerate_omstran.R``).

R's OMstran drops matrix columns without ``drop = FALSE`` (seqdist-OMstran.R:55-58): one sequence gives a 4 x 4
result or "incorrect number of dimensions", and length-1 data without add.column or length-2 data with previous
leave zero columns and stop. By OMstran's definition one sequence is at distance 0 from itself, and data that
recode to empty transition sequences are at distance 0 from each other. ``transindel="prob"`` needs transition
rates, which seqtrate refuses on length-1 data (R stops there too). R's CHI2 does not cap its first overlapping
episode at the sequence length (seqdist-CHI2.R:25-27) while it caps the others (:33). traminex caps it too, its own
correction: length-1 data with step 2 and overlap then have three identical episodes, so the distance is sqrt(3)
times the one-episode distance, and the normalized one equals it. Treating a step wider than the data like R's
step == ncol shortcut (:15-16, one episode) would be another correction (batch-9 review M-10).
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "degenerate_omstran.csv",
    dtype=str,
    keep_default_na=False,
).set_index("case")
ONE = tx.seqdef(np.array([list("AABBC")], dtype=object), alphabet=["A", "B", "C"])
LEN1 = tx.seqdef(
    np.array([["A"], ["B"], ["A"], ["C"], ["B"]], dtype=object), alphabet=["A", "B", "C"]
)
TWO = tx.seqdef(np.array([list("AB"), list("BA"), list("AA")], dtype=object), alphabet=["A", "B"])


def _quiet(*args, **kwargs) -> np.ndarray:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return np.asarray(tx.seqdist(*args, **kwargs), dtype=float)


def _omstran(x, transindel: str, previous: bool, add_column: bool) -> np.ndarray:
    return _quiet(
        x,
        method="OMstran",
        subm="CONSTANT",
        otto=0.3,
        transindel=transindel,
        previous=previous,
        add_column=add_column,
    )


def _r_values(case: str) -> np.ndarray:
    row = GOLDEN.loc[case]
    n = int(row["dim"].split("x")[0])
    return np.array([float(v) for v in row["value"].split()]).reshape(n, n, order="F")


@pytest.mark.parametrize(
    "case, data, previous, expected_shape",
    [
        ("one_constant_noprev_noadd", ONE, False, (1, 1)),
        ("one_constant_prev_noadd", ONE, True, (1, 1)),
        ("len1_constant_noprev_noadd", LEN1, False, (5, 5)),
        ("len1_constant_prev_noadd", LEN1, True, (5, 5)),
        ("two_constant_prev_noadd", TWO, True, (3, 3)),
    ],
)
def test_degenerate_omstran_gives_the_definitions_zeros(
    case, data, previous, expected_shape
) -> None:
    r_row = GOLDEN.loc[case]
    assert r_row["ok"] == "FALSE" or r_row["dim"] != "x".join(
        map(str, expected_shape)
    )  # R's failure
    np.testing.assert_array_equal(
        _omstran(data, "constant", previous, False), np.zeros(expected_shape)
    )


def test_prob_transition_costs_on_length_one_data_are_refused() -> None:
    assert GOLDEN.loc["len1_prob_noprev_add", "ok"] == "FALSE"
    with pytest.raises(ValueError, match="absolute lag must be smaller than the sequence width"):
        _omstran(LEN1, "prob", False, True)


@pytest.mark.parametrize("method", ["CHI2", "EUCLID"])
def test_the_first_overlapping_episode_is_capped(method: str) -> None:
    name = method.lower()
    assert GOLDEN.loc[f"len1_{name}_step2_overlap", "ok"] == "FALSE"
    np.testing.assert_allclose(
        _quiet(LEN1, method=method, step=2, overlap=True, norm="none"),
        np.sqrt(3.0) * _r_values(f"len1_{name}_step1"),
        rtol=0,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        _quiet(LEN1, method=method, step=2, overlap=True, norm="auto"),
        _r_values(f"len1_{name}_step1_auto"),
        rtol=0,
        atol=1e-12,
    )


def test_the_ledger_records_the_capped_episode_as_traminex_s_correction() -> None:
    """Capping the first episode is one correction of R's indexing; treating a step wider than the data like R's
    step == ncol shortcut (seqdist-CHI2.R:15-16, one episode) is another, with a different raw value. The ledger
    says which one traminex chose (batch-9 review M-10)."""
    ledger = pd.read_csv(
        Path(__file__).resolve().parents[1] / "docs" / "API_PARITY_DECISIONS.csv", dtype=str
    )
    (note,) = ledger.loc[ledger["r_name"] == "seqdist", "notes"]
    assert "traminex's correction" in note
    assert "step == ncol shortcut" in note
