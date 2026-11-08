"""seqidegrad reads a logical pow as R does (TX-AUD-20261010-004).

R's public seqidegrad has no spell.integr argument: a logical pow sets it and is
then coerced to numeric (seqidegrad.R:10-15). traminex also offers
``spell_integr``; left at None it follows R's rule, and an explicit value is
kept. Fixture: ``scripts/r_goldens/make_goldens_seqidegrad_logical_pow.R``.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
import pytest
from r_goldens_utils import R_GOLDEN_ATOL

import traminex as tx

GOLDEN = Path(__file__).resolve().parent / "r_goldens" / "seqidegrad_logical_pow.csv"
METHODS = ["RANK", "TRATEDSS", "FREQ"]


def _golden() -> dict[tuple[str, bool], np.ndarray]:
    values: dict[tuple[str, bool], list[tuple[int, float]]] = defaultdict(list)
    with GOLDEN.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            values[(row["method"], row["pow"] == "TRUE")].append(
                (int(row["row"]), float(row["value"]))
            )
    return {key: np.array([value for _, value in sorted(items)]) for key, items in values.items()}


R = _golden()
SEQS = tx.seqdef(
    np.array([list("AABBC"), list("CCBAA"), list("ABCBA"), list("BBBCC")], dtype=object),
    alphabet=["A", "B", "C"],
)
CASES = [(method, pow_value) for method in METHODS for pow_value in (False, True)]


def _degrad(**kwargs) -> np.ndarray:
    return np.asarray(tx.seqidegrad(SEQS, **kwargs), dtype=float).ravel()


def test_every_generated_case_is_exercised():
    assert set(CASES) == set(R)


@pytest.mark.parametrize(("method", "pow_value"), CASES)
def test_logical_pow_matches_r(method: str, pow_value: bool):
    np.testing.assert_allclose(
        _degrad(method=method, pow=pow_value), R[(method, pow_value)], rtol=0, atol=R_GOLDEN_ATOL
    )


@pytest.mark.parametrize(("method", "pow_value"), CASES)
def test_numpy_bool_pow_matches_r(method: str, pow_value: bool):
    np.testing.assert_allclose(
        _degrad(method=method, pow=np.bool_(pow_value)),
        R[(method, pow_value)],
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


@pytest.mark.parametrize(("method", "pow_value"), CASES)
def test_seqindic_degrad_agrees_with_r(method: str, pow_value: bool):
    table = tx.seqindic(SEQS, indic=["degrad"], prec_args={"method": method, "pow": pow_value})
    np.testing.assert_allclose(
        table["Degrad"].to_numpy(dtype=float), R[(method, pow_value)], rtol=0, atol=R_GOLDEN_ATOL
    )


@pytest.mark.parametrize("method", METHODS)
def test_explicit_spell_integr_is_kept(method: str):
    # traminex's own argument wins over the logical-pow rule when given.
    np.testing.assert_allclose(
        _degrad(method=method, pow=False, spell_integr=True),
        _degrad(method=method, pow=0.0),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )
    np.testing.assert_allclose(
        _degrad(method=method, pow=False, spell_integr=False),
        R[(method, False)],
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


@pytest.mark.parametrize("method", METHODS)
def test_numeric_pow_keeps_spell_integration(method: str):
    np.testing.assert_allclose(
        _degrad(method=method, pow=1.0), R[(method, True)], rtol=0, atol=R_GOLDEN_ATOL
    )
