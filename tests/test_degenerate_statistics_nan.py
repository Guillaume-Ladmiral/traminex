"""Degenerate turbulence, duration-variance and entropy cases give R's NaN.

TraMineR lets IEEE arithmetic run on degenerate inputs; traminex substituted a
number in two places (2026-09-22 audit):

* **S3-11, a row with no valid state** (every cell missing and
  ``include_missing=False``). ``seqivardur`` type 1 is the variance and mean of
  an empty vector, NaN, so ``vmax`` is NaN too; type 2 counts the ``K``
  unvisited states as zero-length spells, so the value and mean are 0 and
  ``vmax`` is ``0 / (K - 2)``: 0 when ``K >= 3`` and NaN (0/0) when ``K <= 2``.
  ``seqST`` computes ``log2(phi * (vmax + 1) / (var + 1))`` from those numbers.
  traminex returned 0 for the type-1 turbulence (and every ``K`` for type 2)
  and NaN for the type-2 ``vmax`` whatever ``K`` was. Source of truth:
  ``seqivardur.R`` lines 25-40 and ``seqST.R`` lines 55-63.
* **S3-10, degenerate normalisers.** With a one-state alphabet ``seqient``
  divides by ``log(1) = 0`` (``seqient.R`` lines 22-25) and ``seqST(norm=True)``
  by ``maxT - 1 = 0`` (``seqST.R`` lines 93-99); when every sequence has length
  1, ``maxT`` is 1 as well. R returns the resulting NaN and ``seqici`` and
  ``seqindic`` inherit it; traminex skipped the division and returned the raw
  value (0 or 1).

Expected values: ``scripts/r_goldens/make_goldens_audit20260922_s3b.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import traminex as tx

GOLDEN = Path(__file__).parent / "r_goldens" / "audit20260922_s3b.txt"


def _r() -> dict[str, tuple[str, str]]:
    out = {}
    for line in GOLDEN.read_text().splitlines():
        if line.startswith("#"):
            continue
        case, kind, payload, _warnings = line.split("\t")
        out[case] = (kind, payload)
    return out


R = _r()


def _values(case: str) -> np.ndarray:
    kind, payload = R[case]
    assert kind == "value", (case, kind, payload)
    return np.array([np.nan if v == "NA" else float(v) for v in payload.split(" ")])


def _same(actual, expected) -> None:
    np.testing.assert_allclose(
        np.asarray(actual, dtype=float), expected, rtol=0, atol=1e-12, equal_nan=True
    )


def _sts(rows, alphabet) -> tx.STSList:
    return tx.seqdef(np.array(rows, dtype=object), alphabet=alphabet, right=None)


def _indic_columns(prefix: str) -> set[str]:
    return {case[len(prefix) :] for case in R if case.startswith(prefix)}


N = None

# --------------------------------------------------------------------------- #
# S3-11: rows with no valid state
# --------------------------------------------------------------------------- #
EMPTY_ROWS = {
    "k1": ([["A", "A", "A"], [N, N, N]], ["A"]),
    "k2": ([["A", "B", "B"], [N, N, N]], ["A", "B"]),
    "k3": ([["A", N, "B", "B"], [N, N, N, N], ["B", N, "C", "C"]], ["A", "B", "C"]),
    "k4": ([["A", "B", "C", "D"], [N, N, N, N]], ["A", "B", "C", "D"]),
}


@pytest.mark.parametrize("name", sorted(EMPTY_ROWS))
@pytest.mark.parametrize("type_", [1, 2])
def test_seqivardur_row_without_a_valid_state_matches_r(name, type_):
    result = tx.seqivardur(_sts(*EMPTY_ROWS[name]), type=type_, include_missing=False)
    prefix = f"s311_{name}_type{type_}_seqivardur"
    _same(result, _values(f"{prefix}_value"))
    _same(result.vmax, _values(f"{prefix}_vmax"))
    _same(result.meand, _values(f"{prefix}_meand"))


@pytest.mark.parametrize("name", sorted(EMPTY_ROWS))
@pytest.mark.parametrize("type_", [1, 2])
@pytest.mark.parametrize("norm", [False, True])
def test_seqst_row_without_a_valid_state_matches_r(name, type_, norm):
    result = tx.seqST(_sts(*EMPTY_ROWS[name]), norm=norm, include_missing=False, type=type_)
    suffix = "_norm" if norm else ""
    _same(result, _values(f"s311_{name}_type{type_}_seqST{suffix}"))


@pytest.mark.parametrize("name", sorted(EMPTY_ROWS))
def test_seqindic_turbulence_and_duration_columns_on_empty_rows(name):
    table = tx.seqindic(
        _sts(*EMPTY_ROWS[name]), indic=["turb", "turbn", "turb2", "turb2n", "dustd", "dustd2"]
    )
    prefix = f"s311_{name}_seqindic_"
    assert set(table.columns) == _indic_columns(prefix)
    for column in table.columns:
        _same(table[column], _values(prefix + column))


# --------------------------------------------------------------------------- #
# S3-10: degenerate normalisers
# --------------------------------------------------------------------------- #
DEGENERATE = {
    "k1": ([["A", "A", "A"], ["A", "A", "A"]], ["A"]),
    "k1miss": ([["A", "A", "A"], ["A", N, "A"]], ["A"]),
    "len1": ([["A"], ["B"], ["A"]], ["A", "B"]),
    "len1k1": ([["A"], ["A"]], ["A"]),
    "control": (
        [["A", "B", "B", "C"], ["C", "C", "A", "A"], ["A", "A", "A", "A"]],
        ["A", "B", "C"],
    ),
}
CASES = [(name, wm) for name in sorted(DEGENERATE) for wm in (False, True)]


@pytest.mark.parametrize(("name", "wm"), CASES)
def test_seqient_and_seqici_with_a_degenerate_normaliser_match_r(name, wm):
    x = _sts(*DEGENERATE[name])
    prefix = f"s310_{name}_wm{int(wm)}"
    _same(tx.seqient(x, norm=True, include_missing=wm), _values(f"{prefix}_seqient_norm"))
    _same(tx.seqient(x, norm=False, include_missing=wm), _values(f"{prefix}_seqient_raw"))
    _same(tx.seqici(x, include_missing=wm), _values(f"{prefix}_seqici"))


@pytest.mark.parametrize(("name", "wm"), CASES)
@pytest.mark.parametrize("type_", [1, 2])
def test_normalised_seqst_with_a_degenerate_maximum_matches_r(name, wm, type_):
    x = _sts(*DEGENERATE[name])
    prefix = f"s310_{name}_wm{int(wm)}_type{type_}"
    _same(tx.seqST(x, norm=False, include_missing=wm, type=type_), _values(f"{prefix}_seqST"))
    _same(tx.seqST(x, norm=True, include_missing=wm, type=type_), _values(f"{prefix}_seqST_norm"))


@pytest.mark.parametrize(("name", "wm"), CASES)
def test_seqindic_inherits_the_degenerate_normalisers(name, wm):
    table = tx.seqindic(
        _sts(*DEGENERATE[name]), indic=["entr", "cplx", "turbn", "turb2n"], with_missing=wm
    )
    prefix = f"s310_{name}_wm{int(wm)}_seqindic_"
    assert set(table.columns) == _indic_columns(prefix)
    for column in table.columns:
        _same(table[column], _values(prefix + column))
