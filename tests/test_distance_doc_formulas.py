"""The distance documents state OMloc's normalization and the WHAM, DUR and JSD_TRANS formulas (TX-AUD-20261010-069).

``docs/DISTANCES_NORMALIZATION.md`` had no OMloc entry and did not say that OMloc's raw distance depends on ``norm``
(its local indels use ``maxscost``, which YujianBo sets to two indels, as in R); ``docs/DISTANCES.md`` gave the three
traminex-only methods one line each without a formula. The formulas below are computed independently and compared
with traminex, so the documented formulas are the code's.
"""

from __future__ import annotations

import math
import warnings
from pathlib import Path

import numpy as np
import pytest

import traminex as tx

DOCS = Path(__file__).resolve().parents[1] / "docs"
ROWS = [list("AABBC"), list("ABBCC"), list("CCBAA"), list("BBCCA"), list("AAAAA")]
X = tx.seqdef(np.array(ROWS, dtype=object), alphabet=["A", "B", "C"])
CODES = [[{"A": 0, "B": 1, "C": 2}[s] for s in row] for row in ROWS]


def _shares(row: list[int]) -> np.ndarray:
    counts = np.bincount(row, minlength=3).astype(float)
    return counts / counts.sum()


def _transitions(row: list[int]) -> np.ndarray:
    table = np.zeros((3, 3))
    for a, b in zip(row[:-1], row[1:]):
        table[a, b] += 1
    return (table / table.sum()).ravel()


def _formula(method: str, x: list[int], y: list[int], **options) -> float:
    if method == "WHAM":
        w = np.asarray(options.get("weights", [1.0] * len(x)))
        return float(np.sum(w * (np.asarray(x) != np.asarray(y))) / np.sum(w))
    if method == "DUR":
        u, v = _shares(x), _shares(y)
        metric = options["metric"]
        if metric == "l1":
            return float(np.abs(u - v).sum())
        if metric == "cos":
            return 1.0 - float(u @ v) / (np.linalg.norm(u) * np.linalg.norm(v))
        return float(np.sum((u - v) ** 2 / (u + v + 1e-12)))
    p, q, eps = _transitions(x), _transitions(y), 0.0  # the documented default
    m = (p + q) / 2

    def kl(a: np.ndarray, b: np.ndarray) -> float:
        return sum(
            a[k] * (math.log(a[k] + eps) - math.log(b[k] + eps)) for k in range(a.size) if a[k] > 0
        )

    return 0.5 * kl(p, m) + 0.5 * kl(q, m)


@pytest.mark.parametrize(
    "options",
    [
        dict(method="WHAM"),
        dict(method="WHAM", weights=np.array([1.0, 2.0, 0.5, 3.0, 1.0])),
        dict(method="DUR", metric="l1"),
        dict(method="DUR", metric="cos"),
        dict(method="DUR", metric="chi2"),
        dict(method="JSD_TRANS"),
    ],
    ids=["WHAM", "WHAM-weights", "DUR-l1", "DUR-cos", "DUR-chi2", "JSD_TRANS"],
)
def test_the_documented_formulas_are_the_code_s(options: dict) -> None:
    options = dict(options)
    method = options.pop("method")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        actual = np.asarray(tx.seqdist(X, method=method, **options), dtype=float)
    for i in range(len(ROWS)):
        for j in range(len(ROWS)):
            if i != j:
                assert actual[i, j] == pytest.approx(
                    _formula(method, CODES[i], CODES[j], **options), rel=1e-12
                )


def test_distances_md_states_the_formulas() -> None:
    text = (DOCS / "DISTANCES.md").read_text(encoding="utf-8")
    for fragment in (
        "Σ_t w_t·[x_t ≠ y_t] / Σ_t w_t",
        "Σ_k |u_k − v_k|",
        "1 − u·v / (|u|·|v|)",
        "Σ_k (u_k − v_k)² / (u_k + v_k + 10⁻¹²)",
        "½·KL(p‖m) + ½·KL(q‖m)",
    ):
        assert fragment in text, fragment


def test_normalization_md_has_an_omloc_entry_and_its_norm_dependence() -> None:
    text = (DOCS / "DISTANCES_NORMALIZATION.md").read_text(encoding="utf-8")
    assert "- **OMloc** (`OMVIdistance.cpp`)" in text
    assert "OMloc's raw distance depends on `norm`" in text
    assert "normalizes to `2·2.75 / (2.75 + 3.5) = 0.88`" in text


def test_the_documented_omloc_example_is_the_code_s() -> None:
    # A against BBB, subm=1, expcost=0.25: indel 0.75; maxscost 1 under "none", 1.5 under
    # YujianBo, so the local indel is at most I = 0.25*1.5 + 0.5*1 = 0.875 and M = 2*I + 2*I.
    rows = np.array([["A", None, None], ["B", "B", "B"]], dtype=object)
    pair = tx.seqdef(rows, alphabet=["A", "B"])
    raw = np.asarray(
        tx.seqdist(pair, method="OMloc", subm=1.0, expcost=0.25, norm="none"), dtype=float
    )
    yb = np.asarray(
        tx.seqdist(pair, method="OMloc", subm=1.0, expcost=0.25, norm="YujianBo"), dtype=float
    )
    assert raw[0, 1] == pytest.approx(2.5, rel=1e-12)
    assert yb[0, 1] == pytest.approx(2 * 2.75 / (2.75 + 3.5), rel=1e-12)
