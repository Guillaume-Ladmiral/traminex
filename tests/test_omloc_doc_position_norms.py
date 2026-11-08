"""Which norms stay within [0, 1] for position-dependent OM (TX-AUD-20261010-064).

``docs/OMLOC.md`` said every TraMineR norm of ``seqdist`` stays within ``[0, 1]`` with ``indel_by_pos`` and
``subm_by_pos``. ``gmean``, ``maxdist`` and ``YujianBo`` do: their bound prices each position at the largest
substitution cost or the largest indel pair, and the alignment paying one of the two per position exists.
``maxlength`` (and ``length``, and ``auto`` for OM) divides by the length times the largest indel only, so it can
exceed 1, in TraMineR as here (TraMineR's own OM gives 2 for AAAA against BBBB with CONSTANT costs and indel 1).
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pytest

import traminex as tx

DOC = Path(__file__).resolve().parents[1] / "docs" / "OMLOC.md"
INDEL_BY_POS = np.array([1.0, 0.6, 0.8, 1.0])


def _documented_subm(states: int) -> np.ndarray:
    subm_by_pos = np.full((4, states, states), 2.0)
    for matrix in subm_by_pos:
        np.fill_diagonal(matrix, 0.0)
    return subm_by_pos


def _seqdist(seqmat, norm: str, indel_by_pos=INDEL_BY_POS, subm_by_pos=None) -> np.ndarray:
    states = 3 if subm_by_pos is None else subm_by_pos.shape[1]
    subm_by_pos = _documented_subm(states) if subm_by_pos is None else subm_by_pos
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        D = tx.seqdist(
            seqmat,
            method="OM",
            indel=1.0,
            subm=2.0,
            indel_by_pos=indel_by_pos,
            subm_by_pos=subm_by_pos,
            norm=norm,
        )
    return np.asarray(D, dtype=float)


def test_the_documented_call_exceeds_1_under_maxlength_only() -> None:
    seqmat = np.array([[0, 0, 1, 2], [0, 1, 1, 2], [2, 1, 0, 0]], dtype=np.int32)
    for norm in ("length", "maxlength", "auto"):
        assert _seqdist(seqmat, norm)[1, 2] == pytest.approx(1.2, rel=1e-12)
    for norm in ("gmean", "maxdist", "YujianBo"):
        D = _seqdist(seqmat, norm)
        assert D.min() >= 0.0 and D.max() <= 1.0


def test_the_documented_costs_give_1_7_for_aaaa_against_bbbb() -> None:
    pair = tx.seqdef(np.array([list("AAAA"), list("BBBB")], dtype=object), alphabet=["A", "B"])
    assert _seqdist(pair, "maxlength", subm_by_pos=_documented_subm(2))[0, 1] == pytest.approx(
        1.7, rel=1e-12
    )


def test_gmean_maxdist_and_yujianbo_stay_within_0_and_1_on_random_position_costs() -> None:
    rng = np.random.default_rng(20261022)
    for _ in range(60):
        seqmat = rng.integers(0, 3, size=(6, 5)).astype(np.int32)
        indel_by_pos = rng.uniform(0.1, 3.0, size=5)
        subm_by_pos = rng.uniform(0.0, 6.0, size=(5, 3, 3))
        subm_by_pos = (subm_by_pos + subm_by_pos.transpose(0, 2, 1)) / 2
        for matrix in subm_by_pos:
            np.fill_diagonal(matrix, 0.0)
        for norm in ("gmean", "maxdist", "YujianBo"):
            D = _seqdist(seqmat, norm, indel_by_pos=indel_by_pos, subm_by_pos=subm_by_pos)
            assert D.min() >= -1e-12 and D.max() <= 1.0 + 1e-12, norm


def test_omloc_md_says_which_norms_stay_within_1() -> None:
    text = " ".join(DOC.read_text(encoding="utf-8").split())
    assert "so they stay within `[0, 1]`" not in text
    assert "`gmean`, `maxdist` and `YujianBo` therefore stay within `[0, 1]`" in text
    assert "The example above gives 1.2" in text and "1.7 for `AAAA` against `BBBB`" in text
