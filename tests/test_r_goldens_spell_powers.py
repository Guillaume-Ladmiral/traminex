"""Spell-length powers are R's, bit for bit (TX-AUD-20261019-030).

seqdist raises spell lengths to a power with R's ``^`` (OMslen ``dur.mat ^ (-1 * h)``, seqdist.R:619; OMspell,
NMSMST and SVRspell ``seqdur(...) ^ tpow``, :626), which is the C library's ``pow``. traminex used ``np.power``,
which differs from it in the last bit for 0-76 of 1,000 lengths per exponent, so OMslen and OMspell distances
differed from R's in the last bit. ``scripts/r_goldens/make_goldens_spell_powers.R`` records R's powers (lengths
1-300, seven exponents) and distances in hexadecimal. Its NMSMST, SVRspell and NMS rows are read by
``tests/test_subseq_sum_order.py`` (TX-REP-20261022-001).
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "spell_powers.csv",
    dtype={"kind": str, "key": str, "index": int, "value": str},
)
ROWS = [
    "AABCCDDDDDBBBAAAAAACCCCAABBBCC", "CBBBDDDDDDCCCCCCCCCCCCCCDDDDBB", "ADDDBBBBBBBAAAAAAACCCCCDDDBBBB",
    "AABBBBBBBBBBCCCCCCBDBBBDBBBBBB", "CCCCBAAABBBBBBAAAABBBBBBBBBBAA", "BBBBBAAAAAABBBBBBBCAAAABBBBBBA",
    "AAAAAAAAAAABBAAAAAAAAAAAAAAAAA", "AACCCCCAAAAAAAABDDDDDDDDAAAAAB", "AAAAAAAAAAADDDDDAAAAAACCCAAAAA",
    "BBBBCCCCCCDDDDAAAAAAABBBBBBBBB", "BBBBBBBBCCCCCCCDDCAAADDDDDDBBB", "BBDDDDDDAAAAAACCCCCCCCCCCCCBBD",
    "CCCCCAAAAAAAAAAAADDDDDDDDDDDDD", "DAAAAAAAAAAAAAAAAAAAAAAAAAAABD", "AAAABAAAAAAABBBBBBBCAAAABBDDDD",
    "AAAAADDDDDAAAAACCCBBBBCCCCCCCC", "DDDDDDBBBBBDDDDDDCCBDDDDDBBBBB", "AAAAADDDDDDDDDDDBADDDDDDDDAAAA",
    "BBBBBBBBBBDDDDCDDDDDDBBBBBBBDD", "AABBBBBBBBAABBAAAAAAAAAAAAAACC",
]  # fmt: skip
X = tx.seqdef(np.array([list(row) for row in ROWS], dtype=object), alphabet=["A", "B", "C", "D"])
CASES = {
    "OMslen_h_0.5": dict(method="OMslen", subm="CONSTANT", indel=1.0, h=0.5),
    "OMslen_h_1.3": dict(method="OMslen", subm="CONSTANT", indel=1.0, h=1.3),
    "OMspell_tpow_1.3": dict(method="OMspell", subm="CONSTANT", indel=1.0, expcost=0.5, tpow=1.3),
}


def _r_values(kind: str, key: str) -> np.ndarray:
    rows = GOLDEN[(GOLDEN["kind"] == kind) & (GOLDEN["key"] == key)].sort_values("index")
    return np.array([float.fromhex(value) for value in rows["value"]])


@pytest.mark.parametrize("exponent", ["-2", "-1.3", "-0.5", "0.5", "1.3", "2", "2.5"])
def test_spell_length_powers_are_r_s(exponent: str) -> None:
    from traminex.distances._spells import spell_length_power

    expected = _r_values("power", exponent)
    actual = spell_length_power(np.arange(1, 301, dtype=float), float(exponent))
    assert expected.size == 300
    np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize("case", sorted(CASES))
def test_spell_weighted_distances_are_r_s_bit_for_bit(case: str) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        actual = np.asarray(tx.seqdist(X, **CASES[case]), dtype=float).ravel(order="F")
    np.testing.assert_array_equal(actual, _r_values("dist", case))


def test_r_s_special_cases() -> None:
    """Live R: 3^2 is 9, 5^0 and 1^7.3 are 1, 0^1.3 is 0, and 0^-0.5 and 50^250 are Inf."""
    from traminex.distances._spells import spell_length_power

    assert spell_length_power(np.array([3.0]), 2.0)[0] == 9.0
    assert spell_length_power(np.array([1.0, 5.0]), 0.0).tolist() == [1.0, 1.0]
    assert spell_length_power(np.array([0.0]), -0.5)[0] == np.inf
    assert spell_length_power(np.array([0.0]), 1.3)[0] == 0.0
    assert spell_length_power(np.array([50.0]), 250.0)[0] == np.inf
