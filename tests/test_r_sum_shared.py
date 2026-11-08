"""R's ``sum()`` has one implementation, and its platform dependence is stated (TX-AUD-20261019-069, -11-056).

``disstree.py`` and ``events/_mining.py`` each carried an ``_r_sum`` reproducing R's long-double ``sum()``; they now
share ``traminex._rsum.r_sum``. Its R-identical results need NumPy's ``longdouble`` to be the 80-bit type R's sums use
on x86-64, so the README no longer calls the pure wheel's results portable without that qualification, and the
80-bit golden fixture runs only where ``longdouble`` has R's 63-bit mantissa (a quad-precision ``longdouble``, 112
bits, sums differently).
"""

from __future__ import annotations

import importlib
from pathlib import Path

import numpy as np

from traminex import _rsum

REPO = Path(__file__).resolve().parents[1]


def test_disstree_and_event_mining_share_one_r_sum() -> None:
    disstree_module = importlib.import_module("traminex.disstree")
    mining = importlib.import_module("traminex.events._mining")
    assert disstree_module._r_sum is _rsum.r_sum
    assert mining._r_sum is _rsum.r_sum


def test_r_sum_reads_a_matrix_column_by_column() -> None:
    # Column order adds 1e20 - 1e20 first and keeps both 1s; row order loses the first 1 into 1e20.
    assert _rsum.r_sum(np.array([[1e20, 1.0], [-1e20, 1.0]])) == 2.0
    assert _rsum.r_sum([]) == 0.0


def test_the_80_bit_fixture_runs_only_with_r_s_mantissa() -> None:
    source = (REPO / "tests" / "test_disstree_split_arithmetic.py").read_text(encoding="utf-8")
    assert "np.finfo(np.longdouble).nmant != 63" in source


def test_readme_qualifies_the_pure_wheel() -> None:
    text = " ".join((REPO / "README.md").read_text(encoding="utf-8").split())
    assert "compiler-free and portable" not in text
    # The ICU half follows the owner decision of 2026-10-22 (TX-AUD-20261019-053, option E).
    assert "80-bit `long double`" in text and "pass `alphabet=`" in text
