"""Tests run at the declared matplotlib floor, 3.6.0 (TX-REP-20261025-002).

matplotlib 3.6.0 names a legend's handles ``legendHandles``; ``legend_handles`` arrived in 3.7.
Batch 8 gave two tests a fallback; a later test read ``legend_handles`` alone and failed in the
floors profile. Every line of the tests that reads ``.legend_handles`` names the fallback too.
"""

from __future__ import annotations

from pathlib import Path

TESTS = Path(__file__).resolve().parent


def test_every_legend_handles_read_has_the_floor_fallback() -> None:
    bare = [
        f"{path.name}:{number}"
        for path in sorted(TESTS.glob("*.py"))
        if path.name != Path(__file__).name
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1)
        if ".legend_handles" in line and "legendHandles" not in line
    ]

    assert bare == []
