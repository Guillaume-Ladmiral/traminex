"""The native ICU requirement is declared, and its absence fails loudly.

``traminex/_categorical_order.py`` loads the system ICU C library through
``ctypes.util.find_library("icui18n")`` so that ``seqstatl``, ``seqdiff`` and the
``stats_basic`` grouping path sort state labels the way R does under ``C.UTF-8``.
ICU appeared in no dependency declaration and in no install document, and pip
cannot supply it, so a host without it met a call-time ``RuntimeError`` naming a
library the user had never been told to install (cycle-6 F-011). On the supported
platform (Ubuntu 24.04) ``libicu74`` ships by default, so the realistic exposure
is a minimal container or a non-Ubuntu host.

This module carries no ICU availability guard, unlike
``tests/test_icu_collation_parity.py``: its subject is the *declaration* and the
*absent* case, both of which must be checked wherever the suite runs.
"""

from __future__ import annotations

import ctypes.util
from pathlib import Path

import pytest

from traminex import _categorical_order

REPO = Path(__file__).resolve().parents[1]
DOCUMENTED_IN = ("README.md", "docs/DEV_BUILD.md")

INPUT = ["b", "A", "a", "B", "Z", "é", "e", "_x", "1", "10", "2"]


@pytest.mark.parametrize("relative", DOCUMENTED_IN)
def test_the_native_icu_requirement_is_documented(relative):
    text = (REPO / relative).read_text(encoding="utf-8")
    assert "ICU" in text, f"{relative} does not name the native ICU requirement"
    assert "libicu74" in text, f"{relative} does not name the installable package"


def test_absent_icu_refuses_loudly_rather_than_falling_back(monkeypatch):
    """The recorded decision: hard failure, not a silent portable fallback.

    Any portable fallback (such as the case-folded key traminex once used for
    event and group labels) would permute the state alphabet -- and therefore
    every per-state result -- without telling the user, because no portable key
    reproduces R's collation for case, accents and punctuation. The owner's
    decision of 2026-10-12 keeps ICU a hard requirement; this test exists so a
    fallback cannot appear incidentally.
    """
    _categorical_order._icu_root_collation_api.cache_clear()
    monkeypatch.setattr(ctypes.util, "find_library", lambda name: None)
    try:
        with pytest.raises(RuntimeError) as excinfo:
            _categorical_order.r_character_order(INPUT)
    finally:
        _categorical_order._icu_root_collation_api.cache_clear()

    message = str(excinfo.value)
    assert "ICU" in message
    assert "libicu" in message, "the refusal must name the package that supplies ICU"
