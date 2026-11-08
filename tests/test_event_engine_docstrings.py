"""The event engine's module docstrings describe the port as it is (TX-AUD-20261010-078).

They said the engine reproduces R's C++ "line for line" and that "Everything is faithful", while the parity ledger
records corrections (the seqedist and seqefsub rows); they called the ``-1`` mapping a ``+/-inf`` sentinel, where R's
constraint.cpp sets ``+/-DBL_MAX``; and the mining module said later commits would add functions it already holds.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from traminex.events import _engine, _mining


@pytest.mark.parametrize("phrase", ["line for line", "Everything is faithful"])
def test_the_engine_docstring_makes_no_blanket_fidelity_claim(phrase: str) -> None:
    assert phrase not in " ".join(_engine.__doc__.split())


def test_the_engine_docstring_names_r_s_bound_and_the_recorded_corrections() -> None:
    text = " ".join(_engine.__doc__.split())
    assert "DBL_MAX" in text
    assert "parity ledger" in text


def test_the_mining_docstring_names_only_functions_it_holds() -> None:
    text = " ".join(_mining.__doc__.split())
    assert "Later RP5b commits add" not in text
    named = set(re.findall(r":(?:func|class):`(\w+)`", text))
    assert named and all(hasattr(_mining, name) for name in named), sorted(
        n for n in named if not hasattr(_mining, n)
    )


def test_the_module_index_makes_no_line_for_line_claim() -> None:
    # docs/INDEX.md described _engine.py with the claim the docstring dropped (batch-11 review M-9).
    index = Path(__file__).resolve().parents[1] / "docs" / "INDEX.md"
    text = " ".join(index.read_text(encoding="utf-8").split())
    entry = text[text.index("_engine.py (") :].split(")", 1)[0]
    assert "line-for-line" not in entry and "line for line" not in entry
    assert "parity ledger" in entry
