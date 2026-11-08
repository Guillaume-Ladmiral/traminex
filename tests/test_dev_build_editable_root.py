"""``docs/DEV_BUILD.md`` says the editable install puts the repository root on ``sys.path``
(TX-AUD-20261011-077).

scikit-build-core's editable ``.pth`` file adds each package's parent directory, which for
this flat layout is the repository root, so ``benchmarks``, ``scripts`` and ``tests`` are
importable from anywhere and can mix with other installed top-level names. The layout
stays (owner decision 2026-10-24: escalate, document); the build notes say so and name
the ``src/`` layout that would end it.
"""

from __future__ import annotations

from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def test_dev_build_says_the_editable_install_exposes_the_root() -> None:
    text = " ".join((REPO / "docs" / "DEV_BUILD.md").read_text(encoding="utf-8").split())

    assert "puts the repository root on `sys.path`" in text
    assert "`import benchmarks`" in text
    assert "`src/` layout (the package under `src/traminex/`)" in text
