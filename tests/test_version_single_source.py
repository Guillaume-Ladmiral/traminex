"""Guard tests: the package version is single-sourced and cannot drift.

Phase 3 (version single-sourcing). ``traminex.__version__`` must be *derived*
from ``pyproject.toml`` ``[project].version``: through the installed package
metadata (``importlib.metadata.version``) when that distribution owns the
imported files, and otherwise from the ``pyproject.toml`` beside the imported
tree (TX-AUD-20261011-032) -- one declared number, all readers agree.

These tests fail if anyone reintroduces a hand-maintained literal in
``traminex/__init__.py`` (so ``__version__`` would stop tracking metadata) or
lets the declared sources drift apart (the H1-class defect this phase closes).
"""

from __future__ import annotations

import importlib
import importlib.metadata
import re
from pathlib import Path

import traminex

REPO_ROOT = Path(__file__).resolve().parents[1]


def _pyproject_version() -> str:
    text = (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    m = re.search(r'(?m)^version\s*=\s*"([^"]+)"\s*$', text)
    assert m is not None, "could not find [project].version in pyproject.toml"
    return m.group(1)


def test_version_resolves_from_installed_metadata(monkeypatch):
    """The version is computed from importlib.metadata, not a hand-typed literal.

    ``_resolve_version`` must consult ``importlib.metadata.version``; patching it
    must change the result. A hardcoded literal would ignore the patch -- that is
    the "no single source" failure this guard exists to catch. (Tested via the
    helper rather than by reloading the package, so it is robust to other tests
    that mutate ``sys.modules``.)
    """
    monkeypatch.setattr(
        importlib.metadata,
        "version",
        lambda name: (
            "9.9.9-test" if name == "traminex" else importlib.metadata.distribution(name).version
        ),
    )
    assert traminex._resolve_version() == "9.9.9-test"


def test_all_version_readers_agree():
    """__version__ == installed metadata == pyproject.toml [project].version."""
    metadata_version = importlib.metadata.version("traminex")
    assert traminex.__version__ == metadata_version
    assert traminex.__version__ == _pyproject_version()
