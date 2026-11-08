"""The documents describe the CI workflow as configured, not as having run (TX-AUD-20261019-067).

The repository has no git remote, so the GitHub Actions workflows under ``.github/workflows/``
have never run; several documents said that CI runs, builds or type-checks. The active
documents, and the unreleased part of the changelog, now say what the workflow is configured
to do and that no CI run exists. Released changelog sections keep their history.
"""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

#: "CI" followed by a verb of something done, in the present tense, or CI named as a gate.
_RUN_CLAIM = re.compile(
    r"\b(?:Active CI|CI(?: job)? (?:runs|builds|prepares|type-checks)|runs in CI|CI gates)\b"
)


def _active_texts() -> dict[str, str]:
    texts = {
        name: (ROOT / name).read_text(encoding="utf-8")
        for name in ("README.md", "CONTRIBUTING.md", "ROADMAP.md")
    }
    for path in sorted((ROOT / "docs").glob("*.md")):
        texts[f"docs/{path.name}"] = path.read_text(encoding="utf-8")
    changelog = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    texts["CHANGELOG.md [Unreleased]"] = changelog[: changelog.index("\n## [0.0.8]")]
    return texts


def test_no_active_document_says_ci_has_run() -> None:
    claims = [
        f"{name}: {match.group(0)!r}"
        for name, text in _active_texts().items()
        for match in _RUN_CLAIM.finditer(" ".join(text.split()))
    ]

    assert claims == []


def test_dev_build_says_the_workflows_have_not_run() -> None:
    text = " ".join((ROOT / "docs" / "DEV_BUILD.md").read_text(encoding="utf-8").split())

    assert "has no git remote" in text
    assert "have never run" in text
