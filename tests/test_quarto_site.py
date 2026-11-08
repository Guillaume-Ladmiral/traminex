"""The Quarto site publishes the current documents and opens on INDEX (TX-AUD-20261011-054).

``resources: "*.md"`` matched Markdown at any depth, so a render copied 127 historical files
of ``docs/superpowers/`` (33 of them with local paths) into the site, and with no index page
Quarto redirected ``index.html`` to the alphabetically first page, ``API_CONSISTENCY.html``.
Anchoring the pattern to the top level, excluding ``superpowers/`` and rendering
``INDEX.qmd`` first fixed both (renders of 2026-10-24: 172 files to 45, none under
``superpowers/``; ``index.html`` to ``INDEX.html``).
"""

from __future__ import annotations

from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]


def _project() -> dict:
    return yaml.safe_load((REPO / "docs" / "_quarto.yml").read_text(encoding="utf-8"))["project"]


def test_the_site_publishes_top_level_markdown_only() -> None:
    resources = _project()["resources"]

    assert "!superpowers/" in resources
    assert all(entry.startswith(("/", "!")) for entry in resources)


def test_the_site_opens_on_the_index() -> None:
    assert _project()["render"][0] == "INDEX.qmd"
