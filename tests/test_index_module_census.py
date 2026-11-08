"""``docs/INDEX.md`` names every maintained module of the package.

INDEX.md calls itself the authoritative map of the package and says any
disagreement with it is a bug. Cycle-10 N-16 added the 16 top-level modules it
omitted, but the ``events/`` and ``visualization/`` subpackages were still only
summarized, so 23 of their modules (``visualization/series.py`` among them) were
named nowhere (2026-09-22 audit F-06). A top-level module must be named in the
package section; a subpackage module inside its subpackage's entry.
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]


def _package_section() -> list[str]:
    lines = (REPO / "docs" / "INDEX.md").read_text(encoding="utf-8").splitlines()
    start = lines.index("## Package: traminex/")
    end = next(i for i in range(start + 1, len(lines)) if lines[i].startswith("## "))
    return lines[start + 1 : end]


def _entry_block(section: list[str], directory: str) -> str:
    """The top-level ``- `` entry that introduces ``directory/`` and its continuation."""
    for i, line in enumerate(section):
        if line.startswith("- ") and re.search(rf"(?<![\w/]){re.escape(directory)}/", line):
            block = [line]
            for follower in section[i + 1 :]:
                if follower.startswith("- "):
                    break
                block.append(follower)
            return "\n".join(block)
    raise AssertionError(f"docs/INDEX.md has no package entry introducing {directory}/")


MODULES = [
    path
    for path in subprocess.run(
        ["git", "ls-files", "traminex/*.py", "traminex/**/*.py"],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=True,
    ).stdout.split()
    if Path(path).name != "__init__.py"
]


def test_the_census_covers_the_package():
    assert len(MODULES) > 70
    assert "traminex/visualization/series.py" in MODULES


@pytest.mark.parametrize("module", MODULES)
def test_every_module_is_named_in_the_index(module):
    section = _package_section()
    parts = Path(module).relative_to("traminex").parts
    name = parts[-1]
    text = "\n".join(section) if len(parts) == 1 else _entry_block(section, parts[0])
    assert re.search(
        rf"(?<![\w/]){re.escape(name)}", text
    ), f"docs/INDEX.md does not name {module}" + (
        "" if len(parts) == 1 else f" in its {parts[0]}/ entry"
    )
