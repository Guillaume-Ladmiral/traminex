"""``docs/INDEX.md`` names what the repository tracks (TX-AUD-20261019-064, -076; TX-AUD-20261010-087).

INDEX.md calls itself the authoritative map, yet it omitted governance files, the PROMPT
records at the root, documents, scripts and ``docs/superpowers/``, described ``_spells.py``
and ``omloc.py`` wrongly, left ``auto`` out of the ``TRAMINEX_BINDING`` values, and, with
README.md, presented the 2026-06-12 report as the closeout audit. The census below reads
``git ls-files``: a new or renamed tracked file that INDEX.md does not name fails it.
``tests/test_index_module_census.py`` covers the package's modules.
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
INDEX = (REPO / "docs" / "INDEX.md").read_text(encoding="utf-8")


def _tracked() -> list[str]:
    output = subprocess.run(
        ["git", "-C", str(REPO), "ls-files", "-z"], check=True, capture_output=True
    ).stdout
    return [path.decode("utf-8") for path in output.split(b"\0") if path]


def _section(heading: str) -> str:
    start = INDEX.index(f"\n## {heading}\n")
    end = INDEX.find("\n## ", start + 1)
    return INDEX[start : end if end != -1 else len(INDEX)]


def _subjects(section: str) -> str:
    """The names each bullet introduces: its text before the first " (", continuation lines
    joined. A name in another entry's description does not count (batch-18 review M-2)."""
    bullets: list[str] = []
    for line in section.splitlines():
        stripped = line.strip()
        if stripped.startswith("- "):
            bullets.append(stripped[2:])
        elif bullets and line.startswith(" ") and stripped:
            bullets[-1] += " " + stripped
    return "\n".join(bullet.split(" (", 1)[0] for bullet in bullets)


def _names(name: str, text: str) -> bool:
    """Whether ``text`` names the file ``name`` as a whole name (not inside a path or a longer name)."""
    return re.search(rf"(?<![\w/.-]){re.escape(name)}(?![\w-]|\.\w)", text) is not None


def _names_paired(name: str, text: str) -> bool:
    """A paired file may be written ``STEM.md / .qmd`` (INDEX.md's own convention)."""
    stem, _, suffix = name.rpartition(".")
    paired = rf"(?<![\w/.-]){re.escape(stem)}\.\w+ / \.{re.escape(suffix)}(?![\w])"
    return _names(name, text) or re.search(paired, text) is not None


def _children(directory: str) -> set[str]:
    prefix = f"{directory}/"
    return {
        path[len(prefix) :].split("/")[0] + ("/" if "/" in path[len(prefix) :] else "")
        for path in _tracked()
        if path.startswith(prefix)
    }


def test_index_names_every_tracked_root_file() -> None:
    listed = _subjects(_section("Root")) + _subjects(_section("Prompt and verification records"))
    missing = [path for path in _tracked() if "/" not in path and not _names(path, listed)]

    assert missing == []


def test_index_names_every_tracked_top_level_directory() -> None:
    directories = {path.split("/")[0] for path in _tracked() if "/" in path}
    missing = sorted(
        d for d in directories if re.search(rf"(?<![\w/.-]){re.escape(d)}/", INDEX) is None
    )

    assert missing == []


def test_index_names_every_entry_of_docs_and_scripts() -> None:
    docs = _subjects(_section("Docs: docs/"))
    scripts = _subjects(_section("Scripts: scripts/"))
    missing = [
        f"docs/{name}" for name in sorted(_children("docs")) if not _names_paired(name, docs)
    ]
    missing += [
        f"scripts/{name}" for name in sorted(_children("scripts")) if not _names(name, scripts)
    ]

    assert missing == []


def test_index_describes_the_spell_helpers_and_omloc() -> None:
    def entry(name: str) -> str:
        start = INDEX.index(f"- {name} (")
        return INDEX[start : INDEX.index("\n  - ", start)]

    spells, omloc = entry("_spells.py"), entry("omloc.py")

    assert "kernels" not in spells
    assert "spell_encode" in spells
    assert "`omloc()`" in omloc and "OMLOC.md" in omloc


def _active_documents() -> dict[str, str]:
    names = ["README.md", "CONTRIBUTING.md", "ROADMAP.md"]
    names += [f"docs/{path.name}" for path in sorted((REPO / "docs").glob("*.md"))]
    return {name: (REPO / name).read_text(encoding="utf-8") for name in names}


def test_the_documents_give_the_binding_values_the_code_accepts() -> None:
    # Batch-18 review M-1: docs/DISTANCES.md also gave python|pybind11.
    code = (REPO / "traminex" / "distances" / "__init__.py").read_text(encoding="utf-8")
    accepted = re.search(r"if _backend not in \{([^}]*)\}", code)
    assert accepted is not None
    values = {value.strip().strip('"') for value in accepted.group(1).split(",")}
    written = {
        name: re.findall(r"TRAMINEX_BINDING=([a-z0-9]+(?:\|[a-z0-9]+)+)", text)
        for name, text in _active_documents().items()
    }

    assert written["docs/INDEX.md"]
    wrong = {
        name: found
        for name, found in written.items()
        if any(set(s.split("|")) != values for s in found)
    }
    assert wrong == {}
    assert written["docs/INDEX.md"][0].split("|")[0] == "auto"


def test_the_2026_06_12_audit_is_labelled_historical() -> None:
    # Batch-18 review I-1: docs/README.md also called it the closeout audit.
    readme = (REPO / "README.md").read_text(encoding="utf-8")
    mentions = [
        line
        for text in _active_documents().values()
        for line in text.splitlines()
        if "traminex_audit_report_20260612.md" in line
    ]

    assert mentions
    assert all("closeout audit" not in line.lower() or "historical" in line for line in mentions)
    assert "Prompt and verification records" in readme
