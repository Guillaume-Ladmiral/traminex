#!/usr/bin/env python3
"""Write pyproject.toml's sdist.include list from `git ls-files` and the source policy.

The list stays exact (no globs: scikit-build-core 0.11.0, the declared floor, reads it
with gitignore rules), but it is generated: after adding, removing or renaming a tracked
file, run `make sdist-manifest`. With `--check` the script writes nothing and exits 1
when the list differs, naming each missing or extra entry (TX-AUD-20261019-078).
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path, PurePosixPath

REPO = Path(__file__).resolve().parents[1]

#: Root files in the sdist. 0.11.0 matches a bare name at any depth, so each is
#: written anchored with "/" (TX-AUD-20261010-047).
ROOT_FILES = frozenset(
    {
        ".gitattributes",
        ".pre-commit-config.yaml",
        "CHANGELOG.md",
        "CMakeLists.txt",
        "CONTRIBUTING.md",
        "DUPLICATES.md",
        "ECOSYSTEM.md",
        "LICENSE",
        "MIGRATION_PLAN.md",
        "Makefile",
        "NOTICE",
        "README.md",
        "ROADMAP.md",
        "TraMineR_examples.md",
        "pyproject.toml",
        "traminex_audit_report_20260612.md",
        "verify_audit_improvements.py",
    }
)
#: Every tracked file under these directories is in the sdist.
SOURCE_PREFIXES = (
    ".github/",
    "DATA/",
    "DEMOS/",
    "DOCUMENTATION/",
    "SOURCES/",
    "archive/",
    "benchmarks/",
    "cpp/",
    "examples/",
    "reports/parity/",
    "scripts/",
    "stubs/",
    "tests/",
    "traminex/",
)
#: Of docs/, only these files and the direct children with these suffixes; the
#: subdirectories (plans, prompts, research notes) stay out.
DOCS_FILES = frozenset({"docs/.gitignore", "docs/_quarto.yml"})
DOCS_SUFFIXES = frozenset({".csv", ".md", ".qmd"})

#: Tracked files left out although the rules above take them in: the supervision repository's
#: notes on this project, which hold absolute local paths and are not package documentation.
EXCLUDED_FILES = frozenset({"docs/SUPERVISOR_NOTES.md"})

#: Characters a listed path cannot hold: gitignore pattern and escape characters, the
#: TOML string delimiter, and control characters, which a TOML string cannot hold.
_UNLISTABLE = frozenset('*?[\\"') | frozenset(map(chr, [*range(0x20), 0x7F]))


def _listable(relative: str) -> bool:
    """Whether the path can be an exact entry: no character above, and no trailing space,
    which 0.11.0's gitignore matching drops (batch-17 review M-4)."""
    return not _UNLISTABLE.intersection(relative) and not relative.endswith(" ")


_BLOCK_START = "sdist.include = [\n"
_BLOCK_END = "]\n"


class SdistManifestError(RuntimeError):
    """The tracked files cannot be written as an exact sdist list."""


def is_sdist_source(relative: str) -> bool:
    """Whether the tracked file at ``relative`` (POSIX, from the root) is in the sdist."""
    if relative in EXCLUDED_FILES:
        return False
    if relative in ROOT_FILES or relative in DOCS_FILES:
        return True
    if relative.startswith(SOURCE_PREFIXES):
        return True
    path = PurePosixPath(relative)
    return path.parent == PurePosixPath("docs") and path.suffix in DOCS_SUFFIXES


def _tracked_files(root: Path) -> list[str]:
    output = subprocess.run(
        ["git", "-C", str(root), "ls-files", "-z"], check=True, capture_output=True
    ).stdout
    return [encoded.decode("utf-8") for encoded in output.split(b"\0") if encoded]


def sdist_entries(root: Path = REPO) -> list[str]:
    """The sdist.include entries for the checkout at ``root``, in byte order.

    Raises SdistManifestError for a tracked file missing from the working tree (stage
    its removal with `git rm`) or a path the list cannot hold.
    """
    sources = [relative for relative in _tracked_files(root) if is_sdist_source(relative)]
    unlistable = [relative for relative in sources if not _listable(relative)]
    if unlistable:
        raise SdistManifestError(
            "these tracked paths cannot be listed (a glob, escape, quote or control character,"
            " or a trailing space): " + ", ".join(map(repr, unlistable))
        )
    missing = [relative for relative in sources if not (root / relative).is_file()]
    if missing:
        raise SdistManifestError(
            "tracked but not a file in the working tree (stage a removal with git rm): "
            + ", ".join(missing)
        )
    entries = [relative if "/" in relative else f"/{relative}" for relative in sources]
    return sorted(entries, key=str.encode)


def _block(text: str) -> tuple[int, int]:
    start = text.index(_BLOCK_START) + len(_BLOCK_START)
    return start, text.index(_BLOCK_END, start)


def declared_entries(pyproject_text: str) -> list[str]:
    """The entries of the sdist.include block, one quoted path per line, as written."""
    start, end = _block(pyproject_text)
    return [line.strip().rstrip(",").strip('"') for line in pyproject_text[start:end].splitlines()]


def render(pyproject_text: str, entries: list[str]) -> str:
    """``pyproject_text`` with its sdist.include block replaced by ``entries``."""
    start, end = _block(pyproject_text)
    block = "".join(f'  "{entry}",\n' for entry in entries)
    return pyproject_text[:start] + block + pyproject_text[end:]


def _differences(declared: list[str], expected: list[str]) -> list[str]:
    lines = [f"missing: {entry.lstrip('/')}" for entry in expected if entry not in declared]
    lines += [
        f"not tracked under the policy: {entry.lstrip('/')}"
        for entry in declared
        if entry not in expected
    ]
    if not lines:
        lines.append("the entries are right but not in byte order")
    return lines


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Check only; do not write changes")
    parser.add_argument("--root", type=Path, default=REPO, help="Checkout to read and update")
    args = parser.parse_args(argv)

    pyproject = args.root / "pyproject.toml"
    text = pyproject.read_text(encoding="utf-8")
    expected = sdist_entries(args.root)
    declared = declared_entries(text)
    if declared == expected:
        print(f"sdist.include matches the tracked sources ({len(expected)} paths).")
        return 0
    if args.check:
        print("sdist.include differs from the tracked sources (run `make sdist-manifest`):")
        for line in _differences(declared, expected):
            print(f"  {line}")
        return 1
    pyproject.write_text(render(text, expected), encoding="utf-8")
    print(f"sdist.include now lists {len(expected)} paths.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
