#!/usr/bin/env python3
"""Check that the release gates can run in this interpreter (TX-AUD-20261010-046).

``make release-check`` imports the development tools and an installed
traminex (``scripts/check_docs_consistency.py`` imports it). Without the
documented development install the gates stopped with a bare
``ModuleNotFoundError``; this check runs first and names every missing module
and the install command of docs/DEV_BUILD.md.
"""

from __future__ import annotations

import sys
from importlib.util import find_spec

#: Module -> the gate that needs it.
REQUIRED = {
    "black": "format-check",
    "IPython": "format-check (black's notebook support)",
    "tokenize_rt": "format-check (black's notebook support)",
    "ruff": "lint",
    "mypy": "typecheck, stub-check",
    "yaml": "workflow-check",
    "pytest": "release-check's tests",
    "traminex": "docs-check, stub-check (an installed traminex)",
}
INSTALL = "-m pip install -e '.[plot,perf,extras,tree,dev]'"


def main() -> int:
    missing = {name: purpose for name, purpose in REQUIRED.items() if find_spec(name) is None}
    if missing:
        print("The release gates cannot run in this interpreter; missing modules:")
        for name, purpose in missing.items():
            print(f"- {name} (needed by {purpose})")
        print(
            f"Install the development profile from the repository root: {sys.executable} {INSTALL} "
            "(docs/DEV_BUILD.md)."
        )
        return 1
    print("development environment OK: " + ", ".join(REQUIRED))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
