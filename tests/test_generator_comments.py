"""The R generators point to the documented run command (TX-AUD-20261011-064).

Nine generators said "Requires conda R 4.4.3" and 31 pointed to a private "memory note" for
the run environment. The oracle is R 4.6.2, and scripts/r_goldens/README.md gives the run
command (from the repository root, with the oracle library).
"""

from __future__ import annotations

import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def test_no_generator_names_a_stale_environment() -> None:
    stale = re.compile(r"conda R 4\.4\.3|memory note", re.IGNORECASE)
    found = [
        f"{path.name}:{number}"
        for path in sorted((REPO / "scripts" / "r_goldens").glob("*.R"))
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1)
        if stale.search(line)
    ]

    assert found == []
    assert (REPO / "scripts" / "r_goldens" / "README.md").is_file()
