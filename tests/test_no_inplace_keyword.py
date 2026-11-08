"""No call in the package passes pandas' ``inplace`` keyword (TX-REP-20261025-004).

pandas 3.1 deprecates ``inplace`` (PDEP-8; ``DataFrame.rename`` warns) and pandas 4 removes it.
seqformat's wide-to-long route called ``rename(..., inplace=True)``, so a pandas 3.1 user got a
``Pandas4Warning`` from traminex. traminex supports pandas from 2.1 with no upper bound.
"""

from __future__ import annotations

import ast
from pathlib import Path

PACKAGE = Path(__file__).resolve().parents[1] / "traminex"


def test_no_call_in_the_package_passes_inplace() -> None:
    found = []
    for path in sorted(PACKAGE.rglob("*.py")):
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(node, ast.Call) and any(k.arg == "inplace" for k in node.keywords):
                found.append(f"{path.relative_to(PACKAGE.parent)}:{node.lineno}")

    assert found == []
