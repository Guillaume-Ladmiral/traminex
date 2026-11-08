"""The traminex surface the sibling projects import is documented and resolves (TX-AUD-20261010-048).

weightedcluster and py_seqplot import traminex objects by module path, several of
them outside the 226-name manifest (pairwise, pairwise_between_sets, the DissTree
classes, palette/series/ticks helpers, weighted_pairwise_sum), and no document
said so. docs/CONSUMED_SURFACE.md lists every such import of the siblings'
package sources; this test requires each listed path to resolve, so renaming or
moving one fails here before it breaks a sibling.
"""

from __future__ import annotations

import importlib
import re
from pathlib import Path

import pytest

DOC = (Path(__file__).resolve().parents[1] / "docs" / "CONSUMED_SURFACE.md").read_text(
    encoding="utf-8"
)
ROW = re.compile(r"^\| `(traminex[\w.]*)` \| `([\w()]+)` \| ([^|]+) \|", re.M)
ROWS = ROW.findall(DOC)


def test_the_table_lists_each_path_once() -> None:
    assert len(ROWS) == 19
    assert len({(module, name) for module, name, _ in ROWS}) == len(ROWS)


@pytest.mark.parametrize("module, name, consumers", ROWS, ids=[f"{m}.{n}" for m, n, _ in ROWS])
def test_every_consumed_path_resolves(module: str, name: str, consumers: str) -> None:
    target = importlib.import_module(module)
    assert hasattr(target, name), f"{module}.{name} (used by {consumers.strip()})"
