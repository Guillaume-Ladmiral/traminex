"""The documents describe traminex's plotting as it is (TX-AUD-20261019-074).

ECOSYSTEM.md called it "a lightweight matplotlib/reference layer" and "lightweight reference
plots", and docs/PLOTTING.md and docs/API_PARITY.md said the same, while traminex exports 37
plotting functions under ``traminex.visualization`` with R-checked options (docs/PLOTTING.md).
py_seqplot keeps the ggseqplot-style renderers. The documents then called them "R-faithful",
although the ledger records no plotting row as EXACT (batch-18 review M-8): they re-express the
TraMineR plots.
"""

from __future__ import annotations

import csv
import importlib
import inspect
import pkgutil
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
WORDING = "matplotlib renderers that re-express the TraMineR plots"


def test_no_document_calls_the_plotting_layer_lightweight_or_r_faithful() -> None:
    for name in ("ECOSYSTEM.md", "docs/PLOTTING.md", "docs/API_PARITY.md", "CHANGELOG.md"):
        text = " ".join((REPO / name).read_text(encoding="utf-8").split())

        if name != "CHANGELOG.md":
            assert "lightweight" not in text.lower(), name
        assert "R-faithful matplotlib" not in text, name
        assert "R-faithful TraMineR plots" not in text, name
        assert WORDING in text, name


def test_the_ledger_records_no_plotting_row_as_exact() -> None:
    with open(REPO / "docs" / "API_PARITY_DECISIONS.csv", encoding="utf-8") as handle:
        rows = [row for row in csv.DictReader(handle) if row["scope"] == "PLOTTING"]

    assert rows
    for row in rows:
        expected = {"RE_EXPRESSED", "DIVERGENCE"} if row["owner"] == "traminex" else {"SIBLING"}
        assert row["status"] in expected, row["r_name"]
    plotting = " ".join((REPO / "docs" / "PLOTTING.md").read_text(encoding="utf-8").split())
    assert "records each as RE_EXPRESSED or DIVERGENCE, none as EXACT" in plotting


def test_no_plotting_docstring_calls_the_renderer_faithful() -> None:
    # The renderers' docstrings said "faithful TraMineR" after the documents were corrected
    # (TX-REP-20261025-005).
    import traminex.visualization as viz

    found = []
    for info in pkgutil.walk_packages(viz.__path__, viz.__name__ + "."):
        module = importlib.import_module(info.name)
        docs = [(info.name, module.__doc__)] + [
            (f"{info.name}.{name}", obj.__doc__)
            for name, obj in vars(module).items()
            if inspect.isfunction(obj) and obj.__module__ == info.name
        ]
        found += [
            name for name, doc in docs if doc and "faithful TraMineR" in " ".join(doc.split())
        ]

    assert found == []
