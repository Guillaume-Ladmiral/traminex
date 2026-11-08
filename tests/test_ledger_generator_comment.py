"""The ledger generator's comment states the namespace counts of the oracle (TX-AUD-20261019-057).

``scripts/generate_semantic_parity_ledger.py`` explains why the ledger's denominator is not
"every namespace name": exported names plus registered S3 methods, against every binding in
the two namespaces. It gave 243 and 360 + 66; the ledger's own census gives 167 + 73 = 240,
and live R (TraMineR 2.2-14, TraMineRextras 0.6.9) 358 + 66 bindings.
"""

from __future__ import annotations

import collections
import csv
import re
import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]


def _comment() -> str:
    text = (REPO / "scripts" / "generate_semantic_parity_ledger.py").read_text(encoding="utf-8")
    return " ".join(
        line.strip().lstrip("#") for line in text.splitlines() if line.strip().startswith("#")
    )


def test_the_comment_states_the_ledgers_export_and_s3_census() -> None:
    rows = list(csv.DictReader(open(REPO / "docs" / "API_PARITY_DECISIONS.csv", encoding="utf-8")))
    census = collections.Counter(row["r_kind"] for row in rows)
    found = re.search(
        r"exported \+ registered S3 is (\d+) \((\d+) \+ (\d+)\)", " ".join(_comment().split())
    )

    assert found is not None
    assert tuple(map(int, found.groups())) == (
        census["export"] + census["s3"],
        census["export"],
        census["s3"],
    )


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is not installed")
def test_the_comment_states_the_namespace_binding_counts() -> None:
    found = re.search(
        r"every binding in the two namespaces is (\d+) \+ (\d+)", " ".join(_comment().split())
    )
    expression = (
        'cat(length(ls(asNamespace("TraMineR"), all.names = TRUE)),'
        ' length(ls(asNamespace("TraMineRextras"), all.names = TRUE)))'
    )
    counted = subprocess.run(
        ["Rscript", "-e", expression], capture_output=True, text=True, check=True
    ).stdout.split()

    assert found is not None
    assert list(found.groups()) == counted
