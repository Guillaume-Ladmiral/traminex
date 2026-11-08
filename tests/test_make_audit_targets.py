"""``make audit`` regenerates every artifact DOCS_POLICY.md says it refreshes (TX-AUD-20261010-086).

docs/DOCS_POLICY.md lists ``traminex_api_audit.{csv,json,md}`` and
``traminex_semantic_parity.json`` as "refreshed by ``make audit``", but the target ran only
``audit_traminex.py``, ``generate_parity_table.py`` and ``generate_api_mappings.py``; the JSON's
one writer, ``generate_semantic_parity_ledger.py`` (default outputs docs/API_PARITY_LEDGER.md
and reports/audit/traminex_semantic_parity.json), ran only inside ``make generated-check``.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def test_make_audit_runs_every_generator_of_the_policy_row() -> None:
    policy = (REPO / "docs" / "DOCS_POLICY.md").read_text(encoding="utf-8")
    row = next(line for line in policy.splitlines() if "refreshed by `make audit`" in line)
    dry_run = subprocess.run(
        ["make", "-n", "audit", "PYTHON=python3"],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=True,
    ).stdout

    assert "traminex_semantic_parity.json" in row
    for script in (
        "scripts/audit_traminex.py",
        "scripts/generate_parity_table.py",
        "scripts/generate_api_mappings.py",
        "scripts/generate_semantic_parity_ledger.py",
    ):
        assert f"python3 {script}" in dry_run, script
