#!/usr/bin/env python3
"""Verify that generated audit/parity artifacts are up to date.

This script regenerates the tracked artifacts into a temporary directory and
compares them against the committed copies. It is meant to be strict enough for
CI/release gates while staying deterministic across platforms.
"""

from __future__ import annotations

import argparse
import difflib
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def _run(cmd: list[str], *, cwd: Path) -> None:
    result = subprocess.run(cmd, cwd=cwd, text=True, capture_output=True)
    if result.returncode != 0:
        msg = [
            f"Command failed ({result.returncode}): {' '.join(cmd)}",
            result.stdout.strip(),
            result.stderr.strip(),
        ]
        raise RuntimeError("\n".join(part for part in msg if part))


def _read_normalized(path: Path) -> str:
    return path.read_text(encoding="utf-8").replace("\r\n", "\n").replace("\r", "\n")


def _compare_file(expected: Path, actual: Path) -> str | None:
    exp = _read_normalized(expected)
    act = _read_normalized(actual)
    if expected.name == "traminex_api_audit.md":
        exp = "\n".join(line for line in exp.splitlines() if not line.startswith("Generated on: "))
        act = "\n".join(line for line in act.splitlines() if not line.startswith("Generated on: "))
    elif expected.name == "traminex_api_audit.json":
        exp_obj = json.loads(exp)
        act_obj = json.loads(act)
        exp_obj.pop("generated_at_utc", None)
        act_obj.pop("generated_at_utc", None)
        exp = json.dumps(exp_obj, indent=2, ensure_ascii=False, sort_keys=True)
        act = json.dumps(act_obj, indent=2, ensure_ascii=False, sort_keys=True)
    if exp == act:
        return None
    diff = "".join(
        difflib.unified_diff(
            exp.splitlines(keepends=True),
            act.splitlines(keepends=True),
            fromfile=str(expected),
            tofile=str(actual),
            n=3,
        )
    )
    return diff[:4000]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Check generated audit/parity artifacts for freshness."
    )
    parser.add_argument(
        "--keep-temp",
        action="store_true",
        help="Keep the temporary generation directory for debugging.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    py = sys.executable
    # mkdtemp has no finalizer: --keep-temp transfers ownership to the caller.
    temp_dir = Path(tempfile.mkdtemp(prefix="traminex-generated-check-"))
    try:
        audit_out = temp_dir / "audit"
        parity_md = temp_dir / "API_PARITY_TABLE.md"
        parity_dir = temp_dir / "parity"
        semantic_md = temp_dir / "API_PARITY_LEDGER.md"
        semantic_json = temp_dir / "traminex_semantic_parity.json"

        _run(
            [
                py,
                str(REPO / "scripts" / "audit_traminex.py"),
                "--output-dir",
                str(audit_out),
            ],
            cwd=REPO,
        )
        _run(
            [
                py,
                str(REPO / "scripts" / "generate_semantic_parity_ledger.py"),
                "--decisions",
                str(REPO / "docs" / "API_PARITY_DECISIONS.csv"),
                "--markdown-output",
                str(semantic_md),
                "--json-output",
                str(semantic_json),
            ],
            cwd=REPO,
        )
        _run(
            [
                py,
                str(REPO / "scripts" / "generate_parity_table.py"),
                "--audit-json",
                str(audit_out / "traminex_api_audit.json"),
                "--output",
                str(parity_md),
            ],
            cwd=REPO,
        )
        _run(
            [
                py,
                str(REPO / "scripts" / "generate_api_mappings.py"),
                "--parity-csv",
                str(audit_out / "traminex_api_audit.csv"),
                "--out-dir",
                str(parity_dir),
            ],
            cwd=REPO,
        )

        comparisons = [
            (
                REPO / "reports" / "audit" / "traminex_api_audit.md",
                audit_out / "traminex_api_audit.md",
            ),
            (
                REPO / "reports" / "audit" / "traminex_api_audit.json",
                audit_out / "traminex_api_audit.json",
            ),
            (
                REPO / "reports" / "audit" / "traminex_api_audit.csv",
                audit_out / "traminex_api_audit.csv",
            ),
            (REPO / "docs" / "API_PARITY_TABLE.md", parity_md),
            (REPO / "docs" / "API_PARITY_TABLE.csv", audit_out / "traminex_api_audit.csv"),
            (
                REPO / "reports" / "parity" / "mapping_traminex_vs_TraMineR.csv",
                parity_dir / "mapping_traminex_vs_TraMineR.csv",
            ),
            (
                REPO / "reports" / "parity" / "mapping_traminex_vs_TraMineRextras.csv",
                parity_dir / "mapping_traminex_vs_TraMineRextras.csv",
            ),
            (REPO / "docs" / "API_PARITY_LEDGER.md", semantic_md),
            (REPO / "reports" / "audit" / "traminex_semantic_parity.json", semantic_json),
        ]

        mismatches: list[str] = []
        for tracked, generated in comparisons:
            if not tracked.exists():
                mismatches.append(
                    f"Missing tracked generated artifact: {tracked.relative_to(REPO)}"
                )
                continue
            if not generated.exists():
                mismatches.append(f"Generator did not produce expected file: {generated}")
                continue
            diff = _compare_file(tracked, generated)
            if diff is not None:
                mismatches.append(
                    f"Stale generated artifact: {tracked.relative_to(REPO)}\n" f"{diff}"
                )

        if mismatches:
            print("Generated artifact freshness check FAILED:")
            for item in mismatches:
                print(f"- {item}")
            print(
                "\nRefresh with:\n"
                "  python scripts/audit_traminex.py\n"
                "  python scripts/generate_parity_table.py\n"
                "  python scripts/generate_api_mappings.py\n"
                "  python scripts/generate_semantic_parity_ledger.py"
            )
            return 1

        print("Generated artifact freshness check PASSED.")
        print(
            "- Checked audit reports, parity table, derived mapping CSVs, "
            "and the semantic decision ledger."
        )
        return 0
    finally:
        if args.keep_temp:
            print(f"Temporary files kept at: {temp_dir}")
        else:
            shutil.rmtree(temp_dir)


if __name__ == "__main__":
    raise SystemExit(main())
