"""Audit tooling must not depend on the current working directory (S7-18).

What was wrong (audit S7-18, 2026-09-22): ``scripts/audit_traminex.py`` read
``SOURCES/...``, ``traminex/...`` and ``pyproject.toml`` relative to the
working directory at import time, so even ``--help`` failed with
``FileNotFoundError`` outside the repository root, and its default output
directory was relative too. ``generate_parity_table.py`` and
``calibrate_parity.py`` likewise defaulted to paths relative to the working
directory, so a run from elsewhere read or wrote the wrong place. The scripts
now anchor these paths to the repository that contains them, as their sibling
``generate_api_mappings.py`` already does.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_audit_script_help_runs_outside_the_repository(tmp_path: Path) -> None:
    result = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "audit_traminex.py"), "--help"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "--output-dir" in result.stdout


def test_parity_table_defaults_are_anchored_to_the_repository(monkeypatch) -> None:
    from scripts import generate_parity_table

    monkeypatch.setattr(sys, "argv", ["generate_parity_table.py"])
    args = generate_parity_table.parse_args()

    assert Path(args.audit_json) == ROOT / "reports" / "audit" / "traminex_api_audit.json"
    assert Path(args.output) == ROOT / "docs" / "API_PARITY_TABLE.md"


def test_calibration_defaults_are_anchored_to_the_repository() -> None:
    from scripts import calibrate_parity

    assert calibrate_parity.DEFAULT_REFS == ROOT / "tests" / "references"
    assert calibrate_parity.DEFAULT_OUT == ROOT / "tests" / "references" / "calibrated_params.json"


@pytest.fixture(scope="module")
def checker_subject(tmp_path_factory):
    root = Path(__file__).resolve().parents[1]
    subject = tmp_path_factory.mktemp("generated-subject")
    tracked = subprocess.check_output(["git", "ls-files", "-z"], cwd=root).decode().split("\0")
    for name in filter(None, tracked):
        source = root / name
        target = subject / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
    return subject


@pytest.mark.parametrize("keep", [False, True])
@pytest.mark.parametrize("stale", [False, True])
def test_generated_checker_diagnostic_lifetime_after_exit(checker_subject, tmp_path, keep, stale):
    artifact = checker_subject / "docs/API_PARITY_TABLE.md"
    original = artifact.read_bytes()
    if stale:
        artifact.write_bytes(original + b"\nOwned stale-artifact test.\n")
    command = [sys.executable, "scripts/check_generated_artifacts.py"]
    if keep:
        command.append("--keep-temp")
    try:
        result = subprocess.run(
            command,
            cwd=checker_subject,
            capture_output=True,
            text=True,
            env={**os.environ, "TMPDIR": str(tmp_path), "PYTHONPATH": str(checker_subject)},
        )
        assert result.returncode == int(stale), result.stdout + result.stderr
        match = re.search(r"Temporary files kept at: (.+)", result.stdout)
        if keep:
            assert match, result.stdout
            retained = Path(match.group(1)).resolve()
            assert retained.parent == tmp_path.resolve()
            assert (retained / "API_PARITY_TABLE.md").is_file()
            assert (retained / "audit/traminex_api_audit.json").is_file()
            assert (retained / "traminex_semantic_parity.json").is_file()
            shutil.rmtree(retained)
        else:
            assert match is None
            assert not list(tmp_path.glob("traminex-generated-check-*"))
    finally:
        artifact.write_bytes(original)
