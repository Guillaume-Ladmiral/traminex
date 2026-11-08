"""The type-check gate works with the newest dependencies the floors permit (TX-AUD-20261010-010).

NumPy 2.5 requires Python 3.12 and its stubs use 3.12 syntax, so a mypy target
fixed at Python 3.11 could not check any environment that holds the newest
NumPy (``make typecheck`` exited 2 in NumPy's own stubs). mypy now targets the
interpreter it runs under, and a CI job type-checks with the newest NumPy,
pandas and SciPy on both supported Pythons: on 3.11, where pip resolves the
newest NumPy that supports it, mypy also checks the standard library 3.11
users have (``itertools.batched`` or ``Path.walk``, new in 3.12, would pass a
3.12-only check; batch-5 review I-2). The workflow validator requires both.
"""

from __future__ import annotations

import sys
import tomllib
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import validate_workflows  # noqa: E402

CI = REPO / ".github" / "workflows" / "ci.yml"


def test_mypy_targets_the_running_interpreter() -> None:
    config = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))["tool"]["mypy"]
    assert "python_version" not in config
    lint = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))["tool"]
    assert lint["black"]["target-version"] == ["py311"]
    assert lint["ruff"]["target-version"] == "py311"


def _errors(path: Path) -> list[str]:
    return [
        *validate_workflows.validate_workflow(path),
        *validate_workflows.validate_r5_policy(path),
    ]


def test_ci_typechecks_with_the_newest_dependencies(tmp_path) -> None:
    assert _errors(CI) == []
    text = CI.read_text(encoding="utf-8")
    start = text.index("  typecheck-newest:\n")
    end = text.index("\n  test:\n")
    stripped = tmp_path / "ci.yml"
    stripped.write_text(text[:start] + text[end + 1 :], encoding="utf-8")
    errors = _errors(stripped)
    assert any("typecheck-newest" in error for error in errors), errors


def test_the_newest_dependency_job_upgrades_and_typechecks() -> None:
    text = CI.read_text(encoding="utf-8")
    job = text[text.index("  typecheck-newest:\n") : text.index("\n  test:\n")]
    assert "pip install --upgrade numpy pandas scipy" in job
    assert "make typecheck PYTHON=python" in job


def test_the_newest_dependency_job_runs_on_python_3_11_and_3_12(tmp_path) -> None:
    text = CI.read_text(encoding="utf-8")
    job = text[text.index("  typecheck-newest:\n") : text.index("\n  test:\n")]
    assert 'python-version: ["3.11", "3.12"]' in job
    assert "python-version: ${{ matrix.python-version }}" in job
    one_row = tmp_path / "ci.yml"
    one_row.write_text(
        text.replace('python-version: ["3.11", "3.12"]', 'python-version: ["3.12"]', 1),
        encoding="utf-8",
    )
    assert any("3.11" in error and "typecheck-newest" in error for error in _errors(one_row))
