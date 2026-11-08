"""The maintained stub tree must type-check on its own (audit S7-01, 2026-09-22).

What was wrong: ``stubs/traminex/__init__.pyi`` re-exported ``disstree2dotp``
and ``disstreedisplay`` from ``.disstree`` and ``seqtree2dot`` and
``seqtreedisplay`` from ``.seqtree``, but neither stub module declared them
(the runtime modules ``traminex/disstree.py`` and ``traminex/seqtree.py`` do
re-export all four from ``disstree_viz``). The ``SeqIVarDurResult`` stub also
omitted ``__array_ufunc__``, which ``NDArrayOperatorsMixin`` declares abstract,
so mypy treated the class as abstract. ``MYPYPATH=stubs mypy -p traminex``
reported 5 errors, a user type-checking against the stubs saw those four
names as ``Any``, and nothing in ``make stub-check`` type-checked ``stubs/``.

Source of truth: the runtime modules named above, and mypy itself (the stub
tree is only useful if mypy accepts it). ``make stub-check`` now runs the same
command as this test.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
STUBS = ROOT / "stubs"

pytest.importorskip("mypy")


def test_stub_tree_passes_mypy(tmp_path: Path) -> None:
    env = os.environ.copy()
    env["MYPYPATH"] = str(STUBS)
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "mypy",
            "--config-file",
            str(ROOT / "pyproject.toml"),
            "--cache-dir",
            str(tmp_path / "mypy-cache"),
            "--no-pretty",
            "-p",
            "traminex",
        ],
        cwd=ROOT,
        env=env,
        check=False,
        capture_output=True,
        text=True,
    )
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    # The run must have checked the stub tree, not the runtime package.
    assert "Success: no issues found in" in output, output
    checked = int(output.rsplit("Success: no issues found in", 1)[1].split()[0])
    stub_files = sorted((STUBS / "traminex").rglob("*.pyi"))
    assert checked == len(stub_files), output
