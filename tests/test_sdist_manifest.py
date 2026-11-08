"""The sdist list is generated from git and the policy (TX-AUD-20261019-078).

``scripts/sync_sdist_manifest.py`` writes ``sdist.include`` in ``pyproject.toml`` from
``git ls-files`` and the source policy, so a new tracked file needs one command, not a
hand edit. These tests run it on a small git repository.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from scripts import sync_sdist_manifest

_PYPROJECT = """\
[tool.scikit-build]
sdist.exclude = ["**"]
sdist.include = [
{entries}]

[tool.black]
line-length = 100
"""


def _render(entries: list[str]) -> str:
    return _PYPROJECT.format(entries="".join(f'  "{entry}",\n' for entry in entries))


def _repository(tmp_path: Path, tracked: list[str], untracked: list[str]) -> Path:
    root = tmp_path / "repository"
    root.mkdir()
    subprocess.run(["git", "init", "-q", str(root)], check=True)
    for relative in [*tracked, *untracked]:
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("x\n", encoding="utf-8")
    subprocess.run(["git", "-C", str(root), "add", "--", *tracked], check=True)
    return root


_TRACKED = [
    "README.md",
    "pyproject.toml",
    "docs/INDEX.md",
    "docs/_quarto.yml",
    "docs/superpowers/plans/plan.md",
    "notes/outside.md",
    "tests/test_a.py",
    "traminex/__init__.py",
]
_EXPECTED = [
    "/README.md",
    "/pyproject.toml",
    "docs/INDEX.md",
    "docs/_quarto.yml",
    "tests/test_a.py",
    "traminex/__init__.py",
]


def test_the_list_holds_tracked_files_under_the_policy_only(tmp_path: Path) -> None:
    root = _repository(tmp_path, _TRACKED, ["tests/test_untracked.py", "LICENSE"])

    assert sync_sdist_manifest.sdist_entries(root) == _EXPECTED


def test_check_names_a_missing_and_an_extra_entry(tmp_path: Path, capsys) -> None:
    root = _repository(tmp_path, _TRACKED, [])
    declared = [entry for entry in _EXPECTED if entry != "tests/test_a.py"]
    (root / "pyproject.toml").write_text(_render([*declared, "tests/gone.py"]), encoding="utf-8")

    status = sync_sdist_manifest.main(["--root", str(root), "--check"])

    output = capsys.readouterr().out
    assert status == 1
    assert "missing: tests/test_a.py" in output
    assert "not tracked under the policy: tests/gone.py" in output


def test_check_names_an_order_difference(tmp_path: Path, capsys) -> None:
    root = _repository(tmp_path, _TRACKED, [])
    (root / "pyproject.toml").write_text(_render(_EXPECTED[::-1]), encoding="utf-8")

    status = sync_sdist_manifest.main(["--root", str(root), "--check"])

    assert status == 1
    assert "not in byte order" in capsys.readouterr().out


def test_writing_the_list_replaces_only_the_include_block(tmp_path: Path) -> None:
    root = _repository(tmp_path, _TRACKED, [])
    (root / "pyproject.toml").write_text(_render(["tests/gone.py"]), encoding="utf-8")

    assert sync_sdist_manifest.main(["--root", str(root)]) == 0

    assert (root / "pyproject.toml").read_text(encoding="utf-8") == _render(_EXPECTED)
    assert sync_sdist_manifest.main(["--root", str(root), "--check"]) == 0


def test_a_tracked_file_missing_from_the_tree_is_refused(tmp_path: Path) -> None:
    root = _repository(tmp_path, _TRACKED, [])
    (root / "tests" / "test_a.py").unlink()

    with pytest.raises(sync_sdist_manifest.SdistManifestError, match="tests/test_a.py"):
        sync_sdist_manifest.sdist_entries(root)


# Batch-17 review M-4: 0.11.0's gitignore matching drops a trailing space, and a control
# character makes pyproject.toml invalid TOML.
@pytest.mark.parametrize(
    "name",
    [
        "tests/test_[a].py",
        "tests/test_*.py",
        'tests/q"uote.py',
        "tests/trailing.py ",
        "tests/control\x01.py",
        "tests/tab\t.py",
    ],
)
def test_a_path_the_list_cannot_hold_is_refused(tmp_path: Path, name: str) -> None:
    root = _repository(tmp_path, [*_TRACKED, name], [])

    with pytest.raises(sync_sdist_manifest.SdistManifestError, match="cannot be listed"):
        sync_sdist_manifest.sdist_entries(root)
