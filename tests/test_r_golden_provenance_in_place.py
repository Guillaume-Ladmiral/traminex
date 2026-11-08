"""``--execute`` refuses to rewrite the live checkout's fixtures (TX-AUD-20261010-020).

``--execute`` runs every generator with the repository as its working
directory, so a byte change is written into the tracked fixtures before the
check fails. Without ``--repo`` the repository was the checkout holding the
script. It now refuses that checkout (also through a symlinked path) and any git
work tree with uncommitted changes under tests/r_goldens, unless
``--allow-in-place`` is given; scripts/r_goldens/README.md asks for a disposable
copy. The refusal comes before any R process starts.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import check_r_golden_provenance as checker  # noqa: E402


@pytest.fixture
def no_execution(monkeypatch):
    calls: list[Path] = []

    def record(repo, **kwargs):
        calls.append(repo)
        return {"schema_version": 1}

    monkeypatch.setattr(checker, "execute_provenance_check", record)
    return calls


def test_execute_refuses_the_checkout(no_execution, capsys) -> None:
    assert checker.main(["--execute"]) == 1
    error = capsys.readouterr().err
    assert "--allow-in-place" in error and str(REPO) in error
    assert no_execution == []


def test_execute_refuses_the_checkout_through_a_symlink(no_execution, tmp_path, capsys) -> None:
    link = tmp_path / "checkout"
    link.symlink_to(REPO, target_is_directory=True)
    assert checker.main(["--execute", "--repo", str(link)]) == 1
    assert "--allow-in-place" in capsys.readouterr().err
    assert no_execution == []


def _git(repo: Path, *args: str) -> None:
    subprocess.run(
        ["git", "-c", "user.name=test", "-c", "user.email=test@example.invalid", *args],
        cwd=repo,
        check=True,
        capture_output=True,
    )


def test_a_work_tree_with_uncommitted_fixture_changes_is_refused(tmp_path) -> None:
    fixtures = tmp_path / "tests" / "r_goldens"
    fixtures.mkdir(parents=True)
    (fixtures / "a.csv").write_text("x\n", encoding="utf-8")
    _git(tmp_path, "init", "-q")
    _git(tmp_path, "add", ".")
    _git(tmp_path, "commit", "-q", "-m", "fixtures")
    checker._refuse_in_place(tmp_path)  # clean: a disposable clone is accepted
    (fixtures / "a.csv").write_text("y\n", encoding="utf-8")
    with pytest.raises(checker.ProvenanceError, match="uncommitted"):
        checker._refuse_in_place(tmp_path)


def test_a_plain_copy_is_accepted(tmp_path) -> None:
    (tmp_path / "tests" / "r_goldens").mkdir(parents=True)
    checker._refuse_in_place(tmp_path)


def test_allow_in_place_passes_the_guard(no_execution) -> None:
    assert checker.main(["--execute", "--allow-in-place"]) == 0
    assert no_execution == [REPO]
