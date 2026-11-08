"""Every file under ``DATA/`` is its vendored copy, byte for byte (TX-AUD-20261010-098).

``DATA/<package>/<name>`` mirrors ``SOURCES/<package>/data/<name>``, the CRAN sources the
oracle is pinned to. ``ex1.rda``, ``ex2.rda`` and ``famform.rda`` differed in bytes
while holding identical R objects; they are now the vendored files (owner decision
2026-10-24), and this test keeps the two trees from drifting.
"""

from __future__ import annotations

import subprocess
from pathlib import Path, PurePosixPath

import pytest

REPO = Path(__file__).resolve().parents[1]


def _tracked_data_files() -> list[str]:
    output = subprocess.run(
        ["git", "-C", str(REPO), "ls-files", "-z", "--", "DATA"],
        check=True,
        capture_output=True,
    ).stdout
    return sorted(path.decode("utf-8") for path in output.split(b"\0") if path)


def _vendored_copy(relative: str) -> Path:
    _, package, name = PurePosixPath(relative).parts
    return REPO / "SOURCES" / package / "data" / name


def test_data_holds_tracked_files() -> None:
    assert len(_tracked_data_files()) == 10


@pytest.mark.parametrize("relative", _tracked_data_files())
def test_each_data_file_equals_its_vendored_copy(relative: str) -> None:
    vendored = _vendored_copy(relative)

    assert vendored.is_file(), vendored
    assert (REPO / relative).read_bytes() == vendored.read_bytes()
