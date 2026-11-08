"""The declared scikit-build-core floor can build the sdist (TX-AUD-20261010-047).

pyproject.toml requires ``scikit-build-core>=0.11.0`` but set
``sdist.inclusion-mode``, an option 0.11.0 does not know (build exit 7,
"Unrecognized options"). Without it, 0.11.0 reads the include list with
gitignore rules, under which the bare entry ``README.md`` also matched two
unlisted README files below the root. Without ``inclusion-mode``, though, every
version also reads each nested ``.gitignore`` file, ignored ones included, and
stops at a malformed one (``GitIgnorePatternError``), which a checkout holding
copied tool caches has (batch-5 review I-4). So the option comes back only for
the versions that know it, through an override for scikit-build-core 0.12 and
later, and every root entry is anchored (``/LICENSE``), which 0.11.0 reads as a
root path and the manual mode as the same path. 0.11.0, 0.12.2 and 1.1.0 then
build the same sdist. These checks keep the configuration within 0.11.0's
options; a CI step builds the sdist with scikit-build-core 0.11.0.
"""

from __future__ import annotations

import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PYPROJECT = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
#: tool.scikit-build.sdist options of scikit-build-core 0.11.0 (SDistSettings fields).
SDIST_OPTIONS_0_11 = {"cmake", "exclude", "include", "reproducible"}


def test_the_floor_is_unchanged() -> None:
    assert PYPROJECT["build-system"]["requires"] == ["scikit-build-core>=0.11.0"]


def test_sdist_options_exist_in_the_floor_version() -> None:
    assert set(PYPROJECT["tool"]["scikit-build"]["sdist"]) <= SDIST_OPTIONS_0_11


def test_every_root_include_entry_is_anchored() -> None:
    # A bare name matches at any depth under gitignore rules (0.11.0's reading),
    # so a file of that name anywhere in the tree, ignored ones included, would
    # enter the sdist; an entry with a slash is already relative to the root.
    include = PYPROJECT["tool"]["scikit-build"]["sdist"]["include"]
    assert [entry for entry in include if "/" not in entry] == []


def test_versions_that_know_it_read_no_gitignore_file() -> None:
    overrides = PYPROJECT["tool"]["scikit-build"].get("overrides", [])
    assert {
        "if": {"scikit-build-version": ">=0.12"},
        "sdist": {"inclusion-mode": "manual"},
    } in overrides


def test_ci_builds_the_sdist_with_the_floor_version() -> None:
    ci = (ROOT / ".github" / "workflows" / "ci.yml").read_text(encoding="utf-8")
    assert 'python -m pip install "scikit-build-core==0.11.0"' in ci
