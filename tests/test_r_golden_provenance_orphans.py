"""The static provenance check rejects an orphan that replaces a real fixture (TX-AUD-20261010-056).

Renaming ``dissassoc_r2_cases.csv`` to ``dissZZZ_orphan_written_by_nobody.csv``
kept the count and passed: the literal ``"diss"`` (from ``what = "diss"``)
"named" the orphan through the three-character prefix fallback, and nothing
noticed that a fixture a writer spells in full had gone. Now a fixture name a
writer spells in full must exist, and a prefix names a fixture only when it
looks like a built name: it contains ``_`` or ``-`` and the fixture name goes on
at a boundary (the prefix ends in ``_``/``-``, or the next character is ``_``,
``-``, ``.`` or a digit). Every tracked fixture stays attributed.

What the static check proves is attribution by name: a fixture whose name a
writer builds with ``paste0`` is named only by its prefix, so an orphan with
the same prefix still passes, and only ``--execute`` catches it. The README
says so (batch-5 review M-2), and the last test keeps that statement true.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import check_r_golden_provenance as checker  # noqa: E402


def _copy(tmp_path: Path) -> Path:
    for part in ("scripts/r_goldens", "tests/r_goldens"):
        shutil.copytree(REPO / part, tmp_path / part)
    return tmp_path


def _attribution(root: Path) -> dict[str, set[str]]:
    return checker.static_attribution(root, checker.static_inventory(root))


def test_the_tracked_tree_attributes_every_fixture() -> None:
    attribution = _attribution(REPO)
    assert len(attribution) == checker.EXPECTED_FIXTURE_COUNT
    assert all(attribution.values())


def test_an_orphan_swap_fails(tmp_path) -> None:
    root = _copy(tmp_path)
    fixtures = root / "tests" / "r_goldens"
    (fixtures / "dissassoc_r2_cases.csv").rename(fixtures / "dissZZZ_orphan_written_by_nobody.csv")
    with pytest.raises(checker.ProvenanceError) as error:
        _attribution(root)
    assert "dissassoc_r2_cases.csv" in str(error.value)
    assert "dissZZZ_orphan_written_by_nobody.csv" in str(error.value)


def test_a_fixture_spelled_in_full_must_exist(tmp_path) -> None:
    # The renamed file is still "named" by the old stem as a prefix; only the
    # missing full name gives the swap away.
    root = _copy(tmp_path)
    fixtures = root / "tests" / "r_goldens"
    (fixtures / "dissassoc_r2_cases.csv").rename(fixtures / "dissassoc_r2_cases_renamed.csv")
    with pytest.raises(checker.ProvenanceError, match="dissassoc_r2_cases.csv"):
        _attribution(root)


def test_a_short_word_is_not_a_name_prefix(tmp_path) -> None:
    # md_input_d1.csv is built with paste0("md_input_d", k); replacing it by
    # diss_orphan.csv keeps the count, and no full name goes missing.
    root = _copy(tmp_path)
    fixtures = root / "tests" / "r_goldens"
    (fixtures / "md_input_d1.csv").rename(fixtures / "diss_orphan.csv")
    with pytest.raises(checker.ProvenanceError, match="diss_orphan.csv"):
        _attribution(root)


def test_a_built_name_orphan_passes_only_the_static_check(tmp_path) -> None:
    # md_input_d1.csv is built with paste0("md_input_d", k): md_input_d9.csv
    # keeps the count and is named by the same prefix, so the static check
    # passes it, and the README must not claim otherwise.
    root = _copy(tmp_path)
    fixtures = root / "tests" / "r_goldens"
    (fixtures / "md_input_d1.csv").rename(fixtures / "md_input_d9.csv")
    assert _attribution(root)["md_input_d9.csv"] == {"make_goldens_md_inputs.R"}
    readme = (REPO / "scripts" / "r_goldens" / "README.md").read_text(encoding="utf-8")
    assert "`md_input_d9.csv`" in readme
