"""Every public function that needs the native ICU library is documented as such.

``_categorical_order.r_character_order`` (and ``r_character_rank`` on top of it)
loads the system ICU library to sort labels as R does, and raises
``RuntimeError`` when ICU is absent. ``docs/DEV_BUILD.md`` and ``README.md``
list the functions that can reach it, but the list drifted: ``seqimplic`` and
``seqsurv`` gained R-collation ordering on 2026-09-22 (audit S3-15, S3-22)
without being added (audit F-08). This test finds every call site in the
package and requires its public entry point in both documents.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
# r_factor_levels and stable_sorted_unique order text with r_character_rank
# (TX-AUD-20261010-039, TX-AUD-20261011-015/-017), and r_seqstatl_order with
# r_character_order (the default alphabet of every constructor, owner decision
# D2 of 2026-10-19; TX-AUD-20261010-001, -031), and r_order_keys (R's order() of a
# text sortv, moved from the plotting module; TX-REP-20261021-001) with
# r_character_rank, and r_factor (R's factor(), batch 12) with r_factor_levels, so
# their callers need ICU too.
HELPERS = {
    "r_character_order",
    "r_character_rank",
    "r_factor_levels",
    "r_factor",
    "stable_sorted_unique",
    "r_seqstatl_order",
    "r_order_keys",
}
# Private helpers that call ICU, and the public function they serve.
PRIVATE_ENTRY = {
    "_combined_registry": "seqMD",
    "_ordinary_group": "seqdiff",
    "_seqimplic_levels": "seqimplic",
    "_seqsurv_factor": "seqsurv",
    "_rowmode_scalar": "rowmode",
    "_seqtabstocc_intermediate": "seqtabstocc",
    "resolve_sort_order": "sequence_index_plot",
    "_resolve_axis": "dissrf",
    "_factor": "dissassoc",
    "_renumber": "dissmergegroups",
    "_node_factor": "disstree",
    "_group_factor": "seqecmpgroup",
    "_compare_factor": "dissCompare",
    # The grouped plots' factor levels (batch-2 review M-1): the shared panel
    # helper serves seqedplot and seqHtplot, and the entropy series helper (not
    # a manifest name) serves seqplot_tentrop.
    "_factor_panels": "seqedplot",
    "_r_factor_of": "seqedplot",
    "_hspell_horizontal_sts": "HSPELL_to_STS",
    # traminex/_sequence_builders.py (a private module): the default alphabet of
    # the long constructors (seqdef_from_long, TSE_to_STS).
    "infer_alphabet": "seqdef_from_long",
    # seqformat's STS output: states the input does not declare take R's order.
    "_materialize_sts": "seqformat",
    "grouped_transversal_entropy": "seqplot_tentrop",
}


def _icu_entry_points() -> set[str]:
    entries: set[str] = set()
    for path in sorted((REPO / "traminex").rglob("*.py")):
        if path.name == "_categorical_order.py":
            continue
        for top in ast.parse(path.read_text(encoding="utf-8")).body:
            if not isinstance(top, (ast.FunctionDef, ast.ClassDef)):
                continue
            calls = any(
                isinstance(node, ast.Call) and getattr(node.func, "id", None) in HELPERS
                for node in ast.walk(top)
            )
            if calls:
                name = PRIVATE_ENTRY.get(top.name, top.name)
                assert not name.startswith("_"), (
                    f"{path.relative_to(REPO)}::{top.name} calls ICU; map it to its "
                    "public function in PRIVATE_ENTRY and document that function"
                )
                entries.add(name)
    return entries


ENTRIES = _icu_entry_points()


def test_the_scan_finds_the_known_icu_users():
    assert {"seqstatl", "seqdiff", "seqMD", "FCE_to_TSE", "seqimplic", "seqsurv"} <= ENTRIES


@pytest.mark.parametrize("relative", ["docs/DEV_BUILD.md", "README.md"])
def test_every_icu_entry_point_is_named_in_the_install_documents(relative):
    text = (REPO / relative).read_text(encoding="utf-8")
    missing = sorted(name for name in ENTRIES if f"`{name}`" not in text)
    assert not missing, f"{relative} does not name these ICU-dependent functions: {missing}"


def test_the_scan_sees_the_r_factor_callers():
    """r_factor reaches ICU through r_factor_levels (batch-12 review M-1)."""
    assert {"seqrep_grp", "dissmergegroups", "dissassoc", "disstree", "dissCompare"} <= ENTRIES


def test_every_private_entry_names_a_function_of_the_package():
    """A stale key would map nothing, so a renamed ICU caller would go unchecked (batch-12 review M-1)."""
    defined = set()
    for path in sorted((REPO / "traminex").rglob("*.py")):
        for top in ast.parse(path.read_text(encoding="utf-8")).body:
            if isinstance(top, (ast.FunctionDef, ast.ClassDef)):
                defined.add(top.name)
    assert sorted(set(PRIVATE_ENTRY) - defined) == []
