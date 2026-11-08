"""Removed dead and test-only helpers stay removed (TX-AUD-20261019-063, TX-AUD-20261010-072).

The 2026-10-19 audit's census found them referenced nowhere, only by stubs, or only by tests
that characterized the helpers themselves; none is in traminex/api_manifest.py or
docs/CONSUMED_SURFACE.md, and the sibling projects' own code names none of them.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
REMOVED = {
    "traminex/events/_engine.py": {"has_event", "count_subsequence", "_count_subseq_map"},
    "traminex/format.py": {"_prepare_id_frame", "_build_bounds_maps_horizontal"},
    "traminex/visualization/ticks.py": {"apply_time_axis_one_based"},
    "traminex/stats_complexity.py": {
        "seqivardur_proxy",
        "seqlogp_proxy",
        "seqmaintokens_proxy",
        "seqprecarity_proxy",
        "seqprecstart_proxy",
        "seqibad_proxy",
        "seqintegr_proxy",
        "seqipos_proxy",
        "_state_rank_map",
    },
    "traminex/_sequence_builders.py": {"build_sts_from_encoded", "sts_kwargs_with_codes"},
}


def _defined_names(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    return {
        node.name
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))
    }


@pytest.mark.parametrize("relative", sorted(REMOVED))
def test_removed_helpers_are_not_defined(relative: str) -> None:
    assert not (REMOVED[relative] & _defined_names(REPO / relative))


def test_removed_helpers_are_not_in_the_stubs() -> None:
    names = set().union(*REMOVED.values())
    for stub in sorted((REPO / "stubs").rglob("*.pyi")):
        assert not (names & _defined_names(stub)), stub.relative_to(REPO)
