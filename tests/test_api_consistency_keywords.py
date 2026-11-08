"""docs/API_CONSISTENCY.md states which functions spell the missing-state keyword how (TX-AUD-20261010-029).

Section 1 was headed "Missing-state flags standardized", while the public
signatures use ``include_missing`` in 35 functions, ``with_missing`` in 23 and
both in ``seqtabstocc``, and section 1b listed no ``with.missing`` ->
``include_missing`` rename. Section 1 now lists each group; this test compares
the lists with the signatures of the 226 public names, so a new function or a
renamed keyword must update the document.
"""

from __future__ import annotations

import inspect
import re
from pathlib import Path

import traminex as tx
from traminex.api_manifest import PUBLIC_API

DOC = (Path(__file__).resolve().parents[1] / "docs" / "API_CONSISTENCY.md").read_text(
    encoding="utf-8"
)


def _census() -> dict[str, set[str]]:
    groups: dict[str, set[str]] = {"include_missing": set(), "with_missing": set(), "both": set()}
    for name in PUBLIC_API:
        try:
            params = inspect.signature(getattr(tx, name)).parameters
        except (TypeError, ValueError):
            continue
        has_include, has_with = "include_missing" in params, "with_missing" in params
        if has_include and has_with:
            groups["both"].add(name)
        elif has_include:
            groups["include_missing"].add(name)
        elif has_with:
            groups["with_missing"].add(name)
    return groups


def _documented(label: str) -> tuple[int, set[str]]:
    section = DOC[DOC.index("## 1)") : DOC.index("## 1b)")]
    match = re.search(
        rf"^- {re.escape(label)} \((\d+) functions?\): (.+?)(?=^\S|^- |\Z)", section, re.M | re.S
    )
    assert match, f"section 1 has no list for {label}"
    names = set(re.findall(r"`([^`]+)`", match.group(2).split(" — ")[0]))
    return int(match.group(1)), names


def test_section_1_lists_every_spelling() -> None:
    census = _census()
    for label, key in (
        ("`include_missing`", "include_missing"),
        ("`with_missing`", "with_missing"),
        ("Both", "both"),
    ):
        count, names = _documented(label)
        assert names == census[key], (label, sorted(names ^ census[key]))
        assert count == len(census[key])


def test_section_1_no_longer_claims_one_standard() -> None:
    assert "## 1) Missing-state flags standardized" not in DOC
    section_1b = DOC[DOC.index("## 1b)") : DOC.index("## 1c)")]
    assert "`include_missing`" in section_1b
