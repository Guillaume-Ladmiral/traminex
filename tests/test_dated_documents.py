"""Three dated documents say which of their statements no longer hold (TX-AUD-20261011-030).

DUPLICATES.md said the sibling's ``seqtree.py`` and ``clustassoc.py`` call
``traminex.diss.weighted_pairwise_sum``; at the sibling's e1ca49d only ``_utils.py``,
``_fallback.py`` and ``quality.py`` do. MIGRATION_PLAN.md (written 2026-06-12, completed)
says HAM is length-normalized by default, cites ``seqtransn(normalize=)`` and says
dissassoc's ``permute`` and ``random-sampling`` are not ported; none holds today. The
consumer notice of 2026-09-21 names a re-pin digest that later commits superseded (the
re-pin is owner-deferred with TX-AUD-20261010-099) and cites line numbers of the tree it
read. The dated documents keep their bodies and open with a status note.
"""

from __future__ import annotations

from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def _flat(name: str) -> str:
    return " ".join((REPO / name).read_text(encoding="utf-8").split())


def test_duplicates_names_the_sibling_modules_that_call_the_pairwise_sum() -> None:
    text = _flat("DUPLICATES.md")

    assert "— all now call `traminex.diss.weighted_pairwise_sum`" not in text
    assert (
        "`_utils.py`, `_fallback.py` and `quality.py` call `traminex.diss.weighted_pairwise_sum`"
        in text
    )


def test_the_migration_plan_opens_with_a_dated_status_note() -> None:
    text = _flat("MIGRATION_PLAN.md")
    note = text[text.index("Status note (2026-10-24") :][:1200]

    assert text.index("Status note (2026-10-24") < text.index(
        "This is the written, step-by-step plan"
    )
    for statement in ("HAM", "`seqtransn`", "`norm=`", "`permute`", "`random-sampling`"):
        assert statement in note, statement


def test_the_consumer_notice_opens_with_a_dated_status_note() -> None:
    text = _flat("PYWEIGHTEDCLUSTER_CONSUMER_NOTICE_20260921.md")
    note = text[text.index("Status note (2026-10-24") :][:1200]

    assert text.index("Status note (2026-10-24") < text.index("## 0. CORRECTION")
    assert "TX-AUD-20261010-099" in note
    assert "line numbers" in note
