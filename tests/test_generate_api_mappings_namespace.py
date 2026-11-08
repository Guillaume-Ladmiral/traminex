"""``generate_api_mappings.py`` must read NAMESPACE files with the shared parser.

What was wrong (audit S7-17, 2026-09-22): the script split the
``TraMineRextras`` NAMESPACE with its own regular expression, which does not
strip R comments, so a commented-out ``#seqeordplot`` export became the bogus
name ``##seqeordplot``. No row was misclassified at the time, but a commented
export that matched an R name would have moved it to the wrong mapping file.
``scripts/_r_namespace.parse_namespace`` is the parser the audit and the
semantic ledger already use.
"""

from __future__ import annotations

from pathlib import Path

from scripts._r_namespace import parse_namespace
from scripts.generate_api_mappings import extract_extras_names

ROOT = Path(__file__).resolve().parents[1]


def test_extras_names_come_from_the_shared_namespace_parser() -> None:
    namespace = ROOT / "SOURCES" / "TraMineRextras" / "NAMESPACE"

    names = extract_extras_names(namespace)

    assert names == set(parse_namespace(namespace))
    assert not [name for name in names if name.startswith("#")]


def test_commented_exports_are_ignored(tmp_path: Path) -> None:
    namespace = tmp_path / "NAMESPACE"
    namespace.write_text(
        "export(seqa,\n  # seqold,\n  seqb)\n#export(seqgone)\nS3method(plot, seqa)\n",
        encoding="utf-8",
    )

    assert extract_extras_names(namespace) == {"seqa", "seqb", "plot.seqa"}
