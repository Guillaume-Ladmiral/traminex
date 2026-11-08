from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from scripts.generate_semantic_parity_ledger import LedgerError, parse_namespace

REPO = Path(__file__).resolve().parents[1]


def _namespace_file(tmp_path: Path, text: str) -> Path:
    path = tmp_path / "NAMESPACE"
    path.write_text(text, encoding="utf-8")
    return path


def test_inline_r_comments_inside_export_are_ignored(tmp_path: Path) -> None:
    path = _namespace_file(
        tmp_path,
        """
        export(
            kept,
            also.kept, ## obsolete
            final_kept
        )
        """,
    )

    assert parse_namespace(path) == {
        "kept": "export",
        "also.kept": "export",
        "final_kept": "export",
    }


def test_fully_commented_export_entries_are_ignored(tmp_path: Path) -> None:
    path = _namespace_file(
        tmp_path,
        """
        export(
            visible,
            ## entirely_commented,
            still_visible
        )
        """,
    )

    assert parse_namespace(path) == {
        "visible": "export",
        "still_visible": "export",
    }


def test_r_comments_start_only_outside_quoted_tokens(tmp_path: Path) -> None:
    path = _namespace_file(
        tmp_path,
        """
        export(
            "literal#hash",
            visible # hidden_comment_token
        )
        """,
    )

    assert parse_namespace(path) == {
        "literal#hash": "export",
        "visible": "export",
    }


def test_quoted_s3_generic_and_class_names_are_canonical(tmp_path: Path) -> None:
    path = _namespace_file(
        tmp_path,
        """
        S3method("[", "stslist")
        S3method("as.character", "seqelist")
        """,
    )

    assert parse_namespace(path) == {
        "[.stslist": "s3",
        "as.character.seqelist": "s3",
    }


def test_quoted_setter_exports_include_callable_alias(tmp_path: Path) -> None:
    path = _namespace_file(tmp_path, 'export("alphabet<-", ordinary)\n')

    assert list(parse_namespace(path).items()) == [
        ("alphabet<-", "export"),
        ("alphabet", "export"),
        ("ordinary", "export"),
    ]


def test_export_s3_collisions_remain_errors(tmp_path: Path) -> None:
    path = _namespace_file(
        tmp_path,
        """
        export("[.stslist")
        S3method("[", "stslist")
        """,
    )

    with pytest.raises(LedgerError, match="both export and S3 method"):
        parse_namespace(path)


def test_vendored_namespaces_have_canonical_r_metadata_counts() -> None:
    traminer = parse_namespace(REPO / "SOURCES" / "TraMineR" / "NAMESPACE")
    extras = parse_namespace(REPO / "SOURCES" / "TraMineRextras" / "NAMESPACE")

    traminer_exports = {name for name, kind in traminer.items() if kind == "export"}
    traminer_s3 = {name for name, kind in traminer.items() if kind == "s3"}
    extras_exports = {name for name, kind in extras.items() if kind == "export"}
    extras_s3 = {name for name, kind in extras.items() if kind == "s3"}

    canonical_bracket_methods = {
        "[.seqelist",
        "[.stslist",
        "[.stslist.freq",
        "[.stslist.modst",
        "[.stslist.statd",
        "[.subseqelist",
    }
    incorrectly_quoted_methods = {
        '"[".seqelist',
        '"[".stslist',
        '"[".stslist.freq',
        '"[".stslist.modst',
        '"[".stslist.statd',
        '"[".subseqelist',
    }

    assert "obsolete" not in traminer
    assert "seqeordplot" not in extras
    assert canonical_bracket_methods <= traminer_s3
    assert incorrectly_quoted_methods.isdisjoint(traminer)
    assert (len(traminer_exports), len(traminer_s3), len(traminer)) == (125, 65, 190)
    assert (len(extras_exports), len(extras_s3), len(extras)) == (42, 8, 50)
    assert set(traminer).isdisjoint(extras)
    assert len(set(traminer) | set(extras)) == 240


def test_audit_generation_can_use_a_fixed_timestamp(tmp_path: Path) -> None:
    first = tmp_path / "first"
    second = tmp_path / "second"
    base = [
        sys.executable,
        str(REPO / "scripts" / "audit_traminex.py"),
        "--generated-at",
        "2026-07-13T00:00:00+00:00",
    ]

    subprocess.run(base + ["--output-dir", str(first)], cwd=REPO, check=True)
    subprocess.run(base + ["--output-dir", str(second)], cwd=REPO, check=True)

    assert {path.name: path.read_bytes() for path in first.iterdir()} == {
        path.name: path.read_bytes() for path in second.iterdir()
    }
