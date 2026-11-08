#!/usr/bin/env python3
"""Validate and generate the authoritative semantic R namespace ledger."""

from __future__ import annotations

import argparse
import ast
import csv
import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Iterable

if __package__:
    from ._r_namespace import NamespaceParseError as LedgerError
    from ._r_namespace import parse_namespace, promoted_internals, unrecorded_internal_promotions
else:
    from _r_namespace import NamespaceParseError as LedgerError
    from _r_namespace import parse_namespace, promoted_internals, unrecorded_internal_promotions

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DECISIONS = REPO / "docs" / "API_PARITY_DECISIONS.csv"
DEFAULT_MARKDOWN = REPO / "docs" / "API_PARITY_LEDGER.md"
DEFAULT_JSON = REPO / "reports" / "audit" / "traminex_semantic_parity.json"
MANIFEST = REPO / "traminex" / "api_manifest.py"
SOURCES = REPO / "SOURCES"
NAMESPACES = {
    "TraMineR": REPO / "SOURCES" / "TraMineR" / "NAMESPACE",
    "TraMineRextras": REPO / "SOURCES" / "TraMineRextras" / "NAMESPACE",
}
PACKAGE_ORDER = {name: index for index, name in enumerate(NAMESPACES)}

HEADER = [
    "r_package",
    "r_name",
    "r_kind",
    "python_name",
    "status",
    "scope",
    "option_coverage",
    "owner",
    "evidence",
    "notes",
]
# ``internal`` records an R object that traminex promotes to its public API
# although R defines it without exporting it. The NAMESPACE lists only exports
# and S3 methods, so before this value a promoted internal had no representable
# row and therefore no authoritative disposition (cycle-4b F-004). The
# membership is verified against the vendored source by
# ``scripts/_r_namespace.promoted_internals``.
ALLOWED_KINDS = {"export", "s3", "internal"}
ALLOWED_STATUSES = {
    "EXACT",
    "RE_EXPRESSED",
    "INTERNAL",
    "PLANNED",
    "DIVERGENCE",
    "SIBLING",
    "OUT",
}
ALLOWED_SCOPES = {"ANALYTICAL", "PLOTTING", "DISPLAY", "HELPER", "INTERNAL", "LEGACY"}
# The approved R0 documents require a closed option-coverage vocabulary but do
# not spell it out. These four exhaustive states make completeness explicit.
ALLOWED_OPTION_COVERAGE = {"FULL", "PARTIAL", "NONE", "N/A"}
ALLOWED_OWNERS = {"traminex", "weightedcluster", "py_seqplot", "R-only"}

# These are the historical boilerplate sentences that Phase 4 replaces with
# source-reviewed, row-specific decisions.  Rejecting the known placeholders is
# deliberately narrower than imposing a note-length or keyword heuristic.
PLACEHOLDER_NOTES = {
    "Current public implementation covers the reviewed R contract with behavioral test evidence.",
    "Python protocols result objects or explicit plotting functions replace one-for-one S3 dispatch.",
    "A reference matplotlib renderer re-expresses the base-R graphics contract over verified numerical series.",
    "R internal helper excluded as a public Python symbol by the approved design.",
    "Python mutable STSList metadata replaces the R setter.",
    "Uses a portable Graphviz text representation instead of R device-specific node images.",
    "Mutation on the Python event object replaces the R setter symbol.",
}

# The option audit found source-active R controls with no supported Python
# target.  Keeping this small, explicit registry prevents a known PARTIAL
# contract from being promoted to FULL merely because its existing Python
# arguments have tests.
KNOWN_PARTIAL_OPTIONS = {
    **{name: {"check.costs"} for name in ("seqcost", "seqsubm", "seqdist", "seqMD", "seqdistmc")},
    "plot.stslist": {"title", "cex.plot", "withlegend", "axes", "tlim", "ylas"},
    "disstree2dot": {
        "device",
        "device.args",
        "label.pos",
        "node.pos",
        "split.pos",
        "cex.main",
    },
    "disstree2dotp": {"image.data", "only.leaf", "image.fun", "cex.main", "cex.quality"},
    "disstreedisplay": {
        "image.data",
        "only.leaf",
        "image.fun",
        "cex.main",
        "legend.text",
        "show.tree",
    },
    "seqtree2dot": {"only.leaf", "diss", "cex.main", "xaxis"},
    "seqtreedisplay": {"only.leaf", "diss", "cex.main", "xaxis", "legend.text", "show.tree"},
    # Terms are named predictors: no formula interface, so no interaction terms
    # (TX-AUD-20261010-040).
    "dissmfacw": {"formula"},
}

# Every row here reaches genuine randomness in the reviewed public R contract.
# All mapped Python implementations use a NumPy stream or a deterministic
# plotting re-expression rather than R's global stream.  The validator checks
# the structural consequence (no unqualified EXACT claim); focused tests and
# row-specific notes retain the detailed contract.
DIFFERENT_RANDOM_STREAM_ROWS = {
    "dissassoc",
    "dissCompare",
    "dissmfacw",
    "dissrep",
    "dissrf",
    "disstree",
    "plot.seqpcplot",
    "plot.stslist",
    "seqBIC",
    "seqCompare",
    "seqgen",
    "seqgen.missing",
    "seqLRT",
    "seqpcplot",
    "seqplot.rf",
    "seqpolyads",
    "seqrep",
    "seqrf",
    "seqtree",
}


def current_namespace() -> dict[str, tuple[str, str]]:
    combined: dict[str, tuple[str, str]] = {}
    for package, path in NAMESPACES.items():
        for name, kind in parse_namespace(path).items():
            if name in combined:
                raise LedgerError(f"namespace name appears in both R packages: {name}")
            combined[name] = (package, kind)
    # R objects traminex promotes to its public API although R defines them without
    # exporting them. NAMESPACE lists only exports and S3 methods, so these have no
    # representable row without the ``internal`` kind, and seven public analytical
    # surfaces therefore carried no disposition at all (cycle-4b F-004). Membership
    # is verified against the vendored source, not merely asserted.
    for name, package in promoted_internals(SOURCES).items():
        if name in combined:
            raise LedgerError(f"promoted internal is actually exported: {name}")
        combined[name] = (package, "internal")
    return combined


def public_api() -> list[str]:
    module = ast.parse(MANIFEST.read_text(encoding="utf-8"), filename=str(MANIFEST))
    for node in module.body:
        value: ast.AST | None = None
        if isinstance(node, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == "PUBLIC_API" for target in node.targets
        ):
            value = node.value
        elif (
            isinstance(node, ast.AnnAssign)
            and isinstance(node.target, ast.Name)
            and node.target.id == "PUBLIC_API"
        ):
            value = node.value
        if isinstance(value, (ast.List, ast.Tuple)) and all(
            isinstance(item, ast.Constant) and isinstance(item.value, str) for item in value.elts
        ):
            return [str(item.value) for item in value.elts]
    raise LedgerError(f"could not parse PUBLIC_API from {MANIFEST.relative_to(REPO)}")


def load_decisions(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise LedgerError(f"missing decision source: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != HEADER:
            raise LedgerError(f"decision CSV header must be exactly {HEADER!r}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
    if not rows:
        raise LedgerError("decision source is empty")
    if any(None in row for row in rows):
        raise LedgerError("decision row has columns beyond the fixed schema")
    return rows


def target_phases(notes: str) -> list[str]:
    return re.findall(r"(?:^|;\s*)target=(R[1-4])(?=;|$)", notes)


def evidence_paths(value: str) -> list[Path]:
    return [Path(item.strip()) for item in value.split(";") if item.strip()]


def evidence_categories(row: dict[str, str], paths: list[Path]) -> dict[str, list[Path]]:
    """Classify evidence by what it can establish, independent of filename prose."""
    package_root = Path("SOURCES") / row["r_package"]
    source_roots = (package_root / "R", package_root / "man")
    return {
        "r_contract": [
            path
            for path in paths
            if any(path.is_relative_to(source_root) for source_root in source_roots)
        ],
        "behavioral_test": [
            path
            for path in paths
            if path.parts
            and path.parts[0] == "tests"
            and path.suffix == ".py"
            and path.name.startswith("test_")
        ],
        "r_generator": [
            path
            for path in paths
            if path.parts[:2] == ("scripts", "r_goldens") and path.suffix == ".R"
        ],
    }


def _validate_domain(row: dict[str, str], field: str, allowed: set[str]) -> None:
    if row[field] not in allowed:
        raise LedgerError(
            f"{row['r_name']}: invalid {field}={row[field]!r}; expected one of {sorted(allowed)}"
        )


def _validate_no_unrecorded_promotions(api: Iterable[str]) -> None:
    """Close the validator's one-directional gap (cycle-6 section 8.1).

    Everything below this function enforces R -> Python coverage: every R
    namespace name must have exactly one row. Nothing enforced the reverse, so a
    traminex public name with no row was invisible -- and 64 of the 226
    legitimately have none, being traminex-only extensions such as ``seqRLCP``,
    ``turbulence_proxy`` and the ``*Config`` dataclasses. A future promotion of an
    unexported R internal would therefore land silently in that bucket instead of
    being forced into ``PROMOTED_INTERNALS``.
    """
    unrecorded = unrecorded_internal_promotions(api, SOURCES)
    if unrecorded:
        described = "; ".join(
            f"{python_name} shadows the unexported R internal {r_name} ({packages})"
            for python_name, r_name, packages in sorted(unrecorded)
        )
        raise LedgerError(
            "public names promote unexported R internals without a record: "
            f"{described}. Add each to scripts/_r_namespace.py::PROMOTED_INTERNALS "
            "with its defining source, so it gains a ledger row with r_kind=internal."
        )


def validate_decisions(rows: list[dict[str, str]]) -> None:
    namespace = current_namespace()
    api = set(public_api())
    _validate_no_unrecorded_promotions(api)
    names = [row["r_name"] for row in rows]
    duplicates = sorted(name for name, count in Counter(names).items() if count != 1)
    if duplicates:
        raise LedgerError(f"duplicate decision rows: {duplicates}")
    missing = sorted(set(namespace) - set(names))
    unknown = sorted(set(names) - set(namespace))
    if missing or unknown:
        raise LedgerError(f"decision set mismatch: missing={missing}, unknown={unknown}")

    for row in rows:
        name = row["r_name"]
        expected_package, expected_kind = namespace[name]
        if row["r_package"] != expected_package or row["r_kind"] != expected_kind:
            raise LedgerError(
                f"{name}: expected package/kind {expected_package}/{expected_kind}, got "
                f"{row['r_package']}/{row['r_kind']}"
            )
        _validate_domain(row, "r_kind", ALLOWED_KINDS)
        _validate_domain(row, "status", ALLOWED_STATUSES)
        _validate_domain(row, "scope", ALLOWED_SCOPES)
        _validate_domain(row, "option_coverage", ALLOWED_OPTION_COVERAGE)
        _validate_domain(row, "owner", ALLOWED_OWNERS)

        paths = evidence_paths(row["evidence"])
        if not paths:
            raise LedgerError(f"{name}: evidence is required")
        for relative in paths:
            if relative.is_absolute() or ".." in relative.parts:
                raise LedgerError(
                    f"{name}: evidence must be a repository-relative path: {relative}"
                )
            if not (REPO / relative).is_file():
                raise LedgerError(f"{name}: evidence path does not exist: {relative}")

        categories = evidence_categories(row, paths)
        if not categories["r_contract"] and row["status"] not in {"SIBLING"}:
            expected = f"SOURCES/{row['r_package']}/R or man"
            raise LedgerError(
                f"{name}: evidence must cite reviewed {row['r_package']} source ({expected})"
            )

        if row["notes"] in PLACEHOLDER_NOTES:
            raise LedgerError(f"{name}: generic placeholder note must be replaced")

        executable = row["status"] in {"EXACT", "DIVERGENCE"} or (
            row["status"] == "RE_EXPRESSED" and bool(row["python_name"])
        )
        if executable and row["python_name"] not in api:
            raise LedgerError(
                f"{name}: {row['status']} mapping is not current public API: "
                f"{row['python_name']}"
            )
        if executable and not categories["behavioral_test"]:
            raise LedgerError(f"{name}: {row['status']} requires focused behavioral test evidence")
        if executable and row["python_name"]:
            # A cited test must at least name the mapped function (2026-09-22
            # audit S10-05: plot.stslist cited a test that never called seqplot).
            mentions = re.compile(rf"(?<!\w){re.escape(row['python_name'])}(?!\w)")
            if not any(
                mentions.search((REPO / path).read_text(encoding="utf-8"))
                for path in categories["behavioral_test"]
            ):
                raise LedgerError(
                    f"{name}: no cited behavioral test ever mentions {row['python_name']}"
                )

        if row["status"] == "EXACT":
            if row["option_coverage"] != "FULL":
                raise LedgerError(f"{name}: EXACT requires option_coverage=FULL")
            if not categories["r_contract"] or not categories["behavioral_test"]:
                raise LedgerError(
                    f"{name}: EXACT requires reviewed public-R source or Rd and "
                    "a focused behavioral test"
                )
            if row["scope"] == "ANALYTICAL" and not categories["r_generator"]:
                raise LedgerError(
                    f"{name}: an EXACT analytical row requires the R golden generator "
                    "(scripts/r_goldens/*.R) that records R's result"
                )

        if row["option_coverage"] == "FULL" and executable:
            if not categories["behavioral_test"]:
                raise LedgerError(
                    f"{name}: FULL requires option/default/validation behavioral evidence"
                )
        if name in KNOWN_PARTIAL_OPTIONS and row["option_coverage"] == "FULL":
            unsupported = ", ".join(sorted(KNOWN_PARTIAL_OPTIONS[name]))
            raise LedgerError(
                f"{name}: FULL cannot hide known unsupported effective options: {unsupported}"
            )

        if name in DIFFERENT_RANDOM_STREAM_ROWS:
            if row["status"] == "EXACT":
                raise LedgerError(
                    f"{name}: reviewed random stream differs from R and cannot be EXACT"
                )
            if not categories["behavioral_test"]:
                raise LedgerError(
                    f"{name}: stochastic disposition requires focused behavioral test evidence"
                )

        if row["status"] == "PLANNED" and row["scope"] == "ANALYTICAL":
            phases = target_phases(row["notes"])
            if len(phases) != 1:
                raise LedgerError(
                    f"{name}: PLANNED analytical row requires exactly one target=R1..R4"
                )
            if row["owner"] == "traminex" and not row["python_name"]:
                raise LedgerError(f"{name}: planned traminex analytical row requires python_name")


def _count(rows: Iterable[dict[str, str]], field: str) -> dict[str, int]:
    counts = Counter(row[field] for row in rows)
    return dict(sorted(counts.items()))


def summaries(rows: list[dict[str, str]]) -> dict[str, dict[str, int]]:
    phase_counts = Counter(
        phase for row in rows if row["status"] == "PLANNED" for phase in target_phases(row["notes"])
    )
    return {
        "r_package": _count(rows, "r_package"),
        "r_kind": _count(rows, "r_kind"),
        "status": _count(rows, "status"),
        "scope": _count(rows, "scope"),
        "option_coverage": _count(rows, "option_coverage"),
        "owner": _count(rows, "owner"),
        "target_phase": dict(sorted(phase_counts.items())),
    }


def target_public_api(rows: list[dict[str, str]]) -> dict[str, object]:
    current = public_api()
    current_set = set(current)
    planned = sorted(
        {
            row["python_name"]
            for row in rows
            if row["status"] == "PLANNED"
            and row["scope"] == "ANALYTICAL"
            and row["owner"] == "traminex"
        }
    )
    added = [name for name in planned if name not in current_set]
    target = list(current)
    target.extend(added)
    return {
        "current_count": len(current),
        "planned_analytical_names": planned,
        "added_names": added,
        "added_count": len(target) - len(current),
        "target_count": len(target),
        "target_names": target,
    }


def _sorted_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return sorted(rows, key=lambda row: (PACKAGE_ORDER[row["r_package"]], row["r_name"].lower()))


def render_json(rows: list[dict[str, str]]) -> str:
    ordered = _sorted_rows(rows)
    payload = {
        "schema_version": 1,
        "source": "docs/API_PARITY_DECISIONS.csv",
        "decision_count": len(ordered),
        "summaries": summaries(ordered),
        "target_public_api": target_public_api(ordered),
        "decisions": ordered,
    }
    return json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=False) + "\n"


def _markdown_table(counts: dict[str, int], label: str) -> list[str]:
    lines = [f"### {label}", "", f"| {label} | Count |", "| --- | ---: |"]
    lines.extend(f"| `{name}` | {count} |" for name, count in counts.items())
    lines.append("")
    return lines


def _md_cell(value: str) -> str:
    return value.replace("|", "\\|").replace("\n", " ")


def render_markdown(rows: list[dict[str, str]]) -> str:
    ordered = _sorted_rows(rows)
    summary = summaries(ordered)
    target = target_public_api(ordered)
    # State the denominator the row set actually covers, derived from the rows
    # themselves so the prose cannot drift from the data. "Every namespace name"
    # matched neither reading: exported + registered S3 is 240 (167 + 73), while
    # every binding in the two namespaces is 358 + 66 (TraMineR 2.2-14 and
    # TraMineRextras 0.6.9; cycle-6 F-008, recounted for TX-AUD-20261019-057).
    kinds = summary["r_kind"]
    exports = kinds.get("export", 0)
    s3_methods = kinds.get("s3", 0)
    internals = kinds.get("internal", 0)
    lines = [
        "# Semantic parity ledger",
        "",
        "This file is generated from `docs/API_PARITY_DECISIONS.csv`. Do not edit it by hand.",
        "The ledger accounts exactly once for every name exported by TraMineR 2.2-14 or",
        f"TraMineRextras 0.6.9 ({exports}) and every S3 method their `NAMESPACE` files",
        f"register ({s3_methods}), plus the {internals} verified unexported internals that",
        f"traminex promotes to public names — {exports + s3_methods + internals} rows in all.",
        "`docs/API_PARITY_TABLE.md` publishes the first two groups as",
        f"`R_total_unique={exports + s3_methods}`; the difference between the two figures is",
        "exactly those promoted internals. It does *not* enumerate every binding in the two",
        "namespaces (`ls(asNamespace(...), all.names = TRUE)`), which is far larger.",
        "Semantic decisions are source-reviewed rather than inferred from name similarity.",
        "",
        "## Derived totals",
        "",
        f"- Decisions: {len(ordered)}",
        f"- Current Python public API: {target['current_count']}",
        f"- Ledger-derived target public API: {target['target_count']}",
        f"- New planned analytical public names: {target['added_count']}",
        "",
        "## Summaries",
        "",
    ]
    labels = {
        "r_package": "R package",
        "r_kind": "R kind",
        "status": "Status",
        "scope": "Scope",
        "option_coverage": "Option coverage",
        "owner": "Owner",
        "target_phase": "Target phase",
    }
    for key in (
        "r_package",
        "r_kind",
        "status",
        "scope",
        "option_coverage",
        "owner",
        "target_phase",
    ):
        lines.extend(_markdown_table(summary[key], labels[key]))

    lines.extend(
        [
            "## Planned analytical contracts",
            "",
            ", ".join(f"`{name}`" for name in target["planned_analytical_names"]) or "(none)",
            "",
            "## Target public API additions",
            "",
            ", ".join(f"`{name}`" for name in target["added_names"]) or "(none)",
            "",
            "## Decisions",
            "",
            "| R package | R name | Kind | Python name | Status | Scope | Options | Owner | Evidence | Notes |",
            "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
        ]
    )
    for row in ordered:
        cells = [_md_cell(row[column]) for column in HEADER]
        lines.append("| " + " | ".join(cells) + " |")
    lines.append("")
    return "\n".join(lines)


def _write_or_check(path: Path, content: str, *, check: bool) -> bool:
    if check:
        return path.is_file() and path.read_text(encoding="utf-8") == content
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8", newline="\n")
    return True


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--decisions", type=Path, default=DEFAULT_DECISIONS)
    parser.add_argument("--markdown-output", type=Path, default=DEFAULT_MARKDOWN)
    parser.add_argument("--json-output", type=Path, default=DEFAULT_JSON)
    parser.add_argument(
        "--check", action="store_true", help="check tracked outputs without writing"
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        rows = load_decisions(args.decisions)
        validate_decisions(rows)
        markdown = render_markdown(rows)
        json_text = render_json(rows)
    except LedgerError as exc:
        print(f"Semantic parity ledger validation FAILED: {exc}", file=sys.stderr)
        return 1

    results = {
        args.markdown_output: _write_or_check(args.markdown_output, markdown, check=args.check),
        args.json_output: _write_or_check(args.json_output, json_text, check=args.check),
    }
    stale = [path for path, matches in results.items() if not matches]
    if stale:
        print("Semantic parity ledger generation check FAILED; stale or missing outputs:")
        for path in stale:
            print(f"- {path}")
        return 1
    action = "checked" if args.check else "generated"
    print(f"Semantic parity ledger {action}: {len(rows)} decisions")
    print(json.dumps(summaries(rows), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
