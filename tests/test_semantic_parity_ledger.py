from __future__ import annotations

import ast
import copy
import csv
import importlib
import json
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pytest

from scripts import generate_semantic_parity_ledger as semantic_ledger
from scripts._r_namespace import parse_namespace, promoted_internals
from traminex.api_manifest import PUBLIC_API_SET

REPO = Path(__file__).resolve().parents[1]
DECISIONS = REPO / "docs" / "API_PARITY_DECISIONS.csv"
LEDGER_MD = REPO / "docs" / "API_PARITY_LEDGER.md"
LEDGER_JSON = REPO / "reports" / "audit" / "traminex_semantic_parity.json"
GENERATOR = REPO / "scripts" / "generate_semantic_parity_ledger.py"

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
ALLOWED_OPTION_COVERAGE = {"FULL", "PARTIAL", "NONE", "N/A"}
ALLOWED_OWNERS = {"traminex", "weightedcluster", "py_seqplot", "R-only"}
PACKAGES = {
    "TraMineR": REPO / "SOURCES" / "TraMineR" / "NAMESPACE",
    "TraMineRextras": REPO / "SOURCES" / "TraMineRextras" / "NAMESPACE",
}
R3_ROWS = {
    "seqecreate": "seqecreate",
    "createdatadiscrete": "createdatadiscrete",
    "dissCompare": "dissCompare",
    "group.p": "group_p",
    "seqBIC": "seqBIC",
    "seqCompare": "seqCompare",
    "seqerulesdisc": "seqerulesdisc",
    "seqgen.missing": "seqgen_missing",
    "seqLRT": "seqLRT",
    "seqsamm": "seqsamm",
    "seqsammeha": "seqsammeha",
    "seqsammseq": "seqsammseq",
    "sorti": "sorti",
    "sortv": "sortv",
}
R3_ADDITIONS = {
    "createdatadiscrete",
    "group_p",
    "seqerulesdisc",
    "seqgen_missing",
    "seqsamm",
    "seqsammeha",
    "seqsammseq",
    "sorti",
    "sortv",
}
R4_ROWS = {
    "plot.stslist": "seqplot",
    "plot.SAMM": "plot_SAMM",
    "plot.seqimplic": "plot_seqimplic",
    "seqplot.rf": "seqplot_rf",
    "seqplot.tentrop": "seqplot_tentrop",
    "seqplot.tentrop.m": "seqplot_tentrop_m",
    "seqsplot": "seqsplot",
}
R4_ADDITIONS = {
    "plot_SAMM",
    "plot_seqimplic",
    "seqplot_rf",
    "seqplot_tentrop",
    "seqplot_tentrop_m",
    "seqsplot",
}


def _current_namespace() -> dict[str, tuple[str, str]]:
    combined: dict[str, tuple[str, str]] = {}
    for package, path in PACKAGES.items():
        for name, kind in parse_namespace(path).items():
            assert name not in combined, f"namespace name is exported by both packages: {name}"
            combined[name] = (package, kind)
    # R objects traminex promotes to its public API although R defines them
    # without exporting them. NAMESPACE cannot describe them, so they are carried
    # with r_kind='internal' and verified against the vendored source (F-004).
    for name, package in promoted_internals(REPO / "SOURCES").items():
        assert name not in combined, f"promoted internal is actually exported: {name}"
        combined[name] = (package, "internal")
    return combined


def _load_decisions() -> list[dict[str, str]]:
    assert DECISIONS.is_file(), f"missing authoritative decision source: {DECISIONS}"
    with DECISIONS.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames == HEADER
        rows = list(reader)
    assert rows, "semantic parity decision source is empty"
    assert all(None not in row for row in rows), "decision row has columns beyond the fixed schema"
    return rows


def _evidence_paths(value: str) -> list[Path]:
    return [REPO / item.strip() for item in value.split(";") if item.strip()]


def _target_phases(notes: str) -> list[str]:
    return re.findall(r"(?:^|;\s*)target=(R[1-4])(?=;|$)", notes)


def _rows_by_name() -> dict[str, dict[str, str]]:
    return {row["r_name"]: row for row in _load_decisions()}


def _copied_decisions() -> list[dict[str, str]]:
    return copy.deepcopy(_load_decisions())


def _copied_row(rows: list[dict[str, str]], name: str) -> dict[str, str]:
    return next(row for row in rows if row["r_name"] == name)


def test_validator_rejects_generic_exact_note() -> None:
    rows = _copied_decisions()
    row = _copied_row(rows, "cpal")
    row["notes"] = (
        "Current public implementation covers the reviewed R contract "
        "with behavioral test evidence."
    )

    with pytest.raises(semantic_ledger.LedgerError, match="placeholder"):
        semantic_ledger.validate_decisions(rows)


def test_validator_rejects_exact_without_focused_behavioral_test() -> None:
    rows = _copied_decisions()
    row = _copied_row(rows, "cpal")
    row["evidence"] = "SOURCES/TraMineR/R/cpal.R;docs/API_PARITY.md"
    row["notes"] = (
        "Palette lookup preserves the source state order, missing-state color, "
        "and validation boundary."
    )

    with pytest.raises(semantic_ledger.LedgerError, match="focused behavioral test"):
        semantic_ledger.validate_decisions(rows)


def test_validator_rejects_r_source_from_the_wrong_package() -> None:
    rows = _copied_decisions()
    row = _copied_row(rows, "cpal")
    row["evidence"] = "SOURCES/TraMineRextras/R/rowmode.R;tests/test_palette_legend.py"
    row["notes"] = (
        "Palette lookup preserves the source state order, missing-state color, "
        "and validation boundary."
    )

    with pytest.raises(semantic_ledger.LedgerError, match="TraMineR source"):
        semantic_ledger.validate_decisions(rows)


def test_validator_rejects_stochastic_exact_row_without_stream_parity() -> None:
    rows = _copied_decisions()
    # dissassoc: a stochastic row with FULL coverage (dissmfacw, the earlier example, is
    # PARTIAL since TX-AUD-20261010-040 and would fail the coverage rule first).
    row = _copied_row(rows, "dissassoc")
    row["status"] = "EXACT"
    row["notes"] = "A fixed NumPy seed repeats within Python."

    with pytest.raises(semantic_ledger.LedgerError, match="random stream"):
        semantic_ledger.validate_decisions(rows)


def test_validator_rejects_known_partial_plot_as_full() -> None:
    rows = _copied_decisions()
    row = _copied_row(rows, "plot.stslist")
    row["option_coverage"] = "FULL"

    with pytest.raises(semantic_ledger.LedgerError, match="unsupported effective options"):
        semantic_ledger.validate_decisions(rows)


def test_validator_rejects_executable_difference_without_behavioral_test() -> None:
    rows = _copied_decisions()
    row = _copied_row(rows, "seqhasmiss")
    row["evidence"] = "SOURCES/TraMineR/R/seqhasmiss.R"

    with pytest.raises(semantic_ledger.LedgerError, match="focused behavioral test"):
        semantic_ledger.validate_decisions(rows)


def test_validator_rejects_exact_analytical_row_without_an_r_generator() -> None:
    """An EXACT analytical claim needs the generator that records R's result.

    The validator classified generator evidence but never required it
    (2026-09-22 audit S10-05).
    """
    rows = _copied_decisions()
    # Any EXACT analytical row that cites a generator (seqdur, the earlier
    # example, became a DIVERGENCE on 2026-10-17, TX-AUD-20261010-055).
    row = next(
        candidate
        for candidate in rows
        if candidate["status"] == "EXACT"
        and candidate["scope"] == "ANALYTICAL"
        and "scripts/r_goldens/" in candidate["evidence"]
    )
    row["evidence"] = ";".join(
        path for path in row["evidence"].split(";") if not path.startswith("scripts/r_goldens/")
    )

    with pytest.raises(semantic_ledger.LedgerError, match="R golden generator"):
        semantic_ledger.validate_decisions(rows)


def test_validator_rejects_cited_tests_that_never_name_the_mapped_function() -> None:
    """plot.stslist cited, as its only test for ``seqplot``, one that never calls it."""
    rows = _copied_decisions()
    row = _copied_row(rows, "cpal")
    row["evidence"] = "SOURCES/TraMineR/R/cpal.R;tests/test_seqe_pipeline.py"

    with pytest.raises(semantic_ledger.LedgerError, match="no cited behavioral test ever mentions"):
        semantic_ledger.validate_decisions(rows)


def test_authoritative_decisions_pass_the_direct_validator() -> None:
    semantic_ledger.validate_decisions(_copied_decisions())


def test_decision_source_accounts_for_current_namespace() -> None:
    rows = _load_decisions()
    current = _current_namespace()
    names = [row["r_name"] for row in rows]

    # 240 NAMESPACE names + 7 promoted internals that NAMESPACE cannot describe.
    assert sum(kind != "internal" for _, kind in current.values()) == 240
    assert sum(kind == "internal" for _, kind in current.values()) == 7
    assert len(current) == 247
    assert len(rows) == 247
    assert set(names) == set(current)
    assert not [name for name, count in Counter(names).items() if count != 1]
    for row in rows:
        assert (row["r_package"], row["r_kind"]) == current[row["r_name"]]


def test_every_public_name_shadowing_an_r_internal_is_classified() -> None:
    """cycle-6 section 8.1: close the validator's one-directional gap.

    The validator enforces R -> Python coverage only, so a traminex public name
    with no ledger row is invisible -- and 64 of the 226 legitimately have none,
    being traminex-only extensions. A future promotion of an unexported R
    internal would land silently in that bucket instead of being forced into
    PROMOTED_INTERNALS.

    Scanning the 226 public names against the vendored R sources under exact and
    dot-for-underscore spellings currently yields 14 matches. All 14 must be
    accounted for and none flagged as new: 7 recorded promotions, 6 S3 methods
    the NAMESPACE does register, and ``seqedist``, which TraMineRextras exports
    while TraMineR merely has an internal of the same name.
    """
    from scripts._r_namespace import (
        unrecorded_internal_promotions,
        vendored_object_definitions,
    )

    sources = REPO / "SOURCES"
    definitions = vendored_object_definitions(sources)
    namespace = {}
    for package, path in PACKAGES.items():
        for name, kind in parse_namespace(path).items():
            namespace.setdefault(name, (package, kind))

    matched = {
        name: candidate
        for name in PUBLIC_API_SET
        for candidate in (name, name.replace("_", "."))
        if candidate in definitions
    }
    promotions = {name for name in matched if matched[name] in promoted_internals(sources)}
    registered = {name for name in matched if matched[name] in namespace}

    assert promotions == {
        "STS_to_SPELL",
        "STS_to_SRS",
        "STS_to_TSE",
        "entropy",
        "is_subseqelist",
        "seqasnum",
        "seqerules",
    }
    assert {name for name in registered if namespace[matched[name]][1] == "s3"} >= {
        "plot_SAMM",
        "plot_eseq",
        "plot_seqelist",
        "plot_seqimplic",
        "plot_subseqelist",
        "plot_subseqelistchisq",
    }
    assert namespace["seqedist"] == ("TraMineRextras", "export")
    assert unrecorded_internal_promotions(sorted(PUBLIC_API_SET), sources) == []


def test_an_unrecorded_promotion_is_detected() -> None:
    """The check must discriminate, not merely return an empty list today."""
    from scripts._r_namespace import PROMOTED_INTERNALS, unrecorded_internal_promotions

    sources = REPO / "SOURCES"
    # `seqasnum` is a genuine unexported TraMineR internal; pretend it is not
    # recorded and confirm the scan names it.
    saved = dict(PROMOTED_INTERNALS)
    PROMOTED_INTERNALS.pop("seqasnum")
    try:
        found = unrecorded_internal_promotions(["seqasnum"], sources)
    finally:
        PROMOTED_INTERNALS.clear()
        PROMOTED_INTERNALS.update(saved)

    assert [entry[:2] for entry in found] == [("seqasnum", "seqasnum")]
    assert unrecorded_internal_promotions(["seqasnum"], sources) == []


def test_generated_completeness_sentence_states_the_denominator_it_covers() -> None:
    """cycle-6 F-008: the ledger's prose must describe its own row set.

    It claimed to account for "every current TraMineR and TraMineRextras
    namespace name exactly once" above its decision count. Neither reading of
    that phrase matched: exported names plus registered S3 methods is 240 under
    the 2.2-14 oracle (it was 243 under 2.2-13, which this note first described), while
    every binding in the two namespaces (``ls(asNamespace, all.names = TRUE)``)
    is far larger. The 7 surplus rows are the deliberate promoted internals --
    sound design, but not what the sentence said. The row set itself was correct
    and is not changed here; only the prose was wrong, and no gate compared
    generated prose against the data it summarises.
    """
    rows = _load_decisions()
    kinds = Counter(row["r_kind"] for row in rows)
    text = LEDGER_MD.read_text(encoding="utf-8")
    header = text.split("## Derived totals", 1)[0]

    for count in (kinds["export"], kinds["s3"], kinds["internal"], len(rows)):
        assert str(count) in header, f"the ledger preamble does not state {count}"
    # The other published R denominator must be reconciled, not merely coexist.
    assert f"R_total_unique={kinds['export'] + kinds['s3']}" in header
    assert "promote" in header.lower()
    assert "every current TraMineR 2.2-14 and TraMineRextras 0.6.9" not in header


def test_parity_table_reconciles_its_r_denominator_with_the_ledger() -> None:
    table = (REPO / "docs" / "API_PARITY_TABLE.md").read_text(encoding="utf-8")
    rows = _load_decisions()
    kinds = Counter(row["r_kind"] for row in rows)
    unique = kinds["export"] + kinds["s3"]

    assert f"R_total_unique={unique}" in table
    reconciliation = [
        line for line in table.splitlines() if line.startswith(f"`R_total_unique={unique}`")
    ]
    assert len(reconciliation) == 1, "the R denominator carries no reconciliation line"
    assert "PROMOTED_INTERNALS" in reconciliation[0]


def test_decision_rows_have_valid_domains_mappings_and_evidence() -> None:
    rows = _load_decisions()

    assert {row["r_kind"] for row in rows} <= ALLOWED_KINDS
    assert {row["status"] for row in rows} <= ALLOWED_STATUSES
    assert {row["scope"] for row in rows} <= ALLOWED_SCOPES
    assert {row["option_coverage"] for row in rows} <= ALLOWED_OPTION_COVERAGE
    assert {row["owner"] for row in rows} <= ALLOWED_OWNERS

    for row in rows:
        evidence = _evidence_paths(row["evidence"])
        assert evidence, f"{row['r_name']}: evidence is required"
        assert all(path.is_file() for path in evidence), (
            f"{row['r_name']}: missing evidence path(s): "
            f"{[str(path.relative_to(REPO)) for path in evidence if not path.is_file()]}"
        )
        if row["status"] == "EXACT":
            assert row["option_coverage"] == "FULL"
            assert row["python_name"] in PUBLIC_API_SET
        if row["status"] == "DIVERGENCE":
            assert row["python_name"] in PUBLIC_API_SET
        if row["status"] == "RE_EXPRESSED" and row["python_name"]:
            assert row["python_name"] in PUBLIC_API_SET
        if row["status"] == "PLANNED" and row["scope"] == "ANALYTICAL":
            assert len(_target_phases(row["notes"])) == 1


def test_known_scope_contradictions_stay_explicitly_resolved() -> None:
    rows = _rows_by_name()

    assert "obsolete" not in rows
    assert "seqeordplot" not in rows
    assert "seqeordplot" in PUBLIC_API_SET

    # TraMineR's seqerules is public in Python although R defines it without
    # exporting it. It previously had no row at all -- "intentionally public, no
    # disposition" left seven such surfaces unrecorded (cycle-4b F-004) -- so it is
    # now carried with r_kind='internal'. The TraMineRextras seqerulesdisc API
    # remains distinct.
    assert "seqerules" in PUBLIC_API_SET
    assert rows["seqerules"]["r_kind"] == "internal"
    assert rows["seqerules"]["python_name"] == "seqerules"
    seqerulesdisc = rows["seqerulesdisc"]
    assert (seqerulesdisc["python_name"], seqerulesdisc["status"]) == (
        "seqerulesdisc",
        "DIVERGENCE",
    )
    assert _target_phases(seqerulesdisc["notes"]) == []
    assert "NumPy/SciPy IRLS" in seqerulesdisc["notes"]
    assert "supvars" in seqerulesdisc["notes"]
    assert "empty subsequence" in seqerulesdisc["notes"]

    seqsamm = rows["seqsamm"]
    assert "fractional sublength" in seqsamm["notes"]
    assert "short covariate" in seqsamm["notes"]
    for name in {"seqecreate", "dissCompare", "group.p", "seqsammeha"}:
        assert "collation" in rows[name]["notes"] or "category order" in rows[name]["notes"]

    for name in ("seqcost", "seqsubm", "seqdist", "seqMD", "seqdistmc"):
        assert (rows[name]["status"], rows[name]["option_coverage"]) == ("DIVERGENCE", "PARTIAL")
        assert "check.costs" in rows[name]["notes"]
        # TX-AUD-20261011-041: R's default is accepted, FALSE refused with the reason.
        assert "check_costs=True, R's default, is accepted" in rows[name]["notes"]
        assert "check_costs=False is refused" in rows[name]["notes"]
    assert (
        rows["seqeapplysub"]["status"],
        rows["seqeapplysub"]["option_coverage"],
    ) == (
        "DIVERGENCE",
        "FULL",
    )
    assert "CWIN" in rows["seqeapplysub"]["notes"]
    assert (
        rows["seqtrate"]["status"],
        rows["seqtrate"]["option_coverage"],
        _target_phases(rows["seqtrate"]["notes"]),
    ) == ("DIVERGENCE", "FULL", [])
    assert "time-varying" in rows["seqtrate"]["notes"]
    assert (rows["alphabet"]["status"], rows["alphabet"]["option_coverage"]) == (
        "RE_EXPRESSED",
        "N/A",
    )
    assert "SeqElist" in rows["alphabet"]["notes"]
    assert "PSTf" in rows["alphabet"]["notes"]

    assert all(_target_phases(rows[name]["notes"]) == [] for name in R3_ROWS)
    assert not any(_target_phases(row["notes"]) for row in rows.values())
    assert "isolated NumPy generator" in rows["dissrf"]["notes"]
    assert "wcmdscale-equivalent weighted default MDS" in rows["dissrf"]["notes"]
    assert "wcmdscale-equivalent weighted default MDS" in rows["seqrf"]["notes"]
    assert "corrected weighted boundary" in rows["seqdomassoc"]["notes"]
    assert "state-by-state-by-time" in rows["seqMD"]["notes"]
    assert "overwritten-seed defect" in rows["seqpolyads"]["notes"]
    assert "missing-state survival" in rows["seqsurv"]["notes"]
    assert "intentionally remains quiet" in rows["seqistatd"]["notes"]
    assert "actionable ValueError" in rows["seqient"]["notes"]
    assert "Portable image_fun" in rows["disstree2dot"]["notes"]
    assert rows["disstree2dot"]["option_coverage"] == "PARTIAL"
    assert rows["seqrep"]["status"] == "DIVERGENCE"
    assert "random criterion" in rows["seqrep"]["notes"]
    assert rows["seqhasmiss"]["option_coverage"] == "PARTIAL"
    assert "separate has.miss, has.nr, and has.void vectors" in rows["seqhasmiss"]["notes"]
    assert rows["rowmode"]["status"] == "DIVERGENCE"
    assert "2-D inputs that R flattens" in rows["rowmode"]["notes"]

    expected = {
        "seqlegend": ("state_legend", "RE_EXPRESSED", []),
        "sorti": ("sorti", "DIVERGENCE", []),
        "sortv": ("sortv", "DIVERGENCE", []),
        "seqsplot": ("seqsplot", "RE_EXPRESSED", []),
        "seqplot.rf": ("seqplot_rf", "RE_EXPRESSED", []),
        "seqplot.tentrop": ("seqplot_tentrop", "RE_EXPRESSED", []),
        "seqplot.tentrop.m": ("seqplot_tentrop_m", "RE_EXPRESSED", []),
        "seqintegration": ("seqintegr", "RE_EXPRESSED", []),
        "seqprecorr": ("seqidegrad", "RE_EXPRESSED", []),
        "recodef": ("seqrecode", "RE_EXPRESSED", []),
        "seqfcheck": ("", "INTERNAL", []),
    }
    for r_name, (python_name, status, phases) in expected.items():
        assert (rows[r_name]["python_name"], rows[r_name]["status"]) == (python_name, status)
        assert _target_phases(rows[r_name]["notes"]) == phases


def test_seqdist_deliberate_python_method_default_is_explicitly_decided() -> None:
    row = _rows_by_name()["seqdist"]

    assert row["status"] == "DIVERGENCE"
    assert row["owner"] == "traminex"
    assert "R requires explicit `method`" in row["notes"]
    assert 'Python defaults `method="OM"`' in row["notes"]
    assert "tests/test_seqdist_wrapper.py" in row["evidence"]


def test_summaries_and_phase_schedule_are_derived_from_decisions() -> None:
    rows = _load_decisions()
    summary = semantic_ledger.summaries(rows)

    for field in ("r_package", "r_kind", "status", "scope", "option_coverage", "owner"):
        assert summary[field] == dict(sorted(Counter(row[field] for row in rows).items()))
        assert sum(summary[field].values()) == len(rows)
    assert summary["target_phase"] == {}
    assert not any(row["status"] == "PLANNED" for row in rows)
    partial = {row["r_name"] for row in rows if row["option_coverage"] == "PARTIAL"}
    assert set(semantic_ledger.KNOWN_PARTIAL_OPTIONS) <= partial
    # STS_to_TSE is PARTIAL because R's id and tevent formals have no Python
    # counterpart; see its decision row (cycle-4b F-004).
    assert partial - set(semantic_ledger.KNOWN_PARTIAL_OPTIONS) == {
        "seqhasmiss",
        "STS_to_TSE",
    }


def test_r3_extras_scope_is_closed_on_the_220_name_surface() -> None:
    rows = _rows_by_name()

    assert {name: rows[name]["python_name"] for name in R3_ROWS} == R3_ROWS
    assert not any("target=R3" in row["notes"] for row in rows.values())
    assert not {
        row["r_name"]
        for row in rows.values()
        if row["status"] == "PLANNED" and row["scope"] == "ANALYTICAL"
    }
    r3_baseline = PUBLIC_API_SET - R3_ADDITIONS - R4_ADDITIONS
    assert len(r3_baseline) == 211
    assert R3_ADDITIONS <= PUBLIC_API_SET
    assert len(r3_baseline.union(R3_ADDITIONS)) == 220


def test_r4_visualization_scope_is_closed_on_the_220_to_226_transition() -> None:
    rows = _rows_by_name()

    assert {name: rows[name]["python_name"] for name in R4_ROWS} == R4_ROWS
    assert all(rows[name]["scope"] == "PLOTTING" for name in R4_ROWS)
    assert not any("target=R4" in row["notes"] for row in rows.values())
    assert (rows["plot.stslist"]["status"], rows["plot.stslist"]["option_coverage"]) == (
        "DIVERGENCE",
        "PARTIAL",
    )
    assert "defaults to all rows" in rows["plot.stslist"]["notes"]
    assert "NumPy for random sampling" in rows["plot.stslist"]["notes"]
    assert "missing-state extension" in rows["seqsplot"]["notes"]

    assert R4_ADDITIONS <= PUBLIC_API_SET
    assert len(PUBLIC_API_SET - R4_ADDITIONS) == 220
    assert len(PUBLIC_API_SET) == 226
    assert "seqplot" in PUBLIC_API_SET

    payload = json.loads(LEDGER_JSON.read_text(encoding="utf-8"))
    target = payload["target_public_api"]
    assert target["current_count"] == target["target_count"] == len(PUBLIC_API_SET)
    assert target["planned_analytical_names"] == []
    assert target["added_names"] == []
    assert payload["summaries"]["target_phase"] == {}


def test_manifest_runtime_root_and_stub_share_the_226_name_surface() -> None:
    module = importlib.import_module("traminex")
    stub = ast.parse((REPO / "stubs" / "traminex" / "__init__.pyi").read_text(encoding="utf-8"))
    stub_all: list[str] | None = None
    for node in stub.body:
        if not isinstance(node, ast.Assign):
            continue
        if not any(
            isinstance(target, ast.Name) and target.id == "__all__" for target in node.targets
        ):
            continue
        if isinstance(node.value, (ast.List, ast.Tuple)):
            stub_all = [
                str(item.value)
                for item in node.value.elts
                if isinstance(item, ast.Constant) and isinstance(item.value, str)
            ]
        break

    assert len(PUBLIC_API_SET) == 226
    assert list(module.__all__) == list(semantic_ledger.public_api())
    assert set(module.__all__) == PUBLIC_API_SET
    assert stub_all == list(module.__all__)
    assert all(getattr(module, name) is not None for name in module.__all__)


def test_generated_semantic_ledger_is_tracked_and_deterministic(tmp_path: Path) -> None:
    _load_decisions()
    assert GENERATOR.is_file()
    assert LEDGER_MD.is_file()
    assert LEDGER_JSON.is_file()

    first_md = tmp_path / "first.md"
    first_json = tmp_path / "first.json"
    second_md = tmp_path / "second.md"
    second_json = tmp_path / "second.json"
    base = [sys.executable, str(GENERATOR), "--decisions", str(DECISIONS)]
    subprocess.run(
        base + ["--markdown-output", str(first_md), "--json-output", str(first_json)],
        cwd=REPO,
        check=True,
    )
    subprocess.run(
        base + ["--markdown-output", str(second_md), "--json-output", str(second_json)],
        cwd=REPO,
        check=True,
    )

    assert first_md.read_bytes() == second_md.read_bytes() == LEDGER_MD.read_bytes()
    assert first_json.read_bytes() == second_json.read_bytes() == LEDGER_JSON.read_bytes()
    subprocess.run(
        [sys.executable, str(GENERATOR), "--check"],
        cwd=REPO,
        check=True,
    )
    payload = json.loads(first_json.read_text(encoding="utf-8"))
    # 240 NAMESPACE names + 7 promoted internals (cycle-4b F-004).
    assert payload["decision_count"] == 247
    assert len(payload["decisions"]) == 247
    target = payload["target_public_api"]
    assert target["current_count"] == len(PUBLIC_API_SET)
    assert target["target_count"] == len(PUBLIC_API_SET)
    assert target["planned_analytical_names"] == []
    assert target["added_count"] == 0
    assert target["added_names"] == []
    assert (
        len((PUBLIC_API_SET - R3_ADDITIONS - R4_ADDITIONS) - {"seqcost", "seqeapplysub", "seqsep"})
        == 208
    )
    assert {"seqcost", "seqeapplysub", "seqsep"} <= PUBLIC_API_SET
    assert set(target["target_names"]) == PUBLIC_API_SET
    assert payload["summaries"]["target_phase"] == {}
    assert set(payload["summaries"]) == {
        "r_package",
        "r_kind",
        "status",
        "scope",
        "option_coverage",
        "owner",
        "target_phase",
    }
    assert not any("timestamp" in key.lower() or "generated_at" in key.lower() for key in payload)
