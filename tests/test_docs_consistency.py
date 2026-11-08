import importlib
import inspect
import re
import subprocess
from pathlib import Path
from typing import get_args, get_type_hints

import pytest

import scripts.check_docs_consistency as docs_checker
from scripts.check_docs_consistency import (
    _check_active_command_policy,
    _check_generated_docs_artifacts,
    _check_hardcoded_module_claims,
    _check_python_api_calls,
)
from scripts.check_docs_consistency import (
    main as check_docs,
)


def _write_qmd(path: Path, code: str) -> Path:
    path.write_text(f"```python\n{code}\n```\n", encoding="utf-8")
    return path


def test_docs_checker_accepts_normal_ignored_quarto_state():
    repo = Path(__file__).resolve().parents[1]
    quarto_state = repo / "docs" / ".quarto"
    created = not quarto_state.exists()
    quarto_state.mkdir(exist_ok=True)
    try:
        assert check_docs() == 0
    finally:
        if created:
            quarto_state.rmdir()


def _init_docs_repo(path: Path) -> Path:
    docs = path / "docs"
    docs.mkdir()
    (docs / ".gitignore").write_text("/.quarto/\n/_site/\n/*.html\n", encoding="utf-8")
    subprocess.run(["git", "init", "-q", str(path)], check=True)
    subprocess.run(["git", "-C", str(path), "add", "docs/.gitignore"], check=True)
    return docs


def test_generated_docs_check_accepts_ignored_render_state(tmp_path):
    docs = _init_docs_repo(tmp_path)
    (docs / ".quarto").mkdir()
    (docs / ".quarto" / "cache.json").write_text("{}\n", encoding="utf-8")
    (docs / "_site").mkdir()
    (docs / "_site" / "index.html").write_text("rendered\n", encoding="utf-8")
    (docs / "guide.html").write_text("rendered\n", encoding="utf-8")

    assert _check_generated_docs_artifacts(tmp_path) == []


def test_generated_docs_check_rejects_tracked_render_output(tmp_path):
    docs = _init_docs_repo(tmp_path)
    (docs / "_site").mkdir()
    generated = docs / "_site" / "index.html"
    generated.write_text("rendered\n", encoding="utf-8")
    subprocess.run(
        ["git", "-C", str(tmp_path), "add", "-f", "docs/_site/index.html"],
        check=True,
    )

    assert _check_generated_docs_artifacts(tmp_path) == [
        "Tracked generated docs artifact: docs/_site/index.html"
    ]


def test_generated_docs_check_rejects_unignored_render_output(tmp_path):
    docs = tmp_path / "docs"
    docs.mkdir()
    (docs / "_site").mkdir()

    assert _check_generated_docs_artifacts(tmp_path) == [
        "Generated docs artifact is not ignored: docs/_site"
    ]


def test_active_command_policy_rejects_windows_launch_commands(tmp_path):
    path = tmp_path / "README.md"
    path.write_text(
        "```powershell\n.\\scripts\\project_python.ps1 scripts/audit_traminex.py\n```\n",
        encoding="utf-8",
    )

    problems = _check_active_command_policy(path)

    assert len(problems) == 1
    assert "Windows/PowerShell launch command" in problems[0]


def test_active_command_policy_accepts_datasci_ubuntu_commands(tmp_path):
    path = tmp_path / "README.md"
    path.write_text(
        "```bash\n"
        "PY=~/miniforge3/envs/datasci/bin/python\n"
        "$PY scripts/audit_traminex.py\n"
        "```\n",
        encoding="utf-8",
    )

    assert _check_active_command_policy(path) == []


def test_hardcoded_module_claim_check_rejects_transient_count(tmp_path):
    path = tmp_path / "API_CONSISTENCY.md"
    path.write_text("Bare mypy checks all 78 package modules.\n", encoding="utf-8")

    problems = _check_hardcoded_module_claims(path)

    assert len(problems) == 1
    assert "hard-coded package-module count" in problems[0]


def test_hardcoded_module_claim_check_accepts_derived_wording(tmp_path):
    path = tmp_path / "API_CONSISTENCY.md"
    path.write_text(
        "Bare mypy checks the complete configured package module set.\n",
        encoding="utf-8",
    )

    assert _check_hardcoded_module_claims(path) == []


def test_active_status_check_rejects_stale_not_started_claim(tmp_path):
    path = tmp_path / "ROADMAP.md"
    path.write_text(
        "| Pre-R6 full audit and ecosystem interval | Required before R6; not started |\n",
        encoding="utf-8",
    )
    check = getattr(docs_checker, "_check_stale_active_status", lambda unused: [])

    problems = check(path)

    assert len(problems) == 1
    assert "full-audit status" in problems[0]


def test_active_performance_check_rejects_unqualified_speedup_range(tmp_path):
    path = tmp_path / "demo.py"
    path.write_text("For massive speedups (6-919x), install Numba.\n", encoding="utf-8")
    check = getattr(docs_checker, "_check_unqualified_speedup_claims", lambda unused: [])

    problems = check(path)

    assert len(problems) == 1
    assert "unqualified speedup claim" in problems[0]


def test_active_status_count_and_performance_prose_is_current():
    changelog = (Path(__file__).resolve().parents[1] / "CHANGELOG.md").read_text(encoding="utf-8")
    roadmap = (Path(__file__).resolve().parents[1] / "ROADMAP.md").read_text(encoding="utf-8")
    demo = (Path(__file__).resolve().parents[1] / "examples" / "p0_features_demo.py").read_text(
        encoding="utf-8"
    )

    assert re.search(r"(?i)\ball\s+\d+\s+(?:package\s+|source\s+)?modules\b", changelog) is None
    assert (
        "Pre-R6 full audit and ecosystem interval | Required before R6; not started" not in roadmap
    )
    assert re.search(r"(?i)(?:massive|expected)\s+speedups?", roadmap + demo) is None
    assert "benchmarks.run_benchmarks" in demo
    for qualifier in ("workload", "hardware", "versions"):
        assert qualifier in demo.lower()


def test_constructor_order_and_family_startup_file_disclosures_are_complete():
    repo = Path(__file__).resolve().parents[1]
    parity = (repo / "docs" / "API_PARITY.md").read_text(encoding="utf-8")
    ecosystem = (repo / "ECOSYSTEM.md").read_text(encoding="utf-8")

    # Owner decision of 2026-10-19: every constructor takes R's seqstatl order by default.
    for constructor in (
        "seqdef",
        "seqformat",
        "seqdef_from_long",
        "TSE_to_STS",
        "HSPELL_to_STS",
        "events_to_state_sequences",
    ):
        assert constructor in parity
    assert "`seqstatl` order" in parity
    assert "keeps the caller's order" in parity
    assert "both point here" not in ecosystem
    assert "tool-specific startup files are optional" in ecosystem.lower()


def test_python_api_call_check_rejects_invalid_documented_keywords(tmp_path):
    path = _write_qmd(
        tmp_path / "invalid.qmd",
        "M = tx.seqtransn(X, normalize=False, as_dataframe=True)",
    )

    problems = _check_python_api_calls(path)

    assert len(problems) == 1
    assert "seqtransn" in problems[0]
    assert "unexpected keyword argument 'normalize'" in problems[0]


def test_python_api_call_check_accepts_current_transition_signatures(tmp_path):
    path = _write_qmd(
        tmp_path / "valid.qmd",
        "\n".join(
            [
                "counts_per_sequence = tx.seqtransn(X, norm=False)",
                "count_matrix = tx.seqtrate(X, count=True)",
                "rate_matrix = tx.seqtrate(X, count=False)",
            ]
        ),
    )

    assert _check_python_api_calls(path) == []


def test_python_api_call_check_validates_direct_traminex_imports(tmp_path):
    path = _write_qmd(
        tmp_path / "direct-import.qmd",
        "\n".join(
            [
                "from traminex.diss import dissassoc",
                "result = dissassoc(D, groups, labels=groups)",
            ]
        ),
    )

    problems = _check_python_api_calls(path)

    assert len(problems) == 1
    assert "dissassoc" in problems[0]
    assert "unexpected keyword argument 'labels'" in problems[0]


def test_python_api_call_check_rejects_incomplete_calls(tmp_path):
    path = _write_qmd(
        tmp_path / "missing-argument.qmd",
        "\n".join(
            [
                "from traminex.distances import seqLLCS",
                "length = seqLLCS(first)",
            ]
        ),
    )

    problems = _check_python_api_calls(path)

    assert len(problems) == 1
    assert "seqLLCS" in problems[0]
    assert "missing a required argument: 'seq2'" in problems[0]


def test_python_api_call_check_rejects_removed_root_calls(tmp_path):
    path = _write_qmd(
        tmp_path / "removed-api.qmd",
        "result = tx.definitely_removed_api(X)",
    )

    problems = _check_python_api_calls(path)

    assert len(problems) == 1
    assert "definitely_removed_api" in problems[0]
    assert "not a callable traminex API" in problems[0]


def test_python_api_call_check_rejects_removed_direct_import_calls(tmp_path):
    path = _write_qmd(
        tmp_path / "removed-direct-import.qmd",
        "\n".join(
            [
                "from traminex import definitely_removed_api",
                "result = definitely_removed_api(X)",
            ]
        ),
    )

    problems = _check_python_api_calls(path)

    assert len(problems) == 1
    assert "definitely_removed_api" in problems[0]
    assert "not a callable traminex API" in problems[0]


def test_python_api_call_check_resolves_unaliased_submodule_imports_like_python(tmp_path):
    path = _write_qmd(
        tmp_path / "unaliased-submodule-import.qmd",
        "\n".join(
            [
                "import traminex.diss",
                "counts = traminex.seqtransn(X, norm=False)",
            ]
        ),
    )

    assert _check_python_api_calls(path) == []


def test_python_api_call_check_binds_aliased_submodule_imports_like_python(tmp_path):
    """``import traminex.<m> as alias`` binds ``getattr(traminex, m)`` (S7-06).

    Python compiles it to IMPORT_FROM, the same lookup as
    ``from traminex import <m>``. For the five public functions named like
    their submodule (``disstree``, ``seqdiff``, ``seqemlt``, ``seqformat``,
    ``seqtree``) that is the function, so ``st.seqtree_fit`` fails at run time;
    the checker bound the submodule and accepted it (audit S7-06, 2026-09-22).
    """
    path = _write_qmd(
        tmp_path / "aliased-submodule-import.qmd",
        "\n".join(
            [
                "import traminex.seqtree as st",
                "import traminex.diss as d",
                "fit = st.seqtree_fit(X, y)",
                "tree = st(X, predictors)",
                "total = d.dissvar(D)",
            ]
        ),
    )

    problems = _check_python_api_calls(path)

    assert len(problems) == 1, problems
    assert "traminex.seqtree_fit is not a callable traminex API" in problems[0]


def test_pinned_source_versions_match_vendored_and_ecosystem_contracts():
    repo = Path(__file__).resolve().parents[1]
    duplicates = (repo / "DUPLICATES.md").read_text(encoding="utf-8")
    ecosystem = (repo / "ECOSYSTEM.md").read_text(encoding="utf-8")

    expected = {}
    for package in ("TraMineR", "TraMineRextras"):
        description = (repo / "SOURCES" / package / "DESCRIPTION").read_text(encoding="utf-8")
        match = re.search(r"(?m)^Version:\s*(\S+)\s*$", description)
        assert match is not None
        expected[package] = match.group(1)

    weightedcluster = re.search(
        r"WeightedCluster ([0-9][0-9A-Za-z.-]*) is pyweightedcluster's\s+canonical R source",
        ecosystem,
    )
    assert weightedcluster is not None
    expected["WeightedCluster"] = weightedcluster.group(1)

    for package, version in expected.items():
        assert f"| {package} | **{version}** |" in duplicates


def test_active_parity_doc_uses_evergreen_golden_corpus_reference():
    repo = Path(__file__).resolve().parents[1]
    parity = (repo / "docs" / "API_PARITY.md").read_text(encoding="utf-8")

    assert re.search(r"\b\d+-file current-oracle\s+golden corpus", parity) is None
    assert re.search(
        r"current-oracle golden\s+corpus under `tests/r_goldens`",
        parity,
    )


def test_dissassoc_inventory_lists_every_supported_permutation_mode():
    repo = Path(__file__).resolve().parents[1]
    duplicates = (repo / "DUPLICATES.md").read_text(encoding="utf-8")
    a6 = duplicates.split("- **A6", 1)[1].split("- **A7", 1)[0]
    dissassoc = importlib.import_module("traminex.diss").dissassoc
    supported = get_args(get_type_hints(dissassoc)["weight_permutation"])

    assert "not ported" not in a6
    for mode in supported:
        assert f"`{mode}`" in a6


def test_diss_docstring_test_references_exist():
    repo = Path(__file__).resolve().parents[1]
    diss_module = importlib.import_module("traminex.diss")
    missing = []
    for name, value in vars(diss_module).items():
        if not inspect.isfunction(value) or value.__module__ != diss_module.__name__:
            continue
        for reference in re.findall(r"``(tests/[^`]+\.py)``", inspect.getdoc(value) or ""):
            if not (repo / reference).is_file():
                missing.append((name, reference))

    assert missing == []


def test_every_paired_qmd_is_an_include_wrapper_and_so_cannot_drift(tmp_path):
    """All seven ``docs/`` pairs wrap their ``.md`` instead of duplicating it (W-08).

    Four of them used to be standalone copies, and one -- ``DISTANCES`` -- had
    actually drifted on a user-visible backend switch while every gate passed
    (cycle-6 F-005). The gate added then compares only the backend-switch
    vocabulary, so any other divergence was still possible. A wrapper removes
    the possibility rather than testing for its absence.
    """
    docs = Path(__file__).resolve().parents[1] / "docs"
    paired = [qmd for qmd in sorted(docs.glob("*.qmd")) if qmd.with_suffix(".md").is_file()]

    assert [qmd.stem for qmd in paired] == [
        "API_CONSISTENCY",
        "DISTANCES",
        "DISTANCES_NORMALIZATION",
        "INDEX",
        "PLOTTING",
        "SEQTREE",
        "USAGE_new_apis",
    ]

    for qmd in paired:
        text = qmd.read_text(encoding="utf-8")
        included = {name.strip() for name in docs_checker._INCLUDE_DIRECTIVE.findall(text)}
        assert (
            qmd.with_suffix(".md").name in included
        ), f"docs/{qmd.name} is a standalone duplicate of its .md and can drift"
        # A wrapper is front matter plus the include: no prose of its own.
        body = text.split("---", 2)[-1]
        prose = [
            line
            for line in body.splitlines()
            if line.strip() and not docs_checker._INCLUDE_DIRECTIVE.search(line)
        ]
        assert prose == [], f"docs/{qmd.name} carries its own prose: {prose[:3]}"


def test_the_docs_reader_follows_include_directives():
    """Content checks must see the rendered page, not the wrapper.

    Without this the switch and parameter-naming checks read a wrapper's front
    matter and reported the prose as missing.
    """
    docs = Path(__file__).resolve().parents[1] / "docs"

    wrapper = docs_checker._read(docs / "DISTANCES.qmd")
    rendered = docs_checker._read_rendered(docs / "DISTANCES.qmd")

    assert "TRAMINEX_NUMBA" not in wrapper
    assert "TRAMINEX_NUMBA" in rendered
    assert docs_checker._read(docs / "DISTANCES.md") in rendered


def test_cited_audit_run_paths_must_be_tracked(tmp_path):
    """A document must not cite audit evidence that exists only locally (S7-13).

    ``CHANGELOG.md`` and ``docs/API_PARITY.md`` cited
    ``reports/audit/oracle-2214-substitution-audit-20260920/``, which
    ``.gitignore`` excluded with no negation, so the citation resolved only on
    the author's machine (audit S7-13, 2026-09-22).
    """
    from scripts.check_docs_consistency import _check_cited_audit_runs_are_tracked

    doc = tmp_path / "NOTES.md"
    doc.write_text(
        "Evidence: `reports/audit/run-a/`.\n"
        "Report: `reports/audit/run-b/AUDIT_REPORT.md`.\n"
        "Probe: `reports/audit/run-b/evidence/probe.R`.\n"
        "Generated: `reports/audit/traminex_api_audit.json`.\n",
        encoding="utf-8",
    )
    tracked = {"reports/audit/run-b/AUDIT_REPORT.md", "reports/audit/traminex_api_audit.json"}

    problems = _check_cited_audit_runs_are_tracked([doc], tracked)

    assert len(problems) == 2, problems
    assert "reports/audit/run-a/" in problems[0]
    assert "reports/audit/run-b/evidence/probe.R" in problems[1]

    tracked.add("reports/audit/run-a/REMEDIATION_REPORT.md")
    tracked.add("reports/audit/run-b/evidence/probe.R")
    assert _check_cited_audit_runs_are_tracked([doc], tracked) == []


def test_hardcoded_module_claim_check_rejects_package_source_wording(tmp_path):
    """ "all 78 package source modules" slipped past the pattern (S7-14).

    ``stubs/README.md`` carried that claim, and 403/32/312 stub totals, after the
    code had moved to 79 modules and 402/33/311; the pattern accepted only one of
    "package" or "source", and the file was not among the active documents
    (audit S7-14, 2026-09-22).
    """
    path = tmp_path / "README.md"
    path.write_text(
        "`make typecheck` runs bare mypy across all 78 package source modules.\n",
        encoding="utf-8",
    )

    problems = _check_hardcoded_module_claims(path)

    assert len(problems) == 1, problems


def test_stub_readme_carries_no_transient_counts():
    readme = Path(__file__).resolve().parents[1] / "stubs" / "README.md"

    assert _check_hardcoded_module_claims(readme) == []
    text = readme.read_text(encoding="utf-8")
    assert re.search(r"\b\d+ public functions\b", text) is None
    assert re.search(r"\b\d+ comparable signatures\b", text) is None


def test_python_api_call_check_reads_quarto_python_chunks(tmp_path):
    """Executable ```{python} chunks are checked like ```python fences (S7-21).

    The fence pattern matched only ```python, so the Quarto chunks of
    ``docs/TraMineR_examples.qmd`` were never checked (audit S7-21, 2026-09-22).
    """
    path = tmp_path / "chunks.qmd"
    path.write_text(
        "```{python}\nimport traminex as tx\ntx.seqstatd(X, definitely_not_a_keyword=1)\n```\n",
        encoding="utf-8",
    )

    problems = _check_python_api_calls(path)

    assert len(problems) == 1, problems
    assert "definitely_not_a_keyword" in problems[0]


def test_oracle_version_check_requires_the_vendored_versions(tmp_path):
    """Active parity docs must name the vendored oracle (DOCS_POLICY; S7-21)."""
    from scripts.check_docs_consistency import _check_oracle_versions_named

    repo = tmp_path / "repo"
    for package, version in (("TraMineR", "9.9-9"), ("TraMineRextras", "8.8.8")):
        description = repo / "SOURCES" / package / "DESCRIPTION"
        description.parent.mkdir(parents=True)
        description.write_text(f"Package: {package}\nVersion: {version}\n", encoding="utf-8")
    doc = tmp_path / "API_PARITY.md"
    doc.write_text("Oracle: TraMineR 2.2-14 and TraMineRextras 8.8.8.\n", encoding="utf-8")

    problems = _check_oracle_versions_named(doc, repo)

    assert len(problems) == 1, problems
    assert "TraMineR 9.9-9" in problems[0]

    doc.write_text("Oracle: TraMineR 9.9-9 and TraMineRextras 8.8.8.\n", encoding="utf-8")
    assert _check_oracle_versions_named(doc, repo) == []


def test_plotting_md_names_every_seqplot_type():
    """docs/PLOTTING.md lists every type seqplot routes (TX-AUD-20261019-086: pc was missing)."""
    from traminex.visualization.plots.dispatch import _DISPATCH

    text = (Path(__file__).resolve().parents[1] / "docs" / "PLOTTING.md").read_text(
        encoding="utf-8"
    )
    start = text.index("The `seqplot(x, type=...)` dispatcher")
    listed = set(re.findall(r"^- `([A-Za-z]+)`(?:/`([A-Za-z]+)`)?:", text[start:], flags=re.M))
    names = {name for pair in listed for name in pair if name}
    assert set(_DISPATCH) <= names, sorted(set(_DISPATCH) - names)


def test_unknown_seqdist_keyword_in_a_document_is_reported(tmp_path):
    """seqdist takes its options through ``**kwargs``, which binds any name (TX-AUD-20261019-054).

    The audit's crafted document passed docs-check although seqdist refuses the keyword at run time
    ("Unsupported kwargs for method OM"). The check now uses the keywords seqdist reads itself plus the
    per-method list its run-time check uses.
    """
    path = _write_qmd(
        tmp_path / "crafted.qmd",
        "\n".join(
            [
                "from traminex import seqdist",
                'd = seqdist(x, method="OM", subm=2.0, definitely_not_a_kw=3)',
                # An OM-only traminex option: seqdist refuses it for HAM. (An R argument
                # such as otto is dropped for HAM, as R ignores it: review I-1.)
                'e = tx.seqdist(x, "HAM", trate_smoothing=0.1)',
            ]
        ),
    )

    problems = _check_python_api_calls(path)

    assert len(problems) == 2, problems
    assert "definitely_not_a_kw" in problems[0] and "seqdist" in problems[0]
    assert "trate_smoothing" in problems[1] and "HAM" in problems[1]


def test_valid_seqdist_keywords_pass(tmp_path):
    path = _write_qmd(
        tmp_path / "valid.qmd",
        "\n".join(
            [
                'a = tx.seqdist(x, method="OM", subm="TRATE", indel="auto", with_missing=True)',
                'b = tx.seqdist(x, method="OMstran", subm=2.0, indel=1.0, otto=0.5, '
                'transindel="prob", previous=True, add_column=False, norm="none", weighted=False)',
                'c = tx.seqdist(x, method="om", subm=2.0, refseq=0, full_matrix=False)',
                "d = tx.seqdist(x, method=chosen, anything_goes=1)",
                'e = tx.seqdist(x, method="LCS", **options)',
                'f = tx.seqdist(x, method="DHD", deduplicate=True, condensed=True)',
            ]
        ),
    )

    assert _check_python_api_calls(path) == []


def test_unknown_sts_seqdef_keyword_is_reported(tmp_path):
    """seqdef refuses any extra keyword for the STS format, the default (format.py)."""
    path = _write_qmd(
        tmp_path / "seqdef.qmd",
        "\n".join(
            [
                "x = tx.seqdef(data, definitely_not_a_kw=3)",
                'y = tx.seqdef(data, informat="STS", alphabet=["A", "B"], also_not_a_kw=1)',
                'z = tx.seqdef(spells, informat="SPELL", process=False)',
            ]
        ),
    )

    problems = _check_python_api_calls(path)

    assert len(problems) == 2, problems
    assert "definitely_not_a_kw" in problems[0]
    assert "also_not_a_kw" in problems[1]


def test_every_public_kwargs_callable_is_classified():
    """A public ``**kwargs`` function is either checked or named as a forwarder (TX-AUD-20261019-054)."""
    import traminex as tx

    public_kwargs = {
        name
        for name in tx.__all__
        if callable(getattr(tx, name))
        and any(
            parameter.kind is inspect.Parameter.VAR_KEYWORD
            for parameter in inspect.signature(getattr(tx, name)).parameters.values()
        )
    }

    checked = {"seqdist", "seqdef", *docs_checker._KWARGS_SEQUENCE_BUILDERS}
    assert len(public_kwargs) == 29
    assert set(docs_checker._KWARGS_FORWARDERS) | checked == public_kwargs
    assert not set(docs_checker._KWARGS_FORWARDERS) & checked


def test_parenthesised_traminex_import_gives_its_names(tmp_path):
    """A multi-line ``from traminex import (`` block was read as the API name '(' (TX-AUD-20261019-055).

    The imported names are the refs; an alias is a local name, not an API name.
    """
    path = tmp_path / "imports.md"
    path.write_text(
        "```python\nfrom traminex import (\n    seqdef,\n    seqdist as sd,  # a comment\n)\n```\n"
        "Inline: `from traminex import seqstatd, seqtrate`.\n",
        encoding="utf-8",
    )

    assert docs_checker._extract_doc_api_refs(path) == {"seqdef", "seqdist", "seqstatd", "seqtrate"}


def test_parenthesised_import_of_a_non_public_name_is_reported(tmp_path):
    path = tmp_path / "imports.md"
    path.write_text(
        "```python\nfrom traminex import (seqdef, not_a_public_name)\n```\n", encoding="utf-8"
    )

    assert docs_checker._extract_doc_api_refs(path) == {"seqdef", "not_a_public_name"}


#: Calls seqdist accepts at run time that the first version of the keyword check refused
#: (batch-16 review I-1): R arguments a method does not use, which seqdist drops as R does
#: (TX-AUD-20261010-071), and check_costs, check_max_size and opt_args.
_SEQDIST_CALLS_THE_RUN_TIME_ACCEPTS = (
    'method="HAM", indel=1',
    'method="LCS", indel=1',
    'method="LCP", h=0.5',
    'method="OM", indel=1, subm="CONSTANT", check_costs=True',
    'method="OM", indel=1, subm="CONSTANT", check_max_size=True',
    'method="OM", indel=1, subm="CONSTANT", opt_args=None',
)


@pytest.mark.parametrize("arguments", _SEQDIST_CALLS_THE_RUN_TIME_ACCEPTS)
def test_seqdist_keywords_the_run_time_accepts_pass(tmp_path, arguments):
    import warnings

    import traminex as tx

    x = tx.seqdef([["a", "b", "a", "b"], ["b", "b", "a", "a"]], alphabet=["a", "b"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        eval(f"tx.seqdist(x, {arguments})", {"tx": tx, "x": x})  # the call is valid
    path = _write_qmd(tmp_path / "valid.qmd", f"d = tx.seqdist(x, {arguments})")

    assert _check_python_api_calls(path) == []


@pytest.mark.parametrize(
    "code",
    [
        'opts = {"method": "TWED"}\nd = tx.seqdist(x, **opts, nu=0.5, h=0.5, indel=1, subm=2.0)',
        'args = (x, "TWED")\nd = tx.seqdist(*args, nu=0.5)',
        'opts = {"informat": "SPELL"}\ny = tx.seqdef(spells, **opts, process=False)',
    ],
)
def test_a_method_or_format_given_through_an_expansion_is_not_assumed(tmp_path, code):
    """With ``**options`` or ``*args`` the method (or seqdef's format) is not known (review I-1)."""
    path = _write_qmd(tmp_path / "expansion.qmd", code)

    assert _check_python_api_calls(path) == []


@pytest.mark.parametrize("builder", ["HSPELL_to_STS", "TSE_to_STS", "seqdef_from_long"])
def test_an_unknown_keyword_of_a_sequence_builder_is_reported(tmp_path, builder):
    """These three refuse a keyword outside SEQ_KWARGS (validate_seq_kwargs); they were listed
    as forwarders (batch-16 review M-1)."""
    import traminex as tx

    with pytest.raises(TypeError, match="Unsupported kwargs"):
        getattr(tx, builder)(None, definitely_not_a_kw=3)  # the run time refuses it
    path = _write_qmd(
        tmp_path / "builder.qmd",
        f"x = tx.{builder}(frame, weights=w, labels=names, definitely_not_a_kw=3)",
    )

    problems = _check_python_api_calls(path)

    assert len(problems) == 1, problems
    assert "definitely_not_a_kw" in problems[0] and builder in problems[0]


@pytest.mark.parametrize(
    "statement",
    [
        "from traminex import (\n    seqdef,  # see (R)\n    seqdist,\n)",
        "from traminex import seqdef, \\\n    seqdist",
    ],
)
def test_an_import_with_a_parenthesis_in_a_comment_or_a_continuation_gives_its_names(
    tmp_path, statement
):
    """The reader stopped at the first ')' and at the end of a line (batch-16 review M-2)."""
    path = tmp_path / "imports.md"
    path.write_text(f"```python\n{statement}\n```\n", encoding="utf-8")

    assert docs_checker._extract_doc_api_refs(path) == {"seqdef", "seqdist"}
