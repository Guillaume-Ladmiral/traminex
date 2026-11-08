# Project Index

Authoritative map of directories, modules, and documentation. Other entry
points (`README.md`, `docs/README.md`) are curated views; if they disagree
with this file, this file wins and the disagreement is a bug.

Updated: 2026-10-24 (names every tracked root file, the prompt records, and every
entry of `docs/` and `scripts/`, which `tests/test_index_census.py` checks; corrects
`_spells.py`, `omloc.py` and the `TRAMINEX_BINDING` values; the `graphify-out/` entry
says nothing regenerates or checks it; `DATA/` is checked against the vendored
data). Earlier: 2026-10-22 (adds the tracked `.codebase-memory/`, `.github/`
and `graphify-out/` directories and two documentation pages it left out);
2026-09-20 (parity oracle moved to TraMineR 2.2-14 / TraMineRextras 0.6.9);
2026-07-18 (R5 Ubuntu 24.04 and Python 3.11/3.12 reliability contract;
TraMineR 2.2-13 / TraMineRextras 0.6.9 oracle; current API and module layout).
pyproject.toml is the single declared source of the package version.
`traminex_audit_report_20260612.md` is a historical audit (TraMineR 2.2.12 era);
the findings of the later audit and repair runs are in `docs/FINDINGS_REGISTER.csv`, and
their reports (the latest: the audit-only run of 2026-10-19) are kept in the family
supervision repository's `docs/history/run-reports/`.

## Root

- README.md (overview, quick start, quality-check commands)
- AGENTS.md (instructions for coding agents; read first, with ECOSYSTEM.md)
- CLAUDE.md (points to AGENTS.md and ECOSYSTEM.md)
- ECOSYSTEM.md (the charter of the project family: what traminex,
  pyweightedcluster, py_seqplot and flexscanx each own; owner decisions)
- DUPLICATES.md (the inventory of duplicated functions across the family and the
  convergence decisions)
- MIGRATION_PLAN.md (the convergence plan of June 2026, completed; dated, kept as
  a record)
- ISSUE_SEQDIST_WITH_MISSING_RAW_ARRAY_20260921.md and
  PYWEIGHTEDCLUSTER_CONSUMER_NOTICE_20260921.md (dated reports from a
  pyweightedcluster session; the first is resolved, see its header)
- LICENSE, NOTICE (license and notice files)
- .gitattributes (line-ending policy)
- .gitignore (ignored paths; a tracked file under reports/audit/ needs a negation
  block, docs/DOCS_POLICY.md)
- .pre-commit-config.yaml (the pre-commit hooks: whitespace, end of file, black,
  ruff)
- CHANGELOG.md (release history — single source of truth)
- ROADMAP.md, CONTRIBUTING.md (contributor and roadmap notes)
- pyproject.toml (package metadata, version 0.2.3; mypy/ruff/pytest config)
- CMakeLists.txt (compiler-free default plus opt-in pybind11; see docs/DEV_BUILD.md)
- Makefile (Ubuntu-native lint, typing, workflow, docs, stubs, generated,
  backend, release, and isolated install-profile gates)
- TraMineR_examples.md (reference corpus: R example code collected from the
  TraMineR manual; the Python adaptation lives at docs/TraMineR_examples.qmd —
  the two are companions, not duplicates)
- verify_audit_improvements.py (release-gate verification script)
- traminex_audit_report_20260612.md (historical: the closeout audit of
  2026-06-12 against R 4.4.3 and TraMineR 2.2.12; merge verdicts, R-parity spot
  checks, findings by severity)
- archive/OPEN_QUESTIONS_FOR_AUDIT.md (RESOLVED by the 2026-06-12 audit and
  archived; all items answered in the audit report's open-questions section)

## Prompt and verification records

The task prompts of the project's audit and repair runs and their `.verification.json`
records (39 files, kept at the root until then) were removed on 2026-10-26 by the family
supervision plan (owner decisions D2 and D3). They remain in git history: the last commit
holding them is d2ff400c. The knowledge they governed is in `docs/OWNER_DECISIONS.md` and
`docs/FINDINGS_REGISTER.csv`, and the environment recipe of the repair runs is in
`docs/DEV_BUILD.md`.

## Package: traminex/

- __init__.py (lazy public API resolver; metadata-backed __version__)
- api_manifest.py (canonical top-level public API list; 226 exports)
- sequences.py (STSList core model; SeqMD multidomain container)
- _sequence_builders.py (shared STSList construction helpers; strict
  seq-kwargs validation)
- format.py (seqdef, FCE_to_TSE, TSE_to_STS, HSPELL_to_STS)
- seqformat.py (seqformat, seqdef_from_long, seqstart, toPersonPeriod,
  STS_to_TSE, to_SPS, to_DSS)
- stats.py (facade) with implementations split across stats_basic.py,
  stats_complexity.py, stats_utils.py, stats_p1.py, stats_p2.py
- distances/ (pairwise distance dispatch + costs)
  - __init__.py (seqdist dispatcher, seqLLCS/seqLLCP/seqRLCP, seqsubm,
    om_trate_pairwise)
  - _validation.py (allowed kwargs per method — the authoritative 19-method
    dispatcher list: OM, OMLOC, OMSLEN, OMSPELL, OMSTRAN, HAM, DHD, CHI2,
    EUCLID, LCS, LCP, RLCP, NMS, NMSMST, SVRSPELL, TWED, WHAM, DUR,
    JSD_TRANS. The first 16 are TraMineR 2.2-14 methods, with internal
    case-normalized spellings; WHAM, DUR, and JSD_TRANS are traminex-only
    extensions.)
  - _spells.py (spell helpers the spell-based methods share: spell_encode,
    spell_length_per_position, spell_length_power)
  - _omstran.py (OMstran transition-augmented OM)
  - python_kernels.py (Python kernels incl. NMS/NMSMST/SVRspell)
  - _batched.py (the row engine's dynamic programs for many pairs at once:
    OM with per-state indels, LCS, OMslen, OMspell, OMloc, TWED, and the hard-
    and soft-matching NMS counts; bit-identical to the per-pair kernels)
  - costs.py (TRATE/CONSTANT/INDELS(LOG)/FUTURE/FEATURES substitution costs,
    gower_distance)
  - omloc.py (the separate traminex function `omloc()`; seqdist's OMloc
    method is in the kernels above; see docs/OMLOC.md)
  - _normalization.py (TraMineR DistanceCalculator::normalizeDistance port;
    the norm= vocabulary and per-method maximum-distance bounds)
- diss.py, diss_p1.py (dissimilarity utilities: dissassoc, disscenter,
  dissdomassoc, dissvar, gower_matrix, ...)
- diss_indic.py (dissindic), diss_merge.py (dissmergegroups),
  diss_mfa.py (dissmfacw)
- domassoc.py (seqdomassoc) and group_compare.py (dissCompare, seqCompare,
  seqLRT, seqBIC)
- multidomain.py (seqMD weight grid-search/optimization) and
  multidomain_cat.py (the R-faithful seqMD/seqdistmc CAT engine)
- repseq_rf.py (dissrf, seqrf), samm.py (seqsamm and its companions),
  seqemlt.py (seqemlt), discrete_rules.py (createdatadiscrete,
  seqerulesdisc)
- config.py, config_dataclasses.py (TRATEConfig, DHDConfig, DURConfig, OMConfig, LCSLCPConfig)
- _rformat.py (private: R's fixed-or-scientific text rendering of numbers
  and R's round(),
  used where labels must read as TraMineR's do)
- _rsum.py (private: R's long-double `sum()`, shared by the R-faithful tree
  splits and event-mining tests)
- _rlinpack.py (private: R's LINPACK dqrdc2 rank rule, used for the ranks of
  dissmfacw's designs and the aliased columns of seqerulesdisc's binomial fits)
- _categorical_order.py and _state_order.py (private deterministic ordering:
  ICU root collation for R-compatible label order, and the ranked-state
  contract)
- _deprecated.py (the deprecated-alias resolver shared by public wrappers)
- events/ (event sequences; `__init__.py` keeps the public import surface)
  - _object.py (the R-faithful eseq/seqelist objects: Eseq, SeqElist,
    seqecreate and their accessors)
  - _engine.py (internal port of the algorithms of TraMineR's C++
    event-mining engine, with the corrections the parity ledger records)
  - _mining.py (the R-faithful mining surface over that engine: seqefsub,
    seqeapplysub, seqecmpgroup, rules, constraints)
  - _duration.py (the numeric series behind seqedplot)
  - _helpers.py (traminex-only event helpers with no TraMineR equivalent,
    such as events_transition_matrix and frequent_ngrams)
- seqdiff.py (R-faithful position-wise discrepancy analysis)
- disstree.py and seqtree.py (R-faithful weighted trees), with
  disstree_helpers.py (disstreeleaf, disstree_assign, disstree_get_rules) and
  disstree_viz.py (DOT emission; the system GraphViz `dot` via GVPATH);
  disstree_explore.py and seqtree_explore.py (exploratory sklearn variants)
- patterns.py (frequent_subsequences); `frequent_ngrams` lives in events/
- align_p2.py (alignment helpers)
- cluster.py (cluster_* helpers, seqrep, medoid_index)
- plotting.py (facade) and visualization/ (Matplotlib renderers over series
  the analysis modules compute)
  - series.py (reshape-only adapter from seqstatd/seqmeant/... results to
    plot series; no numeric recomputation)
  - legend.py (state legends), palette.py (state colours), ticks.py (time
    axes), utils.py (layout, figure buffers, save_plot)
  - plots/ (one module per plot family): dispatch.py (the seqplot type
    router), index.py (index plots), distribution.py (state distribution),
    chronogram.py (chronograms), frequency.py (sequence frequency),
    summary.py (mean time, modal state, entropy plots), duration.py
    (duration boxplots), transitions.py (transition heatmap),
    representatives.py (representative sequences), survival.py (spell
    survival), alignment.py (pairwise alignment), discrepancy.py
    (position-wise discrepancy), event_duration.py (seqedplot),
    events_viz.py (event-mining plots), implicative.py (seqimplic),
    multidomain.py (seqplotMD), samm.py (SAMM)

## Compiled core (optional): cpp/

- cpp/src/accel_kernels.h and cpp/bindings/pybind11_module.cpp
  (OpenMP HAM/OM/LCS acceleration plus diagnostics; built via CMakeLists.txt;
  runtime selection with `TRAMINEX_BINDING=auto|python|pybind11`, auto the
  default). The older generic
  kernels remain under cpp/src.

## Docs: docs/

Rendered as a Quarto site (docs/_quarto.yml; output in docs/_site/, not
tracked). For paired files, the `.md` is the single source and the `.qmd` is a
Quarto wrapper that adds front matter and includes it, so the rendered site and
the GitHub-readable file cannot drift. Quarto does not execute plain Markdown
`python` fences. `TraMineR_examples.qmd` is the one standalone Quarto document;
its examples are plain `python` fences as well, so rendering shows them without
running them. `tests/test_doc_examples_execute.py` runs the examples of the fully
runnable documents, including this one. An entry written `X.md / .qmd` below has
such a wrapper.

- README.md (documentation hub / start here)
- _quarto.yml (the Quarto site configuration)
- .gitignore (ignores the rendered site)
- superpowers/ (only specs/2026-07-12-traminex-full-parity-closeout-design.md remains: 11 rows
  of docs/API_PARITY_DECISIONS.csv cite it as evidence; the other plans, prompts and notes of
  earlier work were removed on 2026-10-26 and remain in git history)
- INDEX.md / INDEX.qmd (this file)
- DOCS_POLICY.md (documentation governance)
- USAGE_new_apis.md / .qmd (usage examples)
- API_CONSISTENCY.md / .qmd (missing-state keywords, round-trips, plotting ax)
- CONSUMED_SURFACE.md (traminex objects and module paths the sibling projects import)
- DISTANCES.md / .qmd (core distance guide + hub: HAM, OM, LCS, LCP, RLCP,
  DHD [TraMineR transition-rate default + legacy marginal opt-in], WHAM, DUR,
  JSD_TRANS; cost utilities; backends; links to the parity-ported methods and
  normalization)
- NEW_DISTANCE_METHODS.md (the parity-ported methods added in 0.0.2:
  EUCLID, CHI2 [R-faithful step=1], TWED, OMslen, OMspell, OMstran, NMS,
  NMSMST, SVRspell; plus FEATURES and FUTURE costs. TWED/OMslen/OMspell/OMstran
  support TraMineR-parity norm= normalization (0.0.3); NMS/NMSMST/SVRspell
  refuse norm, matching TraMineR. Default and risk-weighted option routes are pinned
  by live-R goldens.)
- DISTANCES_NORMALIZATION.md / .qmd (norm= options for BOTH the classic
  OM/LCS family and the OM-family spell/warp methods OMloc/TWED/OMslen/
  OMspell/OMstran [maxlength/gmean/maxdist/YujianBo/auto], with per-method
  maximum-distance formulas and R source lines; pinned by
  tests/test_om_normalization.py)
- OMLOC.md (the TraMineR OMloc method, the separate Python `omloc()`
  convenience function, and position-dependent regular OM)
- PLOTTING.md / .qmd (plotting guide)
- SEQTREE.md / .qmd (sequence trees)
- TraMineR_examples.qmd (TraMineR manual examples adapted to Python;
  companion to root TraMineR_examples.md which holds the original R code)
- API_PARITY.md (narrative R-vs-Python coverage overview)
- API_PARITY_TABLE.md / .csv (generated; do not edit by hand)
- API_PARITY_DECISIONS.csv (the authoritative decision table, one row per R name)
- API_PARITY_LEDGER.md (the semantic parity ledger, one row per TraMineR and
  TraMineRextras export, registered S3 method and promoted internal; generated
  from API_PARITY_DECISIONS.csv by scripts/generate_semantic_parity_ledger.py;
  do not edit by hand)
- VISUALIZATION_SERIES_CONTRACT.md (the frozen boundary between traminex and
  downstream visualization libraries such as the sibling py_seqplot)
- DEV_BUILD.md (portable and pybind11/OpenMP build guide)
- REFERENCES.md (reference fixtures policy)
- OWNER_DECISIONS.md (the owner's decisions about traminex, copied from their sources)
- FINDINGS_REGISTER.csv (one row per audit or repair finding identifier, TX-AUD-… and TX-REP-…: title,
  severity, latest status, decision and fixing commits, built from the run records)
- SUPERVISOR_NOTES.md (dated findings about traminex from the family's supervision repository)

## Generated parity/audit artifacts

Source of truth for parity counts. Regenerate with `make audit PYTHON=$PY` (the
four generators below) and verify freshness with `make generated-check PYTHON=$PY`:

```bash
# $PY: the interpreter with the editable install (docs/DEV_BUILD.md, "Environment")
$PY scripts/audit_traminex.py
$PY scripts/generate_parity_table.py
$PY scripts/generate_api_mappings.py
$PY scripts/generate_semantic_parity_ledger.py
$PY scripts/check_generated_artifacts.py
```

- reports/audit/traminex_api_audit.{md,json,csv}, reports/audit/traminex_semantic_parity.json
  and docs/API_PARITY_LEDGER.md
- reports/parity/mapping_traminex_vs_TraMineR.csv,
  mapping_traminex_vs_TraMineRextras.csv
- docs/API_PARITY_TABLE.{md,csv}

## Tests: tests/

- pytest suite (run with the interpreter of docs/DEV_BUILD.md's "Environment";
  compiled-kernel cases skip when
  `_pb_kernels` is not built, which fails the pybind11 gate,
  `scripts/run_pybind11_tests.sh`). One optional smoke test, gated by
  `TRAMINEX_PARITY=1` (`tests/test_parity.py`), reads TraMineR's `actcal` data
  from the local `DATA/` folder and checks only that `seqdef` builds sequences
  from it; it compares nothing with R's output and runs in no gate
- test_phase1_regressions.py, test_phase3_api_cleanup.py (audit-remediation
  regression tests)
- test_api_consistency.py, test_public_api_dupes.py, plus unit tests for
  distances, plotting, patterns, trees, events, multidomain, workflows

## Scripts: scripts/

- audit_traminex.py (Linux audit generator; writes reports/audit/*)
- generate_parity_table.py, generate_api_mappings.py,
  annotate_traminer_list.py, calibrate_parity.py (parity tooling)
- check_docs_consistency.py, check_stub_signatures.py,
  check_generated_artifacts.py, sync_stub_public_api.py, validate_workflows.py
  (consistency gates; see Makefile release-check)
- open_examples.py and disstree_demo.py (examples/developer utilities)
- run_pybind11_tests.sh (temporary compiled-backend build and complete test gate)
- _native_link.sh (the runner's link checks)
- _native_skip_gate.py (the pytest plugin that fails the native gate on a
  compiled-kernel skip)
- check_r_golden_provenance.py (the R-golden inventory and writer-provenance
  proof, static or `--execute`)
- r_goldens/ (the R generators of the tests/r_goldens/ fixtures; see its
  README.md)
- generate_semantic_parity_ledger.py (writes docs/API_PARITY_LEDGER.md and
  reports/audit/traminex_semantic_parity.json from docs/API_PARITY_DECISIONS.csv)
- sync_sdist_manifest.py (writes pyproject.toml's sdist list from `git ls-files`;
  `make sdist-manifest`)
- check_dev_environment.py (names missing development modules before the release
  gates)
- _checkout_imports.py (makes developer commands import this checkout ahead of an
  editable install)
- _r_namespace.py (parses R NAMESPACE files for the audit tooling)
- verify_install_profiles.py (artifact inspection and six fresh-environment
  core/performance/plot/tree/compiled/fallback install profiles)

## Data, examples, references

- DATA/ (sample datasets: copies of the vendored `SOURCES/<package>/data/` files,
  byte for byte, which `tests/test_data_provenance.py` checks), DEMOS/ (demo
  material), examples/ (scripts and
  notebooks), benchmarks/ (the maintained reproducible benchmark CLI and its
  agreement-before-timing contract)
- DOCUMENTATION/ (reference PDFs for TraMineR/TraMineRextras — keep)
- SOURCES/TraMineR (2.2-14), SOURCES/TraMineRextras (0.6.9) — exact CRAN source
  oracles for parity work; see SOURCES/README.md for hashes and refresh steps
- .external/Sequenzo — retained read-only historical/audit evidence used for
  post-hoc provenance comparison only; it is not a current implementation
  input, runtime dependency, build dependency, or distribution input
- PROJECT_PY_SEQPLOT is the sibling dedicated to plotnine/ggseqplot-style
  visualization; traminex owns the numerical series it consumes.
- stubs/ (66-module external review stubs plus the ordered 226-name root;
  bidirectionally checked by scripts/check_stub_signatures.py). Installed
  distributions use inline annotations and `traminex/py.typed`; bare mypy checks
  the complete configured package module set.

## Tooling

- .github/ (GitHub Actions workflows: ci.yml, render-docs.yml, wheels.yml)

## History

- archive/README.md (provenance notes)
- archive/status_reports/ (dated audit reports and milestone summaries;
  snapshots, not source of truth)
- archive/historical_notes/
