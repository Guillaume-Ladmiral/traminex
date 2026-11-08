# Documentation Policy

This document defines which files are authoritative and how to keep docs consistent.

## Source of Truth

- Package metadata/version:
  - `pyproject.toml`
  - `traminex/__init__.py`
- Public API coverage:
  - `traminex/api_manifest.py` (canonical export list)
  - `docs/API_PARITY_DECISIONS.csv` (authoritative semantic R namespace decisions)
  - `docs/API_PARITY_LEDGER.md` (generated human-readable decision ledger)
  - `reports/audit/traminex_semantic_parity.json` (generated machine-readable ledger)
  - `reports/audit/traminex_api_audit.json`
  - `reports/audit/traminex_api_audit.md`
- R oracle versions and artifacts:
  - `SOURCES/README.md`
  - `SOURCES/TraMineR/DESCRIPTION`
  - `SOURCES/TraMineRextras/DESCRIPTION`
- User-facing guides:
  - files under `docs/*.md` and `docs/*.qmd`
- Operational documentation:
  - `scripts/r_goldens/README.md`
  - `reports/parity/README.md`

### The `r_kind` vocabulary in `docs/API_PARITY_DECISIONS.csv`

Every row records what the R object *is*, which determines where its existence can
be verified:

| `r_kind` | Meaning | Verified against |
|---|---|---|
| `export` | R lists the name in `export(...)` | the package `NAMESPACE` |
| `s3` | R registers the method with `S3method(...)` | the package `NAMESPACE` |
| `internal` | R **defines** the object but does **not** export it, and traminex promotes it to its public API | the vendored R source, by `scripts/_r_namespace.promoted_internals` |

`internal` exists because `NAMESPACE` lists only exports and S3 methods, so a
promoted internal had no representable row and therefore no authoritative
disposition at all — seven public analytical surfaces were in that position
(cycle-4b finding F-004). The registry is
`scripts/_r_namespace.PROMOTED_INTERNALS`, and both halves of each entry are
checked: the object is defined in the recorded source file, and the package does
not export it. If R starts exporting one, the check fails and the row must move
to `export`.

Promoted S3 methods (`plot.SAMM`, `plot.seqimplic`, `plot.subseqelist`,
`plot.subseqelistchisq`, `plot.eseq`, `plot.seqelist`) are **not** `internal`:
`NAMESPACE` registers them, so they are carried as `s3`.

`scripts/check_docs_consistency.py` additionally verifies that every public
export whose sanitized name matches a defined R object has a disposition
somewhere, so this class cannot silently reappear.

## Historical vs Active

- Active docs:
  - `README.md`, `CHANGELOG.md` (single source of release history), `CONTRIBUTING.md`,
    `ROADMAP.md`, curated files in `docs/`, `scripts/r_goldens/README.md`, and
    `reports/parity/README.md`
- Historical snapshots:
  - `archive/status_reports/`
- Historical files are retained for traceability and may contain stale metrics/version text.

## Audit and remediation deliverables: what is kept

Audit and repair runs used to write a run directory under `reports/audit/`, with the run's
report and evidence index tracked and its evidence tree ignored. On 2026-10-26 the family
supervision plan (parent repository `PROJECT_SEQANALYSIS_and_SPATSCAN`;
owner decisions D2 and D3) removed those run directories, the loose `full_codebase_*`
reports and the root task prompts. The tracked ones remain in git history; the last commit
holding them is d2ff400c. The knowledge is kept here:

| Knowledge | Where |
|---|---|
| every finding identifier (TX-AUD-…, TX-REP-…): title, severity, latest status, decision, fixing commits | `docs/FINDINGS_REGISTER.csv` |
| the owner's decisions | `docs/OWNER_DECISIONS.md` and `ECOSYSTEM.md` |
| the runs' report texts and finding files | the supervision repository's `docs/history/run-reports/` |
| `traminex_api_audit.{csv,json,md}`, `traminex_semantic_parity.json` | `reports/audit/`: generated artifacts listed under "Source of Truth" above, refreshed by `make audit` and gated by `make generated-check` |

A new audit or repair run keeps its working evidence outside the repository and records
its findings in `docs/FINDINGS_REGISTER.csv` and its decisions in `docs/OWNER_DECISIONS.md`
(the family plan's rule for new evidence). `.gitignore` still excludes `reports/audit/`
except the four generated artifacts; a file to be tracked there needs its own negation,
because a negation cannot reach into an excluded directory and `git add` refuses an
ignored path. `scripts/check_docs_consistency.py` rejects an active document that cites an
untracked path under `reports/audit/<run>/`.

## Generated Artifacts

Do not treat generated docs as source files:

- `docs/_site/`
- `docs/.quarto/`
- `docs/*.html`

The semantic ledger has one hand-reviewed source and two generated views:

- edit `docs/API_PARITY_DECISIONS.csv` when a semantic decision changes;
- regenerate `docs/API_PARITY_LEDGER.md` and
  `reports/audit/traminex_semantic_parity.json` with
  `scripts/generate_semantic_parity_ledger.py`;
- do not infer semantic parity from the fuzzy-name suggestions in the historical
  `traminex_api_audit` reports.

The decision CSV uses the closed option-coverage vocabulary `FULL`, `PARTIAL`,
`NONE`, and `N/A`. `FULL` covers every effective source-active option/default
for the meaningful mapping; `PARTIAL` names unsupported branches; `N/A` is for
protocol, ownership, or excluded surfaces and cannot hide an analytical gap.

The semantic validator requires package-correct public-R source/Rd evidence and
a focused behavioral test for executable mappings, at least one of which names
the mapped Python function. An `EXACT` analytical row must also cite the R
generator (`scripts/r_goldens/*.R`) that records R's result. `DIVERGENCE`
notes name the exact observable boundary; `RE_EXPRESSED` notes name the Python
replacement; stochastic notes separate deterministic sufficient components
from the actual R and NumPy streams. Known placeholder prose and known
unsupported-option `FULL` claims are rejected.

Regenerate when needed from `.qmd` inputs.

## Consistency Rules

- Active docs must not claim a package version different from metadata.
- Active parity docs must use TraMineR 2.2-14 and TraMineRextras 0.6.9;
  `scripts/check_docs_consistency.py` checks that `docs/API_PARITY.md`,
  `docs/INDEX.md` and `scripts/r_goldens/README.md` name the versions vendored
  under `SOURCES/`.
- Active development commands must target Ubuntu 24.04 and the `datasci`
  environment; executable Windows/PowerShell commands belong only in historical
  snapshots and are rejected in active docs.
- Avoid hardcoding fast-changing metrics, including test counts, generator
  counts, fixture counts, and “all N package modules” claims.
- Date active audit/status claims and make completion conditional on the
  corresponding final verification evidence; do not retain obsolete “not
  started” wording after work begins.
- Qualify performance prose by workload, hardware, and versions, and direct
  readers to the maintained agreement-before-timing benchmark instead of
  promising universal numeric speedups.
- Prefer linking to generated reports (`reports/audit/...`) or showing a
  derivation command for parity counts.
- Remove or archive references to deleted aliases and deprecated APIs.
- Keep paired `.qmd` and `.md` sources semantically synchronized; do not track
  Quarto build output.

## Maintenance Workflow

1. Edit authoritative code, decision CSV, or `.qmd`/Markdown sources.
2. Regenerate owned audit, semantic, parity, or R-golden artifacts through
   their scripts; never synthesize R expected values from Python.
3. Inspect every generated diff and run its focused tests.
4. Run `scripts/check_generated_artifacts.py` and
   `scripts/check_docs_consistency.py` with an interpreter that has the editable
   install and the `dev` extra (docs/DEV_BUILD.md, "Environment").
5. Update release notes/changelog only when user-facing behavior changed.
