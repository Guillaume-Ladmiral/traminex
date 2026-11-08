# Traminex Documentation

This page is the entry point for active project documentation. The full
authoritative map of the repository (every doc, module, and artifact) is
`INDEX.md` — this page is a curated subset of it.

## Start Here

- Full project map: `INDEX.md`
- Project overview and install: `../README.md`
- Release history (single source of truth): `../CHANGELOG.md`
- Roadmap: `../ROADMAP.md`
- Audits: every finding of the audit and repair runs is in `FINDINGS_REGISTER.csv`, and the
  runs' reports (the latest: the audit-only run of 2026-10-19) are kept in the family
  supervision repository's `docs/history/run-reports/`;
  `../traminex_audit_report_20260612.md` is a historical closeout audit of the
  TraMineR 2.2.12 era

## User Guides

- Quick usage: `USAGE_new_apis.md`
- API consistency notes: `API_CONSISTENCY.md`
- Distance methods (classic: HAM/OM/LCS/LCP/RLCP/DHD/WHAM/DUR/JSD_TRANS): `DISTANCES.md`
- New distance methods (EUCLID, CHI2, TWED, OMslen, OMspell, OMstran, NMS,
  NMSMST, SVRspell; FEATURES/FUTURE costs): `NEW_DISTANCE_METHODS.md`
- Distance normalization (classic OM/LCS family **and** TWED/OMslen/OMspell/OMstran):
  `DISTANCES_NORMALIZATION.md`
- Plotting: `PLOTTING.md`
- Dedicated plotnine visualization package:
  [`py_seqplot` ownership boundary](INDEX.qmd#data-examples-references)
- Sequence trees: `SEQTREE.md`
- OMloc: `OMLOC.md`
- TraMineR manual examples in Python: `TraMineR_examples.qmd`
  (original R code corpus: `../TraMineR_examples.md`)

## Parity and References

- API parity overview: `API_PARITY.md`
- Function-level parity table: `API_PARITY_TABLE.md`
- Machine-generated parity/audit report: `../reports/audit/traminex_api_audit.md`
- Mapping CSV artifacts: `../reports/parity/`
- References: `REFERENCES.md`

## Developer Docs

- Build notes for compiled kernels: `DEV_BUILD.md`
- Documentation policy and hygiene rules: `DOCS_POLICY.md`
- Public API manifest module: `../traminex/api_manifest.py`
- Type stubs notes: `../stubs/README.md`
- Contributor guide: `../CONTRIBUTING.md`

## Quality Checks

- `python scripts/check_docs_consistency.py`
- `python scripts/check_generated_artifacts.py`
- `python scripts/check_stub_signatures.py`
- `python -m mypy` (the same bare command as `make typecheck`)

## Examples

- Examples directory: `../examples/`
- Open notebooks quickly:
  - `python scripts/open_examples.py`
  - `python scripts/open_examples.py --lab`

## Documentation Health

- Local markdown/qmd links should resolve.
- Version strings in active docs should match package metadata.
- Generated artifacts in `docs/_site/`, `docs/.quarto/`, and `docs/*.html` are build outputs and should not be used as source of truth.
- Historical reports are stored under `../archive/status_reports/`.

Last reviewed: 2026-07-12
