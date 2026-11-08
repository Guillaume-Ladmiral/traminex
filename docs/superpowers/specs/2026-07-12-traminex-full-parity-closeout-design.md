# traminex 1.0 — Full Analytical Parity Closeout Design

Status: approved architecture, 2026-07-12. This design supersedes the release
boundary and environment assumptions in
`2026-06-16-traminex-v1-comprehensive-parity-design.md` while preserving its
R-golden-first discipline, ownership wall, and deliberate Python re-expression
of R's S3 display machinery.

## 1. Goal

Release traminex 1.0 only after it provides full sensible analytical parity with
TraMineR 2.2-13 and TraMineRextras 0.6.9 on Ubuntu 24.04.

"Full sensible analytical parity" means:

- every in-scope analytical R export has an R-derived Python public name
  (`.` becomes `_` where required by Python syntax);
- all documented analytical branches and options are implemented, not only the
  most common/default path;
- every deterministic analytical result is pinned to the current live-R oracle;
- stochastic routines are pinned where reproducible and otherwise receive
  deterministic component goldens plus explicit distribution/structure tests;
- known upstream bugs may remain deliberately corrected only when the divergence
  is named, justified, and regression-tested;
- R S3 printing, summary, subsetting, arithmetic dispatch, and internal helpers
  are not duplicated as one Python symbol per method.

The current 208-symbol public API is a baseline, not a release freeze. The 1.0
surface freezes only after the semantic parity ledger introduced in Phase R0 is
complete and approved.

## 2. Sources of truth

1. Behavioral oracle: vendored CRAN TraMineR 2.2-13 and TraMineRextras 0.6.9
   under `SOURCES/`.
2. Python public surface: `traminex/api_manifest.py`.
3. Semantic parity decisions: the Phase-R0 decision ledger, not the generated
   fuzzy-name audit.
4. Family ownership: `ECOSYSTEM.md` and `DUPLICATES.md`.
5. Active platform: Ubuntu 24.04, using
   `~/miniforge3/envs/datasci` for development and verification.

The generated API audit remains useful for namespace accounting, but its
`suggested` matches are string-similarity candidates and cannot establish
behavioral equivalence.

## 3. Approved scope decisions

### 3.1 Public API rule

Every in-scope analytical R export receives an R-derived public name. Existing
Pythonic entry points may remain as documented aliases when they are useful and
unambiguous. Implemented-but-internal analytical functions, such as the current
`seqeapplysub`, must either become public under the R name or receive an explicit
out-of-scope decision justified by the ledger.

The public API count is derived from approved decisions; no task may preserve
the number 208 merely because it was the historical baseline.

### 3.2 Included analytical work

- missing analytical entry points such as `seqcost`, `seqsep`, and the complete
  `seqerulesdisc` contract;
- completion of composite APIs such as the binary/ranked and duration/subsequence
  branches of `seqindic`;
- complete documented branches of already-exposed functions, including
  representatives, discrepancy permutations, multidomain association, state
  order/equivalence, survival, and linked-polyad analysis;
- the TraMineRextras SAMM family (`seqsamm`, `seqsammseq`, `seqsammeha`) and its
  required analytical data structures;
- TraMineRextras 0.6.9 changes to `dissCompare`, `seqCompare`, and `seqpolyads`;
- complete sequence construction/conversion behavior needed by supported
  analytical workflows;
- analytical support exports required by those workflows, including
  `createdatadiscrete`, `group.p`, `seqgen.missing`, `sorti`, and `sortv`;
- applicable TraMineR 2.2-13 plotting-data controls, including sequence sampling,
  while preserving the visualization ownership boundary.

### 3.3 Re-expressed rather than duplicated

- R `print.*`, `summary.*`, `str.*`, and `as.character.*` behavior is expressed
  through Python result objects, `__repr__`, explicit summaries, and conversion
  functions.
- R `plot.*` methods map to explicit traminex matplotlib/reference functions or
  to the `py_ggplotseq` sibling. A separate top-level alias is added only when it
  improves R workflow compatibility without obscuring ownership.
- setters such as `alphabet<-`, `cpal<-`, and `seqeweight<-` map to Python
  setters/methods plus the callable accessor already defined by the R namespace.

### 3.4 Excluded

- `TraMineRInternal*` and `TraMineR.check.depr.args`;
- legacy `read.tda.mdist`;
- `convert.g`, an ImageMagick document/image conversion utility;
- `pamward`, owned by the weightedcluster sibling;
- one-for-one reproduction of R S3 subset/arithmetic/display machinery;
- a plotnine implementation inside traminex (owned by `py_ggplotseq`);
- any broader flexscanx migration before the later architecture discussion.

`recodef` is re-expressed by `seqrecode`; `seqfcheck` remains internal validation.
Both decisions must be explicit in the Phase-R0 ledger rather than appearing as
unresolved fuzzy matches.

## 4. Findings that require the new plan

The historical plan stack is valuable evidence but is not execution-ready:

- it targets TraMineR 2.2-12 / TraMineRextras 0.6.8 and Windows/PowerShell;
- its Phase 7 plan freezes the public API at 208 before current semantic
  reconciliation;
- it classifies Python `seqerules` as TraMineRextras `seqerulesdisc`, but the
  current code and R sources show that Python `seqerules` faithfully ports
  TraMineR `seqerules`; TraMineRextras `seqerulesdisc` is a separate
  discrete-covariate modeling function;
- it records `seqcost` as an alias of `seqsubm`, while no current Python
  `seqcost` function or public symbol exists;
- it marks several phases complete while documented branches still raise
  `NotImplementedError`;
- it schedules `seqsep`, `sorti`, `sortv`, and `seqsplot` but the closeout status
  does not explain whether they were implemented, re-expressed, or dropped;
- it predates TraMineRextras 0.6.9 multi-sample LRT aggregation, `inc.df`, and
  `seqpolyads(core="auto")`, and TraMineR 2.2-13 plot sampling controls.

The improved plan therefore begins with evidence reconciliation, not new feature
implementation.

## 5. Phase architecture

### R0 — Current semantic ledger and final scope freeze

Account for all 243 current R namespace names exactly once. Replace fuzzy-name
status with approved semantic status, ownership, Python mapping, option coverage,
and evidence links. Derive the target public API and exact subsequent phase task
counts from this ledger.

### R1 — Core contracts, construction, and orchestration

Complete the central APIs on which other parity work depends: `seqcost`, full
`seqindic`, R-name exposure of existing analytical internals, sequence construction
and conversion, and ragged/missing/void distance contracts.

### R2 — Complete existing analytical branches

Remove unsupported documented branches from already-exposed core functions:
representatives, relative-frequency grouping, discrepancy permutations, domain
association, multidomain costs, state ordering/equivalence, sequence likelihood,
survival, linked polyads, and tree output.

### R3 — TraMineRextras completion

Implement full 0.6.9 group comparison, discriminating event rules and their data
preparation, SAMM, and remaining substantive TraMineRextras analytical helpers.

### R4 — Visualization contracts and sibling boundary

Complete current TraMineR plotting-data controls and R-name routing in the
reference matplotlib layer. Define and test the stable analytical-series contract
consumed by `py_ggplotseq`; do not duplicate plotnine rendering.

### R5 — Ubuntu reliability and release hardening

Replace Windows-era gates; declare and test the supported Python/dependency matrix;
type-check all package modules; complete stubs; measure and improve coverage;
consolidate benchmarks; and verify clean installation profiles.

### R6 — Release 1.0

Re-run the namespace ledger, all R generators, all backend and install profiles,
documentation and packaging gates; then freeze the API, bump/version/package/tag,
and publish only with explicit approval.

## 6. Cross-cutting implementation rules

1. Read every vendored R/Rd/C/C++ source that defines or calls the target behavior
   before designing a Python change.
2. Add or update the R generator and failing Python golden test before implementation.
3. Keep one coherent behavior per commit; mechanical manifest/stub/docs updates may
   accompany the behavior they describe.
4. Use pybind11 only for the optional compiled backend on this machine. The default
   installation must remain compiler-free and pure-Python capable.
5. Preserve optional dependency boundaries: core, performance, plotting, tree.
6. No cross-project migration without a separately approved sibling plan.
7. Use the `datasci` environment and native Ubuntu commands; do not add or preserve
   active Windows instructions.
8. Every phase receives a separate bite-sized execution plan after its predecessor's
   outputs fix its exact input scope.

## 7. Verification model

### Deterministic analytics

- R generator under `scripts/r_goldens/`;
- static artifact under `tests/r_goldens/`;
- focused Python test under `tests/test_r_goldens_<area>.py`;
- `rtol=0` and the repository R-golden absolute tolerance unless the R source
  demonstrates a justified numerical tolerance.

### Stochastic analytics

- reproduce the R stream when practical;
- otherwise golden the deterministic sufficient statistics/components;
- test seeded Python reproducibility, output schema, invariants, and edge cases;
- state explicitly which stochastic quantity is not cross-language reproducible.

### Visualization

- golden the producing numerical series;
- verify artist/layer data, ordering, labels, axes, grouping, and sampling;
- avoid pixel baselines except for a small, explicitly justified smoke set.

### Backends and packaging

- pure Python;
- Numba performance extra;
- pybind11/OpenMP compiled extra with byte-identity where promised;
- plotting and graphviz extras;
- core-only clean wheel installation with no compiler or optional dependencies.

## 8. Definition of done

traminex 1.0 is ready only when:

- every current R namespace name appears exactly once in the semantic ledger;
- every in-scope analytical export has an R-derived Python public name;
- every in-scope deterministic analytical branch has current-oracle evidence;
- no documented in-scope analytical branch terminates in `NotImplementedError`;
- every re-expression, exclusion, sibling-owned item, and deliberate divergence is
  explicit; each executable behavior has a test and every decision has documentation;
- the Ubuntu build, type, stub, coverage, docs, benchmark, backend, and clean-install
  gates pass;
- `traminex/api_manifest.py`, stubs, documentation, and built metadata agree on the
  final version and public surface;
- the release artifacts install and run in disposable environments;
- publication receives explicit maintainer approval.
