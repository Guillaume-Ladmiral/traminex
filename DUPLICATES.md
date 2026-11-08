# Duplicate Function Inventory & Convergence Decisions

Companion to `ECOSYSTEM.md`. Master copy lives in `PROJECT_TRAMINEX`.
Last reviewed: 2026-07-18 (reconciled completed Phases 1–3 with the current
Ubuntu/R-oracle and Flexscanx parity contracts). The 2026-06-15 review clarified
that **rflexscan is flexscanx's starting point, not its source of truth**:
flexscanx corrects rflexscan/FleXScan's algorithmic/mathematical errors; only
flexscanx's `sequence/` copies are an R-anchored parity target, via
traminex/TraMineR.

This file lists every place where pyweightedcluster or flexscanx reimplements
something that belongs to traminex's domain, compares the two implementations
**against the original R/C++ packages (the source of truth)**, and records a
recommendation. Possible recommendations:

- **REPLACE** — drop the local copy, use traminex.
- **PROMOTE** — the better implementation currently lives in the *sibling*, not
  in traminex; move/port it up into traminex, then have the sibling use it.
- **KEEP** — deliberate local code; it is genuinely that project's own.
- **DECIDE** — needs the maintainer's judgement before acting.

> The scheduled A1–A8 and B1–B5 migrations are complete; the status sections
> below describe the current implementation. The later inventory preserves the
> evidence and ownership decisions that led to those changes. Broader
> cross-project migration remains unscheduled and requires maintainer approval.

## Phase 1 status: COMPLETE (2026-06-12) — A1, A2, A3, A4 done

pyweightedcluster's distance/cost stack now resolves to traminex, anchored on
TraMineR 2.2-14. One commit per item; both suites green at every step.

- **A3 (costs) — DONE.** `weightedcluster.costs` removed (no internal consumers;
  the clustering pipeline takes caller-supplied `substitution_costs`). traminex
  `seqsubm` is the cost home; new
  `tests/test_phase1_a3_costs_traminex_parity.py` pins CONSTANT and
  weighted/unweighted TRATE to R `TraMineR::seqcost`. The old local matrix was
  verified to agree with R before removal.
- **A1 (distances) — DONE.** `_seqdist.sequences_distance[_between]` delegates
  HAM/OM/LCS/DHD to `traminex.distances.pairwise`; the private 3-tier backend
  (compiled `_wc_core` / numba / pure `_fallback`) is no longer used for
  distances (dead distance kernels stripped from `_fallback.py`/`_numba.py`).
  Normalization is **pinned to R** for the `normalize=True` form clustering uses:
  HAM == R `HAM(sm=unit, norm="maxlength")`; OM == R `OM(sm=2, indel=1,
  norm="maxlength")`; LCS == R `LCS(norm="gmean")`; DHD == R
  `seqdist(method="HAM"/"DHD", sm=...)`. (R-faithful fix: `normalize=False` LCS is
  now R's raw `|a|+|b|-2*lcs`, not half of it; clustering only uses
  `normalize=True`.) Required one traminex change — DHD now honours an explicit
  `subm` (its own commit + R goldens `tests/r_goldens/DHD_subm_*.csv`), since R's
  DHD/HAM honour a user `sm` and traminex previously ignored it.
- **A2 (OMspell) — DONE.** `weightedcluster._omspell` removed;
  `metric="OMSPELL"` delegates to traminex's R-faithful OMspell. `timecost`->R
  `expcost`, `tpow=1`, `normalize=True`/`False` -> R `norm="maxlength"`/`"none"`
  (`tests/test_phase1_a2_omspell_traminex_parity.py`). The old local OMspell
  omitted R's `-1` spell-duration shift and was not R-faithful.
- **A4 (`normalize_metric`) — DONE (KEEP, folded into A1).**
  `weightedcluster._metrics.normalize_metric` stays as a deliberate tiny
  metric-name canonicalizer (charter section 4); it is the front door
  `_seqdist` uses before delegating to traminex, not a distance computation.

## Phase 2 status: COMPLETE (2026-06-13) — A5, A6, A7, A8 done

traminex is now the single faithful, **weighted** home of the discrepancy
analysis and the regression tree, anchored on TraMineR 2.2-14; every
pyweightedcluster inline site consumes it. The durable checks are the
`tests/test_r_goldens_diss*.py` and `tests/test_r_goldens_disstree*.py` families
in traminex plus `tests/test_phase2_a5a8_traminex_parity.py` in
pyweightedcluster.

- **A5 (`dissvar`) — DONE** (traminex 0.0.4). R-faithful weighted discrepancy
  `dissvar(diss, weights=, squared=)`; old simplified dict version relabelled
  `dissvar_summary`. New reusable SS primitives `weighted_pairwise_sum` /
  `weighted_inertia` (ports of `tmrWeightedInertiaDist`). Goldens
  `tests/r_goldens/dissvar_phase2.csv` + shared fixture `phase2_input_*.csv`.
- **A6 (`dissassoc`) — DONE** (traminex 0.0.5). Weighted pseudo-ANOVA +
  [Pseudo F, Pseudo Fbf, Pseudo R2, Bartlett, Levene] with R's permutation
  conventions (variants `none`, `diss`, `replicate`, `rounded-replicate`,
  `group`, `permute`, and `random-sampling`; `random_state` is a Python
  addition). t0/tables R-pinned (`dissassoc_*.csv`); deterministic components
  are R-pinned and seeded p-values are tested reproducibly within Python. Old
  effect size relabelled `dissassoc_effect`.
- **A7 (`disscenter`) — DONE, name clash resolved** (traminex 0.0.6).
  `disscenter` now means R's distance-to-group-centre (incl. `allcenter`,
  `medoids_index="first"/"all"`, weighted, squared; goldens
  `disscenter_*.csv`). The old matrix double-centering helper is renamed
  `double_center_matrix` (deprecated alias `disscenter_double_center`).
- **A8 (`disstree`/`seqtree`) — DONE** (traminex 0.0.7). R-faithful
  weighted-inertia regression tree `traminex.disstree.disstree` (binary level
  partitions, permutation gate, medoids, split schedule, `disstreeleaf`,
  `dissassoc` adjustment) and `traminex.seqtree.seqtree` (weights from
  `seqdef(weights=)`). Goldens `disstree_*.csv` / `seqtree_*.csv` (node
  tables + fitted leaves + adjustment t0, weighted and unweighted). The two
  sklearn stand-ins are relabelled exploratory in
  `disstree_explore.py`/`seqtree_explore.py`; all historical names still
  resolve.
- **Consuming refactor (2.5) — DONE** (pyweightedcluster 0.0.3, dependency
  `traminex>=0.0.7`). Inline weighted maths removed from `seqtree.py`
  (`_weighted_pairwise_sum`), `_utils.py::pairwise_sum` (proptree's
  splitter), `clustassoc.py::_weighted_within_and_total` and `quality.py`'s
  fallback moment loops. At the sibling's commit e1ca49d (read 2026-10-24),
  `_utils.py`, `_fallback.py` and `quality.py` call
  `traminex.diss.weighted_pairwise_sum`; `seqtree.py` and `clustassoc.py` call
  traminex's `dissassoc`, `dissvar`, `weighted_inertia` and `dissmfacw`
  instead. Equivalence proven before deletion
  in `tests/test_phase2_a5a8_traminex_parity.py` (verbatim copies of the
  removed maths as references). The compiled `_wc_core` quality kernel (A9,
  WeightedCluster's own) is untouched.

## Phase 3 status: COMPLETE (2026-06-13) — B1, B2, B3, B4, B5 done

flexscanx's `sequence/` copies are proven behaviour-identical to their traminex
counterparts by parity tests that use traminex at **test time only**; flexscanx
still installs and runs **without** traminex. Run the B1–B5 files under
`FLEXSCANX_PARITY_STRICT=1`; the Flexscanx repository contract pins the exact
traminex commit used rather than treating package version `0.0.8` as a complete
semantic identity. Truth flows traminex -> flexscanx throughout.

- **B1 (`seq_dist_hamming`) — DONE, no behaviour change.** `norm=True` ==
  `seqdist(method="HAM", norm="maxlength")`; `norm=False` == the raw/default
  `seqdist(method="HAM")`; normalized condensed output is the upper triangle of
  the explicitly normalized matrix. Weighted raw == explicit-subm DHD
  (`subm[t][a,b]=w_t*(a!=b)`). Weighted `norm=True` divides by `sum(weights)` —
  a documented flexscanx convention that coincides with traminex for uniform
  weights. `tests/test_traminex_parity_b1_hamming.py`.
- **B2 (`seq_trans`) — DONE, no behaviour change; NOTICE re-attributed.**
  `normalize=True` == `seqtrate(count=False)` and `normalize=False` ==
  `seqtrate(count=True)`, both oriented from-state rows to to-state columns.
  `seqtransn(norm=...)` is a separate per-sequence transition-count statistic,
  not a matrix API. `tests/test_traminex_parity_b2_trans.py` also guards an
  asymmetric orientation case and the absent-source zero row.
- **B3 (`seq_entropy`) — DONE, behaviour fixed (flexscanx 0.0.1 -> 0.0.2).**
  The only drift. `by_time=True` now matches `seqstatd`'s normalized Entropy
  row (divide by `log(K)`, K = alphabet size) via a new `normalize` parameter
  (default True; `normalize=False` keeps the old bits). `by_time=False`
  overall-diversity is unchanged. `tests/test_traminex_parity_b3_entropy.py`.
- **B4 (`seq_freq`) — DONE, no behaviour change.** Counts and proportions ==
  `seqtab`. `tests/test_traminex_parity_b4_freq.py`.
- **B5 (plots) — DONE, no behaviour change; data-level parity.** `seq_dplot`
  distribution == `seqstatd` Frequencies; `seq_fplot` widths == `seq_freq`
  (== seqtab); `seq_iplot` row order == `sequence_index_plot`. Styling differs
  by design. `tests/test_traminex_parity_b5_plots.py`.
- **B6 (`extract_sequences`) — untouched** (flexscanx's own; no parity target).

`sequence/NOTICE.txt` describes the numerical relationships per function.
Behavior is proven under `FLEXSCANX_PARITY_STRICT=1` (all assertions execute)
and standalone install verification keeps traminex out of Flexscanx runtime
requirements. Licensing provenance is tracked separately in Flexscanx's
maintained provenance note and must not be inferred from behavioral parity.
**This completes Phase 3, the final phase of MIGRATION_PLAN.md.**

## Decisions locked (2026-06-12)

- **Weighted discrepancy/tree code (A5–A8): approved.** Make traminex the
  faithful, weighted home. Audit refinement: there is no ready-made
  `dissvar`/`dissassoc`/`disscenter` function in pyweightedcluster to move — the
  faithful weighted maths is only *inline* across its clustering modules. So
  traminex **builds** the canonical primitives from the R source (cross-checked
  against that inline maths), then pyweightedcluster's inline sites are refactored
  to use them.
- **flexscanx (Section B): lightweight but identical.** flexscanx must **install
  and run without traminex** (no runtime dependency). Its `sequence/` copies are
  kept, but they must be made **byte-for-behaviour identical** to traminex's
  canonical versions and guarded by parity tests that use traminex only at
  **test time** (a development/test dependency, never an install requirement).
- **`disscenter` name clash: documented only for now.** No separate rename pass.
  It is resolved as part of executing A7 (when the weighted `disscenter` is
  promoted into traminex). **Resolved 2026-06-13** (see Phase 2 status above):
  `disscenter` = R meaning; double-centering = `double_center_matrix`.
- **Execution: a full written `MIGRATION_PLAN.md` comes first; no code changes
  until the maintainer approves it.**

---

## Pinned source-of-truth versions (updated by WeightedCluster 2.0 Phase 0)

R goldens and parity tests are generated against these exact vendored
snapshots. The WeightedCluster row was upgraded on 2026-07-17; the original
2026-06-12 1.8-1 record remains preserved in pyweightedcluster's historical
materials. If a vendored original is ever upgraded again, regenerate the
goldens and update this table in the same commit.

| R package | Version | Release date | Vendored at |
|---|---|---|---|
| TraMineR | **2.2-14** | 2026-09-15 | `PROJECT_TRAMINEX/SOURCES/TraMineR` (archive `SOURCES/TraMineR_2.2-14.tar.gz`) |
| TraMineRextras | **0.6.9** | 2026-06-20 | `PROJECT_TRAMINEX/SOURCES/TraMineRextras` (archive `SOURCES/TraMineRextras_0.6.9.tar.gz`) |
| WeightedCluster | **2.0** | 2025-12-10 | `PROJECT_PY_WEIGHTEDCLUSTER/orig_R_pkg/WeightedCluster` (archive `orig_R_pkg/WeightedCluster_2.0.tar.gz`) |
| rflexscan | not pinned (partial reference copy, no `DESCRIPTION`; wraps FleXScan 3.1.2) | — | `PROJECT_FLEXSCANX/RFLEXSCAN_ORIGINAL` |

Notes:

- Versions and release dates are taken from each vendored `DESCRIPTION` file
  (`Version:` and `Date/Publication:`, falling back to `Date:`).
- WeightedCluster 2.0 comes from
  `https://cran.r-project.org/src/contrib/WeightedCluster_2.0.tar.gz`; the exact
  archive SHA-256 is
  `6a03fbbc090ac92498a410faab1aeb88d686ff7a81dcfdcb97edacda44f393cc`.
  The historical readable tree is `orig_R_pkg/WeightedCluster-1.8-1`; its exact
  `docs/legacy/WeightedCluster_1.8-1.tar.gz` SHA-256 is
  `daf9987cee5ac9e706a39e72f54a64e49e9fe91b0bc3cfb5e5a94d5076b2df5b`.
- rflexscan is **not** a parity target — and, unlike TraMineR/WeightedCluster,
  **not the source of truth.** flexscanx is a new implementation that *starts
  from* rflexscan and **deliberately corrects its algorithmic/mathematical
  errors** and expands it; correctness is judged against the scan statistics'
  mathematical definitions, not against rflexscan's output (see `ECOSYSTEM.md`
  "Source of truth"). The only R-anchored parity target here is flexscanx's
  `sequence/` sub-package, tested against **traminex/TraMineR** (Section B). The
  vendored rflexscan copy is a reference/starting-point only, so the missing
  version pin is acceptable.
- On the supported Ubuntu 24.04 environment, measure the live oracle instead of
  copying a workstation-specific R version into this document:
  `Rscript -e 'cat(R.version.string, "\n"); cat(as.character(packageVersion("TraMineR")), "\n"); cat(as.character(packageVersion("TraMineRextras")), "\n")'`.
  Golden generators and tests live under `PROJECT_TRAMINEX/scripts/r_goldens/`
  and `PROJECT_TRAMINEX/tests/test_r_goldens_*.py`.

## Source-of-truth facts established by reading the original R/C++ code

These facts (not assumptions) drive the recommendations below. Sources are the
vendored originals:
`PROJECT_TRAMINEX/SOURCES/TraMineR`,
`PROJECT_PY_WEIGHTEDCLUSTER/orig_R_pkg/WeightedCluster`,
`PROJECT_FLEXSCANX/RFLEXSCAN_ORIGINAL`.

1. **WeightedCluster depends on TraMineR.** Its `DESCRIPTION` says
   `Depends: ... TraMineR (>= 2.0-6)`. Its C/C++ `src/` contains **only**
   clustering algorithms (PAM, k-medoids, cluster quality) — **no distance,
   cost, or substitution-matrix code at all.** In R, all sequence distances
   come from TraMineR's `seqdist`. So any distance/cost code inside the Python
   `pyweightedcluster` is **not** a port of WeightedCluster — it is a re-creation
   of TraMineR functionality and belongs in traminex.

2. **Weighting originates in TraMineR, not WeightedCluster.** In the R source,
   per-sequence weights are stored as an attribute on the sequence object
   (`attr(seqdata, "weights")`, set by `seqdef(..., weights=)`). Every weighted
   computation reads that attribute and lives in TraMineR:
   - `seqtrate()` (transition rates → TRATE costs) — reads weights.
   - `dissvar()` (discrepancy / pseudo-variance) — takes `weights`.
   - `dissassoc()` (pseudo-F / pseudo-R², permutation test) — takes `weights`.
   - `disscenter()` (distance of each case to its group centre) — takes `weights`.
   - `disstree()` (regression tree on a distance matrix, weighted inertia
     splits) — takes `weights`.
   WeightedCluster's `as.seqtree.R` does **not** implement weighting itself; it
   reads `attr(seqdata, "weights")` and calls TraMineR's weighted internals
   (`dissvar`, `TraMineRInternalWeightedInertiaDist`, `dissassoc`).
   **Conclusion:** the canonical home for weighted distances, costs, discrepancy
   and trees is traminex. WeightedCluster's genuinely own weighted code is the
   **weighted PAM/k-medoids algorithm, cluster quality, fuzzy and CLARA**
   (its `src/*.cpp` and `R/kmedoids.R`, `clusterquality.R`, `fuzzyfunc.R`).

3. **rflexscan has no sequence-analysis code.** Its source contains spatial scan
   statistics only. So flexscanx's `sequence/` sub-package is a flexscanx-original
   addition, adapted from traminex (as its own `NOTICE.txt` already states). The
   source of truth for those functions is therefore TraMineR.

4. **rflexscan/FleXScan is flexscanx's starting point, not its oracle.** Unlike
   the traminex←TraMineR and pyweightedcluster←WeightedCluster relationships
   (faithful ports, R is the source of truth), flexscanx's spatial-scan core is a
   **new implementation that deliberately corrects critical algorithmic and
   mathematical errors in rflexscan/FleXScan** and may expand beyond it. Its
   acceptance test is **agreement with the scan statistics' mathematical
   definitions and flexscanx's own first-principles tests**, not reproduction of
   rflexscan's numbers. Where the two diverge and rflexscan is wrong, flexscanx's
   behaviour is the deliberate, documented correction — it is never adjusted to
   re-create an rflexscan bug.

---

## A. pyweightedcluster vs traminex

### A-distances (completed: traminex is the faithful implementation)

| # | Former/local concern | Current traminex equivalent | Source-of-truth check | Status |
|---|---|---|---|---|
| A1 | `_seqdist.py` HAM/OM/LCS/DHD/OMSPELL | `traminex.seqdist` / `distances.pairwise` | WeightedCluster R has no distance code; distances are TraMineR's | **REPLACED/delegated** |
| A2 | `_omspell.py` | traminex OMSPELL/OMSLEN | same as A1 | **REPLACED/delegated** |
| A3 | local substitution/indel and transition-rate costs | `traminex.distances.costs`, including weighted TRATE | TraMineR owns costs | **REPLACED/delegated** |
| A4 | `_metrics.py::normalize_metric` aliases | traminex performs the computation after name canonicalization | small helper | **KEPT as canonicalizer** |

**Current HAM normalization contract.** TraMineR and traminex default HAM to
raw `norm="none"`; explicit `norm="maxlength"` returns mismatches divided by
the maximum sequence length. Weightedcluster's normalized clustering interface
deliberately requests the explicit normalized form, while its raw interface
uses the raw form. `tests/test_phase1_a1_distances_traminex_parity.py` pins both
against the R semantics.

### A-discrepancy/trees (completed: canonical weighted primitives in traminex)

The audit found faithful weighted mathematics inline across pyweightedcluster,
not in reusable named functions. Phase 2 built the R-faithful public primitives
in traminex and refactored the consumers. The current inventory is:

| # | Concept | Current faithful version | Status |
|---|---|---|---|
| A5 | `dissvar` (weighted discrepancy / pseudo-variance) | `traminex.dissvar` | **BUILT and delegated** |
| A6 | `dissassoc` (weighted pseudo-F + permutation) | `traminex.dissassoc` | **BUILT and delegated** |
| A7 | `disscenter` (distance-to-group-centre) | `traminex.disscenter`; matrix helper is `double_center_matrix` | **BUILT; name clash resolved** |
| A8 | `seqtree` / `disstree` (weighted-inertia trees) | `traminex.seqtree` / `traminex.disstree`; sklearn helpers are explicitly exploratory | **BUILT and delegated** |

### A-own (genuinely WeightedCluster's → KEEP)

| # | Local code | Why it stays | Recommendation |
|---|---|---|---|
| A9 | `kmedoids.py`, PAM, `quality.py`, `fuzzy.py`, `range.py`, CLARA (`seqclararange.py`), `boot.py` | WeightedCluster's own weighted algorithms (`src/*.cpp`, `R/kmedoids.R` etc.) | **KEEP** — never move to traminex |
| A10 | `seqnull.py` | WeightedCluster's own (`R/seqnull.R`); generates null sequences for cluster significance | **KEEP**; only its inner sequence-generation primitive could later borrow traminex `seqgen` |

## B. flexscanx vs traminex

rflexscan has no sequence code, so everything in `sequence/` is adapted from
traminex (TraMineR). The values handled are tiny (states are `+1/-1/0`).

Per the locked decision, flexscanx keeps its own copies (no runtime traminex
dependency) but each copy must be **identical in behaviour** to traminex's
version and **parity-tested against traminex at test time**.

| # | Local code (flexscanx) | traminex equivalent (parity target) | Status |
|---|---|---|---|
| B1 | `sequence/distances.py` `seq_dist_hamming` | traminex HAM, explicit raw/`maxlength` forms | **KEPT + parity-tested** |
| B2 | `sequence/analysis.py` `seq_trans` (transition matrix) | `traminex.seqtrate(count=False/True)` | **KEPT + parity-tested** |
| B3 | `sequence/analysis.py` `seq_entropy` | `traminex.seqstatd` | **KEPT + parity-tested** |
| B4 | `sequence/analysis.py` `seq_freq` | `traminex.seqtab` / `seqstatf` | **KEPT + parity-tested** |
| B5 | `sequence/plotting.py` `seq_iplot/_dplot/_fplot` | `traminex.sequence_index_plot` etc. | **KEPT + data-level parity-tested** |
| B6 | `sequence/_utils.py` `extract_sequences` | none — turns scan results into sequences | **KEEP** (flexscanx's own job; no parity target) |

**flexscanx rule (decided):** flexscanx installs and runs **without** traminex.
traminex is a **test-only** dependency used to prove these copies stay identical;
if any parity test drifts, fix the flexscanx copy to match traminex (the source
of truth flows traminex → flexscanx, never the reverse).

## C. flexscanx vs pyweightedcluster — examined, no duplication found

The family premise is "flexscanx uses functions from TraMineR **and**
WeightedCluster." The first half is real (Section B). The second half is **not
realised in the code, and currently does not need to be:**

- flexscanx's clustering code (`_consensus.py`, `_stability.py`) operates on
  **spatial** clusters — sets of area IDs compared by Jaccard overlap
  (`match_clusters_jaccard`, `consensus_sweep`). That is flexscanx's own
  scan-statistic domain, **not** WeightedCluster's k-medoids/PAM clustering of
  sequences. No PAM, medoid, silhouette-on-sequences, or cluster-quality code
  from WeightedCluster is duplicated.
- flexscanx computes trajectory **distances** (`seq_dist_hamming`) but never
  **clusters** those trajectories itself, so it has no present need for
  WeightedCluster's algorithms.

**Conclusion:** there is no flexscanx↔pyweightedcluster duplicate to reconcile.
Future agents should not invent a WeightedCluster dependency in flexscanx to
satisfy the layering diagram; the "builds on weightedcluster" relationship is
latent (available if flexscanx ever clusters trajectories) but not active today.

## D. pyweightedcluster modules reviewed and confirmed WeightedCluster-own

For exhaustivity, these were checked against the R source and are **not** traminex
territory (they port WeightedCluster's own R files, no action needed):
`aggregate.py` (aggregatecases.R), `clustassoc.py` (clustassoc.R — covariate↔
clustering association, distinct from TraMineR's `dissassoc`), `cmpcluster.py`,
`range.py` (kmedrange/clustrange), `seqclustname.py`, `proptree.py`
(propclustering.R — but note it carries inline weighted-inertia maths, see A8),
`fuzzy.py`, `boot.py`, `_labels.py` (0-/1-based label conversion utility).

---

## Executed order (historical)

1. **A1 + A2 + A3** delegated pyweightedcluster's distances and costs.
2. **A5–A8** built the faithful weighted discrepancy/tree primitives in
   traminex and refactored the consumers.
3. **B1–B5** retained Flexscanx's standalone copies with exact-revision,
   test-time parity.
4. **A4** remained a deliberate name canonicalizer.

Full sequencing, tests, and rollback are in `MIGRATION_PLAN.md`.

## Packaging reality (affects Phases 0–2)

traminex is **not published on PyPI**; it is a local project. "Add traminex as a
dependency of pyweightedcluster" therefore means a local/path or editable install
(and, for flexscanx, a local test-only install). Establishing how the three
projects reference each other locally (path dependency, editable install, or a
private index) is a Phase 0 prerequisite, not an afterthought.

**Resolved in Phase 0 (2026-06-12):** editable installs into the shared
`datasci` env, in dependency order (traminex first); consumers declare
bare-name `traminex` deps that resolve against the pre-installed local copy. Full
rationale and the dependency-confusion caveat are in `MIGRATION_PLAN.md`
Phase 0.3. The parity harnesses added in Phase 0.1 live at
`PROJECT_TRAMINEX/tests/r_goldens_utils.py`,
`PROJECT_PY_WEIGHTEDCLUSTER/tests/traminex_parity_utils.py`, and
`PROJECT_FLEXSCANX/tests/traminex_parity_utils.py`.

## Migration method used for the completed items

For each item, in order, never as a big batch:

1. **Anchor on the R/C++ source.** The acceptance test is agreement with the
   original package's output, not with either Python copy.
2. Add traminex as a declared dependency of a runtime consumer; Flexscanx is
   the explicit exception and uses it only in parity tests.
3. Write a comparison test: same inputs through the **old local** function, the
   **traminex** function, and (where feasible) the **R** reference; assert they
   match within tolerance.
4. Only when green, switch call sites to traminex and delete the local copy.
5. Keep each project's test suite green at every step. One item per commit.

This enforces "retain the best implementation" by evidence: the winner is the one
that matches the R/C++ source of truth, regardless of which project it lives in
today — which is why A1–A3 move *into* pyweightedcluster's call sites from
traminex, while A5–A8 move *out of* pyweightedcluster up into traminex.
