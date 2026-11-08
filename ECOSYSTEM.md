# Project Family Charter — traminex / pyweightedcluster / py_seqplot / flexscanx

**The master copy of this charter lives in `PROJECT_TRAMINEX`, and since
2026-10-26 it is the only copy: each sibling's `ECOSYSTEM.md` points here and
keeps only its project-specific statements (see "Keeping this file in sync").**

**This file is for any AI assistant or model, not just one tool.** Each project
carries an `AGENTS.md` (read by most coding agents) that points to its applicable
project guidance. Tool-specific startup files are optional; when a project has
one (for example `CLAUDE.md`, `GEMINI.md`, or `.cursorrules`), it should defer to
the same tool-neutral guidance rather than duplicate it. This charter is the
single, tool-neutral source.

Current correction: 2026-10-08 — verified the plotting repository's retained 2026-08-14 rename to PROJECT_PY_SEQPLOT / py_seqplot. Traminex remediation comes first; weightedcluster and flexscanx checks/work are deferred by the owner. This supplies no evidence of successful sibling integration. No sibling synchronization or migration is authorized during this remediation.

Owner decision: 2026-10-12 — **traminex corrects defects of its R sources instead of porting them.** Where TraMineR 2.2-14 or TraMineRextras 0.6.9 behaves wrongly (a result that contradicts the method's definition or R's own documentation, an R error or quirk on valid input), traminex implements the correct behaviour. Each such correction is a deliberate change in the sense of "Source of truth" below: it is recorded as a `DIVERGENCE` row of `docs/API_PARITY_DECISIONS.csv` that states R's behaviour, the correction and the reason, and it is tested. R remains the source of truth for all other behaviour. weightedcluster and flexscanx stay on standby until traminex reaches a stable, release-ready version; this decision authorizes no sibling edit, synchronization or migration, and does not resume R6.

Owner decision: 2026-10-19 — **traminex's default alphabet is R's order.** Every traminex constructor that derives a state alphabet without `alphabet=` (`seqdef`, `seqformat`, `seqdef_from_long`, `TSE_to_STS`, `HSPELL_to_STS`, `events_to_state_sequences`, `seqMD`'s combined alphabet) orders the real states as TraMineR's `seqdef` does through `seqstatl`: by numeric value when every value (a `missing=` code included) is an R number, R's ICU root collation otherwise (ICU stays a hard requirement, decision of 2026-10-12). An explicit `alphabet=` keeps the caller's order, as in R. This reverses the earlier first-appearance default, so per-state results, colours by position, tie rules and state-rank defaults follow R without an explicit alphabet.

Owner decision: 2026-10-26 — **the family supervision plan is approved.** The supervision repository, the four projects' parent folder `PROJECT_SEQANALYSIS_and_SPATSCAN` (its `docs/PLAN.md` and `docs/DECISIONS.md`) supervises the four projects; the standby of weightedcluster and flexscanx is lifted for that plan's work only. Its records cleanup (owner decisions D2 and D3) removed the audit run records; their findings are in `docs/FINDINGS_REGISTER.csv`.

Owner decision: 2026-10-26 — **the supervision repository is the family's supervisor-coordinator, a standing role.** A session started there reads the four projects, runs their documented checks and suites, writes dated findings into each project's `docs/SUPERVISOR_NOTES.md` without asking, and coordinates cross-project work with an ownership table, one project part per row. It changes a project's code, tests, documents, configuration or environment only in a task the owner approved, and it hands no work to other sessions or subagents. Each project's own agent does that project's work, reads an open supervisor entry before working on the area it names, and records the resolution there. Pins, releases, versions, remotes, deletions, backups and conflicts between projects are the owner's decisions. The full text is in the supervision repository's `docs/CHARTER.md`, "Role, duties and authority".

Last reviewed: 2026-07-15 (postponed traminex R6/1.0 while the full traminex
codebase and project are audited and prioritized work continues elsewhere in
the ecosystem). Targeted edit, not a full review: 2026-09-20 (moved the
traminex oracle to TraMineR 2.2-14; TraMineRextras stays 0.6.9). Earlier:
2026-07-12 (added `py_ggplotseq` as the dedicated visualization project,
recorded Ubuntu 24.04 as the sole development platform, upgraded the
traminex oracle to TraMineR 2.2-13 / TraMineRextras 0.6.9, and recorded that
flexscanx should consume more of traminex after a later
architecture decision). Earlier: 2026-06-15 (clarified the source-of-truth
rule: **rflexscan is flexscanx's *starting point*, not its source of truth** —
flexscanx is a new implementation that deliberately corrects
rflexscan/FleXScan's algorithmic and
mathematical errors and expands it; traminex←TraMineR and
pyweightedcluster←WeightedCluster remain faithful, R-anchored ports). Earlier:
2026-06-13 (after completing migration Phase 3 — flexscanx's `sequence/` copies
are now parity-tested identical to traminex, completing the migration; Phase 2
had moved the weighted discrepancy/tree primitives into traminex. See
`DUPLICATES.md` and `MIGRATION_PLAN.md`).

## Domain (read first — this is social science, not biology)

This family is **social-science statistical software**. "Sequence analysis" here
is the social-science methodology (Abbott 1995; Gabadinho, Ritschard, Studer &
Müller 2011; Studer & Ritschard 2016) for comparing **categorical state
trajectories over time** — e.g. a person's month-by-month employment, family, or
housing status across years. A "sequence" is a **discrete time series of states
from a small categorical alphabet**, **not** a DNA/RNA/protein sequence.

Some terms overlap with bioinformatics by historical analogy and can trip naive
content filters; to be explicit: **optimal matching (OM)**, **substitution
costs**, and **indel costs** are parameters of an *edit distance between two
life-course trajectories* (borrowed by analogy from string alignment), operating
on social states, not nucleotides. **TraMineR** (ported by `traminex`) is the
canonical R package for life-course sequence analysis; **WeightedCluster**
(`pyweightedcluster`) clusters trajectories; **rflexscan** is the R package for
spatial-epidemiology scan statistics that `flexscanx` *began from* — but
`flexscanx` is a **new, corrected implementation** of that method, not a faithful
port (see "Source of truth"). Nothing here is biological, genomic, or otherwise
sensitive.

## Source of truth (it differs by project — flexscanx is the deliberate exception)

**traminex, pyweightedcluster, and py_seqplot are faithful ports; their R
originals are the source of truth.** TraMineR 2.2-14 and TraMineRextras 0.6.9
(vendored at `PROJECT_TRAMINEX/SOURCES`) are
the source of truth for `traminex`. WeightedCluster 2.0 is pyweightedcluster's
canonical R source:
its readable tree is `PROJECT_PY_WEIGHTEDCLUSTER/orig_R_pkg/WeightedCluster`,
and its exact CRAN artifact is
`PROJECT_PY_WEIGHTEDCLUSTER/orig_R_pkg/WeightedCluster_2.0.tar.gz` (SHA-256
`6a03fbbc090ac92498a410faab1aeb88d686ff7a81dcfdcb97edacda44f393cc`).
The exact live companion oracle is TraMineR 2.2-14. The prior
WeightedCluster 1.8-1 source is historical evidence at
`PROJECT_PY_WEIGHTEDCLUSTER/orig_R_pkg/WeightedCluster-1.8-1`, with its exact
archive at `PROJECT_PY_WEIGHTEDCLUSTER/docs/legacy/WeightedCluster_1.8-1.tar.gz`
(SHA-256
`daf9987cee5ac9e706a39e72f54a64e49e9fe91b0bc3cfb5e5a94d5076b2df5b`).
When a Python copy disagrees with its R original, the R
original wins — **except** for things we deliberately add or change, which must be
written down as deliberate. When two Python copies disagree, the one that matches
the R source wins; "use traminex's version" is shorthand for "use the version
that faithfully matches the R source" (see `DUPLICATES.md`).

`py_seqplot` is anchored to the vendored ggseqplot source in its own repository.
It consumes traminex's verified sequence objects and analytical series and owns
the plotnine/grammar-of-graphics rendering layer; it must not fork the underlying
sequence numerics.

**flexscanx is NOT a faithful port — rflexscan is its *starting point*, not its
source of truth.** The R package **rflexscan** and the **FleXScan** engine it
wraps (vendored at `PROJECT_FLEXSCANX/RFLEXSCAN_ORIGINAL`) are where flexscanx
*began*, but flexscanx is a **new implementation that deliberately corrects
critical algorithmic and mathematical errors in rflexscan/FleXScan and expands
its capabilities.** The R versions contain genuine bugs that the Python
implementation **must fix**. Consequences:

- **Correctness is judged against the underlying mathematical/statistical
  definitions** of the scan statistics (and flexscanx's own first-principles
  tests), **not** against rflexscan's output. rflexscan is **not a parity
  oracle.**
- **Agreement with rflexscan is a cross-check only where rflexscan is itself
  correct.** Where flexscanx and rflexscan **diverge**, analyse the maths: if
  rflexscan is wrong, flexscanx's behaviour is the **deliberate, documented
  correction**, and flexscanx is **never** "fixed" to reproduce an rflexscan bug.
- **flexscanx may exceed rflexscan's surface** (new models, options, performance);
  rflexscan's feature set is a floor we started from, not a ceiling.
- **One exception inside flexscanx:** its `sequence/` sub-package's source of
  truth **is** traminex/TraMineR (faithful, test-only parity — `DUPLICATES.md`
  Section B), because rflexscan has no sequence-analysis code to begin with.

---

## 1. The four projects and where they live

These four projects form the Python sequence-analysis family. **traminex,
pyweightedcluster, and py_seqplot are faithful Python ports** of their R originals;
**flexscanx is a new, corrected reimplementation** that starts from rflexscan but
fixes its algorithmic/mathematical errors (see "Source of truth"). In R the
originals depend on each other; the Python projects are meant to do the same, but
were each built independently, so the links are not yet visible in the code. This
charter records the intended relationships so that anyone — human or AI assistant
— working in any one folder understands the whole family.

| Project | Folder | Python package | R package & relationship |
|---|---|---|---|
| **traminex** | `PROJECT_SEQANALYSIS_and_SPATSCAN/PROJECT_TRAMINEX` | `traminex` | TraMineR / TraMineRextras (faithful port) |
| **pyweightedcluster** | `PROJECT_SEQANALYSIS_and_SPATSCAN/PROJECT_PY_WEIGHTEDCLUSTER` | `weightedcluster` | WeightedCluster (faithful port) |
| **py_seqplot** | `PROJECT_SEQANALYSIS_and_SPATSCAN/PROJECT_PY_SEQPLOT` | `py_seqplot` | ggseqplot (faithful plotnine port) |
| **flexscanx** | `PROJECT_SEQANALYSIS_and_SPATSCAN/PROJECT_FLEXSCANX` | `flexscanx` | rflexscan (**starting point** — corrected & expanded, **not** a port) |

All active development targets Ubuntu 24.04. Windows paths, Windows Store
Python, and PowerShell-based setup instructions are historical only.

### Current pre-release sequencing decision (2026-07-15)

Traminex R6 (the 1.0 release phase) is **postponed**. Before R6 may be resumed:

1. the full traminex codebase and project must receive a separate, comprehensive
   audit;
2. audit findings must be triaged and any release-blocking remediation must be
   completed or explicitly dispositioned by the maintainer;
3. prioritized work on other components of this ecosystem will proceed under
   those projects' own plans; and
4. the R6 plan must be refreshed against the resulting ecosystem state and the
   maintainer must explicitly authorize its execution.

This sequencing decision does not authorize sibling-repository edits,
cross-project migrations, a traminex version change, tagging, or publication.
On 2026-11-07 the maintainer decided to release below 1.0: traminex 0.2.0, with
the same 226-name public API; on 2026-11-08, the patch release 0.2.1, with that API and
the fixes of 2026-11-07, then the patch release 0.2.2, the same code without this machine's paths in
the published files. The 1.0 release (R6) stays postponed until it is
explicitly resumed.

## 2. Who depends on whom (the layering)

```
                         traminex
          (sequences, distances, costs, trees, stats)
                      /     |      \
                     /      |       \
       pyweightedcluster  py_seqplot  flexscanx
       (clustering)       (plotnine)    (spatial scans)
```

- **traminex is the base.** It owns sequence analysis: building sequences,
  measuring the distance between them, substitution/indel costs, sequence trees,
  and sequence statistics.
- **pyweightedcluster builds on traminex.** Its real job is grouping similar
  sequences (clustering) and judging cluster quality. The distances it clusters
  are traminex's department.
- **py_seqplot builds on traminex.** It is the dedicated user-facing
  grammar-of-graphics visualization library. traminex owns the numerical series
  and the matplotlib renderers that re-express the TraMineR plots (`seqplot` and its
  types, `seqpcplot`, `seqsplot`, the tree DOT output; see docs/PLOTTING.md);
  py_seqplot owns faithful
  ggseqplot-style plotnine renderers and no sequence-analysis numerics.
- **flexscanx builds on both.** Its real job is finding spatial clusters
  (the scan statistic). When it tracks how those clusters evolve over time it
  uses sequence tools (traminex's domain) and may use clustering tools
  (pyweightedcluster's domain).

  **flexscanx dependency rule (decided 2026-06-12):** flexscanx must install and
  run **without** traminex or pyweightedcluster installed. So "builds on" here
  means *follows their definitions*, not *imports them at runtime*. flexscanx
  keeps its own small copies of the shared sequence tools, but those copies must
  stay **identical in behaviour** to traminex's, checked by parity tests that use
  traminex only at test time. The truth flows traminex → flexscanx, never the
  reverse.

  **Planned direction (2026-07-12):** flexscanx should use more of traminex than
  the current small parity-tested surface. The exact runtime/dependency boundary
  will be decided later. Until that discussion occurs, preserve the current
  standalone-install guarantee and do not initiate a migration unilaterally.

  **Note on "builds on pyweightedcluster":** as audited, this link is **latent,
  not active.** flexscanx's own clustering is *spatial* (which areas form a hot
  cluster), which is a different thing from WeightedCluster's clustering of
  sequences; flexscanx does not currently duplicate or need WeightedCluster's
  algorithms. Do not add a WeightedCluster dependency just to match this diagram.

## 3. Ownership map — who is the official home for each kind of function

When you need a function, get it from its owner below. Do **not** write a new
copy in a project that is not the owner.

| Area | Owner | Examples |
|---|---|---|
| Building & formatting sequences | **traminex** | `seqdef`, `seqformat`, `STSList` |
| Distances between sequences | **traminex** | `seqdist`, OM / OMspell / DHD / LCS / Hamming, etc. |
| Substitution & indel costs | **traminex** | `seqsubm`, `costs.*`, transition-rate costs |
| Sequence statistics | **traminex** | `seqstatd`, `seqtransn`, entropy, turbulence |
| Discrepancy analysis & trees | **traminex** | `dissvar`, `dissassoc`, `disscenter`, `disstree`, `seqtree` (R-faithful, weighted); `seqtree_fit`/`disstree_fit` are exploratory sklearn variants |
| Analytical plot series and matplotlib re-expressions of the TraMineR plots | **traminex** | `seqstatd`, `seqrf`, `seqplot` and its types, `seqpcplot`, `seqsplot` |
| Grammar-of-graphics sequence visualization | **py_seqplot** | `ggseqdplot`, `ggseqiplot`, `ggseqrfplot`, all `ggseq*` renderers |
| Clustering of sequences | **pyweightedcluster** | k-medoids (`wcKMedRange`), PAM, CLARA |
| Cluster quality & validation | **pyweightedcluster** | `wcClusterQuality`, bootstrap, fuzzy |
| Null models for cluster testing | **pyweightedcluster** | `seqnull` |
| Spatial scan statistics | **flexscanx** | `scan`, `scan_series`, flexible/circular scan |
| Temporal patterns of spatial clusters | **flexscanx** | `extract_sequences` from scan results |

## 4. The working rule (read before adding sequence code)

> **Before implementing anything sequence-related (a distance, a cost matrix, a
> sequence statistic, a sequence tree) in pyweightedcluster or flexscanx, check
> whether traminex already provides it. If it does, use traminex. If you must
> keep a local copy, write down why in this project's `DUPLICATES.md`.**

Why this rule exists: traminex's distance functions are tested for exact
agreement with the original R package (it has R "golden" parity tests). The
copies in the other two projects are not. Two independent copies of the same
formula will drift apart over time and give different numbers for the same data,
which is the bug this family is trying to avoid.

A deliberate local copy is sometimes the right call (for example, a tiny
self-contained helper where depending on a whole extra package is not worth it).
That is allowed — but it must be a recorded decision, not an accident.

**"Best implementation" is decided by the R source, not by which project it sits
in.** The audit found that for sequence *distances and costs*, traminex has the
better, R-validated, weight-aware version (the siblings should use it). The
*weighted discrepancy and tree* functions (`dissvar`, `dissassoc`, `disscenter`,
`disstree`/`seqtree`) used to be simplified, unweighted stand-ins in traminex
with the faithful weighted maths only inline in pyweightedcluster — **this was
fixed in migration Phase 2 (2026-06-13)**: traminex now hosts the R-faithful,
weighted, R-golden-tested versions, and pyweightedcluster's inline sites call
them. See `DUPLICATES.md` ("Phase 2 status: COMPLETE") for the per-function
record.

**Name-clash warning (RESOLVED 2026-06-13, Phase 2 / A7):** historically,
traminex's `disscenter` double-centred a matrix while R's `disscenter` is
"distance of each case to its group centre". As of traminex 0.0.6,
`disscenter` has the R meaning; the double-centering helper is
`double_center_matrix` (deprecated alias `disscenter_double_center`).

## 5. Current known duplication

The full, function-by-function inventory with a recommendation for each item is
in **`DUPLICATES.md`** (master copy in `PROJECT_TRAMINEX`), where the
distance / cost / discrepancy duplication is now recorded as **migrated**
(Phases 1–2 COMPLETE). Short version:

- **pyweightedcluster** *formerly* carried its own sequence-distance stack;
  migration Phases 1–2 removed it. `_omspell.py` and `costs.py` are deleted and
  `_seqdist.py` is now a thin shim delegating HAM/OM/LCS/DHD/OMSPELL to
  `traminex.distances.pairwise`, with the weighted-inertia maths taken from
  `traminex.diss`. What remains is genuinely WeightedCluster's own: the compiled
  clustering / quality kernel (`_wc_core`), `seqnull.py`, and `seqtree.py`'s
  `as_seqtree` (which now consumes the traminex primitives).
- **flexscanx** carries a small `sequence/` sub-package (Hamming distance,
  transition matrices, entropy, index/distribution/frequency plots). Its own
  `NOTICE.txt` already says these were adapted from traminex.
- **py_seqplot** consumes traminex's public objects and series. It should not
  duplicate distance, statistics, representative, or discrepancy calculations.

## 6. How an AI assistant should use this charter

If you are an AI assistant working in any of these four folders:

1. You learned the family relationships from this file and from the project's
   `CLAUDE.md` / `AGENTS.md`. The sibling projects' absolute paths are in
   section 1 — you can open and read their source directly.
2. Respect the ownership map (section 3) and the working rule (section 4).
3. Do **not** start replacing existing duplicated code with cross-project
   imports on your own initiative. That is a planned, staged migration that
   changes dependencies and needs the maintainer's go-ahead. Propose it; do not
   perform it unprompted.

## 7. Keeping this file in sync

The master copy is `PROJECT_TRAMINEX/ECOSYSTEM.md`, and there is nothing to copy:
since 2026-10-26 each sibling's `ECOSYSTEM.md` points here and keeps only what is
specific to that project (its source of truth, and the statements its own checks
read). Copies made with `cp` had drifted by 69 to 72 lines by 2026-10-25. Edit
this master only. The family's supervision layer (dated decisions, open items,
inventory, the family check) lives in the supervision repository's
`docs/CHARTER.md`.
