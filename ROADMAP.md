# Traminex Roadmap and Issue List

This document tracks parity targets with TraMineR/TraMineRextras, performance goals, and UX improvements.

## Full-parity and reliability closeout (R5 complete; R6 postponed 2026-07-15)

The active scope is derived from
[`docs/API_PARITY_DECISIONS.csv`](docs/API_PARITY_DECISIONS.csv), not from the
historical priority lists below. Its deterministic views are
[`docs/API_PARITY_LEDGER.md`](docs/API_PARITY_LEDGER.md) and
[`reports/audit/traminex_semantic_parity.json`](reports/audit/traminex_semantic_parity.json).

The generated semantic JSON publishes the current namespace, public API,
status, scope, ownership, option-coverage, and target-phase totals. It currently
has no `PLANNED` rows. Derive all counts with the command below rather than
copying fast-changing summaries into this roadmap.

R1 added exactly `seqcost`, `seqeapplysub`, and `seqsep` to the R0 surface of
208 names. R3 added exactly `createdatadiscrete`, `group_p`,
`seqerulesdisc`, `seqgen_missing`, `seqsamm`, `seqsammeha`, `seqsammseq`,
`sorti`, and `sortv` to the 211-name R2 surface. Existing public symbols were
not removed. R4 then added exactly `plot_SAMM`, `plot_seqimplic`,
`seqplot_rf`, `seqplot_tentrop`, `seqplot_tentrop_m`, and `seqsplot` to the
220-name R3 surface. The existing `seqplot` route was extended in place.

| Phase | Remaining planned rows | Analytical rows | R namespace rows |
| --- | ---: | ---: | --- |
| R1 — core contracts (complete) | 0 | 0 | Closed: `alphabet`, `seqcost`, `seqdef`, `seqdist`, `seqeapplysub`, `seqformat`, `seqindic`, `seqivardur`, `seqlength`, `seqsep`, `seqtrate`, `FCE_to_TSE`, `TSE_to_STS`, `toPersonPeriod` |
| R2 — existing analytical branches (complete) | 0 | 0 | Closed: `dissassoc`, `dissdomassoc`, `dissrf`, `disstree2dot`, `seqdomassoc`, `seqient`, `seqfpos`, `seqfposend`, `seqibad`, `seqidegrad`, `seqinsecurity`, `seqistatd`, `seqlogp`, `seqMD`, `seqprecarity`, `seqrep`, `seqrf`, `seqstatd`, `seqpolyads`, `seqsurv` |
| R3 — extras completion (complete) | 0 | 0 | Closed: `createdatadiscrete`, `dissCompare`, `group.p`, `seqBIC`, `seqCompare`, `seqecreate`, `seqerulesdisc`, `seqgen.missing`, `seqLRT`, `seqsamm`, `seqsammeha`, `seqsammseq`, `sorti`, `sortv` |
| R4 — plotting parity (complete) | 0 | 0 | Closed: `plot.stslist`, `plot.SAMM`, `plot.seqimplic`, `seqplot.rf`, `seqplot.tentrop`, `seqplot.tentrop.m`, `seqsplot` |

R5 reliability and R6 release work are quality phases, not additional R
namespace decision rows. R1 was completed through the executable
`2026-07-12-traminex-r1-core-contracts.md`
plan. R2 was completed through
`2026-07-12-traminex-r2-existing-branches.md`.
R3 was completed through
`2026-07-13-traminex-r3-extras-completion.md`.
R4 was completed through
`2026-07-15-traminex-r4-visualization-parity.md`.
R5 completed Ubuntu-native support and build gates, a CI workflow configuration
(never run: the repository has no git remote), whole-package typing and
bidirectional stubs, honest branch coverage, consolidated benchmarks, and six
isolated artifact-install profiles through
`2026-07-15-traminex-r5-ubuntu-reliability.md`.
The source-first R6/1.0.0 release plan was written, but its execution is postponed by
maintainer decision:
`2026-07-15-traminex-r6-release-1.0.md`. These six plans were removed from the tree on
2026-10-26 with the other run records (owner decision D3 of the family supervision plan)
and remain in git history, for example
`git show b03ad054^:docs/superpowers/plans/2026-07-15-traminex-r6-release-1.0.md`.

The first full traminex codebase audit completed on 2026-08-13. Further audits
followed, the latest on 2026-10-19; the repair of their confirmed findings is in
progress, and remediation is complete only when it passes final verification.
The remaining ecosystem interval and any R6 refresh are separate, pending work;
R6 still requires a new, explicit maintainer instruction to execute. Neither
the audit nor remediation authorizes a version bump, tag, publication, sibling
migration, or expansion of traminex's frozen public surface.

| Quality phase | Status | Public/version effect |
| --- | --- | --- |
| R5 — Ubuntu reliability and release hardening | Complete 2026-07-15 | None: 226 names, version 0.0.8 |
| Pre-R6 full audit and remediation | First audit completed 2026-08-13, latest 2026-10-19; repair in progress; remediation completion requires final verification | Ecosystem work remains separately pending; no automatic API/version change |
| R6 — Release 1.0 | Postponed by maintainer 2026-07-15 | No execution until prerequisites and explicit resumption; eventual API remains 226 |

Reproduce these values from the tracked JSON:

```bash
PY=~/miniforge3/envs/datasci/bin/python
$PY scripts/generate_semantic_parity_ledger.py --check
$PY - <<'PY'
import json
from pathlib import Path

p = json.loads(Path("reports/audit/traminex_semantic_parity.json").read_text())
print(p["summaries"]["target_phase"])
print(p["summaries"])
print(p["target_public_api"])
PY
```

## Historical completed roadmap

Priorities: P0 (now), P1 (next minor), P2 (later)

P0 – Parity correctness and immediate performance ✅ COMPLETE
- Stats/Complexity
  - ✅ seqST (exact turbulence) parity with TraMineR
    - Status: COMPLETE - Both type=1 and type=2 implemented
    - type=1 (default): variance over observed runs only
    - type=2: variance denominator includes non-visited states (alphabet size K)
  - ✅ turbulence_exact aggregate
    - Status: COMPLETE - Weighted mean aggregation with full parameter support
- Distances
  - ✅ RLCP accelerator
    - Status: COMPLETE - Numba path verified; historical machine-specific timings are archived
    - Python fallback preserved for systems without Numba
  - ✅ DHD accelerator
    - Status: COMPLETE - Numba path verified; historical machine-specific timings are archived
    - Calibrated profile respected, both mean and sum variants
- Tests/Docs
  - ✅ Golden tests for seqST behavior and normalization
    - Status: COMPLETE - 7 comprehensive tests, all passing
  - ✅ Benchmark infrastructure created
    - `benchmarks/run_benchmarks.py` with agreement-before-timing verification
  - ✅ Documentation updated
    - Detailed completion notes archived under archive/status_reports/

P1 – Feature breadth and robustness ✅ COMPLETE (6/6 items)
- Distances/OM family
  - ✅ Position-dependent substitution matrix support (subm_by_pos [L,K,K])
    - Status: COMPLETE - Verified existing infrastructure, Python kernel supports [L,K,K] arrays
  - ✅ OMloc wrapper and docs (mapping to indel_by_pos/local profiles)
    - Status: COMPLETE - Full implementation with 11 tests, all passing
    - File: traminex/distances/omloc.py (197 lines)
    - Supports multiple indel methods (change, entropy) and substitution modes
    - Canonical public API: omloc
- Events
  - ✅ Extend seqefsub with support/confidence and simple rule generation
    - Status: COMPLETE - Added as_percent parameter, new seqerules function
    - File: traminex/events.py, seqerules() generates association rules with support, confidence, lift
      (since RP5a/RP5b, 2026-06-19: the `traminex/events/` package, `_mining.py`)
    - Tests: 14 comprehensive tests in test_seqerules.py, all passing (since RP5b,
      superseded by the live-R goldens in `tests/test_r_goldens_seqerules.py`)
  - ✅ Gap/window constraints (min_gap, max_gap) and id_weights
    - Status: COMPLETE - Full min_gap/max_gap support with gap_in_time option, id_weights for weighted mining
    - File: traminex/events.py, enhanced SeqEConstraint, seqefsub, seqerules
      (since RP5a/RP5b, 2026-06-19: the `traminex/events/` package, `_mining.py`)
    - Tests: 13 comprehensive tests in test_gap_constraints.py, all passing (since
      RP5b, superseded by the live-R goldens in `tests/test_r_goldens_seqefsub.py`)
- Stats
  - ✅ Exact turbulence type=2 (include non-visited states)
    - Status: COMPLETE (promoted from P1, completed in P0)
  - ✅ seqdiff: per-position state distribution differences between groups
    - Status: COMPLETE - Multiple methods (residuals, chi2, log_odds), bootstrapped CIs, p-value corrections
    - File: traminex/seqdiff.py (358 lines), 15 comprehensive tests all passing
    - Supports 3 statistical methods, bootstrap inference, and multiple testing corrections
- Multi-domain
  - ✅ Domain weight search utilities (grid over weights with silhouette target)
    - Status: COMPLETE - Full grid search and coordinate descent optimization
    - File: traminex/multidomain.py (seqMD_gridsearch_weights, seqMD_optimize_weights)
    - Tests: 21 comprehensive tests in test_multidomain_weights.py, all passing
    - Silhouette-based quality metric with extensible framework

P2 – Performance and developer UX ✅ COMPLETE (4/4 items)
- Performance
  - ✅ Numba paths for OM (scalar/matrix), LCS/LCP, CHI2, JSD_TRANS
    - Status: COMPLETE - All core distance methods now have Numba accelerators
    - File: traminex/distances/python_kernels.py
    - 10 new Numba functions: OM (scalar/matrix), LCS, LCP, CHI2, JSD_TRANS
    - Use the maintained benchmark CLI to measure the relevant workload and environment
    - Automatic dispatch, graceful fallback to Python
  - ✅ Optional pybind11 extension for RLCP/DHD
    - Status: COMPLETE - Added compiled RLCP and DHD dispatch paths when pybind11 is available
    - Python/Numba fallback remains automatic when compiled backends are absent
- Developer ergonomics
  - ✅ Config dataclasses for TRATE/DHD profiles
    - Status: COMPLETE - Type-safe configuration objects for complex parameters
    - File: traminex/config_dataclasses.py (296 lines)
    - 5 dataclasses: OMConfig, TRATEConfig, DHDConfig, LCSLCPConfig, DURConfig
    - Validation, composition, conversion to kwargs
    - Exported in main API: from traminex import OMConfig, etc.
  - ✅ Stronger ndarray shape/dtype hints in core distance docstrings
    - Status: COMPLETE - Added explicit input/output shape and dtype guidance
    - Files: traminex/distances/__init__.py, traminex/distances/omloc.py
- Plotting
  - ✅ Missing/void styling in index/distribution plots (alpha/hatching)
    - Status: COMPLETE - Alpha transparency and optional hatching for special states
    - Files: traminex/visualization/plots/index.py, distribution.py
    - New parameters: missing_alpha, void_alpha, apply_hatching
    - Backward compatible (opt-in)
  - ✅ Save helpers across plot APIs
    - Status: COMPLETE - Unified save_plot() function
    - File: traminex/visualization/utils.py
    - Exported: traminex.plotting.save_plot, traminex.visualization.save_plot
    - Handles (fig, ax) tuples, auto-creates directories, infers format

Distance-metrics parity ✅ COMPLETE (v0.0.2 / v0.0.3)
- Distances/new methods
  - ✅ EUCLID, R-faithful CHI2 (default step=1), TWED, OMslen, OMspell, OMstran,
    NMS, NMSMST, SVRspell — default and risk-weighted option routes are golden-pinned
    against live TraMineR 2.2-14 (tests/test_r_goldens_distances.py,
    tests/test_nms_parity.py)
  - ✅ Substitution costs: seqsubm FEATURES (Gower over state features) and
    FUTURE (chi-squared on common-future transition rates); gower_distance()
- OM-family normalization (v0.0.3)
  - ✅ norm= for TWED/OMslen/OMspell/OMstran to live-R parity
    (maxlength/gmean/maxdist/YujianBo/auto), refused for NMS/NMSMST/SVRspell —
    mirrors DistanceCalculator::normalizeDistance (tests/test_om_normalization.py)
- Parity remediation (closeout audit 2026-06-12)
  - ✅ DHD reconciled to TraMineR's transition-rate Dynamic Hamming
    (time-varying TRATE, cval=4); golden-pinned against live R
  - ✅ Version identity 0.0.3 synchronized; consistency gate green

Tracking
- Link issues to this roadmap and update statuses as PRs merge.
- Documentation hygiene and archival policy live in docs/DOCS_POLICY.md.
