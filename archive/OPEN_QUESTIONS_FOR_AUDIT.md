# Open Questions for the Independent Audit

Prepared during a tidy/index/graph pass on branch `feat/distance-metrics-parity`
(2026-06-11). This pass organized and documented what the code **is**; it did
**not** verify correctness. The items below are the claims that need an
independent check. R is available for goldens at:

    %USERPROFILE%\miniforge3\envs\datasciencex\Scripts\Rscript

Vendored R ground truth: `SOURCES/TraMineR` (2.2-12),
`SOURCES/TraMineRextras` (0.6.8). Do not modify these trees.

---

## 1. CHI2 default-behavior change (R-faithful `step=1`) — needs a live-R golden

- **What changed.** The new CHI2 implementation reproduces the
  Studer–Ritschard chi-squared and sets the **default `step=1`** (one interval
  per position), described as "matching TraMineR" in
  `docs/NEW_DISTANCE_METHODS.md`. This is a behavior change from any prior
  CHI2/aggregate-interval expectation.
- **Why it is only weakly guarded.** The parity test
  `tests/test_chi2_euclid_parity.py` validates the algebra against a Python
  re-implementation, and its reference cases use a **single interval
  (`step = L`)**, not the runtime default `step=1`. So the default path is
  covered by property/algebraic tests, **not** by a TraMineR golden.
- **Check.** Generate `seqdist(..., method="CHI2")` goldens from R TraMineR
  (`seqdist-CHI2.R`, `chisq.cpp`) at the default settings and compare to
  `traminex.seqdist(x, method="CHI2")` with no `step` argument. Confirm the
  `step=1` default actually matches TraMineR's default, and that `norm` (a
  boolean here, distinct from the OM-style `norm="length"/"max"`) matches.
- Pointers: `traminex/distances/python_kernels.py` (`_chi2_*`),
  `traminex/distances/_validation.py` (CHI2 allowed kwargs).

## 2. SVRspell / NMSMST — validated by brute force only, no R at build time

- **What was done.** `tests/test_nmsmst_svrspell_parity.py` validates the
  C++/Python ports (`NMSMSTdistance.cpp`, `NMSDURSoftdistance.cpp`) against an
  **independent brute-force embedding enumeration**, explicitly because R was
  not available when the methods were written.
- **Check (now that Rscript exists).** Produce R goldens for `NMSMST` and
  `SVRspell` (and `NMS`, `NMSdistance.cpp`) from `SOURCES/TraMineR` and compare
  to `traminex.seqdist(..., method="NMSMST"|"SVRspell"|"NMS")`. Brute force
  confirms the spec the author *implemented*; it does not confirm that spec
  equals TraMineR's, e.g. soft-match (`prox`), `kweights`, and `tpow`
  conventions.
- Pointers: `traminex/distances/python_kernels.py` (NMS/NMSMST/SVRspell
  branches), `docs/NEW_DISTANCE_METHODS.md` (R-reference column).

## 3. OM-family normalization contradiction — which document is true?

- **The contradiction.** `docs/NEW_DISTANCE_METHODS.md` states that the
  spell/subsequence methods **return raw, un-normalized distances**
  ("OM-family normalisation … can be layered later"), and the code backs this:
  `traminex/distances/_validation.py` does **not** list `norm` among the
  allowed kwargs for `TWED`, `OMSLEN`, `OMSPELL`, `OMSTRAN`, `NMS`, `NMSMST`,
  `SVRSPELL`. Meanwhile `docs/DISTANCES_NORMALIZATION.md` documents
  `norm=length|max|None` **only for the classic OM/LCS family** and says
  nothing about the new methods.
- **Status after this pass.** This is most likely **not** a true contradiction
  but a scoping gap: DISTANCES_NORMALIZATION predates the new methods and only
  ever covered the classic family. During this pass `docs/INDEX.md` and
  `docs/README.md` were annotated to scope DISTANCES_NORMALIZATION to the
  classic family and to state that the new methods are raw. The knowledge graph
  (`graphify-out/GRAPH_REPORT.md`) confirms the "raw un-normalized" grouping as
  an explicit hyperedge over exactly those seven methods.
- **Check.** Confirm that (a) the new methods genuinely expose no normalization
  at runtime (try `seqdist(x, method="TWED", norm="length")` → should raise),
  and (b) decide whether leaving them raw is acceptable for release or whether
  TraMineR parity requires `maxdist`/`YujianBo` normalization. If normalization
  is required, DISTANCES_NORMALIZATION must be extended, not just scoped.
- **RESOLVED (branch `feat/om-normalization`).** Reading the R source confirmed
  TraMineR *does* accept `norm` for OM/OMloc/OMstran/OMspell/OMslen/TWED but
  *refuses* it for NMS/NMSMST/SVRspell (`seqdist.R` line 205). Normalization was
  therefore genuinely missing for TWED/OMslen/OMspell/OMstran and is now
  implemented to live-R parity (`tests/test_om_normalization.py`,
  `traminex/distances/_normalization.py`), with `norm` refused for
  NMS/NMSMST/SVRspell. `docs/DISTANCES_NORMALIZATION.md/.qmd` were extended to
  cover the new methods (not just scoped).

## 4. The feature branch is unmerged and ahead of `master`

- `feat/distance-metrics-parity` is **not merged**. As of this pass it carries
  the five original distance-method commits **plus** the commits added during
  this tidy pass (audit-remediation layer split into phase-1 / phase-3 /
  generators / checks, one **WIP** CI commit, repo cleanup, doc unification,
  and the graph). It deliberately remains unmerged pending this audit.
- New distance methods on the branch: Gower/FEATURES/FUTURE, EUCLID,
  R-faithful CHI2, TWED, OMslen, OMspell, OMstran, NMS, NMSMST, SVRspell.
- **Version note (do not act on without audit sign-off).** `pyproject.toml`,
  `traminex/__init__.py`, README, RELEASE_NOTES, and the CHANGELOG `[0.0.1]`
  section all agree on **0.0.1**; the CHANGELOG `[Unreleased]` section holds P1
  additions. The substantial new distance surface on this branch **likely
  warrants a future version bump**, but that is the audit/maintainer's call —
  it was intentionally **not** bumped here.

---

## Secondary items surfaced during this pass (lower priority)

- **WIP CI job is intentionally red.** The `consistency` job added to
  `.github/workflows/ci.yml` (commit subject begins `WIP(ci): … DO NOT ENABLE
  AS-IS`) runs a bare `python -m mypy` that follows imports into unannotated
  modules and fails (36 errors across 13 files, mostly pre-existing in
  `traminex/distances/*`, `stats_p1/p2`, `visualization/*`). The repo-blessed
  type path `scripts/typecheck.ps1` (`mypy --follow-imports=skip`) is green.
  Resolution: either run mypy with `--follow-imports=skip` in CI or annotate
  the followed modules. All other steps of that job pass locally.
- **`check_stub_signatures.py` import requirement.** The new stub-signature
  gate needs `traminex` importable. The local `.venv` has an **older**
  `traminex` (0.1.0 metadata) registered but not importable as a package, so
  the script must run from the repo root or with `PYTHONPATH` set to the repo.
  CI installs `-e .`, which resolves this. Auditors running locally: set
  `PYTHONPATH` to the repo root.
- **Observed local state (observation, not endorsement):** pytest **312
  passed, 2 skipped**; `ruff check .`, `scripts/typecheck.ps1`,
  `check_docs_consistency.py`, `check_stub_signatures.py` (with PYTHONPATH),
  `check_generated_artifacts.py`, and `verify_audit_improvements.py` all pass.
  These confirm internal consistency, **not** parity with R.
