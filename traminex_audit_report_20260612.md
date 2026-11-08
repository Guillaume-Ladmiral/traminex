# Traminex Independent Closeout Audit — 2026-06-12

**Auditor:** Independent closeout audit (did not author any of the code under review).
**Repo state:** branch `master` @ `64582423`; working tree clean at audit start.
**Ground truth:** live R **4.4.3 (ucrt)** + **TraMineR 2.2.12** (compiled from
`SOURCES/TraMineR_2.2-12.tar.gz`); vendored `SOURCES/TraMineR` (2.2-12) /
`SOURCES/TraMineRextras` (0.6.8). Python **3.12.12**, `.venv` interpreter.
**Method:** trust nothing — every claim re-derived from live R / SOURCES / a clean run.

> **Scope note (mandate):** the audit itself was *report, don't fix.* A
> subsequent **remediation pass (same day, at maintainer request)** implemented
> all findings — see **§0. Remediation status** below. The findings in §8 are
> retained as the original audit record.

---

## 0. Remediation status (fixes applied 2026-06-12)

All audit findings were fixed in a follow-up pass. Fresh verification:
**`pytest` 436 passed / 1 skipped** (was 432/2); all gates green
(`check_docs_consistency`, `check_stub_signatures`, `check_generated_artifacts`,
`verify_audit_improvements`, `mypy --follow-imports=skip`, `pre-commit`).

| Finding | Fix | Evidence |
|---|---|---|
| **H1** version mismatch / CI red | Bumped `__init__.py`, `README.md`, `docs/INDEX.md` to **0.0.3**; regenerated audit artifacts | `check_docs_consistency` PASS (Version: 0.0.3); `check_generated_artifacts` PASS |
| **M1** DHD not TraMineR's | `method="DHD"` now computes TraMineR's transition-rate Dynamic Hamming (time-varying TRATE, `cval=4`); legacy marginal metric kept as opt-in (`variant=`/`smooth_window=`/`dhd_profile=`). New golden `DHD_default`. | live R match ≤ **3.6e-14** on 3 datasets incl. a zero-denominator edge case; `tests/test_r_goldens_distances.py::DHD_default` passes |
| **M1** doc error (API_PARITY:99) | Corrected during the audit (commit `9e04dc69`) | — |
| **M2** orphaned/absent release tags | Annotated tag **`v0.0.3`** created on the remediated HEAD (orphaned `v0.1.1` left for the maintainer to delete) | `git tag` |
| **L1** `test_seqtree_smoke` always-skipped | Replaced broken guard with `pytest.importorskip("sklearn")` | test now runs & passes (skip count 2→1) |
| **L2** INDEX staleness | Updated version/test-count/date | — |
| **L3** ROADMAP frozen | Added distance-metrics + OM-normalization + remediation milestones | — |
| **L5** OMslen_gmean golden gap | Pinned `gmean` + `auto` goldens | `test_om_normalization.py` 93→95 pass |
| **L6** stray `#pragma once` | Removed from `cpp/bindings/pybind11_module.cpp` | pybind11 still builds |
| **L4** nanobind build | Upstream (nanobind 2.9.2 + g++ 14); **not fixed** — documented as an external limitation | — |

The DHD default change is **BREAKING** (documented in `CHANGELOG.md [0.0.3]`):
callers relying on the old marginal-frequency DHD must pass `variant="mean"`.

---

## 1. Merge verdict — `feat/om-normalization`

**The task framed this branch as "unmerged." It is not.** As of the audited
tree it is **already fully merged** into `master`:

```
git merge-base --is-ancestor feat/om-normalization master  → YES (0 commits ahead)
d31e16d5  Merge feat/om-normalization: OM-family norm (TWED/OMslen/OMspell/OMstran)
64582423  chore(release): bump version to 0.0.3 (OM-family normalization)
```

The merge commit `d31e16d5` carries an **audit-claim** in its message ("93 tests
pass… all 23 pinned goldens reproduced bit-for-bit from live R… spot-checks
match", co-authored "Claude Fable 5"). Per mandate I treated that claim as
unverified and re-derived it independently.

**Verdict: the already-merged work is SOUND. Re-verified independently.** No
re-merge is required (nothing to merge). Concretely:

| Check | Result |
|---|---|
| `tests/test_om_normalization.py` | **93 passed** (1.45 s) |
| Goldens **non-circular** (literal R values, not `tx.seqdist` output) | **Confirmed** — `GOLDENS` are hard-coded float matrices; the test compares `tx.seqdist(...)` *against* them |
| All 23 pinned goldens reproduced from **live R** | **Confirmed** — re-ran the committed `scripts/r_goldens/make_goldens_norm.R`; `max|liveR − pinnedGolden| = 0.000e+00` over all 23 |
| `maxdist` spot-check (requested) | e.g. `OMspell_exp05_maxdist[0,1] = 0.46875` — matches live R |
| `YujianBo` spot-check (requested) | e.g. `OMstran_otto05_YujianBo[0,1] = 0.51162790697674421` — matches live R |
| `norm` **refused** for NMS/NMSMST/SVRspell | **Confirmed in live R**: `'norm' isn't implemented yet for 'method = "NMS"'` (and NMSMST/SVRspell); traminex raises `ValueError` to match (`_validation.py:31,69`) |

The transcribed transform (`traminex/distances/_normalization.py`,
`normalizeDistance` codes maxlength/gmean/maxdist/YujianBo) reproduces TraMineR's
`src/distancecalculator.h` exactly on the pinned cases.

---

## 2. Whole-repo ground truth (clean tree)

| Item | Result |
|---|---|
| **Full pytest suite** (`pytest -q`) | **432 passed, 2 skipped** (434 collected), 33 s |
| Skips | `test_parity.py` (needs `TRAMINEX_PARITY=1`); `test_seqtree_smoke.py` (sklearn guard — but see Finding L1: the guard never fires) |
| **pre-commit** (`pre-commit run --all-files`) | trailing-whitespace ✓ · end-of-file-fixer ✓ · black ✓ · ruff ✓ — **all pass** |
| **mypy** (`--follow-imports=skip`, the blessed gate) | **Success, no issues in 14 files** |
| `check_stub_signatures.py` | **PASS** (36 stubs) |
| `check_generated_artifacts.py` | **PASS** (see §5) |
| `verify_audit_improvements.py` | **PASS** |
| `check_docs_consistency.py` | **FAIL** → version mismatch (Finding **H1**) |

**Which path the tests exercise:** the **pure-Python kernel path**. No compiled
binding is installed (`traminex.distances._impl_name() → "python"`, `_pb=None`,
`_nb=None`); there is no `*.pyd` in the package. The default `pyproject.toml`
build sets both bindings `OFF`.

**C++ core build (`CMakeLists.txt`, `cpp/`):**
- **pybind11 target builds cleanly** with the conda g++ 14.3 + Ninja →
  `_pb_kernels.cp312-win_amd64.pyd` (~2 MB). Placing it and re-running the
  distance/R-golden suites with `TRAMINEX_BINDING=pybind11` →
  **34 passed**, proving the native path is R-faithful and equivalent to Python
  (HAM/OM/OM-mat/LCS/LCP/RLCP/DHD dispatch). Artifacts were then removed; tree
  left clean.
- **nanobind target fails to build** — a template error *inside* nanobind
  2.9.2's own headers (`analyze_method` incomplete type) under g++ 14.3. This is
  an upstream nanobind/compiler incompatibility, **not** a traminex defect
  (Finding L4). The nanobind path is optional and unused by default.

---

## 3. R-parity spot checks (3 random methods, none from §1)

Random draw (seeded) from the previously-implemented pool → **DHD, LCP, RLCP**.
Verified on a **fresh, independent 5×7 / 4-letter dataset** (deliberately *not*
the repo's 6×8 fixture), live R vs `tx.seqdist`:

| Method | `max|tx − liveR|` | Verdict |
|---|---|---|
| **LCP** | `0.000e+00` | ✅ **Confirmed** (exact) |
| **RLCP** | `0.000e+00` | ✅ **Confirmed** (exact) |
| **DHD** | `2.210e+01` | ❌ **Discrepancy** (Finding **H2/M1**) |

DHD detail: traminex returns values in `[0,1]` (`~0.37–0.73`); R returns raw
`~11–22`; the ratios are **not constant**, so it is a genuine algorithmic
divergence, not a scaling. traminex's `method="DHD"` computes a
**marginal-frequency-weighted Hamming** averaged over L
(`python_kernels.py:827` `1 − 0.5(p_t[x]+p_t[y])`, ÷L). TraMineR's DHD is
**Lesnard's transition-rate dynamic Hamming**: `seqdist.R:451` builds
`sm <- seqsubm(..., "TRATE", cval=4, time.varying=TRUE)` (four transition-rate
terms, summed, not averaged). No R-golden guards DHD anywhere in `tests/`.

> As corroboration of the broader surface I also re-ran the committed
> CHI2/EUCLID/NMS/NMSMST/SVRspell goldens against live R (§7, OQ#1/#2):
> **all 12 reproduce bit-for-bit.** DHD is the lone parity outlier found.

---

## 4. Verdict per claim

| # | Claim (source) | Verdict | Evidence |
|---|---|---|---|
| C1 | om-norm goldens come from live R, not traminex (non-circular) | **Confirmed** | 23/23 reproduce from live R, 0.000 diff |
| C2 | "93 tests pass on the branch" (merge msg) | **Confirmed** | `93 passed` |
| C3 | `norm` refused for NMS/NMSMST/SVRspell (seqdist.R:205) | **Confirmed** | live R errors; traminex raises |
| C4 | `norm="auto"` → YujianBo for OM-family (seqdist.R:316-322) | **Confirmed** | `*_auto` goldens == `*_YujianBo`, live R |
| C5 | CHI2 default `step=1` matches TraMineR default (OQ#1) | **Confirmed** | `CHI2_default` golden == live R; traminex matches |
| C6 | NMS/NMSMST/SVRspell match TraMineR (OQ#2) | **Confirmed** | live-R goldens reproduce; traminex matches |
| C7 | LCP / RLCP parity | **Confirmed** | fresh dataset, exact |
| C8 | `method="DHD"` ports TraMineR's DHD | **Discrepancy** | differs from live R; different algorithm |
| C9 | "extra metrics (DHD, DUR, JSD_TRANS) not present in TraMineR" (API_PARITY.md:99) | **Discrepancy** | **DHD *is* a TraMineR method**; statement is factually wrong |
| C10 | Package version is 0.0.3 (release commit) | **Discrepancy** | only pyproject+CHANGELOG bumped; code/README/INDEX still 0.0.2 |
| C11 | Generated parity artifacts are fresh | **Confirmed** | `check_generated_artifacts.py` PASS |
| C12 | pybind11/nanobind compiled backends build (`docs/DEV_BUILD.md`) | **Partly** | pybind11 ✓; nanobind ✗ (upstream) |
| C13 | INDEX.md lists every doc & method | **Confirmed** | 22/22 docs listed; 18/18 methods match `_validation.py` |

---

## 5. Regenerated-artifact comparison

`check_generated_artifacts.py` regenerates into a temp dir (via
`audit_traminex.ps1` + `generate_parity_table.py` + `generate_api_mappings.py`)
and diffs against the committed copies. **Result: PASS** — every tracked
generated artifact is up to date:

| Artifact | Status |
|---|---|
| `reports/audit/traminex_api_audit.{md,json,csv}` | ✅ fresh |
| `docs/API_PARITY_TABLE.{md,csv}` | ✅ fresh |
| `reports/parity/mapping_traminex_vs_TraMineR.csv` | ✅ fresh |
| `reports/parity/mapping_traminex_vs_TraMineRextras.csv` | ✅ fresh |

No stale-generated-artifact finding. (Note: `docs/API_PARITY.md` is
**hand-maintained**, not in this generated set — so its DHD error in C9 is not a
generator-staleness issue.)

**INDEX.md coverage:** all 22 `docs/*.md|.qmd` files are listed; the method list
(INDEX.md:46-47) equals `_validation.py`'s 18 keys exactly
(OM, HAM, WHAM, LCS, LCP, RLCP, DUR, JSD_TRANS, DHD, CHI2, EUCLID, TWED, OMslen,
OMspell, OMstran, NMS, NMSMST, SVRspell).

---

## 6. Documentation freshness

| Surface | Verdict | Specific notes |
|---|---|---|
| `README.md` | **Current** (minor stale) | Commands **verified working** (`project_python.ps1` wrapper, quick-usage snippet, tests, quality-checks). Stale: line 5 "version `0.0.2`" → should be 0.0.3 (H1). |
| `CHANGELOG.md` | **Current** | `[0.0.3] – 2026-06-12` entry is accurate. (The `[0.0.2]` "known gap" note for OM-family norm is correctly superseded by `[0.0.3]`; historical, not a contradiction.) |
| `ROADMAP.md` | **Outdated** | Frozen at P0/P1/P2 "✅ COMPLETE"; says nothing about the EUCLID/CHI2/TWED/OM-family distance work or OM-normalization (0.0.2/0.0.3). Nothing false, but incomplete (L3). |
| `WARP.md` | **Absent (removed)** | Deleted in `1d170092` ("Warp terminal no longer used"). The audit task's premise that WARP.md holds agent build/test instructions is **moot**. Those instructions now live in `README.md` + `docs/INDEX.md`, both **verified to work**. (Remaining "WARP" hits in docs refer to TWED time-**warp**.) |
| `docs/INDEX.md` | **Current + complete** (minor stale) | Doc & method coverage complete (§5). Stale: line 18 "version 0.0.2", line 125 "312 passed, 2 skipped" (now 432/2), header "Updated 2026-06-11", and it still describes the branch as the audit's open item (now merged). (L2) |
| `docs/API_PARITY.md` | **Contradictory** | Line 99: "extra metrics (DHD, DUR, JSD_TRANS) **not present in TraMineR**" — DHD *is* a TraMineR method (M1). |
| `docs/DISTANCES_NORMALIZATION.md/.qmd` | **Current** | Extended to cover TWED/OMslen/OMspell/OMstran norms; golden-pinned. |
| `docs/NEW_DISTANCE_METHODS.md`, `DISTANCES.md/.qmd`, others | **Current** | Consistent with code & tests. |
| `DEMOS/` (TraMineR, TraMineRextras) | **N/A** | Vendored R reference corpus; stable. |
| `examples/` | **Current** | All `examples/*.py` compile clean. |
| vignettes | **N/A** | None (R concept; no `vignettes/`). |

---

## 7. Version identity

The two merges (distance-metrics + OM-norm) landed a release bump that was
**applied incompletely**. Release commit `64582423` touched **only**
`CHANGELOG.md` and `pyproject.toml`.

| Location | Declared version | Status |
|---|---|---|
| `pyproject.toml:11` | **0.0.3** | ✅ canonical |
| `CHANGELOG.md` (latest section) | **0.0.3** | ✅ agrees |
| `traminex/__init__.py:15` `__version__` | **0.0.2** | ❌ **stale** |
| `README.md:5` | **0.0.2** | ❌ **stale** |
| `docs/INDEX.md:18` | **0.0.2** | ❌ **stale** |
| `RELEASE_NOTES_v0.0.1.md` | 0.0.1 | ➖ historical (the 0.0.1 notes; OK) |
| git tags | **`v0.1.1`** → `fa35229f` | ❌ **orphaned** — abandoned numbering scheme, points to an unrelated old commit; **no `v0.0.1/0.0.2/0.0.3` tags exist** |

**Consequence (elevates severity):** the repo's own `check_docs_consistency.py`
fails on the `__init__.py` vs `pyproject.toml` mismatch, and **CI's
`consistency` job runs that script** (`.github/workflows/ci.yml:31`) — so **CI is
red on `master`**. The mypy/stub/generated gates are all green; the version
mismatch is the sole blocker.

**Recommended fix (maintainer, not auditor):** bump `traminex/__init__.py`,
`README.md:5`, `docs/INDEX.md:18` to `0.0.3`; tag `v0.0.3` on HEAD; consider
deleting/retargeting the orphaned `v0.1.1`.

---

## 8. Findings by severity

### HIGH
- **H1 — Version identity broken; CI red on `master`.** `pyproject`/`CHANGELOG`
  say `0.0.3`; `traminex/__init__.py`, `README.md`, `docs/INDEX.md` say `0.0.2`.
  `check_docs_consistency.py` fails and is wired into CI → master CI fails.
  (Release commit `64582423` bumped only 2 of the version sites.) See §7.

### MEDIUM
- **M1 — DHD is not TraMineR's DHD, and the doc misstates it.** `method="DHD"`
  returns a marginal-frequency Hamming, not Lesnard's transition-rate dynamic
  Hamming (`seqdist.R:451`); confirmed numerically against live R (§3). It is
  **unguarded by any R-golden** and reuses TraMineR's exact method name, so a
  user porting R code gets silently different numbers. `docs/API_PARITY.md:99`
  compounds this by wrongly asserting DHD is "not present in TraMineR."
  *Recommend:* add an R-golden for DHD (it will fail until reconciled) **or**
  document the divergence explicitly and stop implying parity.
- **M2 — git release tags absent/orphaned.** Only `v0.1.1` exists, on an
  abandoned numbering, pointing to an unrelated commit; 0.0.x releases are
  untagged. See §7.

### LOW
- **L1 — `tests/test_seqtree_smoke.py` never runs.** Guard
  `@pytest.mark.skipif("sklearn" not in globals(), …)` checks the module's
  globals, which never contains `sklearn` (it's never imported there) → the
  smoke test is **always skipped**, even though sklearn is installed. Dead test.
- **L2 — `docs/INDEX.md` minor staleness:** version `0.0.2`, "312 passed",
  "Updated 2026-06-11", and an open-item description that predates the merge.
- **L3 — `ROADMAP.md` frozen** at P0/P1/P2; omits the distance-metrics-parity and
  OM-normalization milestones.
- **L4 — nanobind backend won't compile** with the bundled nanobind 2.9.2 + g++
  14.3 (upstream header bug). pybind11 builds fine; nanobind path is optional.
- **L5 — OMslen_gmean golden coverage asymmetry:** the om-norm test pins 3 of 5
  norms for `OMslen_gmean` (`gmean`/`auto` produced by R but not pinned). Both
  reproduce in live R; purely a coverage gap.
- **L6 — `cpp/bindings/pybind11_module.cpp:1`** has `#pragma once` in a main file
  (harmless compiler warning).

### OBSERVATION (not a defect)
- The om-normalization work was **already merged before this audit ran**; the
  merge message asserts a prior audit. Every asserted check was independently
  re-confirmed here (§1). The "circular fixtures" critical-risk was checked and
  **does not exist** — goldens are literal live-R values.

---

## 9. Resolved open questions (`OPEN_QUESTIONS_FOR_AUDIT.md`)

| # | Item | Status |
|---|---|---|
| 1 | CHI2 default `step=1` needs a live-R golden | ✅ **Answered/Resolved** — `CHI2_default` (no `step` arg) reproduces bit-for-bit from live R; `tx.seqdist` matches. step=1 == TraMineR default. |
| 2 | SVRspell/NMSMST (+NMS) need R goldens (brute-force only) | ✅ **Answered/Resolved** — live-R goldens for NMS/NMSMST/SVRspell reproduce bit-for-bit; `tx.seqdist` matches. No longer brute-force-only. |
| 3 | OM-family normalization contradiction / scoping | ✅ **Resolved** — `norm` implemented for TWED/OMslen/OMspell/OMstran to live-R parity (23 goldens), refused for NMS/NMSMST/SVRspell (live R confirms); `DISTANCES_NORMALIZATION` extended. Independently re-verified. |
| 4 | Feature branch unmerged; version bump deferred | ✅ **Answered** — both `feat/distance-metrics-parity` and `feat/om-normalization` are **merged** into `master`. Version **was** bumped to 0.0.3 but **incompletely** → Finding **H1**. |
| S-a | WIP CI job's bare-mypy step is red | ✅ **Resolved** — CI now runs `mypy --follow-imports=skip` (green). ⚠️ But the same CI job still fails via `check_docs_consistency.py` (H1). |
| S-b | `check_stub_signatures.py` needs `traminex` importable | ✅ **Resolved** — passes (36 stubs); `.venv` traminex imports fine. |
| S-c | Observed local state "312 passed, 2 skipped" | ✅ **Updated** — now **432 passed, 2 skipped**. |

`OPEN_QUESTIONS_FOR_AUDIT.md` is fully resolved and is archived to
`archive/OPEN_QUESTIONS_FOR_AUDIT.md` as part of this audit.

---

## 10. Knowledge graph

The committed graph (`graphify-out/`, dated **2026-06-11**) **predates** the
OM-normalization module: `graph.json` has **0** references to
`traminex/distances/_normalization.py` or its symbols (`resolve_norm_code`,
`normalize_value`, `_maxscost`). Per the audit instruction the graph is
refreshed in this audit's doc-reconciliation commit. The new module adds a
normalization layer (`resolve_norm_code` → integer `normalizeDistance` codes;
`normalize_value` implementing maxlength/gmean/maxdist/YujianBo; `_maxscost`
Yujian–Bo cap) wired into the OM-family dispatch and `_validation` (norm-refusal
set NMS/NMSMST/SVRspell). See the refreshed `graphify-out/GRAPH_REPORT.md`.

---

## 11. Verification before completion

- Full suite re-run from clean tree: **432 passed, 2 skipped**.
- All live-R reproductions: **0.000e+00** max abs diff (23 norm goldens + 12
  main-distance goldens).
- Pre-commit, mypy, stub, generated-artifact, verify_audit gates: **green**;
  `check_docs_consistency`: **red on version only** (H1).
- Working tree confirmed clean of edits to `traminex/`, `cpp/`, `tests/`,
  `SOURCES/`. Scratch `_audit_*` files used during the audit were removed.

**Bottom line:** the distance-metrics and OM-normalization work is
**R-faithful and soundly merged**, with parity independently re-confirmed
against live TraMineR 2.2-12 — **except DHD** (M1), which diverges and is
unguarded. The one release-blocking issue is the **incomplete 0.0.3 version
bump (H1)**, which currently makes CI red on `master`. Both are reported, not
fixed, per mandate.
