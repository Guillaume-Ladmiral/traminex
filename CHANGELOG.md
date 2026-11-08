# Changelog

All notable changes to this project will be documented in this file.

The format is based on Keep a Changelog, and this project adheres to Semantic Versioning.

## [Unreleased]

## [0.2.2] - 2026-11-08

A patch release of 0.2.1, by the maintainer's decision of 2026-11-08. The package code and the
public API (226 names) are those of 0.2.1.

### Changed

- The published files no longer name this machine's folders, so that the repository can be made public
  (maintainer's decision of 2026-11-08). Home-folder paths are written `~/` in commands, `$HOME` in the
  helper-script text of `docs/DEV_BUILD.md` and `$(HOME)` in the `Makefile`; the family's folders are named
  inside their parent folder `PROJECT_SEQANALYSIS_and_SPATSCAN`; a Windows path in
  `archive/OPEN_QUESTIONS_FOR_AUDIT.md` is written from `%USERPROFILE%`. On this machine every default
  resolves to the same location as before. The test that compares the commit hook's black with the
  Makefile's `PYTHON` now expands `$(HOME)` as make does; without that it skipped. Sixteen files, no
  package code.

## [0.2.1] - 2026-11-08

A patch release of 0.2.0, by the maintainer's decision of 2026-11-08. The public API keeps its 226 names.

The findings of the supervision repository's study of 2026-10-25, fixed at the maintainer's request of
2026-11-07 ("fix all known issues"). The public API keeps its 226 names. R's values and behaviour for
these cases come from `scripts/r_goldens/make_goldens_study_findings.R` and are pinned by
`tests/test_r_goldens_study_findings.py` (209 generators, 819 fixtures, 205 writers).

- **`seqdef(states=...)` no longer merges distinct states** (finding F-01). A new name that is also an old
  one was renamed again, so states `1 2 3` renamed `2 3 4` gave `4-4-4` for every row, and swapping two
  names made both rows the same. Every state is now renamed once. R renames one state at a time and
  merges them too, against `seqdef.Rd`, which defines `states` as one label per alphabet state: a
  correction of an R defect, recorded as a `DIVERGENCE` in `docs/API_PARITY_DECISIONS.csv`. Names
  disjoint from the alphabet give the same result as before, and as R.
- **`seqdist(method="NMS")` no longer raises `IndexError`** when every pair of a slice has an empty
  sequence on one side (finding F-04). It returns R's matrix.
- **`dissrep` and `seqrep` refuse weights without a positive total with a `ValueError`** that names them,
  where Python's float division raised `ZeroDivisionError` (finding F-05). R stops too, with "missing
  value where TRUE/FALSE needed".
- **`TSE_to_STS` in event mode is linear in the number of rows** (finding F-02): each person's rows were
  found by scanning every row. 4,000 persons with 4 events each: 24.90 s before, 0.57 s after, the same
  output; R takes 1.70 s.
- **NMSMST and SVRspell run compiled when Numba is dispatched** (finding F-06). On sequences with many
  spells their Python kernels were 40 to 140 times slower than R; the compiled copies of the same loops
  give the same bytes, and on the four workloads `docs/DISTANCES.md` gives took 0.03 to 0.42 s against R's
  0.14 to 0.31 s, after a one-time compilation of about 1 s.
  `docs/DISTANCES.md` gave them as "close to TraMineR" from a low-spell workload; it now gives both.
- **`dissassoc`'s permutation test and `disstree`'s weighted inertia run compiled when Numba is dispatched**
  (TX-AUD-20261010-083). The weighted pair sum in R's order took most of the time; its compiled copy gives the
  same bytes. `dissassoc` on `mvad`'s 712 weighted cases with `R=200` and `weight_permutation="diss"`: 0.88 s
  before, 0.115 s after, the same output; R takes about 0.08 s.
- **The ICU line shows under `make test-python` and `make test-pybind11`** (finding F-03): their `-q`
  hid the header that names the loaded ICU, so the line is also written at the end of a quiet run.
- **With one sequence, a `sortv` string that names no method is refused, as with several**, and a
  one-element list is the sort vector, which gives R's order. R takes any such string as a sort vector of
  length 1 when there is one sequence; recorded as a `DIVERGENCE` of `plot.stslist`.
- Documents: `seqmeant`'s docstring says that the result has a row, zero when uncounted, for each marker
  the object stores (observation O-2); `docs/CONSUMED_SURFACE.md` lists three changes of 2026-09-22 to
  2026-10-17 that reach the other libraries; `ECOSYSTEM.md` names the family folders where they are now;
  `ROADMAP.md` points at git history for the six plans removed on 2026-10-26; `MIGRATION_PLAN.md` no longer
  calls the other projects' copies identical, which they are not; `docs/DEV_BUILD.md` names the
  three release steps no workflow runs (TX-AUD-20261010-091); `scripts/r_goldens/README.md` asserts
  R 4.6.2. The CI workflow stops when its R is not the oracle's (TX-AUD-20261011-079, option (b)).
  `docs/FINDINGS_REGISTER.csv`: the 72 `SUSPECTED` rows each re-registered under a later id carry that
  id's outcome (70 fixed, 2 escalated).

## [0.2.0] - 2026-11-07

The first release since 0.0.8 (2026-06-13), numbered below 1.0 by the maintainer's decision of
2026-11-07; the 1.0 procedure (R6) stays postponed. The public API is unchanged: the same 226 names.
The entries below are the changes made since 0.0.8. An older tag in this history, on a 2025 commit, has
no entry in this changelog (TX-AUD-20261010-089); this release is numbered above it.

- **The oracle's R is 4.6.2** (2026-11-02, owner decision; it was 4.6.1). The machine's R moved through
  `apt upgrade`, and `scripts/check_r_golden_provenance.py` asserts the exact version, so two tests of
  `tests/test_r_golden_provenance.py` failed. Under R 4.6.2 the full two-pass proof, with both packages
  compiled from `SOURCES/`, left all 816 fixtures byte-identical. So no expected value changes. Only
  the header line of 11 fixtures, which names R's version, is rewritten by its generator. Statements of what was
  measured under R 4.6.1, in tests and in comments, stay as records.
- **Tests against R for what the other libraries use and no test here ran** (2026-11-02, owner decision; no
  change to `traminex/`). The weighted branches of `disstree`'s split permutation test, which
  `disstree(D, predictors, weights=w)` enters by default: R's permutations are replayed from
  `tests/r_goldens/disstree_weighted_gate.txt`, and every p-value and the whole tree must equal R's.
  `squared`, `first` and the rounding of weights, in `disstree` and `seqtree`. `resolve_sort_order`'s refusals,
  `_exact_near_ties`' early return, and direct tests of `apply_time_axis`, which were only run through
  other functions (`tests/test_r_goldens_consumed_items.py`). The golden inventory counts follow (208 generators,
  818 fixtures, 204 writers).

- Prevent the chronogram example's main title from overlapping panel headings
  and keep its time labels and shared colorbar visible.

### Added

- A `SeqElist` (the result of `seqecreate`) can be indexed with a list or
  array of positions or with a boolean mask, as R's subset takes any index
  vector; the result is a `SeqElist` of the selected sequences with their
  weights. Positions are 0-based and follow NumPy's rules: a negative position
  counts from the end, and a mask of the wrong length raises `IndexError`.
  It raised `TypeError` (TX-AUD-20261011-058).
- `make sdist-manifest` writes the sdist list (`sdist.include` in `pyproject.toml`) from
  `git ls-files` and the source policy in `scripts/sync_sdist_manifest.py`; `--check`
  names each missing or extra entry. The list stays exact, with no globs, for
  scikit-build-core 0.11.0, but a new tracked file no longer needs a hand edit, and the
  packaging test takes its policy from the script. It refuses a tracked file missing from
  the working tree and a path no exact entry can hold: a glob, escape, quote or control
  character, or a trailing space, which 0.11.0's gitignore matching would drop
  (TX-AUD-20261019-078).

### Changed

- **`seqrecode`'s `otherwise` leaves the missing and void markers alone
  (2026-10-23, owner decision; TX-AUD-20261019-033).** It recodes the real states
  `recodes` does not list; missing and void cells stay markers unless `recodes`
  names them, so no sequence changes length. R, and traminex before, turned them
  into the `otherwise` state too: `a-c-%-%` became `A-O-O-O`, its void padding
  counted as time spent in O. An `otherwise` naming a marker itself (as
  `seqipos` uses it), not a real state spelled like an unused marker, still
  sends every unlisted cell to that marker. A recode named like a marker in use
  combines with `otherwise` as it does without it.
- **The default alphabet is R's order (2026-10-19, owner decision;
  TX-AUD-20261019-071, -038).** Without `alphabet=`, `seqdef` (STS, SPS and
  SPELL input), `seqformat` (to `"STS"`), `seqdef_from_long`, `TSE_to_STS`
  (state and event mode), `events_to_state_sequences` and `seqMD`'s combined
  alphabet now order the states as TraMineR's `seqdef` does through `seqstatl`:
  by numeric value when every value (a `missing=` code included) is an R
  number, R's collation otherwise. They used to keep the order of first
  appearance (`events_to_state_sequences` used code points; `seqMD` sorted
  numeric-looking labels such as `10.1` as text). Consequences: state codes,
  per-state rows and columns (`seqstatd`, `seqtrate`, costs), colours by
  position, the tie rules of `seqmodst` and `seqgranularity(method="mostfreq")`
  and the state-rank defaults of `seqibad`, `seqidegrad`, `seqprecarity` and
  `seqinsecurity` now equal R's without an explicit alphabet, and `seqpolyads`
  no longer refuses default role objects whose states first appear in
  different orders. `states=`, `labels=` and `cpal=` keep pairing with R's
  order. A pandas categorical SPELL status keeps its declared categories in
  `seqformat`'s own output, as R's factor status does, and an STSList input
  keeps its alphabet's order there (a traminex choice: R returns a plain
  matrix). **Migration:** pass `alphabet=` to `seqdef`. `seqformat`,
  `seqdef_from_long`, `TSE_to_STS`, `HSPELL_to_STS` and `seqMD` take no
  `alphabet=`: pass their result to `seqdef(x, alphabet=[...])`, which keeps
  every cell and marker. `events_to_state_sequences` takes `states=`.
- Removed helpers that no package code used (TX-AUD-20261019-063,
  TX-AUD-20261010-072): the eight `traminex.stats_complexity.*_proxy`
  functions (`seqivardur_proxy`, `seqlogp_proxy`, `seqmaintokens_proxy`,
  `seqprecarity_proxy`, `seqprecstart_proxy`, `seqibad_proxy`,
  `seqintegr_proxy`, `seqipos_proxy`; use the public functions of the same
  name without the suffix), `traminex.visualization.ticks.apply_time_axis_one_based`
  and several private ones. None was in the public API.
- **The parity oracle is now TraMineR 2.2-14 (2026-09-20).** `SOURCES/`
  vendors the 2.2-14 CRAN source release (published 2026-09-15; CRAN MD5
  `609ae6049321405695e5ed6aa7a7c6e7` matches the `PACKAGES` index, internal
  `MD5` manifest 347/347); TraMineRextras stays at 0.6.9. `src/` is
  byte-identical between the two releases, and 23 of 23 reference quantities
  recomputed on the same corpus are byte-identical, so **no distance, cost,
  statistic or tree result changes**. The 2.2-14 `NAMESPACE` drops
  `Math`/`Ops`/`Summary` for `eseq`/`seqelist` and adds `summary.seqelist`,
  `summary.eseq` and `print.eseqsum`, so the semantic ledger moves from 250 to
  247 decisions and the R denominator from 243 to 240. All 108 golden
  generators, the provenance checker's version guard, the sdist allow-list and
  the prose now name 2.2-14; dated historical records keep their original
  version. The complete 2.2-13 → 2.2-14 delta is in the report of the oracle-2214
  substitution audit of 2026-09-20, kept in the family supervision repository's
  `docs/history/run-reports/` (its evidence was removed on 2026-10-26).
- **`om_trate_pairwise(indel=None)` now resolves its automatic indel as R
  does** — `max(sm)/2` (`seqdist.R:359-365`), the same rule
  `seqdist(subm="TRATE", indel="auto")` uses. It previously used
  `0.5 * mean(off-diagonal)`, so the package's two TRATE entry points returned
  different numbers for one corpus (`[0, 0.3125, 0.625]` against
  `[0, 1.25, 4.0]`). `indel_from_subm` remains re-exported from
  `traminex.distances` for the historical rule.
- **traminex deliberately does not reproduce 2.2-14's `INDELS` regression.**
  2.2-14's new internal `checkcost()` call aborts on any `sm` containing
  `Inf`, which `INDELS` produces when two or more states are absent from the
  data (for time-varying costs, from one time slice), so
  `seqsubm(method="INDELS")`, `seqdist(sm="INDELS")` and the time-varying
  `INDELS` costs fail upstream on such data (corrected 2026-10-13, audit
  TX-AUD-20261010-025: this entry first said they always fail; with at most
  one absent state 2.2-14 computes them). traminex keeps computing them and
  still refuses an explicit non-finite matrix with a named error. Recorded in
  `docs/API_PARITY.md`.
- **`OMConfig`'s default indel is now `"auto"`, as in `seqdist` (audit
  TX-AUD-20261011-043).** It was `1.0`. `"auto"` resolves as `seqdist`
  resolves it: half the largest substitution cost for a number or a matrix
  (`seqdist.R:359-365`), and the cost method's own indel for a method string
  such as `"TRATE"`. `OMConfig()` gives the same distances as before (its
  cost 2.0 gives an indel of 1.0), but other costs do not: `OMConfig(subm=4.0)`
  now uses an indel of 2.0 and a matrix whose largest cost is 5 an indel of
  2.5, so distances change (an un-normalized one can only grow). Pass
  `indel=1.0` to keep the previous values.

### Fixed

- `seqformat`'s wide-to-long route no longer passes `inplace=True` to
  `DataFrame.rename`, which pandas 3.1 deprecates (PDEP-8): with pandas 3.1 every such call
  raised a `Pandas4Warning` from traminex (TX-REP-20261025-004).
- The install-profile smoke test checks the distances it computes, [[0, 2], [2, 0]] for
  A-B-B and A-A-B, not only their shape (TX-AUD-20261019-058). `make test-pybind11` builds
  and tests with the Makefile's `PYTHON`; the runner, started directly, falls back to the
  `python3` on `PATH` instead of one user's interpreter (TX-AUD-20261011-078).
- The documents no longer claim a CI run, build or type check: the repository has no
  git remote, so the workflows under `.github/workflows/` have never run.
  `docs/DEV_BUILD.md`, `ROADMAP.md` and five earlier entries of this section now describe
  the configured workflow, and a test keeps the active documents from claiming a CI run
  (TX-AUD-20261019-067).
- The ICU loader reads the version through a library symlink: an absolute path, as
  `ctypes.util.find_library` gives on macOS (often the unversioned `libicui18n.dylib`),
  is resolved to its file before the version is read. The `icuinN.dll` form is removed:
  on Windows `find_library("icui18n")` can give only an `icui18n[.dll]` name, so that
  form was never reached; such a name is refused with the supported form named. Only
  the Linux form has been checked on a real host (TX-REP-20261019-005).
- `docs/DEV_BUILD.md` says that the editable install puts the repository root on
  `sys.path` (scikit-build-core's `.pth` file adds each package's parent directory, the
  root for this flat layout), so `benchmarks`, `scripts` and `tests` are importable from
  anywhere and can mix with other installed names, and that a `src/` layout would end it.
  The layout is unchanged (TX-AUD-20261011-077).
- The R-golden provenance proof (`--execute`) records the version of every R package
  the oracle reaches: the recursive Depends/Imports closure, computed by R at run time,
  of the packages the generators load or call with `::` (40 packages under the
  maintained oracle, `cluster` among them). It recorded only R, TraMineR,
  TraMineRextras and `survival`, although generators also load WeightedCluster, vegan,
  colorspace and RColorBrewer. A package of the closure that R cannot find stops the
  proof (TX-AUD-20261011-063).
- `traminex.__version__` is the version of the files actually imported. It returned the
  installed distribution's version whenever one was installed, so a copied tree or a
  second worktree imported ahead of it reported the installed version, not its own. The
  installed version now counts only when that distribution owns the imported
  `__init__.py` (a wheel install puts it beside its metadata, with or without a RECORD;
  an editable install names its tree in `direct_url.json`); any other tree reads the
  `pyproject.toml` beside it
  (TX-AUD-20261011-032).
- A benchmark report names what it timed: `config.kernel_options` holds the options the
  kernel ran with (`{"variant": "mean"}` for DHD, the frequency variant, not R's default
  DHD; `{"metric": "chi2"}` for DUR; indel 1 and substitution cost 2 for OM), and
  `metadata.package` holds the imported traminex's file, `__version__` and a SHA-256 of
  its Python sources, where it held only the installed distribution's version. The
  `--method` help says the same (TX-AUD-20261011-007).
- `DATA/TraMineR/ex1.rda`, `ex2.rda` and `famform.rda` are now the vendored files of
  `SOURCES/TraMineR/data/`: their bytes differed, while every R object they hold was
  `identical()` to the vendored one, so no data value changes. A test requires each of
  the 10 `DATA/` files to equal its vendored copy byte for byte (TX-AUD-20261010-098).
- `docs/INDEX.md` says that no Make target, workflow or script regenerates or checks
  `graphify-out/`, the knowledge-graph snapshot of 2026-07-03, which is kept as it is
  (TX-AUD-20261010-092).
- `docs/INDEX.md` names every tracked root file, top-level directory and entry of
  `docs/` and `scripts/` (a test reads `git ls-files`), lists the task prompts and
  verification records at the root as active or historical, describes `_spells.py`
  (spell helpers) and `omloc.py` (the separate `omloc()`) as they are, and gives
  `TRAMINEX_BINDING=auto|python|pybind11` with auto the default. It and README.md label
  the 2026-06-12 report as a historical audit and name the latest one
  (TX-AUD-20261019-064, TX-AUD-20261019-076, TX-AUD-20261010-087).
- `ECOSYSTEM.md`, `docs/PLOTTING.md` and `docs/API_PARITY.md` describe traminex's plotting
  as matplotlib renderers that re-express the TraMineR plots (the ledger records none as
  EXACT), where they called it a lightweight reference layer (TX-AUD-20261019-074).
- `docs/DEV_BUILD.md` no longer says the `datasci` environment has traminex installed (it
  has neither traminex nor `tokenize_rt`), and gives the private virtual environment,
  layered on `datasci` with `--system-site-packages`, that the gates run in and that leaves
  `datasci` unchanged; `docs/DOCS_POLICY.md`, `docs/INDEX.md` and `CONTRIBUTING.md` stop
  sending the checks to the `datasci` interpreter (TX-AUD-20261010-084).
- `make audit` also runs `scripts/generate_semantic_parity_ledger.py`, so it refreshes
  `reports/audit/traminex_semantic_parity.json` and `docs/API_PARITY_LEDGER.md`, as
  `docs/DOCS_POLICY.md` says it does (TX-AUD-20261010-086).
- `gower_matrix` takes any weight total, as R's `gower_matrix` does: a zero total gives R's
  NaN matrix and a negative one R's finite values (the shares `w / sum(w)` still sum to 1),
  without a warning; it refused both with "weights must have a positive sum". `dissmfacw`
  still refuses such weights before it centres (TX-AUD-20261011-005).
- The parity ledger records four departures from R that its rows called EXACT, with a
  test of each and R's values beside them: `is_stslist` is True for a `seqdef` object
  with `id=` and `weights=`, which R's `is.stslist` rejects; `seqLLCS` and `seqLLCP`
  compare the non-void cells, where R compares the first `seqlength()` cells and so
  differs for an interior void; `seqecontain` accepts a single event sequence, which R
  refuses (kept by the owner's decision). No behaviour changes (TX-AUD-20261011-005).
- `DUPLICATES.md` names the sibling modules that call `weighted_pairwise_sum` today, and
  `MIGRATION_PLAN.md` and the 2026-09-21 consumer notice open with a dated status note
  listing their statements that no longer hold (HAM's default, `seqtransn`'s `norm=`,
  `dissassoc`'s permutations; a superseded re-pin digest) (TX-AUD-20261011-030).
- The Quarto documentation site publishes the top-level documents only and opens on the
  index: it copied the 127 historical files of `docs/superpowers/` and redirected its
  landing page to `API_CONSISTENCY.html` (TX-AUD-20261011-054).
- The R golden generators name R 4.6.1, the oracle, and point to
  `scripts/r_goldens/README.md` for the run command, where nine named a conda R 4.4.3 and
  31 a private note; only comments changed (TX-AUD-20261011-064).
- `make install-profiles` writes its report to `build/install-profiles.json`, which git
  ignores, instead of `install-profiles.json` at the repository root. `make test-pybind11`
  removes, before it builds, the `traminex/_pb_kernels*.so` link a killed run left behind,
  which every later run refused. Such a link stays live while its build directory exists
  (until `/tmp` is cleared): the runner now records its process ID in that directory, and
  a link into one of its build directories is removed, with the directory, when the
  recorded process is gone and no process has the module loaded, or when the target is
  gone. A link of a run in progress, a regular file or a link pointing elsewhere is still
  refused (TX-AUD-20261010-088).
- `scripts/calibrate_parity.py` records only parameters that reproduce the R reference,
  every value within an absolute 1e-6 (TRATE) or 2e-3 (DHD), the tolerances of
  `tests/references/README.md` without a relative part, and otherwise records none and
  says so. It recorded the closest candidate
  whatever its error: against R's default DHD (5.5 for two sequences) it wrote the
  frequency variant `sum` with `smooth_window` 5 (0.55), which a `dhd_profile="calibrated"`
  call would then use. No frequency variant reproduces R's DHD; the script now says that
  R's DHD is the public `seqdist(method="DHD")` and writes no `"dhd"` section. When
  nothing reproduces it writes an empty calibration (`{}`), which a calibrated profile
  refuses, so an earlier file's sections do not survive (TX-AUD-20261019-087).
- Two warnings have R's text: TWED with one indel cost per state warns " [!!] indel vector
  not supported by the chosen method, max(indel) used instead!" (`seqdist.R:301-305`), and
  `seqrf` / `dissrf` with no degrees of freedom left warn " [!!] Negative  or zero df! k
  = 3  wsum =  3", as `dissrf.R:262` builds it (the weight sum as R's `format()` writes
  it). traminex had its own wording for both (TX-REP-20261023-008).
- The documentation check (`make docs-check`) reports an unknown keyword in a documented
  `seqdist` call with a literal method, and in a `seqdef` call with the STS format: both
  take their options through `**kwargs`, which accepts any name, so such a call passed
  the check although it fails when run. A keyword is accepted when `seqdist` accepts it at
  run time: the keywords it reads itself (read from its source), for one of R's methods
  every R argument (one the method does not use is ignored, as in R) and `check_costs`,
  `check_max_size` and `opt_args`, and the per-method list of its run-time check. A method
  or format passed through `**options` or `*args` is not assumed. `HSPELL_to_STS`,
  `TSE_to_STS` and `seqdef_from_long` are checked against the names they accept; the 24
  other public functions with `**kwargs` pass their keywords on and are named as
  unchecked (TX-AUD-20261019-054).
- The documentation check reads a `from traminex import (...)` statement spread over
  several lines as the names it imports; it read the parenthesis as an API name, so a
  valid document failed, and it took an alias (`name as alias`) for a name
  (TX-AUD-20261019-055).
- The test that runs the documented examples compiles each one without the test
  module's `from __future__ import annotations`, so an annotation naming an undefined
  name fails there as it fails in a reader's script (TX-AUD-20261010-063).
- The R-golden provenance check (`check_r_golden_provenance.py --execute`) asserts R
  4.6.1, the documented oracle, with the package versions; it printed R's version
  but did not check it. The suite's check that every generator pins TraMineR now
  requires code comparing `packageVersion("TraMineR")` with `"2.2.14"`; the bare text
  passed in a comment or a printed version (TX-AUD-20261010-090).
- `make release-check` passes its `PYTHON` to the check that the pre-commit formatter
  revisions are not older than the maintained environment's tools, so the check also
  runs where that path does not exist, as in the configured CI workflow
  (`PYTHON=python`; the workflow has never run, the repository having no git
  remote); it read only the Makefile's default path, which exists on
  one machine, and skipped elsewhere (TX-AUD-20261011-076).
- The stub check (`make stub-check`) checks a renamed re-export (`from .x import a as b`)
  in a stub without a literal `__all__` when `b` is public at run time: such a stub
  describes the runtime `b` as `a`, and only the `b as b` form was checked there
  (TX-REP-20261019-004).
- The MDS orientation test checks that the oriented axis is R's raw axis up to its
  sign; it asserted R's raw sign, which follows the numerical library, so it failed on
  an oracle whose LAPACK returns the other sign (TX-REP-20261019-006).
- Two tests that compared with R values copied by hand now read them from the R
  goldens: the 13 OMstran option values (already recorded by
  `make_goldens_distances_void.R`) and `seqdef` on a numeric factor (a new case of
  `make_goldens_default_alphabet.R`). A test requires every test module whose test names
  claim agreement with R to read a golden or run R, apart from two named modules whose
  claims rest on R's source (TX-AUD-20261019-088).
- Every docstring example runs as written: the examples of the distance configuration
  classes, `save_plot`, `seqient` and `seqici` used names they never defined or showed
  `True` where NumPy 2 prints `np.True_` (14 of 114 failed). A test runs doctest over every
  traminex module (TX-AUD-20261010-097).
- Every comparison in the R-golden tests states its tolerance: 108 calls relied on a
  library default (`assert_allclose`'s rtol 1e-7, `np.allclose`'s 1e-5, `pytest.approx`'s
  1e-6, pandas' 1e-5) and now apply the documented contract, absolute 1e-10
  (`R_GOLDEN_ATOL`), which all of them meet; the three looser tolerances of the discrete
  sequence-rule tests (rtol up to 2e-6) are tightened to it as well. A check that two sets
  of values differ requires a difference above 1e-6 (`R_DIFFERENCE_ATOL`; the differences
  asserted are 0.2 or more), since with the equality tolerance it would pass on 1e-10. A
  test fails on a default tolerance in a golden test, and on a "values differ" check with
  the equality tolerance (TX-AUD-20261019-075, TX-AUD-20261011-065).
- The pybind11 gate (`scripts/run_pybind11_tests.sh`) fails when a test skips because
  the compiled kernels are unbuilt or not in use, as `docs/DEV_BUILD.md` requires; pytest
  exited 0 with every compiled test skipped. Other skips (Graphviz, Rscript) do not fail
  it. The optional `TRAMINEX_PARITY=1` test is described as what it is: it checks that
  `seqdef` builds TraMineR's `actcal` data from the local `DATA/` folder and compares
  nothing with R's output (TX-AUD-20261019-077, TX-AUD-20261011-075).
- The test suite treats a warning (`UserWarning`, including the cost diagnostics, and
  `DeprecationWarning`) as an error. It emitted 269 warnings, R's diagnostics met by
  tests that feed unusual costs or data on purpose, among which a new, unexpected one
  would have gone unnoticed; each of the 52 tests concerned now names the diagnostics it
  expects. Blocks inside tests that switch warnings off still hide them: 81 calls ignore
  every warning and 7 a whole category (TX-AUD-20261019-068). A `DeprecationWarning`
  raised inside matplotlib's own modules is not an error: with matplotlib 3.6.0, the
  declared floor, pyparsing 3.3's deprecations of the names it calls stopped the suite
  from loading 72 test modules (TX-REP-20261025-001).
- The index plot (`seqIplot`, through `select_index_rows`) takes a unique prefix of
  `sample_meth` (`"sys"`, `"r"`), as R's `plot.stslist` reads it with `match.arg`; it took
  exact names only (TX-REP-20261023-007).
- `seqemlt` stops with R's message, "subscript out of bounds", when no state is observed
  (no sequence, or only missing cells); it raised Python's `ZeroDivisionError`
  (TX-REP-20261023-005).
- `seqdiff` stops with R's message, "invalid 'ncol' value (too large or NA)", when a
  comparison window has no complete case (no sequences, only void or only missing cells),
  and `dissassoc` (and so `seqdiff`) gives R's NaN results for an all-zero weight vector;
  both raised Python's `ZeroDivisionError` (TX-AUD-20261010-096).
- `seqalign` refuses indices that are not two integer positions (a boolean, a fraction,
  one or three positions) with a `TypeError` naming them; they raised raw `TypeError`,
  `IndexError` or unpacking errors (TX-AUD-20261011-074).
- `seqimplic` refuses a group without one entry per sequence; a shorter group raised a
  raw `IndexError` and a longer one was computed silently. R stops on a longer group and
  returns no finite value for a shorter one (TX-AUD-20261011-073).
- `seqsha` stops, as R does, on a missing time or event (R's message) and on a time
  below 1, fractional or infinite; it computed with a time of 0 (writing the first person's event at the end)
  and raised numpy's errors otherwise. `event` must be logical (booleans, or 0 and 1):
  text events were read as `True` ("no" counted as an event) and are now refused
  (TX-AUD-20261011-069).
- `seqtransn_count` counts and normalizes transitions as R's `seqtransn` does (it now
  equals `seqtransn(x, norm=normalized)`): missing cells are left out of the spells
  and the count is divided by each sequence's own length minus 1. It counted changes
  into and out of the missing marker and divided by the padded width, so `B-B-C`
  padded to 6 positions gave 0.2 instead of 0.5; `turbulence_proxy` averages the new
  values. docs/SEQTREE.md no longer says that `seqtree_explore`'s summary feature is
  `seqtransn_count` (TX-AUD-20261011-068).
- `seqpcfilter` takes a unique prefix of the method (`"cum"` for `"cumfreq"`), as R's
  `match.arg` does; it took exact names only. `level` is checked for the matched method,
  so `seqpcfilter("cum", None)` is refused where R lets it through (TX-AUD-20261011-002).
- `seqgen` reads `p` as probability weights and normalizes them, as R's `sample()`
  does (`p=[5, 3, 2]` draws as `[0.5, 0.3, 0.2]`); it refused weights that did not sum
  to 1. A negative, missing or infinite weight, or no positive weight, stops with
  `sample()`'s message (TX-AUD-20261011-072).
- `seqgen` keeps the alphabet's type, as R does: a numeric alphabet gives numbers (it gave
  text) and a boolean one booleans; a text alphabet still gives an object array of
  states. A one-value alphabet, also given as a scalar (`alphabet=5`), draws that value;
  R draws 1 to 5 there, a quirk of `sample()` (TX-AUD-20261011-001).
- `seqdiff` refuses a group mixing text and numbers with the TypeError of the other
  group functions (it failed while sorting the levels). That refusal, and the plots'
  refusal of such a `sortv`, name the kinds of value found ("text and numbers",
  "numbers and booleans", "text, numbers and booleans"); a group of numbers and
  booleans was said to mix text and numbers. `seqindic_dyn` on an object
  without sequences stops with R's message, "invalid 'ncol' value (too large or NA)",
  where numpy's "zero-size array" error came out (TX-AUD-20261010-096).
- `entropy` of uniform data is exactly 1 (normalized) and `ln K` (raw): a zero frequency
  contributes 0, where `1e-12` was added inside the log (0.9999999999971146). The same
  holds for `seqtree_explore`'s summary entropy and `local_indel_by_pos(method="entropy")`.
  `entropy`'s docstring says that the missing state is one of the `K` states
  (TX-AUD-20261010-095). `seqdist(method="JSD_TRANS")`'s logarithm constant `eps` now
  defaults to 0 (it was `1e-12`), so disjoint transition distributions are exactly `ln 2`
  apart; only positive shares enter the logarithm (batch-15 review M-2).
- `seqdomassoc` takes `assoc` as one string (`assoc="LRT"`), as R does; a string was read
  character by character and refused (TX-AUD-20261010-062).
- `seqstart` stops with R's messages on an empty window: with no sequences and no
  `tmax` ("There are no data in the specified new time frame"), with no sequences
  and no `tmin` ("invalid 'ncol' value (too large or NA)"), and when `tmin` is more
  than one position after `tmax` ("invalid 'ncol' value (< 0)"); it raised numpy's
  "zero-size array" and "negative dimensions" errors. The ledger row records that
  it takes an `STSList` where R's `seqdata` may be any data frame or matrix
  (TX-AUD-20261011-062).
- `seqdef_from_long` and `TSE_to_STS`'s state mode read numeric states as R prints
  them, as `seqdef` does: `1.0` and `2.0` are the states `1` and `2`; they were
  `1.0` and `2.0` (TX-AUD-20261010-073).
- `seqfpos`'s per-state table (`state=None`) no longer has a void column with
  `include_missing=True`: void is not a state (it gave the position where a
  sequence ends). `docs/API_CONSISTENCY.md` states the three per-state width
  rules, where it said every per-state result is marker-padded
  (TX-AUD-20261019-070). `seqfposend`'s table follows the same rule; it too had a
  void column (TX-REP-20261023-004).
- Both marker representations of an object give the same results in `alphabet`,
  `seqstatf` and `seqtrate`: an object declaring the `nr` marker's text in its
  alphabet beside an external `nr_code` got an extra `*` state (and an `IndexError`
  from `seqstatf(include_missing=True)`) where its canonical twin did not
  (TX-AUD-20261010-094, TX-AUD-20261019-042). An `STSList` refuses `nr` text equal
  to `void` text, as `seqdef` does (TX-AUD-20261019-042). The same holds for
  `seqetm` (a 4 by 4 matrix with `*` twice) and for `seqformat`'s STS, DSS and SPELL
  outputs from an `STSList` and `STS_to_SPELL` with `right` (the missing marker
  written `<NR>`) (TX-REP-20261023-003).
- `seqMD_gridsearch_weights` refuses a `weight_step` above 1 (the grid runs from 0 to
  1, so a larger step left no vector and ended in `RuntimeError: No valid weight
  vector evaluated`); it and `seqMD_optimize_weights` refuse group labels with a
  missing value; `SeqMD` refuses `names` that do not give one name per channel.
  The docstring examples run, with their real output, and `docs/INDEX.md` names
  the five Config classes (TX-AUD-20261019-041, TX-AUD-20261010-075).
- `seqdef_from_long` and `TSE_to_STS`'s state mode keep the persons in order of
  first appearance, as R's route for long data (reshape, then `seqdef`) and
  `TSE_to_STS`'s event mode do; they sorted them (`p1, p10, p2` for `p10, p2, p1`)
  (TX-AUD-20261019-039). They refuse a missing id with R's message ("missing values
  in 'row.names' are not allowed"), as R's route stops at `seqdef`; it was kept as
  a person with no id.
- `HSPELL_to_STS`'s per-row bounds `tmin_col` and `tmax_col` (a traminex
  extension) leave the periods outside each person's window unfilled, as R's
  `endObs` does: they become missing (before) or void (after). `fillblanks`
  filled them, so a person's unobserved periods read as the filler state
  (TX-AUD-20261019-037).
- `seqdef_from_long`, `TSE_to_STS` (both modes) and `seqformat`'s STS output from
  a table give the states `seqdef`'s labels and colours, the missing marker
  `missing_color` and void white, as R's route (reshape, then `seqdef`) does; they
  had none, so plots used another palette and coloured the void marker like a
  state (TX-AUD-20261019-035). The long builders take `labels=` and `cpal=` with one
  entry per real state, as `seqdef` and R do; they refused them when a marker was
  needed (TX-REP-20261019-013).
- `seqrplot` (and `seqplot(type="r")`) titles its y axis as R does,
  "N representative(s) (n=N)" ("weighted n=N" for weighted sequences), and
  labels none of the bands; it read "Coverage (na%)" with a percentage per band
  (TX-REP-20261023-001).
- `sequence_index_plot` and `state_distribution_plot` draw missing cells in
  the missing colour, opaque, as R does: `missing_alpha` now defaults to 1
  (it was 0.3 and 0.5, which washed the missing colour out with white or
  transparency). Pass the old value to keep the previous look
  (TX-AUD-20261011-067).
- `seqsurvplot`'s legend names the states by their labels, as `seqsplot`'s
  legend does; it showed the state codes (TX-AUD-20261011-067).
- `seqplot` takes R's `cpal`, `ylab` and `with_legend` for every type (and
  with `group`); `ylab`, `with_legend`, and `cpal` for the index, distribution
  and entropy routes raised `TypeError`. As in R, `cpal` is ignored for
  `type="Ht"`, and `with_legend="auto"` (the default) draws one state legend
  for the figure, except for `Ht` and `pc`: `seqplot` now shows a state legend
  for the index, frequency, modal-state, mean-time, representative and
  relative-frequency plots, which it drew without one. Pass `with_legend=False`
  for the previous look. `seqplot` leaves alone the legends it did not draw (in
  a figure given through `ax`); `seqfplot`, `seqmsplot`, `seqmtplot`,
  `seqrplot` and `seqrfplot` keep drawing no legend when given `main` or
  `group`; `plot_SAMM` draws one state legend for the figure, as R's
  `plot.SAMM` (TX-AUD-20261011-067).
- `seqfplot`, `seqmsplot`, `seqrplot` and `seqrfplot` put each time tick at
  the centre of its cell, as R does; the ticks sat on the cells' left edges
  (TX-AUD-20261011-067).
- Time ticks follow R's grid in every plot that draws them: a tick every
  `xtstep` positions from the first, plus the last with `tick_last`. As R's
  `axis()`, the labels are drawn left to right and one that would come closer
  than the width of an "m" to the last label drawn is left out, on the axis'
  size when it is drawn. Plots showed at most 10 evenly spaced ticks, ignoring
  the `xtstep` grid (`xtstep=6` on 72 positions gave ticks at 1, 9, 17, ... in
  place of R's 1, 7, 13, ...). `apply_time_axis`'s `max_labels` no longer
  limits the labels (TX-AUD-20261011-067).
- `seqfplot` labels its y axis as R does: 0% at the bottom and the cumulative
  percentage of the shown sequences at the top, titled "Cum. % freq. (n=N)"
  ("weighted n=N" for weighted sequences). It labelled each band with its own
  percentage under the title "Cum. % freq." (TX-AUD-20261011-067).
- Every `DeprecationWarning` (deprecated arguments of `seqdist`, `seqdiff`,
  `seqformat`, `seqtab`, `seqtabstocc`, `seqipos`, `seqtrate` and
  `seqtree_leaf_summary`; the deprecated `disscenter_double_center` and
  `seqeordplot`) is attributed to the caller's line, so Python's default
  warning filters show it, as R always shows its warning. Several were
  attributed to a line inside traminex (a decorator, or a traminex function
  that forwards the argument) and so were hidden (TX-AUD-20261011-031).
- `cluster_chronograms` and `cluster_chronograms_with_representatives` refuse
  `labels` that do not hold one label per sequence, and
  `sequence_index_plot_grouped` refuses an `order` that is not a permutation of
  the rows, with a `ValueError` and no figure left open. Shorter labels and
  such orders were accepted silently, and longer labels raised a raw
  `IndexError` (TX-AUD-20261011-029).
- `chronogram` (both drawing routes), `cluster_chronograms` and
  `cluster_chronograms_with_representatives` draw one row per state of
  `seqstatd`: the real states, plus the missing state when
  `include_missing=True`. They drew a row for every internal code, so the
  missing and void markers appeared as states (empty rows by default, and
  the void row always). The representative markers sit on the new rows; a
  state without a row is not marked (TX-AUD-20261019-052).
- `seqsurvplot`'s docstring says that a state without events is drawn as a
  flat line at 1, as the function draws it; it said such a state got no curve
  (TX-AUD-20261019-051).
- `seqedplot` recycles a `cpal` shorter than the events, so every event's curve
  and legend entry has a colour, in every group panel; such a palette raised
  `KeyError` naming the first event without a colour. R draws those events in
  no colour. An empty `cpal` raises `ValueError` (TX-AUD-20261019-050). A
  string `cpal` is one colour, as R's `cpal = "red"`; it was split into
  characters and raised a Matplotlib error (TX-REP-20261023-002).
- `plot_subseqelistchisq` gives every group's panel the same frequency range,
  from 0 (or a smaller frequency) to the largest frequency of any group among
  the rows shown (the first `top`), as R's default `ylim = "uniform"`; each panel was scaled to its own bars, so bars of
  different panels could not be compared by length (TX-AUD-20261019-049).
- A weighted sequence index plot (`sequence_index_plot`, `seqplot(type="i")`)
  draws each row with a height equal to its sequence's weight, as R's
  `barplot(width = weights)` does: the rows touch, from 0 up to the total
  weight, and the row labels and highlights follow them. All rows were drawn
  with the same height. Equal weights and `weighted=False` keep equal rows
  (TX-AUD-20261019-048).
- `docs/PLOTTING.md` lists `seqplot`'s `pc` route (`seqpcplot`), which the
  dispatcher already took (TX-AUD-20261019-086).
- `disstreedisplay` and `seqtreedisplay` find `dot` when `GVPATH` names
  GraphViz's `bin` directory, as well as its installation directory or the
  executable; a `bin` directory was looked up as `bin/bin/dot` and missed.
  When `dot` cannot be found the error names `GVPATH`, the paths tried and
  the `.dot` file left for rendering by hand; it said only "GraphViz 'dot'
  was not found" (TX-AUD-20261010-082).
- `seqtree_fit`'s features (and so `seqtree_cv_fit` and `seqtree_predict`)
  are over the alphabet's states: the missing and void markers have no
  duration or transition column, a marker cell joins no transition (the
  transition count included), and the entropy is over the states alone
  (normalised by the log of their number).
  The features had marker columns and counted marker cells, so a padded
  `B B % %` had the transition `tr_B_%` and an entropy of 0.5. The
  exploratory helpers' documentation now says what they are, an MDS
  embedding or per-sequence features fed to a scikit-learn decision tree,
  instead of calling them an "analog" of `disstree` (TX-AUD-20261010-081).
- `disstree_plot` and `seqtree_plot` (the exploratory tree plots) raise
  `ImportError: matplotlib is required for ...; install the 'plot' extra`
  without matplotlib, as the package's other plotting functions do; they
  let a bare `ModuleNotFoundError` through (TX-AUD-20261011-033).
- Importing `traminex.disstree`, `traminex.seqtree` or `traminex.cluster`
  no longer imports scikit-learn, which set `KMP_DUPLICATE_LIB_OK` and
  `KMP_INIT_AT_FORK` in the process environment. scikit-learn is imported
  when an exploratory helper (`disstree_fit`, `seqtree_fit` and their
  plotting, prediction and DOT helpers) or the scikit-learn clustering
  helper runs; without it they raise the same "scikit-learn is required
  for ..." error (TX-AUD-20261011-034).
- `weighted_pairwise_sum` and `weighted_inertia` check their `indices`:
  zero-based integer positions, in range and not repeated. A float or
  boolean array raises `TypeError` and a negative, out-of-range or repeated
  position `ValueError`; they were truncated, counted from the end or
  counted twice, and `weighted_inertia(D, w)`, with the weights in the
  indices' place, returned 0.0. The signatures are unchanged
  (TX-AUD-20261011-036).
- `dissindic` accepts a zero weight, as TraMineRextras does: that case's
  marginality and gain are NaN (its residual is 0/0) and the other cases
  are computed. It raised `ValueError: weights must be strictly positive`; a
  negative weight is still refused. A group whose weights are all zero takes
  no part: its cases' indicators are NaN and the others' are those of the
  data without it; R agrees unless that group is the first level, where its
  QR adds a direction outside the design (recorded as a correction of R). A
  character `group`, which R's `dissindic` refuses, is computed as R computes
  the same group as a factor (recorded as a correction of R)
  (TX-AUD-20261011-066).
- `seqtree` passes `seqdist_args` to `seqdist` as given: only `None` (the
  default) takes R's default, LCS with `norm="auto"`. An empty dict calls
  `seqdist` with no arguments, which raises `ValueError`, as R's
  `seqdist.args = list()` stops; it built an LCS tree (TX-AUD-20261010-080).
- `dissrf` and `seqrf` accept a text or categorical `sortv`, ordered as R's
  `order(sortv)` orders it: text in R's collation (ICU), a pandas
  categorical by its categories' order, missing values last. They raised
  `ValueError: sortv must contain numeric values` (TX-REP-20261021-001).
- `dissrf` and `seqrf` sort cases on the same default MDS axis whatever
  the scale of the distances, as R does. traminex's test for a repeated
  leading eigenvalue was absolute, so distances around 1e-10 (eigenvalues
  near 1e-20) counted as a repeated root and the cases kept their input
  order: an R2 of 0.167 where R gives 0.623. The test is now relative
  (TX-AUD-20261011-052).
- Medoids break ties between identical cases as R does: each case's inertia
  contribution is summed in the order of R's C kernel
  (`tmrWeightedInertiaContrib`), so identical cases tie exactly and the
  first is the medoid. On the weighted mvad tree one node's medoid was case
  139 where R's is 66 (both at distance 0). `dissrf` and `seqrf` take that
  exact first minimum; their near-tie tolerance (relative and absolute
  1e-12), which chose the wrong medoid when distances are around 1e-14 or
  smaller, is removed. The weighted mean `disscenter` subtracts from the
  contributions is R's `weighted.mean`, with R's long-double sum, so
  `disscenter`'s values equal R's to the last bit. Affects `disscenter`,
  `disstree`, `seqtree`, `dissrf` and `seqrf` (TX-AUD-20261011-023).
- `weighted_inertia`, `weighted_pairwise_sum` and weighted `dissvar` add
  their terms and the weight total in the order of R's C kernel
  (`tmrWeightedInertiaDist`), so they equal R's values to the last bit; they
  summed with NumPy's pairwise sums, which differ in the last bits. The sum
  is vectorised instead of a Python loop over rows: `dissassoc`'s
  permutation test on 712 weighted cases (R = 200) runs about 1.8 times as
  fast; it is still several times slower than R (TX-AUD-20261010-083).
- `dissassoc` computes a group whose weights sum to a negative total, as R
  does; only a total of exactly 0 gives discrepancy 0. Its Levene residuals
  use R's `weighted.mean` (`sum((x * w)[w != 0]) / sum(w)`) without a guard,
  so a zero total gives NaN, as in R. traminex skipped negative totals and
  returned 0 residuals for them, which changed the explained and residual
  sums of squares and the Levene statistic (TX-AUD-20261010-079).
- The repr of a `DissTree` is a short summary, as R's `print.disstree`
  summarises a tree: node, leaf and depth counts, n and the global R2, for
  example `DissTree(n=711.57, nodes=7, leaves=4, depth=3, R2=0.028542)`. A
  `DissTreeNode` shows its id, n and its children's ids. Both printed every
  node with its index arrays and the predictor frame (about 21,000
  characters for a two-predictor tree on 712 cases). The fields are
  unchanged (TX-AUD-20261019-065).
- `disstree_assign` reads each rule's levels as R's parser does: `' & '`
  inside a quoted level is part of the level, so the rules of a tree split
  on a level such as `u & v` assign as R's do (they raised `ValueError`).
  A level may be written in single or double quotes, and its backslash
  escapes are read as R reads them (`\n` is a newline, `\x41` is `A`, and
  likewise octal, `\u` and `\U`). Text R cannot parse, such as a doubled
  quote `''`, an escape R does not know (`\q`) or the rule R's
  `disstree.get.rules` writes for a level that holds a quote, raises
  `ValueError`, as R's `disstree.assign` stops; before, some of it was read
  silently. Checked against R (TX-AUD-20261010-065).
- A root-only tree (no split, for example with `max_depth=1`) has R's `NA`
  rule: `disstreeleaf(label=True)` gives `None` for every case (it gave
  `''`), `disstree_get_rules` gives `[None]` (it gave `['']`), and
  `disstree_assign` assigns no profile to a `None` rule and returns NaN, as
  R returns `NA` (it raised `ValueError`). Checked against R
  (TX-AUD-20261019-047).
- `disstree` and `seqtree` read a predictor's levels as R's `factor()` does:
  values equal to 15 significant digits, such as 0.3 and 0.1 + 0.2, are one
  level and no longer split, and a numeric split's threshold is the level as
  R writes it, read back (1/3 splits at 0.333333333333333). The grouped
  plots (`seqplot` with `group`, `seqHtplot`, `seqedplot`,
  `seqplot_tentrop`) draw such values as one group, named as R names it
  (`1`, not `1.0`); they drew two groups with the same title. Checked against
  R (TX-AUD-20261011-055, TX-REP-20261022-006).
- Numeric groups are named and grouped as R's `factor()` does them in
  `dissassoc`, `disscenter`, `dissvar_grp`, `dissCompare` (and the
  `seqCompare` family's sets), `seqrep_grp` and `dissmergegroups`: a group
  1.0 is named `1` (it read `1.0`), 100000 `1e+05`, and values equal to 15
  significant digits, such as 0.3 and 0.1 + 0.2, form one group (they formed
  two). A whole number beyond R's integer range (2147483647) is a double, as
  in R, so two 64-bit integers equal as doubles form one group. `dissvar_grp` without a group gives one value named `1`, as R; it
  raised TypeError. `seqrep_grp`'s keys are the levels themselves, where R
  prefixes a space. Checked against R (TX-AUD-20261019-046,
  TX-AUD-20261011-004).
- The logit fits of `seqerulesdisc` use R's binomial family beyond a linear
  predictor of 30 in size: there R sets the fitted probability to machine
  precision from 0 or 1 and its derivative to machine precision
  (TX-REP-20261022-004). On separated data traminex's values were up to about
  155 times larger, so the estimates drifted from R's (by 2.8e-6 relative on
  one tested slope); they now agree with R to about 1e-9 or better. Checked
  against R.
- When the binomial fit of a `seqerulesdisc` rule does not converge within 25
  iterations (separated data), it warns "glm.fit: algorithm did not converge"
  and keeps the last estimate, as R's `glm` does; the whole call failed with
  `RuntimeError` (TX-AUD-20261019-085). Checked against R.
- The binomial fits of `seqerulesdisc` also warn "glm.fit: fitted
  probabilities numerically 0 or 1 occurred" when a fitted probability is
  within 10 times machine precision of 0 or 1, as R's `glm` does, converged
  or not (TX-REP-20261022-005). Checked against R.
- Ranks follow R's rule (LINPACK `dqrdc2`, TX-REP-20261019-002): a column is
  set aside when what remains of it, after the columns before it, falls below
  a tolerance times its own length. `dissmfacw` uses R's `qr()` tolerance
  (1e-7) on the weighted design; it compared each column with the largest, so
  a predictor in very large units (1e9) made the others count as redundant
  (rank 1 where R has 3, and NaN statistics). The binomial fits of
  `seqerulesdisc` set columns aside as R's `glm` does (tolerance 1e-11, on the
  weighted design at every step); a nearly collinear design was refused as
  "rank deficient". Their standard errors are now computed from the final QR
  factor, as R's `summary.glm` does, which stays accurate on a nearly
  collinear design.
- `STS_to_TSE` writes no record for a void cell (the end of a shorter
  sequence), so with `state_changes_only=False` (one record per cell)
  `TSE_to_STS` reads a table of sequences of different lengths back with
  their void cells; it read the void marker `%` back as a real state. A
  sequence of void cells only has no record (and is not read back), and a
  table without records keeps its columns. A missing cell keeps its record,
  with the `*` marker. With `state_changes_only=True` the state after an
  internal void starts a new run (TX-AUD-20261019-082).
- The parity ledger now records how `seqefsub` differs from R
  (TX-AUD-20261010-076): a `min_support` that no subsequence reaches gives an
  empty result where R stops with an error, and `weighted=False` leaves the
  caller's weights unchanged where R's `str.subseq` route leaves them all at 1.
  It also records that a `SeqElist` slice holds copies of its sequences, while
  R's slices share them. No behaviour changed; R's results are now recorded
  as goldens.
- Two event sequences (`Eseq`, the elements of a `seqecreate` result) compare
  equal when their id, times, events, dictionary, weight and length are equal
  and, with fractional times, their events were given in an order that
  rounds the gaps alike (the order decides the printed gaps and the counts);
  comparing two sequences of more than one event raised `ValueError`. An
  `Eseq` stays unhashable, as it can be changed in place (TX-AUD-20261011-059).
- `seqeapplysub` uses a numeric `method` only when it equals 1 to 5, as R's
  `method %in% 1:5` does: 2.5 now leaves the constraint's count method (it
  was truncated to 2), and a number written as text, such as `"2"`, now
  counts as that method (it was ignored). Checked against R
  (TX-AUD-20261011-057).
- `seqeconstraint` refuses a fractional `count_method` such as 2.5 with R's
  message " [!] unknown count.method input"; it truncated it to 2. Checked
  against R (TX-AUD-20261010-077). The parity ledger now also records two
  existing differences: for CWIN without a window size R's message names
  'CMINWIN' (a slip) where traminex's names 'CWIN', and a count method given
  as text ("3") is refused, where R accepts it and `seqefsub` then fails.
- `seqerules(sortv="rule")` orders the rules as R does, in R's collation:
  `(a) => (B)` before `(B) => (a)`; it ordered code points, upper case first.
  Numeric sorts are unchanged. Checked against R (TX-AUD-20261019-083).
- `seqecontain` reads the requested events as R does: `1.0` names the event
  `1` (it matched nothing), and a single value may be given without a list (it
  raised TypeError); any other collection, such as a set or a generator, is
  read as its values. Checked against R (TX-AUD-20261019-044).
- `seqerulesdisc` writes the numbers of its rule labels as R prints them, so a
  p-value that rounds to 0.0002 reads `p=2e-04`, not `p=0.0002`; R builds the
  label with `paste(round(...))`. The numbers are rounded as R's `round()`
  rounds, which differs from Python's at decimal ties (6.165 gives 6.16, and a
  p-value of 5e-05 gives `p=0`). Checked against R (TX-AUD-20261019-043).
- `seqdist`'s acceptance of `method`, `subm` and `norm` in any letter case and
  of the `norm` aliases `None`, `'length'` and `'max'`, and `seqMD`'s acceptance
  of `method` and `norm` in the same way, are recorded as a traminex
  convenience: R refuses each of these spellings, and traminex gives R's result
  for the canonical one. The seqdist and seqMD ledger rows and
  `docs/DISTANCES_NORMALIZATION.md` list them, and a test checks them against R.
  Owner decision (TX-AUD-20261019-061, TX-AUD-20261010-074).
- The documents (README, `docs/DEV_BUILD.md`, `seqdef`'s help) give a
  version-proof rule for the state order: without `alphabet=`, labels are
  ordered by the ICU library that is loaded, as R orders them with its own, and
  ICU versions change over time on both sides (a few characters, such as `₨`
  and `₧`, order differently between versions); pass `alphabet=` explicitly for
  an order that matches R's under any ICU version. The test run prints the ICU
  version traminex loaded; a version other than R's is shown, not refused.
  Owner decision (TX-AUD-20261019-053).
- Direct calls to the private compiled module `traminex._pb_kernels` are
  validated: a non-finite or negative indel, a non-finite or negative
  substitution cost between occurring states, a DHD frequency outside [0, 1]
  where the kernel reads it, and input that is not two-dimensional raise
  `ValueError` (an infinite, NaN or negative indel gave 4, NaN and -8; a NaN
  cost was ignored; 1-D input raised `RuntimeError`). The kernels release the
  GIL, so other Python threads run during a long call (2.7% of their idle
  rate before, 79% after with one OpenMP thread). The dispatcher's results
  are unchanged; a zero indel stays valid, as the dispatcher passes one for an
  all-zero TRATE matrix (TX-AUD-20261011-045).
- The raw-array Gower extension (`gower_distance` and the raw-array `FEATURES`
  costs) no longer depends on the platform's `long double`. It relied on an
  80-bit one for weights near the largest double, ranges spanning the whole
  double line and subnormal weights; where `long double` is double it gave
  NaN, NaN and 0.5 for hand values 1, 1 and 0.75. It now computes in double
  with exact power-of-two scaling of each numeric column and of each pair's
  weights, giving the hand values on every platform. Values on x86-64 Linux
  can differ from the earlier ones in the last bits (at most 3 units in the
  last place, 90% of cells identical, on 400 random tables). The STSList
  `FEATURES` path, which follows `cluster::daisy`, is unchanged
  (TX-AUD-20261011-035).
- `docs/DISTANCES.md` and the `JSD_TRANS` kernel's docstring say that a
  sequence with fewer than two positions has no transition and that a
  distance involving one is 0 (TX-AUD-20261011-046).
- `docs/DISTANCES.md` no longer gives fixed figures for the last-bit
  differences between the Numba and reference `JSD_TRANS` and `DUR` chi2
  kernels: the difference grows with the input (8.88e-16 for `JSD_TRANS` on
  20 states and 60 sequences of length 60, above the 4.4e-16 it stated). The
  agreement test gains that input (TX-AUD-20261010-068).
- `docs/OMLOC.md` no longer says every TraMineR norm stays within `[0, 1]` for
  position-dependent OM: `gmean`, `maxdist` and `YujianBo` do, but `maxlength`
  (and `length`, and `auto` for OM) can exceed 1, as in TraMineR; its own
  example gives 1.2 (TX-AUD-20261010-064).
- `docs/NEW_DISTANCE_METHODS.md` says that `EUCLID` ignores `weights`, as R
  does (`seqdist-CHI2.R:6-8`), that the `CHI2`/`EUCLID` `norm` is `"none"` or
  `"auto"` (a logical `norm` is R's deprecated spelling, and its example now
  uses `"auto"`), and what `"auto"` computes: each interval's term divided by
  its largest possible value, then the distance by the square root of the
  number of intervals (it said only the latter). The `norm` refusal for
  `NMS`, `NMSMST` and `SVRspell`, its message included, cites `seqdist.R`
  lines 206-208 instead of 205; `ROADMAP.md` names the latest audit
  (2026-10-19); `docs/INDEX.md` lists the tracked `.codebase-memory/`,
  `.github/` and `graphify-out/` directories and the two pages it left out
  (TX-AUD-20261019-027, TX-AUD-20261011-053).
- `docs/DISTANCES_NORMALIZATION.md` documents `OMloc`'s normalization (R's
  `indel`, the lengths and the corrected bound) and that its raw distance
  depends on `norm`, as in R, since `YujianBo` sets the `maxscost` its local
  indels use; `docs/DISTANCES.md` gives the `WHAM`, `DUR` and `JSD_TRANS`
  formulas, checked against the code by a test (TX-AUD-20261010-069).
- `gmean` of `OMloc`, `OMspell` and `TWED` keeps R's bound, so it is R's value
  and can exceed 1 where R's does. Elzinga's `gmean` stays at least 0 only for a
  bound close to the true maximum, and with the corrected bounds below and R's
  lengths it fell below 0 where R's value was within [0, 1] (`OMloc`, `A`
  against `AAB`, `sm` 4, `expcost` 0, `context` 0.5: -2.46 against R's 1). It
  takes R's normalization indel and R's `maxscost` (the upper triangle of `sm`,
  capped at two indels), so for asymmetric costs the whole-matrix maximum of
  TX-AUD-20261019-008 no longer applies to it. `maxdist` and `YujianBo` keep the
  corrected bounds. Owner decision (TX-REP-20261022-002).
- Normalized `TWED` distances (`maxdist`, `YujianBo` and so `auto`) no longer
  exceed 1, up to rounding in the last place where the bound is reached.
  R's bound caps the largest substitution cost at two indels,
  but TWED pays the costs directly, so with a small `indel` its normalized
  values exceeded 1 (`maxdist` 1.58, `gmean` 3.91); against an empty sequence
  R's distance `n * indel` exceeded its bound (1.667). The bound now takes the
  uncapped largest cost, and is the empty-sequence distance itself there. Raw
  distances are R's; `gmean` keeps R's bound (TX-REP-20261022-002). Recorded as
  a divergence from R (TX-REP-20261021-005).
- `OMspell` refuses data in which a spell's indel cost
  `indel + expcost * (t^tpow - 1)` is negative, which happens with `tpow`
  below 0 once `expcost * (1 - t^tpow)` exceeds the indel. R computes on and
  returns negative distances there (`AABBAABB` against `BBAABBAA`, `sm` 4,
  `indel` 0.3, `expcost` 3, `tpow` -2: -15.6, `YujianBo` 2.89), which no
  normalization bounds; traminex refuses them as it refuses negative indel and
  substitution costs. A zero spell cost stays valid. Recorded as a divergence
  from R (TX-REP-20261022-003).
- Normalized `OMspell` distances (`maxdist`, `YujianBo` and so `auto`) no
  longer exceed 1, up to rounding in the last place where the bound is
  reached. R's bound leaves out the spell-duration costs, so its
  normalized values exceeded 1 (4.86 at `expcost=1, tpow=2`; 9 at `expcost=0`
  when the substitution cost is below the indel). The bound is now the largest
  cost, over every spell structure of the two lengths, of pairing spells in
  order and indelling the rest; with R's defaults it is R's. Raw distances are
  R's; `gmean` keeps R's bound (TX-REP-20261022-002). Recorded as a divergence
  from R (TX-REP-20261021-004).
- Normalized `OMloc` distances (`maxdist`, `YujianBo` and so `auto`) no
  longer exceed 1, up to rounding in the last place where the bound is
  reached. R's bound prices an indel at `expcost * max(sm) + context`,
  but OMloc's local indel costs up to `expcost * maxscost + context * max(sm)`,
  so with `expcost` below 0.5 normalized values exceeded 1, with R's `CONSTANT`
  costs too (1.25 for A against BBB at `expcost=0.25`). The bound is now the cost
  of an alignment that always exists, priced with those largest costs; raw
  distances are R's; `gmean` keeps R's bound (TX-REP-20261022-002). Recorded
  as a divergence from R (TX-REP-20261021-003).
- `NMS`, `NMSMST` and `SVRspell` give identical sequences a distance of 0 when
  their counts overflow, as R does; they returned NaN. Between different
  sequences, where R's double sums overflow to NaN or Inf (one-state sequences
  of 515 positions, `tpow=150` on spells of 50), the distance is recomputed
  with counts rescaled by powers of two and is finite unless it is itself
  beyond the double range; NumPy's overflow warnings no longer reach the caller.
  Values that R computes finitely are unchanged (TX-AUD-20261019-029,
  TX-AUD-20261011-047).
- `NMS`, `NMSMST` and `SVRspell` add their weighted subsequence counts in R's
  order, one length after another and skipping zero weights; `np.sum` rounded
  differently, so the distances differed from R's in the last bit
  (TX-REP-20261022-001).
- Spell-length powers (`OMslen`'s `h`, the `tpow` of `OMspell`, `NMSMST` and
  `SVRspell`) are computed as R's `^` computes them, the C library's `pow`;
  `np.power` differed from it in the last bit for some lengths, so `OMslen` and
  `OMspell` distances differed from R's in the last bit (TX-AUD-20261019-030).
- NMS and SVRspell with a proximity matrix holding negative entries are tested
  against R (the full matrix and `refseq=2`, new generator
  `make_goldens_negative_prox.R`), and every row-engine method is tested to
  evaluate at most the pairs a `refseq` returns, with values equal to the full
  matrix's or, with asymmetric costs, to the per-pair reference's
  (TX-REP-20261019-001).
- `docs/DISTANCES.md` says that only the row engine evaluates just the pairs a
  `refseq` call returns; the dense Numba OM, LCS, LCP and RLCP kernels, DHD and
  CHI2/EUCLID compute the full matrix (TX-AUD-20261019-081).
- `seqdist(method="OMstran", transindel="prob")` on data without a transition
  (sequences of length 1) is refused, as transition-rate costs are; it computed
  with empty rates. The other degenerate OMstran and CHI2/EUCLID inputs on which
  R stops or returns a wrong shape are now recorded as divergences, with
  traminex's values tested: OMstran's from its definition, and CHI2/EUCLID's
  with the first overlapping episode capped as R caps the others, one of two
  possible corrections (TX-AUD-20261019-080).
- `seqMD`'s combined sequence object uses seqdef's default markers (`*` and `%`)
  as R's does, not channel 1's, and the deprecated `what="sm"`/`"seqmc"` warn
  with R's text as a `DeprecationWarning` (they were a `UserWarning`). Zero
  channel weights with `link="mean"`, which R turns into NaN costs, stay refused
  and are recorded as a divergence (TX-AUD-20261010-074).
- `seqdiff` no longer prints seqdist's " [!!] " notices ("with.missing set as
  FALSE", "Empty sequence(s)", and for CHI2/EUCLID with `breaks` the notice
  that `step` or `overlap` is ignored) once per sliding window; R's seqdiff
  silences them. seqdist's notice that `step` or `overlap` is ignored when
  `breaks` is given now has R's text (" [!!] 'step' has been ignored because
  'breaks' is specified") (TX-AUD-20261019-059).
- An invalid `TRAMINEX_NUMBA` value is refused at import whether or not numba
  is installed; without numba it was accepted silently
  (TX-AUD-20261011-049).
- `TRAMINEX_BINDING` is read as `TRAMINEX_NUMBA` is: surrounding spaces are
  ignored and an empty value means `auto`; `" python "` and `""` were refused
  at import (TX-AUD-20261019-031).
- A NumPy scalar or 0-d array substitution cost (`subm=np.float32(2)`,
  `np.int64(2)`, `np.array(2.0)`) gives the distances of the Python number; it
  raised IndexError on most routes. A boolean or complex `subm` is refused, as
  R refuses a logical or complex `sm`, at R's `sm` check (after `norm` and
  `indel`); a boolean was used as 1 and a complex number as its real part
  (TX-AUD-20261019-028).
- The configuration objects agree with `seqdist`: `DHDConfig` accepts the
  `variant` spellings `seqdist` accepts (`"MEAN"`), and an option nothing would
  read is refused with `TypeError` instead of being ignored: a
  `calibrated_params_path` without a calibrated profile (in `seqdist` and the
  low-level `pairwise` too), and `OMConfig(trate_config=...)` with a `subm`
  other than `"TRATE"` (TX-AUD-20261019-026).
- `seqRLCP` no longer prints seqdist's " [!!] seqdist: 'with.missing' set as
  FALSE" notice on every call with data that hold no missing cells; it passes
  `with_missing=True` only when they do (TX-AUD-20261019-025).
- The low-level `pairwise(..., method="OM", subm="TRATE")`, which
  `docs/DISTANCES.md` documents, also works on rows with void cells and with a
  per-state indel vector, with the costs `seqdist` resolves for the same data;
  those routes refused it with the message meant for the other methods
  (TX-AUD-20261019-024).
- `seqdist(method="DHD", subm="CONSTANT")` is refused with R's message
  ("'sm = \"CONSTANT\"' is not relevant for DHD, consider HAM instead"), and an
  invalid `norm` for HAM no longer says that `gmean` and `maxdist` are not
  implemented: both compute, as in R. `seqdist(method="OMstran")` with generated
  costs that are not finite (INDELS with two states the data never show) says
  that the substitution costs are not finite; it said the indel cost should be
  positive (TX-AUD-20261019-079).
- `check_costs` is read as R's `if (check.costs)` reads it, in `seqdist`,
  `seqcost`, `seqsubm`, `seqMD` and `seqdistmc`: the strings `"FALSE"`, `"F"`,
  `"false"` and the number 0 are R's FALSE and refused as `False` is; other
  strings, NaN (R's NA) and `None` (R's NULL) are refused with R's messages.
  `"FALSE"` and `"no"` used to run the checks as if true, and `None` was
  reported as `check_costs=False`. A 0-d array is read as the value it holds
  (`np.array(True)` is true). R reads the value only where it checks costs, and
  computes whatever it is for LCS, CHI2/EUCLID, LCP/RLCP, NMS/NMSMST/SVRspell,
  HAM with a matrix, HAM's and DHD's default costs and `seqMD(what="MDseq")`;
  traminex reads it on every route, as it refuses `False` on every route
  (TX-REP-20261019-009).
- `seqdist` checks its arguments in R's order (norm, indel, sm, prox, link,
  step), so with two invalid arguments it names the one R names; it named
  `link` where R names `norm`, and `prox` where R names `sm`. A `norm` spelled
  with surrounding spaces is refused at the `norm` check, as R refuses it
  (TX-REP-20261019-010).
- At the declared matplotlib floor (3.6.0) every legend keeps a label that
  starts with an underscore, such as a state or group named `_x`: matplotlib
  3.6.0 drops such a label even with explicit handles, so the legends of
  `state_distribution_plot`, `seqsurvplot`, `seqdiffplot`, `seqedplot`,
  `seqHtplot` and `seqplot_tentrop` lost it there. Every legend is now made
  with stand-in labels and then given the real text (TX-AUD-20261019-060).
- Sixteen public functions that take a further required argument now check
  their object argument first and raise `TypeError` naming the expected type
  where a plain array raised `AttributeError`: `seqMD_gridsearch_weights`,
  `seqMD_optimize_weights`, `seqrecode`, `seqimplic`, `seqsha`, `seqstart`,
  `cluster_representatives`, `cluster_summary`,
  `frequent_event_pairs_timegap`, `seqtree_cv_fit`, `seqtree_fit`,
  `seqtree_predict`, `seqtree2dot`, `disstree2dot_explore`, `disstree2dotp`
  and `disstree_predict` (TX-AUD-20261019-036).
- Without Numba (the core install, or `TRAMINEX_NUMBA=off`), dense LCS runs
  through the batched row engine, as LCS on sequences of unequal length
  already did, instead of one pure-Python dynamic program per pair: 100 `mvad`
  sequences of 70 positions took 0.48 s instead of 15.1 s, with identical
  results. The test guarding this now runs with Numba dispatch both on and off
  in the ordinary suite; it had failed only under `TRAMINEX_NUMBA=off`
  (TX-AUD-20261019-022).
- The workflow validator (`scripts/validate_workflows.py`) reads
  backslash-continued lines joined, so a `|| true` on a continuation line no
  longer passes, and counts `-m`, `--ignore`, `--ignore-glob` and test paths as
  test selection in a gate's test command, as it counted `-k` and
  `--deselect`; a gate writes an option's value attached (`--option=value`).
  The test command is recognised whatever runs it (`$PYTHON -m pytest`,
  `coverage run -m pytest`, `env VAR=value pytest`), and `PYTEST_ADDOPTS`,
  options with which pytest runs no test (`--setup-plan`, `--fixtures`,
  `--version`, ...) and options that rewrite its configuration (`-o`,
  `--override-ini`, `-c`, `--config-file`) are refused (TX-AUD-20261019-021).
- The legends of `sequence_index_plot(..., legend=True)` and
  `sequence_index_plot_grouped` use the colours the cells are drawn with, as
  R's `seqIplot` legend takes a user `cpal`; with a user `cmap` they showed the
  states in the object's palette (TX-REP-20261021-002).
- A user colour list (`cpal`, or `cmap` for the index and distribution plots)
  with one colour per real state, R's convention, now draws missing cells in the
  object's missing colour, as R appends `missing.color` after a user `cpal`; a
  list with one more colour draws them in that colour, and void cells stay
  white. The list was recycled over every code, so `seqfplot`,
  `sequence_index_plot`, `state_distribution_plot` (cells and legend) and
  `seqmsplot` drew missing cells in the first state's colour
  (TX-AUD-20261019-020).
- A categorical `sortv` given as a pandas Series of category dtype orders the
  index plot by its categories, as R orders a factor by its levels, on every
  route, `seqplot(..., group=)` included; it was ordered as text, and the grouped
  route also lost a `pd.Categorical`'s order (TX-AUD-20261019-019).
- `seqrfplot` (and `seqplot(type="rf")`) draws each box of the right-hand panel
  inside its group's band, with R's widths (`0.8 * h / max(h)` on R's axis); the
  boxes were several bands tall and overlapped (TX-AUD-20261019-018).
- `dissrep` and `seqrep` order candidates whose scores are equal in exact arithmetic
  by row, as R's `order()` does on a tie, on every build: such scores are summed
  again with correct rounding, so the summation order of the linear-algebra
  library no longer decides which comes first (TX-AUD-20261019-017).
- `seqentrans(avg_occ=True)` after `seqefsub(..., weighted=False)` on weighted
  sequences divides `Count` by the number of sequences, so `Avg.occ` is the
  average number of occurrences by sequence; it divided by the total weight, as
  R does (a tenth of the average with weights 10). `seqentrans` also returns a
  copy and leaves its argument unchanged, as R does; it added the columns to the
  argument. On a DataFrame without a positive `seqe_weight_total` attr, its message
  says the caller sets that attr to the total the counts were taken on; it named
  `seqefsub`, which sets no such attr (TX-AUD-20261019-016).
- `createdatadiscrete` accepts fractional ages as TraMineRextras does (rows
  `trunc(end - start + 1)`, ages `start, start + 1, ...`, an independent event
  from row `trunc(age - start + 1)`); it refused them (TX-AUD-20261019-015).
- `STS_to_SPELL` reads the object as it is, as TraMineR's internal
  `STS_to_SPELL` does, so trailing nonresponse cells (`seqdef(..., right=None)`)
  form a spell; it re-ran `seqdef(right="DEL")` first and dropped them. Pass
  `right="DEL"` for the previous output, which `seqformat(to="SPELL")` and
  `seqsurv` keep, as in R (TX-AUD-20261019-014).
- `seqintegr`, `seqibad`, `seqipos(index="integr")` with `include_missing=True` and
  `seqindic`'s `Integr` and `Bad` columns with `with_missing=True` leave void positions out of the
  integration weights: a short sequence's value no longer depends on the length
  of the other sequences (a-b padded with two void cells gave 0.1 and 0.2, now
  1/3 and 2/3). R counts the padding; recorded as a divergence
  (TX-AUD-20261019-012).
- The low-level `pairwise` refuses a NaN or infinite substitution cost between
  observed states on every route (data with void cells, a state-vector indel,
  `OMslen`, `TWED`), as its dense `OM` route already did; those routes computed
  silently (TX-AUD-20261019-010).
- With an asymmetric substitution-cost matrix, `pairwise(ref_sets=)` on data
  with void cells, and `pairwise(ref_sets=)` and `pairwise_between_sets` when the
  second set repeats a sequence of the first, return the distance with the first
  set's sequence first, as R's `refseq=list(a, b)` and `seqdist` do; they
  returned the other direction (TX-AUD-20261019-009).
- With an asymmetric substitution-cost matrix whose largest cost lies below the
  diagonal, the normalization bound of `OM`, `OMslen`, `OMspell` and `TWED`
  distances (`maxdist`, `gmean`, and `pairwise`'s `norm="max"`) takes the
  largest off-diagonal cost of the whole matrix. It took the largest cost of the
  upper triangle only, as R does, so it was not the maximum possible distance it
  stands for and normalized distances exceeded 1. That cause is removed; the
  other causes, in `OMloc`, `OMspell` and `TWED`, are corrected separately
  (TX-REP-20261021-003 to -005); `gmean` of `OMspell` and `TWED` keeps R's
  upper-triangle bound (TX-REP-20261022-002). Recorded as a divergence from R
  (TX-AUD-20261019-008).
- `seqtab`'s `idxf` attribute is R's one-based row number, beside the
  one-based ranks in `idxs`; it was zero-based. `source_indices`, which has no R
  counterpart, still gives zero-based rows of the input for indexing
  (TX-AUD-20261019-034).
- `seqindic` names its rows by the sequence ids, as R names them by
  `rownames(seqdata)` (`1..n` for an object without ids). They were numbered
  `0..n-1`, so the result could not be joined back to the sequences by id
  (TX-AUD-20261010-093).
- An event sequence whose id is -1 is counted by `seqefsub` (and the other
  miners) and printed with its gaps. The engine used -1 as its "no sequence
  yet" and "generic pattern" markers, as TraMineR's C++ does, so such a
  sequence was never counted (an empty result) and printed like a pattern; R
  stops on this input. Recorded as a divergence (TX-AUD-20261019-006).
- A substitution matrix given as nested lists or tuples is read as the
  array it spells on every route. On sequences of different lengths an
  asymmetric matrix so given was applied in the reverse direction and without
  R's "Non-symmetric substitution cost matrix!" warning (TX-AUD-20261019-003).
- With `subm="INDELS"`, an `alphabet=` state that no compared sequence
  contains no longer changes any distance. Its stand-in cost (`DBL_MAX`, as in
  R) set the normalization bound, so every normalized OM, HAM, DHD, OMloc,
  OMslen, OMspell and TWED distance was 0 or NaN, OMloc's and TWED's raw
  distances were inflated, and OMstran's raw and normalized distances changed
  (it rescales its costs by the largest indel cost); R does the same, and the
  correction is recorded as a divergence. The costs `seqcost` returns are unchanged (TX-AUD-20261019-002).
- `seqdist(refseq=0)` chooses among equally frequent sequences as R does:
  by the text of R's `seqtab(seqnum(...))`, with zero-based state codes and a
  missing cell written `*`. The key used one-based codes, so with ten or more
  states, or a missing cell, another reference sequence could be chosen and
  every returned distance changed (TX-AUD-20261019-001).
- `seqdef_from_long` and `TSE_to_STS` (state mode) keep a real state spelled
  like the missing marker (`*`) as a state. They read it as the marker, so it
  merged with genuinely missing cells and sequence lengths changed; the
  marker now gets another spelling when a state uses its text, as on
  `seqdef`'s SPS and SPELL routes (TX-AUD-20261019-005).
- `seqdef_from_long` and `TSE_to_STS` (both modes) refuse `nr` equal to
  `void` (including `nr="%"` with the default `void`), with `seqdef`'s
  message "'nr' and 'void' must be distinct character markers". They
  accepted it and silently wrote one of the two markers as `<VOID>`
  (TX-REP-20261019-012).
- `seqpolyads` compares its roles as R does: by `alphabet(x, with_missing)`
  only (seqpolyads.R:52-57), recoding each role onto one codebook before
  stacking. It refused a role list in which only some roles held a
  missing cell, with `with_missing=True`, where R computes
  (TX-AUD-20261019-004).
- **OMslen's `link` is checked before computing (audit TX-AUD-20261011-051).**
  An invalid `link` (and, through `seqdist`, `link="none"`, which R refuses for
  OMslen) was accepted when there was one sequence or every row was identical,
  because it was checked only when a pair was computed. It is now checked first,
  at every layer, as R checks it.
- **Low-level OM takes a per-state indel without void (audit
  TX-AUD-20261011-050).** `traminex.distances.pairwise(..., method="OM",
  indel=<vector>)` raised `TypeError: only length-1 arrays can be converted to
  Python scalars` unless the matrix held a void code, although that layer
  documents per-state indels and its void route and `seqdist` computed them.
  The dense route now computes them too.
- **The distance configuration dataclasses take `seqdist`'s values (audit
  TX-AUD-20261011-043).** `OMConfig` refused `indel="auto"` and per-state
  indels with raw NumPy errors, accepted a negative cost matrix, and fixed
  `indel=1.0` even for `subm="TRATE"`, where `seqdist`'s automatic indel
  differs; `OMConfig` and `LCSLCPConfig` refused norms `seqdist` accepts
  (`maxlength`, `gmean`, ...), and `TRATEConfig(smoothing="x")` raised a raw
  NumPy error. They now take `seqdist`'s values (`OMConfig`'s indel defaults to
  `"auto"`; see Changed), refuse what `seqdist` refuses, and name the field in
  each refusal, also for a cost that is not a number (`OMConfig(subm=None)`).
- **CHI2/EUCLID windows must be integer positions (audit TX-AUD-20261011-048).**
  `seqdist(method="CHI2")` truncated `step=2.5` and window positions such as
  `breaks=[(2.7, 5)]` with `int()` (R refuses the step and, for `(2.7, 5)`,
  silently uses positions 2-4). A fractional step and any `breaks` pair that is
  not two integer positions are now refused, as unordered and out-of-range pairs
  already were; R's documentation defines breaks as ordered pairs of integers.
  The low-level `pairwise` applies the same step rule (it truncated `step=2.5`
  to 2), and `step=0` without `breaks`, which both layers read as the whole
  sequence, is refused, as R stops on it.
- **Reference distances with asymmetric costs follow R's direction (repair
  finding TX-REP-20261017-001).** With an asymmetric substitution matrix,
  `seqdist(..., refseq=k)` and an external `refseq` gave d(reference, x) for
  the sequences after the reference, and index sets whose first set came after
  the second gave d(set-2, set-1); R gives d(x, reference) and d(set-1,
  set-2). The reference routes now evaluate each pair in R's direction, also
  with `deduplicate=True`, and a normalized reference distance divides that
  oriented distance by R's bound.
- **`seqdist` reads `refseq` as R does (audit TX-AUD-20261011-044).**
  `refseq=[1, 2]` was read as two singleton index sets and returned a 1x1
  block, where R stops on the vector `c(1, 2)`. A flat list, tuple or array of
  numbers is now refused; two index sets are two collections (R's
  `list(1, 2)` is `([1], [2])`). An index out of range or negative raised
  `IndexError`, while index sets raised `ValueError`, and a string, set,
  mapping or other object `TypeError`; every `refseq` refusal is now a
  `ValueError`. A NumPy integer such as `np.int64(2)` was taken for a
  sequence object and failed; it is now one index.
- **`seqdist` refuses transition-rate costs it cannot compute (audit
  TX-AUD-20261011-039, TX-AUD-20261011-040).** With a one-state alphabet or
  sequences of length 1, `seqcost(..., "TRATE")` and `seqtrate` already
  refused (and R stops), but `seqdist(method="OM", subm="TRATE")` returned
  zeros and `seqdist(method="DHD")` returned numbers. Both now refuse with
  the same message as `seqcost` and `seqtrate`.
- **`seqdist` gives R's empty-sequence and `with.missing` messages (audit
  TX-AUD-20261011-042).** R warns "Empty sequence(s) ... may cause inconsistent
  results!" for a sequence with nothing to compare, and says when it resets
  `with.missing` to FALSE because nothing is missing; traminex was silent.
  Both are now `UserWarning`s with R's text, given once per call on every
  `refseq` route, and OMloc's refusal of an empty sequence names it, as R's
  does.
- **Generated cost matrices get R's cost messages (audit TX-AUD-20261011-080).**
  TraMineR 2.2-14 runs `checkcost` on every matrix `seqcost` generates; traminex
  gave its messages only for a matrix passed explicitly. `seqcost` and `seqsubm`
  now print them (a triangle-inequality violation of TRATE costs, a largest cost
  of 0 or below, zero costs, a cost above twice the indel), and `seqdist` prints
  them once per call, on every `refseq` route too, for a cost-method string, for
  LCS's and HAM's constant costs and for DHD's generated costs, as R does (an
  index or index-set reference printed every message twice, batch-6 review);
  `seqMD` silences them for each channel, as R does, and keeps its check of the
  combined costs.
- **`check_costs=True` is accepted (audit TX-AUD-20261011-041).** `seqdist`,
  `seqcost`, `seqsubm`, `seqMD` and `seqdistmc` refused `check_costs` at every
  value, while the ledger said only `False` was refused. `True`, R's default,
  is now accepted; `False` is still refused, with the reason: R's
  `check.costs = FALSE` computes with a negative indel or negative costs and
  returns negative distances, so traminex keeps its cost checks mandatory.
- **`seqdist` accepts R's arguments for every R method (audit
  TX-AUD-20261010-071).** R's `seqdist` takes all of its arguments for every
  method and ignores the ones a method does not use; traminex refused them
  (`TypeError: Unsupported kwargs`), so a call giving R's defaults explicitly,
  such as `tpow=1` or `link="mean"` with OM, failed for every method. Each R
  argument (traminex spelling) is now accepted for each of R's 16 methods: the
  checks R makes for every method (norm, indel, sm, prox, link, step) apply to
  every method, and an argument the method does not use is ignored, as in R.
  A number or vector `subm` is refused for a method without costs, as R
  refuses it (a number stays traminex's extension for the methods with costs).
  Unknown keywords still raise. `opt_args` with a `tokdep.coeff` entry,
  TraMineR's undocumented OMtspell route, is refused; any other `opt_args`,
  such as a number, is ignored, as R ignores it.
- **`seqsubm` shows its real `method` default (audit TX-AUD-20261010-061).**
  The signature and the stub showed `method = "TRATE"`, but a sequence object
  without `method` raised, as R's `seqsubm` does (R's `seqcost` has no default
  for it); only a raw integer array used `"TRATE"`. The default is now `None`:
  required for a sequence object, `"TRATE"` for a raw array.
- **The relative-frequency MDS axis has a fixed orientation (repair finding
  TX-REP-20261017-006).** `seqrf`, `dissrf`, `seqplot_rf`, `seqrfplot`,
  `seqplot(type="rf")` and `seqdomassoc(rep_method="eq.group")` sort sequences
  by the first MDS coordinate, whose sign R takes from LAPACK, so their groups,
  medoids and R2 followed the numerical library (on one input R and traminex had
  opposite signs). The axis is now oriented so that its coordinate of largest
  absolute value is positive; coordinates within a relative 1e-10 of it count
  as tied and the first decides, so symmetric data, whose tied coordinates
  differ only in rounding, give the same groups on every build. The results
  equal R's `seqrf` given the axis oriented by that rule.
- **`seqdim` follows R's rule for compressed sequences (audit
  TX-AUD-20261011-071).** It looked only at the first element:
  `seqdim(["AAB", "A-B-C"])` gave `[1, 2]` where R gives `2 3`, a leading
  missing value hid the compressed rows, a single `"A-B-C"` raised, and a
  2 x 2 character matrix with a dash gave `[4, 2]` where R gives `2 2`. It now
  applies R's `seqfcheck` over every element, as R does, and reads numbers
  mixed with text as text, as R's `c()` makes them (`[1, "A-B"]` gives `2 2`).
- **Generator-authored R-golden records are labelled (audit
  TX-AUD-20261011-020).** Some fixture records are written by the generator
  rather than taken from an R call: R statements re-implemented because the
  function only draws (plot.stslist sampling, seqedplot series, seqsplot
  panels), hand-typed messages, traminex-side metadata, a per-group
  computation, mapped error codes and the MDS axis oriented by traminex's rule.
  `scripts/r_goldens/AUTHORED_RECORDS.csv`
  lists each with its reason, the README rule points to it, and a test keeps
  the list true to the generators and fixtures.
- **The surface the sibling projects import is documented (audit
  TX-AUD-20261010-048).** weightedcluster and py_seqplot import 19 traminex
  objects by module path, 9 of them outside the public manifest
  (`pairwise`, `pairwise_between_sets`, the `DissTree` classes, palette, series
  and ticks helpers), and no document said so. `docs/CONSUMED_SURFACE.md` lists
  them, and a test requires every listed path to resolve. No sibling is changed.
- **The declared scikit-build-core floor builds the sdist (audit
  TX-AUD-20261010-047).** `pyproject.toml` requires `scikit-build-core>=0.11.0`
  but set `sdist.inclusion-mode`, which 0.11.0 does not know (exit 7). The option
  is now set only for scikit-build-core 0.12 and later, through an override,
  and every root entry is anchored (`/LICENSE`, which 0.11.0 would otherwise
  match at any depth); 0.11.0, 0.12.2 and 1.1.0 build the same sdist, and the CI
  workflow is configured to build it with 0.11.0 (never run: no git remote).
  With 0.12 and later, which read no `.gitignore` file,
  the sdist also builds from a checkout holding a malformed ignored
  `.gitignore`; 0.11.0 reads them, so docs/DEV_BUILD.md says to build with it
  from a clean clone. The floor is unchanged.
- **`make release-check` names its prerequisites (audit TX-AUD-20261010-046).**
  In an interpreter without the development install the gate stopped with a
  bare `ModuleNotFoundError` (`tokenize_rt`, then `traminex`). It now first runs
  `scripts/check_dev_environment.py`, which lists each missing module, the gate
  that needs it and the documented install command; docs/DEV_BUILD.md states
  the prerequisite.
- **ICU documentation and portable library names (audit TX-AUD-20261010-030).**
  The documents listed `seqtab` among the ICU users only for tied
  frequencies; it needs ICU whenever it orders two or more labels. They now
  also say that the order follows the loaded ICU's major version (conda's
  ICU 78 orders `₧ ₨` where ICU 74 and R give `₨ ₧`) and how to see which
  library is loaded. The loader reads the ICU version from `.so.N` and
  `.N.dylib` names, and names the library it cannot read (an `icuinN.dll` form,
  first added here, could not be reached and was removed; see
  TX-REP-20261019-005 above).
- **A wrong-type input raises `TypeError` (audit TX-AUD-20261010-029,
  TX-REP-20261018-001).** Called with a plain array, sixteen public functions
  raised `AttributeError` from deep inside (`seqidegrad`, `seqibad`, `seqipos`,
  `seqMD_distances`, the event helpers `events_to_state_sequences`,
  `frequent_events`, `frequent_event_pairs`, `events_transition_matrix`,
  `events_transition_degrees` and `seqentrans`, and the tree functions
  `disstreedisplay`, `seqtreedisplay`, `disstree_plot`, `seqtree_plot`,
  `disstree_leaf_summary` and `disstree_feature_importance`). They now raise
  `TypeError` naming the object they expect, as the other functions do.
- **`seqsubsn` returns a 1-D array (audit TX-AUD-20261010-029).** It returned
  shape `(n, 1)` while every other per-sequence indicator returns a length-n
  array (docs/API_CONSISTENCY.md section 1c; R returns `n x 1` for all of them).
  Code that indexed `result[i, 0]` now uses `result[i]`.
- **docs/API_CONSISTENCY.md states the missing-state keyword rule (audit
  TX-AUD-20261010-029).** Section 1 was headed "standardized", while 35 public
  functions take `include_missing`, 23 take `with_missing` and `seqtabstocc`
  takes both. It now lists each group, records the rename in section 1b, and a
  test compares the lists with the signatures.
- **`seqnum` has an R golden (audit TX-AUD-20261010-021).** Its ledger row
  claimed exact parity, but no generator called R's `seqnum`. A generator now
  records R's result for a non-sorted alphabet, numeric states, and missing and
  void cells with `with.missing` FALSE and TRUE; traminex matches it.
- **The workflow validator rejects a neutralised CI (audit
  TX-AUD-20261010-057).** A `ci.yml` left with only a manual trigger and its test
  commands turned into `--collect-only || true` passed validation. ci.yml must
  now run on push and pull_request with filters that keep the default branch,
  and a gate step (a make target, pytest, `run_pybind11_tests.sh`) may not
  follow its command with `||` or `; true`, feed a pipe without pipefail, run
  `set +e` or under another shell, only collect or select tests (`-k`,
  `--deselect`), set `continue-on-error` to anything but false, or carry a
  condition other than the matrix backend, nor may its job carry one.
- **The static R-golden provenance check rejects an orphan swap (audit
  TX-AUD-20261010-056).** Replacing a real fixture with an orphan file of the
  same count passed, because a three-character literal such as `"diss"` (from
  `what = "diss"`) "named" the orphan as a prefix. A fixture name a writer
  spells in full must now exist, and a prefix names a fixture only when it looks
  like the start of a built name; every fixture stays attributed. An orphan
  that keeps a built name's prefix still passes the static check, which the
  README now says; `--execute` rejects it.
- **The R-golden provenance proof no longer rewrites the checkout (audit
  TX-AUD-20261010-020).** `scripts/check_r_golden_provenance.py --execute`
  ran every generator in the checkout holding the script when `--repo` was
  omitted, so a byte change overwrote tracked fixtures before the check failed.
  It now refuses that checkout (also through a symlinked path) and any git
  work tree with uncommitted fixture changes, before any R process starts,
  unless `--allow-in-place` is given; the README runs it on a clone with
  `--repo`.
- **R-golden generators pin their locale (audit TX-AUD-20261010-019).** Run
  from a shell in the C locale, the generators changed 24 tracked fixtures:
  text sorted by bytes, non-ASCII literals read as bytes even where the ICU
  root collation was set, and R's messages quoted with `'` instead of `‘ ’`.
  Every generator now starts by switching to `C.UTF-8` with the ICU root
  collation and R's messages in English (stopping if that locale or ICU is
  missing); the locale alone had left R's messages in the shell's language, so
  a `LANGUAGE=de` shell still wrote 45 fixtures in German. `--execute` runs the
  scripts with `LC_ALL=C.UTF-8` and `LANGUAGE=en` and records them. All 175
  generators now write byte-identical fixtures from shells with `LC_ALL=C`,
  `LC_ALL=C.UTF-8`, `LANG=fr_FR.UTF-8` and `LANGUAGE=de`.
- **The type-check gate works with the newest permitted NumPy (audit
  TX-AUD-20261010-010).** With NumPy 2.5 (which requires Python 3.12), `make
  typecheck` exited 2 in NumPy's own stubs, because mypy was fixed at a Python
  3.11 target, and a 3.12 target found two annotation errors. mypy now checks
  for the interpreter it runs under, the two errors are fixed, and a CI job is
  configured to type-check with the newest NumPy, pandas and SciPy the floors
  permit, on Python 3.11 and on 3.12 (never run: no git remote); the 3.11 row still rejects a standard-library call
  new in 3.12, such as `itertools.batched`. The workflow validator requires the
  job and both rows. No dependency floor changes.
- **The stub gate checks re-exports (audit TX-AUD-20261010-009).** A root stub
  that re-exported another function under a public name (`seqtransn as seqtab`)
  passed every stub check, because only names and signatures were compared.
  `scripts/check_stub_signatures.py` now resolves each public re-export to its
  source module and original name and requires the runtime object to be that
  object; it reports how many re-exports it compared.
- **`seqpcplot` drops the missing cells of a state sequence object (batch-4
  review of audit TX-AUD-20261011-027).** As R's `seqdss` does, missing cells
  are now removed before runs are collapsed, so `A*ABC` is drawn as A-B-C
  rather than A-*-A-B-C, and no missing state appears on the y axis, which now
  shows the state labels instead of the alphabet codes (the alphabet for an object
  built without labels, as R's `seqdef` labels default to the states).
- **Grouped representative and per-type plots take R's options (batch-4 review
  of audit TX-AUD-20261010-044).** `seqrfplot(x, diss, group=g)` raised on its
  default call, because its two-panel default reached every panel. With `group`
  each level now gets its medoid panel, as R draws grouped plots, and
  `which_plot="both"` with `group` is refused as in R; the default is now
  `None`, which still draws both panels without `group`. `seqfplot`,
  `seqmtplot`, `seqmsplot`, `seqdHplot` and `seqrplot` take R's `main`
  (default `"auto"`), which raised `TypeError` before (for `seqrplot`, whenever
  `diss` was given). A per-sequence `weights` vector is split per panel, where
  `seqdHplot(group=g, weights=w)` raised.
- **Legends keep states, groups and events whose label starts with "_"
  (TX-REP-20261017-003).** Matplotlib's automatic legend leaves out such
  labels, so a state named `_x` vanished from the legends of the distribution
  plot (`seqplot(type="d")`), `seqsurvplot` and `seqdiffplot`, and an event
  named `_x` from `seqedplot`'s. They now pass their entries explicitly, as
  `seqHtplot` does, and list every state, group or event.
- **`transition_heatmap` draws `seqtrate`'s rates (audit TX-AUD-20261011-028).**
  The heatmap counted transitions over every code, the missing and void markers
  included, and was always weighted, although its documentation called the
  rates `seqtrate`'s (row B `0.24/0.45/0.12/0/0.18` where `seqtrate` gives
  `0.3/0.55/0.15`). It now draws `seqtrate(x, weighted=weighted)`, or its
  counts with `normalize=False`, over the states `seqtrate` reports, and takes
  a new `weighted=True` keyword.
- **`dissCompare` draws the larger group first (audit TX-AUD-20261011-022).**
  With sampling (`s > 0`) and the smaller group as the first level, the first
  group's sample took cases of the other group, so the result depended on the
  order of the group levels; R's `dissCompare` does the same (mean sampled LRT
  5.8 against 9.6 for the same split), while R's `seqCompare` puts the larger
  group first. `dissCompare` now draws the larger group first, so the same seed
  gives the same result for either order. This is a recorded correction of R.
- **`seqMD_optimize_weights` never returns all-zero weights (audit
  TX-AUD-20261011-018).** The search could set every channel weight to 0,
  report convergence, and return weights that `seqMD_distances` then refused
  ("weights must have positive sum"). A zero-sum weight vector is no longer a
  candidate, and the weights and score are returned as Python floats.
- **`seqerulesdisc` keeps the valid rules when one rule's predictor is constant
  (audit TX-AUD-20261011-010).** A constant predictor made the whole call fail
  ("binomial model design matrix is rank deficient"). R returns every rule but
  reports the age coefficient under the constant predictor's name. The
  predictor now has no estimate, as in R's `glm`, so that rule is left out with
  a warning, and the other rules are returned with R's values. This is a
  recorded correction of R.
- **`seqdomassoc(with_missing=True)` keeps the missing level for every domain
  (audit TX-AUD-20261011-013).** The missing row or column was added only to a
  domain that had missing cells; R keeps it for every domain, so its empty cells
  enter the degrees-of-freedom adjustment (df 4 where traminex gave 1 on the
  recorded example). It is now kept as in R.
- **`seqdomassoc`'s `seqrep_args` and `seqrf_args` replace the defaults
  (audit TX-AUD-20261011-012).** A given mapping was merged into the defaults
  (`coverage=0.8, pradius=0.1`; `k=20`), so `seqrep_args={"pradius": 0.2}`
  kept coverage 0.8. As in R, a given mapping now replaces them, and `seqrep`'s
  and `seqrf`'s own defaults apply to what it leaves out.
- **`seqdomassoc(wrange=...)` counts co-occurrences within each person (audit
  TX-AUD-20261011-009).** The windowed cross-table copied R's `seqxtabs.win`,
  which windows all cases' states end to end and so pairs states of different
  people, ignores the weights and leaves out the missing state with
  `with_missing=True`. It now counts, for each person, the state at each
  position against the states inside the window in the other domain, weighted
  by that person's weight. On two people with constant sequences `A`/`X` and
  `B`/`Y`, `wrange=(0, 1)` gives `A-X 5, B-Y 5` where R gives `3, 2, 2, 3`.
  This is a recorded correction of R.
- **Index plots draw the first sequence at the bottom (audit
  TX-AUD-20261010-059).** R's index plot is a horizontal bar plot whose first
  bar, the first sequence, is the lowest; `sequence_index_plot` (also
  `seqplot(type="i")`) and `sequence_index_plot_grouped` drew the first
  sequence at the top. They now draw it at the bottom, as R does.
- **`seqpcplot` places simultaneous events together and draws state sequences
  (audits TX-AUD-20261010-045, TX-AUD-20261011-027).** Events were placed at
  1, 2, 3, ... in each sequence, so two events at the same time stood one step
  apart, and a state sequence object was refused. As in R, an event's position
  is now the rank of its distinct time within its sequence (simultaneous events
  share it), and a state sequence object is drawn through its distinct
  successive states, also by `seqplot(type="pc")`.
- **`plot_subseqelistchisq` colours bars as R does (audit
  TX-AUD-20261010-045).** The bars were coloured on a continuous red-blue scale
  with positive residuals red, and the panels titled "Group g1". As in R, each
  bar now takes the colour of its residual's bin, cut at the 5% and 1% normal
  quantiles (Bonferroni-adjusted when `seqecmpgroup` used it), from the five
  Brewer "RdBu" colours, so a positive residual is blue, and each panel is
  titled by its group label.
- **`seqmtplot` draws R's bars and no error bars by default (audits
  TX-AUD-20261010-043, -058).** R's `seqmtplot` draws no error bars unless asked
  (`seqmeant`'s `serr` defaults to `FALSE`); traminex drew 95% whiskers by
  default. The default is now `errbar=False`. The bars are R's `seqmeant` rows:
  the states, plus the missing state with `include_missing=True`; ragged data
  no longer gets zero bars for the missing and void markers (`A, B, *, %` became
  `A, B`).
- **Survival curves run to the last follow-up time (audits TX-AUD-20261010-042,
  TX-AUD-20261011-026).** `seqsurvplot`, `seqsplot` and `seqedplot(type=
  "survival")` drew each curve only to its last event, and a state whose spells
  are all censored had no curve. As R's `plot.survfit` draws them, each curve now
  runs on to its stratum's last follow-up time, censored or not (mvad 1:150:
  employment to 72 rather than 48), a state without events is a flat line at
  1, and `seqedplot` stops at the right end of its `ages` window. `SeqSurv` has
  a new field, `follow_up`, holding each stratum's last time.
- **`seqplot_rf` prints R's statistic (audit TX-AUD-20261010-041).**
  TraMineRextras' `seqplot.rf` prints, and writes into its subtitle, its own
  quality statistic from squared distances to the group medoids and to the
  global medoid; `seqplot_rf` printed `seqrf`'s instead (for the same groups,
  `0.2` where R prints `0.2307692`). It now prints R's statistic with R's seven
  significant digits and uses it in the quality line (`NaN` as R spells it when
  there is one group or one sequence per group).
- **`dissmfacw` uses the rank of a rank-deficient design (audit
  TX-AUD-20261010-040).** When the predictors are collinear (five region dummies
  that sum to one, a numeric copy of another predictor, a constant), R's
  `dissmfacw` projects onto one arbitrary extra direction and counts degrees of
  freedom in columns; traminex copied that with a different arbitrary direction,
  so the two disagreed and both were wrong. The projection is now onto the
  column space and the degrees of freedom are ranks: the model and every
  estimable term give the values of the same design without its redundant
  columns, and a term that adds nothing given the others has no test (NaN, with
  a warning). The ledger also records that terms are named predictors, without a
  formula interface or interaction terms.
- **`seqplot` draws one panel per group level (audit TX-AUD-20261010-044).**
  `seqplot(group=...)` raised `TypeError` for every type but `Ht`. It now draws
  one panel per level of `factor(group)`, each the plot of that level's
  sequences (with `diss` and a `sortv` vector taken per panel), titled
  as R does: by the level, by `"<main> - <level>"` for one `main`, or by `main`
  itself when it gives one title per level. As in R, types `"r"` and `"rf"`
  first compute on the whole set what the panels share: the distances (from
  the distance keywords when `diss` is not given) and `dmax` for `"r"`, and the
  default MDS sort variable for `"rf"`; computed per panel, they gave other
  representatives and medoids than R's grouped plot. Type `"r"` refuses a level
  with fewer than two sequences, as R does. `seqfplot`, `seqmtplot`,
  `seqmsplot`, `seqdHplot`, `seqrplot` and `seqrfplot` take `group=` the same
  way. As for `seqHtplot`, numeric groups keep `factor()`'s numeric order and
  the titles have no leading space (R's internal `group()` adds one).
- **`seqrep` honours `random_state` and refuses keywords it cannot use (audit
  TX-AUD-20261010-015).** With a dissimilarity matrix given, `seqrep` dropped
  every extra keyword, including `random_state`, so `criterion="random"`
  followed the global NumPy stream and a misspelled option passed silently, also
  through `seqrplot`. `seqrep` (and `seqrep_grp`) now take `random_state`; with
  `diss` given, a keyword meant for `seqdist` raises `TypeError` (R drops it).
- **Tree displays show the node text (audit TX-AUD-20261010-014).**
  `disstreedisplay`, `seqtreedisplay` and the `*2dot` writers drew empty node
  boxes unless a node-image function was given. R always draws node images and
  writes the node text into their titles: the branch label (`Root`, `<= 3.1`,
  `[ b,c ]`), `n: .. s2: ..` and `Split: .. R2:..`. When no node image is drawn,
  that text is now the box's label; with images at the leaves only, the inner
  boxes stay blank, as in R.
- **`seqrecode` keeps TraMineR's colours and takes vector labels (audit
  TX-AUD-20261011-024).** States kept as they were, and the `otherwise` state,
  were coloured grey (`#808080`); as in R a retained state now keeps its colour,
  a recoded state takes the colour of its first source state and the `otherwise`
  state that of the first original state no recode uses; a recode to the missing
  or void code no longer shifts the other colours. `labels` and `cpal` also
  accept R's plain vectors (one entry per new state); a list raised
  `AttributeError`.
- **`seqpm` matches POSIX character classes as R does (audit
  TX-AUD-20261010-054).** R's `seqpm` uses R's regular expressions, where
  bracket classes such as `[[:upper:]]` work; Python's `re` read them as plain
  character sets, so `[[:upper:]]B` silently matched nothing. The classes are
  now translated (Unicode classes, as R's under C.UTF-8, so `[[:alpha:]]` also
  matches `é`), and a backslash inside a bracket is an ordinary character, as
  in R (`[\]A]` is the bracket `[\]` followed by `A]`); the rest of a pattern
  is still read by Python's `re`.
- **`seqivardur`'s result can be copied and pickled (audit
  TX-AUD-20261010-053).** `copy.copy` and `pickle` raised `RecursionError`, and
  `copy.deepcopy` returned a bare array without `vmax` and `meand`, because the
  result forwarded every unknown attribute, the copy protocol's included, to its
  value array. Private and special names are no longer forwarded.
- **`seqindic_dyn` checks `indic` as TraMineRextras does (audit
  TX-AUD-20261010-051).** An indicator group name (`"basic"`, `"all"`,
  `"complexity"`, ...) silently gave the first indicator's values, a list raised
  `TypeError`, and a name `seqindic` does not compute raised `IndexError`. As
  in R, a group name now stops with "Bad indic value, group name not
  supported!", a list uses its first value with R's message as a `UserWarning`,
  and an unknown name is refused by name.
- **`seqtab(include_missing=False)` keeps shorter sequences (audit
  TX-AUD-20261010-027).** The option (a traminex extension) dropped every
  sequence with a void cell as if void were missing, so the shorter sequences of
  a ragged object vanished and `ex1` gave an empty table. Only sequences with a
  missing cell are now left out; `seqfplot(include_missing=False)` follows.
- **`seqstatl` reads number-like states as R's `as.numeric` does (audit
  TX-AUD-20261010-050).** It used Python's `float()`, which accepts digit-group
  underscores and refuses hexadecimal: states `1_1, 1_2, 2_1` came back as the
  numbers 11, 12, 21, and `0x10, 2` as text. They now follow R: `1_1` is text,
  `0x10` is 16, and `" 2"`, `1e2` and `Inf` are numbers.
- **Group labels read as R writes them (audit TX-AUD-20261010-026, -049).**
  `seqdiff`'s discrepancy columns, `group_p`'s levels and `rowmode` used Python's
  text form (`False`, `100000`, `1000000000000000`, `0.30000000000000004`) where
  TraMineR uses `as.character` (`FALSE`, `1e+05`, `1e+15`, `0.3`), and `seqdiff`
  kept an unused level of a categorical group as an empty column; it now drops
  it, as R's `factor()` does.
- **A real state spelled like an external marker stays a real state (audit
  TX-AUD-20261010-034).** When an object kept its missing (or void) marker
  outside the alphabet and a real state used the same text (alphabet `A, *` with
  `nr="*"`), the internal canonical form merged the marker into the real state's
  slot, so operations that canonicalize (for example `seqgen_missing(p_cases=0)`)
  turned the real `*` cells into missing ones. The marker now gets its own
  spelling (`<NR>`, `<VOID>`), as `seqdef` gives it.
- **`HSPELL_to_STS` checks its column selectors (audit TX-AUD-20261010-036).**
  A numeric selector outside the table raised a raw `IndexError`, a missing name
  a raw `KeyError`, and an R-style one-based selector that landed on a text
  column raised pandas' "Unable to parse string" without naming it. They now
  stop with R's `undefined columns selected`, or name the column that cannot be
  read as times. The other column names (`id_col`, `pvar`, `endObs`, the
  per-person bounds `tmin_col` and `tmax_col`, and the vertical layout's
  columns) are checked the same way; a misspelled bound was silently ignored
  and the others raised a raw `KeyError`. `docs/API_PARITY.md` now says which
  functions take zero-based and which one-based numeric selectors; it claimed
  zero-based everywhere.
- **Vertical spells give `HSPELL_to_STS`'s horizontal result (audit
  TX-AUD-20261011-019).** The vertical layout (one row per spell, a traminex
  extension) had its own rules: uncovered periods became a real `"_"` state,
  `endObs` meant the last observed period, `pvar` realigned through `seqstart`,
  there were no time labels, a missing end raised a raw `ValueError`, and ids
  were sorted as text (`10, 2, 3`). The spells are now regrouped to one row per
  person (first-appearance order) and converted by the R-faithful horizontal
  route, so both layouts give the same object; `endObs`, `pvar` and per-person
  bounds must hold one value per person, and may have any name (a column named
  like a regrouped column, such as `id` or `b1`, no longer overwrites it).
- **`HSPELL_to_STS` keeps the period before the end of observation
  (TX-REP-20261017-005).** With `endObs`, TraMineRextras blanks every cell from
  time `endObs - 1` on, because it shifts the spell times by `1 - tmin` but
  `endObs` only by `- tmin`: a spell X 2000-2002 with `endObs = 2003` lost 2002.
  traminex reproduced that; it now blanks from `endObs` on, reading `endObs`
  as the first unobserved period, a recorded correction of R
  (`docs/API_PARITY.md`, `HSPELL_to_STS` section).
- **`HSPELL_to_STS` gives the states the alphabet R's `seqdef` gives them (audit
  TX-AUD-20261010-031).** It sorted the states by code point (`1, 10, 11, 12,
  2, 3, 9`; `A, B, C, _x, a, b, c`), whereas `TraMineRextras::HSPELL_to_STS`
  returns a character table that `seqdef` sorts as `seqstatl` does: numbers by
  value (`1, 2, 3, 9, 10, 11, 12`), R's collation otherwise (`_x, a, A, b, B,
  c, C`). Two or more text states now need ICU, as for `seqdef`.
- **`seqdef_from_long` keeps its time axis and accepts `time_grid=` (audit
  TX-AUD-20261010-018).** It pivoted the person-period rows to one column per
  observed time and then dropped those column names, so the object had no time
  labels and a period nobody reports vanished without trace. The time values are
  now the `xtlab` labels (`2001, 2002, 2003, 2005`), and the new keyword
  `time_grid` (as in `HSPELL_to_STS`) places the sequences on a fixed grid
  (`range(2001, 2006)` gives five positions); without it no period is invented,
  so wave-coded times such as 1, 3, 5 keep three positions. A grid entry given
  twice is refused (it repeated the period). Date and date-time labels read as
  R's `as.character` gives them: the date alone at midnight (`2001-01-01`, not
  `2001-01-01 00:00:00`), the clock time otherwise, never the time zone.
- **`STS_to_TSE` keeps the object's ids, and `to_SPS`/`to_DSS` leave void
  cells out (audit TX-AUD-20261010-033).** `STS_to_TSE` numbered the sequences
  1..n even when the object had ids, and sorted its rows by id; it now carries
  the ids in the input order (R numbers the sequences only when no `id` is given,
  and then warns that it replaces the original ids). The string helpers
  `traminex.seqformat.to_SPS` and `to_DSS` printed the void marker as a state
  (`A/1-B/1-%/2`); like R's SPS and DSS forms they now end a shorter sequence
  where its void cells start (`A/1-B/1`), and a missing cell stays a state.
- **SPELL frames work with pandas operations that compare attributes
  (TX-REP-20261017-004).** `STS_to_SPELL` and `seqformat(to="SPELL")` kept the
  per-spell state codes as a NumPy array in `DataFrame.attrs`, so `astype` with
  a mapping and `pd.concat` of two such frames raised "The truth value of an
  array ... is ambiguous". The codes are now a tuple of integers.
- **`STS_to_SPELL` labels the spells with the object's ids (audit
  TX-AUD-20261010-017).** It numbered them 1..n even when `seqdef(id=...)` had
  set ids, while TraMineR (and traminex's own `seqformat(to="SPELL")`) uses the
  sequence object's row names: `id=["p1", "p2"]` now gives `p1, p1, p2, p2`
  instead of `1, 1, 2, 2`. An object without ids keeps the integer ids.
- **The row engine of `seqdist` evaluates many pairs at once, and `refseq`
  computes only the pairs it returns (audit TX-AUD-20261010-022).** OM with
  per-state indel costs or an asymmetric matrix, OMstran, OMslen, OMspell,
  OMloc, TWED, LCS on sequences of unequal length and NMS ran one pure-Python
  dynamic program per pair: on 60 sequences of length 72 all but OMspell were
  30 to 320 times slower than TraMineR, and a `refseq` call computed the whole
  matrix to return one column. Each dynamic program now runs for many pairs at
  once with the same operations in the same order, so results are unchanged
  bit for bit, and those methods are 30 to 60 times faster; the measured times
  are in `docs/DISTANCES.md`. The pairs are generated in bounded blocks, so the
  memory these routes need beyond the distance matrix stays small at any number
  of sequences, and a call that evaluates only a handful of pairs (fewer than 5
  to 9, by method) runs the per-pair reference, which is faster there. NMS
  stays about four to seven times slower than TraMineR.
- **The grouped plots order text groups as TraMineR does (batch-2 review of the
  collation repairs; TX-REP-20261017-002).** `seqHtplot`, `seqedplot` and
  `seqplot_tentrop` sorted text group levels by code point or in pandas' order
  (`10, 9, A>%, A>*, B, _c, ...`) where TraMineR's `factor()` and `seqplot` use
  collation (`_c, 10, 9, a, A>*, A>%, ...`); they now follow R, and a group whose
  name starts with `_` is no longer left out of the legend. Text mixed with numbers
  is refused, as an R vector holds one type. Numeric groups keep `factor()`'s
  order (`1, 2, 10`) in `seqHtplot`, where R's `seqplot` orders them as text
  (`1, 10, 2`) with a space before each title: a recorded correction of R.
- **Event labels and text groups are ordered with TraMineR's collation in the
  event and group-comparison functions (audit TX-AUD-20261011-015, -017).**
  `seqecreate`'s event dictionary,
  `seqecmpgroup`, `group_p`, `seqsammeha` and the group comparisons (`dissCompare`,
  `seqCompare`, `seqLRT`, `seqBIC`) used a portable case-folded order that differs
  from R for punctuation and accents: a transition to missing (`A>*`) sorted after
  one to void (`A>%`), which changed `seqefsub`'s tie order, and `élevé` after
  `faible`, which changed the group `dissCompare` samples first and so its sampled
  statistic. They now use R's ICU collation (ICU becomes required for two or more
  text labels there); the portable key is gone.
- **`sorti`, `sortv`, the index-plot `sortv="from.start"|"from.end"` and
  `seqrf` rank missing and void cells after the states however a sequence
  object stores them (audit TX-AUD-20261011-006, -014).** They ranked cells by
  their integer codes, so an `STSList` whose missing or void code sat inside or
  before the real states (as `events_to_state_sequences(align_missing=True)`
  can produce) sorted differently from the same cells as `seqdef` stores them,
  and `seqrf`'s groups and R² changed. Cells now rank as TraMineR's factor
  levels: the states, then missing, then void (`sorti`/`sortv` rank void as
  missing, as TraMineRextras does).
- **A text `sortv` orders index-plot rows as R's `order()` does (audit
  TX-AUD-20261010-060).** It was sorted by code point (`['b', 'A', 'a', 'B',
  '_x', 'c']` gave rows 2, 4, 5, 3, 1, 6 where TraMineR gives 5, 3, 2, 1, 4, 6),
  and a missing value raised `TypeError`. Text now follows R's collation (ICU)
  with missing values last; a missing value (`None` or `NaN`) among numbers or
  logicals also sorts last, after `Inf`, as R's `order()` puts `NA`; text mixed
  with numbers is refused.
- **Text groups are ordered as R's `factor()` orders them in the discrepancy and
  tree functions (audit TX-AUD-20261010-039).** `dissassoc`, `dissvar_grp`,
  `disscenter`, `disstree`/`seqtree`, `seqrep_grp` and `dissmergegroups` sorted
  text levels by code point (`B, _c, a, é`) where TraMineR uses collation
  (`_c, a, B, é`), so group rows, tree rules, representative groups and merged
  group numbers came in a different order (statistics unchanged). They now use
  R's order, which needs ICU for two or more text levels; numbers, booleans and
  declared pandas categories keep their order.
- **Concurrent `seqefsub` calls no longer stop each other's mining early (audit
  TX-AUD-20261011-021).** The miner stops when a pass adds no node to its prefix
  tree; traminex counted the nodes of every tree in one process-wide counter, as
  TraMineR's single-threaded C++ does, so a mining call in another thread could
  end a pass and drop subsequences (6 of 7 in a reproduced interleaving). Each
  tree now counts its own nodes.
- **Subsetting a `seqefsub` result gives a subsequence list, as in TraMineR
  (audit TX-AUD-20261010-038).** `fs[0:2]` returned the label strings, so
  `seqecmpgroup(fs[0:2], group)` raised "subseq should be a subseqelist", and a
  list of positions or a boolean mask raised a raw `TypeError`. A slice, a
  sequence of positions or a mask now returns a list of the same class with its
  statistics rows (a `seqecmpgroup` result keeps its group labels); an integer
  still returns the label.
- **`seqedist` gives event sequences without events their defined distance and
  records both departures from TraMineR (audit TX-AUD-20261010-037).** TraMineR
  returns 0 between an empty event sequence and every later one, so its matrix
  depends on the row order; traminex already returned the cost of inserting the
  other sequence's events, without recording the difference (the ledger row said
  EXACT). For a raw distance of 0 whose bound is 0 too (two empty sequences, or
  every `idcost` 0) the default `norm="YujianBo"` returned NaN, as TraMineR
  does; it now returns 0, as TraMineR's own distance normaliser does for a zero
  distance. Both are recorded as divergences in the semantic ledger.
- **Calibrated distance profiles fail loudly and share one precedence rule
  (audit TX-AUD-20261010-032).** `trate_profile="calibrated"` without a
  loadable calibration file silently used the default TRATE settings, and
  `dhd_profile="calibrated"` silently switched to the legacy DHD with
  `variant="mean"`; a file without the profile's section or with a misspelt
  key was ignored the same way. These now raise a `ValueError`. An explicitly
  passed value now overrides the file's for both profiles (the TRATE profile
  let the file win, the DHD profile the argument).
- **The low-level `pairwise(method="CHI2"|"EUCLID", weights=)` checks the
  weights and EUCLID ignores them, as TraMineR does (audit TX-AUD-20261011-008,
  TX-AUD-20261010-067).** Only the length was checked: NaN or infinite weights
  gave an all-zero CHI2 matrix and a negative weight altered it, where `seqdist`
  refuses such weights. They now get `seqdist`'s case-weight contract (finite,
  non-negative; all-zero weights still give TraMineR's all-zero CHI2 matrix).
  EUCLID read the weights to decide which states are present, so a zero weight
  changed its distances; TraMineR never weights EUCLID, and neither does traminex.
- **`seqdist` refuses a NaN TWED `nu`/`h`, OMslen `h`, OMspell `expcost` or
  OMloc `expcost`/`context` (audit TX-AUD-20261010-013).** Their sign checks
  were false for NaN, so the call returned a matrix of NaN where TraMineR
  stops. They now raise a `ValueError` naming the parameter, in `seqdist` and
  in the low-level `pairwise`. As in TraMineR, `tpow=NaN` and `nu=Inf` still
  compute (NaN and infinite distances).
- **The low-level `pairwise(method="DHD")` checks explicit costs and applies
  `norm` as `seqdist` does (audit TX-AUD-20261010-012, -023).** With explicit
  costs it accepted negative, NaN and non-zero-diagonal values (a negative
  matrix gave negative distances) and raised a raw `IndexError` for an
  undersized array, where TraMineR refuses all four; it divided by the
  sequence length for `norm="maxlength"` and refused the other norms. Without
  explicit costs it ignored `norm` altogether. Both routes now use TraMineR's
  bound, the sum over positions of the largest cost, for every norm. The
  legacy marginal-frequency DHD (`variant`, `smooth_window`, calibrated
  profile), which has no normalisation, now refuses a `norm` other than
  `None`/`"none"` instead of ignoring it, in `pairwise` and in `seqdist`.
- **`seqdist(method="TWED", indel="auto")` resolves TWED's automatic indel
  (audit TX-AUD-20261010-011).** `"auto"` is TraMineR's default `indel`; for
  TWED R resolves it to `2 * max(sm) + nu + h` over the whole substitution
  matrix. traminex refused the string ("indel='auto' is not a number", or
  "could not convert string to float" with an explicit matrix). It now gives
  R's distances, as an omitted `indel` already did; the low-level `pairwise`
  still takes only numbers.
- **`seqdef_from_long` and `TSE_to_STS` read `labels=` and `cpal=` in the same
  order as `seqdef` (repair finding TX-REP-20261012-001).** Both document
  `seqdef`'s alphabet order but paired supplied labels and colours with their
  first-appearance alphabet, so after this release's `seqdef` correction the same
  data in long form gave the states different labels than in wide form.
- **`seqplotMD(type="i", sortv="from.start"|"from.end")` keeps each person on
  the same row in every domain (audit TX-AUD-20261010-008).** Each panel
  sorted its own domain, so the panels showed different people in a row. As in
  TraMineR, the multidomain sequences are now sorted once and every panel
  follows that order.
- **`seqdiffplot` plots `seqdiff`'s default series (audit
  TX-AUD-20261010-007).** It declared `norm=None` and always passed it on, so
  by default it drew the discrepancy of unnormalised LCS distances, four times
  `seqdiff`'s and TraMineR's default (`norm="auto"`), and refused
  `seqdist_args=` unless `norm="auto"` was also given. Its default is now
  `seqdiff`'s.
- **`disstree_assign` reads rules for covariates whose names contain dots,
  spaces or R backquotes (audit TX-AUD-20261010-006).** The rule parser took
  the covariate name as a run of word characters, so rules for R-style names
  such as mvad's `N.Eastern` (including the rules `disstree_get_rules` writes
  for them) raised "cannot parse rule term". A name now runs up to the
  operator. For a non-syntactic name such as `my var`, R's own
  `disstree.assign` cannot read the rules its `disstree.get.rules` writes;
  traminex assigns them (recorded as a divergence in the semantic ledger).
- **`DHDConfig()` and `default_configs()["DHD"]` select TraMineR's DHD (audit
  TX-AUD-20261010-005).** `DHDConfig.variant` defaulted to `"mean"` and was
  always passed on, which switched `seqdist(method="DHD")` to traminex's legacy
  marginal-frequency metric (row 1 of a four-sequence example: 0, 0.25, 0.525,
  0.575 instead of TraMineR's 0, 5.17, 14.67, 13.33). `variant` now defaults to
  `None` and is passed only when set; `"mean"`/`"sum"` still select the legacy
  metric explicitly. As a consequence, `DHDConfig(profile="calibrated")` now
  uses the variant its calibration file names instead of forcing `"mean"`
  (repair finding TX-REP-20261012-002).
- **`seqidegrad(pow=False)` computes R's index without spell integration
  (audit TX-AUD-20261010-004).** R's public `seqidegrad` reads a logical `pow`
  as its spell-integration switch; traminex used `False` only as the number 0
  and kept spell integration on, so `pow=False` gave different values (RANK:
  0.3, -0.3, 0, 0.2 instead of R's 1, -1, 0, 1). `spell_integr` now defaults to
  `None`, meaning R's rule; an explicit value still applies. `seqindic`'s
  degradation indicator uses the same rule.
- **`seqdist` accepts `refseq` for CHI2 and OMstran on a weighted sequence
  object (audit TX-AUD-20261010-003).** Two index sets, an index and `0` (most
  frequent sequence) raised "weights= conflicts with the weights the sequence
  object carries" for both `weighted=True` and `weighted=False`, because the
  weights the call had taken from the object were passed back into its own
  full-matrix computation. The results now equal TraMineR's.
- **`seqecreate` labels numeric and logical events as TraMineR does and
  honours a numeric `end_event` (audit TX-AUD-20261010-002, -016).** Event
  labels are R's `as.character` of the values (`1e+05`, `TRUE`, `1` for the
  double `1.0`) instead of Python's `str` (`100000.0`, `True`), and
  `end_event` is matched after the same rendering: `end_event=99` was silently
  ignored, leaving the end event in the sequences and every length at -1.
- **`seqdef` gives every state the name, label and colour TraMineR gives it
  (audit TX-AUD-20261010-001).** Without `alphabet=`, R pairs `states=`,
  `labels=`, `cpal=` and its default palette with the `seqstatl` order of the
  observed states (numeric order when every state is an R number, R's
  collation otherwise). traminex paired them by position with its
  first-appearance alphabet, so `seqdef(biofam[, 10:25], labels=..., states=...)`
  silently decoded the first case as `...-L-M-...` where R decodes
  `...-LM-LMC-...`, and default colours went to different states than in R.
  The stored alphabet keeps its documented first-appearance order. A state
  spelled as a non-canonical number (`01`, `0x10`) is kept and ordered by its
  value; R turns such cells into `NA` (recorded in the seqdef ledger row).
  As in R, a non-numeric `missing=` code present in the data makes the order
  alphabetical even when the states are numbers. With two or more non-numeric
  states and no `alphabet=`, `seqdef` now needs ICU, like the other
  R-collation routes (`docs/DEV_BUILD.md`). A `states=` vector of the wrong
  length is still refused, now with a message that gives both lengths.
- **Verified remediation (2026-10-08):** validate each dense between-set input
  before dtype promotion; preserve real states named `*` or `%` in costs,
  missing detection, DSS and plot styling; compare scalar LLCS/LLCP missing
  states independently of storage codes. Refuse malformed condensed vectors,
  symmetrize finite extreme and subnormal distances safely, and handle empty
  OM references and normalized Hamming corpora. Raw TRATE uses its documented
  formula without transitions; time-varying raw costs retain the empty time
  axis and validate options. Native LCP/RLCP promotes arithmetic before
  multiplication. Generated-artifact diagnostics survive `--keep-temp`, docs
  wrap long prose/references, and the CI workflow is configured to prepare the
  exact vendored R oracle (never run: no git remote).
  Rendered navigation retains its accessible button role; code blocks support
  keyboard scrolling, and linked code has sufficient color contrast.
  The sequence-tree guide renders its API and feature lists as lists, and
  per-state semantic notes describe markers at their actual code positions.
  Tree rendering resolves callback images beside the DOT file; safe adjacent
  image aliases prevent Graphviz from producing invalid SVG image references
  for filenames containing XML punctuation.

- Preserve explicit real-state identities through modes, recoding, granularity,
  polyads and group comparison; retain modal time metadata and count real dash
  states separately from alignment gaps. Widen transition indexing in heatmaps
  and exploratory tree features. Sum weighted distance pairs once to prevent
  finite overflow, retaining the R kernel's triangle order. Normalize numeric
  person IDs consistently in vertical spell bounds and alignment. Large
  subsequence counts return float64 infinity, with raw/normalized turbulence
  retaining IEEE results; the semantic ledger documents the exact-integer
  correction to R's repeated-state overflow failure. Native LCS conversion and
  inclusive OM counters now avoid signed integer overflow.

- **Remaining-issue repair (2026-10-08):** preserve raw Gower distances
  with extreme or subnormal feature weights and nearby large values;
  exclude missing values across feature types and reject invalid numeric
  and asymmetric binary observations. Raw sequence-cost builders enforce
  the integer state-code contract, avoid narrow-integer transition
  overflow, and accept valid uint64 frequency inputs on NumPy 1.26.
  The corrected raw Gower calculation uses more memory and measured
  3.38 times slower on the recorded 64-row, eight-feature workload.

- Numeric sequence conversion, integration and ranked indices use real-state ranks
  independently of missing-marker position; main tokens include external missing
  codes. Normalized multidomain and weighted-Hamming calculations handle large
  finite weights, and tree image export creates new output directories first.

- Make exploratory sequence-tree features and predictions independent of whether
  declared missing/void codes are stored inside or outside the alphabet.

- **Remediation (2026-10-08):** refuse oversized state codes before narrowing on
  ragged and void routes; align EMLT and fitted Markov likelihood with factory
  missing-state codebooks; honor vertical spell overlap controls; exclude
  weighted self-distances before arithmetic; restore group representative
  defaults and forward selection options. Explicit survival layout disabling
  works at the Matplotlib floor, DOT identifiers support filename punctuation,
  zero-smoothing JSD avoids undefined log terms, and calibration accepts a
  basename output. Reference checks report absent and partial expectations
  honestly. Raw Gower rejects invalid feature weights. Numerical examples,
  semantic option boundaries, plotting ownership and rendered links are corrected.

- **`seqcost` no longer downgrades `with_missing` behind `seqdist`'s back
  (2026-09-20; adopted from TraMineR 2.2-14).** 2.2-13 reset `with.missing` to
  `FALSE` whenever the *source* carried no missing cell, silently undoing the
  decision `seqdist.R:112-113` had already taken across both the source and an
  external `refseq`. The undersized `sm` then produced an out-of-bounds read in
  R — 2.2-13 returned `5.09888213631624e-310, 5.09888181284790e-310, 2` on one
  run and `0, 0, 2` on the next for the same input. `seqdist` with a
  missing-bearing `refseq` previously raised `IndexError` on
  `CONSTANT`/`INDELS`/`INDELSLOG` (`TRATE` was already correct) and now matches
  live 2.2-14 exactly. An explicit substitution matrix too small for the
  reference is refused by name instead of reaching the kernel.
- **`seqdist` applies R's `with.missing` reset to plain integer matrices too
  (2026-09-22).** `seqdist.R:111-116` resets `with.missing` to `FALSE` when
  neither the data nor a sequence `refseq` contains a missing value, but
  traminex applied that rule only to sequence objects. Once `seqcost` stopped
  downgrading the flag (above), a plain array with `with_missing=True` had its
  generated costs widened by a missing state it could never reach: OM, OMslen,
  OMspell and OMstran raised `ValueError: … indel vector length must equal the
  active state count` for `INDELS`/`INDELSLOG` (and OMstran for any vector
  indel), and OMslen, OMspell, OMloc, OMstran and TWED silently returned
  different distances. A plain array now gives TraMineR 2.2-14's result for
  every method and cost string, as the sequence object already did. With an
  external `refseq`, a plain array must also use only codes that are real
  states of the reference: a reference whose missing code equalled a state the
  array used was silently scored as missing, where R refuses the pair ("same
  alphabet"). Reported downstream by `pyweightedcluster`; see
  `ISSUE_SEQDIST_WITH_MISSING_RAW_ARRAY_20260921.md`.
- **Validation gaps in the low-level distance API are closed (2026-09-22).**
  `pairwise(ref_sets=...)` validates its indices like `seqdist` (a negative
  one wrapped, a fractional one was truncated); `seqLLCS`/`seqLLCP` refuse a
  multi-row array as they refuse a multi-row object (row 0 was used
  silently); a cost-method string for OMloc/OMslen/OMspell/OMstran/TWED at
  the `pairwise` layer is refused by name (it failed converting `"TRATE"` to
  a float); state codes beyond 32 bits are refused instead of wrapping to 0;
  the dense default DHD refuses a misspelt `norm` instead of ignoring it;
  and legacy DHD smoothing no longer crashes on zero-width input.
- **Three traminex-only distance helpers agree with the R-facing paths
  (2026-09-22).** `seqRLCP` on a ragged object compared the void padding as a
  state; it now returns `seqdist(method="RLCP", with_missing=True)`'s
  distances, and its `norm="length"` divides by the sum of the two sequence
  lengths instead of twice the padded width. `seqsubm`/`seqcost` on a plain
  integer matrix gave a state that never occurs an `INDELS` cost of `1e12`
  (and nudged `INDELSLOG` off `log 2`) where the sequence-object path gives
  R's largest-double clamp; both now agree. `seqdist` on a weighted object
  refuses an explicit sequence-weight keyword (`weights=` for CHI2/OMstran,
  `trate_weights=`), which it used to drop silently in favour of the object's
  own weights; set weights on the object and choose them with `weighted=`.
- **`seqdist` accepts TraMineR's deprecated argument forms (2026-09-22).**
  `method="OMopt"`/`"LCSopt"` map to `OM`/`LCS` and a logical `norm` to
  `"auto"`/`"none"`, each with a `DeprecationWarning`; a `refseq` list of more
  than two index sets uses the first two with a warning; all were refused.
  `refseq=True`, which R refuses, is refused instead of being read as row 1.
- **Three numeric edges of `seqdist` now match TraMineR (2026-09-22).** CHI2/
  EUCLID with `step` beyond the sequence length and `overlap=True` build R's
  three intervals instead of one (the result was all zeros); OM with an
  infinite indel and `norm="maxdist"` gives R's NaN off the diagonal on dense
  input as it already did on ragged input; and OMslen halves a cost before
  weighting it, so an `INDELS` cost at the largest double no longer overflows
  to `inf`. TWED's finite border stand-in is recorded as a deliberate
  divergence (`docs/API_PARITY.md`).
- **DHD's time-first cost array is documented and a TraMineR-layout array is
  named (2026-09-22).** `seqdist(method="DHD", subm=array)` takes `(T, K, K)`
  while `seqcost(..., time_varying=True)["sm"]` returns R's `(K, K, T)`, so
  passing one to the other failed with a shape error -- or, when `T == K`,
  "diagonals must be zero". The refusal now names the layout and the
  conversion, `np.moveaxis(sm, -1, 0)`; the layout itself is unchanged.
- **`seqdist`'s TraMineR norms account for position-dependent OM costs
  (2026-09-22).** With `subm_by_pos`/`indel_by_pos` (a traminex extension) the
  normalization used the plain scalar `indel`/`subm`, so "normalized" OM
  distances reached 2.5; they now use the position maxima, as `docs/OMLOC.md`
  already specified for the low-level `'max'`.
- **The low-level `pairwise` keeps one meaning for `norm="length"` (2026-09-22).**
  As soon as any row contained the void symbol -- and in
  `pairwise_between_sets` whenever a set had a shorter row -- the row engine
  read `'length'` as TraMineR's `maxlength`, so the same equal-length pair got
  LCS `0.667` or `1.333`, LCP/RLCP `1` or `2`, OM with `indel=2` `1.333` or
  `0.667`. The documented low-level meaning (`l1 + l2` for LCS/LCP/RLCP, the
  longer length for OM) now holds on every path; `seqdist` was not affected.
- **`kweights` for NMS, NMSMST and SVRspell follow R (2026-09-22).** A scalar
  crashed with `IndexError` where R recycles it; a negative or NaN weight was
  accepted where R refuses it; and an infinite weight gave 0 instead of R's NaN
  because a NaN squared distance was clamped to 0.
- **Asymmetric costs give identical sequences identical distances in OMloc,
  OMslen, OMspell, OMstran and TWED (2026-09-22).** With asymmetric substitution
  costs the direction of a comparison matters; TraMineR deduplicates first, so
  the first-occurrence order decides it, while these engines compared pairs in
  row order -- two copies of one sequence could get different distances to a
  third (OMspell `2` against `2.5`). They now use R's order and, like OM and
  HAM, warn that the matrix is not symmetric.
- **Caller-supplied costs are checked on every OM-family and TWED path
  (2026-09-22).** TraMineR's `checkcost` refuses a wrong-sized substitution
  matrix, a negative cost and a non-zero diagonal for OM, OMloc, OMslen,
  OMspell, OMstran and TWED; traminex refused them only on the dense OM path,
  so every other method, and OM on ragged data, returned negative distances
  and non-zero distances between identical sequences. A sequence object's
  matrix must now be exactly one row and column per alphabet state (plus the
  missing state under `with_missing`), and a non-positive scalar indel is
  refused on the ragged OM path as it already was on the dense one.
- **`seqcost` and `seqsubm` give the same costs for both marker
  representations (2026-09-22).** They did not canonicalize, so an object with
  declared external `nr`/`void` codes made `INDELS`/`INDELSLOG` fail (broadcast
  `ValueError`, or `IndexError` when time-varying) and `FUTURE` with
  `with_missing=True` return different costs when `void_code` equalled the
  alphabet length.
- **CHI2/EUCLID and NMS/NMSMST/SVRspell are sized on the alphabet
  (2026-09-22).** On dense input the state count came from the codes a corpus
  happened to use, and `global_pdotj="obs"` was read from `seqmeant` rows that
  include the void marker. So an object whose last state is unobserved, a
  ragged object, and a plain integer matrix all refused `"obs"` (the plain
  matrix with `TypeError: x must be an STSList`), a full-alphabet vector or
  `prox` was refused, and one sized to the observed codes -- which R refuses --
  was accepted. All now match TraMineR 2.2-14; a plain matrix's alphabet is
  `0..max`, as for its cost strings.
- **`HSPELL_to_STS` horizontal spells reproduce TraMineRextras (2026-09-22).**
  Rows were keyed by id text, so they came out sorted as text (`1, 10, 11, 12,
  2, ...`) and misaligned with anything matched by position, and repeated ids
  were silently merged; now there is one row per input row in input order and
  a repeated id is refused by name. A spell with a missing end was dropped
  (and the time grid shortened) where R runs it to `tmax`; uncovered cells
  were filled with `"_"` where R leaves them NA (`fill` now defaults to
  `None`; vertical input keeps `"_"`); `endObs` censored one period late and
  wrote the filler, where R blanks from `endObs - 1`; one `fixed_status` was
  applied to every spell where R takes one per spell; integer statuses read
  `"1.0"` whenever any column was float (the same `iterrows` coercion crashed
  vertical input); `overwrite=False` now paints only uncovered cells; and the
  time labels are R's `a<t>`. One deliberate divergence: an open spell with
  `tmin < 1` runs to `tmax`, where R stops it early (`docs/API_PARITY.md`).
- **Index plots draw every state in its own palette colour (2026-09-22).**
  `sequence_index_plot` (and so `seqplot(type="i"/"I")`, `seqplotMD`,
  `plot_SAMM(type="i")` and the representative overlay) and
  `sequence_index_plot_grouped` passed state codes to `imshow` without a value
  range, so matplotlib stretched the palette over the codes on screen. Whenever
  the displayed rows did not span the whole codebook -- a sampled subset, an
  unobserved last state, missing values elsewhere in the object -- states were
  drawn in a neighbour's colour (e.g. `C` in the missing-value grey), out of
  step with the legend. The colour scale is now fixed to the codebook.
- **`seqistatd`, `seqstatd`, `seqmeant` and `seqmodst` give the same result for
  both marker representations (2026-09-22).** They read their per-state axis
  from the codebook without canonicalizing, so an object with declared external
  `nr`/`void` codes got R's unpadded columns while its in-alphabet twin got the
  documented marker-padded ones (`seqistatd` 3 columns against 5, `seqmeant` 3
  rows against 5), and `seqmodst` returned its modal sequence on the external
  alphabet. Cycle 10 reported this (N-14); its 2026-09-20 erratum judged the
  three statistics unaffected because the twin test covered only the one
  configuration — nr only, `include_missing=True` — in which the shapes
  coincide. They now canonicalize at entry, so an **external-marker input
  returns the padded shape** its in-alphabet twin always had; `seqdef`-built
  objects are unchanged.
- **`seqrecode(otherwise=…)` folds marker cells in both marker
  representations.** R's `recodeFunc` starts every factor cell — nr and void
  included — at `otherwise`; the external-marker regime kept `*`/`%` as
  ordinary states because the mapping loop iterated alphabet text.
- **`cluster_summary` emits the same per-state columns in both marker
  representations** (its column set was keyed on `len(x.alphabet)`).
- **`seqtree_fit`'s scikit-learn refusal named `seqtree`**, the R-faithful
  weighted tree, which needs no scikit-learn at all; three refusals also named
  the import root `sklearn` instead of the installable distribution.
- **The triangle-inequality scan no longer leaks
  `RuntimeWarning: invalid value encountered in scalar subtract`** on
  non-finite cost matrices (`Inf - Inf` is a silent NaN in the R original's C
  loop).
- **Every golden generator asserts its oracle version.** 77 of 108 asserted
  nothing, which lets a regeneration during an oracle substitution run against
  the wrong TraMineR and record a provenance it does not have; two gate tests
  now pin both the assertion's presence and that no generator names a
  superseded version.
- **`docs/INDEX.md`'s self-declared authoritative module map** listed 16 fewer
  modules than the package contains.
- **`EventSequence.time` is read-only once validated.** `__setattr__`
  re-validates the whole triple on assignment, so `e.time = [0.0, nan]` is
  refused — but the same NaN written *through* the array (`e.time[1] = nan`)
  was accepted, leaving the object holding a value its own validator rejects.
  The stored vector is already a fresh copy, so freezing it does not affect the
  caller's array; reassignment remains the supported way to change it.
- **`seqdist`'s "Allowed kwargs" line no longer advertises a refused keyword.**
  `DUR`'s `eps` was removed as unused but left in the allowed set so its
  targeted refusal would fire first, with the side effect that an unrelated bad
  keyword reported `Allowed kwargs: ['eps', 'metric']` — offering the parameter
  the next call rejects. Both messages now come from one table
  (`ALWAYS_REFUSED_PAIRWISE_KWARGS`) and cannot disagree.
- **`docs/DEV_BUILD.md` records the stray `_pb_kernels` hazard.** Backend
  selection is import-based, so a forgotten copy or a stale symlink silently
  switches later runs onto the compiled kernels; nothing in the package warns
  about it (cycle-10 recommendation 4).

- **The remaining external-marker entry points are canonicalized
  (2026-09-05 cycle-9 audit remediation; public API count unchanged).** The
  cycle-8 sweep missed nine public surfaces that still crashed, refused, or
  drew wrongly on the declared-external representation; all now fold markers
  at entry: `seqdist` (and an `STSList` `refseq`) — whose state-count bound
  refused external codes the kernels bound by value, breaking the public
  `seqCompare(a, b, with_missing=True)` route on mixed-width missing input —
  `seqMD`/`seqdistmc` channels, `seqsplot`, `seqstart`, `seqasnum`,
  `seqgen_missing` (which also no longer collides the appended nr with an
  external `void_code == len(alphabet)` and now extends `cpal` alongside
  `labels`), and the renderer set `sequence_index_plot`(+grouped/overlay),
  `seqfplot`, `state_legend`(+figure), `seqmtplot`, `seqmsplot`,
  `cluster_chronograms`(+representatives), `seqrplot`/`seqrfplot`/
  `seqplot_rf`. `STSList._canonicalized()` now merges an external marker into
  an existing alphabet slot bearing the marker's text instead of appending a
  duplicate symbol. `seqalign` builds an explicit code→row map (nr last,
  R's `sm` layout) and refuses a missing state under `include_missing=False`
  with a clear `ValueError` where it previously crashed with an opaque
  `IndexError`; `seqalignplot` inherits the repair. `seqsamm`'s `spell`
  metadata now holds `alphabet(seqdata)` — the real states only — matching
  `seqsamm.R`. `TRATEConfig.smoothing` defaults to the R-faithful `0.0`,
  closing the last non-zero TRATE default (the config-dataclass route).
  The two-representation equivalence suite grew from the cycle-8 anchors to
  cover every repaired entry point, including byte-identical PNG rendering
  comparisons. The stub gate's deliberate annotation-drift blindness is now
  documented in `scripts/check_stub_signatures.py`; four informational
  divergences are recorded in `docs/API_PARITY.md`.
- **Both marker representations of the same data are now analytically
  equivalent (2026-09-02 audit remediation; public API count unchanged).**
  `STSList` gained `_canonicalized()` — identity for in-alphabet markers; for
  declared external `nr_code`/`void_code` an equivalent copy with the markers
  folded into the alphabet — and every consumer the cycle-8 audit proved wrong
  on the external representation canonicalizes at entry (`seqient`, `entropy`,
  `transition_heatmap`, `seqtransn`, `seqsubsn`, `seqfpos`/`seqfposend`,
  `seqST`, `seqivolatility`, `seqivardur`, `seqlogp`, `seqdss`,
  `seqprecarity`, `seqinsecurity`, `seqibad`, `seqidegrad`, `seqintegr`,
  `seqipos`, `seqindic`, `seqcta`, `seqsurv`, `seqpolyads`, `seqdomassoc`,
  the `*_proxy` helpers, and the distribution/chronogram/transition/survival
  renderers). Silent divergences (e.g. `seqient` 0 vs 0.579, normalized
  `seqST` `[3, 1]` vs `[1, 0]`) and the `seqdss`-induced crashes of six public
  functions are gone; `tests/test_external_marker_equivalence.py` pins the
  equivalence with live-R anchors.
- **`seqmodst(weighted=False)` no longer returns the weighted modal
  sequence** when weights are attached (live R: `A-A`/0.5, not `B-B`/0.909).
- **`seqdist` option handling matches TraMineR**: explicit `subm` is required
  for OM/OMloc/OMslen/OMspell/OMstran/TWED exactly where R stops with
  `'sm' must be specified` (HAM/DHD/LCS keep R's auto-generated costs);
  OMspell rejects a negative `expcost`; CHI2/EUCLID refuse an odd `step` with
  `overlap=True` and accept only R's `auto`/`none` norm spellings; the
  friendly `norm` aliases are now true aliases (`'length'` = `maxlength`,
  `'max'` = `maxdist`, upper-triangle `maxscost` bound) — previously
  `'length'` silently halved LCS/LCP/RLCP distances; TRATE smoothing defaults
  are the R-faithful `0.0` at every entry point (previously `1e-9` on the
  low-level paths, which disagreed with `seqdist` by ~1e-10);
  `TRAMINEX_BINDING` is validated at import time like `TRAMINEX_NUMBA`.
- **Tree/diss helpers align with R**: `disstree_get_rules` returns rules in
  R's factor-level (parent-expansion) order instead of a lexicographic sort;
  `_fit_statistics` warns on degenerate df like `dissrf`; `_stars` matches
  R's `cut()` domain; DOT quality lines keep R's separator space;
  `disstree2dot`/`disstree2dotp`/`seqtree2dot` default `use_title=True` like
  R; `seqmtplot` whiskers use the weighted `nbseq`; `seqcta` validates
  `subseq` against the sequence length with a clear error.
- **Event layer**: fractional event times reach the state grid
  (`floor(min)`..`ceil(max)`); `SeqElist` flat views exclude `end_event`
  rows and group subsets rebuild their flat views, so the duck-typed
  counting helpers no longer miscount; `EventSequence` re-validates on every
  field assignment; the three `support_fraction` denominators are documented.
- **Gates**: `scripts/check_stub_signatures.py` now compares literal default
  values, which caught and fixed the `seqedist` stub drift plus stale
  smoothing/`use_title` defaults in five more stubs.
- **Cleanups**: the dead `traminex_core` OM/LCS/mean-Hamming kernels and the
  never-dispatched `pairwise_ham` binding are removed (live `traminex_accel`
  kernels untouched); one shared `_state_color` replaces five copies;
  `ensure_cmap` warns before falling back over a short `cpal`; `seqrfplot`
  accepts a fitted `rf=` instead of smuggling results through a figure
  attribute; `dissrf` is routed from `repseq_rf`; duplicated DHD calibration
  resolution and the dead `_cmdscale`/overwritten `seqrep` score paths are
  gone. Deliberate divergences are documented in place and in the parity
  ledger (`disscenter(medoids_index=True)`, `disstree_assign` first-match,
  GVPATH on all platforms, `seqplot` `i`/`I`, the load-bearing 1e-12 medoid
  near-tie tolerance).

- **Reference-set cross-blocks no longer re-estimate dataset-fitted costs
  (2026-08-23; public API count unchanged).** `traminex.distances.pairwise(...,
  ref_sets=(a, b))` rebuilt a corpus from the two reference sets alone and
  delegated to `pairwise_between_sets`, so any method whose parameters are
  estimated from the data was re-fitted on `A ∪ B` — up to 13.5% away from the
  same package's full-matrix result and from TraMineR on the corpus measured.
  `pairwise(ref_sets=...)` now computes the fit on the whole corpus and returns
  the block, matching `seqdist(..., refseq=...)` and R exactly. Given only two
  sets, `pairwise_between_sets` cannot know the corpus to fit from, so it now
  raises `ValueError` for a dataset-fitted configuration (`CHI2`, `EUCLID`,
  `OMSTRAN`, `DHD` without an explicit `subm`, `OM` with `subm="TRATE"`) and
  names the two supported routes. Pair-local configurations — including every
  shape `weightedcluster._seqdist` delegates — are unchanged and keep the
  efficient rectangular path. The public `seqdist(refseq=...)` and
  `seqdist(ref_sets=...)` routes were already correct and are unaffected.
- **`seqLLCS`/`seqLLCP` refuse a multi-sequence object (2026-08-23).** Both
  silently narrowed to row 0 and returned a plausible scalar. `TraMineR`'s
  `seqLLCS.R`/`seqLLCP.R` reach `.C(C_cLCS, as.integer(seq1), as.integer(seq2))`,
  which a multi-row `stslist` cannot satisfy, so live R stops. They now raise
  `ValueError: sequence objects must contain exactly one sequence` — the message
  the already-repaired sibling `seqcomp` uses. Single-row objects and bare integer
  arrays are unchanged. **Behaviour change:** code relying on the row-0 narrowing
  must slice explicitly. `seqmpos` deliberately still accepts a whole object,
  which is exact R parity.
- **R's non-fatal cost diagnostics are emitted (2026-08-23).** Of R `checkcost`'s
  four non-fatal diagnostics traminex emitted only the symmetry warning.
  `checktriangleineq` is now ported, and a caller-supplied cost matrix can raise
  three further warnings with R's wording: a substitution cost greater than twice
  a scalar indel, a substitution triangle-inequality violation, and an
  indel-vector triangle-inequality violation. A triangle-inequality violation
  means the matrix is not a metric, which silently invalidates MDS, medoid
  selection and clustering downstream. R signals the substitution case through
  `message()`, so it carries the distinct
  `traminex.distances.python_kernels.CostDiagnosticMessage` category (a
  `UserWarning` subclass) and `seqdiff` silences exactly that category around its
  internal `seqdist` call, mirroring R. Diagnostics fire only for a
  caller-supplied matrix, never for one generated from `"TRATE"`/`"CONSTANT"`/
  `"INDELS"`/`"INDELSLOG"`, matching `seqdist.R`.
- **A non-finite automatic indel is refused (2026-08-23).** An explicit `Inf`
  substitution cost — even in a cell no observed state can reach — made
  `indel = max(subm)/2` infinite and every distance `inf`, with only numpy
  overflow warnings. `seqdist` now raises. R's whole-matrix `max(sm)/2` scope is
  **retained deliberately** (it is observable: changing only an unobserved cost
  changes R's normalized distances), and the `subm="INDELS"` + external `refseq`
  route, where `[inf, inf, inf]` is exact R parity, is unchanged.
- **Layer-aware rejection messages (2026-08-23).** `seqdist(sm=...)` named only
  `'TRATE'` although `seqdist` accepts all four `seqdist.R` cost strings; it now
  names all four, while the low-level `pairwise` keeps its accurate `'TRATE'`-only
  text. `seqdist` also validates `norm` against its true accepted set
  (`None`, `'none'`, `'length'`, `'max'`, `'maxlength'`, `'gmean'`, `'maxdist'`,
  `'YujianBo'`, `'auto'`, case-insensitively; `'max'`/`'maxdist'` are `OM`-only)
  instead of surfacing the three-value low-level message. No working alias was
  removed.
- **A missing optional dependency names its extra (2026-08-23).**
  `traminex.seqplot`, `sequence_index_plot`, `state_legend` and `chronogram`
  raised a bare `ImportError: No module named 'matplotlib'`. The lazy facade
  resolver now re-raises naming the extra, matching the exemplary
  `scikit-learn is required for disstree_fit` wording.
- **The R-golden provenance proof is executable again (2026-08-23).** One fixture
  recorded the resolved R `survival` version as *data*, so
  `check_r_golden_provenance.py --execute --passes 2` aborted on any machine whose
  copy differed. A third-party version is now provenance metadata: the generator
  asserts a documented floor (`survival >= 3.8.9`) where a mismatch stops
  generation, and the checker records the resolved version in its evidence JSON.
  No analytic golden value changed.
- **`seqcomp` and `seqmpos` compare states by label (2026-09-22, audit
  S3-01).** R compares the two objects' factor columns, which compares state
  labels; traminex compared the integer codes. Two objects whose alphabets
  differed in content or order were therefore compared code by code:
  `A-B-B` against `B-A-A` was reported identical with 3 matching positions
  (R: `FALSE`, 0), and `A-B-B` coded under the alphabet `["B", "A"]` did not
  match itself. Both functions now decode each object's states before
  comparing, and `seqmpos` excludes each object's void and missing cells by
  that object's own markers. Results for two objects that share one codebook
  are unchanged.
- **The degradation indicators honour all of TraMineR's `method` and
  `penalized` options (2026-09-22, audit S3-02).** `seqidegrad`,
  `seqprecarity`, `seqinsecurity` and the ranked `seqindic` columns share one
  degradation term, which computed the `TRATE` weighting for `method="ONE"`,
  for the `"+"` methods (`"RANK+"`, `"FREQ+"`, `"TRATE+"`, `"TRATEDSS+"`) and
  for any unknown or misspelled method, returned zeros for `penalized=True`
  (R reads it as `"BOTH"`) and for an unknown `penalized`, and accepted any
  `weight_type` and a `border_effect` of one or less. `"ONE"` now gives every
  transition weight one, the `"+"` methods multiply the result by the mean
  transition weight per spell, a logical `penalized` is read as in R, and the
  values R refuses raise `ValueError`. `penalized="NO"` (or `False`) still
  returns zeros from `seqidegrad`, where R's own `seqidegrad` fails; this is
  recorded as a divergence.
- **The `TRATE` weighting of the degradation indicators no longer bridges
  missing cells (2026-09-22, audit S3-03).** With `method="TRATE"` and
  missing values excluded, `seqidegrad`, `seqprecarity` and `seqinsecurity`
  dropped the missing cells before pairing states, so `A * A B` counted an
  `A -> A` and an `A -> B` transition. TraMineR weighs transitions with
  `seqtrate`, which pairs neighbouring positions only, so a missing cell
  breaks the transition. The rates now match `seqtrate` (and TraMineR): on a
  three-sequence example `seqidegrad(method="TRATE")` gave
  `0.35, -0.672, 0.156` where TraMineR gives `0.30, -0.672, 0.10`. Data
  without internal missing cells are unaffected.
- **`seqstatd(include_missing=True)` normalises entropy as TraMineR does on
  data without missing values (2026-09-22, audit S3-04).** R always adds the
  missing state to the state list when missing values are included, so the
  normalised entropy is divided by `log(K + 1)` and the palette gains the
  missing colour. traminex counted the missing state only when the object
  contained a missing cell, so on missing-free data it divided by `log(K)`
  (0.6309 where R gives 0.5 on a three-state alphabet, `NaN` where R gives 0
  on a one-state alphabet) and left the missing colour out of `cpal`. Both now
  match R; `seqient` already handled this case. The frequency matrix keeps its
  documented marker-padded shape.
- **`seqalign(include_missing=True)` takes the cost matrix TraMineR takes on
  data without missing values (2026-09-22, audit S3-05).** R sizes the
  substitution-cost matrix with the missing state whenever missing values are
  included, so it accepts the `(K + 1) x (K + 1)` matrix that
  `seqsubm(..., with_missing=True)` builds and refuses a `K x K` one. traminex
  added the missing state only when a cell was missing, so on missing-free
  data it rejected its own `seqsubm(..., with_missing=True)` output ("sm must
  be 3×3 matrix") and accepted the matrix R refuses. It now follows R.
- **`seqtransn(pweight=True, with_missing=True)` returns NaN where TraMineR
  returns `NA` (2026-09-22, audit S3-06).** The transition weights come from
  `seqtrate`, whose rate matrix has no missing state, so TraMineR gives `NA`
  for any sequence with a transition to or from the missing state. traminex
  skipped those transitions and returned a finite count (`0 0 0` where R gives
  `NA 0 NA`), silently treating them as weightless. Such sequences now get
  NaN; sequences without such a transition are unchanged.
- **`seqmodst` no longer invents a modal state where no cell is valid
  (2026-09-22, audit S3-07).** At a position where every cell is excluded
  (missing with `include_missing=False`, or void) or the total weight is
  zero, TraMineR finds no modal state and leaves `NA`. traminex took the
  first state of the alphabet, reporting `A-A-B-B` where R has `A-NA-B-B`,
  and `A-A-A-A` for all-zero weights where R has `NA-NA-NA-NA`. Such a
  position now holds the missing state (a sequence object cannot hold `NA`)
  with frequency NaN, and, as in R, the modal sequence then counts no
  occurrence.
- **`seqmeant(prop=True)` is NaN when a sequence has no valid cell
  (2026-09-22, audit S3-08).** TraMineR averages each sequence's state
  proportions, and a sequence with no valid cell (for example, entirely
  missing with `include_missing=False`) has undefined proportions, so every
  mean is NaN. traminex treated that sequence's proportions as zeros and
  returned finite means (0.111, 0.333, 0.222 where R gives NaN), which also
  contradicted its own `seqistatd(prop=True)`. It now returns NaN as R does;
  the documented marker-padded rows stay zero.
- **`seqstatf` keeps TraMineR's missing-state row and NaN percentages
  (2026-09-22, audit S3-09).** With `include_missing=True`, R's table always
  has a row for the missing state, at zero when no cell is missing; traminex
  added it only when a cell was missing, so it returned 3 rows where R returns
  4. With a zero total weight, R's `Percent` is `0/0 = NaN`; traminex
  reported 0. Both now match R, and the two marker representations of the
  same data give the same table.
- **Degenerate turbulence and entropy cases give TraMineR's NaN (2026-09-22;
  audit S3-10, S3-11).** A sequence whose every cell is missing has no spell
  when missing states are excluded. For it `seqST` returned a turbulence of 0
  where TraMineR returns NaN (type 1, and type 2 with fewer than three states),
  and `seqivardur(type=2)` returned a `vmax` of NaN where TraMineR returns 0
  once there are three or more states. With a one-state alphabet the
  normalised `seqient` and `seqici` returned 0, and `seqST(norm=True)`
  returned the unnormalised value (1) there and when every sequence has
  length 1; TraMineR divides by a zero maximum in these cases and returns
  NaN. All of these now match TraMineR, including the corresponding `seqindic`
  columns; ordinary data are unaffected.
- **`seqauto` builds TraMineRextras' contingency table (2026-09-22; audit
  S3-13).** R's `table()` keeps a state that is seen only next to a missing
  lag (for instance one that occurs only at the first position) as an
  all-zero row or column, and `chisq.test` then returns NaN for both the
  auto-association and the p-value; traminex dropped that state and returned
  finite values. A table with a single row or column now gets `chisq.test`'s
  goodness-of-fit p-value and an infinite Cramer's V, as in R (traminex gave
  NaN and a p-value of 1), and a 1 x 1 table is refused with R's message
  instead of returning a p-value of 1.
- **`seqsurv` ignores zero case weights as TraMineRextras does (2026-09-22;
  audit S3-14).** R's survival summary keeps only the times at which the
  weighted number of ending spells is positive. traminex also emitted a row
  for every time at which only zero-weight spells ended (seven rows where R
  has five in the audit example), and a NaN survival row when such a time came
  after every weighted spell had ended. Those rows are gone; the survival
  estimates and standard errors of the remaining rows are unchanged.
- **`seqimplic` orders groups and adds the missing state as TraMineRextras
  does (2026-09-22; audit S3-15).** Group levels were sorted as text, so
  numeric groups `1, 2, 10` came out as `1, 10, 2`, `a, A, b` as `A, a, b`,
  and a pandas categorical lost its declared order; the first axis of
  `indices` was permuted accordingly. Levels now follow R's `factor()`.
  With `with_missing=True` the missing state is always added to the
  alphabet (label "Missing", the missing colour, an all-NaN slice), as in R,
  even when no cell is missing. With `na_rm=False` a missing group value is no
  longer a level of its own, and the positions it enters are NaN, as in R.
- **The `weighted` flag of `seqmodst` and `seqmeant` follows TraMineR
  (2026-09-22; audit S3-16).** TraMineR records an object as unweighted when
  it has no weights, when `weighted=False` is requested, or when every weight
  is 1. `seqmodst(...).weighted_` simply echoed the `weighted` argument (so it
  was True for an object without weights), and `seqmeant(...).attrs["weighted"]`
  was True for all-one weights. Both now agree with R and with `seqstatd`;
  no computed value changes.
- **Four statistics functions check their arguments as TraMineR does
  (2026-09-22; audit S3-17).** `seqST` refuses a `type` other than 1 or 2
  (it silently computed type 1). `seqivolatility` refuses a non-boolean
  `adjust`. `seqintegr(state=...)` accepts a single state of the alphabet,
  plus the missing state only with `include_missing=True`; it used to accept
  the missing or void marker without it and return a proportion above 1, and
  it now returns R's zero column for the missing state when no cell is
  missing instead of raising. `seqipos` accepts R's obsolete
  `index="integration"` as `"integr"` with a `DeprecationWarning` instead of
  refusing it.
- **`seqipos` volatility across an internal gap is recorded as a deliberate
  difference from TraMineR (2026-09-22; audit S3-12).** When the binary
  recoding in `seqipos(index="volatility")` leaves void inside a sequence
  (a missing cell with `include_missing=False`, or a state outside explicit
  `pos_states`/`neg_states`), TraMineR's `seqdss` stops at the gap and drops
  the later spells. traminex already counted every observed cell; this is
  now documented in `docs/API_PARITY.md` and the parity ledger, and pinned by
  a test. `seqindic`'s `Vpos` column carries the same difference. No value
  changed.
- **traminex-only summaries no longer count the void padding as a state
  (2026-09-22; audit S3-19).** On sequences of unequal length, `entropy`
  counted the padding cells as a state (0.64 instead of 0 when a single state
  is observed), `seqtransn_count` and `turbulence_proxy` counted the step into
  the padding as a change of state (which `seqtree_fit` used as a feature),
  and the matrix behind `transition_heatmap` had transitions into and out of
  the padding. All three now ignore it; equal-length data are unaffected.
- **`seqtab` no longer slows down quadratically with the number of distinct
  sequences (2026-09-22; audit S3-18).** Finding each pattern's first
  sequence scanned all sequences once per pattern; it now comes from the same
  sort that finds the patterns. With 16,000 distinct sequences
  `seqtab(idxs=0)` took 4.6 s and now takes 0.19 s; the results are
  unchanged.
- **`seqtrate` counts transitions without a per-sequence Python loop
  (2026-09-22; audit S3-18).** The time-varying form on 2,000 sequences of
  length 20 over 8 states took 7.2 s and now takes 0.04 s; `seqlogp`, whose
  default fits time-varying rates, drops from 7.1 s to 0.17 s. Sums are still
  taken sequence by sequence, so every result is bit-for-bit unchanged.
- **Event sequences with fractional times now match TraMineR (2026-09-22).**
  TraMineR stores an event sequence as a chain of gaps built by inserting the
  events one at a time in the order they were given, and reads every count,
  age and printed gap back from that chain. traminex rebuilt the gaps from the
  sorted times, which rounds differently once times are not whole numbers: an
  event at 2.0 after events at 0.1 and 0.2 has age 1.9999999999999998 in R, so
  `seqeapplysub(method="age")` gave 2 where R gives 1, and `seqefsub`,
  `seqeapplysub` (CWIN, CDIST), `seqeage`, `seqedist` and the printed sequences
  could differ as well. `seqecreate` now records the input order and the
  engine replays it, so all of these match R; results for whole-number times
  are unchanged (audit S5-01).
- **`seqecmpgroup` forms its groups as R's `factor(group)` does
  (2026-09-22).** Every group value was turned into text and the text sorted,
  so numeric groups 2 and 10 came out as `10 2` with their `Freq.*` and
  `Resid.*` columns swapped, a pandas categorical lost its category order,
  `b B a` sorted as `B a b` instead of R's `a b B`, `True`/`False` groups were
  not named `TRUE`/`FALSE`, and a missing value became a group called `"None"`
  — a different test from R's. Numbers now sort by value and are named as R
  prints them, a categorical keeps its order and drops unused categories, and
  a missing value (`None`, `NaN`, `pd.NA`) leaves its sequence out of the
  group tables as R leaves out `NA` (audit S5-02).
- **Slicing a `SeqElist` keeps exactly the events of the selected sequences
  (2026-09-22).** The flat `id`/`time`/`label` views of a slice were rebuilt
  by assuming one run of identical ids per sequence. A sequence without events
  (an id whose only row is its end event) broke that: `el[1:3]` raised
  `IndexError` and `el[0:2]` returned the next sequence's events; ids that
  merge into one sequence (1.2 and 1.7 both become 1) returned the wrong
  event. Each sequence's flat rows are now recorded when the list is built; a
  hand-built list whose rows cannot be matched to its sequences is refused
  with a message instead of being guessed (audit S5-03).
- **Stored event vectors refuse in-place writes (2026-09-22).** Following
  `EventSequence.time` (frozen in cycle 10), `EventSequence.id` and `.label`
  now refuse `append`, item assignment and the other in-place list edits,
  which slipped past the validation that assignment applies —
  `ev.label.append("C")` left four labels against three times. The flat
  `id`/`time`/`label` views of a `SeqElist`, which have no setter and must
  describe the same events as its sequences, are read-only too; building a
  `SeqElist` copies the vectors it is given, so the caller's stay writable.
  Reading the vectors, comparing them with lists and assigning a new
  `EventSequence` vector work as before (audit S5-05).
- **`seqecontain` accepts a `seqefsub` result (2026-09-22).** As in R, where
  `seqecontain` starts by replacing a subsequence list with its subsequences,
  each subsequence is now tested for the given events, with the same any/all
  rule; traminex raised `TypeError` (audit S5-06).
- **`seqecmpgroup`'s `pvalue_limit` difference from R is now documented
  (2026-09-22).** When no subsequence passes the limit, traminex returns an
  empty result, while R keeps `1:0` rows and so still returns the top
  subsequence although its p-value fails the limit. The behaviour is
  unchanged; the docstring and the parity ledger now record it as a
  deliberate difference (audit S5-07).
- **Event-sequence gaps are printed as R prints them (2026-09-22).** R renders
  each gap with its `format()`, which keeps seven significant digits and uses
  scientific notation when that is shorter; traminex printed `100000` and
  `0.0001` where R prints `1e+05` and `1e-04`, so `as_character` and the row
  labels of `seqeapplysub` differed. One helper now applies R's rule to
  every gap, including the trailing gap to the end of observation (audit
  S5-08).
- **`seqefsub` warns about count method 6 as R does (2026-09-22).** A
  constraint with `count_method=6` (the age of first occurrence, meant for
  internal use) now raises R's warning "count method 6 is only for internal
  use" before mining; traminex mined silently (audit S5-09).
- **A single number sets the length or weight of a one-sequence list
  (2026-09-22).** As in R, `eseq.lengths = 5` and `eseq.weights = 2.5` now
  work for a `SeqElist` holding one sequence, and are refused with R's
  "should be of the same size" message for a longer one; both setters
  raised `TypeError: len() of unsized object` (audit S5-10).
- **`seqeapplysub` names the missing window size for CWIN (2026-09-22).**
  Counting with the CWIN method, which needs a window size, failed with
  `ValueError: cannot convert float NaN to integer` from inside the counting
  engine when none was set (R aborts the session). It now raises the same
  message as `seqefsub`: "count method 3 or 'CWIN' requires a window.size"
  (audit S5-11).
- **`seqsammeha` names a numeric typology as R does (2026-09-22).** A
  typology given as floats (as a NumPy or pandas column usually is) produced
  indicator columns `SAMM1.0`/`SAMM2.0` and levels `1.0`/`2.0`; R writes the
  values with `as.character`, giving `SAMM1`/`SAMM2` and `1`/`2`. Values are
  now written the same way. The parity ledger also records that R fills a
  numeric typology's indicator columns by position, so a typology such as
  2, 1 swaps them in R and 10, 20 fails there, while traminex matches by
  value (audit S5-12).
- **The parity ledger describes TraMineR 2.2-14's new summary methods
  accurately (2026-09-22).** The `summary.seqelist`, `summary.eseq` and
  `print.eseqsum` rows said the per-sequence transition count was available
  through `seqentrans`, which in fact only augments a `seqefsub` result and
  rejects an event-sequence list. The rows now say that `len()`, the event
  dictionary and `seqelength` give the other components and that the
  transition count has no accessor of its own (audit S5-13).
- **`seqdist`'s cost diagnostics follow TraMineR 2.2-14 (2026-09-22).** The
  warnings for a caller-supplied cost matrix still used TraMineR 2.2-13's
  wording, and a test claimed they had been checked against 2.2-14. They now
  reproduce 2.2-14's `checkcost` word for word: the asymmetry warning reads
  `Non-symmetric substitution cost matrix!`; the "greater than twice the indel"
  and indel triangle-inequality diagnostics are `CostDiagnosticMessage`s, as
  they are `message()`s in R, so `seqdiff` now silences them as R does; two
  missing diagnostics were added, for zero costs between distinct states and
  for a largest cost of 0 or below (an all-zero matrix used to produce a
  spurious triangle-inequality message instead); states are named by their
  labels and numbers printed as R prints them; and a DHD time-varying cost
  array now gets the same diagnostics, one asymmetry warning per asymmetric
  position. No distance value changes. Code that matched the old asymmetry
  text must match the new one. Audit S10-01, S1a-11.
- **Distance documentation now matches behaviour (2026-09-22).** The
  `seqdist` docstring and `docs/API_CONSISTENCY.md` said the `norm="max"`
  alias was OM-only; `seqdist` accepts `"max"` (as `"maxdist"`) and `"length"`
  (as `"maxlength"`) for every method that takes a TraMineR norm. The
  documentation of the low-level `norm="length"` (the default of
  `om_trate_pairwise`) now says it omits TraMineR's indel factor, and
  `docs/API_PARITY.md` no longer claims that every explicit non-finite cost
  matrix is refused: a non-finite cost that no observed state can reach is
  computed when the indel is given. No behaviour changes. Audit S1a-13.
- **The Numba backend's agreement is documented as measured (2026-09-22).**
  `docs/DISTANCES.md` said only `JSD_TRANS` could differ between the Numba and
  reference kernels; `DUR` with `metric="chi2"` also differs in the last bits
  (at most 4.4e-16 measured, seen with 9 states). The text now gives both
  bounds and names `TRAMINEX_NUMBA=off` for bit-for-bit reproduction, and the
  backend-agreement test covers a 9-state `DUR` corpus. No result changes.
  Audit S1b-14.
- **`seqtab(idxs=[0])` selects every rank, as TraMineR's `idxs = c(0)` does
  (2026-09-22; audit S3-21).** R has no scalars, so a one-element vector is
  the same as the number: `0` means every distinct sequence. traminex
  accepted `idxs=0` but refused `[0]`, `(0,)` and `np.array([0])`; a
  one-element `idxs` now behaves exactly like the scalar. Vectors of two or
  more ranks still refuse any rank below 1, as in R.
- **`seqimplic` labels numeric and logical groups as TraMineRextras does
  (2026-09-22; audit S3-20).** Group levels were Python's text for the
  values: `2.0`, `10.0`, `100000.0`, `False`/`True`. They are now R's `2`,
  `10`, `1e+05`, `FALSE`/`TRUE`, and two values that R prints identically
  (such as `0.1 + 0.2` and `0.3`) share one level, as in R.
- **`seqsurv` orders and labels groups as TraMineRextras does (2026-09-22;
  audit S3-22).** Text group levels were sorted by character code (`A`, `a`,
  `b`), so the per-group curves, their labels and the default colours came
  in a different order from R's (`a`, `A`, `b`). They now follow R's
  collation, as `seqimplic`, `seqdiff` and `rowmode` already do, and numeric
  and logical groups are labelled `2`, `10`, `TRUE` instead of `2.0`, `10.0`,
  `True`. `seqsplot` and `seqsurvplot`, which build their groups the same
  way, follow suit.
- **`dissassoc(weight_permutation="diss")` now permutes the Levene test
  (2026-09-22, audit S4-01).** The permuted Levene statistic was computed from
  the unpermuted group members, so it never changed and its p-value was always
  1. TraMineR takes the residuals of the permuted members
  (`dissassocweighted.R:119-125`), and traminex now does the same. The Levene
  p-value of `dissassoc` in this mode, and of the final adjustment of `disstree`
  and `seqtree` fitted with `weight_permutation="diss"`, changes accordingly;
  the observed statistics and the other four p-values are unchanged.
- **`disstree` handles categorical and boolean predictors as TraMineR does
  (2026-09-22, audit S4-02).** An ordered pandas categorical with text
  categories raised `ValueError: could not convert string to float`; one with
  numeric categories was split in numeric order instead of the declared order;
  an unordered categorical listed its levels alphabetically; and a boolean
  column produced a threshold rule (`b <= 0`). As in R (`disstree.R:404`,
  `:543`), a categorical now keeps its declared category order, an ordered
  categorical splits between consecutive categories and reports level labels,
  only a numeric column reports a threshold, and a boolean column is a factor
  with levels `FALSE`/`TRUE` (rule `b %in% c('FALSE')`). `disstree_assign`
  compares profile values with rule labels the way R's `%in%` does, so a
  boolean profile matches those rules.
- **Permutation p-values are NaN where TraMineR's are NA (2026-09-22, audits
  S4-03 and S4-07).** A statistic that could not be computed (for example the
  Bartlett statistic of `dissassoc` when a group has a single case) received a
  p-value of 0, which reads as maximally significant; a NaN among the permuted
  statistics lowered the p-value; and `dissmfacw(nperm=1)` reported p-values of
  1. As in R's `TraMineR.permutation`, a p-value is now NaN when its statistic
  is NaN in the observed data or in any permutation, and every p-value is NaN
  when there is no permutation (`R <= 1`, `nperm=1`), in `dissassoc` (all
  weighting modes, and so in `disstree`/`seqtree` adjustments) and in
  `dissmfacw`. `dissmfacw` also no longer raises `ZeroDivisionError` when all
  dissimilarities are zero (identical sequences): like R it returns NaN
  statistics.
- **`disstree` puts the same child on the left as TraMineR (2026-09-22, audit
  S4-04).** A level subset and its complement often explain exactly the same
  discrepancy, and then rounding in the last bit decides which becomes the left
  child, and so the node numbers, the fitted leaf ids and the order of the
  rules. traminex summed differently from R and could put the other child on
  the left. Candidate splits are now evaluated with R's arithmetic (R's `sum()`
  accumulates in long double; its C routines add in case order), so the tree
  matches R bit for bit where NumPy's long double matches R's, as on Linux
  x86-64. The partitions and statistics were already correct. The shared
  phase-2 test inputs are now stored with 17 significant digits, so that the
  tests read exactly the numbers R's goldens were computed from; no golden value
  changed.
- **`disstree` rules print numeric thresholds as TraMineR does (2026-09-22,
  audit S4-05).** A threshold of 123456 printed as `inc <= 1.2346e+05` instead of
  R's `inc <= 123456`, and since `disstree_assign` evaluates the printed rule, a
  case at 123457 was assigned to the wrong leaf. Thresholds now use R's
  `format(x, digits = 5)`, which keeps fixed notation unless it is wider than
  scientific notation, in `disstree_get_rules` and in leaf labels.
- **`dissassoc` refuses `weights` with `weight_permutation="none"` (2026-09-22,
  audit S4-06).** The combination returned statistics that mixed weighted sums
  of squares with unweighted group sizes (a negative Bartlett statistic, for
  example). TraMineR has no test for it and stops; traminex now raises
  `ValueError` naming the weighted alternatives. Without `weights`, `"none"` is
  the unweighted test as before.
- **`dissrep` on identical sequences follows TraMineR (2026-09-22, audit
  S4-08).** When every dissimilarity is zero, traminex reported no similar
  sequences (`nb = 0`) and quality statistics of 0. As in R, `nb` is now the
  total weight (every sequence is identical to the representative), and the
  `Q` statistics and the overall quality are NaN, because there is no
  discrepancy to explain. The same plain division gives NaN or `-Inf` for any
  representative whose distance-to-centre total is zero.
- **`seqrep_grp` leaves out sequences whose group is missing (2026-09-22, audit
  S4-10).** A missing group value (`None`, `NaN`, `pd.NA`) raised a `TypeError`
  or the misleading "Each group must have 2 or more cases". As in
  TraMineRextras, such a sequence now belongs to no group.
- **`dissassoc`, `disscenter` and `dissvar_grp` keep a categorical's level
  order (2026-09-22, audit S4-11).** A pandas categorical group with categories
  `low < mid < high` produced rows, columns or keys in alphabetical order
  (`high, low, mid`). As in R's `factor`, and as `seqdiff` already did, the
  declared category order is now kept and unused categories are dropped. The
  values themselves were already correct.
- **`dissmergegroups` no longer ignores an undefined quality (2026-09-22,
  audit S4-12).** When a group's weighted size is below 1 the average
  silhouette width cannot be computed, and traminex silently returned the input
  grouping with quality NaN. As in R, merging is now refused with a
  `ValueError`; when no merge is attempted (`min_group` already reached) the
  grouping is returned with quality NaN and a `RuntimeWarning`, as
  WeightedCluster warns.
- **`disstree` explains why it refuses two of R's weight permutations
  (2026-09-22, audit S4-13).** With weights, `weight_permutation="permute"` and
  `"random-sampling"` were refused with a message claiming that R accepts only
  the five other values. R does run them, but its split test then replicates
  cases by their weights truncated to integers (a weight of 0.5 counts as 0).
  traminex still refuses them, now saying so, and the parity ledger records
  the difference for `disstree` and `seqtree`.
- **`cluster_summary` and `cluster_representatives` keep the cluster labels
  (2026-09-22, audit S4-14).** Labels were converted with `int()`, so text
  labels raised `ValueError` and fractional labels were truncated (clusters 1.5
  and 1.7 were both reported as `1`). Labels are now returned as given;
  integer labels are unchanged.
- **`disscenter` documents how it treats missing group values (2026-09-22,
  audit S4-15).** TraMineR stops when the group holds a missing value; traminex
  leaves such a case out of every group (distance 0 in the default mode, its
  distance to each centre with `allcenter=True`, never a medoid). This extension
  was undocumented; the behaviour is unchanged and is now described in the
  docstring and the parity ledger.
- **`dissdomassoc` reads an asymmetric matrix as TraMineR does (2026-09-22,
  audit S4-09, in part).** TraMineR takes each domain matrix and the joint
  matrix from their lower triangle (`as.dist`); traminex took the upper one.
  Results for symmetric matrices are unchanged; for an asymmetric one (possible
  with asymmetric substitution costs) the correlations now match R.
- **Missing cells in `seqdef_from_long` and `TSE_to_STS` are markers, not
  states (2026-09-22, audit S2-05).** A period with no row for a person, or a
  missing state value, became a real state spelled `"nan"`, `"None"` or
  `"<NA>"`, and in such an unbalanced panel integer states read `"1.0"`.
  `TSE_to_STS` in event mode did the same with the cells an NA transition in
  `stm` leaves. They now follow `seqdef`'s default rules, as R does when
  `seqdef` reads `TSE_to_STS`'s result: leading and internal missing cells
  become the nonresponse marker and trailing ones void, and integer states keep
  their integer text. Such cells are refused by name when `nr=None` or
  `void=None` disables the marker they need.
- **The long-to-wide pivot honours an ordered categorical time column
  (2026-09-22, audit S2-06).** `seqformat(from_="long", to="wide")`, and so
  `seqdef_from_long` and `TSE_to_STS` in state mode, sorted the time labels
  with Python's `sorted()`, which put `t10` before `t2` even when the time
  column was an ordered categorical declaring `t1, ..., t12`. Columns now
  follow the time column's own order: categories for a categorical, value
  order for numbers and dates (unchanged). Plain string labels still sort as
  text; the `seqdef_from_long` docstring now says so and how to avoid it.
- **Tied events and person ids follow R's collation (2026-09-22, audit
  S2-07).** `FCE_to_TSE` and `TSE_to_STS` ordered events that share a time,
  and `toPersonPeriod` ordered person ids, by Python code points, so `C`
  came before `b` and `_` after the capitals. R's `order()` sorts these
  character keys by the locale's collation (ICU: `_x, a, b, B, C`). With an
  explicit transition matrix the order of tied events decides the resulting
  state, so `TSE_to_STS` could return a different state than R (or refuse a
  transition R makes); `FCE_to_TSE` event rows (and character times) and
  `toPersonPeriod` rows came back in a different order. All three now use the
  ICU root collation traminex already uses for `seqstatl`, which makes the
  system ICU library (`libicu74`) a requirement for them too.
- **`seqMD` orders its combined alphabet as R does (2026-09-22, audit
  S2-08).** The combined states were sorted by Python code points, so with
  mixed-case state names `B+Y` came before `a+x` in the combined sequence
  object and in the rows and columns of `what="cost"`. R takes `seqdef`'s
  default alphabet, sorted in the locale's collation (`a+x, B+x, B+Y`);
  `seqMD` and `seqdistmc` now use the same ICU collation. Distances are
  unchanged.
- **A real state spelled `"*"` is no longer turned into a missing value
  (2026-09-22, audit S2-09).** The SPS, SPELL and TSE routes keep a literal
  `"*"` in the data as an ordinary state, but three functions then matched the
  missing-value marker by its text. `STS_to_SPELL` and
  `seqformat(to="SPELL")` dropped the real `"*"` spells (`with_missing=False`)
  or merged them into the adjacent missing spell (`with_missing=True`);
  `seqstart` turned the real `"*"` cells into missing cells; and
  `seqgen_missing` added a second `"*"` to the alphabet. The markers are now
  told apart by their codes, and a marker that must be added takes a spelling
  no state uses (`"<NR>"`), as `seqdef` already does; R's `seqstart` likewise
  fills with `NA`, distinct from a literal `"*"`.
- **`seqdef` documents how it reads a literal marker spelling (2026-09-22,
  audit S2-10).** On the STS route a cell or `alphabet` entry spelled like
  `nr` or `void` is read as that marker (so `alphabet=["A", "*"]` with
  `nr="*"` declares the marker), while the SPS and SPELL routes keep it as a
  real state. R lists such a spelling as a state but then also counts its
  cells as missing. The behaviour is unchanged; the docstring and the parity
  record now state it and how to keep a literal `"*"` state.
- **`seqMD` pads a shorter channel with void, as R does (2026-09-22, audit
  S2-11).** With channels of unequal width, the columns beyond a shorter
  channel's end were filled with that channel's missing-value marker whenever
  the channel had missing values, giving combined states such as `*+y` (and
  `*+%` where R has a void cell). R pads with void and applies
  `fill.with.miss` only to the channel's own columns, so these are now `%+y`
  and void.
- **`seqMD(with_missing=False)` refuses a channel that has missing values
  (2026-09-22, audit S2-12).** For `what="cost"` and `what="diss"` it only
  warned and went on to build costs that include the missing state, where R
  stops ("found missing values in channel i"). It now raises `ValueError`
  with that message; `what="MDseq"` still only warns, as in R.
- **Sequence objects keep their own copy of the weights (2026-09-22, audit
  S2-13).** `seqdef(weights=...)`, the `STSList` constructor and assigning
  `x.weights` stored a float64 NumPy array by reference, so changing the
  caller's array afterwards silently changed the object's weights. They now
  copy it, as R does.
- **`stlab` and `cpal` agree between the two marker representations
  (2026-09-22, audit S2-14).** For an object whose `nr`/`void` codes are
  declared outside the alphabet they returned only the real states' entries,
  while the same data with in-alphabet markers returned the documented
  codebook-padded list (markers included). Both now return the padded list.
- **The `STSList` constructor refuses misaligned time labels and a duplicated
  alphabet (2026-09-22, audit S2-15).** It checked the length of `labels`
  but accepted an `xtlab` whose length differs from the number of positions
  and an alphabet naming the same state twice, which gave shifted time labels
  or two codes no function could tell apart. Both now raise `ValueError`, as
  `seqdef` already does. A `cpal` covering only the real states is still
  accepted; plotting pads it with a warning.
- **`seqMD_gridsearch_weights(return_all=True)` reports every weight vector
  on the scale it tested (2026-09-22, audit S2-16).** Each vector is rescaled
  to sum to the number of channels before it is scored, and `best_weights` was
  reported on that scale, but `all_weights` listed the grid as given (summing
  to 1), so the best vector was not among the reported ones. `all_weights` is
  now the rescaled vectors; scores are unchanged.
- **Calibrated-parameter files and `OMConfig` validate their input
  (2026-09-22, audit S2-17).** `load_calibrated_params` raised
  "'PosixPath' object is not iterable" for a single `pathlib.Path`; it now
  accepts one. A `trate.smoothing` in the file that is negative, NaN,
  infinite or not a number surfaced later as an unrelated cost error ("OM
  indel cost must be strictly positive"); it is now refused by name,
  naming the file. `OMConfig(indel=True)` and `OMConfig(subm=True)` were
  silently used as 1; they now raise `TypeError`.
- **`seqdef` documents where its time labels come from (2026-09-22, audit
  S2-18).** `seqdef` takes `xtlab` only from `cnames`; R keeps the source
  column names. The behaviour is unchanged; the docstring and the parity
  record now say so, and how to keep the column names
  (`cnames=list(frame.columns)`).
- **`seqMD(what="cost")` returns R's costs when a combined state has a void
  component (2026-09-22, audit S2-19).** A combined state such as `%+y` --
  from the padding of a shorter channel, or a trailing void in a channel
  without missing values -- has no CAT cost. R returns the cost matrix with
  `NA` for every cost pairing that state with another (and 0 on its
  diagonal); traminex raised "combined state component '%' is not in channel
  1 costs". It now returns R's matrix, with `NaN` for `NA`. `what="diss"`
  still refuses such data, now with R's message ("NA indel and/or
  substitution CAT costs, cannot compute MD distances!").
- **`seqecmpgroup` orders tied subsequences as R does (2026-09-22).** R's
  `chisq.test` applies one Yates correction, `min(0.5, abs(x - E))`, to the
  whole 2 x 2 table and adds up the table, its margins and the statistic's
  terms in long double, in R's order. traminex corrected each cell separately
  and summed in double precision, so statistics that are equal in exact
  arithmetic came out equal in traminex but a few last bits apart in R (tiny
  residues, about 1e-31, stand in for zero there), and tied subsequences were
  listed in a different order. The group tables, statistics, expected counts
  and residuals now match R bit for bit, and so does the order. Only the
  p-values, which come from SciPy rather than R's `pchisq`, still agree only to
  rounding (audit S5-16).
- **`seqpcplot` filters whole trajectories, as TraMineR does (2026-09-22,
  audit S6-02).** A `seqpcfilter` was evaluated on one value per *sequence*,
  so `cumfreq` kept only as many sequences as it took to reach the level: with
  a trajectory holding 60% of the sequences, level 0.5 drew one of its three
  sequences, and level 0.7 drew two of them and dropped the next trajectory.
  The filter now sees, as in `seqpcplot.R`, one weighted share per distinct
  trajectory (simultaneous events forming one step) plus the share of
  sequences without events, and draws every sequence of each trajectory it
  keeps. `minfreq` results change only where simultaneous events now make a
  trajectory of their own; all three filter methods match live TraMineR
  2.2-14.
- **Tree DOT labels format numbers as TraMineR does (2026-09-22, audit
  S6-03).** `disstree2dot`, `disstree2dotp`, `seqtree2dot` and the display
  functions formatted node sizes, discrepancies and R2 with Python's `%g`
  rule, which turns to scientific notation as soon as a number has more
  integer digits than `digits`: a node of 1,200 cases read `n: 1.2e+03`
  where TraMineR writes `n: 1200`, and `1234.5` read `1.23e+03` instead of
  `1234`. Labels now follow R's `format(x, digits = digits)`, which keeps
  fixed notation unless it is wider than scientific notation. The
  global-quality lines now show five significant digits
  (`getOption("digits") - 2`, as in `disstree2dot.R`) whatever `digits` is.
- **`duration_boxplot(by_state=True)` splits by state without seaborn too
  (2026-09-22, audit S6-04).** With `seaborn=False` the per-state split and
  the cluster `labels` were silently ignored: the plot showed a single box
  pooling every sequence-and-state cell, most of them zeros from states a
  sequence never visits and from the missing and void marker columns. The
  matplotlib path now draws one box per state (per state and cluster when
  `labels` is given). Both paths now leave out the missing and void markers,
  which the seaborn violins listed as if they were states.
- **`seqalignplot` draws both marker representations alike (2026-09-22,
  audit S6-12).** It did not canonicalize its input, so for an object with a
  declared external missing code the missing cells were painted in the
  alignment-gap grey instead of the missing state's palette colour. It now
  canonicalizes at entry, as `docs/PLOTTING.md` promises for every renderer;
  `duration_boxplot`, the other renderer the audit named, draws both
  representations alike since the `by_state` fix above.
- **`seqmsplot` reports the modal sequence's weighted frequency
  (2026-09-22, audit S6-05).** The "Modal state sequence (k occurrences,
  freq=…%)" annotation divided by the number of rows, so weighted data
  showed 25.0% where TraMineR shows 17.8%. It now divides by the sum of the
  weights, as `plot.stslist.modst.R` does, and prints the rounded value as
  R does (`25%`, not `25.0%`).
- **`seqedplot(group=...)` orders and fills its panels as TraMineRextras
  does (2026-09-22, audit S6-06).** Group values were compared as text, so
  numeric groups 1, 2 and 10 gave panels in the order 1, 10, 2; missing
  values became panels titled `None` or `nan`; and supplying fewer `axes`
  than groups silently dropped panels. Panels now follow R's
  `factor(group)`: numbers in numeric order and printed as R prints them,
  text in sorted order, a categorical in its declared order, and sequences
  with a missing group value in no panel. A `group` of the wrong length or
  an `axes` list of the wrong size is refused.
- **`seqHtplot(group=...)` leaves sequences with a missing group out
  (2026-09-22, audit S6-07).** A `NaN` group value added an all-`NaN` line
  labelled `nan`, a `None` value raised `TypeError`, and a group of the
  wrong length was not refused. As in TraMineR, sequences with a missing
  group now belong to no group line; the group levels follow R's `factor`
  order and labels (`1`, not `1.0`), and a group of the wrong length is
  refused.
- **A failed plotting call no longer leaves an open figure behind
  (2026-09-22, audit S6-09).** Several renderers created their figure
  before checking their arguments or applying caller styles, so an error
  left the figure in pyplot's registry: twenty failing `plot_seqimplic`
  calls left twenty figures open. Every public plotting function now
  closes the figures it opened when it raises; figures the caller
  supplied are left open.
- **Legends show the curves' colours when no palette is set (2026-09-22,
  audit S6-10).** For a sequence object without `cpal` (or a `SeqImplic`
  without one), `seqsplot` and `plot_seqimplic` let matplotlib colour the
  curves but drew every legend entry in the first default colour; in
  `seqsplot` a state could also change colour between panels. State
  curves now use the default state palette the other renderers use, and
  `seqsplot(per_state=True)` gives more than eight groups explicit
  colours, so each legend entry matches its curve.
- **State legends name states by their labels, as `seqlegend` does
  (2026-09-22, audits S6-11 and S6-14).** `state_legend`,
  `state_legend_figure` (and so index-plot legends) and the
  `state_distribution_plot` legend listed the internal codes followed by
  the missing and void markers (`A, B, C, D, *, %`). They now list each
  state's label, add a `missing` entry only when missing values are
  present (for the distribution plot, only when they are plotted with
  `include_missing=True`) and never list void, matching live TraMineR
  2.2-14. `state_legend` returns those displayed labels. With a palette
  shorter than the alphabet, the legend helpers raised `IndexError` or
  `ValueError`; they now warn and fall back to the default palette exactly
  as the plots do.
- **A colormap gives every plot the same, distinct state colours
  (2026-09-22, audit S6-13).** A named matplotlib colormap such as
  `cmap="viridis"` coloured state *k* with entry *k* of its 256-colour
  table in the bar and band plots (`state_distribution_plot`, `seqfplot`,
  `seqmtplot`, `seqmsplot`, `seqrplot`, `seqrfplot`, `seqsurvplot`,
  `seqalignplot`), so six states came out in six nearly identical dark
  purples, while the index plot spread them over the whole map. A colormap
  (by name or as an object) is now sampled once per state, at the colours
  the index plot uses, and a colour list colours state *k* with entry *k*
  (recycled when shorter) in every renderer.
- **Plotting functions refuse a wrong input type with `TypeError`
  (2026-09-22, audit S6-15).** `seqfplot`, `seqrfplot`, `seqplot_rf`,
  `seqsurvplot`, `seqalignplot`, `seqplot_tentrop`, `cluster_chronograms`,
  `cluster_chronograms_with_representatives`,
  `sequence_index_plot_grouped` and `overlay_representatives_on_index_plot`
  raised `AttributeError` when given a plain array instead of a sequence
  object; `seqedplot` did the same for a non-event-sequence input and
  `seqplotMD` for a channel list holding plain arrays. They now raise the
  same `TypeError` as the other entry points.
- **Large index plots render quickly (2026-09-22, audit S6-16).**
  `sequence_index_plot` labelled every displayed row on the y axis, and
  laying out those labels dominated the run time: 5,000 rows took about 53
  seconds to render and save, against 0.6 seconds now. At most 50 evenly
  spaced rows (the first and last included) are labelled, each with its own
  identifier; plots of up to 50 rows keep every label.
- **Survival plots accept a short palette as the other plots do
  (2026-09-22, audit S6-14).** A sequence object whose `cpal` has fewer
  colours than states made `seqsurvplot` and `seqsplot` fail with
  `IndexError` inside `seqsurv`, while every other renderer warns and falls
  back to the default palette. `seqsurv` now leaves such a palette
  unresolved (`cpal=None` in its result) and both plots warn and fall back
  like the rest, so each state keeps the same colour across plots.
- **Tree display copes with unusual file and column names (2026-09-22,
  audit S6-17).** `disstreedisplay` and `seqtreedisplay` passed relative
  paths to Graphviz as given, so a file name starting with `-` (for example
  `"-Tsvg"`) was read as a `dot` option and rendering failed; the paths are
  now passed in absolute form. A double quote in a predictor name, a
  sequence-legend label or an image path ended the DOT string early and
  made the file invalid (TraMineR writes such names raw too); it is now
  escaped, and DOT output for all other names is unchanged.
- **`disstreedisplay` and `seqtreedisplay` accept an image file name, as
  in TraMineR (2026-09-22, audit S6-18).** `filename` was always treated as
  a stem, so the R-style call `disstreedisplay(tree, "tree.png")` wrote
  `tree.png.png` and `tree.png.dot`. A `filename` that ends with the
  extension of `image_format` is now the image path itself (with the `.dot`
  file beside it as `tree.dot`); other names keep the stem behaviour. The
  parity ledger now also records that traminex asks Graphviz for
  `image_format` directly, where R renders PNG or JPEG and converts other
  formats with ImageMagick.
- **`save_plot(bbox_inches="standard")` works as documented (2026-09-22,
  audit S6-19).** The helper documented `"standard"` for default margins but
  passed it straight to matplotlib, which raised `AttributeError`. It now
  keeps matplotlib's default margins for `"standard"` (and for `None`).
- **`chronogram(seaborn=True)` centres its time labels on the columns
  (2026-09-22, audit S6-21).** The seaborn heatmap draws each time position
  as a cell one unit wide, but the time labels were placed on the cell
  edges, half a cell to the left of their columns. They now sit on the cell
  centres, as on the default (non-seaborn) path.
- **`seqplot` honours `weighted` for every plot type (2026-09-22, audit
  S6-20).** `seqplot(type="d", "Ht", "dH" or "f", weighted=False)` raised
  `TypeError` because those renderers have no `weighted` option, while
  TraMineR passes it to `seqstatd` or `seqtab`. `seqplot` now draws the
  unweighted series for `weighted=False` and the weighted one otherwise,
  without changing the caller's object.
- **The maintained type stubs pass mypy again, and the stub gate now checks
  them (2026-09-22, audit S7-01).** `stubs/traminex/__init__.pyi` re-exported
  `disstree2dotp`, `disstreedisplay`, `seqtree2dot` and `seqtreedisplay` from
  stub modules that did not declare them, and the `SeqIVarDurResult` stub
  lacked `__array_ufunc__`, which made the class abstract. Type-checking the
  stub tree (`MYPYPATH=stubs mypy -p traminex`) reported 5 errors and treated
  the four names as untyped. The stubs now declare them, and `make stub-check`
  runs that mypy command. The installed package, which relies on its inline
  annotations, was not affected.
- **The maintained stubs declare the real class constructors (2026-09-22,
  audit S7-02).** The `DURConfig` stub still listed the `eps` field that the
  runtime class no longer has, so `DURConfig(metric="chi2", eps=0.1)` passed a
  type check and then raised `TypeError`. Nine result classes (`DissAssocResult`,
  `DissTree`, `DissTreeNode`, `DissTreeSplit`, `SeqElist`, `SeqStatDResult`,
  `SeqIVarDurResult`, `SeqImplic`, `SeqPcFilter`) had stubs with no constructor,
  the `SeqSurv` stub did not declare its `max_time` attribute, and the
  `SeqImplic` stub declared a `max_time` attribute the runtime class does not
  have. The stubs now match the runtime classes, and
  `scripts/check_stub_signatures.py` compares class constructors and dataclass
  fields as well as function signatures, so this drift fails `make stub-check`.
- **Three documented examples no longer raise (2026-09-22, audit S7-04).**
  `docs/NEW_DISTANCE_METHODS.md` (`OMstran`) and `docs/USAGE_new_apis.md`
  (`seqrep` without a distance matrix, and `seqdist(method="OM")`) omitted the
  substitution cost that `seqdist` requires for the OM family, just as TraMineR
  requires `sm`, so each raised `ValueError`. The examples now pass
  `indel=1.0, subm=2.0`, and `tests/test_doc_examples_execute.py` runs every
  Python example in the six fully runnable documents so this cannot recur
  unnoticed.
- **The development guide states that the test suite runs from a git
  checkout (2026-09-22, audit S7-03).** The source distribution ships `tests/`
  and `scripts/`, but some tests call `git ls-files` and the semantic-ledger
  and documentation checks read repository-only evidence that the sdist
  deliberately excludes, so running `pytest` in an unpacked sdist fails.
  `docs/DEV_BUILD.md` now says so; no packaged file changed.
- **`help(traminex)` and `inspect.getmembers(traminex)` work without
  matplotlib (2026-09-22, audit S7-05).** Both call `getattr` on every name
  `dir(traminex)` lists and skip only `AttributeError`, so on a core install the
  informative `ImportError` raised by the 37 plotting names made both fail.
  `dir(traminex)` now leaves out the plotting names while matplotlib cannot be
  found. `__all__` still lists all 226 names and accessing a plotting name
  still raises the `ImportError` naming the `plot` extra, so
  `from traminex import *` continues to raise it on a core install
  (documented in `docs/API_CONSISTENCY.md`).
- **`import traminex.<name> as alias` is documented for the five names shared
  by a function and its submodule (2026-09-22, audit S7-06).** For
  `disstree`, `seqdiff`, `seqemlt`, `seqformat` and `seqtree`, Python resolves
  `import traminex.seqtree as st` with the same attribute lookup as
  `from traminex import seqtree`, so `st` is the public function, not the
  submodule. That is unavoidable while `traminex.seqtree` stays the function,
  and it was undocumented; `docs/API_CONSISTENCY.md` now explains it and the
  working alternatives. The documentation checker also bound such an alias to
  the submodule, so it would have accepted a documented `st.seqtree_fit(...)`
  call that fails at run time; it now binds it as Python does.
- **`seqrep`'s deprecated-argument warnings are shown again (2026-09-22,
  audit S7-07).** `seqrep` is wrapped by an input-validation decorator, so the
  fixed stack level of the shared deprecation helper attributed
  `'trep' is deprecated for seqrep()` (and the `tsim` / `dist_matrix`
  warnings) to a file inside the package, where Python's default warning
  filters hide a `DeprecationWarning`. TraMineR always shows the equivalent
  warning. The helper now attributes every deprecated-argument warning to the
  first caller outside the package.
- **Documentation that had drifted from the code is corrected (2026-09-22,
  audits S7-08 and S7-15).** `TRAMINEX_BINDING` and `TRAMINEX_NUMBA` are read
  when `traminex.distances` is first imported, not by `import traminex`
  (`docs/USAGE_new_apis.md`, `docs/DISTANCES.md`). `om_trate_pairwise`
  defaults to `indel=None` (`docs/DISTANCES.md`), and `seqtree` takes
  `(seqdata, predictors, ...)` rather than R's `(formula, data, ...)`
  (`docs/SEQTREE.md`). `docs/README.md` now lists the same bare `mypy` command
  as `make typecheck`; `docs/INDEX.md` places `frequent_ngrams` in `events/`,
  describes the `.qmd` files as include wrappers, and dates its oracle
  reference. `ROADMAP.md` points two completed items at the current
  `traminex/events/` package and tests; `CONTRIBUTING.md` no longer calls the
  100-character line length black's default or refers to a missing
  `CONTRIBUTORS.md`; `ECOSYSTEM.md` records the 2026-09-20 oracle change.
- **The wheel workflow's verification step installs the wheel's dependencies
  (2026-09-22, audit S7-10).** `.github/workflows/wheels.yml` installed the
  built wheel with `--no-deps` into a runner that had only build tools, then
  imported `traminex.distances`, which needs NumPy, pandas and SciPy; the same
  sequence fails with `ModuleNotFoundError: No module named 'numpy'` in a fresh
  environment. The step now lets pip install the declared dependencies, and
  `scripts/validate_workflows.py` rejects a wheel installed with `--no-deps`.
  Checked statically and in a local environment; no CI run was made.
- **The workflow policy check reads the job structure (2026-09-22, audit
  S7-09).** `scripts/validate_workflows.py` accepted a workflow whenever the
  expected words appeared anywhere in its text, so a `windows-2022` job, a
  Python 3.13 matrix entry or a make target that does not exist passed. It now
  requires every job's `runs-on` to be an Ubuntu runner, every configured
  Python version (matrix or `setup-python`) to be 3.11 or 3.12, every `make`
  target named in a `run` step to exist in the `Makefile`, and, for `ci.yml`
  and `wheels.yml`, both supported versions and (in `ci.yml`) an actual
  `make release-check` command.
- **`duration_boxplot(by_state=True)` width-normalises its violins on seaborn
  0.12 too (2026-09-22, audit S7-11).** The plot passed
  `density_norm="width"`, a keyword seaborn introduced in 0.13. seaborn 0.12,
  which the `plot` extra still allows, silently ignored it and drew
  area-scaled violins (measured widths 0.29 / 0.8 / 0.27 instead of 0.8 each).
  The plot now passes `density_norm` or the older `scale` keyword, whichever
  the installed `seaborn.violinplot` declares, so no deprecated argument is
  used on seaborn 0.13 or later.
- **`check_stub_signatures.py --runtime-root` checks the tree it is given
  (2026-09-22, audit S7-12).** It parsed the source files under
  `--runtime-root` but imported the package through Python's import system,
  where an editable install's finder outranks `sys.path`, so a copy of the
  package was compared against the installed checkout's objects and could pass
  while drifting from its stubs. The checker now puts ordinary path lookup
  first, as `make test-pybind11` already does, and reports any module that
  still comes from outside `--runtime-root` instead of comparing it.
- **Cited audit evidence is tracked, and the documentation check enforces it
  (2026-09-22, audit S7-13).** `CHANGELOG.md` and `docs/API_PARITY.md` cited
  the 2.2-14 substitution audit's directory, which `.gitignore` excluded, so
  the citation resolved only on the author's machine; and the twelve audit and
  remediation runs from 2026-08-26 to 2026-09-20 had not tracked their reports
  as `docs/DOCS_POLICY.md` requires. Their reports and evidence indexes (22
  small Markdown files) are now tracked, the two citations name the tracked
  report, and `scripts/check_docs_consistency.py` rejects an active document
  that cites an untracked path under `reports/audit/<run>/`.
- **`stubs/README.md` no longer carries stale counts (2026-09-22, audit
  S7-14).** It stated 403 public functions, 32 classes, 312 signatures and 78
  package source modules after the code had moved on (402, 33, 311 and 80 at
  this change). It now points to the totals the check prints. The
  documentation check's pattern for hard-coded module counts also missed the
  "package source" wording, and the file was not among the checked documents;
  both are fixed.
- **A seaborn-backed plot without seaborn says how to fix it (2026-09-22,
  audit S7-16).** `seqfrequency_plot` (seaborn by default), and
  `duration_boxplot(by_state=True)`, `chronogram(seaborn=True)` and
  `transition_heatmap(seaborn=True)`, import seaborn inside the function,
  where the package's "install the extra" message can never see the failure,
  so a user with matplotlib but not seaborn got a bare
  `No module named 'seaborn'`. They now raise an `ImportError` that names the
  `plot` extra and the `seaborn=False` opt-out. The unreachable scikit-learn
  and numba entries were removed from the package's extras map: those imports
  are guarded, and their features either fail with their own message
  (`scikit-learn is required for ...`) or fall back to the reference kernels.
- **The parity-mapping generator reads R NAMESPACE files with the shared
  parser (2026-09-22, audit S7-17).** `scripts/generate_api_mappings.py`
  split the TraMineRextras NAMESPACE with its own regular expression, which
  kept R comments and produced the bogus name `##seqeordplot`. It now uses
  `scripts/_r_namespace.parse_namespace`, as the API audit and the semantic
  ledger do. The generated mapping files are unchanged.
- **The audit tooling no longer depends on the working directory
  (2026-09-22, audit S7-18).** `scripts/audit_traminex.py` read its inputs
  relative to the current directory as soon as it was loaded, so even
  `--help` failed with `FileNotFoundError` outside the repository root, and
  its default output directory was relative too. `generate_parity_table.py`
  and `calibrate_parity.py` defaulted to paths relative to the current
  directory as well. All three now anchor their inputs and default outputs to
  the repository that contains them.
- **`scripts/annotate_traminer_list.py` writes LF line endings (2026-09-22,
  audit S7-19).** It wrote `reports/parity/traminer_API_list_with_traminex.csv`
  with the `csv` module's default CRLF endings although `.gitattributes`
  declares LF for CSV files, so regenerating the file rewrote every line. The
  regenerated file is now byte-identical to the tracked one.
- **The benchmark names the real reason the Numba backend is unavailable
  (2026-09-22, audit S7-20).** With `TRAMINEX_NUMBA=off`,
  `benchmarks/run_benchmarks.py --backend numba` told users to install the
  `perf` extra even when numba was installed. It now says that numba is
  installed but dispatch is disabled by `TRAMINEX_NUMBA=off`. An unused
  constant, whose explanation of the Numba comparison tolerance now sits next
  to the tolerance itself, was removed.
- **The documentation check covers Quarto chunks and the oracle version, and
  the Quarto examples are executed (2026-09-22, audit S7-21).**
  `scripts/check_docs_consistency.py` checked only ```` ```python ```` fences,
  never executable ```` ```{python} ```` chunks, and did not enforce
  `docs/DOCS_POLICY.md`'s rule that active parity documents name the current
  oracle. It now reads both kinds of block and requires `docs/API_PARITY.md`,
  `docs/INDEX.md` and `scripts/r_goldens/README.md` to name the TraMineR and
  TraMineRextras versions vendored under `SOURCES/`. The 28 examples in
  `docs/TraMineR_examples.qmd` are plain `python` fences, which Quarto shows
  without running, so `tests/test_doc_examples_execute.py` now runs them; the
  sentence in `docs/INDEX.md` that said they run at render time is corrected,
  and a stale comment about the paired `.qmd` files is updated.
- **`make lint` checks every tracked Python file, and a stale lint
  exception is gone (2026-09-22, audit S7-22).** `make lint` passed
  `traminex tests scripts benchmarks` while `make format-check` also covered
  `examples`, `docs` and `verify_audit_improvements.py`, so ruff never saw
  those files; it now uses the same paths. The `F401`/`E402` exception for
  `traminex/__init__.py` in `pyproject.toml` was no longer needed (ruff passes
  without it) and was removed. `README.md` and `docs/DEV_BUILD.md` now show
  `make lint`.
- **The install documents name every function that needs the ICU library
  (2026-09-22, audit F-08).** `README.md` and `docs/DEV_BUILD.md` listed
  `seqstatl`, `seqdiff`, a vaguely named "grouping path", the event
  converters and `seqMD`; `seqtab` (tied frequencies), `seqtabstocc`,
  `rowmode`, and since today `seqimplic` and `seqsurv` (so `seqsplot` and
  `seqsurvplot`) also order labels with ICU and raise without it. Both
  documents now name them, and a test keeps the list complete.
- **`disstree` stops on missing values where TraMineR does (2026-09-22, audit
  S4-16).** A missing dissimilarity between two cases, a missing weight, or a
  predictor level whose cases all weigh 0 makes a node's discrepancy
  undefined. TraMineR then stops, but traminex skipped the predictor and
  silently returned a tree with the root alone. `disstree` (and `seqtree`) now
  raise `ValueError` in that case. As in R, the root is still returned alone
  when no split is attempted, and the diagonal of the matrix is ignored: a
  missing value there no longer turns the whole tree into NaN.
- **`dissrep` reports an undefined mean distance as NaN (2026-09-22, audit
  S4-17).** For a representative whose assigned sequences all have weight 0,
  the mean distance `MD` was reported as 0. As in TraMineR (`MD = SD/na`), it
  is now NaN, since there is no weight to average over.
- **`disstree` no longer warns on a root-only tree with a NaN discrepancy
  (2026-09-22, audit S4-18).** The split schedule averaged tied ranks over
  every distinct value, NaN included, and NumPy warned "Mean of empty slice";
  R prints nothing. The schedule itself was already R's.
- **`seqMD(what="cost")` ignores `method`, as TraMineR does (2026-09-22,
  audit S10-12).** With `what="cost"` and a cost-method string,
  `method="DHD"` made the costs time-varying, so `seqMD` returned a
  states-by-states-by-time array where TraMineR returns the stationary
  matrix, and an explicit two-dimensional `sm` was refused. Only
  `what="diss"` with `method="DHD"` is time-varying now; for `what="cost"`
  the dimensions of `sm` decide, as in R.
- **`seqMD`'s combined costs are symmetric for an asymmetric channel matrix
  (2026-09-22, audit S10-13).** TraMineR computes each combined cost once for
  a pair of states and mirrors it; traminex took the transposed channel
  costs below the diagonal, so with an asymmetric channel matrix its
  `what="cost"` costs and `what="diss"` distances differed from R's. They now
  match R, for stationary and time-varying costs.
- **`seqMD` reports TraMineR's cost diagnostics for each channel (2026-09-22,
  audit S10-11).** With `what="cost"` or `"diss"`, TraMineR's `seqMD` checks
  every channel's costs -- explicit, or generated from `"TRATE"` and the other
  method strings -- before combining them, and reports zero costs, a largest
  cost of 0 or below, triangle-inequality violations, costs above twice the
  indel and asymmetry. traminex reported none of these for the channels; it
  now gives the same diagnostics in R's order and wording, with R's
  `message()` diagnostics as `CostDiagnosticMessage` warnings. No cost or
  distance value changes.

### Added

- **`TRAMINEX_NUMBA=auto|off` (2026-08-23).** Numba is a *third* distance backend,
  active by default whenever `numba` imports and orthogonal to `TRAMINEX_BINDING`,
  which selects only between the Python and compiled implementations. There was
  no supported way to select the reference kernels. `TRAMINEX_NUMBA=off` now does
  so; like `TRAMINEX_BINDING` it is read **once, at import time**. A new
  agreement suite compares both backends across all seven production dispatch
  sites: seven of eight routes are byte-identical and only `JSD_TRANS` differs, at
  ~1e-16.
- **Parameter-naming and result-orientation conventions**, in
  `docs/API_CONSISTENCY.md`: how R argument names map to Python (dots,
  camelCase, keyword collisions, R's leading `seqdata` spelled `x`) plus a table
  of the semantic renames no rule predicts; and that per-sequence indicators are
  1-D here where R returns an `n x 1` matrix.
- **`r_kind = internal` in `docs/API_PARITY_DECISIONS.csv`**, for R objects
  traminex promotes to its public API although R defines them without exporting
  them. Seven such surfaces previously carried no authoritative disposition
  because the vocabulary could not express them. Membership is verified against
  the vendored R source.

### Changed

- **Parity dispositions corrected (2026-08-23; behaviour unchanged).** `seqmeant`,
  `cpal` and `stlab` move `EXACT` -> `DIVERGENCE`: their length follows the
  internal codebook, which appends the `nr` and `void` markers, where R returns
  one entry per real state. Values agree exactly at 15 digits. The padding is
  **not** uniform — `seqstatf` and `seqtrate` follow the real alphabet exactly as
  R does — so it is recorded per function. `[.stslist` and `rbind.stslist` no
  longer claim an `STSList` indexing or container protocol: none exists, and row
  selection is done by slicing `.data` and reconstructing.
- **Self-contained R-source reimplementation (2026-07-22; public API and
  analytical contracts unchanged).** Ten formerly Sequenzo-adapted Python/C++
  implementation units were independently reimplemented from the exact
  vendored TraMineR 2.2-13 and TraMineRextras 0.6.9 sources, plus the exact
  WeightedCluster 2.0 ASW authority used by `dissmergegroups`. The compiled
  HAM/LCS/OM kernels now use direct width-bounded dynamic programming and a
  project-original static OpenMP row partition, with no prefix/suffix trim,
  bit-parallel masks, or Sequenzo scheduling helpers. Historical notices and
  dated entries remain retained as provenance for earlier revisions.
- **Integrity-remediation follow-through (2026-07-18; public API unchanged).**
  `omloc` now rejects invalid weights and localized-cost options consistently;
  even-width smoothing preserves sequence length; and duration boxplots use the
  Matplotlib-compatible tick-label path across the supported dependency range.
  Active transition, indicator, event-mining, and multidomain examples now use
  the current public contracts. Lazy root-API resolution now also preserves
  callable exports when a dependency imports a same-named submodule. The
  documentation gate binds explicit traminex examples against their runtime
  signatures so stale or incomplete calls fail the release check.
- **Stable `traminex.distances` submodule boundary (2026-07-15; root API count
  unchanged).** The supported dispatcher callables now have an explicit
  submodule `__all__` inventory and matching stub/documentation contracts.
  `pairwise` and `pairwise_between_sets` retain their existing signatures and
  dense/ragged, normalization, deduplication, and rectangular behavior. They
  intentionally remain submodule-only rather than expanding the 226-name root
  manifest.
- **`dissmfacw` mixed-predictor parity (2026-07-15; public signature and root
  API count unchanged).** Numeric predictor arrays now remain one continuous
  model-matrix column instead of being incorrectly dummy-coded as categories;
  string/object, boolean, and pandas categorical inputs retain R treatment
  coding. Mixed categorical/numeric term order, weighted and unweighted
  PseudoF/PseudoR2/explained-SS results, complete-case subsetting, constants,
  and weight validation are pinned to live TraMineR 2.2-13 goldens, including
  the exact design required by pyweightedcluster `clustassoc`. Numeric category
  labels should be passed as strings or a pandas categorical.
- **R5 Ubuntu reliability and release hardening (2026-07-15; no analytical,
  API, ledger, golden, or version change).** The supported matrix is Ubuntu
  24.04 with Python 3.11/3.12. Default sdists/wheels are compiler-free and the
  sole compiled profile is opt-in pybind11/OpenMP. The Make gates and the
  configured CI workflow (never run: no git remote) are Ubuntu-native; obsolete PowerShell and nanobind development paths were
  removed.
- Whole-package mypy now checks the complete configured package module set with zero errors.
  Installed distributions ship `traminex/py.typed`; the maintained external
  stubs are checked bidirectionally across 66 declared modules and against the
  ordered 226-name root manifest.
- Branch coverage now records statements and branches explicitly, with focused
  tests for configuration, backend dispatch, input validation, optional tree
  dependencies, and shared visualization utilities. The four ad hoc timing
  scripts were replaced by one reproducible benchmark CLI that checks backend
  agreement before timing and records environment metadata without imposing a
  performance threshold.
- A maintained artifact verifier hashes and inspects one default sdist, one
  pure wheel, and one compiled wheel, then reuses them across six fresh
  core/performance/plot/tree/compiled/fallback environments. It proves Python
  3.11 dependency floors, Python 3.12 current profiles, import isolation,
  optional boundaries, Graphviz rendering, pybind11/OpenMP identity, automatic
  fallback, and clear explicit-backend failure.

- **Linux and R-oracle refresh (2026-07-12).** Ubuntu 24.04 is now the sole
  development platform. The vendored behavioral references are TraMineR
  2.2-13 and TraMineRextras 0.6.9, verified from the official CRAN source
  archives and compiled with R 4.6.1 on Linux. All 75 R-golden generators were
  rerun. `dissrf` now follows 2.2-13's normalized-weight and `k=1` behavior,
  and `seqstatl` follows the Linux oracle's missing/void ordering.
- Active setup, build, API parity, family, plotting, and source-provenance docs
  were refreshed. `py_ggplotseq` is recorded as the dedicated plotnine
  visualization sibling, and flexscanx's planned broader use of traminex is
  recorded as a later architecture decision.

### Performance

- `seqdist(method="DHD")` with generated costs (the default and `subm="TRATE"`)
  builds its time-varying costs once per call; the costs checked for seqcost's
  messages were built again inside the kernel (TX-REP-20261019-007).
- **RP9 — Acceleration: a compiled, OpenMP-parallel C++ distance kernel behind
  `seqdist` for OM, LCS and HAM (2026-06-25; speed only, no behaviour change, public
  API unchanged at 208). Adapted from Sequenzo's BSD-3 C++ (attribution header on
  `cpp/src/accel_kernels.h`, entry in `NOTICE`); the OMspell sub-family and the `xsimd`
  paths are deliberately not borrowed.**
  - Held to a strict **byte-identity gate**: every accelerated result is **bit-for-bit
    identical** (`np.array_equal`, `atol=0`) to the pure-Python / numba kernel, proven by
    `tests/test_accel_byte_identity.py` across a randomised corpus + edge cases and a
    clean re-pass of the live-R golden suite under the compiled backend. The OM DP border
    accumulates by repeated addition (matching the reference exactly, not `i*indel`); the
    LCS is bit-parallel (integer); HAM is the raw mismatch count; matrix-OM trimming is
    used only for zero-diagonal substitution matrices.
  - Measured speedup vs the numba path (n=1500, L=30): **HAM ≈ 368×, LCS ≈ 8.3×, OM ≈
    1.4×**; vs the previous single-threaded kernel, OM ≈ 31× and LCS ≈ 138×. The compiled
    path needs no numba.
  - **Opt-in build (default unchanged):** the default `pip install` stays pure-Python; the
    kernel is built by enabling `-DTRAMINEX_BUILD_PYBIND11=ON` and selected via
    `TRAMINEX_BINDING=pybind11`. `CMakeLists.txt` now links OpenMP when available (the
    OpenMP build adds a `libgomp` runtime dependency; see `docs/DEV_BUILD.md`).

### Added

- **Event-duration plots — completing the one live event-sequence plot RP7a left
  out (2026-06-24; BUILD-from-R, no Sequenzo/ggseqplot borrow; structure tests
  keyed to live-R-golden'd series, NOT pixel baselines). +2 public symbols
  (206 → 208):**
  - **`seqedplot`** (TraMineRextras event-duration plot): one coloured curve per
    non-ignored event — `type="survival"` (the unweighted Kaplan-Meier curve of
    the time to each event's first occurrence, censored at `seqelength`) or
    `type="hazard"` (the weighted mean number of occurrences per equal-width age
    period) — with `group` faceting, `ignore`, an `ages` window and `breaks`.
  - **`seqeage`** (events-package internal, not a public symbol): the event-age
    table, a pure-Python port of TraMineR's unexported `seqeage` / C
    `tmreventinseq`, golden-pinned to live R. It is `seqedplot`'s only missing
    dependency.
  - The two producing series are golden-pinned to live R (survival vs
    `summary(survfit(...))`, unweighted — proven to ignore weights; hazard vs
    `xtabs` over R's `cut()` break-points, weighted/binned — reproducing R's
    `dig.lab` break-point rounding exactly) and the renderer is structure-tested
    and mutation-checked.
  - **`seqeordplot`**: superseded upstream (suppressed in TraMineRextras 0.4.3,
    2018-08-22, in favour of `seqpcplot`); traminex adds a deprecation shim that
    warns and forwards to the already-provided `seqpcplot` — no analytical port.
- **RP7a — Visualization, matplotlib core (roadmap Phase 7, spec P5; structure
  tests keyed to the already-R-golden'd series, NOT pixel baselines; all
  BUILD-from-R-design — no Sequenzo rendering borrow, so no BSD-3/NOTICE).**
  VERIFY first proved the existing plots draw the R series and surfaced one
  wrong-target (seqfrequency_plot is a simplified bar chart of seqtab, not R's
  seqfplot); the maintainer chose to ADD a faithful seqfplot and keep
  seqfrequency_plot. Each plot its own commit + structure test vs the golden
  series. +17 public symbols (189 → 206):
  - Foundation: `traminex/visualization/series.py` (internal plot-series adapter:
    the `sortv` resolver, reshape-only) and a faithful `sortv=` on
    `sequence_index_plot` (seqIplot "from.start"/"from.end"/vector).
  - Base-R / TraMineRextras: **`seqfplot`** (coloured state-segment frequency
    bands, width∝Percent), **`seqmtplot`** (mean-time bars + 95% t-CI),
    **`seqmsplot`** (modal-state bars), **`seqHtplot`** / **`seqdHplot`**
    (transversal entropy / distribution+entropy; `seqHtplot(group=)` folds the
    TraMineRextras `seqplot.tentrop` grouped entropy), **`seqrplot`**
    (representative bands ∝ coverage), **`seqrfplot`** (the Phase-4 `seqrf` edge:
    medoid index + per-group residual boxplot), **`seqsurvplot`** (KM spell-survival
    step curves), **`seqalignplot`** (alignment op grid), **`seqdiffplot`**
    (per-position discrepancy lines), **`seqplotMD`** (side-by-side multidomain
    panels).
  - Event-sequence plots (consolidated from RP5/RP5b): **`plot_subseqelist`**
    (seqefsub Support/Count bars), **`plot_subseqelistchisq`** (seqecmpgroup
    per-group Freq bars coloured by Pearson residual), **`seqpcplot`** +
    **`plot_eseq`** / **`plot_seqelist`** (parallel-coordinate event plot over the
    faithful Eseq/SeqElist; one polyline per sequence).
  - **`seqplot(type=...)`** dispatcher mirroring `seqplot.R` (i/I/d/Ht/dH/f/mt/ms/
    r/rf; event `pc` via `seqpcplot`).
  - **Out of RP7a:** `seqedplot`/`seqeordplot` (no ported producing series at the
    time — the TraMineRextras event-age/survfit extraction was unported); both
    were completed on 2026-06-24 (see the event-duration plots entry above). The
    `plotnine` `ggseqplot` port is RP7b.
  - **RP7a audit fix (2026-06-21):** `seqrfplot`'s right panel now draws R's
    **weighted** Tukey-hinge box (`wtd.boxplot.tmr`) instead of matplotlib's
    unweighted `np.percentile` — the box hinges/median/whiskers reproduce R's
    weighted five-number summary (golden `seqrf_uw_box`); diverged by up to 1.5
    distance units when the relative-frequency boundary weights are non-uniform.
- **RP5a — Event-sequence framework, R-faithful half (TraMineR 2.2-12 /
  TraMineRextras 0.6.8 goldens).** A gate-2 source audit found the prior event
  framework wrong-target at every level (and the intended Sequenzo borrow itself
  non-faithful — a placeholder miner with fake count methods). The maintainer
  chose "do the R-faithful half now, defer the engine to RP5b". Landed, each its
  own commit + live-R golden:
  - **`seqecreate`** now returns the faithful TraMineR `eseq`/`seqelist` object
    (`Eseq`/`SeqElist`): event dictionary (= R factor levels), per-sequence events
    sorted by (time, code), unit weights, an unset (-1) stored length, and the
    C-exact `as_character`/`str` label strings (`"(A)-2-(B,C)"`, `"1-(A)-2-(B)"`,
    `"(A,B)"`; simultaneous events comma-grouped & code-sorted).
  - **`seqetm`** (state→event transition matrix) and **`seqe2stm`** (event→state
    machine) — BUILD-from-R.
  - Accessors **`seqeid`**, **`seqeweight`** (+ `SeqElist.weights` setter =
    `seqeweight<-`), **`seqecontain`** (flat event-list any/all), and the
    RE-EXPRESSED type guards **`is_eseq`/`is_seqelist`/`is_stslist`**.
  - `events.py` promoted to an **`events/` package** (`_object` faithful /
    `_approx` approximate) with the public import surface unchanged.
- **RP5b — Event-mining engine, the R-faithful half RP5a deferred (TraMineR
  2.2-12 / TraMineRextras 0.6.8 goldens).** A pure-Python port of TraMineR's
  vendored ~2,345-line C++ event engine (`events/_engine.py`:
  `Sequence`/`SequenceEventNode` with the six count methods, the prefix-tree
  miner, `tmrfindsubsequences`/`tmrmatrixsubseqinseq`/`tmrseqedist`), with the
  R-faithful surface built on it (`events/_mining.py`), each landed with a live-R
  golden:
  - **`seqeconstraint`** — R fields `max.gap`/`window.size`/`age.*`/`count.method`
    ∈ 1..6 (COBJ/CDIST_O/CWIN/CMINWIN/CDIST), with R's validation.
  - **`seqefsub`** — the prefix-tree frequent-subsequence miner returning a
    faithful **`SubSeqElist`** (R `subseqelist`: generic `(A)-(B)` labels, a
    `data` frame of `Support`/`Count`); **`is_subseqelist`** type guard.
  - the per-sequence subsequence-count matrix (events-package public) over all
    **six count methods** + a `rules` mode.
  - **`seqerules`** — full statistics `Support`/`Conf`/`Lift`/`Standardlift`/
    `JMeasure`/`ImplicStat`/`p-value`/`p.valueB1`/`p.valueB2` (ported internal
    `implicativestat`; degenerate rules reproduce R's exact `NaN`/`±Inf`).
  - **`seqecmpgroup`** — discriminating-subsequence chi² (`SubSeqElistChisq` with
    per-group `Freq`/`Resid`, Yates' correction for 2×2, Bonferroni).
  - **`seqedist`** — event-sequence edit distance (per-event indel `idcost`,
    time-warp `vparam`, `interval` 1/2/3, `norm` 0/1/2).
  - **`seqecreate`** now honours a factor / pandas `Categorical` event column's
    level order as the event dictionary (load-bearing for the miner's
    comma-grouping and Count-tie ordering).
  - **`seqentrans`** re-pointed onto the faithful `SubSeqElist`.
- **RP6a — Indicators & advanced statistics, first half (sequence operations +
  VERIFY reconciliations; TraMineR 2.2-12 / TraMineRextras 0.6.8 goldens).** RP6
  was split 6a/6b; 6a lands the lower-risk, mostly-deterministic items, each its
  own commit + live-R golden:
  - **`seqconc`** — per-sequence string of non-void states joined by `sep`.
  - **`seqdecomp`** — split joined state strings into a padded matrix
    (`miss` → `None`).
  - **`seqdim`** — `[Rows, Columns]` (decompresses `"-"`/`":"`-joined strings).
  - **`seqmpos`** — number of matching positions between two sequences (void
    excluded always; `nr` unless `with_missing`).
  - **`seqcomp`** — whether two sequences are identical (non-void length + states).
  - **`seqsubsn`** — number of distinct subsequences per sequence (`nsubs` DP, DSS
    by default).
  - **`seqpm`** — regex pattern match over `seqconc` strings → `nbocc` + 1-based
    match indices.
  - **`seqrep_grp`** (R `seqrep.grp`) — per-group representatives with a shared
    `dmax = max(diss)`; per-group quality `Statistics`.
  - **`seqrecode`** / **`seqgen`** confirmed R-faithful (VERIFY; `seqgen` via a
    structure/seeded test — stochastic, no numeric golden).
- **RP6b — Indicators & advanced statistics, second half (TraMineR 2.2-12 /
  TraMineRextras 0.6.8 goldens).** The advanced statistics plus the faithful R
  `seqformat` STS→SPELL/SRS hub, each its own commit + live-R golden. Phase 6 is
  now complete.
  - **`STS_to_SPELL`** / **`STS_to_SRS`** — the faithful R `seqformat` spell /
    shifted-replication converters (closes RP6a's deferred T18). Gate-2 surprise:
    R's `seqformat` re-runs `seqdef(right="DEL")` internally and unifies `nr`+void
    into NA, so a trailing missing run is dropped while leading/interior `nr` is
    kept. The minimal hub needed by 6b was built; the rest (TSE/SPS/DSS round-trips)
    stays labelled-simplified.
  - **`seqidegrad`** — degradation index (`"Degrad"`); generalized the existing
    correction with the `spell.integr=TRUE` branch (the destination spell's
    integration weights, computed from `seqdur`).
  - **`seqinsecurity`** — insecurity index (type=2 of `seqprecarity.private`):
    `stprec[first]*integr1 + (ici + correction)`, optional `bound` clamp.
  - **`seqimplic`** — implicative statistics of a group over positions (a
    `levels × alphabet × positions` array); ports the TraMineRextras
    `implicativestat` (weighted, multi-level, continuity).
  - **`seqindic_dyn`** (R `seqindic.dyn`) — `seqindic` over a sliding /
    incremental window → the `dynin` matrix.
  - **`seqsha`** — sequence-history person-period converter.
  - **`seqgranularity`** — time-granularity coarsening → a new `STSList`
    (`first`/`last`/`first.valid`/`last.valid`/`mostfreq`).
  - **`seqpcfilter`** — the parallel-coordinate line-filter constructor
    (`minfreq`/`cumfreq`/`linear`); the `seqpcplot` plot is RP7.
  - **`seqauto`** — state autocorrelation by lag (SRS + Cramer's V).
  - **`seqcta`** — competing-trajectory-analysis converter (via SRS).
  - **`seqsurv`** — weighted Kaplan-Meier spell survival per state. BUILD-from-R,
    not the Sequenzo borrow: the Sequenzo KM guards `at_risk > d` and so yields 0
    std.err at `surv=0` where R's `survfit` yields NaN (reproduced R exactly).
  - **`seqpolyads`** — linked-polyad analysis; the deterministic `observed.dist` /
    `mean.dist['Obs']` are golden-pinned, the Monte-Carlo part is structure/seeded-
    tested (R's RNG is not reproducible from numpy).
  All are BUILD-from-R (TraMineR/TraMineRextras GPL source) — docstring provenance
  only, no BSD-3 / NOTICE. SAMM stays DEFER (D1); all visualization is RP7.

### Changed

- **`seqLLCS` / `seqLLCP` re-pointed to the R scalar contract (BREAKING, RP6a
  T12).** Live R's `seqLLCS(seq1, seq2)` / `seqLLCP(seq1, seq2)` return the
  *scalar* length of the longest common subsequence / common prefix between two
  sequences (compared over each sequence's non-void length, `nr` kept). The prior
  traminex functions returned an n×n *distance matrix* — a wrong-target name
  clash. The LCS/LCP distance matrices remain available via
  `seqdist(method="LCS"/"LCP")`. Pinned by `tests/test_r_goldens_seqLLCS.py`.
- **`seqalign` now aligns over the non-void length, with a faithful cost alphabet
  (RP6a T17 VERIFY fix).** R aligns over `seqlength()` positions, so trailing
  structural void is dropped before the DP (`nr` kept); and the substitution-cost
  alphabet is the real states (+`nr` iff `with.missing`), not the padded
  `x.alphabet`. Clean rectangular alignments are byte-unchanged. Pinned by
  `tests/test_r_goldens_seqalign.py`.
- **`seqformat` TSE/SPELL converters labelled as simplified, not R's hub (RP6a
  T18).** `seqformat` (wide↔long), `TSE_to_STS`/`STS_to_TSE`/`HSPELL_to_STS`/
  `FCE_to_TSE` are simplified "state-at-time" helpers, **not** R's event/tevent-
  driven `seqformat` STS/SPELL/SPS/DSS/SRS/TSE hub. The faithful hub (STS→SPELL/
  SRS, needed by `seqsurv`/`seqauto`/`seqcta`) is deferred to RP6b.
- **`seqelength` re-pointed to the R-faithful stored length (BREAKING).**
  TraMineR `seqelength` is the per-sequence stored observation time (`-1` for a
  TSE-created eseq; set via the `SeqElist.lengths` property = `seqelength<-`),
  **not** the number of events per id (the prior, non-R behaviour). A set length
  shows in `as_character` as a trailing gap (`"(A)-2-(B)-3"`). Pinned by
  `tests/test_r_goldens_eseq_length.py`.

- **The event miner `seqeconstraint` / `seqefsub` / `seqerules` is now R-faithful
  (BREAKING, RP5b — replaces the prior approximate miner).** `seqeconstraint` takes
  R's `max.gap`/`window.size`/`age.*`/`count.method` ∈ 1..6 (not the approximate
  `contiguous`/`min_gap`/`gap_in_time`/`count_method` ∈ {OCC, SEQS}); `seqefsub`
  returns a faithful `SubSeqElist` (not a flat DataFrame) and takes
  `min_support`/`pmin_support`/`max_k`/`weighted` (not `max_len`/`as_percent`/
  `id_weights`); `seqerules` takes a `SubSeqElist` and returns R's full rule
  statistics (not the approximate `antecedent`/`consequent`/support/conf/lift).
  Their approximate behaviour tests were rewritten/removed in favour of live-R
  goldens. The `seqerules`→`seqerulesdisc` rename remains dropped (traminex
  `seqerules` is R `seqerules`). The genuinely-non-R helpers
  (`events_transition_matrix`, `frequent_*`, `events_to_state_sequences`) are
  unchanged. See `docs/API_PARITY.md`.

- **`seqdomassoc` weighted path: documented + pinned the intentional divergence
  from R (RP3/RP4 audit).** The unweighted path is byte-faithful to R. R's
  `seqdomassoc(weighted=TRUE)` mis-assigns case weights (column-major `unlist`
  vs case-major weight-rep — a latent R bug); traminex keeps the sensible
  per-case weighting (consistent with `dissvar_grp`). The weighted golden now
  asserts traminex's result via R's own `loglin` on the correctly-weighted
  cross-table, and the divergence is documented in `docs/API_PARITY.md`. Pinned
  by `tests/test_r_goldens_seqdomassoc.py` (weighted + unweighted).

### Fixed

- **RP5 + RP6 independent audit-and-fix (live TraMineR 2.2-12 / TraMineRextras
  0.6.8 / survival 3.8.6).** An adversarial re-verification of the event framework
  (RP5a/RP5b) and the indicators/advanced-statistics (RP6a/RP6b) against new
  corners the shipped fixtures never exercised found six real port bugs, each fixed
  with its own live-R golden:
  - **`seqetm`** no longer double-appends `nr`/`void` to the alphabet (traminex's
    `x.alphabet` already folds them in) — it produced an N+2 matrix with a
    duplicated final row/col on any data carrying a missing or void state.
  - **`seqeid`** honors the id-column dtype: integer ids keep their literal value
    and factor (pandas Categorical) ids use their declared level codes, matching
    R's `if (is.character(id)) id <- factor(id); id <- as.integer(id)` (only
    character ids get a sorted-unique rank). Integer ids were being string-sorted.
  - **`seqeapplysub(method="age")`** truncates the first-occurrence age toward zero
    (C++ accumulates into an `int`), matching R (e.g. `1.9999 → 1`).
  - **`seqidegrad` / `seqinsecurity` / `seqprecarity`** weight the FREQ/TRATE/
    TRATEDSS transition counts by the sequence weights (R's `seqtrate` default is
    `weighted=TRUE`); the rate-based methods previously ignored weights.
  - **`STS_to_SPELL`** re-derives interior/leading `void` as `nr` (kept under
    `with_missing`), matching R's `seqformat` nr+void→NA unification; only the
    trailing missing run is dropped. Interior void was being silently dropped.
  - **`seqsurv`** weighted `std.err` is now the infinitesimal-jackknife robust SE
    (R sets `robust=TRUE` for fractional weights via `has.rwt`), reproduced to
    machine precision; integer weights keep Greenwood. The weighted golden pins
    `std.err` (was curve-only).

  Documented R-divergences ratified (R crashes/aborts on degenerate input where
  traminex returns the sensible value): a degenerate `CWIN` window, `seqsha`
  start-align `ma<L`, `seqpcfilter` unreachable `cumfreq` level. Findings:
  `docs/superpowers/research/2026-06-20-rp5-rp6-audit/findings.md`.

- **`seqdist(method="HAM")` is now R-faithful (RP5/RP6 audit, P1; BREAKING).** It
  returns the **raw** per-position substitution-cost sum by default (R `norm="none"`;
  unit cost per mismatch, or `sum_t sm[a,b]` with an explicit `sm`), matching R.
  Previously it was length-normalized (`mismatches / L`) with no opt-out and rejected
  the `norm=`/`sm=` arguments. The `/L` proportion is now available via
  `norm="maxlength"` (== `"auto"` for fixed-length data, dividing by `L * max(sm)`).
  This was the one cross-domain item the audit had only *proposed*; on review the
  blast radius was contained (no test pinned the `/L` value; pyweightedcluster computes
  HAM via the DHD-unit-sm path and normalizes itself, so it is unaffected), so it was
  fixed. `seqpolyads(method="HAM")` (R's default method) now matches R exactly. Pinned
  by `tests/test_r_goldens_distances.py` (HAM raw / `norm="maxlength"` / custom `sm`).

- **`seqdiff` `with_missing=True` is now faithful on ragged sequences (RP3/RP4
  audit).** `seqdiff` now passes the structural `void` code to the within-window
  `seqdist`, so void is deleted (not treated as a state) under `with_missing=True`,
  matching R's `with.missing=TRUE`. The shipped golden used a no-missing fixture;
  the void path was unexercised and diverged from R. Interior-missing (`nr`) was
  already faithful. Pinned by `tests/test_r_goldens_seqdiff_void.py`.

- **`seqdist` now deletes the structural `void` symbol on ragged sequences
  (RP3/RP4 audit).** TraMineR's `seqdist` removes void (the right padding of
  ragged sequences) before computing an order-based distance, keeping the
  missing state `nr` as a regular state (`with.missing=TRUE`); traminex's
  distance kernels previously treated `void` as an ordinary state, inflating
  common subsequences and diverging from R by up to the void-run length. When a
  state-sequence object's `void_code` is present, `seqdist`/`pairwise` now route
  the order-based methods (`LCS`/`OM`/`LCP`/`RLCP`, `norm=None`) through the
  per-pair kernels on the void-deleted codes. The dense fast-path is byte-
  unchanged when no void cells are present. This makes `seqdiff`,
  `seqCompare`/`seqLRT`/`seqBIC` faithful on ragged data with `with_missing=True`.
  R's `norm="auto"` (geometric-mean) on ragged sequences remains unsupported
  (not in the distance layer). Pinned by `tests/test_r_goldens_distances_void.py`.

### Added

- **`seqefsub(str_subseq=)` — specific (non-mined) subsequences (RP5/RP6 audit).**
  Ports R's `seqecreatesub` label parser so `seqefsub` can score given generic
  subsequence labels (e.g. `"(A)-(B,C)"`) instead of mining: `Count = Σ wᵢ·count`,
  `Support = (Σ wᵢ·COBJ-count)/totseq`, in the given order. Closes the prior
  `NotImplementedError`. Pinned by `tests/test_r_goldens_seqefsub.py`.

- **`seqrf` — relative-frequency representative sequences (roadmap Phase 4 —
  Dissimilarity & group, ★; ADAPT from Sequenzo).** New public symbol porting
  TraMineR 2.2-12 `seqrf`: resolves the `sortv` (`"mds"` / `"from.start"` /
  `"from.end"` / vector), runs `dissrf`, and returns the representative
  sub-sequences (`seqtoplot`) weighted by the group `heights`, plus the full
  `dissrf` result. This is the numeric backbone of the upcoming Phase-7
  `seqrfplot` (`plot.seqrf`). Adapted from
  `sequenzo/representative_sequences.py` (BSD-3, see NOTICE). Pinned by
  `tests/test_r_goldens_diss_rf.py`.

- **`seqdomassoc` — cross-domain association of multidomain sequences (roadmap
  Phase 4 — Dissimilarity & group; ADAPT from Sequenzo).** New public symbol
  porting TraMineR 2.2-12 `seqdomassoc` (`rep.method="overall"`): for each domain
  pair, cross-tabulates the position-wise states and reports the independence
  likelihood-ratio (`LRT`) and/or Cramer's `V` (with `loglin`'s structural-zero
  df adjustment). Adapted from
  `sequenzo/multidomain/association_between_domains.py` (BSD-3, see NOTICE).
  `rep.method` `"rep"` / `"eq.group"` raise `NotImplementedError`. Pinned by
  `tests/test_r_goldens_seqdomassoc.py`.

- **`dissmergegroups` — greedy ASW group merging (roadmap Phase 4 —
  Dissimilarity & group; ADAPT from Sequenzo).** New public symbol porting
  TraMineR 2.2-12 `dissmergegroups`: iteratively merges the group pair whose
  merge least degrades the partition's average silhouette width (ASW), stopping
  at `min_group` / a `crit` quality-loss threshold. The greedy-loop shape is
  adapted from `sequenzo/.../merge_cluster_groups.py` (BSD-3, see NOTICE), but
  the **ASW follows `WeightedCluster::wcClusterQuality` / R directly** —
  Sequenzo's own `_asw` is *not* faithful to WeightedCluster (wrong a(i)/b(i)
  denominators) and was not used. Verified byte-for-byte against live
  WeightedCluster 2.0; pinned by `tests/test_r_goldens_diss_mfa.py`.

- **`dissindic` — individual diss-based indicators (roadmap Phase 4 —
  Dissimilarity & group; ADAPT from Sequenzo).** New public symbol porting
  TraMineRextras 0.6.8 `dissindic`: per-observation `marginality` (residual to
  the group model) and `gain` (discrepancy explained by the group). Adapted from
  `sequenzo/discrepancy_analysis/stats/individual_indicators.py` (BSD-3, see
  NOTICE); reuses traminex's `gower_matrix` + weighted hat matrix. (R's wrapper
  errors for a single grouping variable — a `data[seqok, ]` one-column-drop bug
  — so the golden replicates the documented `dissmfacwresid` math with live R's
  own `gower_matrix` + `qr`.) Pinned by `tests/test_r_goldens_diss_mfa.py`.

- **`seqCompare` / `seqLRT` / `seqBIC` — LRT / BIC group comparison from
  sequences (roadmap Phase 4 — Dissimilarity & group, B-7).** New public symbols
  porting TraMineRextras 0.6.8 `seqCompare`/`seqLRT`/`seqBIC`: compute the
  sequence `seqdist` and run the same LRT/BIC two-group split as `dissCompare`.
  Only the deterministic `s=0`, `group`-based, single-set path is supported.
  Pinned by `tests/test_r_goldens_group_compare.py` vs live R.

- **`dissCompare` — LRT / BIC group comparison from dissimilarities (roadmap
  Phase 4 — Dissimilarity & group).** New public symbol porting TraMineRextras
  0.6.8 `dissCompare`: the likelihood-ratio test, ΔBIC and Bayes factor of a
  two-group split of a dissimilarity matrix. Only the deterministic `s=0`
  (no-sampling) path is exposed. Pinned by
  `tests/test_r_goldens_group_compare.py` vs live R.

- **`dissvar_grp` — per-group discrepancy (roadmap Phase 4 — Dissimilarity &
  group).** New public symbol porting TraMineRextras 0.6.8 `dissvar.grp`: the
  weighted discrepancy `dissvar` within each group level. Pinned by
  `tests/test_r_goldens_diss_mfa.py` vs live R.

- **`gower_matrix` gains R's `squared` / `weights` parameters (roadmap Phase 4 —
  Dissimilarity & group, T15).** Now a faithful port of TraMineR 2.2-12
  `gower_matrix(diss, squared=TRUE, weights=NULL)` — `(I - 1 s^T)(D'/-2)(I - s
  1^T)` with `s = w/sum(w)` and `D' = D^2` when `squared`. The defaults
  (`squared=True`, unweighted) keep the classical `-0.5 J D^2 J` MDS Gram matrix,
  so existing callers are unaffected. Pinned by `tests/test_r_goldens_gower.py`.

- **Core-stats live-R golden coverage (roadmap Phase 3 — Core-stats VERIFY,
  T3-T11).** Pins the core statistics to TraMineR 2.2-12 across the weighted /
  `with.missing` / ragged-void paths (shared A/B/C fixtures in
  `scripts/r_goldens/make_goldens_corestats.R` + `tests/corestats_fixtures.py`):
  `seqistatd` (per-sequence state-duration counts, `tests/test_r_goldens_seqistatd.py`);
  `seqmeant` (mean time per state incl. Var/Stdev/SE and `prop`, `tests/test_r_goldens_seqmeant.py`);
  `seqmodst` (modal-state sequence + frequencies + occurrences, `tests/test_r_goldens_seqmodst.py`);
  `seqtab` / `seqtabstocc` (weighted unique-sequence frequency table; state
  co-occurrence frequencies, `tests/test_r_goldens_seqtab.py`); `seqtrate`
  (public transition-rate matrix, `tests/test_r_goldens_costs_trate.py`).

- **Tree visualization design layer — `disstree2dotp` / `seqtree2dot` /
  `disstreedisplay` / `seqtreedisplay` (roadmap Phase 2 — Trees, B-12; decision
  D4).** `disstree2dotp` writes `disstree2dot` plus a global-fit *quality* node
  (Pseudo F / Pseudo R2 / Levene with R significance stars); `seqtree2dot` adds an
  alphabet *legend* node for sequence trees; `disstreedisplay` / `seqtreedisplay`
  build the `.dot` and render it to an image by shelling to GraphViz `dot`, an
  **optional** extra (`pip install 'traminex[tree]'` + a system GraphViz) with
  graceful fallback — when `dot` is absent the `.dot` is still written and a clear,
  actionable error is raised, and the test suite passes with GraphViz not
  installed. Pinned by **structure tests**, not numeric goldens
  (`tests/test_disstree_display.py`). Documented divergence: TraMineR renders a
  jpeg quality/legend panel and a per-node sequence plot in every node; that rich
  per-node imagery is deferred to the Phase-7 viz layer, so here the quality/legend
  panels are text DOT nodes and the per-node images are omitted (topology, fit
  statistics and split structure are faithful). Adapts the optional-graphviz /
  graceful-fallback mechanism of Sequenzo's
  `discrepancy_analysis/trees/tree_visualization.py` (BSD-3, see `NOTICE`); the DOT
  bytes follow TraMineR.

- **Faithful `disstree2dot` — GraphViz DOT export of the R-faithful tree (roadmap
  Phase 2 — Trees, T20).** `disstree2dot(tree, filename, use_title=, show_depth=,
  digits=)` now ports TraMineR 2.2-12 `disstree2dot` over the fitted `disstree`
  object: it writes a `<filename>.dot` whose node names, `n`/`s2`/`Split`/`R2`
  labels, the `show_depth` rank-schedule block and edges byte-match R's text path
  (`image.fun=NULL`). Pinned by a **structure golden** (decision D4), not a numeric
  one (`tests/r_goldens/disstree_dot_*.dot`; generator
  `scripts/r_goldens/make_goldens_disstree_dot.R`, test
  `tests/test_r_goldens_disstree_dot.py`). The previous public `disstree2dot` was
  the sklearn/MDS *exploratory* exporter — it is retained, renamed
  `disstree2dot_explore`. Per-node image rendering is not ported (text/topology
  only; render the `.dot` with graphviz `dot`).

- **`disstreeleaf` / `disstree_get_rules` / `disstree_assign` — tree leaf helpers
  (roadmap Phase 2 — Trees, B-5).** Three TraMineR functions over a fitted
  `disstree`: `disstreeleaf(tree, label=True, collapse=)` returns each case's
  root→leaf path label (e.g. `"X > 3.1 & A [b; c] & A [b]"`); `disstree_get_rules`
  returns the per-leaf classification rules as **R-syntax** condition strings
  (`A %in% c('b', 'c') & X <= 3.1`), byte-identical to R's `disstree.get.rules`;
  `disstree_assign(rules, profile)` assigns covariate profiles to leaves, returning
  the 1-based leaf index per row. Pinned to live TraMineR 2.2-12 goldens
  (`tests/r_goldens/disstree_full_*`; generator
  `scripts/r_goldens/make_goldens_disstree_full.R`, test
  `tests/test_r_goldens_disstree_full.py`). Adapts the recursion shape of Sequenzo's
  `discrepancy_analysis/trees/tree_leaf_helpers.py` (BSD-3, see `NOTICE`); the label /
  rule strings follow TraMineR directly (Sequenzo emits Python-syntax conditions).

- **`seqemlt` — EMLT analysis of time-stamped situations (roadmap Phase 1, B-2).**
  New `seqemlt(seqdata, a=1.0, b=1.0, weighted=True)` ports TraMineRextras
  `seqemlt`: per-situation (time-stamped state) future transition profiles
  (time-discounted by `1/(a*lag + b)`), squared-Euclidean profile distances,
  their Benzécri double-centred covariance, principal-component sequence
  coordinates (via a `princomp(cor=TRUE)` replication), and the situation
  correlation matrix. The sign-stable quantities (correlation, profile distances,
  component standard deviations) are pinned to live TraMineRextras 0.6.8 goldens
  (`tests/r_goldens/seqemlt_*.csv`; generator
  `scripts/r_goldens/make_goldens_seqemlt.R`, test
  `tests/test_r_goldens_seqemlt.py`); the sign-/rotation-arbitrary coordinates are
  pinned via the invariant Gram matrix `coord @ coord.T` (plan decision 0.2.3).
  Adapts the structure of Sequenzo's `emlt/seqemlt.py` (BSD-3, see `NOTICE`),
  re-verified against live R.

- **`dissmfacw` — weighted multi-factor analysis of a distance matrix (roadmap
  Phase 1, B-1; the family unlock).** New `dissmfacw(diss, factors, weights=None,
  nperm=1000, seed=None, gower=False, squared=False)` ports TraMineR
  `dissmfacw`: a distance-based weighted pseudo-ANOVA (McArdle & Anderson 2001) —
  weighted Gower double-centring, QR-based weighted hat matrices, and Type-III
  "backward" sums of squares — returning per-term and `Total` PseudoF, PseudoR2,
  explained SS, and permutation p-values. The deterministic decomposition is
  pinned to a live TraMineR 2.2-12 golden (`tests/r_goldens/dissmfacw_det.csv`;
  generator `scripts/r_goldens/make_goldens_diss_mfa.R`, test
  `tests/test_r_goldens_diss_mfa.py`); the stochastic p-values are verified
  separately with a fixed seed (plan decision 0.2.2). **This unblocks the
  downstream pyweightedcluster `clustassoc` rewrite (sibling follow-on).** Adapts
  the structure of Sequenzo's `discrepancy_analysis/stats/multifactor_association.py`
  (BSD-3, see `NOTICE`), re-verified against live R.

### Changed

- **BREAKING: `dissrf` is now the R-faithful weighted relative-frequency grouping
  (roadmap Phase 4 — Dissimilarity & group, T13).** Replaces the simplified
  gmedoid-sort stand-in with TraMineR 2.2-12 `dissrf` for `grp.meth="prop"`
  (Gilbert Ritschard's weight-aware grouping into `k` equal-weight groups with
  fractional boundary sequences). The return dict now exposes R's fields
  (`medoids`, `R2`, `Fstat`, `pvalue`, `heights`, `index_list`, `weights_list`,
  `dist_list`, `sizes`, …); the old `group_assignment`/`F_stat`/`p_value` keys are
  gone. New params `sortv` / `grp_meth` / `squared` / `pow`; only `grp_meth="prop"`
  is implemented. The impl moved to `traminex/repseq_rf.py` (ADAPT from Sequenzo,
  BSD-3, see NOTICE) and is re-exported from `traminex.align_p2`. Pinned by
  `tests/test_r_goldens_diss_rf.py`.

- **BREAKING: `seqrep` now extracts representative sequences via `dissrep`, not
  an ad-hoc seqtab/medoid heuristic (roadmap Phase 4 — Dissimilarity & group,
  T13).** Following TraMineR 2.2-12 `seqrep`, it builds (or accepts via `diss=`)
  the pairwise dissimilarity matrix and delegates to `dissrep`, returning a dict
  with `rep_indices`, the `representatives` code rows, and the dissrep
  `statistics` / `quality` / `distances` / `rep_group`. Defaults follow R
  (`coverage=0.25`, `pradius=0.10`); the `k`/`sep`/`coverage`-key fields are
  gone (use `nrep` for a fixed count). `criterion="prob"/"mscore"` are not yet
  implemented. Pinned by `tests/test_r_goldens_dissrep.py`.

- **BREAKING: `seqdiff` is now the R-faithful sliding-window discrepancy
  analysis, not per-position chi²/residuals (roadmap Phase 4 — Dissimilarity &
  group, T16).** Following TraMineR 2.2-12 `seqdiff`, each window position runs
  `dissassoc` on the within-window `seqdist` and the result is a dict
  `{"stat", "discrepancy"}` of tables — `stat` columns Pseudo F / Pseudo Fbf /
  Pseudo R2 / Bartlett / Levene, `discrepancy` columns the group levels plus
  Total. New params `cmprange` / `method` / `norm` / `with_missing` / `weighted`
  / `squared`; the old `method="residuals"/"prop_diff"/"chi2"/"log_odds"` /
  `bootstrap` API is removed. `norm` defaults to `None` (unnormalised LCS), since
  R's `norm="auto"` (gmean) is not yet in the distance layer. Pinned by
  `tests/test_r_goldens_seqdiff.py`.

- **BREAKING: `dissdomassoc` is now the R-faithful inter-domain association, not
  a per-domain variance summary (roadmap Phase 4 — Dissimilarity & group, T14).**
  Following TraMineR 2.2-12 `dissdomassoc`, it returns weighted Pearson /
  Spearman correlation matrices between the domains' dissimilarity vectors (with
  the slope t-test p-values from the vendored `wtd.cor`), their R2, and
  Cronbach's alpha — selectable via `what=`, with `weights=` / `dnames=` /
  `jointdiss=`. The old `{domain_i_var, joint_overall_var}` summary is removed.
  Pinned by `tests/test_r_goldens_dissdomassoc.py`.

- **BREAKING: `seqindic` is now the composite per-sequence indicator table, not
  a one-hot state-indicator matrix (roadmap Phase 3 — Core-stats VERIFY, T11).**
  R's `seqindic(indic=, with.missing=)` returns one column per indicator; the
  default `c("visited","trans","entr","cplx","turb2n")` is reproduced, along with
  the basic/diversity/complexity groups (`lgth, nonm, dlgth, visited, visitp,
  recu, trans, transp, entr, cplx, turb, turbn, turb2, turb2n`). The binary /
  ranked indicators (`ppos/integr/prec/...`) and `meand/dustd/nsubs` are deferred
  (raise `NotImplementedError`). Signature `seqindic(x, *, indic=..., with_missing=,
  w=)`; the old `by`/`as_dataframe` parameters are gone. Pinned by
  `tests/test_r_goldens_seqindic.py` (generator
  `scripts/r_goldens/make_goldens_seqindic.R`).

- **BREAKING: `seqtransn` is now the per-sequence transition count, not the
  dataset K×K transition matrix (roadmap Phase 3 — Core-stats VERIFY, T4).** R's
  `seqtransn` returns `length(seqdss) - 1` for each sequence (number of state
  changes), with `with_missing`, `norm` (divide by the non-void length − 1) and
  `pweight` (weight each transition by `1 − rate`). traminex previously returned
  the dataset-level transition *matrix*; that matrix is now the private
  `_transition_count_matrix` (the transition heatmap uses it / `seqtrate`). For
  transition *rates* use `seqtrate`. Pinned by `tests/test_r_goldens_seqtransn.py`
  (generator `scripts/r_goldens/make_goldens_seqtransn.R`).

- **BREAKING: `seqnum` is now the numeric recoding of TraMineR `seqnum`, and
  `seqasnum` is newly exposed (roadmap Phase 3 — Core-stats VERIFY, T11).**
  `seqnum` previously returned a per-sequence count of distinct states; it now
  recodes the sequences over a numeric `0..K-1` alphabet (R `seqnum`). The new
  `seqasnum` (R `seqasnum`) returns the numeric code matrix (`NaN` for the
  structural `void` always, and the missing code `nr` unless `with_missing`).
  Pinned by `tests/test_r_goldens_seqasnum.py` (generator
  `scripts/r_goldens/make_goldens_seqasnum.R`).

- **BREAKING: `seqstatl` now returns the list of distinct states present, not a
  last-position frequency table (roadmap Phase 3 — Core-stats VERIFY, T10).** R's
  `seqstatl` is `sort(unique(as.vector(seqdata)))` — the sorted observed states
  (numeric-aware; including the missing `nr` and structural `void` symbols when
  present). traminex previously returned a last-position `Freq`/`Percent`
  DataFrame (wrong target). Now returns a `list[str]`; the `weighted` /
  `include_missing` parameters are dropped (R has none). Pinned by
  `tests/test_r_goldens_seqstatl.py` (generator
  `scripts/r_goldens/make_goldens_seqstatl.R`).

- **BREAKING: `seqstatf` is now the overall state frequency, not the
  first-position distribution (roadmap Phase 3 — Core-stats VERIFY, T10).** R's
  `seqstatf` returns `colSums(seqistatd(seqdata, with.missing) * weights)` — the
  weighted total time spent in each state across *all* positions. traminex
  previously counted only the first position. The result is now indexed by
  `statl` (the alphabet, plus the missing code `nr` when `include_missing`;
  never the structural `void`), so the `nr`/`void` rows are no longer emitted
  unless requested. Pinned by a live-R golden,
  `tests/test_r_goldens_seqstatf.py` (generator
  `scripts/r_goldens/make_goldens_seqstatf.R`).

- **`seqprecarity` now supports a custom `state_order` (full alphabet permutation)
  — roadmap Phase 2 — Trees, item 4.** Previously `state_order`/`state_equiv` raised
  `NotImplementedError`; a full-permutation `state_order` now drives both the start
  cost (`seqprecstart`) and the degradation correction's penalty signs (R
  `seqidegrad.private` rank reordering), pinned to live TraMineR 2.2-12 goldens
  (`tests/r_goldens/seqprecarity_A_ord_*`; generator
  `make_goldens_precarity_order.R`). `state_equiv` (equivalence classes) and a
  *partial* `state_order` (non-comparable states) still raise a clear
  `NotImplementedError` (they need R's transition-collapse path). The
  `disstree(weight_permutation=…)` guard now raises `ValueError` for an unknown
  value (input validation; all five R `weight.permutation` values are ported — not
  a gap).

- **BREAKING — `disstree2dot` now means the R-faithful tree exporter** (port of
  TraMineR `disstree2dot` over the fitted `disstree`; see Added). The previous
  sklearn/MDS exploratory exporter is renamed **`disstree2dot_explore`**; callers
  that passed a `disstree_fit` result to `disstree2dot` must switch to
  `disstree2dot_explore`.

- **BREAKING — `seqentrans` is now the faithful TraMineRextras subsequence-table
  augmentation (roadmap Phase 1, T1).** `seqentrans(fsubseq, avg_occ=False)` no
  longer returns an event-label transition matrix; matching
  `TraMineRextras::seqentrans`, it now takes a `seqefsub` result and returns it
  augmented with `ntrans` (number of `-`-separated transitions in each
  subsequence label), `nevent` (`ntrans - 1 + ` number of `,`-separated tokens),
  and — when `avg_occ=True` — `avg_occ` (`Count / sum(seqeweight(eseq))`, R's
  `Avg.occ`). `seqefsub` now carries the total event-sequence weight in
  `result.attrs["seqe_weight_total"]`. The event-label adjacency matrix that
  `seqentrans` previously returned is now reached **only** via
  `events_transition_matrix` (unchanged). Pinned to live TraMineRextras 0.6.8
  goldens on TraMineR's own subsequence labels (`tests/r_goldens/seqentrans*.*`;
  generator `scripts/r_goldens/make_goldens_seqentrans.R`, test
  `tests/test_r_goldens_events_full.py`).

- **BREAKING — `seqMD` is now the faithful TraMineR CAT multidomain engine; new
  `seqdistmc` (roadmap Phase 1, D2).** `seqMD(channels, ..., what=)` no longer
  returns a lightweight `SeqMD` container; it now returns the TraMineR
  "additive-trick" (CAT) result — the combined-state multidomain sequence
  (`what="MDseq"`, default), the CAT substitution/indel cost structure
  (`what="cost"`), or the pairwise distance matrix (`what="diss"`). The new
  `seqdistmc(channels, ...)` is the `what="diss"` alias (TraMineR
  `seqdistmc`/`seqMD`). Both are pinned to live TraMineR 2.2-12 goldens
  (`tests/r_goldens/seqMD_*.csv`; generator
  `scripts/r_goldens/make_goldens_seqMD.R`, test `tests/test_r_goldens_seqMD.py`).
  The previous **additive** per-channel combination is preserved as the named
  **non-R extension** `seqMD_distances` (unchanged behaviour;
  `D = sum_c w_c * seqdist(channel_c)`), and the `SeqMD` container dataclass is
  retained for it and the `seqMD_gridsearch_weights` / `seqMD_optimize_weights`
  weight-search utilities (now documented as building on `seqMD_distances`, not
  the CAT path). Construct the channel container directly with `SeqMD([ch1, ch2])`
  (the old `seqMD([ch1, ch2], names=...)` builder is removed). The CAT engine
  adapts the structure of Sequenzo's `multidomain/cat.py` (BSD-3, see `NOTICE`),
  re-verified against live R.

- **Version single-sourcing (roadmap Phase 3, release-integrity).**
  `traminex.__version__` is now read from the installed package metadata
  (`importlib.metadata.version("traminex")`, single-sourced from
  `pyproject.toml` `[project].version`), with a graceful fallback (parse
  `pyproject.toml`, then `"0+unknown"`) for an uninstalled source tree; the
  hand-maintained literal in `traminex/__init__.py` is removed. A new guard test
  `tests/test_version_single_source.py` asserts all readers agree
  (`__version__` == installed metadata == `pyproject.toml` version), and
  `scripts/check_docs_consistency.py` now (a) verifies the version is
  single-sourced and (b) enforces that the README / CONTRIBUTING / `docs/INDEX.md`
  current-version strings match `pyproject.toml`. Stale `0.0.3` version strings in
  those docs are corrected to `0.0.8`, and `docs/API_PARITY.md` is reconciled to
  the current 136-name public API. No runtime/API behaviour change; closes the
  closeout audit's **H1** version-drift finding. Version stays **0.0.8**
  (changes remain under `[Unreleased]`).

- **BREAKING — eight complexity/precarity indicators are now R-faithful**
  (roadmap Phase 2). `seqintegr`, `seqprecstart`, `seqibad`, `seqivardur`,
  `seqmaintokens`, `seqipos`, `seqlogp`, and `seqprecarity` are rewritten to match
  TraMineR 2.2-12 exactly and pinned to live-R goldens
  (`tests/r_goldens/seq{integr,precstart,ibad,ivardur,maintokens,ipos,logp,precarity}_*.csv`;
  generator `scripts/r_goldens/make_goldens_precarity.R`, test
  `tests/test_r_goldens_precarity.py`). Numerical output and signatures differ from
  the previous proxies, e.g.: `seqlogp` is now a **non-negative Markov surprisal**
  (`-log P` under time-varying transition rates) and raises when missing states are
  present with `include_missing=False`; `seqprecstart` returns a **per-state cost
  vector** over the alphabet (was per-sequence); `seqivardur` returns the
  **per-sequence variance of spell durations** (was per-state); `seqipos` defaults
  to the **DSS positive-spell share**; `seqintegr` is the position-weighted state
  presence; `seqmaintokens` selects states by mean time (`seqmeant`). The seven
  newly-public names join `api_manifest` (`seqibad` was already exported). The old
  proxies are retained privately as `*_proxy` (e.g. `seqlogp_proxy`), not exported.
  `seqprecarity` ports the default `method="TRATEDSS"` path (type=1) plus a private
  `_seqidegrad_correction`; custom `state_order`/`state_equiv` raise
  `NotImplementedError`.

### Fixed

- **`dissrep` quality statistics are now R-faithful (roadmap Phase 4 —
  Dissimilarity & group, T13).** The representative *selection* already matched
  TraMineR, but the quality table used an unweighted row-mean for the distance to
  the global centre (R uses `disscenter`), an unweighted `na` count (R weights
  it), and `V = DC / group_weight` (R uses the within-group
  `mean(disscenter(diss[sel, sel], weights))`). Candidate ordering now also uses a
  stable sort matching R's `order(score, decreasing=)` tie-break. Pinned by
  `tests/test_r_goldens_dissrep.py` vs live TraMineR 2.2-12.

- **`seqient` (and therefore `seqici`) no longer count structural `void` in the
  longitudinal entropy (roadmap Phase 3 — Core-stats VERIFY, surfaced by T11
  `seqindic`).** With `include_missing=True` on ragged/void-padded data the
  per-sequence state distribution included the structural `void` column, so the
  entropy (and the complexity index built on it) was inflated — e.g. a constant
  void-padded sequence reported nonzero entropy. The distribution is now taken
  over `statl` only (real alphabet + nr when `with.missing`), matching the
  function's own normalizer and TraMineR 2.2-12.

- **Public `seqtrate` is now R-faithful for `with.missing` and gains weighting
  (roadmap Phase 3 — Core-stats VERIFY, T9).** `seqtrate(include_missing=True)`
  previously added the structural `void` symbol to the state set on ragged data;
  R's `seqtrate` only adds the missing code (`nr`) when `with.missing=TRUE` and
  never treats `void` as a state. A `weighted=True` parameter (R's default) was
  added so the rates honour sequence weights. Empty source rows stay zero (R's
  `PA == 0`). Pinned by a live-R golden across the weighted / `with.missing` /
  ragged paths in `tests/test_r_goldens_costs_trate.py` (generator
  `scripts/r_goldens/make_goldens_seqtrate.R`).

- **`seqmeant(prop=True)` now averages per-sequence proportions R-faithfully
  (roadmap Phase 3 — Core-stats VERIFY, T6).** R normalizes each sequence's
  durations by its own total over `statl` (valid length) *before* the weighted
  average; traminex previously divided the mean by the rectangular sequence
  length `L`, biasing every value on ragged / `with.missing` data. The
  `prop=False` Mean/Var/Stdev/SE path was already faithful. Pinned by
  `tests/test_r_goldens_seqmeant.py`.

- **`seqstatd` no longer counts structural `void` as a state (roadmap Phase 3 —
  Core-stats VERIFY, T3).** With `include_missing=True` on ragged/void-padded
  data the per-position frequencies, `ValidStates` and normalized `Entropy`
  previously divided by a denominator that included the structural `void`
  symbol, biasing every value. R's `seqstatd` adds only the missing (`nr`) code
  to `statl` when `with.missing=TRUE` and never treats `void` as a state; the
  entropy normalizer is the count of `statl` columns. Now matches TraMineR
  2.2-12 across the weighted / `with.missing` / ragged paths, pinned by
  `tests/test_r_goldens_seqstatd.py` (generator
  `scripts/r_goldens/make_goldens_seqstatd.R`).

- **Precarity/complexity indicators are now R-faithful in the void / empty /
  missing-alphabet corners** (Phase 2 executed-audit follow-up; pinned by new R
  goldens in `tests/r_goldens/seq*_{Kr,J,A_m,V,W}*.csv`). Each was confirmed
  against live TraMineR 2.2-12:
  - **`seqlogp` on void-padded (ragged) sequences.** The time-varying transition
    rates and first-position frequencies counted the structural `void` symbol as a
    state, so void padding contaminated the *shared* rate array and produced wrong
    surprisals — even for full-length rows and on the default `include_missing=False`
    path. `_timevarying_trate` gained an additive `ignore_code` (DHD callers
    unchanged) and `seqlogp` now excludes `void` from both the rates and the begin
    frequencies. (audit F1)
  - **`seqmaintokens(mint=…, k=…)` when more than `k` states pass the threshold.**
    Now reproduces TraMineR's exact result — R returns positions *into the
    thresholded subset* (an upstream `order()` quirk), not alphabet indices. (audit F2)
  - **`include_missing=True` on data with no missing values.** `seqintegr`,
    `seqprecstart`, `seqibad`, `seqprecarity` (and `seqient`'s normalization) now
    append the phantom `nr` state that R's `alphabet(with.missing)` always carries,
    so the effective alphabet size is `k_real + 1` (consistent with `seqivardur` /
    `seqivolatility` / `seqST`). (audit F3)
  - **Fully-empty (all-`void`) rows now return `NaN`** (not `0`) from `seqintegr`,
    `seqibad`, `seqivardur`, and `seqprecarity`, matching R's division-by-zero. (audit F4)
  - **`seqici` no longer counts `void` transitions.** Its normalized transition term
    is now `(DSS length − 1)/(non-void seqlength − 1)` over void-stripped rows
    (= R `seqtransn(norm=TRUE)`), so `seqici` and `seqprecarity` are correct on
    ragged data. (audit F5)
- **`seqST` / `turbulence_exact` no longer count structural `void` as a state.**
  Under `include_missing=True` with void-padded (ragged) input, `seqST` previously
  treated right-padding `void` as a visited state (and the `norm=True` reference
  alphabet leaked `void`/`nr` and used `ncol` instead of the max non-void length),
  diverging from TraMineR. `seqST` now strips `void` per row and uses the
  effective alphabet (`k_real + include_missing`, `void` always excluded) for both
  the per-sequence turbulence and the normalization reference — matching R for
  rectangular **and** ragged inputs. Pinned to R goldens
  (`tests/r_goldens/seqST_{A,B,C}.csv`, generator
  `scripts/r_goldens/make_goldens_seqst.R`, test `tests/test_r_goldens_seqst.py`).
  Rectangular void-free results are unchanged.

- **`seqST` `type=2` turbulence now matches R's `seqivardur` type=2.** The
  duration-variance and its maximum used an alphabet-size (`K`) denominator and a
  type-1-style `var.max`, and fell back to the type-1 formula whenever the number
  of spells reached the alphabet size — diverging from TraMineR even for
  rectangular void-free input. `_dur_var_and_max` type=2 now uses R's denominator
  `d + nnvisit` (with `nnvisit = K − #visited states`) and R's `maxnnv`-based
  `var.max`. Pinned to R goldens (`tests/r_goldens/seqST_{A,B,C}_t2.csv`). `type=1`
  (the default) is unchanged.

### Added

- **Core representation transforms `seqdss`, `seqdur`, `seqivolatility`**
  (roadmap Phase 1). R-faithful ports of TraMineR 2.2-12 `seqdss` (distinct
  successive states, returning an `STSList`), `seqdur` (spell durations,
  NaN-padded), and `seqivolatility` (within-sequence volatility,
  `w*visit + (1-w)*transitions`, `adjust`/`include_missing` honoured). Faithful
  for rectangular **and ragged (void-padded)** sequences. Pinned to R goldens for
  three fixtures incl. a ragged one (`tests/r_goldens/seqdss_*.csv`,
  `seqdur_*.csv`, `seqivol_*.csv`; generator
  `scripts/r_goldens/make_goldens_seqdss_seqdur.R`,
  test `tests/test_r_goldens_seqdss_seqdur.py`). Exported via `api_manifest`.

## [0.0.8] - 2026-06-13

### Fixed

- **`seqsubm(method="TRATE")` is now bit-faithful to R `TraMineR::seqcost(method="TRATE")`**
  (Audit Phase B). Two divergences from the R source were corrected in
  `traminex.distances.costs.trate_substitution_matrix` (and its time-varying
  sibling):
  - The default additive `smoothing` was `1e-9`; R's `seqtrate` uses **no**
    smoothing. The default is now `0.0`, so the cost matrix matches R to
    floating-point precision. The `smoothing` parameter is retained for callers
    who want Laplace-style regularization (the OM-TRATE convenience path keeps
    its small explicit default).
  - A **source state with no observed transitions** (R's `if (PA == 0)
    tmat[x, ] <- 0`) was being turned into a uniform `1/K` rate row by the
    smoothing, biasing every cost involving that state by up to `1/K` (~0.33 for
    K=3). Empty source rows are now zeroed exactly as R does, independent of the
    `smoothing` value. The OM-TRATE convenience path (which passes its own
    `trate_smoothing`) benefits from the same correction.
  Pinned inside the oracle by `tests/test_r_goldens_costs_trate.py`, which
  replicates R's exact `seqcost`/`seqtrate` TRATE formula (weighted and
  unweighted, with an absent-source fixture) and checks CONSTANT, with no live-R
  dependency. The downstream pyweightedcluster A3 parity test (live R, atol 1e-8)
  is unaffected and now agrees to machine precision.

## [0.0.7] - 2026-06-13

### Added

- **R-faithful weighted regression trees `disstree` and `seqtree`**
  (migration plan Phase 2.4 / DUPLICATES A8). `traminex.disstree.disstree(
  diss, predictors, weights=, min_size=, max_depth=, R=, pval=,
  weight_permutation=, squared=, first=, random_state=)` ports TraMineR
  2.2-12 `disstree.R`: recursive binary partitioning of predictor levels
  minimising residual weighted inertia (ordered/numeric thresholds,
  categorical level subsets, NA group handling), the permutation gate
  (`DTNdissassocweighted`; variants none/diss/replicate/group), node
  medoids, the R2-based split schedule, leaf memberships
  (`disstreeleaf`) and the `dissassoc` adjustment. `traminex.seqtree.
  seqtree(seqdata, predictors, diss=None, weighted=True, ...)` ports
  `seqtree.R`, reading per-sequence weights from `seqdef(weights=)` like
  R's `attr(seqdata, "weights")`. Anchored to live R (node tables, fitted
  leaves, adjustment t0; weighted and unweighted) by
  `tests/test_r_goldens_disstree.py` (goldens
  `tests/r_goldens/disstree_*.csv` / `seqtree_*.csv`, generator
  `scripts/r_goldens/make_goldens_disstree.R`).

### Changed

- **The sklearn tree stand-ins are relabelled as exploratory** and moved to
  `traminex.disstree_explore` / `traminex.seqtree_explore` (`disstree_fit`,
  `seqtree_fit`, and their helpers; docstrings now state they are NOT the
  TraMineR trees). All historical names keep resolving — top-level imports
  are unchanged and `traminex.disstree` / `traminex.seqtree` re-export them.

## [0.0.6] - 2026-06-13

### Changed

- **`disscenter` now has the R TraMineR meaning; the name clash is resolved**
  (migration plan Phase 2.3 / DUPLICATES A7). `disscenter(diss, group=None,
  medoids_index=None, allcenter=False, weights=None, squared=False)` ports
  TraMineR 2.2-12 `disscenter`: the (weighted) distance of each case to its
  group centre, with `allcenter=True` (distance of every case to each group
  centre) and `medoids_index="first"/"all"` (group medoids, 0-based indices)
  supported. Anchored to live R by `tests/test_r_goldens_disscenter.py`
  (goldens `tests/r_goldens/disscenter_*.csv`, generator
  `scripts/r_goldens/make_goldens_disscenter.R`).
  **Breaking rename:** the previous traminex `disscenter` was an unrelated
  matrix double-centering helper; it is now `double_center_matrix`, with the
  deprecated alias `disscenter_double_center` (emits DeprecationWarning) as
  a greppable migration target. Family-wide, `disscenter` is now unambiguous
  (ECOSYSTEM.md name-clash warning closed).

## [0.0.5] - 2026-06-13

### Changed

- **`dissassoc` is now the R-faithful weighted discrepancy analysis**
  (migration plan Phase 2.2 / DUPLICATES A6). `dissassoc(diss, group,
  weights=None, R=1000, weight_permutation="replicate", squared=False,
  random_state=None)` ports TraMineR 2.2-12 `dissassoc` /
  `dissassocweighted`: weighted pseudo-ANOVA table, per-group discrepancies,
  and the five test statistics (Pseudo F, Pseudo Fbf, Pseudo R2, Bartlett,
  Levene) with R's permutation test (`pval = mean(t >= t0)` including the
  observed permutation). Supported permutation variants: `none`, `diss`,
  `replicate`/`rounded-replicate`, `group` (R's `permute` and
  `random-sampling` are not ported). `random_state` is a deliberate Python
  addition for reproducibility. Observed statistics are anchored to live R
  by `tests/test_r_goldens_dissassoc.py` (goldens
  `tests/r_goldens/dissassoc_*.csv`, generator
  `scripts/r_goldens/make_goldens_dissassoc.R`); p-values are stochastic and
  tested structurally. The previous simplified effect size is kept as
  **`dissassoc_effect`** (relabelled exploratory helper; `dissassoc_test`
  still wraps it).

## [0.0.4] - 2026-06-13

### Changed

- **`dissvar` is now the R-faithful weighted discrepancy** (migration plan
  Phase 2.1 / DUPLICATES A5). `dissvar(diss, weights=None, squared=False)`
  ports TraMineR 2.2-12 `dissvar` (Studer & Ritschard discrepancy analysis):
  unweighted `sum(diss)/(2n^2)`, weighted `tmrWeightedInertiaDist(var=TRUE)`.
  Returns a float, anchored to live R by `tests/test_r_goldens_dissvar.py`
  (goldens `tests/r_goldens/dissvar_phase2.csv` + shared Phase 2 fixture
  `phase2_input_*.csv`, generators `scripts/r_goldens/make_goldens_dissvar.R`
  / `make_goldens_phase2_inputs.R`). The previous simplified, unweighted
  dict-returning `dissvar` is kept as **`dissvar_summary`** (relabelled
  exploratory helper).

### Added

- **Weighted-inertia primitives** `weighted_pairwise_sum(diss, indices,
  weights)` and `weighted_inertia(...)` (ports of TraMineR's
  `tmrWeightedInertiaDist` building blocks) — the reusable SS quantities the
  rest of Phase 2 (`dissassoc`, `disscenter`, `disstree`) and
  pyweightedcluster's Phase 2.5 refactor consume.

### Added

- **DHD honours a caller-supplied substitution matrix**
  (`pairwise(method="DHD", subm=...)`). A static `(K, K)` matrix reproduces R
  `seqdist(method="HAM", sm=...)` and a time-varying `(L, K, K)` array reproduces
  R `seqdist(method="DHD", sm=...)` — the position-wise substitution sum
  `sum_t subm_t[a_t, b_t]`, with `norm="maxlength"` dividing by length. Previously
  traminex's DHD always rebuilt transition-rate costs and silently ignored a
  user `sm`, a real faithfulness gap vs TraMineR. Bare `DHD` (no `subm`) is
  unchanged. New R goldens (`tests/r_goldens/DHD_subm_*.csv`, generator
  `scripts/r_goldens/make_goldens_dhd_subm.R`) anchor it to TraMineR 2.2-12.
  Exposes the option weightedcluster's Phase 1 / A1 migration needs to route its
  position-cost-weighted Hamming distances through traminex (see
  `MIGRATION_PLAN.md`).

### Fixed

- **Editable/source installs no longer fail.** Removed an invalid
  `[tool.scikit-build.metadata]` block in `pyproject.toml` whose
  `requires-dist = []` (a list where scikit-build-core expects a per-field
  dict) made every `pip install -e .` / source build abort with
  `TypeError: Expected typing.Dict[str, typing.Any], got list`. The block had
  no effect; removing it unblocks installing traminex into the shared
  environment as the family's local dependency (see `MIGRATION_PLAN.md`
  Phase 0.3). No runtime behaviour change; test suite unchanged
  (442 passed, 1 skipped).

### Added

- **Cross-project parity harness** (`tests/r_goldens_utils.py`,
  `tests/r_goldens/*.csv`): a file-based loader for per-golden R reference
  matrices, shared by the family migration's later phases (see
  `MIGRATION_PLAN.md` / `DUPLICATES.md`). Same `%.15g` / `atol=1e-10`
  convention as the existing inline goldens.

## [0.0.3] - 2026-06-12

### Added

- **OM-family distance normalization** (TraMineR parity, live-R golden-pinned in
  `tests/test_om_normalization.py`): `TWED`, `OMslen`, `OMspell`, `OMstran` now accept
  `norm=` with the TraMineR vocabulary `maxlength` / `gmean` / `maxdist` / `YujianBo` /
  `auto` (plus `length` / `max` aliases; `auto` resolves to `YujianBo`). The shared
  transform mirrors `DistanceCalculator::normalizeDistance`
  (`SOURCES/TraMineR/src/distancecalculator.h`); each method's maximum-distance / length
  inputs are transcribed from its C++ `distance()` and documented in
  `traminex/distances/_normalization.py` and `docs/DISTANCES_NORMALIZATION.md`.
  `norm=None` leaves the raw distances unchanged.
- `NMS` / `NMSMST` / `SVRspell` explicitly **refuse** `norm` (raising `ValueError`),
  matching TraMineR's `seqdist()` (`seqdist.R` line 205).

### Changed (closeout-audit remediation, 2026-06-12)

- **BREAKING — `method="DHD"` now reproduces TraMineR's Dynamic Hamming Distance.**
  The default DHD is now Lesnard's transition-rate metric (time-varying TRATE
  substitution costs, `cval=4`, `tratecostBoth`; `seqdist.R` line 451 /
  `seqcost.R`), golden-pinned against live TraMineR 2.2-12
  (`tests/test_r_goldens_distances.py::DHD_default`). The previous marginal
  state-frequency metric was **not** TraMineR's DHD; it remains available as an
  explicit opt-in via `variant=` (`"mean"`/`"sum"`), `smooth_window=`, or
  `dhd_profile="calibrated"`.
- Version identity synchronized to **0.0.3** across `traminex/__init__.py`,
  `README.md`, and `docs/INDEX.md` (the 0.0.3 release commit had bumped only
  `pyproject.toml` + this changelog, leaving the consistency gate / CI red).

### Fixed (closeout-audit remediation, 2026-06-12)

- `tests/test_seqtree_smoke.py` used a guard (`"sklearn" not in globals()`) that
  always skipped; it now runs via `pytest.importorskip("sklearn")`.
- Pinned the two missing `OMslen` (`link="gmean"`) normalization goldens
  (`gmean`, `auto`) in `tests/test_om_normalization.py`.
- Removed a stray `#pragma once` from `cpp/bindings/pybind11_module.cpp` (warning
  in a translation-unit main file).

## [0.0.2] - 2026-06-11

### BREAKING

- **CHI2 rewritten to be R-faithful** (TraMineR `seqdist-CHI2.R` / `chisq.cpp`):
  pdotj-weighted chi-squared on per-interval state distributions with the
  default `step=1` (one interval per position), matching TraMineR's default.
  Results differ from any earlier CHI2 behavior. `norm` for CHI2/EUCLID is a
  boolean (TraMineR `norm="auto"`), distinct from the OM-style
  `norm="length"/"max"`. Verified against live TraMineR 2.2-12 goldens
  (`tests/test_r_goldens_distances.py`).

### Added

- Distance methods (TraMineR parity, all golden-pinned against live
  TraMineR 2.2-12): `EUCLID`, `CHI2`, `TWED`, `OMslen` (with `h`, default 0.5),
  `OMspell` (with `tpow`, default 1), `OMstran` (`otto`, `transindel`,
  `previous`, `add_column`), `NMS`, `NMSMST`, `SVRspell`.
- Substitution costs: `seqsubm(method="FEATURES")` (Gower over per-state
  features), `seqsubm(method="FUTURE")` (chi-squared on common-future
  transition rates), and `gower_distance()`.
- P1 additions: `seqdist` methods `WHAM`, `DHD`, `DUR` (l1/cos/chi2),
  `JSD_TRANS`; `seqsubm` `CONSTANT` and `INDELSLOG` helpers; entropy overlays
  via `state_distribution_plot(..., mode="Ht"|"dH")`; `seqfind`; Quarto docs
  site.
- R golden regression suite (`tests/test_r_goldens_distances.py`, 27 cases)
  with regeneration scripts under `scripts/r_goldens/`.

### Fixed (independent verification against live TraMineR 2.2-12)

- `NMS`/`NMSMST`/`SVRspell` returned the squared distance; TraMineR applies a
  final `sqrt` in the R wrapper (`seqdist.R`) on top of the C++ kernels.
- `OMslen` weighted positions by raw spell length and halved the
  `link="mean"` cost; R uses weight `dur^(-h)` (h default 0.5) and the sum
  `sm*(w_i+w_j)`.
- `OMspell` used raw spell durations; R uses `dur^tpow - 1`.
- TWED border sentinel `1000*maxscost` replaced with `inf` (could be exceeded
  under large `nu`/`h`); `INDELSLOG` uses `log(2/(1+p))`; OM `norm="max"`
  min-bound — all confirmed correct against R goldens.

### Changed

- Normalization standardized around `norm`; RLCP normalization added; TRATE
  weighting passthrough; small dispatcher/plotting fixes.
- New spell/subsequence methods expose raw (un-normalized) distances; `norm`
  is rejected for them. TraMineR-side `norm` support for
  OMslen/OMspell/OMstran/TWED remains a known gap (R has no `norm` for
  NMS/NMSMST/SVRspell either).
- CI consistency job runs mypy with `--follow-imports=skip` (matches
  `scripts/typecheck.ps1`).

## [0.0.1] - 2025-09-11
### Added
- Initial M1 baseline and audit integration.
- STS_to_TSE function for round-tripping state sequences to time-stamped records.
- seqtree_predict to mirror disstree_predict for feature-based trees.
- toPersonPeriod accepts time_col to harmonize with TSE_to_STS pipelines.
- Plotting functions (sequence_index_plot, state_distribution_plot, seqfrequency_plot, chronogram, transition_heatmap) now accept ax and return (fig, ax).
- seqstatd accepts weights (defaults to x.weights) for weighted state distributions; state_distribution_plot can pass weights.
- Documentation: API_CONSISTENCY.md and DISTANCES_NORMALIZATION.md with examples and guidance.
- Tests: public API has no duplicates; round-trip STS_to_TSE->TSE_to_STS; seqtree_predict; events_to_state_sequences(align_missing=True); plotting ax passthrough; seqstatd(weights).

### Changed
- Standardized missing-state flags across stats/plotting to `include_missing`.
- Distances documentation clarified to prefer norm over normalize for public API.
- __all__ is now de-duplicated while preserving order.

### Fixed
- Minor formatting and EOF fixes via pre-commit hooks.
