# API parity with TraMineR and TraMineRextras

The current R oracle is TraMineR 2.2-14 plus TraMineRextras 0.6.9. Exact CRAN
archive checksums, Ubuntu build instructions, and the oracle refresh policy are
recorded in [`SOURCES/README.md`](../SOURCES/README.md).

## Authoritative and generated sources

- [`traminex/api_manifest.py`](../traminex/api_manifest.py) is the canonical
  ordered Python public API.
- [`API_PARITY_DECISIONS.csv`](API_PARITY_DECISIONS.csv) is the hand-reviewed
  semantic disposition of every exported/S3 name parsed from the two pinned R
  namespaces.
- [`API_PARITY_LEDGER.md`](API_PARITY_LEDGER.md) and the
  [machine-readable semantic JSON](../reports/audit/traminex_semantic_parity.json)
  are deterministic generated views. Their summaries publish current status,
  scope, option coverage, owner, target-phase, and public-API totals.

The older [name audit](../reports/audit/traminex_api_audit.md), its
[JSON](../reports/audit/traminex_api_audit.json), the
[name-parity table](API_PARITY_TABLE.md), and the files under
[`reports/parity`](../reports/parity/) are discovery artifacts. Exact or fuzzy
name similarity is not behavioral parity.

Focused tests load the current-oracle golden corpus under `tests/r_goldens`;
the operational regeneration contract is documented in
[`scripts/r_goldens/README.md`](../scripts/r_goldens/README.md).

Reproduce the current totals instead of copying them into prose:

```bash
PY=~/miniforge3/envs/datasci/bin/python
$PY scripts/generate_semantic_parity_ledger.py --check
$PY - <<'PY'
import json
from pathlib import Path

report = json.loads(
    Path("reports/audit/traminex_semantic_parity.json").read_text(encoding="utf-8")
)
print(report["decision_count"])
print(report["summaries"])
print(report["target_public_api"])
PY
```

## Evidence and validation contract

The CSV schema is fixed. Validation enforces namespace/package/kind identity,
closed categorical domains, repository-relative existing evidence, and current
public mappings for executable statuses. Evidence must cite the defining
public-R source (or Rd), not merely an existing file from the wrong package.

For executable rows:

- `EXACT` requires `FULL`, reviewed R source, and focused behavioral evidence.
  Nontrivial deterministic analytics normally cite their public-R generator
  and static fixture as well.
- `DIVERGENCE` names the observable R behavior, the tested Python behavior,
  and the exact differing boundary.
- `RE_EXPRESSED` names the Python object, protocol, setter, or plotting route
  replacing the R representation.
- `FULL` means every effective documented/source-active option and default in
  the meaningful mapped contract is supported and tested. Known unsupported
  options force `PARTIAL`; `N/A` cannot hide an analytical option gap.
- Known boilerplate notes are rejected. Generated summaries, including option
  coverage, are derived only from the authoritative CSV.

The direct negative and determinism tests live in
[`tests/test_semantic_parity_ledger.py`](../tests/test_semantic_parity_ledger.py).

## Sequence alphabet order

Without `alphabet=`, every constructor orders the real states as TraMineR's
`seqdef` does (owner decision, 2026-10-19): the `seqstatl` order of the values
it sees -- numeric order when every value, a `missing=` code included, is an R
number; R's collation otherwise (ICU root collation, see
[`DEV_BUILD.md`](DEV_BUILD.md)) -- followed by the replacement states of
`left`, `gaps` and `right`. This holds for `seqdef` (STS, SPS and SPELL input),
`seqformat` (STS output), `seqdef_from_long`, `TSE_to_STS` (state and event
mode), `HSPELL_to_STS`, `events_to_state_sequences` and the combined alphabet
of `seqMD`. Objects built from another object (`seqrecode`, `seqdss`,
`seqgranularity`, subsets) keep their input's order, as in R. `states=`,
`labels=`, `cpal=` and the default palette are paired with that order, so the
stored alphabet, each state's name, label and colour, per-state rows and
columns, tie rules (`seqmodst`, `seqgranularity(method="mostfreq")`) and
state-rank defaults (`seqibad`, `seqidegrad`, `seqprecarity`,
`seqinsecurity`) equal R's. `labels=` and `cpal=` give one entry per real
state, as in R, or one per codebook entry; by default the labels are the
states and the colours `seqdef`'s palette, and a marker takes its own text
and `missing_color` (nr) or white (void) (2026-10-23, TX-AUD-20261019-035,
TX-REP-20261019-013). An explicit `alphabet=` keeps the caller's order,
as in R. `seqformat`'s own output keeps the declared categories of a
categorical SPELL status, as R's factor status keeps its levels, and the
alphabet order of an STSList input (a traminex choice: R returns a plain
character matrix there, and `seqdef` of either result re-sorts). The cases
are recorded from R by
[`make_goldens_default_alphabet.R`](../scripts/r_goldens/make_goldens_default_alphabet.R).
Until 2026-10-19 the default was the order of first appearance. To obtain a
specific order, pass `alphabet=` to `seqdef`; `seqformat`, `seqdef_from_long`,
`TSE_to_STS`, `HSPELL_to_STS` and `seqMD` take no `alphabet=`, so pass their
result to `seqdef(x, alphabet=[...])`, which keeps every cell and marker;
`events_to_state_sequences` takes `states=`.

A data frame mixing numeric and text columns is read cell by cell with R's
`as.character`, so `1` and `"1"` are one state. R's `seqdef` runs `as.matrix` on
the frame, which pads a numeric column's numbers (`" 1"` beside `"1"`) or makes NA
cells and a duplicated alphabet: a recorded correction of R (TX-AUD-20261019-062;
R's results in `scripts/r_goldens/make_goldens_construction_items.R`).

One deliberate difference: R re-renders a numeric `seqstatl` result with
`as.character` and then matches cells by text, so a state spelled `01` or
`0x10` becomes `NA` in R; `seqdef` keeps such a cell as its state and orders it
by its numeric value. A `states=` vector must have one entry per state of the
alphabet, as the Rd requires; R does not check and pairs a longer vector's
first entries with the observed states.

A missing cell is decided by `seqdef`'s default rules (`left = NA`,
`gaps = NA`, `right = "DEL"`) in `seqdef` itself when those arguments are not
given, and always in the other three constructors: leading and internal
missing cells become the nonresponse marker and trailing ones void, the
markers following the real states. For `seqdef_from_long` and the state mode
of `TSE_to_STS` that covers a period with no row for a person as well as an
explicit missing state; for the event mode of `TSE_to_STS` it covers the NA
cells an NA transition in `stm` leaves, exactly as `seqdef` does on R's result.

Three `seqdef` inputs that R accepts are refused (2026-10-17, TX-AUD-20261010-035;
R's results are recorded by `scripts/r_goldens/make_goldens_seqdef_nr_void_neutral.R`,
and for the third by `scripts/r_goldens/make_goldens_construction_items.R`):

- **`nr` equal to `void`.** R then handles a missing cell as void, so the
  missing observation disappears (the row `NA, A, B` becomes `A-B`, length 2),
  although `seqdef.Rd` gives the two codes separate roles. traminex refuses the
  ambiguous coding: a recorded correction of R.
- **`left`, `gaps` or `right` = `"NEUTRAL"`.** R keeps a plain `NA` cell that is
  neither the missing nor the void code (`NA-A-B`, length 3). An STSList has no
  such cell; supporting it would change the data model, so it is refused (owner
  decision, 2026-10-18).
- **`left`, `gaps` or `right` equal to the `nr` marker (`"*"`).** R lists the
  marker as a state of the alphabet (`A B *`), twice when two of them name it
  (`A B * *`), so one text is both a state and the missing code. traminex
  refuses a replacement state spelled like the `nr` marker: a recorded
  correction of R (TX-AUD-20261019-040).

## Phase 4 corrections

The source-wide review removed unqualified exact claims from stochastic
operations whose deterministic components match R but whose random streams do
not. It also corrected container/helper and plotting-option claims:

| R row(s) | Final correction |
| --- | --- |
| `seqgen`, `dissmfacw`, `dissrep`, `seqrep`, `disstree`, `seqtree` | `DIVERGENCE/FULL`; deterministic contracts remain evidenced, but RNG-dependent results are not cross-language stream exact |
| `rowmode` | `DIVERGENCE/FULL`; scalar behavior is R-golden backed, while Python returns `None` for R zero-length results and applies a row-wise 2-D extension where R flattens |
| `seqhasmiss` | `DIVERGENCE/PARTIAL`; Python returns silent combined booleans rather than R's printed counts and separate missing/nonresponse/void vectors |
| `plot.stslist` | `DIVERGENCE/PARTIAL`; tested selection behavior remains, but identified R graphics controls have no supported target |
| `disstree2dot`, `disstree2dotp`, `disstreedisplay`, `seqtree2dot`, `seqtreedisplay` | `DIVERGENCE/PARTIAL`; portable DOT/text behavior remains, while named R device/image/layout controls have no supported target |

The R1–R4 scheduling rows remain closed with no `PLANNED` decision. The
historical public-name additions remain unchanged; Phase 4 corrects semantic
claims and one existing `seqrep` typing/result path rather than adding a new
top-level name. [`ROADMAP.md`](../ROADMAP.md) retains the phase inventory and
links back to the generated summaries.

## Recorded informational divergences (cycle 9)

The 2026-09-04 informational notes are historical observations. Current
boundaries and corrections are recorded below; the earlier description did
not establish that every supported codebook leaves numerical results unchanged.

- **Exploratory sequence-tree features canonicalize markers (2026-10-08).**
  `seqtree_fit`, `seqtree_cv_fit` and `seqtree_predict` resolve declared external
  missing/void codes before building features. Equivalent codebooks now produce
  the same duration/transition columns, feature names and entropy normalization.
  Previously external markers could cause reshape errors or change predictions;
  that was a defect, not an informational padding difference. These helpers
  remain exploratory sklearn models. `cluster_representatives` selects medoids from supplied distances and decodes
  sequence strings through the declared codebook; it does not build these
  features. `cluster_summary` already canonicalizes at entry.
- **EMLT and fitted Markov likelihood use actual state codes (2026-10-08).**
  Factory-generated missing-first and interleaved codebooks map real states to
  situation indicators and fitted transition axes by code. Missing/void markers
  do not become ordinary EMLT situation columns. Existing hand-built objects
  remain subject to the STSList input contract; this correction adds no new
  factory or layout restriction.
- **Plot group-panel level order uses codepoint sort.** `seqedplot` and
  `seqHtplot(group=)` panels and the grouped transversal-entropy helper order
  text levels by Python/pandas codepoint sort where R factor levels use locale
  collation, so panel order can differ even for ASCII labels that differ in
  case or start with punctuation: for `b, A, a, B, _x` traminex gives
  `A, B, _x, a, b` and R `_x, a, A, b, B` (2026-09-22 audit S3-23; the note
  said "identical for ASCII"). `plot_SAMM` origin groups use the package's
  case-folding order, which agrees with R on those labels; accents can still
  differ. Numeric levels sort numerically, as in R.
- **CHI2/EUCLID `step=0` is refused with a named error.** R's global check
  lets 0 through (`is.a.positive.integer(0)` is TRUE) and CHI2 then stops
  inside `seq()` ("invalid '(to - from)/by'"), unless `breaks` replaces the
  step; traminex refuses it as "'step' must be a positive integer", in
  `seqdist` and in the low-level `pairwise` (2026-10-19, batch-6 review; it
  used to read 0 as the whole sequence length).

## `seqipos` volatility across an internal gap (2026-09-22)

`seqipos(index="volatility")` recodes the states outside `pos.states` and
`neg.states` to void (`seqipos.R:51`): the missing state when
`with.missing = FALSE`, and any state left out of explicit `pos.states` and
`neg.states`. With `dss = FALSE` (the default for this index) that void can
sit inside a sequence, but `seqdss.R` walks only the first `seqlength` cells,
as if void were always trailing padding, so R loses every spell after the
gap. On `A-B-*-*-C-A` with `pos.states = "A"`, R reads `p-n` instead of
`p-n-n-p` and returns 0.667; `seqindic`'s `Vpos` inherits the same value.

**traminex does not reproduce this.** It strips the void wherever it lies and
counts every observed cell, which gives R's own `seqivolatility` on the
gap-free binary sequence `p-n-n-p` (0.833), the value R's `dss = TRUE` route
also returns on this example. The R values are recorded by
`scripts/r_goldens/make_goldens_audit20260922_s3b.R` and asserted in
`tests/test_seqipos_volatility_gaps.py`, which fails if a future oracle stops
dropping the spells (audit S3-12).

## `HSPELL_to_STS` horizontal spells (2026-09-22)

Horizontal input reproduces `TraMineRextras::HSPELL_to_STS` followed by R's
`seqdef` defaults, cell by cell (`tests/test_r_goldens_hspell_matrix.py`): one
row per input row in input order, a missing begin skips a spell and a missing
end runs it to `tmax`, uncovered cells are NA unless `fillblanks`,
`fixed_status` is one status per spell, `endObs` blanks every cell from time
`endObs` on, and time labels are `a<t>`. Four differences are deliberate:

- **`endObs` blanks from `endObs` on (2026-10-17, TX-REP-20261017-005).** R
  shifts the spell times by `1 - tmin` but `endObs` only by `- tmin`
  (`HSPELL_to_STS.R:25,67-68`), so `seqindex >= endObs` also blanks the period
  `endObs - 1`, before the end of observation (a spell X 2000-2002 with
  `endObs = 2003` loses 2002). traminex applies the shift the spell times get,
  so `endObs` is the first unobserved period (owner decision, 2026-10-18);
  R's rows stay recorded in `tests/r_goldens/hspell_horizontal_matrix_cases.txt`
  (`endobs_scalar`, `endobs_column`).

- **An open spell with `tmin < 1` runs to `tmax`.** R sets a missing end to
  the *unshifted* `tmax` after shifting every time by `1 - tmin`
  (`HSPELL_to_STS.R:37`), so with, say, ages starting at 0 the spell stops one
  period early. traminex ends it at `tmax`; the R value is recorded in
  `tests/r_goldens/hspell_horizontal_matrix_cases.txt` (`open_end_tmin_zero`).
- **Numeric column selectors are zero-based** here (negative positions count
  from the end), as in `seqdef`; R's are one-based. traminex is not uniform:
  TraMineR's `seqformat` (`var=`, for a matrix and a data frame alike) and
  TraMineRextras' `FCE_to_TSE` (`cols=`) and `TSE_to_STS` in event mode keep R's
  one-based selectors. A selector outside the table, or a
  name it does not have, stops with R's `undefined columns selected`; this holds
  for every column argument (`begin`, `end`, `status`, `id_col`, `pvar`, `endObs`,
  `tmin_col`, `tmax_col`, and the vertical layout's `start_col`, `end_col` and
  `state_col`). A column that cannot be read as times is named in the error
  (TX-AUD-20261010-036).
- **Duplicate ids are refused** with a named `ValueError`; R fails later, when
  it sets duplicate row names.

## Deliberate divergence from the 2.2-14 oracle (2026-09-20)

TraMineR 2.2-14 added an internal `checkcost()` call inside `seqcost()`. That
call aborts on any substitution matrix containing `Inf`, because it computes
`tol <- 1e-07 * max(sma)`, which is `Inf`, and then evaluates
`any((sm - t(sm)) > tol)` — that is `any(NaN > Inf)`, which is `NA`, and `if
(NA)` stops with `missing value where TRUE/FALSE needed`. The symmetry test is
unchanged from 2.2-13; what is new is that `seqcost()` now reaches it.

`INDELS` legitimately produces such a matrix: a state absent from the data (or
from one time slice) gets an indel of `.Machine$double.xmax`, so the cost
between two absent states, `xmax + xmax`, overflows to `Inf` (one absent state
alone gives `xmax + 1/f`, which stays finite). The observable consequence in
2.2-14 is that `seqsubm(seqdata, method = "INDELS")`, `seqdist(..., sm =
"INDELS")` and the time-varying `INDELS` costs abort whenever two or more
states of the alphabet are absent from the data — for time-varying costs, from
any one time slice — where 2.2-13 returned values. With at most one absent
state 2.2-14 computes them, and `INDELSLOG` (`log(2 / (1 + f))`, always finite)
never aborts. On `mvad[1:100, ]`, for example, OM with `sm = "INDELS"` succeeds
(every state occurs) while DHD with `sm = "INDELS"` aborts (some months lack
two states).

**traminex does not reproduce this.** It computes those results as before. An
*explicitly supplied* matrix with a non-finite cost that an observed state can
reach, or one from which an automatic indel would be infinite, is refused with
a named `ValueError`; a non-finite cost between states the data never use is
accepted and computed when the indel is given, because no distance can reach
it. Both R versions stop on all of these with only `missing value where
TRUE/FALSE needed` (`tests/test_seqdist_infinite_cost_contract.py`). Reproducing the abort would be a regression in
working, tested functionality, so it is recorded here as a deliberate
divergence rather than inherited.

The affected golden generators pass 2.2-14's own `check.costs = FALSE` switch,
which controls only the diagnostics and changes no recorded `indel` or `sm`
value. The full 2.2-13 → 2.2-14 delta is recorded in the report of the oracle-2214
substitution audit of 2026-09-20, kept in the family supervision repository's
`docs/history/run-reports/`; its evidence tree was removed on 2026-10-26.

A second, narrower divergence follows from the same `checkcost` rewrite, which
demoted a zero indel from an error to a warning and added division-by-zero
guards to `OMstran`. Three degenerate inputs that 2.2-13 refused are now
accepted:

| case | what 2.2-14 returns | traminex |
|---|---|---|
| `OMstran`, scalar `indel=Inf` | `NaN` off the diagonal | refuses: *indel cost must be finite and strictly positive* |
| `OMstran`, `indel=[Inf, 1]` | `NaN` off the diagonal | refuses, same message |
| `OM`, all-zero `sm`, automatic indel | an all-zero matrix | refuses: *indel cost must be strictly positive* |
| `OM` with a cost-method string (`sm="CONSTANT"` etc.) and a scalar `indel` of 0 or below | negative distances, or an all-zero matrix: R checks the indel only for an explicit matrix | refuses on every path, as the dense path always did (the ragged one returned negative distances even between identical sequences until 2026-09-22) |
| `TWED` with costs small relative to `nu`/`h` | `TWEDdistance.cpp` fills the impossible border cells with the finite stand-in `1000 * maxscost`, which then binds and caps distances (e.g. `20` for a length-1 sequence against longer ones with costs of 0.02 and `nu = 5`) | keeps the infinite border of the TWED definition, so `R <= traminex` wherever the stand-in binds and the two agree otherwise |
| `refseq` index sets containing `0` | silently drops the `0` (`is.a.positive.integer` accepts 0 and R indexing skips it) | refuses by name |
| `HAM` with an invalid explicit `sm` (negative, non-zero diagonal) | computes: a dangling `else` in `seqdist.R:377-386` means R never calls `checkcost` for HAM | refuses by name, and gives the triangle-inequality diagnostic and the asymmetry warning |
| An explicit `sm` whose largest cost is 0 or below | `checkcost` reports *max substitution cost less or equal to 0!* and returns before its refusals, so a negative matrix gives negative distances and a wrong-sized all-zero one gives meaningless tiny values | an all-zero matrix of the right size gets the same message; a negative or wrong-sized one is refused as any other is |
| `LCS` with `norm="gmean"` or `"maxdist"` on a one-state alphabet | 1 between every two sequences, from its 1x1 cost matrix (with a two-state alphabet, the definition) | the definition: `1 - LCS/sqrt(l1 l2)` and `raw/(l1 + l2)` (TX-AUD-20261011-038) |
| `CHI2` or `EUCLID` with one sequence or one state | stops: *'x' must be an array of at least two dimensions* | 0, a sequence's distance to itself and between identical state distributions (TX-AUD-20261011-040) |
| `OMstran` whose transition alphabet has one element (one state, or no transition between two) | stops: its loop `1:(alphabet_size-1)` runs `1:0` (`seqdist-OMstran.R:126`) | the value R's algorithm gives (TX-AUD-20261011-039) |
| Transition-rate costs (OM with `sm="TRATE"`, DHD) on a one-state alphabet or sequences of length 1 | stops (*subscript out of bounds*, *undefined columns selected*) | refuses, as its `seqcost` and `seqtrate` do (TX-AUD-20261011-039, -040) |

Since 2026-09-22 the refusals `checkcost` does make -- a wrong-sized explicit
matrix, a negative cost, a non-zero diagonal -- apply to every OM-family and
TWED path, regular or ragged, not only to dense OM
(`tests/test_seqdist_cost_checks.py`).

`checkcost`'s non-fatal diagnostics are reproduced word for word for a
caller-supplied matrix (OM family and TWED) and a DHD time-varying array
(`tests/test_checkcost_diagnostics_parity.py`, against live 2.2-14 via
`scripts/r_goldens/make_goldens_audit20260922_s10.R`): a largest cost of 0 or
below, zero costs between distinct states, a triangle-inequality violation, an
indel triangle-inequality violation or a cost above twice the indel, and
asymmetry. R's `message()` diagnostics are
`traminex.distances.python_kernels.CostDiagnosticMessage` warnings and its
`warning()` a `UserWarning`. 2.2-14 also runs `checkcost` inside `seqcost()`
(`seqcost.R:286`), so it prints these diagnostics for every matrix it
generates, and since 2026-10-18 traminex does too (owner decision; audit
TX-AUD-20261011-080, formerly the maintainer question S10-01): for the matrices
`seqcost` and `seqsubm` return, once per `seqdist` call (with or without
`refseq`) with a cost-method string,
for the constant costs LCS and HAM without `sm` generate
(`seqdist.R:370-374, 452-455`) and for DHD's generated time-varying costs, with
`seqcost`'s own indel (`tests/test_r_goldens_generated_cost_messages.py`).
`seqMD` silences them for each channel, as R's
`suppressMessages(seqcost(..., check.costs = FALSE))` does.

`seqMD` (and `seqdistmc`) with `what="cost"` or `"diss"` gives the same
diagnostics for every channel's costs, in R's order, before the channels are
combined: `seqMD.R:266-267` calls `checkcost` itself, on explicit costs and
equally on the costs it generates from a method string, which it builds with
`check.costs = FALSE` precisely because it checks them there (`seqMD.R:241-243`).
That call is seqMD's own and is in 2.2-13 too, so it is not the `seqcost`
question above, and traminex follows it (2026-09-22 audit S10-11;
`tests/test_seqmd_cost_diagnostics.py`, against live 2.2-14 via
`scripts/r_goldens/make_goldens_audit20260922_s10b.R`). With `what="diss"`, the
combined costs then get `seqdist`'s own diagnostics, as in R.

The R goldens record what R actually does, so the oracle artifact stays
truthful; the divergence is asserted in
`tests/test_r_goldens_distances_void.py::DELIBERATE_REFUSALS_OF_R_ACCEPTED_INPUT`,
guarded by a test that fails if a future oracle refuses them again. The
rationale is the one already stated in
`tests/test_seqdist_infinite_cost_contract.py`: traminex need not reproduce
R's internal error text, but it must not return a silently meaningless matrix.

## Stochastic contracts

No mapped Phase 4 stochastic contract claims exact seeded R/NumPy stream
identity. The approved evidence model records the R stream where practical,
pins deterministic sufficient components, then tests the actual Python API for
reproducibility, schema, invariants, edge behavior, and nonmutation.

| Public R row(s) | Deterministic evidence | Exact cross-language limit |
| --- | --- | --- |
| `seqgen` | shape, alphabet/probability mechanics, R global-RNG advancement | Python uses a local NumPy PCG64 generator; sampled states differ from R Mersenne-Twister |
| `seqgen.missing` | rounded counts, precedence, marker metadata, boundaries | local NumPy `random_state` draws replace R global `sample`/`runif` draws |
| `seqpolyads` | resampling roles, weights, replacement, draw-plan invariants | Python honors the supplied local seed and avoids R's overwritten-seed defect; sampled polyadic sequences can differ |
| `dissrf`, `seqrf`, `seqplot.rf` | group sizes, medoids/distances, deterministic ordering, weighted MDS | random group allocation uses an isolated NumPy generator rather than R global sampling and can alter groups/plots |
| `plot.stslist` | systematic selection, weights, factor order, missing-last behavior | random row selection uses NumPy and can choose different rows; Python also has a different default |
| `dissassoc` | observed weighted pseudo-statistics and every supported permutation mode | permutation/bootstrap p-values use a local NumPy stream and are not R-stream exact |
| `dissmfacw` | formula-order deterministic decomposition and complete-case behavior | permutation p-values use a local NumPy stream and are not R-stream exact |
| `disstree`, `seqtree` | R=1 structure, weighted statistics, assignments and rules | split/final-adjustment permutations use a local NumPy stream and can change the returned tree |
| `dissrep`, `seqrep` | deterministic criteria, scores, representatives, quality and metadata | random criterion uses NumPy global draws; random scores, ordering, and representatives can differ from R |
| `dissCompare`, `seqCompare`, `seqLRT`, `seqBIC` | deterministic aggregation, balanced-plan invariants, schema and collation boundary | sampled comparisons use NumPy draw plans rather than R `set.seed`/`sample`, so sampled statistics can differ |
| `seqpcplot`, `plot.seqpcplot` | Python path/filter/artist structure | R may consume/set global RNG for overlap placement; Python uses a deterministic Matplotlib re-expression |

The asymmetric R fixtures for `seqgen` and random representatives are produced
by
[`make_goldens_stochastic_phase4.R`](../scripts/r_goldens/make_goldens_stochastic_phase4.R)
and loaded by
[`test_r_goldens_stochastic_phase4.py`](../tests/test_r_goldens_stochastic_phase4.py).
Other rows cite their focused generator/test evidence directly in the ledger.

## Ownership boundary

traminex owns sequence construction, distances/costs, statistics, discrepancy
analysis, R-faithful trees, and its matplotlib renderers that re-express the
TraMineR plots. `disstree_fit` and `seqtree_fit` are exploratory scikit-learn
helpers (an MDS embedding or per-sequence features fed to a decision tree), not
TraMineR tree mappings.

Dedicated grammar-of-graphics sequence rendering belongs to the
`py_seqplot` sibling, and Ward-constrained PAM belongs to
`weightedcluster`, as defined by [`ECOSYSTEM.md`](../ECOSYSTEM.md). The semantic
ledger records those rows as `SIBLING`; this phase does not edit or migrate
sibling projects.

## Regeneration on Ubuntu

```bash
PY=~/miniforge3/envs/datasci/bin/python
$PY scripts/audit_traminex.py
$PY scripts/generate_parity_table.py
$PY scripts/generate_api_mappings.py
cp reports/audit/traminex_api_audit.csv docs/API_PARITY_TABLE.csv
$PY scripts/generate_semantic_parity_ledger.py
$PY scripts/check_generated_artifacts.py
$PY scripts/check_docs_consistency.py
$PY -m pytest -q tests/test_semantic_parity_ledger.py
```

Edit only `docs/API_PARITY_DECISIONS.csv` for semantic decision changes. Do not
hand-edit the generated semantic Markdown/JSON, audit reports, name tables, or
R fixtures.
