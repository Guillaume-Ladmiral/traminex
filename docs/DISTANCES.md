# Distance guide

Overview of `traminex.seqdist()` distance methods and the supporting cost
utilities. This page covers the **classic / core** methods; the TraMineR
parity-ported methods added in 0.0.2 are documented in
[`NEW_DISTANCE_METHODS.md`](NEW_DISTANCE_METHODS.md), and `norm=` normalization
in [`DISTANCES_NORMALIZATION.md`](DISTANCES_NORMALIZATION.qmd).

All methods are reached through `tx.seqdist(x, method=..., **kwargs)`. The
authoritative per-method kwarg list lives in
`traminex/distances/_validation.py`.

## Dispatcher inventory

The 16 TraMineR 2.2-14 methods are `OM`, `OMloc`, `OMslen`, `OMspell`,
`OMstran`, `HAM`, `DHD`, `CHI2`, `EUCLID`, `LCS`, `LCP`, `RLCP`, `NMS`,
`NMSMST`, `SVRspell`, and `TWED`. traminex adds three dispatcher extensions:
`WHAM`, `DUR`, and `JSD_TRANS`. The validation table uses uppercase,
case-normalized keys such as `OMLOC` and `SVRSPELL`.

## Core methods

- **HAM** — raw per-position mismatch/substitution-cost sum, matching TraMineR
  `norm="none"`. Use `norm="maxlength"` (or `"auto"` for fixed-length data)
  for a length-normalized result.
- **OM** — Optimal Matching (edit distance).
  - Costs: `indel` (scalar or per-position via `indel_by_pos`) and substitution
    `subm` (scalar, K×K matrix, or per-position `subm_by_pos`; see
    [`OMLOC.md`](OMLOC.md)).
  - Scalar indel costs are strictly positive (positive infinity retains the
    TraMineR sentinel behavior). Fixed substitution costs are finite and
    non-negative; matrix diagonals are zero. Position-dependent costs are
    finite, strictly positive for indels, and zero-diagonal/non-negative for
    substitution arrays.
  - `norm` accepts `None`, the full TraMineR vocabulary (`"none"`,
    `"maxlength"`, `"gmean"`, `"maxdist"`, `"YujianBo"`, `"auto"`), and the
    traminex aliases `"length"` and `"max"`.
- **R cost-method strings** — methods that consume substitution costs accept
  `subm="TRATE" | "CONSTANT" | "INDELS" | "INDELSLOG"`. `seqdist` resolves
  these through `seqcost` before method dispatch, including vector indels and
  declared-but-unobserved states; method-specific refusal and infinity behavior
  follows TraMineR.
- **`subm` is required where R requires `sm`** — `seqdist` stops for `OM`,
  `OMloc`, `OMslen`, `OMspell`, `OMstran`, and `TWED` when `subm` is missing,
  exactly as TraMineR stops with `'sm' must be specified`; only `HAM`
  (CONSTANT 1), `DHD` (TRATE, `cval=4`) and `LCS` (CONSTANT 2) auto-generate
  costs. The low-level `pairwise` keeps its documented scalar defaults.
- **OMloc** — TraMineR's context-sensitive localized OM method. It is distinct
  from the Python `omloc()` convenience function and position-dependent costs
  on regular OM; see [`OMLOC.md`](OMLOC.md).
- **LCS** — `d = |a| + |b| − 2·LCS(a,b)`.
- **LCP** — longest common prefix: `d = |a| + |b| − 2·LCP(a,b)`.
- **RLCP** — reverse LCP (longest common suffix).
- **DHD** — Dynamic Hamming Distance.
  - **Default = TraMineR's transition-rate metric** (Lesnard 2010): time-varying
    TRATE substitution costs (`cval=4`), summed position-by-position. Pinned
    against live TraMineR 2.2-14 (`tests/test_r_goldens_distances.py::DHD_default`).
  - A legacy marginal state-frequency metric is available as an explicit opt-in
    via `variant="mean"|"sum"`, `smooth_window=`, or `dhd_profile="calibrated"`.
    It has no TraMineR counterpart and no normalisation: any `norm` other than
    `None`/`"none"` is refused.
  - Calibrated profiles (`dhd_profile="calibrated"`, and
    `trate_profile="calibrated"` for OM with `subm="TRATE"`) read the first
    existing file of `calibrated_params_path`/`calibrated_params_paths`, or of
    the `TRAMINEX_CALIBRATED_PARAMS` environment variable. An explicitly passed
    value (`variant`, `smooth_window`, `trate_smoothing`, `trate_mode`,
    `trate_normalize`) overrides the file's. When a value must come from the
    file, the call raises `ValueError` if no file loads, if the file has no
    `"dhd"`/`"trate"` section, or if that section holds an unknown key; it
    never falls back to the defaults silently. A calibration path given without
    a calibrated profile is refused with `TypeError`, since nothing would read it.
    `scripts/calibrate_parity.py` writes a section only for parameters that
    reproduce its R reference (every value within an absolute 1e-6 for TRATE and
    2e-3 for DHD, the tolerances of `tests/references/README.md` without a
    relative part). No frequency variant reproduces
    R's default DHD, so against that reference it writes no `"dhd"` section.
    When nothing reproduces, it writes an empty calibration (`{}`), which a
    calibrated profile refuses, so an earlier file's sections do not survive.
  - `norm` divides by TraMineR's bound, the sum over positions of the largest
    cost, in `seqdist` and in the low-level `pairwise` alike (`"auto"` is
    `"maxlength"`; with equal lengths `maxlength`, `maxdist` and `gmean` use
    that bound). The low-level `pairwise` checks an explicit cost array or
    static `(K, K)` matrix as `seqdist` does: finite, non-negative, zero
    diagonal, covering every observed state.
  - An explicit time-varying cost array is **time-first**, `(T, K, K)`. TraMineR
    and `seqcost(..., time_varying=True)["sm"]` use `(K, K, T)`; convert with
    `np.moveaxis(sm, -1, 0)` (the refusal says so when it sees that layout).
- **WHAM** — weighted Hamming with per-position `weights` (traminex extension):
  `d(x, y) = Σ_t w_t·[x_t ≠ y_t] / Σ_t w_t`, the weighted share of positions
  where the two sequences differ (without `weights`, every `w_t` is 1 and `d`
  is the share of mismatched positions). Weights must be a one-dimensional,
  finite, non-negative vector with one value per sequence position and a
  positive sum.
- **DUR** — duration-distribution distance, `metric="l1"|"cos"|"chi2"`
  (default `"l1"`; traminex extension). With `u` and `v` the shares of
  positions each sequence spends in each state: `"l1"` is `Σ_k |u_k − v_k|`;
  `"cos"` is `1 − u·v / (|u|·|v|)` (0 when either vector is zero); `"chi2"` is
  `Σ_k (u_k − v_k)² / (u_k + v_k + 10⁻¹²)`.
- **JSD_TRANS** — Jensen–Shannon divergence on transition distributions
  (traminex extension). With `p` and `q` the two sequences' transition
  distributions (each ordered pair of consecutive states, as a share of the
  sequence's transitions) and `m = (p + q)/2`, the distance is
  `½·KL(p‖m) + ½·KL(q‖m)`, where `KL(a‖b) = Σ_{a_k > 0} a_k·(ln(a_k + eps) − ln(b_k + eps))`
  with the natural logarithm. A sequence with fewer than two positions has no
  transition, and a distance involving one is 0 (so every distance is 0 for
  sequences of length 1). Its optional logarithm-smoothing constant `eps`
  (default 0) must be finite and non-negative; only positive shares enter the
  logarithm, so with the default the distance is the Jensen–Shannon divergence
  itself (`ln 2` for disjoint transition distributions).

```python
import numpy as np
import traminex as tx

mat = np.array([[0, 1, 2], [0, 2, 2], [1, 1, 0]], dtype=int)

D_om   = tx.seqdist(mat, method="OM", indel=1.0, subm=2.0, norm="length")
D_ham  = tx.seqdist(mat, method="HAM")
D_dhd  = tx.seqdist(mat, method="DHD")                  # TraMineR transition-rate DHD
D_dhdm = tx.seqdist(mat, method="DHD", variant="mean")  # legacy marginal metric
D_wham = tx.seqdist(mat, method="WHAM", weights=np.array([1.0, 2.0, 1.0]))
```

### Distribution-based and traminex-only methods

```python
# CHI2 / EUCLID: `step` is the interval width; overlap=True slides the window
# by half a step, overlap=False uses disjoint windows. `breaks` supplies the
# intervals directly instead, as 1-based inclusive pairs, and then wins over
# both (see the seqdist Rd).
D_chi2    = tx.seqdist(mat, method="CHI2", step=2, overlap=True)
D_chi2_br = tx.seqdist(mat, method="CHI2", breaks=[(1, 2), (2, 3)])
D_euclid  = tx.seqdist(mat, method="EUCLID", step=2)

# DUR: distance between per-sequence duration distributions
D_dur_l1   = tx.seqdist(mat, method="DUR", metric="l1")
D_dur_cos  = tx.seqdist(mat, method="DUR", metric="cos")
D_dur_chi2 = tx.seqdist(mat, method="DUR", metric="chi2")

# JSD_TRANS: Jensen-Shannon divergence between transition distributions
D_jsd = tx.seqdist(mat, method="JSD_TRANS")
```

## Cost utilities

- `trate_substitution_matrix(seqmat, smoothing=0.0, mode="symmetric", normalize="none|mean_offdiag|clip_0_2")`
  — defaults follow R's `seqcost(method="TRATE")` formula (no smoothing;
  states without outgoing transitions have zero transition-rate rows).
  `smoothing>0` is an opt-in Laplace knob. The raw-array extension also applies
  this formula to singleton trajectories: distinct states have cost 2 with the
  default symmetric mode. Pinned R singleton TRATE itself raises an error.
- `trate_substitution_matrix_timevarying(seqmat, ...)` returns `[L, K, K]`,
  including `[L, 0, 0]` for an empty corpus, and validates its options on empty
  input. This raw cost utility has its own position-cost convention.
- `indel_from_subm(subm, strategy="mean_offdiag_half|mean_offdiag")`
- `local_indel_by_pos(seqmat, method="change|entropy", base=1.0, scale=1.0, smooth_window=None)`
- `seqsubm(seqmat, method="TRATE|CONSTANT|INDELS|INDELSLOG|FEATURES|FUTURE", ...)`
  — see [`NEW_DISTANCE_METHODS.md`](NEW_DISTANCE_METHODS.md) for FEATURES/FUTURE.

```python
from traminex.distances import seqsubm
from traminex.distances.costs import indel_from_subm

subm  = seqsubm(mat, method="TRATE")  # default smoothing=0; == TraMineR seqcost(TRATE)
subm_const = seqsubm(mat, method="CONSTANT", cval=2.0)
subm_log   = seqsubm(mat, method="INDELSLOG")
indel = indel_from_subm(subm_log)
```

## Dispatcher and convenience

- `pairwise(seqmat, method="OM", subm="TRATE", trate_smoothing=..., trate_mode=..., trate_normalize=..., norm="length|max")`
- `pairwise_between_sets(A, B, *, method="OM", deduplicate=False, **kwargs)`
  returns the rectangular `A × B` cross-block. Dense and canonical ragged rows
  are supported; method options and normalization are the same as `pairwise`.
- `om_trate_pairwise(seqmat, *, indel=None, trate_smoothing=..., trate_mode=..., trate_normalize=..., norm="length")`;
  `indel=None` resolves the automatic indel as R does, `max(sm)/2`

### Reference sets and dataset-fitted methods

TraMineR estimates substitution costs and state frequencies from the **whole**
`seqdata` before it selects a reference set, so `seqdist(s, sm="TRATE",
refseq=list(1:2, 3:4))` equals the corresponding block of the full-corpus
matrix exactly. traminex reproduces that: `pairwise(corpus, method=...,
ref_sets=(a0, b0))` computes the fit on `corpus` and returns the block, so it
agrees with `pairwise(corpus, ...)`, with `seqdist(..., refseq=[a, b])`, and with
R. The two number the sequences differently: `ref_sets` takes zero-based row
positions, `seqdist`'s `refseq` R's one-based indices, so `a0` is `a` minus one
(`a0 = [i - 1 for i in a]`).

A configuration is **dataset-fitted** when its parameters are estimated from the
corpus: `CHI2`, `EUCLID`, `OMSTRAN`, `DHD` without an explicit `subm`, and `OM`
with `subm="TRATE"`. For these, `pairwise_between_sets(A, B, ...)` raises
`ValueError`: given only the two sets it cannot know the corpus the parameters
should be fitted from, and silently fitting on `A ∪ B` would return distances
inconsistent with both the full-corpus matrix and TraMineR. Pass the whole
corpus and select the sets instead:

```python
# dataset-fitted: give the dispatcher the corpus it must fit on
a0, b0 = [i - 1 for i in a], [i - 1 for i in b]                # ref_sets: zero-based
block = pairwise(corpus, method="CHI2", ref_sets=(a0, b0))
block = seqdist(x, method="OM", subm="TRATE", refseq=[a, b])   # R-facing route, one-based
```

Every other method is **pair-local** — its distances depend only on the pair
being compared — so `pairwise_between_sets(A, B, ...)` keeps the efficient
rectangular path and returns the same values as the corresponding block of any
larger corpus.

`traminex.distances` is a supported submodule API with an explicit `__all__`
inventory. In particular, `pairwise` and `pairwise_between_sets` are stable
submodule functions used by sibling packages. They deliberately remain absent
from the root `traminex` manifest: importing them from `traminex.distances`
keeps the distance-engine boundary explicit and leaves the 226-name root API
unchanged. Backend implementation modules and underscored helpers are not part
of this export inventory.

Dense `pairwise`, raw-matrix `seqdist` and raw-array `seqsubm` inputs must be
two-dimensional integer arrays with non-negative state codes fitting a 32-bit
integer. The raw substitution-cost builders enforce the same contract.
Invalid numeric dtypes are not
silently truncated before dispatch: floats and booleans are not coerced to
state integers.

## Backends

The Python/Numba path is always available. The supported optional compiled
backend is pybind11, selected with `TRAMINEX_BINDING=auto|python|pybind11` (auto
the default); dispatch falls back to Python only under automatic selection. An explicit unavailable or
unknown binding is an error, so a requested compiled gate cannot silently pass
on Python. See [`DEV_BUILD.md`](DEV_BUILD.md).

Numba accelerates several pure-Python kernels automatically when installed. It
is a **third** backend, orthogonal to `TRAMINEX_BINDING`, and is selected by
`TRAMINEX_NUMBA=auto|off` (default `auto`). Like `TRAMINEX_BINDING`, it is read
once, when `traminex.distances` is first imported -- not by `import traminex`,
which resolves public names lazily -- so assigning to `os.environ` after the
first distance computation has no effect; set it before the process starts. An
unknown value is an error. `off` uses the reference NumPy
kernels; it is the supported way to reproduce a result on the reference path,
and the way `tests/test_numba_backend_agreement.py` compares the two backends
across all seven production dispatch sites. The two backends are not
byte-identical everywhere: `JSD_TRANS` and `DUR` with `metric="chi2"` can differ
in the last bits, from floating-point reassociation inside the Numba kernels.
The difference is a few units of double-precision rounding and grows with the
number of states and the sequence length, so no fixed bound applies; for `DUR`
it appeared only with more than four states in the measurements. The other
routes were byte-identical in every case measured. Use `TRAMINEX_NUMBA=off` when a result
must be reproduced bit for bit on the reference path.

The maintained benchmark CLI checks agreement before warmup or timing and
applies no performance threshold.

### Speed of the row engine

OM with per-state indel costs (`INDELS`, `INDELSLOG`, an indel vector) or an
asymmetric substitution matrix, OMstran, OMslen, OMspell, OMloc, TWED, LCS on
sequences of unequal length (and on equal lengths when Numba is not installed
or is switched off), and NMS with hard or soft matching run in the Python row
engine. It evaluates each dynamic program for many pairs at once,
cell by cell with the same floating-point operations in the same order as the
per-pair reference, so its results equal that reference bit for bit
(`tests/test_distance_batched_engine.py`). In the row engine a `refseq` call
evaluates only the pairs it returns; the dense Numba OM, LCS, LCP and RLCP
kernels, DHD and CHI2/EUCLID compute the full matrix and return the requested
part, with the same values. The batched sweep costs about as much for one pair as for a
few dozen, so when a call evaluates fewer pairs than the method's measured
crossover (5 for OM and OMslen, 6 for OMspell and OMstran, 7 for OM on equal
lengths, 8 for LCS and OMloc, 9 for TWED) the engine runs the per-pair
reference instead, with the same values; NMS is batched at any number of
pairs, because its batched counts are faster even for one pair. Before this
change the engine ran one pure-Python dynamic program per pair.

Measured on the first 60 `mvad` sequences (months 15 to 86, so 72 positions,
weighted), `seqdist` with default settings unless stated, on an Intel Xeon
Silver 4214R with Python 3.12.14, NumPy 2.3.5 and Numba installed, no compiled
backend, against R 4.6.1 with TraMineR 2.2-14 on the same machine. "Before" and
TraMineR are single runs of 2026-10-13; "After" is the faster of two runs on
2026-10-17 at a load average of 9 to 10, so other work shared the machine:

| Method | Before | After | TraMineR |
| --- | ---: | ---: | ---: |
| OM, `subm="INDELSLOG"` | 15.8 s | 0.29 s | 0.22 s |
| OMslen | 20.2 s | 0.44 s | 0.21 s |
| OMspell | 0.15 s | 0.02 s | 0.18 s |
| OMloc, `expcost=0.5` | 20.4 s | 0.66 s | 0.23 s |
| OMstran, `otto=0.1` | 15.6 s | 0.29 s | 0.48 s |
| TWED, `nu=0.5`, `h=0.5` | 19.2 s | 0.62 s | 0.21 s |
| NMS | 199 s | 3.7 s | 0.96 s |
| NMS, soft (`prox` 1 on the diagonal, 0.3 elsewhere) | 221 s | 4.7 s | 0.70 s |
| OMslen, `refseq=0` | 19.7 s | 0.07 s | 0.18 s |

Every result agrees with TraMineR to a relative difference of at most
2.6e-14, as it did before. NMS remains about four to seven times slower than
TraMineR: to reproduce TraMineR's results its counts must be added in
TraMineR's order, which leaves NumPy a long sequence of small steps per pair.
On this workload, with 3.6 spells per sequence, NMSMST (0.39 s against 0.20 s)
and SVRspell (0.43 s against 0.18 s) were close to TraMineR. Both count over
spells, though, and on sequences with many spells their Python kernels were 40
to 140 times slower than TraMineR (2026-10-25 study, finding F-06). Since
2026-11-07, when Numba is dispatched, both run compiled copies of the same
loops, whose counts equal the Python kernels' bit for bit
(`tests/test_nmsmst_svrspell_parity.py`). Measured on 2026-11-07 on the same
machine, with Python 3.12.16, NumPy 2.3.5, Numba 0.67.0 and R 4.6.2, at a load
average of about 15, single runs, not counting Numba's compilation of each
method, about 0.6 to 1.3 s once per process:

| Workload | NMSMST, Numba off | NMSMST, Numba | TraMineR | SVRspell, Numba off | SVRspell, Numba | TraMineR |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 30 random sequences of length 40, 4 states | 9.90 s | 0.03 s | 0.16 s | 23.16 s | 0.05 s | 0.17 s |
| 60 random sequences of length 20, 4 states | 6.17 s | 0.05 s | 0.14 s | 18.63 s | 0.07 s | 0.15 s |
| first 60 `mvad` sequences, columns 17 to 86 (70 months) | 0.18 s | 0.04 s | 0.18 s | 0.40 s | 0.04 s | 0.17 s |
| first 200 `mvad` sequences, columns 17 to 86 (70 months) | 1.87 s | 0.40 s | 0.26 s | 3.68 s | 0.42 s | 0.31 s |

`benchmarks/run_benchmarks.py`
times only the common methods (DHD, DUR, HAM, JSD_TRANS, LCP, LCS, OM, RLCP);
the figures above come from timing `seqdist` directly on this workload.

## Parity

Classic and parity-ported methods are pinned against live TraMineR 2.2-14 goldens
(`tests/test_r_goldens_distances.py`, `tests/test_om_normalization.py`). See
[`API_PARITY.md`](API_PARITY.md) and the optional reference fixtures described in
[`REFERENCES.md`](REFERENCES.md).
