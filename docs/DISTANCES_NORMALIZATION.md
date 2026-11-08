# Distance normalization reference

This note documents the `norm=` post-normalization options for sequence distances in
traminex, for both the classic OM/LCS family and the OM-family spell/warp methods
(OMloc, TWED, OMslen, OMspell, OMstran). All behavior is a port of TraMineR 2.2-14
(`SOURCES/TraMineR`); the per-method maxima are pinned against live-R goldens in
`tests/test_om_normalization.py`.

Normalization is backend-independent. The pure-Python/Numba path and maintained
pybind11/OpenMP kernels are compared on raw distance contracts before benchmark timing;
normalization remains in the shared dispatcher and is pinned by the same R
goldens.

## Classic OM/LCS family

For `OM`, `LCS`, `LCP`, `RLCP`, `seqdist` accepts **the full TraMineR
normalization vocabulary as well as** traminex's two friendly aliases. Every
value below is verified R-exact against live TraMineR 2.2-14.

traminex-friendly aliases (true aliases — same numbers, verified by
`tests/test_audit_fixes_distances.py`):

- `'length'` — exactly TraMineR `norm="maxlength"` (divide by the longer
  length, which for `OM` is the sequence length times the indel cost), for
  `OM`, `LCS`, `LCP`, and `RLCP` alike.
- `'max'` — exactly `norm="maxdist"` (divide by the maximum possible
  distance), for `OM`, `LCS`, `LCP`, and `RLCP` alike.

`seqdist` accepts both aliases for every other method that takes a TraMineR
`norm` too (`OMloc`, `OMslen`, `OMspell`, `OMstran`, `TWED`, `HAM`, `DHD`),
with the same meaning (`tests/test_distance_wrappers_norm.py`).
- `None` — no normalization; the same as `'none'`.

TraMineR's own six names, all accepted here too:

- `'none'`, `'maxlength'`, `'gmean'`, `'maxdist'`, `'YujianBo'`, `'auto'`.

Spellings are case-insensitive. `'auto'` resolves as TraMineR does —
`maxlength` for `OM`, `gmean` for `LCS`/`LCP`/`RLCP`.

`seqdist`'s `method` and `subm` are read in any letter case too (`"om"`, `"trate"`), and so
are `seqMD`'s `method` and `norm` (its `norm` goes through `seqdist`, aliases included). TraMineR accepts only its exact spellings and stops on the others,
on `norm=NULL` and on the aliases `'length'` and `'max'`. These spellings are a traminex
convenience (TX-AUD-20261019-061): each gives TraMineR's result for the canonical
spelling, but code that uses them does not run in R unchanged.

> The **low-level** `traminex.distances.pairwise` deliberately keeps the
> narrower classic vocabulary (`None`, `'length'`, `'max'`) **with its own
> historical semantics**: its `'length'` divides LCS/LCP/RLCP distances by
> `l1 + l2` (not `max(l1, l2)`) and an OM distance by the longer length
> (without TraMineR's indel factor, so with an indel other than 1 it is not
> R's `maxlength`; `om_trate_pairwise` uses it by default), on regular and
> variable-length rows alike,
> while its `'max'` divides by the maximum possible distance, as `maxdist` does.
> `seqdist` resolves the aliases to the
> canonical TraMineR names before dispatch, so the low-level spelling
> difference never reaches `seqdist` callers.

```python
from traminex import seqdef, seqdist
import numpy as np

X = seqdef(np.array([["A","B","A"],["B","B","A"]], dtype=object))

D_len = seqdist(X, method="OM", subm="TRATE", norm="length")
D_max = seqdist(X, method="OM", subm="TRATE", norm="max")
D_raw = seqdist(X, method="OM", subm="TRATE", norm=None)
```

## OM-family spell/warp methods (OMloc, TWED, OMslen, OMspell, OMstran)

These five methods accept the full TraMineR normalization vocabulary. They use the same
underlying transform as the classic family (TraMineR
`DistanceCalculator::normalizeDistance`, `src/distancecalculator.h`), selected by name:

| `norm` | Formula (raw `d`, per-pair max `M`, lengths `l1,l2`) | Notes |
|---|---|---|
| `None` / `'none'` | `d` | raw distance, unchanged |
| `'maxlength'` (alias `'length'`) | `d / max(l1, l2)` | can exceed 1 when `d` can exceed the longer length's indel budget |
| `'gmean'` | `1 - (M - d) / (2·√l1·√l2)` | Elzinga-style; not generally in `[0,1]` |
| `'maxdist'` (alias `'max'`) | `d / M` | in `[0, 1]` when `d` cannot exceed `M` |
| `'YujianBo'` | `2d / (d + M)` | in `[0, 1]` when `d` cannot exceed `M` (Yujian & Bo 2007) |
| `'auto'` | resolves to `'YujianBo'` | TraMineR's automatic choice for these methods |

`d == 0` always maps to `0`, so identical sequences normalize to `0` under every option.

The per-pair maximum `M` and lengths `l1,l2` are method-specific, transcribed from each
method's C++ `distance()` (see `traminex/distances/_normalization.py` for the exact
formulas and R source lines):

- **OMloc** (`OMVIdistance.cpp`): R's normalization indel is
  `indel = expcost·max(sm) + context` (`seqdist.R` line 723), `l = len·indel`, and
  `maxscost` is computed from that `indel` as below. An indel costs
  `expcost·maxscost + context·(sm[prev, state] + sm[next, state])/2`, at most
  `I = expcost·maxscost + context·max(sm)`; `M = |n−m|·I + min(m,n)·min(S, 2I)` (`2I`
  under `YujianBo`), `S` the largest off-diagonal cost. R's `M` prices the indels at
  `indel` (TX-REP-20261021-003, below). Because `maxscost` sets the local indels,
  OMloc's raw distance depends on `norm`, as in R: `YujianBo` fixes `maxscost` at
  `2·indel`, the other options use `min(S, 2·indel)`. Example: `subm=1`,
  `expcost=0.25` (so `context=0.5` and `indel=0.75`): `A` against `BBB` is 2.5 under
  `norm="none"`, but 2.75 behind `YujianBo` (`maxscost` 1.5 instead of 1), which
  normalizes to `2·2.75 / (2.75 + 3.5) = 0.88`. `gmean` keeps R's bound
  (TX-REP-20261022-002, below).
- **TWED** (`TWEDdistance.cpp`): `indel = 2·max(sm) + nu + h` (when not given);
  `M = |n−m|·(nu+h+c) + 2·(c+nu)·min(m,n)`, `c` the largest off-diagonal cost of the whole
  matrix, uncapped (at least `2·indel` under `YujianBo`); with an empty sequence
  `M = |n−m|·indel`, the distance R gives there; `l = len·indel`. R's `c` is `maxscost`,
  capped at `2·indel` (TX-REP-20261021-005, below). `gmean` keeps R's bound
  (TX-REP-20261022-002, below).
- **OMslen** (`OMvdistance.cpp`): the substitution cost matrix is `sm/2` for `link="mean"`
  (and the spell term is doubled) or `sm` for `link="gmean"`;
  `M = |n−m|·indel + f·maxscost·min(m,n)` with `f=2` for mean else `1`; `l = len·indel`.
- **OMspell** (`OMPerdistance.cpp`): normalization uses the **full** sequence lengths
  (not the DSS spell count); `l = len·indel`. R's `M = |n−m|·indel + maxscost·min(m,n)`
  leaves out the spell-duration costs (a spell of length `t` costs
  `indel + expcost·(t^tpow − 1)` to insert or delete), so R's `maxdist`, `gmean` and
  `YujianBo` values exceed 1, even at `expcost=0` when the substitution cost is below the
  indel. traminex's `M` is the largest cost, over the spell structures of the two
  lengths, of pairing spells in order and indelling the rest (TX-REP-20261021-004,
  below); with R's defaults (`CONSTANT` costs, `indel=1`, `expcost=1`, `tpow=1`) it is
  R's. Example: four sequences of length 5 over three states (`AABBC`, `ABBCC`,
  `CCBAA`, `BBCCA`), `subm="CONSTANT"`, `indel=1`, `expcost=10`: R's largest `maxdist`
  is 2.4 and `YujianBo` 1.411765; traminex's are 0.292683 and 0.45283, the bound 82
  being reached by `AAAAA` against `BBBBB`. `gmean` keeps R's bound
  (TX-REP-20261022-002, below). With `tpow` below 0 a spell's indel cost
  `indel + expcost·(t^tpow − 1)` falls as the spell lengthens and is negative once
  `expcost·(1 − t^tpow)` exceeds the indel; R then returns negative distances, which
  no bound covers, and traminex refuses such data (TX-REP-20261022-003). A zero spell
  cost stays valid.
- **OMstran** (`OMVIdistance.cpp`): after the transition augmentation, OM runs with a
  per-token indel vector; the normalization `indel` is `max(indel_vector)` and `maxscost`
  comes from the transformed substitution matrix; lengths use the augmented length.

`maxscost` is `min(max-off-diagonal substitution cost, 2·indel)`, except for `'YujianBo'`
where it is fixed at `2·indel` (TraMineR `OMdistance.cpp` lines 23-34).

```python
import numpy as np
import traminex as tx

X = tx.seqdef(np.array([["A","A","B","C"],["A","B","B","C"]], dtype=object))

tx.seqdist(X, method="OMspell", indel=1.0, subm=2.0, expcost=0.5, norm="YujianBo")
tx.seqdist(X, method="TWED", nu=0.5, h=1.0, subm=2.0, norm="maxdist")
tx.seqdist(X, method="OMslen", indel=1.0, subm=2.0, link="mean", norm="auto")  # -> YujianBo
```

## Asymmetric costs and the maximum possible distance

One deliberate difference from TraMineR (TX-AUD-20261019-008). `maxdist` divides a
distance by the maximum possible distance between sequences of the two lengths, and
`gmean` uses it too (`YujianBo` uses two indels in its place, as in R). For `OM`,
`OMslen`, `OMspell` and `TWED` that
bound is built from the largest substitution cost; R reads it from the upper triangle
of `sm` only (`OMdistance.cpp`), so with a larger cost below the diagonal the bound is
too small and normalized distances exceed 1 (`maxdist` reaches 3 on a four-state
example). traminex takes the largest off-diagonal cost of the whole matrix, still
capped at two indels as R caps it; with a symmetric matrix nothing changes. `gmean` of
`OMspell` and `TWED` keeps R's upper-triangle bound (TX-REP-20261022-002, below). Raw
distances are R's, and so are `OMstran`'s normalized values and TWED's border
value (`OMloc`'s bound is corrected for another cause, below). R's values are recorded by
`scripts/r_goldens/make_goldens_asym_bound.R`.

This removes one cause of normalized values above 1; the others, in `OMloc`, `OMspell`
and `TWED`, are corrected below. Normalized `OM` and `OMslen` distances stay within 1
over 200 random asymmetric matrices with indels from 0.25 to 3
(`tests/test_asymmetric_norm_bound.py`).

`TWED` (TX-REP-20261021-005). R caps the largest cost of its bound at two indels
(`TWEDdistance.cpp` line 153 with `OMdistance.cpp`'s `maxscost`), but TWED's matches
and warps pay `sm` directly, so with a small `indel` R's normalized values exceed 1
(`maxdist` 1.58 and `gmean` 3.91 on two six-sequence examples); and against an empty
sequence R sets the distance to `n·indel` (lines 144-149), above its bound (`maxdist`
1.667). traminex keeps R's form with the uncapped largest cost, which bounds the path
that matches along the shorter length and warps the rest, and uses `|n−m|·indel` with
an empty sequence. Raw distances, TWED's border value and the lengths `l` stay R's.

`OMspell` (TX-REP-20261021-004). R's bound leaves out the spell-duration costs
(`OMPerdistance.cpp` line 108), so R gives `OMspell` a `maxdist` of 4.86 with
`indel=0.5, expcost=1, tpow=2`, and 9 at `expcost=0` for `AAAAAAAAAA` against
`ABABABABAB` with substitution costs 0.1 and `indel=1`. traminex's bound is the largest
cost, over every spell structure of the two lengths, of pairing spells in order
(substituting, or deleting and inserting) and indelling the rest: with `p` and `q`
spells, `expcost·(g_p(m) + g_q(n)) + min(p, q)·min(S, 2·indel) + |p − q|·indel`, `g` the
largest sum of `t^tpow − 1` over `p` spell lengths (`2·indel` in place of
`min(S, 2·indel)` under `YujianBo`; for `tpow < 0` the duration terms are left out and
a pair costs at most `min(max(S, expcost), 2·indel)`). Raw distances and the lengths
`l` stay R's; R's values are recorded by `scripts/r_goldens/make_goldens_norm_bounds.R`.

`OMloc` (TX-REP-20261021-003). R's bound prices an indel at
`expcost * max(sm) + context` (`seqdist.R` line 723), but OMloc's local indel
costs `expcost * maxscost + context * (sm[prev, state] + sm[next, state]) / 2`
(`OMVIdistance.h`), up to `I = expcost * maxscost + context * max(sm)`. With
`expcost` below 0.5 R's normalized values therefore exceed 1, with its
`CONSTANT` costs too (`maxdist` 1.25 for A against BBB at `expcost=0.25`).
traminex's bound is `|n - m| * I + min(m, n) * min(S, 2I)` (`2I` under
`YujianBo`), `S` the largest off-diagonal cost: the cost of substituting along
the shorter length (or deleting and inserting) and indelling the rest, which is
always possible. Raw distances and the `gmean` lengths stay R's; R's values are
recorded by `scripts/r_goldens/make_goldens_norm_bounds.R`.

`gmean` of `OMloc`, `OMspell` and `TWED` keeps R's bound (TX-REP-20261022-002, owner
decision of 2026-10-22). Elzinga's `gmean`, `1 − (M − d)/(2·√l1·√l2)`, stays at least 0
only while `M − d ≤ 2·√(l1·l2)`; R's lengths `l = len·indel` ensure that for a bound close
to the true maximum (seqdist.Rd: "to avoid negative outcomes"), and the larger corrected
bounds above broke it: `OMloc`, `A` against `AAB` with `sm` 4, `expcost` 0 and `context`
0.5, gave −2.46 against R's 1. So `gmean` takes R's own bound, built from R's
normalization indel and R's `maxscost` (the upper triangle of `sm`, capped at two
indels): it is R's value given R's raw distance, also with asymmetric costs, and it can
exceed 1 where R's does. `maxdist` and `YujianBo` (and so `auto`) keep the corrected
bounds.

## INDELS costs and states no sequence contains

One deliberate difference from TraMineR (TX-AUD-20261019-002). `subm="INDELS"`
gives a state of frequency 0 -- an `alphabet=` state that no sequence contains --
the stand-in indel `.Machine$double.xmax`. That cost enters no distance between
the data's sequences, but R lets it set the normalization bound (and OMloc's
local costs), so every normalized OM, HAM, DHD, OMloc, OMslen, OMspell or TWED
distance becomes 0 or NaN (TWED's raw distances become infinite), and OMstran's
raw and normalized distances change, because OMstran rescales its costs by the
largest indel and substitution costs (`seqdist-OMstran.R` lines 73-77).
traminex keeps such a state out of the scale: its distances equal R's with that
state left out of the alphabet (for DHD, per position: the bound sums the largest cost between
the states present at each position). The costs `seqcost` and `seqsubm` return
stay R's, and a state that a compared sequence (an external `refseq` included)
contains keeps its cost. The cases are recorded from R by
`scripts/r_goldens/make_goldens_indels_unobserved.R`.

## Methods that refuse `norm`

`NMS`, `NMSMST`, and `SVRspell` do **not** accept `norm` — exactly as TraMineR's
`seqdist()` refuses it for them (`seqdist.R` lines 206-208). traminex raises a `ValueError`
in that case. Their only transform is the final `sqrt` applied in the R wrapper, which
traminex already replicates. `CHI2`/`EUCLID` accept only `norm="none"` (the default; also
`None`) or `"auto"` (`seqdist.R:314-315`), distinct from the OM-style names above; a
logical `norm` is TraMineR's deprecated spelling (`True` for `"auto"`, `False` for
`"none"`) and gives a `DeprecationWarning`. `docs/NEW_DISTANCE_METHODS.md` describes what
`"auto"` computes.
