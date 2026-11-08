# New distance & cost methods (TraMineR parity)

These methods close the gap between `traminex` and the TraMineR `seqdist()` / `seqcost()`
catalogue. The implementations below follow the TraMineR C++/R reference. **Default and
risk-weighted option routes are pinned** against live TraMineR 2.2-14 goldens
(`tests/test_r_goldens_distances.py`, `tests/test_nms_parity.py`, and
`tests/test_om_normalization.py`); this does not mean every possible parameter combination
has been exhaustively golden-tested. Use them through `traminex.seqdist(x, method=..., **kwargs)`
and `traminex.seqsubm(seqmat, method=..., **kwargs)`.

## Distance methods (`seqdist`)

| Method | Description | Key kwargs | R reference |
|---|---|---|---|
| `EUCLID` | Euclidean distance on per-interval state-distribution vectors | `step`, `overlap`, `norm` (`"none"`/`"auto"`), `global_pdotj`; `weights` is accepted and ignored, as in R | `seqdist-CHI2.R` (`euclid=TRUE`) |
| `CHI2` | pdotj-weighted chi-squared on per-interval state distributions | `step` (default 1), `overlap`, `norm` (`"none"`/`"auto"`), `weights`, `global_pdotj` | `seqdist-CHI2.R`, `chisq.cpp` |
| `TWED` | Time-Warp Edit Distance (Marteau) | `nu` (required, >0), `h` (≥0, default 0.5), `subm` | `TWEDdistance.cpp` |
| `OMslen` | OM weighted by spell lengths (position weight = dur^(−h)) | `indel`, `subm`, `link` (`mean`/`gmean`), `h` (default 0.5) | `OMvdistance.cpp` |
| `OMspell` | OM on the DSS with spell-duration costs (durations dur^tpow − 1) | `indel`, `subm`, `expcost` (timecost), `tpow` (default 1) | `OMPerdistance.cpp` |
| `OMstran` | OM on transition-augmented sequences | `otto` (required, in [0,1]), `transindel` (`constant`/`prob`/`subcost`), `previous`, `add_column`, `indel`, `subm` | `seqdist-OMstran.R` |
| `NMS` | Number of Matching Subsequences (Elzinga) | `kweights`, `prox` (soft-match matrix) | `NMSdistance.cpp` (identity/default); `NMSMSTSoftdistanceII.cpp` (non-identity `prox`) |
| `NMSMST` | NMS weighted by Minimum Shared Time | `kweights`, `tpow` | `NMSMSTdistance.cpp` |
| `SVRspell` | Subsequence Vectorial Representation (spell, duration-weighted, soft) | `kweights`, `tpow`, `prox` | `NMSDURSoftdistance.cpp` |

Notes:

- **CHI2** reproduces the Studer–Ritschard chi-squared exactly: with `a` and `b` the two
  sequences' state *proportions* in an interval and `p` the overall (weighted) state
  proportions `pdotj` there, the distance is `√(Σ_intervals Σ_states (a − b)² / p)`;
  **EUCLID** is the same with `p = 1`. The default `step` is **1** (one interval per
  position), matching TraMineR. `norm` is `"none"` (the default; also `None`) or
  `"auto"` (`seqdist.R:314-315`), distinct from the OM-style `norm="length"/"max"`.
  `"auto"` divides each interval's term by its largest possible value — `1/p₁ + 1/p₂`
  for CHI2, `p₁` and `p₂` the interval's two smallest positive overall proportions (no
  division when only one state occurs there), and 2 for EUCLID — and the distance by
  `√(number of intervals)` (`seqdist-CHI2.R:98-115, 173`). A logical `norm` is
  TraMineR's deprecated spelling (`True` for `"auto"`, `False` for `"none"`) and gives a
  `DeprecationWarning`. EUCLID ignores `weights`, as TraMineR does: its `seqdist` sets
  `weighted <- FALSE` for EUCLID (`seqdist-CHI2.R:6-8`).
- **OMslen/OMspell/OMstran/TWED** now support TraMineR-parity normalisation via
  `norm=` (`maxlength`/`gmean`/`maxdist`/`YujianBo`/`auto`, with `length`/`max` accepted as
  aliases; `auto` resolves to `YujianBo`). See `docs/DISTANCES_NORMALIZATION.md` for the
  per-method maximum-distance formulas and R source lines, and
  `tests/test_om_normalization.py` for the live-R goldens. `norm=None` returns the raw
  distance unchanged.
- **NMS/NMSMST/SVRspell** do **not** accept `norm` — TraMineR's `seqdist()` refuses it for
  them (`seqdist.R` lines 206-208); traminex raises a `ValueError` to match. Their only transform
  is the final `sqrt` applied in the R wrapper, which traminex replicates.
- Default and risk-weighted option routes are pinned against live TraMineR 2.2-14 goldens,
  including non-identity NMS proximity, in `tests/test_r_goldens_distances.py`,
  `tests/test_nms_parity.py`, and `tests/test_om_normalization.py`.

## Substitution-cost methods (`seqsubm`)

| Method | Description | Key kwargs | R reference |
|---|---|---|---|
| `FEATURES` | Gower dissimilarity over a per-state feature table | `state_features` (K×F), `feature_weights`, `feature_types`, `ordinal_levels` | `seqcost.R` (`daisy(metric="gower")`) |
| `FUTURE` | Chi-squared distance on common-future transition rates | `weights` | `seqcost.R` (`chisqdista`) |

`traminex.distances.costs.gower_distance(features, *, weights, types, ordinal_levels)` is also
available directly. Supported feature types: `numeric`/`interval`, `ordinal`, `nominal`/`factor`,
and `asymm` (asymmetric binary, joint-absence excluded).

The raw-array FEATURES extension drops constant numeric columns and returns
zero for pairs with no comparable positive-weight feature. The STSList
`seqcost`/`seqsubm` FEATURES path follows `cluster::daisy`: constant columns
remain in the denominator, and incomparable pairs return NaN. Raw-array Gower
calculations are done in double precision on every platform: each numeric
column, and each pair's weights, are scaled by exact powers of two, which
preserves very small effective weights, differences between nearby large
values and ranges spanning the whole double line without depending on the
platform's `long double`.
Missing feature values are excluded for every supported feature type. Observed
numeric values must be finite, and asymmetric binary observations must be 0 or 1.
These raw-array validation rules do not change the separate R-faithful STSList
FEATURES implementation.

## Examples

```python
import numpy as np
import traminex as tx

seqmat = np.array([[0, 0, 1, 1, 2], [0, 1, 1, 2, 2], [2, 2, 1, 0, 0]], dtype=np.int32)

# Time-warp edit distance
D_twed = tx.seqdist(seqmat, method="TWED", nu=0.5, h=1.0, subm=2.0)

# Spell-based OM
D_slen = tx.seqdist(seqmat, method="OMslen", indel=1.0, subm=2.0, link="mean")
D_spell = tx.seqdist(seqmat, method="OMspell", indel=1.0, subm=2.0, expcost=0.5)

# Transition-augmented OM
D_stran = tx.seqdist(
    seqmat, method="OMstran", indel=1.0, subm=2.0, otto=0.5, transindel="constant"
)

# Subsequence-vector family
D_nms = tx.seqdist(seqmat, method="NMS")
D_mst = tx.seqdist(seqmat, method="NMSMST")
D_svr = tx.seqdist(seqmat, method="SVRspell", tpow=1.0)

# Chi-squared / Euclidean on state distributions
D_chi2 = tx.seqdist(seqmat, method="CHI2", step=1, norm="auto")
D_euc = tx.seqdist(seqmat, method="EUCLID", step=2)

# Feature-based (Gower) substitution costs
state_features = np.array([[0.0, "x"], [1.0, "y"], [2.0, "x"]], dtype=object)
sm = tx.seqsubm(seqmat, method="FEATURES", state_features=state_features,
                feature_types=["numeric", "nominal"])
```
