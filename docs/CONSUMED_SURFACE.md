# Surface consumed by the sibling projects

The public API is the 226 names of `traminex/api_manifest.py`. The sibling
projects of the family (`ECOSYSTEM.md`) also import traminex objects by
module path, some of them outside the manifest. This page lists every such
import in their package sources (tests, benchmarks and scripts are left out),
read on 2026-10-18 (audit TX-AUD-20261010-048). They are a contract in
practice: renaming, moving or changing one of them breaks the consumer, so a
change to any of them is decided with the siblings in view (the siblings are
on standby; no sibling is changed from this repository).
`tests/test_consumed_surface.py` requires every path below to resolve.

flexscanx's package source imports nothing from traminex (it adapts code and
names traminex only in comments). A name imported from the root package that
is in the manifest is not listed: it is the public API.

| Module | Name | Consumer | In the manifest |
|---|---|---|---|
| `traminex.format` | `STS_to_SPELL` | py_seqplot: `py_seqplot/_series.py` | yes (imported by module path) |
| `traminex.sequences` | `STSList` | py_seqplot: `py_seqplot/_series.py` | yes (imported by module path) |
| `traminex.visualization.series` | `resolve_sort_order` | py_seqplot: `py_seqplot/_series.py` | no |
| `traminex.diss` | `dissassoc` | weightedcluster: `weightedcluster/proptree.py` | yes (imported by module path) |
| `traminex.diss` | `weighted_pairwise_sum` | weightedcluster: `weightedcluster/_fallback.py`, `weightedcluster/_utils.py`, `weightedcluster/quality.py` | yes (imported by module path) |
| `traminex.disstree` | `DissTree` | weightedcluster: `weightedcluster/seqtree.py` | no |
| `traminex.disstree` | `DissTreeNode` | weightedcluster: `weightedcluster/seqtree.py` | no |
| `traminex.disstree` | `DissTreeSplit` | weightedcluster: `weightedcluster/seqtree.py` | no |
| `traminex.disstree` | `disstree` | weightedcluster: `weightedcluster/proptree.py` | yes (imported by module path) |
| `traminex.distances` | `pairwise` | weightedcluster: `weightedcluster/_seqdist.py` | no |
| `traminex.distances` | `pairwise_between_sets` | weightedcluster: `weightedcluster/_seqdist.py` | no |
| `traminex.format` | `STSList` | weightedcluster: `weightedcluster/proptree.py` | yes (imported by module path) |
| `traminex.format` | `seqdef` | weightedcluster: `weightedcluster/proptree.py` | yes (imported by module path) |
| `traminex.seqformat` | `seqformat` | weightedcluster: `weightedcluster/seqclustname.py` | yes (imported by module path) |
| `traminex.seqtree` | `seqtree` | weightedcluster: `weightedcluster/proptree.py` | yes (imported by module path) |
| `traminex.visualization` | `state_legend` | weightedcluster: `weightedcluster/fuzzy.py` | yes (imported by module path) |
| `traminex.visualization.palette` | `ensure_cmap` | weightedcluster: `weightedcluster/fuzzy.py` | no |
| `traminex.visualization.series` | `select_index_rows` | weightedcluster: `weightedcluster/fuzzy.py` | no |
| `traminex.visualization.ticks` | `apply_time_axis` | weightedcluster: `weightedcluster/fuzzy.py` | no |

## Behaviour changes to listed objects

- **2026-10-23 — `pairwise`'s TWED warning for an indel vector has R's text (batch 16;
  TX-REP-20261023-008).** With `method="TWED"` and one indel cost per state,
  `distances.pairwise` (and `seqdist`) still use the largest cost and warn; the
  `UserWarning` now reads " [!!] indel vector not supported by the chosen method,
  max(indel) used instead!", R's message (`seqdist.R:301-305`), instead of "TWED does
  not support an indel vector; max(indel) is used". Recorded in `CHANGELOG.md`. Known
  consumer: weightedcluster's `_seqdist.py` passes a numeric indel only for OM and
  OMSPELL, so it does not meet this warning, and none of its warning filters names the
  old text. The sibling is not changed from this repository.
- **2026-10-23 — `select_index_rows` takes a unique prefix of `sample_meth` (batch 15;
  TX-REP-20261023-007).** `visualization.series.select_index_rows` keeps its name and
  signature; `sample_meth` also accepts a unique prefix of `"systematic"` or `"random"`
  (`"sys"`, `"r"`), as R's `plot.stslist` reads it with `match.arg`; other values are
  refused as before. Recorded in `CHANGELOG.md`. Known consumer: weightedcluster's
  `fuzzy.py` passes its own `sample_meth` through; full names behave as before. The
  sibling is not changed from this repository.
- **2026-10-23 — `JSD_TRANS`'s `eps` defaults to 0 (batch 15; review M-2 of
  TX-AUD-20261010-095).** `distances.pairwise` (and `seqdist`) keep their names and
  signatures; `method="JSD_TRANS"` without `eps` now gives the Jensen–Shannon
  divergence itself, where the default `1e-12` inside the logarithm made it up to about
  2e-12 smaller (`ln 2` minus 2.0e-12 for disjoint transitions). Recorded in
  `CHANGELOG.md`. Known consumers: none uses `JSD_TRANS` (weightedcluster's and
  py_seqplot's package sources name it nowhere).
- **2026-10-23 — `dissassoc` with a zero total weight (batch 15; TX-AUD-20261010-096).**
  `diss.dissassoc` keeps its name and signature. With weights that sum to 0 it returns
  R's floating-point results (the Total row's discrepancy and the statistics NaN), where
  Python's float division raised `ZeroDivisionError`. Recorded in `CHANGELOG.md`.
  Known consumer: weightedcluster's `proptree.py` calls `dissassoc` for its tree
  tests and adjustment; a call with all-zero weights now gives NaN statistics instead
  of an exception. Other inputs are unaffected. The sibling is not changed from this
  repository.
- **2026-10-23 — `STSList` refuses equal missing and void texts; `seqformat`'s
  STS output from a table has `seqdef`'s labels and colours (batch 14;
  TX-AUD-20261019-042, TX-AUD-20261019-035).** `sequences.STSList` (also
  imported as `format.STSList`) refuses an `nr` text equal to its `void` text
  ("'nr' and 'void' must be distinct character markers"), as `seqdef` does;
  it built such an object, whose missing and void cells could not be told
  apart. `seqformat` from a table to STS gives the states `seqdef`'s labels
  and colours, the missing marker `missing_color` and void white, where the
  result had none. Recorded in `CHANGELOG.md`. Known consumers: py_seqplot's
  `_series.py` rebuilds a derived object with `nr=str(seqdata.nr)` and
  `void=str(seqdata.void)`, which is unchanged for objects from `seqdef` and
  the builders (their markers are `*` and `%` by default); for a hand-built
  object whose markers are both `None`, both become the text `"None"`, and the
  rebuild that was accepted is now refused. weightedcluster's `proptree.py`
  checks the type and calls `seqdef`, which already refused equal markers, and
  its `seqclustname.py` converts an `STSList` to SPS, which the labels change
  does not touch. Both marker representations of an object now give the same
  `seqformat` output (STS, DSS and SPELL from an `STSList`) and the same
  `STS_to_SPELL` output with `right` (TX-REP-20261023-003): an object declaring the
  missing marker's text in its alphabet beside an external code had that marker
  written `<NR>`. py_seqplot calls `STS_to_SPELL` without `right`, which read the
  object this way already, and weightedcluster's SPS output was already the same.
  The siblings are not changed from this repository.
- **2026-10-23 — time ticks and labels (batch 13; TX-AUD-20261011-067).**
  `visualization.ticks.apply_time_axis` keeps its name and signature but now
  puts a tick at every `xtstep`-th position from the first (plus the last with
  `tick_last`), as R's plots do, where it showed at most 10 evenly spaced ticks;
  each tick gets its label, and, as R's `axis()`, a label that would come closer
  than the width of an "m" to the last one drawn is left out when the axis is
  drawn, so the labels shown depend on the figure's size. `max_labels` no longer
  limits the labels. Recorded in `CHANGELOG.md`. Known consumer: weightedcluster's
  `fuzzy.py` calls `apply_time_axis(ax, seqdata, length)` for its fuzzy index
  plot, which now shows R's tick grid (for example 1, 7, ..., 67 with
  `xtstep=6` on 72 positions, and every position with the default step) and
  only the labels that fit. The sibling is not changed from this repository.
- **2026-10-22 — normalized OMloc, OMspell and TWED distances, OMspell's spell costs
  and the NMS family (batch 10 and the owner's decision on gmean;
  TX-REP-20261021-003 to -005, TX-REP-20261022-001 to -003, TX-AUD-20261019-029,
  -030).** `pairwise` gives OMloc, OMspell and TWED `maxdist` and `YujianBo`
  values with a corrected bound, so they no longer exceed 1 and differ from R's
  where R's bound was too small; `gmean` keeps R's bound (owner decision).
  OMspell refuses data whose spell indel cost would be negative (`tpow` below 0
  with a large `expcost`). Raw distances change only in the last bits, now R's
  (spell-length powers, the NMS-family sums), and where R's NMS-family counts
  overflow, which now give finite distances. Recorded in `CHANGELOG.md`. Known
  consumer: weightedcluster's `_seqdist.py` maps its OMSPELL metric to
  `norm="maxlength"` or `"none"`, which the bound changes do not touch; it can
  see the last-bit raw changes and the new refusal.
- **2026-10-21 — `STS_to_SPELL` reads the object as it is (TX-AUD-20261019-014).**
  As TraMineR's internal `STS_to_SPELL`, it keeps a spell of trailing
  nonresponse cells (an object built with `seqdef(..., right=None)`); it used to
  re-run `seqdef(right="DEL")` first and drop them. `STS_to_SPELL(x, right="DEL")`
  and `seqformat(x, to="SPELL")` give the previous output. Objects without
  trailing nonresponse are unaffected. Recorded in `CHANGELOG.md`. Known consumer
  affected: py_seqplot's `_series.py` (the ggseqiplot data builder) calls
  `STS_to_SPELL(seqdata, with_missing=True)` to reproduce R's
  `seqformat(seqdata, to = "SPELL")`, which drops trailing missing cells; for
  such an object it now receives the trailing missing spell, and adding
  `right="DEL"` to that call keeps R's output. The sibling is not changed from
  this repository; the owner is told (batch-8 review I-4).
- **2026-10-21 — reference direction with asymmetric costs (TX-AUD-20261019-009).**
  `pairwise(ref_sets=)` and `pairwise_between_sets` return d(first-set row,
  second-set row) with an asymmetric substitution matrix on every route, as R's
  `refseq=list(a, b)` and `seqdist` do; the void route of `pairwise(ref_sets=)`
  and a second set repeating a row of the first gave the other direction.
  Symmetric costs are unaffected. Recorded in `CHANGELOG.md`.
- **2026-10-19 — default alphabet in R's order (owner decision).** `seqdef` and
  `seqformat` (and every other constructor) now order the states of an
  alphabet they derive as R's `seqdef` does: by numeric value when every state
  is an R number, R's collation otherwise. Before, the states kept their order
  of first appearance. A consumer that relied on that order (state codes,
  per-state rows, colours by position) gets R's order instead; passing
  `alphabet=` to `seqdef` keeps any specific order. `seqformat` takes no
  `alphabet=`: pass its result to `seqdef(x, alphabet=[...])`, which keeps
  every cell and marker. Recorded in `CHANGELOG.md`.
- **2026-10-17 — `STS_to_SPELL` labels spells with the object's ids (TX-AUD-20261010-017,
  `5b3938a0`; added here 2026-11-07).** As R's `STS_to_SPELL`, the `id` column holds the
  object's ids (R's row names) where it held row positions. Known consumer: py_seqplot's
  `_series.py` matched spells to rows by position, so against this traminex `ggseqiplot` of an
  object with ids drew a blank plot without error and `ggseqrplot`'s index panel collapsed
  (measured by the supervision repository on 2026-11-02). py_seqplot matches the ids since
  its release commit of 2026-11-07 (`f370c13f`).
- **2026-10-13 — `pairwise(method="DHD")` validates explicit costs (`9ac62e0d`; added here
  2026-11-07).** With explicit costs, DHD refuses non-finite values, as R's DHD stops.
  `pairwise(method="HAM")` and `seqdist` still accept an `Inf` cost for an unused state, as R's
  HAM does. Known consumer: weightedcluster's `_seqdist.py` routed HAM through DHD, so four of
  its tests failed against this traminex (measured 2026-11-02).
- **2026-09-22 — R's text for the asymmetric-cost warning (`df708c99`; added here
  2026-11-07).** The `UserWarning` for an asymmetric substitution matrix reads
  "Non-symmetric substitution cost matrix!", R's message. Known consumer: six of
  weightedcluster's tests matched the older wording (measured 2026-11-02).

These three were made after the siblings' last tested commit (`6e76440b`) and before this list
began; the supervision repository found them by running the siblings against the current
traminex. One more change reaches a consumer through the public API, so it is outside the
table: `dissmfacw` gives NaN for a term that adds no rank to the design (`3a3a5149`,
2026-10-17; a `DIVERGENCE` row of `docs/API_PARITY_DECISIONS.csv`), where R gives a finite
value of an arbitrary direction; weightedcluster's `clustassoc` is built on it.
