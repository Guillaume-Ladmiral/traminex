# Visualization-series consumer contract

This document freezes the supported boundary between `traminex` and downstream
visualization libraries, especially the sibling `py_seqplot`. It describes
analytical data and metadata, not a renderer ABI. `traminex` remains the source
of numerical truth; consumers may reshape these outputs but must not recompute
sequence statistics.

## Dependency direction and public boundary

The dependency is one-way: `py_seqplot` may depend on `traminex`; `traminex`
does not import or require `py_seqplot`. Supported consumers use top-level
public functions, `traminex.STSList`, and the public
`traminex.visualization.series.resolve_sort_order` helper. Names beginning with
an underscore are private and are not covered by compatibility guarantees.

Read-only inspection on 2026-07-15 found these production imports in
`py_ggplotseq` (now `py_seqplot`): `STSList`, `seqdef`, `STS_to_SPELL`, `resolve_sort_order`,
`seqdss`, `seqtrate`, `seqstatd`, `seqmeant`, `seqmodst`, `seqtab`, `seqrep`, and
`seqrf`. Private visualization helpers, including `_wtd_boxplot_stats`, are not
part of this contract and must not be imported by downstream production code.
External grammar-of-graphics renderers own renderer-owned drawing summaries
over the public `seqrf` box inputs described below.

## `STSList` contract

An `STSList` is a categorical state-trajectory matrix:

- `data` is an integer NumPy array of shape `(N, L)`: rows are sequences and
  columns are successive time positions. Row and time order are stable.
- Codes are zero-based indexes into `alphabet`. `labels` and `cpal`, when
  present, use the same state order. Consumers must not alphabetically reorder
  any of them.
- `nr_code` identifies nonresponse/missing and `void_code` identifies structural
  void. `missing_color` styles missing values; structural void is never a
  substantive state.
- `weights`, when present, has length `N`; `ids`, when present, has length `N`.
- `start`, `xtlab`, `xtstep`, and `tick_last` define the time-axis metadata.
  Derived sequence objects preserve applicable metadata unless their time axis
  necessarily changes.
- Public analytical functions and series helpers do not mutate the caller.

Consumers should require an actual `STSList`, not duck-type arbitrary matrices.
Malformed matrix dimensions, IDs, weights, palettes, and analytical arguments
must be refused rather than silently repaired.

## Verified analytical outputs

### Sequence definition, spells, and durations

`seqdef` is the canonical constructor. `STS_to_SPELL` returns spell rows without
changing sequence order. `seqdss` returns an `STSList` with the same `N`, state
order, weights, IDs, and special codes; its columns are distinct-successive
states and right padding is structural void.

### Transition rates

`seqtrate(..., as_dataframe=False)` returns a `(K, K)` float array over the
included states. Row `i` is the origin and column `j` the destination. Both axes
use state order. Rows are denominator-normalized over eligible transitions;
zero-denominator rows remain defined by the documented analytical routine.

### State distribution and entropy

`seqstatd` returns a `SeqStatDResult` mapping:

- `Frequencies`: `(L, K)` — time on rows, state on columns;
- `ValidStates`: `(L,)` — eligible weighted denominator at each time;
- `Entropy`: `(L,)` — normalized transversal entropy.

Structural void is always excluded from denominators. Missing is excluded unless
`include_missing=True`. With weights enabled, attached or explicit weights are
used; `weighted=False` uses unit weights. Excluded state columns remain aligned
with the `STSList` codebook and carry zero frequency, so consumers select columns
by `nr_code`/`void_code` rather than assuming a shortened matrix.

### Mean, modal, and frequency summaries

- `seqmeant(..., serr=True)` is a state-ordered DataFrame with columns
  `Mean`, `Var`, `Stdev`, and `SE`.
- `seqmodst` returns a one-row `STSList` of length `L`; `frequencies_` is the
  aligned modal share, while `weighted_` and `nbseq_` describe its denominator.
- `seqtab` returns `sequence`, `count`, `prop`, and `cumprop`, ordered by the
  routine's stable frequency rules. Its R-compatible default is the ten highest
  frequency ranks in SPS form; `idxs=0` requests every distinct sequence and
  `format="STS"` requests full state strings. `DataFrame.attrs` carries aligned
  `freq`, `idxs`, `idxf`, `source_indices`, and `representatives` metadata.
  Consumers that need state codes must use `representatives`/`source_indices`
  instead of parsing or matching the display string.

### Representatives, relative-frequency groups, and box inputs

`seqrep` returns representative code rows in `(R, L)` orientation plus their
statistics. `seqrf` retains its legacy array-valued `seqtoplot` in `(G, L)`,
zero-based `medoids`, group `heights`, and the full public `rf` mapping. It also
returns exactly one metadata-complete sequence field, `medoid_sequences`, as an
`STSList`. This field contains the medoid rows in output order, preserves the
source codebook, labels, palette, IDs, time-axis, missing/void, issue, and other
rendering metadata, and carries the RF group heights as its row weights. Its
arrays and mutable metadata are copies: mutating it does not mutate the source
`STSList`, `seqtoplot`, or the public RF arrays.

For each group, `index_list`, `dist_list`, and `weights_list` have equal length
and corresponding order. Those three lists are the supported weighted-box
inputs; a renderer may reshape them or compute renderer-owned drawing summaries
in its own R-golden-tested layer. The private matplotlib helper
`_wtd_boxplot_stats` is not a public cross-project API.

Distance inputs must be finite square matrices representing the same `N` rows
as the sequence object. Unsupported grouping methods, sort keys, invalid
weights, and malformed distance shapes must raise instead of producing an
invented visualization series.

## Change control

Changing any orientation, ordering, denominator, missing/void rule, metadata
field, or refusal behavior above requires an explicit compatibility decision and
consumer-contract test update. New plot styles do not belong here. Proposed
sibling migrations remain proposals until the maintainer schedules them.
