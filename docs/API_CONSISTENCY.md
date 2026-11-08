# API consistency and examples

This document summarizes recent API consistency improvements and shows examples for new/updated functions.

## 1) Missing-state keywords: two spellings

R's missing-state switch is `with.missing`. traminex spells it in one of two
ways, chosen function by function when each was ported rather than by one
rule. Use the keyword the function's signature shows: the other spelling
raises `TypeError` (section 1b), except in `seqsplot`, which hands extra
keywords to matplotlib, where it fails with `AttributeError`. The lists are the
complete public set and are checked against the signatures by
`tests/test_api_consistency_keywords.py`.

- `include_missing` (35 functions): `chronogram`, `cluster_chronograms`,
  `cluster_chronograms_with_representatives`, `seqalign`, `seqalignplot`,
  `seqdHplot`, `seqdss`, `seqdur`, `seqfplot`, `seqfpos`, `seqfposend`, `seqHtplot`,
  `seqibad`, `seqici`, `seqient`, `seqintegr`, `seqipos`, `seqistatd`, `seqivardur`,
  `seqivolatility`, `seqlogp`, `seqmaintokens`, `seqmeant`, `seqmodst`, `seqmsplot`,
  `seqmtplot`, `seqprecarity`, `seqprecstart`, `seqST`, `seqstatd`, `seqstatf`,
  `seqtab`, `seqtrate`, `state_distribution_plot`, `turbulence_exact`.
- `with_missing` (23 functions): `alphabet`, `seqasnum`, `seqcost`, `seqdiff`,
  `seqdomassoc`, `seqformat`, `seqidegrad`, `seqimplic`, `seqindic`, `seqindic_dyn`,
  `seqinsecurity`, `seqlength`, `seqMD`, `seqmpos`, `seqnum`, `seqpolyads`,
  `seqrep_grp`, `seqsplot`, `seqsubsn`, `seqsurv`, `seqsurvplot`, `seqtransn`,
  `STS_to_SPELL`.
- Both (1 function): `seqtabstocc` — it takes `include_missing` and accepts R's
  `with_missing` as a compatibility spelling; giving both is an error.

Examples:

- seqistatd(x, include_missing=False)
- seqstatd(x, include_missing=False, weights=...)
- seqtabstocc(x, include_missing=False)
- seqindic(x, with_missing=False)
- seqnum(x, with_missing=False)
- seqtrate(x, include_missing=False)
- seqdist(x, with_missing=True) — accepted through `**kwargs`, and
  **mandatory** when the sequences contain missing values: otherwise seqdist
  raises `ValueError: with_missing must be true when sequences contain missing
  values`. It appears in no signature or stub because seqdist is a `**kwargs`
  dispatcher; see its docstring's "Common options" list.

## 1b) Parameter naming: how R argument names map to Python

- These rules are mechanical and cover the large majority of ported arguments:
  - dots become underscores: R `with.missing` -> `with_missing`, `w.rank` -> `w_rank`.
  - camelCase becomes snake_case: R `dropMatrix` -> `drop_matrix`,
    `firstState` -> `first_state`.
  - a name that collides with a Python keyword gains a trailing underscore:
    R `except` -> `except_`.
  - R's leading `seqdata` argument is spelled `x` (87 of the 167 R exported
    functions with at least one formal lead with `seqdata`; 107 traminex exports
    lead with `x`). Twelve traminex functions keep the name `seqdata`.
  - R's replacement-function `value` argument has no Python equivalent by design.
- The renames below are *semantic* -- no rule above predicts them -- so they are
  listed individually. R name, the R formal, and the traminex parameter:

  | R function | R formal | traminex parameter |
  |---|---|---|
  | `seqfind` | `x` | `ref` |
  | `seqfind` | `y` | `target` |
  | `seqsubsn` | `DSS` | `dss` |
  | `gower_matrix` | `diss` | `D` |
  | `dissdomassoc` | `domdiss` | `domain_dists` |
  | `rowmode` | `except` | `except_` |
  | the functions of section 1's `include_missing` list | `with.missing` | `include_missing` |

- Plotting arguments `main` (10 functions) and `group` (8) are re-expressed by
  the plotting layer and are recorded as `RE_EXPRESSED` in
  `docs/API_PARITY_DECISIONS.csv`; they are owned there, not here.
- Passing an R spelling that traminex does not accept raises `TypeError`, so the
  failure is loud rather than silent.

## 1c) Result orientation: per-sequence indicators are 1-D

- R returns most per-sequence indicators as an `n x 1` **matrix** with a named
  column; traminex returns a **1-D array of length n**. Values, order and
  missing-state handling are unchanged — only the trailing axis differs.
- Measured on a 6-sequence fixture: `seqient`, `seqtransn`, `seqlength`,
  `seqici` and `seqST` are `6x1` in R and length-6 in traminex. `seqsubsn` is
  `n x 1` in R (column "Subseq.") and length-n in traminex; it returned `(n, 1)`
  until 2026-10-18 (TX-AUD-20261010-029), and
  `tests/test_per_sequence_indicator_shape.py` checks the shape of each. `seqivolatility`
  and `seqivardur` return a plain vector in R too, so those already agree.
- Consequence for porting: `result[i, 1]` in R is `result[i]` here, and
  `nrow(result)` is `len(result)`. Use `np.asarray(result).reshape(-1, 1)` if a
  column matrix is required downstream.
- This is a package-wide convention, recorded here once rather than in each
  decision row. Per-**state** results follow one of three rules, recorded in each
  function's row of `docs/API_PARITY_DECISIONS.csv` and checked by
  `tests/test_per_state_shape_and_selection.py` (the table modes by
  `tests/test_construction_items.py`; TX-AUD-20261019-070):
  - codebook padding: `seqstatd`'s Frequencies, `seqmeant` and `seqistatd` give one
    entry per code of the object's codebook, the missing and void markers included
    (zero when unused), where R gives one per state of `statl`;
  - R's state axis: `seqstatf` and `seqtrate` give the real states, plus the missing
    state when `include_missing=True`, never void, as R;
  - the table modes of `seqfpos` and `seqfposend` (traminex extensions; R's
    functions take one state) give the real states, plus the missing state when
    `include_missing=True`, never void (TX-REP-20261023-004 for `seqfposend`).

## 2) New conversions and round-trip symmetry

- STS_to_TSE(x, state_changes_only=True|False, time_col='time'): converts state sequences to a long table of (id, time, state), with times x.start + position. A void cell (the end of a shorter sequence) gets no record, so with state_changes_only=False (one record per cell) TSE_to_STS reads it back as void, except for a position no sequence reaches and a sequence of void cells only, which has no record; a change table does not read back. A missing cell keeps its record with the nr marker, which TSE_to_STS reads as a state of that name.
- toPersonPeriod(x, time_col='timestamp'): extract a person-period table; time_col lets you harmonize with TSE_to_STS.
- events_to_state_sequences(..., align_missing=False): when True, treat the fill state as missing (nr), aligning downstream stats/plots.

### Examples

- STS to TSE (per-time records) and back:

```python
from traminex import seqdef, STS_to_TSE, TSE_to_STS
import numpy as np
import pandas as pd

arr = np.array([["A","A","B","B"],["A","B","B","B"]], dtype=object)
x = seqdef(arr, start=1)
tse = STS_to_TSE(x, state_changes_only=False, time_col="time")
x2 = TSE_to_STS(tse, id_col="id", time_col="time", state_col="state")
```

- events to STS with blanks treated as missing:

```python
from traminex import seqecreate, events_to_state_sequences
import pandas as pd

df = pd.DataFrame({"id":["p1","p1","p1"],"timestamp":[2,4,5],"event":["X","Y","Z"]})
events = seqecreate(df, id_col="id", time_col="timestamp", event_col="event")
sts = events_to_state_sequences(events, fill="_", align_missing=True)
```

- toPersonPeriod with custom time column:

```python
from traminex import toPersonPeriod
pp = toPersonPeriod(x, time_col="time")
```

## 3) Distances normalization parameters

- For classic OM/LCS/LCP/RLCP distances via `seqdist`:
  - `norm` accepts the full TraMineR vocabulary plus two traminex aliases,
    case-insensitively: `None | 'none' | 'length' | 'max' | 'maxlength' |
    'gmean' | 'maxdist' | 'YujianBo' | 'auto'`. Every value is R-exact.
    - `'length'`: alias of `'maxlength'` (divide by the longer sequence
      length, for OM times the indel cost)
    - `'max'`: alias of `'maxdist'` (divide by the maximum possible distance
      given indel/subm)
    - `None` / `'none'`: no post-normalization
    - `'gmean'`, `'YujianBo'`, `'auto'`: TraMineR's own names, same semantics
  - All four classic methods accept the whole set, and `seqdist` accepts both
    aliases for every other method that takes a TraMineR `norm` as well.
  - The low-level `traminex.distances.pairwise` keeps the narrower classic
    vocabulary (`None | 'length' | 'max'`); see
    [DISTANCES_NORMALIZATION.md](DISTANCES_NORMALIZATION.qmd).

### Examples

- TRATE + OM with normalization:

```python
from traminex import seqdist
D_len = seqdist(x, method="OM", subm="TRATE", norm="length")
D_max = seqdist(x, method="OM", subm="TRATE", norm="max")
D_raw = seqdist(x, method="OM", subm="TRATE", norm=None)
```

## 4) Plotting ax passthrough

- All main plot functions accept an optional ax and return (fig, ax):
  - sequence_index_plot, state_distribution_plot, seqfrequency_plot, chronogram, transition_heatmap.
- This allows composing plots into existing figures/subplots.

## 5) Trees: prediction helper

- seqtree_predict(result, x, features='durations+transitions+summary') mirrors disstree_predict for feature-based trees.
- If feature names differ in order, they are aligned; if required features are missing, an error is raised.

## 6) Package, typing, and API identity

- `traminex.api_manifest.PUBLIC_API`, runtime `traminex.__all__`, and the root
  review stub are the same ordered 226-name list.
- The package version is 0.2.2 and is derived from installed metadata.
- Bare `mypy` checks the complete configured package module set. Installed wheels expose inline
  annotations through `traminex/py.typed`; the separate `stubs/` tree is a
  bidirectional review gate and is not installed as a competing package.
- Without the `plot` extra (matplotlib), the plotting names stay in
  `traminex.__all__`, because the manifest does not vary by environment, but
  they are left out of `dir(traminex)`, so `help(traminex)` and
  `inspect.getmembers(traminex)` work on a core install. Accessing one raises
  `ImportError` naming the extra to install. `from traminex import *` walks
  `__all__` and therefore raises that `ImportError` on a core install; import
  the names you need explicitly instead.
- Five public functions share their name with the submodule that defines
  them: `disstree`, `seqdiff`, `seqemlt`, `seqformat` and `seqtree`.
  `traminex.seqtree` is always the function, even after the submodule has been
  imported. Python compiles `import traminex.seqtree as st` to the same
  attribute lookup as `from traminex import seqtree`, so `st` is the
  function too, and `st.seqtree_fit` raises `AttributeError`. To reach other
  names in such a submodule, import them from it
  (`from traminex.seqtree import seqtree_fit`) or use
  `importlib.import_module("traminex.seqtree")`.
