# Localized optimal matching: three distinct interfaces

traminex has three related interfaces whose names and cost models are easy to
confuse. Only the first is TraMineR's `OMloc` method.

## 1. TraMineR-style `seqdist(method="OMloc")`

Use this interface for parity with TraMineR 2.2-14:

```python
import numpy as np
import traminex as tx

x = tx.seqdef(
    np.array(
        [
            ["A", "A", "B", "C"],
            ["A", "B", "B", "C"],
            ["C", "B", "A", "A"],
        ]
    )
)

d = tx.seqdist(
    x,
    method="OMloc",
    subm=2.0,
    expcost=0.25,
    context=0.5,
    norm="auto",
)
```

`subm` is the state substitution cost, `expcost` controls the base insertion
and deletion cost, and `context` controls the adjacent-state contribution.
Both `expcost` and `context` must be non-negative. When `context` is omitted it
defaults to `1 - 2 * expcost`, as in TraMineR. `norm` accepts the TraMineR
normalization modes supported by `seqdist`, including `none`, `auto`,
`maxlength`, `gmean`, `maxdist`, and `YujianBo`.

This path supports `STSList` missing/void handling, ragged sequences, and the
`seqdist` reference forms. It refuses empty sequences, matching the R method.
Its full ragged, missing, normalization, refusal, and external-reference
contract is generated from the R oracle and executed in
`tests/test_r_goldens_distances_void.py`.

Implementation: `traminex/distances/__init__.py` dispatches to the TraMineR
recurrence in `traminex/distances/python_kernels.py`.

## 2. Python convenience function `omloc()`

`traminex.omloc()` is a separate Python API. It derives position-specific indel
costs from a dense integer-coded sequence matrix and optionally derives
position-specific substitution costs, then calls regular OM:

```python
seqmat = np.array(
    [
        [0, 0, 1, 2],
        [0, 1, 1, 2],
        [2, 1, 0, 0],
    ],
    dtype=np.int32,
)

d = tx.omloc(
    seqmat,
    indel_method="change",       # or "entropy"
    indel_base=1.0,
    indel_scale=0.8,
    indel_smooth=3,
    subm_timevarying=True,
    subm_base="TRATE",           # also "CONSTANT", scalar, 2-D, or 3-D array
    norm="length",
    condensed=False,
)
```

This convenience function is not an alias for
`seqdist(method="OMloc")`: it uses the regular OM recurrence with
`indel_by_pos` and, when requested, `subm_by_pos`. It accepts a dense,
equal-width matrix with at least one row, positive width, and non-negative
integer state codes; the inferred state count is `max(code) + 1`, so unused
intermediate codes are allowed. Boolean/float codes, negative codes, and empty
dimensions are refused. A global substitution array must be `(K, K)`; a
time-varying call accepts `(K, K)` for replication or `(L, K, K)` explicitly.
`trate_weights` is consumed only with `subm_base="TRATE"`.

This interface does not provide `STSList` void/ragged normalization or
`seqdist` reference-set forms. Its focused behavior and validation tests are in
`tests/test_omloc.py`.

Implementation: `traminex/distances/omloc.py`, with cost construction in
`traminex/distances/costs.py`.

## 3. Position-dependent costs on regular OM

Call regular OM directly when the position profiles are already available:

```python
length = seqmat.shape[1]
states = int(seqmat.max()) + 1

indel_by_pos = np.array([1.0, 0.6, 0.8, 1.0])
subm_by_pos = np.full((length, states, states), 2.0)
for matrix in subm_by_pos:
    np.fill_diagonal(matrix, 0.0)

d = tx.seqdist(
    seqmat,
    method="OM",
    indel=1.0,
    subm=2.0,
    indel_by_pos=indel_by_pos,
    subm_by_pos=subm_by_pos,
    norm="length",
)
```

`subm_by_pos` has shape `(L, K, K)`. At dynamic-programming cell `(i, j)`,
the substitution slice is indexed by `min(i - 1, j - 1)`.
`indel_by_pos` supplies the position cost used by the regular OM recurrence.
Its values must be finite and strictly positive. Every `subm_by_pos` slice
must be finite, non-negative, and zero-diagonal.
For `norm="max"`, the denominator uses the maximum across all position-specific
substitution costs and the maximum position-specific indel cost. `seqdist`'s
TraMineR norms (`maxlength`, `gmean`, `maxdist`, `YujianBo`) use the same
maxima wherever TraMineR's formulas use the indel and the largest
substitution cost. `gmean`, `maxdist` and `YujianBo` therefore stay within `[0, 1]`:
their bound prices each position at the largest substitution cost or the
largest pair of indels, and an alignment paying one of the two at each position
always exists. `maxlength` (and `length`, and `auto`, which resolves to
`maxlength` for OM) divides by the length times the largest indel only, so it
can exceed 1 (up to 2), in TraMineR as here: TraMineR's own OM gives 2 for
`AAAA` against `BBBB` with `CONSTANT` costs and `indel=1`. The example above
gives 1.2 for its second and third sequences, and its costs give 1.7 for `AAAA`
against `BBBB`.
Position-dependent regular-OM costs are a traminex extension and are refused
for ragged rows because there is no single unambiguous shared position axis.

Constant-slice equivalence and local-indel behavior are guarded by
`tests/test_om_subm_by_pos.py`, `tests/test_om_extensions.py`, and
`tests/test_local_indel.py`. The recurrence lives in
`traminex/distances/python_kernels.py`.
