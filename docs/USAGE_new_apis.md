# traminex: Quick usage for core APIs

This guide summarizes the main utilities currently exposed in traminex.
It also covers discrepancy analysis, domain association, the exploratory
`SeqMD` container, and representative sequences.

## Distance computation (seqdist)

```python
import pandas as pd
import traminex as tx

x = tx.seqdef(pd.DataFrame({
  "t1": ["A","A","B"],
  "t2": ["B","A","B"],
  "t3": ["B","B","A"],
}))

# Optimal matching (OM)
D = tx.seqdist(x, method="OM", indel=1.0, subm=2.0, norm="length")

# Hamming
Dh = tx.seqdist(x, method="HAM")
```

Backend selection is process-wide. Leave `TRAMINEX_BINDING` unset for automatic
pybind11-or-Python selection, set it to `python` for the portable path, or set it
to `pybind11` to require the compiled backend. Those three values are its whole
domain. A required backend that is not installed raises an error instead of
silently falling back.

The Numba kernels are a **third** backend and are orthogonal to
`TRAMINEX_BINDING`: they activate on their own whenever `numba` is importable,
including under `TRAMINEX_BINDING=python`. Select them with `TRAMINEX_NUMBA`,
which accepts `auto` (the default: use them when available) or `off` (use the
reference NumPy kernels). `TRAMINEX_NUMBA=off` is the supported way to reproduce
a result on the reference path or to compare the two backends.

Both variables are read **once per process, when the distance module
`traminex.distances` is first imported**. `import traminex` alone does not
import it, because public names resolve lazily; the first use of a distance
function, or of any function that computes distances, does. From then on,
assigning to `os.environ` has no effect. Set both variables in the environment
before the process starts (or in the command that launches it) rather than
relying on when that first import happens.

Core usage needs only NumPy, pandas, and SciPy. Install `traminex[perf]` for
Numba, `traminex[plot]` for the matplotlib layer, or `traminex[tree]` plus the
system Graphviz `dot` executable for tree rendering.

## Spell conversions

- TSE_to_STS(long): strict id/time/state -> STS
  Accepts only sequence metadata kwargs supported by `seqdef_from_long`; spell-specific kwargs such as `fill` are rejected.
- HSPELL_to_STS: supports vertical and horizontal spell inputs
- FCE_to_TSE: fixed-column events to TSE
- toPersonPeriod: STS -> person-period long

```python
import pandas as pd
import traminex as tx

# Vertical spells
spells = pd.DataFrame({
  "id": [1,1,2],
  "start": [1,3,2],
  "end": [2,4,3],
  "state": ["A","B","A"],
})
X = tx.HSPELL_to_STS(spells, id_col="id", start_col="start", end_col="end", state_col="state")

# Horizontal spells
hspell = pd.DataFrame({
  "id": ["p1","p2"],
  "begin1": [1, 1], "end1": [3, 2], "status1": ["A","A"],
  "begin2": [2, 3], "end2": [4, 4], "status2": ["B","B"],
})
Xh = tx.HSPELL_to_STS(hspell, id_col="id",
  begin=["begin1","begin2"], end=["end1","end2"], status=["status1","status2"],
  overwrite=True, fillblanks="_"
)

# Align to age if you have a birth column
hspell_age = hspell.assign(birth=[1990, 1992])
Xa = tx.HSPELL_to_STS(hspell_age, id_col="id",
  begin=["begin1","begin2"], end=["end1","end2"], status=["status1","status2"],
  pvar="birth"
)
```

### Overlap options and interval closure

When horizontal spells overlap, you can choose how overlaps resolve with `overlap_priority` and how intervals are closed with `interval_closed`:
- overlap_priority: `"later"` (default), `"earlier"`, `"longer"`, or `"shorter"`
- interval_closed: `"both"` (default), `"left"` (end excluded), or `"right"` (start excluded)

```python
import pandas as pd
import traminex as tx

h = pd.DataFrame({
  "id": ["p1"],
  "begin1": [1], "end1": [3], "status1": ["A"],
  "begin2": [2], "end2": [4], "status2": ["B"],
})
# Later spell B overwrites overlap, both ends inclusive
X_later = tx.HSPELL_to_STS(h, id_col="id",
  begin=["begin1","begin2"], end=["end1","end2"], status=["status1","status2"],
  overlap_priority="later", interval_closed="both"
)
# Earlier spell A keeps priority
X_earlier = tx.HSPELL_to_STS(h, id_col="id",
  begin=["begin1","begin2"], end=["end1","end2"], status=["status1","status2"],
  overlap_priority="earlier"
)
# Visualize
_ = tx.sequence_index_plot(X_later)
_ = tx.sequence_index_plot(X_earlier)
```

## Event subsequences (minimal pipeline)

```python
import pandas as pd
import traminex as tx

# TSE events
tse = pd.DataFrame({
  "id": ["i1","i1","i2"],
  "timestamp": [1,2,1],
  "event": ["A","B","A"],
})
e = tx.seqecreate(tse)
out = tx.seqefsub(e, min_support=1, max_k=2)
ct = tx.seqeconstraint(max_gap=2, count_method="COBJ")
out2 = tx.seqefsub(e, min_support=1, max_k=3, constraint=ct)
out.data.head(), out2.data.head()
```

## Multi-domain sequences (SeqMD)

```python
A = tx.seqdef(pd.DataFrame({"t1": ["A","A","B"], "t2": ["B","A","B"]}))
B = tx.seqdef(pd.DataFrame({"t1": ["X","Y","X"], "t2": ["Y","Y","X"]}))
obj = tx.SeqMD([A,B], names=["edu","work"])  # exploratory container
# Additive combination of per-channel distances
D_md = tx.seqMD_distances(obj, method="HAM", weights=[0.5, 1.5])
```

## Discrepancy analysis: seqdiff

```python
import numpy as np
from traminex import seqdiff
X = tx.seqdef(pd.DataFrame({
  "t1": ["A","A","B","B"],
  "t2": ["B","B","A","A"],
  "t3": ["A","A","A","A"],
}))
group = np.array([0,1,0,1])
out = seqdiff(X, group, method="LCS")
out["stat"].head(), out["discrepancy"].head()
```

## Domain association summary

```python
from traminex.diss import dissdomassoc
D1 = tx.seqdist(A, method="HAM")
D2 = tx.seqdist(B, method="HAM")
res = dissdomassoc([D1, D2], jointdiss=D_md)
```

## Representative sequences (seqrep)

```python
from traminex.cluster import seqrep
# Neighbourhood-density criterion. Without `diss`, seqrep computes the distance
# matrix itself from the seqdist arguments it forwards, as R's seqrep does; the
# OM family needs an explicit substitution cost (`subm`), as R needs `sm`.
res1 = seqrep(X, criterion="density", coverage=0.5, method="OM", indel=1.0, subm=2.0)
# Distance-based medoids (three representatives)
D = tx.seqdist(X, method="OM", indel=1.0, subm=2.0)
res2 = seqrep(X, diss=D, nrep=3, criterion="dist")
```

## Plotting and time axis

- sequence_index_plot: index plot
- state_distribution_plot: state distribution plot
- chronogram: chronogram

```python
import matplotlib.pyplot as plt
import pandas as pd
import traminex as tx

df = pd.DataFrame({
  "t1": ["A","B"],
  "t2": ["B","B"],
  "t3": ["A","B"],
  "t4": ["A","A"],
  "t5": ["B","A"],
})
X = tx.seqdef(df, xtstep=2, start=2000)
fig, ax = tx.sequence_index_plot(X)
fig2, ax2 = tx.state_distribution_plot(X)
fig3, ax3 = tx.chronogram(X)
plt.close('all')
```

## Dissimilarity helpers (minimal)

```python
import numpy as np
from traminex.diss import dissassoc, disscenter, dissvar

D = np.array([[0,1,2],[1,0,1],[2,1,0]], dtype=float)
labels = np.array([0,0,1])

print(dissassoc(D, labels))
C = disscenter(D)
print(dissvar(D, labels))
```
