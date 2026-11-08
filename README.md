# traminex

Python package for sequence analysis, inspired by TraMineR and TraMineRextras.

Current package version: `0.2.3`.

> **Project family.** traminex is the **foundation** of a four-project family.
> `pyweightedcluster` owns sequence clustering, `py_seqplot` owns the dedicated
> plotnine visualization layer, and `flexscanx` owns corrected spatial scan
> statistics while consuming more of the foundation over time. See
> [ECOSYSTEM.md](ECOSYSTEM.md) for the charter, [DUPLICATES.md](DUPLICATES.md)
> for the duplicate inventory, and [MIGRATION_PLAN.md](MIGRATION_PLAN.md) for the
> convergence plan. Any AI assistant should read [AGENTS.md](AGENTS.md) first.

## Install a release

traminex is not on PyPI. Each release is published on the project's GitHub
[releases page](https://github.com/Guillaume-Ladmiral/traminex/releases), with
its packages and a `SHA256SUMS` file. Python 3.11 or 3.12 is required. Install
the pure-Python wheel of the current release; pip takes NumPy, pandas and SciPy
from PyPI:

```bash
python -m pip install https://github.com/Guillaume-Ladmiral/traminex/releases/download/v0.2.3/traminex-0.2.3-py3-none-any.whl
```

To add an optional profile (`plot`, `perf`, `extras` or `tree`, described in the
next section), put its name before the address:

```bash
python -m pip install "traminex[plot] @ https://github.com/Guillaume-Ladmiral/traminex/releases/download/v0.2.3/traminex-0.2.3-py3-none-any.whl"
```

On Linux x86-64 with CPython 3.12 you can install the compiled wheel instead,
`traminex-0.2.3-cp312-cp312-linux_x86_64.whl` from the same release; it adds the
optional pybind11/OpenMP distance backend (see "Packaging and optional compiled
kernels" below). Both wheels need the ICU library, and tree plots need GraphViz;
pip cannot install either, and the next section says how to get them. To check
the files you downloaded, save `SHA256SUMS` from the release beside them and run
`sha256sum -c --ignore-missing SHA256SUMS`.

## Install (development)

Development targets Ubuntu 24.04 and supports Python 3.11 and 3.12. Use the
family-standard Python 3.12 miniforge environment:

```bash
~/miniforge3/envs/datasci/bin/python -m pip install -e '.[plot,perf,extras,tree,dev]'
```

The core install needs only NumPy, pandas, and SciPy. Optional profiles are
`perf` (Numba), `plot` (matplotlib/seaborn), and `extras` (scikit-learn).

The `tree` profile installs nothing. Tree rendering
(`disstreedisplay`/`seqtreedisplay`) runs the **system GraphViz `dot`
executable** as a subprocess, which pip cannot install; supply it from your
package manager (`sudo apt install graphviz` on Ubuntu 24.04) and either put it
on `PATH` or point `GVPATH` at it. traminex does not use the Python `graphviz`
distribution.

A second native library is required and pip cannot install it either: **ICU**
(`libicu74` on Ubuntu 24.04, present in a default install). Wherever two or
more labels must be ordered as R orders them -- the default alphabet, R's
`seqdef` order, that every sequence constructor derives when no `alphabet=` is
given and two or more states are not numbers (`seqdef`, `seqformat`,
`seqdef_from_long`, `TSE_to_STS`, `HSPELL_to_STS` and
`events_to_state_sequences`); state labels in `seqstatl`,
`seqdiff`, `seqtab`, `seqtabstocc` and `rowmode`; group
levels in `seqimplic` and `seqsurv` (so also `seqsplot` and `seqsurvplot`);
tied events and person ids in `FCE_to_TSE`, `TSE_to_STS` (event mode) and
`toPersonPeriod`; the combined alphabet of `seqMD`; a text `sortv` in
`sequence_index_plot`, `seqplotMD`, `dissrf` and `seqrf`; text group levels in the grouped plots
`seqHtplot`, `seqedplot`, `seqplot_tentrop` and `seqplot` (and its per-type plots)
with `group`; event labels in `seqecreate` (and so the tie order of `seqefsub`); the rule text of `seqerules` with `sortv="rule"`; text group levels in `dissassoc`, `dissvar_grp`,
`disscenter`, `disstree`/`seqtree`, `seqrep_grp`, `dissmergegroups`,
`seqecmpgroup`, `group_p`, `seqsammeha` and the comparison functions
`dissCompare`, `seqCompare`, `seqLRT` and `seqBIC` -- traminex uses ICU's root
collation, because that is what R resolves under `C.UTF-8` and Python's code-point order
differs from it for case, accents and punctuation. If ICU cannot be loaded
these raise `RuntimeError` naming it; traminex deliberately
does **not** fall back to a portable ordering, because a silent fallback would
change which state is which in every per-state result without saying so. On a
minimal container install ICU explicitly:

```bash
sudo apt install libicu74     # Ubuntu 24.04
```

The order is ICU's root collation of the ICU that is loaded:
`ctypes.util.find_library("icui18n")` can find a conda environment's own ICU
(the `datasci` interpreter finds that environment's ICU 78) or the system's (ICU
74 on Ubuntu 24.04). R orders with the ICU it was built against (R 4.6.2 on
Ubuntu 24.04: ICU 74.2). ICU versions change over time, in R and in Python
environments alike, and a few characters order differently between major
versions -- `₨ ₧` under ICU 74 and in R 4.6.2 on Ubuntu 24.04, `₧ ₨` under ICU
78 -- so an automatic order can differ from R's when the two versions differ.
For an order that matches R's whatever the versions, pass `alphabet=`
explicitly: R and traminex keep the given order under any ICU version. This rule
fixes the order of the states; other text orders (group levels, event labels,
the combined alphabet of `seqMD`) follow the loaded ICU as well. The test run
prints the ICU version traminex loaded. To see which library is loaded: `python
-c "import ctypes.util; print(ctypes.util.find_library('icui18n'))"`.

## Quick usage

```python
import numpy as np
from traminex import seqdef, seqdist, seqstatd

X = np.array([
    ["A", "A", "B", "C"],
    ["A", "B", "B", "C"],
    ["C", "B", "A", "A"],
], dtype=object)

sts = seqdef(X, nr="*")
D = seqdist(sts, method="OM", subm="TRATE", norm="length")
st = seqstatd(sts)
```

## Tests

```bash
make test-python
make test-pybind11
```

## Quality checks

```bash
make release-check
make lint
make typecheck
make stub-check
```

## Documentation

- Authoritative project map: `docs/INDEX.md`
- Documentation hub: `docs/README.md`
- API usage examples: `docs/USAGE_new_apis.md`
- Distances guide: `docs/DISTANCES.md` (new methods: `docs/NEW_DISTANCE_METHODS.md`)
- Plotting guide: `docs/PLOTTING.md`
- API parity overview: `docs/API_PARITY.md`
- Public API manifest: `traminex/api_manifest.py`
- Release history: `CHANGELOG.md`
- Audits: every finding of the audit and repair runs, with its status and decision, is in
  [`docs/FINDINGS_REGISTER.csv`](docs/FINDINGS_REGISTER.csv); the owner's decisions are in
  [`docs/OWNER_DECISIONS.md`](docs/OWNER_DECISIONS.md). The runs' reports, the latest being the
  audit-only run of 2026-10-19, are kept in the family supervision repository
  (the parent folder `PROJECT_SEQANALYSIS_and_SPATSCAN`, folder `docs/history/run-reports/`).
  [`traminex_audit_report_20260612.md`](traminex_audit_report_20260612.md) is a historical
  closeout audit of the TraMineR 2.2.12 era
- The task prompts and verification records of those runs were removed from the root on
  2026-10-26; see `docs/INDEX.md` under "Prompt and verification records"

## Packaging and optional compiled kernels

The default wheel is compiler-free (`py3-none-any`) and installs on any
platform, but two R-faithful details depend on the platform. Several sums
reproduce R's `sum()`, which accumulates in an 80-bit `long double` on x86-64,
so they match R to the last bit only where NumPy's `long double` is that type
(as on x86-64 Linux). Without `alphabet=`, label order follows the root
collation of the ICU library that is loaded, as R's follows its own; ICU
versions change over time, so for an order that matches R's whatever the
versions, pass `alphabet=` explicitly (`docs/DEV_BUILD.md`).
The only supported compiled distance backend is the opt-in pybind11/OpenMP profile. See
`docs/DEV_BUILD.md` for the verified Ubuntu build, benchmark, and six isolated
install-profile instructions. Automatic dispatch may fall back to Python;
explicit `TRAMINEX_BINDING=pybind11` fails clearly when that backend is absent.

## License

GPL-3.0-or-later
