# Developer build notes

Ubuntu 24.04 with Python 3.11 or 3.12 is the supported development matrix. The
package works without a compiler; NumPy/Numba provide the portable fallback.
The only supported compiled backend is pybind11.

## Environment

Use the family-standard interpreter and the system GCC/CMake/Ninja toolchain:

```bash
PY=~/miniforge3/envs/datasci/bin/python
$PY --version
gcc --version
g++ --version
cmake --version
ninja --version
```

Install the editable package and build requirements in that environment:

```bash
$PY -m pip install -e '.[plot,perf,extras,tree,dev]'
```

That command changes the shared `datasci` environment, which as maintained has
neither traminex nor `tokenize_rt` (which `black[jupyter]`, and so
`make docs-check`, needs). To leave `datasci` unchanged, use a private virtual
environment layered on it, as the maintained gates of this checkout do:

```bash
B=$HOME/.local/share/traminex-repair
~/miniforge3/envs/datasci/bin/python -m venv --system-site-packages "$B/venv"
"$B/venv/bin/python" -m pip install --no-deps --no-build-isolation -e .
"$B/venv/bin/python" -m pip install tokenize_rt==6.2.0
PY=$B/venv/bin/python
```

`--system-site-packages` reuses `datasci`'s NumPy, pandas, SciPy and tools, and
`--no-deps --no-build-isolation` installs only the editable traminex; pass
`PYTHON=$PY` to the Make gates.

### The private environment's helper scripts

The maintained gates of this checkout run with `source $B/env.sh`, and long runs are
logged with `$B/run_logged.sh`. Their full text follows, with the home folder written `$HOME`; the files themselves spell
the folder out, and their SHA-256 values (`ba1def00…` and `08e34582…`, recorded for repair
run `repair-20261012T025723Z`) are those of that spelled-out text. Either can be recreated:

```bash
# $B/env.sh
# Environment for every traminex repair command: `source $HOME/.local/share/traminex-repair/env.sh`
export B=$HOME/.local/share/traminex-repair
export PY=$B/venv/bin/python
export TRAMINEX_PYTHON=$PY
export PYTHONPYCACHEPREFIX=$B/pycache
export GVPATH=$HOME/miniforge3/envs/datasci
export MPLBACKEND=Agg
export R_LIBS=$B/rlib
```

```bash
# $B/run_logged.sh
#!/bin/bash
# Run one command with env.sh sourced, keep its full output, and append a JSON line
# (start, end, exit, cwd, output, argv) to a log that COMMANDS.json is built from.
# Usage: run_logged.sh <log.jsonl> <output-file> -- <command> [args...]
# Prints the last 5 lines of output; exits with the command's status.
set -u
if [ $# -lt 4 ] || [ "$3" != "--" ]; then
  echo "usage: run_logged.sh LOG.jsonl OUTPUT -- COMMAND [ARGS...]" >&2; exit 2
fi
log=$1; out=$2; shift 3
source $HOME/.local/share/traminex-repair/env.sh
mkdir -p "$(dirname "$log")" "$(dirname "$out")"
start=$(date -u +%Y-%m-%dT%H:%M:%SZ); cwd=$(pwd -P)
"$@" > "$out" 2>&1; rc=$?
end=$(date -u +%Y-%m-%dT%H:%M:%SZ)
/usr/bin/python3 - "$log" "$start" "$end" "$rc" "$cwd" "$out" "$@" <<'PY'
import json, sys
log, start, end, rc, cwd, out, *argv = sys.argv[1:]
with open(log, "a", encoding="utf-8") as f:
    f.write(json.dumps({"start": start, "end": end, "exit": int(rc), "cwd": cwd, "output": out, "argv": argv}) + "\n")
PY
tail -n 5 "$out"
exit $rc
```

The environment's R oracle library is `$B/rlib`, built as `scripts/r_goldens/README.md`
describes, with `RLIB=$B/rlib` in place of the temporary folder.

The editable install also puts the repository root on `sys.path`. Its `.pth`
file in `site-packages` (`_editable_skbc_traminex.pth`) holds, besides the import
hook, the parent directory of each package, as scikit-build-core writes it; with
this flat layout (`traminex/` at the root) that parent is the root itself. So
`import benchmarks` works from any directory and loads this checkout's
`benchmarks/`, and the root's `scripts/` and `tests/` folders can mix with, or
lose to, an installed top-level package of the same name. Run the gates and the
scripts from the repository root. A `src/` layout (the package under
`src/traminex/`) would end this, since the parent directory would then hold only
the package; that move is a separate decision, not taken.

### Native requirements pip cannot install

Two system libraries are required and appear in no dependency declaration,
because pip has no way to supply either. Both ship in a default Ubuntu 24.04
install; a minimal container may not have them.

| Requirement | Ubuntu 24.04 package | Used by | Behaviour when absent |
|---|---|---|---|
| ICU | `libicu74` | R-compatible label collation wherever two or more labels are ordered: the default alphabet, R's `seqdef` order, that every sequence constructor derives when no `alphabet=` is given and two or more states are not numbers (`seqdef`, `seqformat`, `seqdef_from_long`, `TSE_to_STS`, `HSPELL_to_STS` and `events_to_state_sequences`; owner decision of 2026-10-19); state labels in `seqstatl`, `seqdiff`, `seqtab`, `seqtabstocc` and `rowmode`; group levels in `seqimplic` and `seqsurv` (so also `seqsplot` and `seqsurvplot`); tied events and person ids in `FCE_to_TSE`, `TSE_to_STS` (event mode) and `toPersonPeriod`; the combined alphabet of `seqMD`/`seqdistmc`; two or more text group levels in `dissassoc`, `dissvar_grp`, `disscenter`, `disstree`/`seqtree` (text predictors), `seqrep_grp` and `dissmergegroups`; a text `sortv` in `sequence_index_plot` (also `seqplot(type="i")`), `seqplotMD`, `dissrf` and `seqrf`; text group levels in the grouped plots `seqHtplot`, `seqedplot`, `seqplot_tentrop` and `seqplot` (and its per-type plots) with `group`; event labels in `seqecreate` (so the tie order of `seqefsub`); the rule text of `seqerules` with `sortv="rule"`; text group levels in `seqecmpgroup`, `group_p`, `seqsammeha`, `dissCompare`, `seqCompare`, `seqLRT` and `seqBIC` | `RuntimeError: R-compatible character collation requires ICU (…)` at call time |
| GraphViz `dot` | `graphviz` | `disstreedisplay`, `seqtreedisplay` | `RuntimeError: GraphViz 'dot' was not found` at render time |

```bash
sudo apt install libicu74 graphviz
```

`traminex/_categorical_order.py` loads ICU through its stable C API with
`ctypes.util.find_library("icui18n")`; R itself is never invoked. The `dot`
locator honours `GVPATH` first and then `PATH`
(`traminex/disstree_viz.py::_locate_dot`), so an off-`PATH` `dot` can be
selected without changing `PATH`.

The label order is ICU's root collation of the ICU that is loaded:
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

**Both failures are deliberately loud.** In particular, traminex has no portable
fallback ordering for a missing ICU: no portable key reproduces R's collation for
case, punctuation and accents, so a fallback would silently permute the state
alphabet — and therefore every per-state result — without the user being told.
The owner confirmed this policy on 2026-10-12 for every function that needs R's
collation; the case-folded key traminex once used for event and group labels was
replaced by ICU collation on 2026-10-13 (TX-AUD-20261011-015, -017). Changing this
to a fallback is a maintainer decision that must be recorded as a new divergence,
not made incidentally.

### Testing a second checkout (a git worktree)

The editable install registers a scikit-build-core import hook
(`_editable_skbc_traminex` in `site-packages`) in `sys.meta_path`, ahead of the
ordinary path search. It maps every `traminex` module to the checkout that was
installed, whatever the current directory, `PYTHONPATH` or
`sys.path.insert(0, ...)` say. So in a second checkout of the repository, such
as a `git worktree`, a plain `python` silently imports the *installed*
checkout's code.

Child interpreters make this worse, because a fix applied inside one process
does not reach them. Several tests start `python -c` children that import
traminex -- `tests/test_accel_thread_determinism.py` and
`tests/test_numba_backend_agreement.py` among them. Run from a second checkout,
those children test the installed checkout's code instead, and under
`make test-pybind11` they fail at setup because the installed checkout has no
`_pb_kernels` link (recorded by the 2026-09-19 audit).
`scripts/_checkout_imports.prefer_checkout_imports()` reorders the finders
only in the process that calls it, and `python -S` is not passed on to
children either.

The mitigation that reaches every child is a `sitecustomize.py`, kept in a
directory **outside** the repository, that removes traminex's hook when each
interpreter starts:

```python
# sitecustomize.py
import sys

sys.meta_path[:] = [
    finder
    for finder in sys.meta_path
    if not (
        type(finder).__name__.startswith("ScikitBuild")
        and "traminex" in getattr(finder, "pkgs", ())
    )
]
```

Put its directory, then the second checkout's root, on `PYTHONPATH` for every
command. Child processes inherit the variable, so they run the hook too:

```bash
export PYTHONPATH=/path/to/shim-dir:/path/to/second-checkout
$PY -c "import subprocess, sys, traminex; print(traminex.__file__); \
subprocess.run([sys.executable, '-c', 'import traminex; print(traminex.__file__)'])"
```

Both printed paths must lie in the second checkout; treat a run as evidence
for that checkout only after confirming it. Only traminex's hook is removed,
so the other editable packages in the environment keep working. A
`sitecustomize.py` placed first on `PYTHONPATH` takes the place of any other
`sitecustomize` module in the environment. Leave the shim out when you work
in the installed checkout itself.

## Local pybind11 build

```bash
PY=~/miniforge3/envs/datasci/bin/python
PYBIND11_DIR=$($PY -m pybind11 --cmakedir)
cmake -S . -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DTRAMINEX_BUILD_PYBIND11=ON \
  -DPython_EXECUTABLE="$PY" \
  -Dpybind11_DIR="$PYBIND11_DIR"
cmake --build build --parallel
```

The resulting `_pb_kernels*.so` must be importable as `traminex._pb_kernels`.
For a local source-tree test, copy or symlink the module into `traminex/`; do not
commit the generated shared object.

**Remove it when you are done.** Backend selection is *import-based*:
`TRAMINEX_BINDING=auto` (the default) resolves to `pybind11` the moment
`traminex._pb_kernels` imports, so a forgotten copy or a stale symlink silently
switches every later run -- including test runs you believe are pure Python --
onto the compiled kernels, and a symlink pointing at a deleted build directory
fails at import instead. Nothing in the package warns about this. `.gitignore`
keeps the object untracked and the sdist allow-list excludes it, so the hazard
never reaches a distribution; `make test-pybind11` links and unlinks its own
copy and refuses to overwrite one that is already there.

The hazard does reach anything that imports this working tree, though: an
environment with traminex installed editable from it (such as the private one
under "Environment") imports the link, and the downstream
`pyweightedcluster` verifies the package tree on every import and refuses *any*
symlink in it. So for as long as `make test-pybind11` runs, that project cannot
import traminex. The script's `EXIT` trap removes the link on a normal exit, on a
failure, and on `SIGINT`, `SIGTERM`, `SIGHUP` and `SIGQUIT`; only an
uncatchable `SIGKILL` (or a power loss) leaves it behind, together with the build
directory it points into, so the link stays live until `/tmp` is cleared. The
runner writes its process ID into its build directory (`runner.pid`), and the
next run, before it builds (`scripts/_native_link.sh`), removes a link into one of
its own `/tmp/traminex-r1-pybind11.*` build directories when the target is gone,
or when the recorded process no longer exists and no process has the module
loaded; it then removes that build directory too. A link of a run in progress, a
build directory without `runner.pid`, a regular file or a link pointing elsewhere
is left alone and refused. After any interrupted
run, confirm that `find traminex -type l` prints nothing. A copied module would
not be safer: the consumer does not refuse a regular file, and `auto` dispatch
would then use a stale compiled kernel silently.

The canonical gate performs that temporary linking safely, verifies pybind11
and OpenMP, runs the complete suite, and cleans up:

```bash
make test-pybind11
```

It builds and tests with the Makefile's `PYTHON` (`make test-pybind11 PYTHON=...`),
which it passes to `scripts/run_pybind11_tests.sh` as `TRAMINEX_PYTHON`; run
directly, the script uses `TRAMINEX_PYTHON` or the `python3` on `PATH`.

## Artifacts

The raw CMake option defaults to OFF. A normal build therefore produces a
compiler-free `py3-none-any` wheel plus the source distribution:

```bash
$PY -m build --sdist --wheel
```

The source distribution carries `tests/` and `scripts/` for reference, but the
test suite and the release gates must be run from a git checkout, not from an
unpacked sdist. Some tests list tracked files with `git ls-files`, and the
semantic-ledger and documentation-consistency checks read repository-only
evidence that the sdist allow-list deliberately leaves out, such as
`reports/audit/traminex_semantic_parity.json` and the design note under
`docs/superpowers/specs/` that rows of `docs/API_PARITY_DECISIONS.csv` cite.

The sdist is an explicit allow-list (`sdist.include` in `pyproject.toml`).
`make sdist-manifest` writes it from `git ls-files` and the source policy in
`scripts/sync_sdist_manifest.py`; run it after staging a new, removed or renamed
file (`scripts/sync_sdist_manifest.py --check` writes nothing and names each
missing or extra entry). scikit-build-core 0.12 and later, which an isolated `python -m build` installs
unless a version is pinned, read it in manual mode and read no `.gitignore`
file. The declared floor, 0.11.0, does not know that mode: it reads the list
with gitignore rules, so every root entry is anchored with `/`, and it reads
each nested `.gitignore` file while walking the tree, ignored directories
included. With 0.11.x, build the sdist from a clean clone or a `git archive`
copy: in a checkout that holds a copied tool cache with a malformed
`.gitignore` (a bare `!` line, for example) the build stops with
`GitIgnorePatternError`.

The opt-in compiled wheel passes `-DTRAMINEX_BUILD_PYBIND11=ON` **as a CMake
argument through scikit-build config settings**, together with
`wheel.platlib=true` and an empty `wheel.py-api`, which produces a platlib
wheel with the interpreter's native ABI tag:

```bash
PYBIND11_DIR="$($PY -m pybind11 --cmakedir)"
$PY -m build --wheel \
  -Ccmake.args="-DTRAMINEX_BUILD_PYBIND11=ON;-Dpybind11_DIR=$PYBIND11_DIR" \
  -Cwheel.platlib=true -Cwheel.py-api=
```

Setting `TRAMINEX_BUILD_PYBIND11=ON` only as an environment variable does
**nothing** — the build silently produces the pure `py3-none-any` wheel. The
maintained implementation and installed-wheel smoke are in
`.github/workflows/wheels.yml`.

Default artifact construction must remain compiler-free. A deliberately
invalid `CC`/`CXX` is therefore a valid audit probe for the default wheel, but
not for the explicit compiled profile.

## Isolated install profiles

The maintained verifier reuses exactly one default sdist, one default pure
wheel, and one explicit pybind11/OpenMP wheel across six fresh environments:
core, performance, plotting, tree rendering, compiled identity, and automatic
fallback. Imports must resolve from the temporary environment, never this
checkout.

```bash
make install-profiles \
  SDIST=/absolute/path/to/traminex-0.2.2.tar.gz \
  PURE_WHEEL=/absolute/path/to/traminex-0.2.2-py3-none-any.whl \
  COMPILED_WHEEL=/absolute/path/to/traminex-0.2.2-cp312-cp312-linux_x86_64.whl \
  PYTHON311=/absolute/path/to/python3.11 \
  UV=/absolute/path/to/uv \
  DOT=/absolute/path/to/dot \
  INSTALL_PROFILE_EVIDENCE=/absolute/path/to/install-profiles.json
```

Without `INSTALL_PROFILE_EVIDENCE`, the report is written to
`build/install-profiles.json`, under the git-ignored `build/` directory.

The verifier hashes and inspects all artifacts before installation, records
every command and `pip check`/`pip freeze`, and deletes only a temporary root
carrying its unpredictable ownership marker. The Python 3.11 core/performance
profiles use the declared dependency floors; Python 3.12 exercises the current
plot, tree, compiled, and fallback profiles.

## OpenMP

`cpp/src/accel_kernels.h` contains independently implemented, R-source-derived
OpenMP-parallel OM, LCS, and unit-HAM kernels. LCS and OM use direct one-row
dynamic programming. Exact duplicate rows are grouped before evaluation, and
the remaining upper-triangle pairs are divided into contiguous static ranges;
each worker owns its workspace. CMake links `OpenMP::OpenMP_CXX` when the
Ubuntu toolchain provides it. Check the active backend and thread visibility
with:

```bash
TRAMINEX_BINDING=pybind11 $PY - <<'PY'
from traminex import _pb_kernels
from traminex import distances

assert distances._impl_name() == "pybind11", distances._impl_name()
assert hasattr(_pb_kernels, "accel_openmp_enabled")
assert _pb_kernels.accel_openmp_enabled() is True
print("verified backend:", distances._impl_name())
print("OpenMP threads:", _pb_kernels.accel_max_threads())
PY
```

A Python fallback or skipped compiled differential module is a failed pybind11
gate even when pytest exits zero. `scripts/run_pybind11_tests.sh` enforces the
second half: a test that skips because the compiled kernels are unbuilt or not
in use (the reasons listed in `scripts/_native_skip_gate.py`) makes the run
fail and is named in the summary; other skips, such as a missing Graphviz or
Rscript, do not.

The compiled kernels must remain byte-identical to the Python/Numba reference
on their supported raw-distance contracts. Run:

```bash
TRAMINEX_BINDING=pybind11 $PY -m pytest -q \
  tests/test_accel_byte_identity.py \
  tests/test_compiled_dispatch_rlcp_dhd.py
```

## Python fallback and benchmarks

Force the portable backend with `TRAMINEX_BINDING=python`:

```bash
TRAMINEX_BINDING=python $PY -m pytest -q
$PY -m benchmarks.run_benchmarks --method OM --backend python \
  --backend numba --N 8 --L 6 --K 3 --warmups 1 --repeats 2
```

## R oracle build on Ubuntu

The parity oracle is TraMineR 2.2-14 plus TraMineRextras 0.6.9. Build the exact
vendored CRAN archives into an isolated R library:

```bash
mkdir -p /tmp/traminex-r-lib
R CMD INSTALL --library=/tmp/traminex-r-lib SOURCES/TraMineR_2.2-14.tar.gz
R_LIBS=/tmp/traminex-r-lib R CMD INSTALL \
  --library=/tmp/traminex-r-lib SOURCES/TraMineRextras_0.6.9.tar.gz
R_LIBS=/tmp/traminex-r-lib Rscript -e \
  'stopifnot(packageVersion("TraMineR") == "2.2.14", packageVersion("TraMineRextras") == "0.6.9")'
```

See `SOURCES/README.md` and `scripts/r_goldens/README.md` for provenance and
golden regeneration.

## Release checks

```bash
make release-check
make lint
make typecheck
make stub-check
```

These Make targets are the gates that have run. The GitHub Actions workflows
under `.github/workflows/` (`ci.yml`, `render-docs.yml`, `wheels.yml`) are
configured, but the repository has no git remote, so they have never run; a
statement below about a workflow job describes its configuration.

Three release steps are run by hand before a release, because no workflow job
runs them (TX-AUD-20261010-091): the six isolated install profiles
(`make install-profiles`, under "Isolated install profiles" above), the
two-pass provenance proof of the R goldens in a disposable clone
(`scripts/check_r_golden_provenance.py --execute --passes 2`, as
`scripts/r_goldens/README.md` gives it), and rendering the documentation
(`make docs-render`, which needs Quarto).

The gates need the editable install with the `dev` extra shown under
"Environment" (`black[jupyter]` brings `tokenize-rt`, and `make docs-check`
imports the installed traminex). `make release-check` first runs
`scripts/check_dev_environment.py`, which names every missing module and that
install command instead of stopping with a bare `ModuleNotFoundError`.

`make typecheck` runs mypy for the interpreter it runs under: the configuration
fixes no target version, because NumPy 2.5 and later require Python 3.12 and
ship stubs in 3.12 syntax. Python 3.11 syntax stays checked by black and ruff
(`target-version = py311`). The CI workflow is configured to run
`make typecheck` on Python 3.11 and on 3.12 after upgrading NumPy, pandas and
SciPy to the newest versions the dependency floors permit (job
`typecheck-newest`). The 3.11 row is the one that rejects a
standard-library call new in 3.12, such as `itertools.batched`; on 3.11, pip
resolves the newest NumPy that still supports 3.11. Locally, `make typecheck`
under a 3.12 interpreter does not catch such a call; run it under 3.11 too.
