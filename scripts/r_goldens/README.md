# R-golden operations

The files in this directory generate the static expected values under
`tests/r_goldens/`. The supported oracle is Ubuntu 24.04 with R 4.6.2,
TraMineR 2.2-14, and TraMineRextras 0.6.9. Archive hashes, extracted-tree
integrity, and the oracle refresh policy live in
[`SOURCES/README.md`](../../SOURCES/README.md).

## Isolated oracle library

Run from the repository root. Install the retained CRAN source archives into an
isolated library and verify every resolved version before generating evidence:

```bash
R=/usr/bin/Rscript
RLIB=/tmp/traminex-r-lib

mkdir -p "$RLIB"
R CMD INSTALL --library="$RLIB" SOURCES/TraMineR_2.2-14.tar.gz
R_LIBS="$RLIB" R CMD INSTALL \
  --library="$RLIB" SOURCES/TraMineRextras_0.6.9.tar.gz
R_LIBS="$RLIB" "$R" -e \
  'cat(R.version.string, as.character(packageVersion("TraMineR")), as.character(packageVersion("TraMineRextras")), sep="\n"); stopifnot(getRversion() == "4.6.2", packageVersion("TraMineR") == "2.2.14", packageVersion("TraMineRextras") == "0.6.9")'
```

Every new or changed generator must also assert the package version inside the
script, as code that stops on another version: `stopifnot()` of
`packageVersion("TraMineR") == "2.2.14"`, or an `if` on the inequality whose body calls
`stop()` (a comment, or a printed or stored comparison, is not an assertion;
`tests/test_r_golden_provenance.py` checks every generator), and call the public package namespace (`TraMineR::name` or
`TraMineRextras::name`) where practical. Private calls require an explicit
evidence reason.

This library pins **only** the two retained archives. Packages that ship with
the R installation — `survival`, used by `TraMineRextras::seqsurv` — still
resolve from the system library even under `R_LIBS`. Never write such a version
into a fixture as data: assert a documented compatible range in the generator
with `stopifnot`, and let the provenance checker record the resolved version in
its evidence JSON. See
[`SOURCES/README.md`](../../SOURCES/README.md) ("What the isolated library does
*not* pin") for the rationale and the maintained floors.

## Targeted regeneration

Regenerate only the generator owned by the contract being changed:

```bash
R_LIBS="$RLIB" "$R" scripts/r_goldens/make_goldens_stochastic_phase4.R
git diff -- tests/r_goldens
PY=~/miniforge3/envs/datasci/bin/python
$PY -m pytest -q tests/test_r_goldens_stochastic_phase4.py
```

Inspect every fixture diff. A generator completing successfully does not prove
that a changed expected value is correct.

## Complete sorted regeneration

The inventory has three explicit roles:

- **canonical fixture writers** create or refresh reviewed files in
  `tests/r_goldens/`;
- `make_goldens_corestats.R` is a **sourced helper** that is also executed as a
  standalone health check; and
- `make_goldens_distances.R`, `make_goldens_norm.R`, and
  `make_goldens_norm_costs.R` are **stdout-only diagnostics** retained for
  interactive oracle inspection.

Do not infer ownership from filenames or execute an unclassified glob. The
provenance checker contains the maintained role declarations, executes every
script in byte-sorted order, detects which fixtures each canonical writer
directly rewrites, and rejects missing or undeclared multiple ownership.

First validate the static inventory:

```bash
git status --short
PY=~/miniforge3/envs/datasci/bin/python
$PY scripts/check_r_golden_provenance.py
```

Besides the pinned counts, static mode requires every fixture to be named by
at least one canonical writer: a string literal in the writer (or its last
path component) must be the fixture's relative path, file name or stem, or,
for names built with `paste0`, a name prefix of its stem: a literal that
contains `_` or `-` after which the name goes on at a boundary (the literal ends
in `_` or `-`, or the next character is `_`, `-`, `.` or a digit). A word such
as `"diss"` from `what = "diss"` is not a name prefix. Every fixture name a
writer spells in full must also exist, so replacing such a fixture by an
orphan of the same count fails. A fixture whose name a writer builds is named
by its prefix only, so an orphan with that prefix passes the static check:
`md_input_d1.csv` (from `paste0("md_input_d", k)`) renamed to
`md_input_d9.csv` keeps the count and stays attributed. The static check thus
attributes fixtures by name only (it is also asserted by
`tests/test_r_golden_provenance.py` and
`tests/test_r_golden_provenance_orphans.py`); only the execute mode below
proves which writer actually writes each fixture. It rejects that orphan: the
writer writes `md_input_d1.csv` again, and the fixture count it checks after
every script no longer matches.

Then run the complete two-pass proof in a disposable repository copy. The
checker verifies the exact oracle versions (R 4.6.2 and the two packages), invokes all declared script roles,
requires every fixture to have a direct canonical writer, compares the full
fixture tree with its reviewed starting bytes after each pass, and requires
the ownership map to be identical between passes:

```bash
COPY="$(mktemp -d)/traminex"
git clone --quiet . "$COPY"
R_LIBS="$RLIB" "$PY" scripts/check_r_golden_provenance.py \
  --repo "$COPY" --execute --passes 2 --rscript "$R" \
  --evidence /tmp/traminex-r-golden-provenance.json
```

The generators run with the copy as their working directory, so they write
only the copy's fixtures. `--execute` refuses the checkout that holds the
script (also through a symlinked path) and any git work tree with uncommitted
changes under `tests/r_goldens`, before any R process starts; a clone of
committed work, as above, is accepted. `--allow-in-place` overrides the
refusal.

The evidence JSON records the complete inventory, role for every script,
direct fixture ownership, oracle versions, per-script output hashes, and the
fixture-tree hash after each pass. The oracle versions include `dependencies`:
the version R loads for every package of the recursive Depends/Imports closure
of the packages the generators load or call with `::` (TraMineR,
TraMineRextras, WeightedCluster, vegan, survival, colorspace, RColorBrewer;
found by a static scan of the scripts), computed by R at run time. A package of
that closure that R cannot find stops the proof. Only the two archives and the
`survival` floor are asserted; the other versions are recorded, not pinned.
An unexplained byte change or second-run
ownership difference is a failed evidence design, not an acceptable refresh.

## Evidence rules

- Serialize numeric values with enough precision for the focused comparison;
  preserve integer/double/logical/character distinctions where observable.
- Preserve string bytes, ordering, dimensions, row/column names, factor levels,
  missing values, warnings, messages, and errors needed by the contract.
- Avoid timestamps, absolute paths, locale accidents, and unclassified random
  output.
- Every script's first code line pins the locale and R's message language:
  `stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")`.
  Without it, a run from a C-locale shell sorted text by bytes, read non-ASCII
  literals as bytes (even with `icuSetCollate`) and quoted messages with `'`
  instead of `‘ ’`, changing 24 tracked fixtures; and the locale alone leaves
  R's messages in the shell's language (a `LANGUAGE=de` shell wrote 45 fixtures
  in German). Without ICU, `icuSetCollate` only warns, so the line also
  requires it. `--execute` also runs every
  script with `LC_ALL=C.UTF-8` and `LANGUAGE=en` and records them in its
  evidence; `tests/test_r_golden_locale_pin.py` requires the line in every
  script. A script does not switch the locale again: a later
  `Sys.setlocale` to a UTF-8 collation resets R's ICU collator (it undid the
  pin in three scripts when R started in a C shell). Byte collation,
  `Sys.setlocale("LC_COLLATE", "C")`, the same from any shell, is the one
  allowed later switch.
- Deterministic contracts compare public R values directly.
- Stochastic contracts either reproduce the R stream with an asymmetric
  fixture or record deterministic sufficient components plus R global-RNG
  behavior. Python tests then prove its actual seeded reproducibility, schema,
  invariants, edge cases, and the precise quantity that is not cross-language
  reproducible.
- Never use Python output to create an expected R value.
- A record a generator writes itself rather than taking from an R call (an R
  statement re-implemented because the function only draws, a hand-typed
  message, traminex-side metadata, a mapped error code) is listed in
  `scripts/r_goldens/AUTHORED_RECORDS.csv` with the reason;
  `tests/test_r_golden_authored_register.py` keeps the list true to the
  generators and fixtures it names.

Each generator owns one or more files in `tests/r_goldens/`; a focused Python
test must load those files and prove the documented public contract. The
semantic decision row cites the defining R source, generator when applicable,
and focused test.
