# R oracle provenance

The parity oracle is the exact CRAN source release of:

- TraMineR 2.2-14, published 2026-09-15
- TraMineRextras 0.6.9, published 2026-06-20

Official package pages:

- https://cran.r-project.org/web/packages/TraMineR/index.html
- https://cran.r-project.org/web/packages/TraMineRextras/index.html

The extracted source trees are `SOURCES/TraMineR/` and
`SOURCES/TraMineRextras/`. The downloaded CRAN archives are **tracked** in this
repository as reproducibility inputs, and are part of the sdist allow-list. The
blanket `*.tar.gz` rule in `.gitignore` carries two explicit negations for them
(`!SOURCES/TraMineR_2.2-14.tar.gz`, `!SOURCES/TraMineRextras_0.6.9.tar.gz`), so
neither archive is ignored:

| Archive | CRAN MD5 | SHA-256 |
|---|---|---|
| `TraMineR_2.2-14.tar.gz` | `609ae6049321405695e5ed6aa7a7c6e7` | `8638b25d65d49ce78a78cc0b7ef3a1a4647d2f18d2d07b747b6da7e13aa71f43` |
| `TraMineRextras_0.6.9.tar.gz` | `834d4bfdc8f3bb0804d5ca225064e781` | `6bf0215bca2063669c5775b3656e0b3b05b5f820e343968042b588d490958b7b` |

The MD5 values come from CRAN's `src/contrib/PACKAGES` index. Every file in
both extracted trees was also verified against the package's internal `MD5`
manifest: 347 TraMineR entries and 89 TraMineRextras entries passed.

## Verified Ubuntu build

Verified on Ubuntu 24.04, Linux x86_64, R 4.6.2, GCC/G++ 13.3:

```bash
mkdir -p /tmp/traminex-r-lib
R CMD INSTALL --library=/tmp/traminex-r-lib SOURCES/TraMineR_2.2-14.tar.gz
R_LIBS=/tmp/traminex-r-lib R CMD INSTALL \
  --library=/tmp/traminex-r-lib SOURCES/TraMineRextras_0.6.9.tar.gz
R_LIBS=/tmp/traminex-r-lib Rscript -e \
  'stopifnot(packageVersion("TraMineR") == "2.2.14", packageVersion("TraMineRextras") == "0.6.9")'
```

CRAN does not provide a general Ubuntu binary on the package pages. The source
archives are therefore the portable reference and are compiled locally. Do not
introduce Windows binaries or Windows-specific oracle instructions.

## What the isolated library does *not* pin

The isolated library above contains exactly the two retained archives. It does
**not** pin the packages those two depend on that ship with the R installation
itself. `TraMineRextras::seqsurv` calls `survival::survfit`, and `survival` is a
*recommended* R package: with `R_LIBS` pointing at the isolated library it still
resolves from the system library (verified — `/usr/lib/R/site-library/survival`).

Consequently a third-party version must never be written into a fixture under
`tests/r_goldens/` as data. Doing so makes that fixture reproducible only on a
machine whose copy of the package matches exactly, which breaks the documented
two-pass provenance proof for everyone else (cycle-4b finding F-002). Instead:

- assert a **documented compatible range** inside the generator with
  `stopifnot`, so a mismatch stops generation rather than silently changing a
  recorded value; and
- let `scripts/check_r_golden_provenance.py` record the resolved version in its
  evidence JSON, where a change is informative but is not a fixture byte-diff.

The maintained floors live in `UNPINNED_DEPENDENCY_FLOORS` in that script.
Currently: `survival >= 3.8.9`.

## Refresh procedure

1. Download the new source archives and record CRAN MD5 plus SHA-256 values.
2. Extract them without modifying their contents and replace the two vendored
   trees.
3. Refresh `DATA/`, `DEMOS/`, and the CRAN reference PDFs in `DOCUMENTATION/`.
4. Compile both packages into an isolated R library on Ubuntu.
5. Run every script under `scripts/r_goldens/` with `R_LIBS` pointing at that
   library.
6. Run the complete Python suite and reconcile every behavioral delta before
   claiming the new release as the oracle.
7. Regenerate the API audit and parity tables.

Historical audit documents may continue to name the oracle that was current
when they were written. Active docs and tests must name the versions above.
