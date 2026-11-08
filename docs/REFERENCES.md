# Optional pinned-R reference fixtures

The consumer is `tests/test_parity_refs.py`. Its default absence is a skip;
input alone cannot verify parity. A complete input/expected pair executes a
comparison. Missing additional expectations are reported as partial coverage.

Use TraMineR 2.2-14 (R reports 2.2.14). `input_wide.csv` has time columns and
sequence rows, with nonmissing character states. Both sides use the default
alphabet, which is R's `seqdef` order in traminex too (owner decision of
2026-10-19), so no alphabet is passed. Frequencies are
transposed explicitly from R state-by-time to Python time-by-state; no marker
columns are included when no missing or structural void cells occur. Missing
or ragged fixtures need explicit marker/name alignment rather than a blind
transpose or positional padding.

Run from the repository root with a pinned isolated R library. The expected
files are CSV matrices without headers or row names; R generates every value.

```r
library(TraMineR)
stopifnot(packageVersion("TraMineR") == "2.2.14")
base <- "tests/references"
df <- read.csv(file.path(base, "input_wide.csv"), check.names=FALSE,
               stringsAsFactors=FALSE)
seq <- seqdef(df)
write_reference <- function(value, name) {
  write.table(value, file.path(base, name), sep=",", row.names=FALSE,
              col.names=FALSE)
}
write_reference(t(seqstatd(seq)$Frequencies), "expected_seqstatd_freq.csv")
sm <- seqsubm(seq, method="CONSTANT", cval=2)
write_reference(as.matrix(seqdist(seq[1:2,], method="OM", indel=1, sm=sm)),
                "expected_dist_om.csv")
write_reference(seqsubm(seq, method="TRATE"), "expected_trate_subm.csv")
write_reference(as.matrix(seqdist(seq[1:2,], method="DHD")), "expected_dhd.csv")
```

The TRATE comparison uses the documented low-level symmetric cost extension
with smoothing 1e-6 and tolerance 1e-6; the R fixture above is unsmoothed. The
DHD comparison uses public `seqdist` with R's default transition costs and the
maintained tolerance 2e-3. The low-level `variant="mean"` frequency extension
is a different distance and must not be compared with this R default. These
four comparisons do not certify every public option or backend.

Run the consumer in the isolated development interpreter:

```bash
TRAMINEX_BINDING=python TRAMINEX_NUMBA=off "$REPAIR_PY" -m pytest -q tests/test_parity_refs.py
```
