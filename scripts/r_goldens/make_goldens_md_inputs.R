# Shared multidomain input fixtures for RP4 (dissdomassoc T14, seqdomassoc).
#
# Three deterministic 12x8 state-sequence channels (alphabet 1:3), their OM
# distance matrices (sm="CONSTANT", indel=1), and a non-integer weight vector.
# Dumped to tests/r_goldens/ as md_input_* so the Python tests load the exact
# same numbers (no cross-language seqdist coupling).
#
# Requires R 4.6.2 + TraMineR 2.2-14; see scripts/r_goldens/README.md.

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
source("scripts/r_goldens/make_goldens_phase2_inputs.R")  # write_phase2_csv()
suppressMessages(library(TraMineR))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

md_inputs <- function() {
  n <- 12
  L <- 8
  # Deterministic but non-degenerate channels (seeded sample; Python loads the
  # dumped matrices, so R's RNG need not be replicated).
  set.seed(20260618)
  mk <- function() matrix(sample(1:3, n * L, replace = TRUE), nrow = n, ncol = L)
  M1 <- mk()
  M2 <- mk()
  M3 <- mk()
  s1 <- suppressMessages(seqdef(M1, alphabet = 1:3))
  s2 <- suppressMessages(seqdef(M2, alphabet = 1:3))
  s3 <- suppressMessages(seqdef(M3, alphabet = 1:3))
  d1 <- suppressMessages(seqdist(s1, method = "OM", sm = "CONSTANT", indel = 1))
  d2 <- suppressMessages(seqdist(s2, method = "OM", sm = "CONSTANT", indel = 1))
  d3 <- suppressMessages(seqdist(s3, method = "OM", sm = "CONSTANT", indel = 1))
  w <- 0.5 + 0.1 * (1:n)
  list(
    n = n, L = L,
    s = list(s1, s2, s3),
    smat = list(M1 - 1L, M2 - 1L, M3 - 1L),  # 0-based code matrices
    d = list(d1, d2, d3),
    w = w
  )
}

if (sys.nframe() == 0L || identical(environment(), globalenv())) {
  inp <- md_inputs()
  for (k in 1:3) {
    write_phase2_csv(
      paste0("md_input_d", k), inp$d[[k]],
      paste0("# multidomain input: OM(sm=CONSTANT) distance matrix, channel ", k, "\n")
    )
    write_phase2_csv(
      paste0("md_input_s", k), matrix(as.integer(inp$smat[[k]]), nrow = inp$n),
      paste0("# multidomain input: 0-based state-code matrix, channel ", k, "\n")
    )
  }
  write_phase2_csv(
    "md_input_weights", matrix(inp$w, nrow = 1),
    "# multidomain input: non-integer weights 0.5 + 0.1*i, i=1..12\n"
  )
  cat("MD_INPUTS_DONE\n")
}
