# Phase 2 shared input fixtures (migration plan Phase 2 / DUPLICATES A5-A8).
#
# Generates the deterministic dissimilarity matrix and weight vectors that all
# Phase 2 goldens (dissvar / dissassoc / disscenter / disstree) are computed
# from. The matrix is dumped to tests/r_goldens/ so the Python tests load the
# exact same numbers (no cross-language RNG coupling).
#
# Construction: 15 cases in 3 latent groups of 5 (groups rep(1:3, each=5)),
# 2-d coordinates with a deterministic (non-random) pattern, Euclidean
# distances. Weights: a non-integer vector (for weight.permutation="diss")
# and an integer vector (for "replicate").
#
# Requires R; see scripts/r_goldens/README.md for the run command.

# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
phase2_inputs <- function() {
  n <- 15
  i <- 1:n
  grp <- rep(1:3, each = 5)
  # Deterministic 2-d configuration: group centres far apart, within-group
  # spread irregular (no ties in the distance matrix).
  cx <- c(0, 10, 20)[grp] + (i %% 5) * 0.7 + (i %% 3) * 0.31
  cy <- c(0, 5, -5)[grp] + (i %% 4) * 1.13 + (i %% 7) * 0.19
  D <- as.matrix(dist(cbind(cx, cy)))
  dimnames(D) <- NULL
  wreal <- 0.5 + 0.1 * i           # 0.6, 0.7, ..., 2.0 (non-integer)
  wint <- rep(c(1, 2, 3), 5)       # integer weights for "replicate"
  list(D = D, grp = grp, wreal = wreal, wint = wint, n = n)
}

write_phase2_csv <- function(name, M, header, fmt = "%.15g") {
  out_dir <- file.path("tests", "r_goldens")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  M <- as.matrix(M)
  path <- file.path(out_dir, paste0(name, ".csv"))
  con <- file(path, "w")
  cat(header, file = con)
  for (r in seq_len(nrow(M))) {
    cat(paste(sprintf(fmt, M[r, ]), collapse = ","), "\n", file = con, sep = "")
  }
  close(con)
}

if (sys.nframe() == 0L || identical(environment(), globalenv())) {
  inp <- phase2_inputs()
  # The inputs are written with 17 significant digits, which round-trip R's
  # doubles exactly: the goldens are computed from R's in-memory values, and a
  # result decided by the last bit (disstree's choice between a level subset
  # and its complement, audit S4-04) is reproducible only from exact inputs.
  write_phase2_csv(
    "phase2_input_D", inp$D,
    "# Phase 2 shared input: 15x15 Euclidean distance matrix (3 latent groups of 5), see make_goldens_phase2_inputs.R\n",
    fmt = "%.17g"
  )
  write_phase2_csv(
    "phase2_input_weights_real", matrix(inp$wreal, nrow = 1),
    "# Phase 2 shared input: non-integer weights 0.5 + 0.1*i, i=1..15\n",
    fmt = "%.17g"
  )
  write_phase2_csv(
    "phase2_input_weights_int", matrix(inp$wint, nrow = 1),
    "# Phase 2 shared input: integer weights rep(c(1,2,3), 5)\n",
    fmt = "%.17g"
  )
  cat("PHASE2_INPUTS_DONE\n")
}
