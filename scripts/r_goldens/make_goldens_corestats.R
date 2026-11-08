# Shared input fixtures for the RP3 core-stats goldens (T3-T11).
#
# Reuses the Phase-1 / seqST A/B/C fixtures (see make_goldens_seqst.R and
# make_goldens_seqdss_seqdur.R) plus a per-sequence weights vector, so every
# core-stats golden (seqstatd/seqistatd/seqmeant/seqmodst/seqtab/seqtabstocc/
# seqtrate/seqstatf/seqstatl/seqasnum/seqtransn/seqindic) is computed from the
# same numbers across the weighted, with.missing and ragged (void) paths.
#
#   A  : void-free, no missing (alphabet A,B,C,D)
#   B  : interior missing "*" (right=NA -> trailing "*" stay nr, void-free)
#   C  : ragged, trailing "*" -> void (default right="DEL"), interior "*" -> nr
#   w  : per-sequence weights c(1.5, 2.0, 0.5, 3.0)
#
# Run from the repository root with the oracle library (scripts/r_goldens/README.md):
#   Rscript scripts/r_goldens/make_goldens_corestats.R   (sourced by the per-symbol generators)
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# %.15g numeric matrix writer (matches r_goldens_utils.load_golden).
write_num_matrix <- function(name, m, header) {
  m <- as.matrix(m)
  con <- file(file.path(out_dir, paste0(name, ".csv")), open = "wt")
  cat(header, file = con)
  txt <- apply(m, 1, function(r) paste(formatC(r, format = "g", digits = 15), collapse = ","))
  writeLines(txt, con)
  close(con)
}

corestats_weights <- function() c(1.5, 2.0, 0.5, 3.0)

corestats_fixtures <- function() {
  w <- corestats_weights()

  A <- rbind(
    c("A", "A", "B", "C", "C", "C"),
    c("A", "B", "B", "B", "C", "A"),
    c("D", "D", "D", "D", "D", "D"),
    c("A", "B", "A", "B", "A", "B")
  )
  sdA <- seqdef(A, alphabet = c("A", "B", "C", "D"))
  sdAw <- seqdef(A, alphabet = c("A", "B", "C", "D"), weights = w)

  B <- rbind(
    c("A", "A", "*", "B", "B", "C"),
    c("*", "A", "A", "A", "*", "*"),
    c("A", "*", "A", "B", "B", "B"),
    c("C", "C", "C", "C", "C", "C")
  )
  sdB <- seqdef(B, alphabet = c("A", "B", "C"), missing = "*", right = NA)
  sdBw <- seqdef(B, alphabet = c("A", "B", "C"), missing = "*", right = NA, weights = w)

  C <- rbind(
    c("A", "A", "B", "C", "C", "C"),
    c("A", "B", "B", "*", "*", "*"),
    c("A", "*", "A", "B", "*", "*"),
    c("C", "C", "*", "*", "*", "*")
  )
  sdC <- seqdef(C, alphabet = c("A", "B", "C"), missing = "*")
  sdCw <- seqdef(C, alphabet = c("A", "B", "C"), missing = "*", weights = w)

  list(A = sdA, Aw = sdAw, B = sdB, Bw = sdBw, C = sdC, Cw = sdCw, w = w)
}

if (sys.nframe() == 0L || identical(environment(), globalenv())) {
  invisible(corestats_fixtures())
  cat("CORESTATS_FIXTURES_OK\n")
}
