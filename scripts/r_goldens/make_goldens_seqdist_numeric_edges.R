# Golden generation: numeric edges of seqdist (2026-09-22 audit S1a-05,
# S1b-08, S1b-10, S1b-12).
#
# - TWED with small costs and a large stiffness: TWEDdistance.cpp's border
#   stand-in 1000 * maxscost is chosen, capping the distance (traminex keeps
#   the infinite border of the TWED definition; a deliberate divergence).
# - CHI2 with step > ncol and overlap = TRUE: three intervals (seqdist-CHI2.R),
#   not one.
# - OM with an infinite indel and norm = "maxdist": maxpossiblecost is
#   abs(n - m) * indel + ..., so 0 * Inf makes equal-length distances NaN.
# - OMslen with INDELS costs and a missing-bearing reference: R halves the
#   cost before weighting it (seqdist.R:737), so a cost at the largest double
#   does not overflow.
#
# Consumed by tests/test_seqdist_numeric_edges.py.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

rows <- character()
record <- function(case, d) {
  v <- as.vector(as.matrix(d))
  rows <<- c(rows, paste(case, paste(sprintf("%.17g", v), collapse = " "), sep = "\t"))
}

ragged <- suppressMessages(seqdef(data.frame(
  V1 = c("A", "C", "A", "C"), V2 = c("B", NA, "B", "C"), V3 = c("C", NA, NA, "C"),
  V4 = c("A", NA, NA, "C"), V5 = c("B", NA, NA, "C"), V6 = c("C", NA, NA, "C"),
  V7 = c("A", NA, NA, "C"), V8 = c(NA, NA, NA, "C"), V9 = c(NA, NA, NA, "C"),
  V10 = c(NA, NA, NA, "C"), stringsAsFactors = FALSE
), alphabet = c("A", "B", "C")))
small <- matrix(c(0, 0.01, 0.02, 0.01, 0, 0.015, 0.02, 0.015, 0), nrow = 3, byrow = TRUE)
record("twed_stand_in_nu1", suppressMessages(seqdist(ragged, method = "TWED", sm = small, nu = 1, h = 0.5, norm = "none")))
record("twed_stand_in_nu5", suppressMessages(seqdist(ragged, method = "TWED", sm = small, nu = 5, h = 0.5, norm = "none")))

five <- suppressMessages(seqdef(matrix(c(
  "A", "B", "B", "C", "A",
  "B", "B", "A", "A", "C",
  "C", "A", "A", "B", "B"
), ncol = 5, byrow = TRUE), alphabet = c("A", "B", "C")))
record("chi2_step6_overlap", suppressWarnings(suppressMessages(seqdist(five, method = "CHI2", step = 6, overlap = TRUE))))
record("chi2_step6_no_overlap", suppressWarnings(suppressMessages(seqdist(five, method = "CHI2", step = 6))))
record("om_inf_indel_maxdist", suppressMessages(seqdist(five, method = "OM", sm = "CONSTANT", indel = Inf, norm = "maxdist")))

source_seqs <- suppressMessages(seqdef(matrix(c(
  "A", "B", "B", "C", "A", "A", "C", "C", "B", "B", "A", "A"
), ncol = 4, byrow = TRUE), alphabet = c("A", "B", "C")))
ref <- suppressMessages(seqdef(matrix(c("A", "B", NA, "C"), nrow = 1), alphabet = c("A", "B", "C")))
record("omslen_indels_missing_ref", suppressMessages(seqdist(
  source_seqs, method = "OMslen", sm = "INDELS", refseq = ref, with.missing = TRUE
)))

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
writeLines(c("# case\tR values, column-major", rows), file.path(out_dir, "seqdist_numeric_edges.txt"))
cat("SEQDIST_NUMERIC_EDGES_GOLDENS_DONE\n")
