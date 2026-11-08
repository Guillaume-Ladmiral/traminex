# Golden generation: accepted asymmetric substitution costs with duplicated rows.
#
# seqdist removes duplicated sequences before its kernels run (seqdist.R:524)
# and checkcost warns "Non-symmetric substitution cost matrix!".  With
# asymmetric costs the direction of each comparison is observable, so the
# deduplicated, first-occurrence order decides every value: identical
# sequences always get identical distances.  traminex evaluated OMloc, OMslen,
# OMspell, OMstran and TWED pairs in row order, so identical sequences could
# get different distances (2026-09-22 audit S1b-03).
#
# Consumed by tests/test_seqdist_asymmetric_duplicates.py.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

sequences <- suppressMessages(seqdef(
  matrix(c("A", "A", "B", "C",
           "C", "C", "B", "A",
           "A", "A", "B", "C",
           "B", "A", "C", "C",
           "C", "C", "B", "A"), ncol = 4, byrow = TRUE),
  alphabet = c("A", "B", "C")
))
sm <- matrix(c(0, 1, 3,
               2, 0, 1,
               1, 4, 0), 3, byrow = TRUE)
method_args <- list(
  OMloc = list(),
  OMslen = list(indel = 1.5),
  OMspell = list(indel = 1.5),
  OMstran = list(indel = 1.5, otto = 0.5),
  TWED = list(indel = 1.5, nu = 0.5, h = 0.5)
)
rows <- character()
warned <- character()
for (m in names(method_args)) {
  seen <- character()
  d <- withCallingHandlers(
    suppressMessages(do.call(TraMineR::seqdist, c(list(seqdata = sequences, method = m, sm = sm), method_args[[m]]))),
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  d <- as.matrix(d)
  for (i in seq_len(nrow(d))) {
    for (j in seq_len(ncol(d))) {
      rows <- c(rows, paste(m, i - 1, j - 1, sprintf("%.17g", d[i, j]), sep = ","))
    }
  }
  warned <- c(warned, paste(m, paste(unique(trimws(seen)), collapse = " | "), sep = "\t"))
}

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
writeLines(c("method,row,col,value", rows), file.path(out_dir, "seqdist_asymmetric_duplicates.csv"))
writeLines(c("# method\tR warnings", warned), file.path(out_dir, "seqdist_asymmetric_duplicates_warnings.txt"))
cat("SEQDIST_ASYMMETRIC_DUPLICATES_GOLDENS_DONE\n")
