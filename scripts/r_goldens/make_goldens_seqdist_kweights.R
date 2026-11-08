# Golden generation: seqdist's kweights contract for NMS, NMSMST and SVRspell.
#
# seqdist.R:265-273 refuses a non-numeric or negative kweights (a NaN makes the
# check itself fail), recycles a scalar with a warning, and warns when the
# vector is shorter than the sequences.  An infinite weight gives NaN
# distances.  traminex crashed on a scalar, accepted negative weights, and
# clamped NaN squared distances to 0 (2026-09-22 audit S1b-06).
#
# Consumed by tests/test_seqdist_kweights.py.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

sequences <- suppressMessages(seqdef(
  matrix(c("A", "A", "B", "C", "C",
           "A", "B", "B", "C", "A",
           "C", "C", "A", "A", "B",
           "B", "A", "C", "C", "C"), ncol = 5, byrow = TRUE),
  alphabet = c("A", "B", "C")
))
cases <- list(
  scalar = 0.5,
  short = c(1, 0.5),
  full = c(1, 0.8, 0.6, 0.4, 0.2),
  infinite = c(Inf, 1, 1, 1, 1),
  negative = c(1, -0.5, 1, 1, 1),
  not_a_number = c(NaN, 1, 1, 1, 1)
)
rows <- character()
outcomes <- character()
for (m in c("NMS", "NMSMST", "SVRspell")) {
  for (case in names(cases)) {
    seen <- character()
    d <- tryCatch(
      withCallingHandlers(
        suppressMessages(TraMineR::seqdist(sequences, method = m, kweights = cases[[case]])),
        warning = function(w) {
          seen <<- c(seen, trimws(conditionMessage(w)))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) e
    )
    if (inherits(d, "error")) {
      outcomes <- c(outcomes, paste(m, case, "error", trimws(conditionMessage(d)), sep = "\t"))
      next
    }
    outcomes <- c(outcomes, paste(m, case, "warnings", if (length(seen)) paste(unique(seen), collapse = " | ") else "none", sep = "\t"))
    d <- as.matrix(d)
    for (i in seq_len(nrow(d))) {
      for (j in seq_len(ncol(d))) {
        rows <- c(rows, paste(m, case, i - 1, j - 1, sprintf("%.17g", d[i, j]), sep = ","))
      }
    }
  }
}

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
writeLines(c("method,case,row,col,value", rows), file.path(out_dir, "seqdist_kweights_cases.csv"))
writeLines(c("# method\tcase\tkind\tR messages", outcomes), file.path(out_dir, "seqdist_kweights_outcomes.txt"))
cat("SEQDIST_KWEIGHTS_GOLDENS_DONE\n")
