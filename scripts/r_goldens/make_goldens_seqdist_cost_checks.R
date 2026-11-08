# Golden generation: seqdist's checks on a caller-supplied substitution matrix.
#
# checkcost (TraMineR-checkcost.R, called from seqdist.R:377-381 for the OM
# family and TWED) stops on a wrong-sized matrix, a negative cost and a
# non-zero diagonal, for regular and ragged data alike.  traminex applied the
# same checks only on its dense OM path (2026-09-22 audit S1b-01).  This
# fixture records R's outcome for every (method, data shape, defect).  It also
# records what R returns for a non-positive scalar indel with a cost-method
# string, which R never checks (seqdist.R: checkcost runs for sm.type
# "matrix" only): negative distances, or an all-zero matrix.
#
# Consumed by tests/test_seqdist_cost_checks.py.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

X <- matrix(c(0, 1, 2, 1, 1, 1, 0, 2, 2, 0, 1, 1), ncol = 4, byrow = TRUE)
regular <- suppressMessages(seqdef(X, alphabet = 0:2))
Xr <- X
Xr[3, 4] <- NA
ragged <- suppressMessages(seqdef(Xr, alphabet = 0:2, right = "DEL"))

base <- matrix(2, 3, 3)
diag(base) <- 0
defects <- list(
  diagonal = { m <- base; m[2, 2] <- 1; m },
  negative = { m <- base; m[1, 3] <- -1; m[3, 1] <- -1; m },
  oversized = { m <- matrix(2, 4, 4); diag(m) <- 0; m },
  undersized = { m <- matrix(2, 2, 2); diag(m) <- 0; m }
)
method_args <- list(
  OM = list(indel = 1),
  OMloc = list(),
  OMslen = list(indel = 1),
  OMspell = list(indel = 1),
  OMstran = list(indel = 1, otto = 0.5),
  TWED = list(indel = 1, nu = 0.5, h = 0.5)
)
lines <- character()
for (shape in c("regular", "ragged")) {
  data <- if (shape == "regular") regular else ragged
  for (m in names(method_args)) {
    for (defect in names(defects)) {
      args <- c(list(seqdata = data, method = m, sm = defects[[defect]]), method_args[[m]])
      outcome <- tryCatch(
        { suppressWarnings(suppressMessages(do.call(TraMineR::seqdist, args))); "computed" },
        error = function(e) paste("error:", trimws(conditionMessage(e)))
      )
      lines <- c(lines, paste(shape, m, defect, outcome, sep = "\t"))
    }
  }
}

indel_lines <- character()
for (shape in c("regular", "ragged")) {
  data <- if (shape == "regular") regular else ragged
  for (indel in c(-1, 0)) {
    d <- suppressMessages(TraMineR::seqdist(data, method = "OM", sm = "CONSTANT", indel = indel))
    indel_lines <- c(indel_lines, paste(shape, indel, paste(sprintf("%.17g", as.vector(as.matrix(d))), collapse = ","), sep = "\t"))
  }
}

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
writeLines(c("# shape\tmethod\tdefect\tR outcome", lines), file.path(out_dir, "seqdist_cost_checks.txt"))
writeLines(
  c("# shape\tindel\tR distance matrix, column-major", indel_lines),
  file.path(out_dir, "seqdist_cost_checks_indel.txt")
)
cat("SEQDIST_COST_CHECKS_GOLDENS_DONE\n")
