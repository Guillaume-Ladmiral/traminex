# Golden generation: 2026-09-22 whole-codebase audit, package s4b (follow-up
# findings S4-16 and S4-17 on dissimilarity trees and representatives).
#
# Each section reproduces one finding against TraMineR 2.2-14 and writes one
# fixture under tests/r_goldens/; the consuming test is named in the section
# header. Values are written with sprintf("%.17g"); R's NA is written as "NA",
# and an R error is recorded as its message.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

num <- function(x) paste(sprintf("%.17g", as.vector(x)), collapse = " ")
write_fixture <- function(name, header, lines) {
  writeLines(c(paste("#", header), lines), file.path(out_dir, name))
}
one_line <- function(text) gsub("\\s+", " ", trimws(text))

## --------------------------------------------------------------------------
## S4-16: disstree compares candidate splits with if (SCres > currentSCres)
## (disstree.R:476); a missing dissimilarity off the diagonal, a missing weight
## or a predictor level of weight 0 makes that comparison NA and R stops. With
## no comparison (max.depth = 1, one level only) R returns the root alone; the
## diagonal is never read by the C kernels, so a missing value there is ignored.
## Consumed by tests/test_disstree_missing_values.py.
## --------------------------------------------------------------------------
pos <- c(0.0, 0.1, 0.2, 0.5, 0.6, 0.7, 5.0, 5.1, 5.2, 5.5, 5.6, 5.7)
D0 <- as.matrix(dist(pos))
w0 <- c(1.5, 0.5, 2.0, 1.0, 0.7, 1.3, 2.2, 0.9, 1.1, 1.6, 0.4, 1.8)
two_levels <- rep(c("a", "b"), each = 6)
three_levels <- c(rep("a", 6), rep("b", 4), rep("c", 2))
with_value <- function(D, i, j, value) {
  D[i, j] <- value
  D[j, i] <- value
  D
}
nan_weight <- w0
nan_weight[2] <- NaN
zero_c <- w0
zero_c[11:12] <- 0
cases <- list(
  nan_off_diagonal = list(D = with_value(D0, 2, 9, NaN), v = two_levels, w = NULL, depth = 5),
  na_off_diagonal = list(D = with_value(D0, 2, 9, NA), v = two_levels, w = NULL, depth = 5),
  nan_off_diagonal_weighted = list(D = with_value(D0, 2, 9, NaN), v = two_levels, w = w0, depth = 5),
  nan_weight = list(D = D0, v = two_levels, w = nan_weight, depth = 5),
  zero_weight_level = list(D = D0, v = three_levels, w = zero_c, depth = 5),
  nan_off_diagonal_depth_one = list(D = with_value(D0, 2, 9, NaN), v = two_levels, w = NULL, depth = 1),
  nan_off_diagonal_one_level = list(D = with_value(D0, 2, 9, NaN), v = rep("a", 12), w = NULL, depth = 5),
  nan_diagonal = list(D = with_value(D0, 3, 3, NaN), v = two_levels, w = NULL, depth = 5),
  inf_off_diagonal = list(D = with_value(D0, 2, 9, Inf), v = two_levels, w = NULL, depth = 5)
)
lines <- paste("pos", "values", num(pos), sep = "\t")
for (case in names(cases)) {
  k <- cases[[case]]
  lines <- c(lines,
    paste(case, "D", num(t(k$D)), sep = "\t"),
    paste(case, "v", paste(k$v, collapse = " "), sep = "\t"),
    paste(case, "weights", if (is.null(k$w)) "none" else num(k$w), sep = "\t"),
    paste(case, "max_depth", k$depth, sep = "\t"))
  outcome <- tryCatch({
    tr <- disstree(k$D ~ v, data = data.frame(v = k$v), weights = k$w, R = 1, min.size = 0,
      max.depth = k$depth, weight.permutation = "diss")
    c(paste(case, "outcome", "tree", sep = "\t"),
      paste(case, "fitted", paste(tr$fitted[, 1], collapse = " "), sep = "\t"),
      paste(case, "root_vardis", num(tr$root$info$vardis), sep = "\t"),
      paste(case, "adjustment_t0", num(tr$info$adjustment$stat[, 1]), sep = "\t"))
  }, error = function(e) paste(case, "outcome", "error", one_line(conditionMessage(e)), sep = "\t"))
  lines <- c(lines, outcome)
}
write_fixture(
  "disstree_missing_values.txt",
  "disstree of D on v with R = 1 and min.size = 0; missing or infinite dissimilarities and weights",
  lines
)

## --------------------------------------------------------------------------
## S4-17: dissrep's mean distance is MD <- SD/na (dissrep.R:151-153), a plain
## division, so a representative whose assigned cases all weigh 0 has MD = 0/0.
## Consumed by tests/test_dissrep_zero_weight_representative.py.
## --------------------------------------------------------------------------
rep_pos <- c(0, 0.05, 0.1, 10)
rep_w <- c(1, 1, 1, 0)
r <- suppressMessages(dissrep(as.matrix(dist(rep_pos)), criterion = "density", weights = rep_w, nrep = 2))
write_fixture(
  "dissrep_zero_weight_representative.txt",
  "dissrep with nrep = 2 where the second representative carries weight 0; statistics column-major",
  c(paste("positions", num(rep_pos), sep = "\t"),
    paste("weights", num(rep_w), sep = "\t"),
    paste("representatives", paste(as.vector(r), collapse = " "), sep = "\t"),
    paste("statistics", num(as.matrix(attr(r, "Statistics"))), sep = "\t"),
    paste("quality", num(attr(r, "Quality")), sep = "\t"))
)

cat("AUDIT20260922_S4B_GOLDENS_DONE\n")
