# Live-R goldens for how group levels are rendered as text (repair of audit
# findings TX-AUD-20261010-026 and -049).
#
# seqdiff names its discrepancy columns after factor(group) levels, and
# TraMineRextras' group.p and rowmode label with as.character: FALSE/TRUE,
# 1e+05, 1e+15, 0.3 for 0.1 + 0.2, 0 for -0, 1e-20, 123456; factor() drops an
# unused level of a factor.
#
# factor() keeps NaN as a level of its own, so seqdiff gets a NaN column and
# group.p a NaN level, while table() (rowmode) leaves NaN out. traminex reads a
# NaN group as missing; the nan_as_na cases are the same groups with NA for NaN.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_group_label_rendering.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages({library(TraMineR); library(TraMineRextras)})
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14",
          as.character(packageVersion("TraMineRextras")) == "0.6.9")
icuSetCollate(locale = "root")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

strings <- c("AABB", "ABBB", "AAAB", "BBBA", "ABAB", "BAAA", "AABA", "BBAB")
seqs <- suppressMessages(TraMineR::seqdef(as.data.frame(do.call(rbind, strsplit(strings, ""))), alphabet = c("A", "B")))
groups <- list(
  logical = rep(c(TRUE, FALSE), 4),
  doubles = rep(c(100000, 2, 0.1, 1e15), 2),
  float_sum = rep(c(0.1 + 0.2, 2), 4),
  unused_level = factor(rep(c("b", "a"), 4), levels = c("a", "b", "zzz_unused")),
  edge_numbers = rep(c(-0, 1e-20, 123456, 2), 2)
)
rows <- list()
for (name in names(groups)) {
  columns <- colnames(suppressMessages(TraMineR::seqdiff(seqs, group = groups[[name]]))$discrepancy)
  rows[[length(rows) + 1L]] <- data.frame(function_name = "seqdiff", case = name, position = seq_along(columns),
                                          label = columns, stringsAsFactors = FALSE)
}
nan_groups <- list(nan_level = c(NaN, 1, 2, NA, 1, 2, NaN, 1), nan_as_na = c(NA, 1, 2, NA, 1, 2, NA, 1))
for (name in names(nan_groups)) {
  columns <- colnames(suppressMessages(TraMineR::seqdiff(seqs, group = nan_groups[[name]]))$discrepancy)
  rows[[length(rows) + 1L]] <- data.frame(function_name = "seqdiff", case = name, position = seq_along(columns),
                                          label = columns, stringsAsFactors = FALSE)
  levels_p <- levels(TraMineRextras::group.p(nan_groups[[name]]))
  rows[[length(rows) + 1L]] <- data.frame(function_name = "group.p", case = name, position = seq_along(levels_p),
                                          label = levels_p, stringsAsFactors = FALSE)
}
for (name in c("doubles", "edge_numbers", "logical")) {
  levels_p <- levels(TraMineRextras::group.p(groups[[name]]))
  rows[[length(rows) + 1L]] <- data.frame(function_name = "group.p", case = name, position = seq_along(levels_p),
                                          label = levels_p, stringsAsFactors = FALSE)
}
modes <- list(big = c(1e5, 1e5, 2), huge = c(1e15, 1e15, 3), tiny = c(1e-20, 1e-20, 2), plain = c(123456, 123456, 1),
              truth = c(TRUE, TRUE, FALSE), nan = c(NaN, NaN, 1))
for (name in names(modes)) {
  rows[[length(rows) + 1L]] <- data.frame(function_name = "rowmode", case = name, position = 1L,
                                          label = TraMineRextras::rowmode(modes[[name]]), stringsAsFactors = FALSE)
}
write.csv(do.call(rbind, rows), file.path(out_dir, "group_label_rendering.csv"), row.names = FALSE, quote = TRUE)
