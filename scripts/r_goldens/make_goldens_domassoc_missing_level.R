# Live-R goldens for seqdomassoc's missing level (repair of audit finding
# TX-AUD-20261011-013).
#
# seqdomassoc's position-wise table has every level of the stslist factors and
# then removes void and, without with.missing, the missing state
# (seqdomassoc.R:206-208), so with with.missing = TRUE the missing row and column
# are kept for every domain, also one without missing cells; with struct.zero,
# its zero cells then enter the degrees-of-freedom adjustment.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_domassoc_missing_level.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
icuSetCollate(locale = "root")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

with_na <- suppressMessages(TraMineR::seqdef(
  data.frame(a = c("A", "B", NA, "A", "B", "A"), b = c("A", NA, "B", "B", "B", "A"),
             c = c("B", "B", "A", "A", "B", "B")), alphabet = c("A", "B")))
without_na <- suppressMessages(TraMineR::seqdef(
  data.frame(a = c("X", "Y", "Y", "X", "Y", "X"), b = c("Y", "Y", "X", "X", "Y", "X"),
             c = c("X", "Y", "Y", "Y", "X", "X")), alphabet = c("X", "Y")))
also_na <- suppressMessages(TraMineR::seqdef(
  data.frame(a = c("X", NA, "Y", "X", "Y", "X"), b = c("Y", "Y", "X", NA, "Y", "X"),
             c = c("X", "Y", "Y", "Y", "X", "X")), alphabet = c("X", "Y")))
cases <- list(one_domain_without_missing = list(with_na, without_na),
              both_domains_with_missing = list(with_na, also_na))
rows <- list()
for (name in names(cases)) {
  for (struct in c(TRUE, FALSE)) {
    result <- suppressMessages(TraMineR::seqdomassoc(cases[[name]], with.missing = TRUE, cross.table = TRUE,
                                                     struct.zero = struct))
    stats <- unclass(result)
    for (column in colnames(stats)) {
      rows[[length(rows) + 1L]] <- data.frame(case = name, struct_zero = struct, section = "stat", name = column,
                                              value = sprintf("%.17g", stats[1, column]), stringsAsFactors = FALSE)
    }
    table <- attr(result, "cross.tables")[[1]]
    for (i in seq_len(nrow(table))) for (j in seq_len(ncol(table))) {
      rows[[length(rows) + 1L]] <- data.frame(case = name, struct_zero = struct, section = "table",
                                              name = paste(rownames(table)[i], colnames(table)[j], sep = "|"),
                                              value = sprintf("%.17g", table[i, j]), stringsAsFactors = FALSE)
    }
  }
}
write.csv(do.call(rbind, rows), file.path(out_dir, "domassoc_missing_level.csv"), row.names = FALSE, quote = TRUE)
