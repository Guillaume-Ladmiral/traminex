# Live-R goldens for seqidegrad with a logical pow (repair of audit finding
# TX-AUD-20261010-004).
#
# R's public seqidegrad has no spell.integr argument: a logical pow sets it and
# is then coerced to numeric (seqidegrad.R:10-15), so pow=FALSE computes the
# index without spell integration.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_seqidegrad_logical_pow.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rows_in <- lapply(c("AABBC", "CCBAA", "ABCBA", "BBBCC"), function(s) strsplit(s, "")[[1]])
seqs <- suppressMessages(TraMineR::seqdef(
  as.data.frame(do.call(rbind, rows_in), stringsAsFactors = FALSE), alphabet = c("A", "B", "C")
))

rows <- list()
for (method in c("RANK", "TRATEDSS", "FREQ")) {
  for (pow in c(FALSE, TRUE)) {
    value <- as.vector(TraMineR::seqidegrad(seqs, method = method, pow = pow))
    rows[[length(rows) + 1L]] <- data.frame(
      method = method, pow = if (pow) "TRUE" else "FALSE", row = seq_along(value),
      value = sprintf("%.15g", value), stringsAsFactors = FALSE
    )
  }
}
write.csv(do.call(rbind, rows), file.path(out_dir, "seqidegrad_logical_pow.csv"),
          row.names = FALSE, quote = TRUE)
