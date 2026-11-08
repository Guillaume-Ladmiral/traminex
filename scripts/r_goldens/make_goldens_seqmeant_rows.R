# Live-R goldens for the bars of seqmtplot (repair of audit findings
# TX-AUD-20261010-043 and -058).
#
# seqmtplot is seqplot(type = "mt"), which calls seqmeant (seqplot.R:345-346).
# seqmeant's serr defaults to FALSE (seqmeant.R:5), so no error bars are drawn
# unless asked (plot.stslist.meant.R:59-65), and its rows, one bar each, are the
# alphabet, plus the missing state with with.missing = TRUE. The seqmeant values
# themselves are R-golden tested elsewhere (make_goldens_seqmeant*.R).
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_seqmeant_rows.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
icuSetCollate(locale = "root")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ragged <- suppressMessages(TraMineR::seqdef(data.frame(a = c("A", "A"), b = c(NA, "B"), c = c("B", NA),
                                                       d = c("B", NA))))
complete <- suppressMessages(TraMineR::seqdef(data.frame(a = c("A", "B"), b = c("B", "B"), c = c("C", "A"))))
rows <- list()
for (name in c("ragged", "complete")) {
  seqs <- get(name)
  for (with_missing in c(FALSE, TRUE)) {
    result <- TraMineR::seqmeant(seqs, with.missing = with_missing)
    rows[[length(rows) + 1L]] <- data.frame(case = name, with_missing = with_missing,
                                            position = seq_len(nrow(result)), row = rownames(result),
                                            serr_default = as.character(formals(TraMineR::seqmeant)$serr),
                                            se_attribute = as.character(attr(result, "se")),
                                            stringsAsFactors = FALSE)
  }
}
write.csv(do.call(rbind, rows), file.path(out_dir, "seqmeant_rows.csv"), row.names = FALSE, quote = TRUE)
