# Live-R goldens for where an index plot draws each sequence (repair of audit
# finding TX-AUD-20261010-059).
#
# plot.stslist's index plot is barplot(seqbar, horiz = TRUE) with one column per
# sequence (plot.stslist.R:131,149-151), which draws the first sequence at the
# bottom. The bar midpoints barplot returns are recorded (traced on exit), in the
# order of the plotted sequences, for seqIplot with all sequences shown.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_index_orientation.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
pdf(NULL)

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

seqs <- suppressMessages(TraMineR::seqdef(
  data.frame(t(sapply(c("AABBC", "ABBBC", "BBAAC", "CCCBA"), function(s) strsplit(s, "")[[1]]))),
  alphabet = c("A", "B", "C")))
midpoints <- NULL
trace(graphics::barplot, print = FALSE, where = asNamespace("TraMineR"),
      exit = quote(midpoints <<- returnValue()))
TraMineR::seqIplot(seqs, with.legend = FALSE)
untrace(graphics::barplot, where = asNamespace("TraMineR"))
stopifnot(length(midpoints) == nrow(seqs))
write.csv(data.frame(sequence = seq_len(nrow(seqs)), midpoint = sprintf("%.17g", as.vector(midpoints))),
          file.path(out_dir, "index_orientation.csv"), row.names = FALSE, quote = TRUE)
