# Live-R goldens for sorting sequences that hold missing and void cells (repair
# of audit findings TX-AUD-20261011-006 and -014).
#
# R orders cells over the factor levels c(alphabet, nr, void) (seqdef.R:158):
# plot.stslist (plot.stslist.R:62-67) and seqrf (seqrf.R:11-16) sort with
# do.call(order, ...) over the state columns; TraMineRextras' sorti/sortv first
# replace void by nr (sortv.R:5). The index-plot order is the permutation
# plot.stslist draws, captured by an exit tracer (private access: the drawn
# order is observable only inside plot.stslist).
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_sort_marker_layouts.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages({library(TraMineR); library(TraMineRextras)})
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14",
          as.character(packageVersion("TraMineRextras")) == "0.6.9")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cells <- rbind(c("B", NA, "A"), c(NA, "A", "C"), c("A", "C", NA), c("C", NA, NA), c("A", "A", "B"), c(NA, NA, "B"))
seqs <- suppressMessages(TraMineR::seqdef(cells, alphabet = c("A", "B", "C")))
stopifnot(identical(levels(seqs[, 1]), c("A", "B", "C", "*", "%")))

captured <- NULL
trace(TraMineR:::plot.stslist, where = asNamespace("TraMineR"), print = FALSE,
      exit = quote(captured <<- sortv))
drawn <- list()
for (sortv in c("from.start", "from.end")) {
  captured <- NULL
  pdf(NULL)
  suppressMessages(TraMineR::seqIplot(seqs, sortv = sortv, with.legend = FALSE))
  invisible(dev.off())
  stopifnot(length(captured) == nrow(seqs))
  drawn[[sortv]] <- captured
}
untrace(TraMineR:::plot.stslist, where = asNamespace("TraMineR"))

record <- function(case, values) {
  data.frame(case = case, position = seq_along(values), value = as.character(values), stringsAsFactors = FALSE)
}
rows <- list(
  record("cells", as.vector(t(as.matrix(seqs)))),
  record("sorti_beg", TraMineRextras::sorti(seqs, start = "beg")),
  record("sorti_end", TraMineRextras::sorti(seqs, start = "end")),
  record("sortv_beg", TraMineRextras::sortv(seqs, start = "beg")),
  record("sortv_end", TraMineRextras::sortv(seqs, start = "end")),
  record("plot_from_start", drawn[["from.start"]]),
  record("plot_from_end", drawn[["from.end"]])
)
write.csv(do.call(rbind, rows), file.path(out_dir, "sort_marker_layouts.csv"), row.names = FALSE, quote = TRUE)

# seqrf from.start / from.end on mvad[1:60, 15:86] with missing cells (no void).
data(mvad, package = "TraMineR")
mv <- as.matrix(mvad[1:60, 15:86])
for (i in seq_len(nrow(mv))) for (j in seq_len(ncol(mv))) if ((i * 7 + j * 3) %% 29 == 0) mv[i, j] <- NA
alphabet <- c("employment", "FE", "HE", "joblessness", "school", "training")
smv <- suppressMessages(TraMineR::seqdef(mv, alphabet = alphabet, right = NA))
dmv <- suppressMessages(TraMineR::seqdist(smv, method = "HAM", with.missing = TRUE))
rf_rows <- list()
for (sortv in c("from.start", "from.end")) {
  rf <- suppressMessages(TraMineR::seqrf(smv, diss = dmv, k = 10, sortv = sortv, weighted = FALSE))
  rf_rows[[length(rf_rows) + 1L]] <- data.frame(case = sortv, kind = "R2", position = 1L,
                                                value = sprintf("%.15g", rf$rf$R2), stringsAsFactors = FALSE)
  rf_rows[[length(rf_rows) + 1L]] <- data.frame(case = sortv, kind = "medoid", position = seq_along(rf$rf$medoids),
                                                value = as.character(rf$rf$medoids), stringsAsFactors = FALSE)
}
write.csv(do.call(rbind, rf_rows), file.path(out_dir, "sort_marker_layouts_seqrf.csv"), row.names = FALSE, quote = TRUE)
# The seqrf inputs as R holds them: the cells ("" for NA) and R's HAM matrix.
write.csv(mv, file.path(out_dir, "sort_marker_layouts_seqrf_cells.csv"), row.names = FALSE, na = "", quote = TRUE)
write.csv(dmv, file.path(out_dir, "sort_marker_layouts_seqrf_diss.csv"), row.names = FALSE)
