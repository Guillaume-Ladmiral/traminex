# Live-R goldens for the row order seqplotMD draws in each domain panel with a
# string sortv (repair of audit finding TX-AUD-20261010-008).
#
# With sortv "from.start"/"from.end", the default type = "i" and the default
# dom.crit = 0, seqplotMD orders the multidomain sequences
# seqMD(channels, what = "MDseq") once and passes that order to every panel's
# plot.stslist (seqplotMD.R:205-258, 486-489, 589-598), which draws
# x[order(sortv), ] (plot.stslist.R:63-75). The
# drawn order exists only inside that call, so it is recorded by tracing the
# non-exported TraMineR:::plot.stslist and taking order(sortv) at entry; a
# private access, needed because no public value carries the order.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_seqplotmd_sort_order.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
stopifnot(capabilities("ICU"))
invisible(icuSetCollate(locale = "root"))
stopifnot(identical(icuGetCollate(), "root"))

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

domain1 <- suppressMessages(TraMineR::seqdef(matrix(
  c("A", "A", "B", "B", "A", "A", "A", "B", "B", "B", "B", "A", "A", "A", "A"), ncol = 3, byrow = TRUE
)))
domain2 <- suppressMessages(TraMineR::seqdef(matrix(
  c("y", "x", "x", "x", "x", "y", "x", "y", "y", "y", "y", "x", "x", "x", "y"), ncol = 3, byrow = TRUE
)))

captured <- list()
trace(TraMineR:::plot.stslist, where = asNamespace("TraMineR"), print = FALSE,
      tracer = quote(captured[[length(captured) + 1L]] <<- order(sortv)))

rows <- list()
for (sortv in c("from.start", "from.end")) {
  captured <- list()
  pdf(NULL)
  suppressMessages(TraMineR::seqplotMD(list(domain1, domain2), sortv = sortv,
                                       with.legend = FALSE))
  invisible(dev.off())
  stopifnot(length(captured) == 2L)
  for (panel in seq_along(captured)) {
    rows[[length(rows) + 1L]] <- data.frame(
      case = sortv, panel = panel, position = seq_along(captured[[panel]]),
      row = captured[[panel]], stringsAsFactors = FALSE
    )
  }
}
untrace(TraMineR:::plot.stslist, where = asNamespace("TraMineR"))

write.csv(do.call(rbind, rows), file.path(out_dir, "seqplotmd_sort_order.csv"),
          row.names = FALSE, quote = TRUE)
