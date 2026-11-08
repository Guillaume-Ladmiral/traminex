# Live-R record of seqedist on event sequences that have no events (repair of
# audit finding TX-AUD-20261010-037; R defects corrected by traminex).
#
# TraMineR's seqedist.cpp resets its indices at :156-157 and reads the result at
# :222; when the row sequence is empty its loop never runs, so the distance
# between an empty sequence and any later sequence is 0, and the matrix depends
# on the row order. For two empty sequences YujianBo divides 0 by 0 (:231). The
# fixture keeps R's values so the divergence stays visible; traminex's tests
# assert the values derived from the definition.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_seqedist_empty.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages({library(TraMineR); library(TraMineRextras)})
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14",
          as.character(packageVersion("TraMineRextras")) == "0.6.9")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

record_case <- TraMineR::seqecreate(id = c(1, 2, 2, 2, 3, 3, 4), timestamp = c(5, 0, 2, 6, 1, 4, 3),
                                    event = c("END", "X", "Y", "END", "Y", "X", "END"), end.event = "END")
fresh_case <- TraMineR::seqecreate(id = c(1, 1, 1, 2, 3, 3, 3, 4, 5, 5, 5),
                                   timestamp = c(0, 1, 3, 2, 0, 2, 4, 1, 1, 1, 3),
                                   event = c("A", "B", "END", "END", "B", "A", "END", "END", "A", "B", "END"),
                                   end.event = "END")
reordered <- TraMineR::seqecreate(id = c(1, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
                                  timestamp = c(2, 1, 0, 1, 3, 0, 2, 4, 1, 1, 3),
                                  event = c("END", "END", "A", "B", "END", "B", "A", "END", "A", "B", "END"),
                                  end.event = "END")
cases <- list(
  record = list(eseq = record_case, idcost = c(1, 1, 1), vparam = 0.1),
  fresh = list(eseq = fresh_case, idcost = c(1, 2, 0.5), vparam = 0.2),
  reordered = list(eseq = reordered, idcost = c(1, 2, 0.5), vparam = 0.2)
)

rows <- list()
for (case_name in names(cases)) {
  for (nm in c("none", "max", "YujianBo")) {
    args <- cases[[case_name]]
    value <- TraMineRextras::seqedist(args$eseq, idcost = args$idcost, vparam = args$vparam, norm = nm)
    rows[[length(rows) + 1L]] <- data.frame(
      case = case_name, norm = nm,
      row = rep(seq_len(nrow(value)), times = ncol(value)),
      col = rep(seq_len(ncol(value)), each = nrow(value)),
      value = sprintf("%.15g", as.vector(value)),
      stringsAsFactors = FALSE
    )
  }
  rows[[length(rows) + 1L]] <- data.frame(
    case = case_name, norm = "sequences", row = seq_along(cases[[case_name]]$eseq), col = 0L,
    value = as.character(cases[[case_name]]$eseq), stringsAsFactors = FALSE
  )
}
write.csv(do.call(rbind, rows), file.path(out_dir, "seqedist_empty.csv"), row.names = FALSE, quote = TRUE)
