# Live-R goldens for where seqpcplot places events (repair of audit findings
# TX-AUD-20261010-045 and TX-AUD-20261011-027).
#
# seqpcplot (order.align = "first") places each event at the rank of its
# distinct time within its sequence (seqpcplot.R:209-211, ordering() :985-990),
# so simultaneous events share an x, and the y position is the event's level.
# A state sequence object is drawn through its DSS (:131-140); seqdss drops
# missing cells before it collapses runs (with.missing = FALSE), so no missing
# level remains (:149-157), and the y levels are the state labels (:142;
# batch-4 review I-3). The points are those of the returned object (pts: x, y
# and the first id of each distinct trajectory); every sequence below has its
# own trajectory.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_seqpcplot_positions.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
icuSetCollate(locale = "root")
pdf(NULL)

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tse <- data.frame(id = c(1, 1, 1, 2, 2, 3, 3, 3, 4),
                  time = c(1, 1, 2, 1, 3, 2, 2, 5, 1),
                  event = c("a", "b", "c", "a", "c", "b", "c", "a", "c"))
eseq <- TraMineR::seqecreate(id = tse$id, timestamp = tse$time, event = tse$event)
sts <- suppressMessages(TraMineR::seqdef(
  data.frame(t(sapply(c("AABBC", "CCAAB", "BBBBB", "ABABA"), function(s) strsplit(s, "")[[1]]))),
  alphabet = c("A", "B", "C")))
# Internal, leading and trailing missing cells (the trailing one becomes void).
sts_missing <- suppressMessages(TraMineR::seqdef(
  data.frame(t(sapply(c("A*ABC", "*CAAB", "B*B*B", "ABAB*"), function(s) strsplit(s, "")[[1]]))),
  alphabet = c("A", "B", "C"), missing = "*", labels = c("lab-A", "lab-B", "lab-C")))
inputs <- list(events = eseq, states = sts, states_missing = sts_missing)

rows <- list()
for (case in names(inputs)) {
  result <- TraMineR::seqpcplot(inputs[[case]], plot = FALSE)
  pts <- result$pts
  stopifnot(!anyDuplicated(unique(pts[, c("traj", "id")])$id))
  for (traj in unique(pts$traj)) {
    p <- pts[pts$traj == traj, ]
    rows[[length(rows) + 1L]] <- data.frame(case = case, id = p$id[1], position = seq_len(nrow(p)),
                                            x = p$x, y = p$y,
                                            ylevel = result$ylevs[p$y], stringsAsFactors = FALSE)
  }
}
write.csv(do.call(rbind, rows), file.path(out_dir, "seqpcplot_positions.csv"), row.names = FALSE, quote = TRUE)
write.csv(tse, file.path(out_dir, "seqpcplot_positions_tse.csv"), row.names = FALSE, quote = TRUE)
