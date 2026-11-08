# Live-R goldens for how far R draws a survival curve (repair of audit findings
# TX-AUD-20261010-042 and TX-AUD-20261011-026).
#
# seqsplot (plot.stslist.surv.R:28-30) and seqedplot's survival type
# (seqedplot.R:201-235) draw survfit objects with plot.survfit: each curve starts
# at (0, 1), steps at its event times and runs on to its last time, censored or
# not, so a state whose spells are all censored is a flat line at 1. The
# coordinates are recorded as drawn, by tracing graphics::lines while the plots
# are made; a tracer on plot.stslist.surv numbers the group panels.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_survival_extent.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages({library(TraMineR); library(TraMineRextras)})
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14",
          as.character(packageVersion("TraMineRextras")) == "0.6.9")
icuSetCollate(locale = "root")
pdf(NULL)

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

data(mvad)
alph <- sort(unique(as.vector(as.matrix(mvad[, 15:86]))))
mv <- suppressMessages(TraMineR::seqdef(mvad[1:150, 15:86], alphabet = alph))
toy <- suppressMessages(TraMineR::seqdef(
  data.frame(t(sapply(c("AABBC", "ABBBC", "BBAAC", "AABBB", "BBBAA"), function(s) strsplit(s, "")[[1]]))),
  alphabet = c("A", "B", "C")))
tse <- data.frame(id = c("s1", "s1", "s2", "s3", "s3", "s4", "s4", "s5"),
                  timestamp = c(1, 3, 2, 1, 4, 1, 3, 2), event = c("A", "B", "A", "B", "A", "A", "A", "B"))
eseq <- TraMineR::seqecreate(id = tse$id, timestamp = tse$timestamp, event = tse$event)
TraMineR::seqelength(eseq) <- c(5, 5, 6, 4, 5)

drawn <- list()
panel <- 0L
# plot.survfit draws each curve with lines() of an x/y list; lines() of a survfit
# object (from lines.survfit) is skipped, its own inner lines() call is recorded.
trace(graphics::lines, print = FALSE, where = asNamespace("survival"),
      tracer = quote({
        xy <- if (is.list(x)) x else list(x = x, y = y)
        if (is.numeric(xy$x)) drawn[[length(drawn) + 1L]] <<- list(panel = panel, x = xy$x, y = xy$y)
      }))
trace("plot.stslist.surv", print = FALSE, where = asNamespace("TraMineRextras"),
      tracer = quote(panel <<- panel + 1L))
cases <- list(
  toy = quote(TraMineRextras::seqsplot(toy)),
  mvad150 = quote(TraMineRextras::seqsplot(mv)),
  mvad150_male = quote(TraMineRextras::seqsplot(mv, group = mvad$male[1:150])),
  seqedplot = quote(TraMineRextras::seqedplot(eseq, type = "survival"))
)
rows <- list()
for (name in names(cases)) {
  drawn <- list()
  panel <- 0L
  eval(cases[[name]])
  for (i in seq_along(drawn)) {
    rows[[length(rows) + 1L]] <- data.frame(case = name, panel = max(drawn[[i]]$panel, 1L), curve = i,
                                            x = paste(sprintf("%.17g", drawn[[i]]$x), collapse = "|"),
                                            y = paste(sprintf("%.17g", drawn[[i]]$y), collapse = "|"),
                                            stringsAsFactors = FALSE)
  }
}
untrace(graphics::lines, where = asNamespace("survival"))
untrace("plot.stslist.surv", where = asNamespace("TraMineRextras"))
write.csv(do.call(rbind, rows), file.path(out_dir, "survival_extent.csv"), row.names = FALSE, quote = TRUE)
write.csv(mvad[1:150, c("male", names(mvad)[15:86])], file.path(out_dir, "survival_extent_mvad.csv"),
          row.names = FALSE, quote = TRUE)
