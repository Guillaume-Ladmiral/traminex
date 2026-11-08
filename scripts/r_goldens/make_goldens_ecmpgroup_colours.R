# Live-R goldens for the bar colours of plot.subseqelistchisq (repair of audit
# finding TX-AUD-20261010-045).
#
# plot.subseqelistchisq colours each bar by cut(residual, c(-Inf, -l2, -l1, l1,
# l2, Inf), labels = brewer.pal(5, "RdBu")), with l = -qnorm(c(0.05, 0.01)),
# Bonferroni-adjusted when seqecmpgroup used it, so a positive residual is blue,
# and titles each panel with its group label (seqecmpgroup.R:105-142). The
# colours and titles are recorded as drawn, by tracing graphics::barplot.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_ecmpgroup_colours.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
icuSetCollate(locale = "root")
pdf(NULL)

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 40 event sequences: in group g1 most have a then b, in g2 most have c then b.
ids <- integer(0); times <- numeric(0); events <- character(0)
for (i in 1:40) {
  first <- if ((i <= 20 && i %% 10 != 0) || (i > 20 && i %% 10 == 0)) "a" else "c"
  ids <- c(ids, i, i); times <- c(times, 1, 2); events <- c(events, first, "b")
  if (i %% 3 == 0) { ids <- c(ids, i); times <- c(times, 3); events <- c(events, "d") }
}
tse <- data.frame(id = ids, time = times, event = events)
group <- rep(c("g1", "g2"), each = 20)
eseq <- TraMineR::seqecreate(id = tse$id, timestamp = tse$time, event = tse$event)
fsub <- TraMineR::seqefsub(eseq, pmin.support = 0.05)

drawn <- list()
trace(graphics::barplot, print = FALSE, where = asNamespace("TraMineR"),
      tracer = quote(drawn[[length(drawn) + 1L]] <<- list(main = list(...)$main, col = list(...)$col)))
rows <- list()
for (method in c("chisq", "bonferroni")) {
  drawn <- list()
  cmp <- TraMineR::seqecmpgroup(fsub, group = group, method = method)
  plot(cmp)
  for (i in seq_along(drawn)) {
    rows[[length(rows) + 1L]] <- data.frame(method = method, panel = i, title = drawn[[i]]$main,
                                            position = seq_along(drawn[[i]]$col),
                                            subseq = as.character(cmp$subseq), colour = drawn[[i]]$col,
                                            stringsAsFactors = FALSE)
  }
}
untrace(graphics::barplot, where = asNamespace("TraMineR"))
write.csv(do.call(rbind, rows), file.path(out_dir, "ecmpgroup_colours.csv"), row.names = FALSE, quote = TRUE)
write.csv(tse, file.path(out_dir, "ecmpgroup_colours_tse.csv"), row.names = FALSE, quote = TRUE)
