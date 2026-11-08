# Golden generation: plotting and display semantics (2026-09-22 audit, slice S6).
#
# Four small oracles for the plotting findings of the 2026-09-22 audit:
#
# * seqpcplot_filter_trajectories.csv -- which event sequences seqpcplot
#   colours under a seqpcfilter.  seqpcplot.R:633-639 evaluates the filter on
#   propid.tot, one frequency per DISTINCT trajectory (events grouped by
#   timestamp rank) plus a final "no observation" share (seqpcplot.R:517);
#   traminex filtered a per-sequence vector (S6-02).  A sequence is "drawn"
#   when its trajectory keeps its colour, i.e. is not painted with hide.col.
#   The cases are built so that no tie at a cumfreq cut-off is decided by
#   punctuation, where R's order depends on the session collation; LC_COLLATE
#   is pinned to C anyway so the file does not depend on the machine.
# * r_format_digits.csv -- format(x, digits = d), the rule disstree2dot.R
#   uses for DOT node labels (formd, digits = 3) and for the quality lines
#   (getOption("digits") - 2 = 5).  traminex used "%.{d}g" (S6-03).
# * seqmsplot_info_text.tsv -- the info text plot.stslist.modst.R draws:
#   nbrep / attr(x, "nbseq"), where nbseq is the SUM OF WEIGHTS (seqmodst.R:65)
#   and R pastes round(., 1) as text (S6-05).
# * plot_group_levels.tsv -- the panel order and missing-value handling of the
#   group factors the plots build: factor() in seqedplot.R:53 (S6-06) and
#   TraMineR:::group() in seqplot.R/seqsplot.R (S6-07, S6-08).
#
# Consumed by tests/test_seqpcplot_trajectory_filter.py,
# tests/test_disstree_dot_number_format.py, tests/test_seqmsplot_weighted_info.py
# and tests/test_seqedplot_group_panels.py.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
invisible(Sys.setlocale("LC_COLLATE", "C"))

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --------------------------------------------------------------------------- #
# 1. seqpcplot filter over distinct trajectories                              #
# --------------------------------------------------------------------------- #
pc_cases <- list(
  shares_602020 = list(
    id = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5),
    timestamp = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1),
    event = c("A", "B", "A", "B", "A", "B", "A", "C", "B", "C"),
    weights = c(1, 1, 1, 1, 1)
  ),
  first_appearance_tie = list(
    id = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5),
    timestamp = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1),
    event = c("B", "C", "A", "B", "A", "B", "A", "B", "A", "C"),
    weights = c(1, 1, 1, 1, 1)
  ),
  weighted = list(
    id = c(1, 1, 2, 2, 3, 3, 4, 4),
    timestamp = c(0, 1, 0, 1, 0, 1, 0, 1),
    event = c("A", "B", "A", "B", "B", "C", "A", "C"),
    weights = c(1, 1, 3, 0.5)
  ),
  simultaneous = list(
    id = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5),
    timestamp = c(0, 0, 0, 1, 0, 1, 0, 0, 0, 1),
    event = c("A", "B", "A", "B", "A", "B", "A", "B", "A", "B"),
    weights = c(1, 1, 1, 1, 1)
  )
)
pc_filters <- list(
  list(case = "shares_602020", method = "cumfreq", level = 0.5),
  list(case = "shares_602020", method = "cumfreq", level = 0.7),
  list(case = "shares_602020", method = "linear", level = NA),
  list(case = "first_appearance_tie", method = "cumfreq", level = 0.7),
  list(case = "first_appearance_tie", method = "minfreq", level = 0.2),
  list(case = "weighted", method = "cumfreq", level = 0.5),
  list(case = "weighted", method = "cumfreq", level = 0.8),
  list(case = "simultaneous", method = "cumfreq", level = 0.5),
  list(case = "simultaneous", method = "minfreq", level = 0.5)
)

# Trajectory key: event labels at each timestamp rank joined by "+", ranks
# joined by ">" -- the identity seqpcplot.R:350-364 builds (order.align =
# "first" ranks distinct timestamps; simultaneous events share a rank).
key_from_points <- function(x, y, labels) {
  ranks <- sort(unique(x))
  paste(vapply(ranks, function(r) {
    paste(labels[sort(y[x == r])], collapse = "+")
  }, character(1)), collapse = ">")
}

hide <- "#CCCCCC"
pc_rows <- character()
for (spec in pc_filters) {
  inp <- pc_cases[[spec$case]]
  seqe <- seqecreate(id = inp$id, timestamp = inp$timestamp, event = inp$event)
  seqeweight(seqe) <- inp$weights
  filter <- if (is.na(spec$level)) {
    seqpcfilter(spec$method)
  } else {
    seqpcfilter(spec$method, level = spec$level)
  }
  res <- seqpcplot(seqe, filter = filter, plot = FALSE, hide.col = hide)
  labels <- res$ylevs
  pts <- res$pts
  traj_key <- character()
  traj_drawn <- logical()
  for (t in sort(unique(pts$traj))) {
    sub <- pts[pts$traj == t, ]
    traj_key[t] <- key_from_points(sub$x, sub$y, labels)
    traj_drawn[t] <- any(toupper(sub$col.1) != hide)
  }
  tse <- TraMineR:::seqe2TSE(seqe)
  seq_ids <- unique(tse$id)
  for (position in seq_along(seq_ids)) {
    rows <- tse[tse$id == seq_ids[position], ]
    ranks <- match(rows$timestamp, sort(unique(rows$timestamp)))
    codes <- match(as.character(rows$event), labels)
    key <- key_from_points(ranks, codes, labels)
    hit <- which(traj_key == key)
    stopifnot(length(hit) == 1)
    pc_rows <- c(pc_rows, paste(
      spec$case, spec$method,
      if (is.na(spec$level)) "none" else sprintf("%.17g", spec$level),
      position, key, if (traj_drawn[hit]) "TRUE" else "FALSE",
      sep = ","
    ))
  }
}
writeLines(c("case,method,level,sequence,trajectory,drawn", pc_rows),
  file.path(out_dir, "seqpcplot_filter_trajectories.csv"))

# --------------------------------------------------------------------------- #
# 2. format(x, digits = d)                                                    #
# --------------------------------------------------------------------------- #
values <- c(
  0, 1, 2, 15, 100, 999, 1000, 1200, 1234.5, 9996, 99999, 100000, 123456,
  1e6, 1234567, 1e15, 0.5, 0.1, 0.3, 2 / 3, 1 / 3, 0.123456, 0.0123456,
  0.00123456, 0.000123456, 1.23456e-5, 1e-10, 2.5, 12.345, 123.45, 12.34,
  0.567, 1.23, 3.14159265358979, 0.00174977, 1.5e10, 0.001, 0.0001, 1e-5,
  99.995, 9.9999, 0.99999, 5e-324, 1.7976931348623157e308, -0.5, -1234.5,
  -1e-5, -2 / 3, NaN, Inf, -Inf
)
fmt_rows <- character()
for (digits in c(3L, 5L)) {
  for (value in values) {
    fmt_rows <- c(fmt_rows, paste(
      sprintf("%.17g", value), digits, format(value, digits = digits), sep = ","
    ))
  }
}
writeLines(c(
  paste0("# getOption(\"digits\") = ", getOption("digits")),
  "value,digits,formatted", fmt_rows
), file.path(out_dir, "r_format_digits.csv"))

# --------------------------------------------------------------------------- #
# 3. seqmsplot info text                                                      #
# --------------------------------------------------------------------------- #
rows <- list(
  c("A", "A", "B", "B", "C", "C"),
  c("A", "B", "B", "C", "C", "C"),
  c("B", "B", NA, "C", "D", "D"),
  c("A", "A", "A", "A", "B", "B"),
  c("C", "C", "C", "D", "D", NA),
  c("A", NA, "B", "B", "C", "C"),
  c("D", "D", "D", "D", "D", "D"),
  c("A", "A", "B", "B", "C", "C")
)
modst_data <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
modst_cases <- list(
  weighted = c(1.0, 2.0, 0.5, 3.0, 1.5, 1.0, 0.25, 2.0),
  unweighted = NULL
)
ms_rows <- character()
for (name in names(modst_cases)) {
  s <- suppressMessages(seqdef(modst_data,
    weights = modst_cases[[name]], right = "DEL", alphabet = c("A", "B", "C", "D")
  ))
  m <- seqmodst(s)
  nbrep <- attr(m, "Occurrences")
  n <- attr(m, "nbseq")
  # plot.stslist.modst.R: ctfreq and txt, verbatim.
  ctfreq <- round((nbrep / n) * 100, 1)
  txt <- paste("Modal state sequence (", nbrep, " occurrences, freq=", ctfreq, "%)", sep = "")
  ms_rows <- c(ms_rows, paste(name, nbrep, sprintf("%.17g", n), txt, sep = "\t"))
}
writeLines(c("case\toccurrences\tnbseq\ttext", ms_rows),
  file.path(out_dir, "seqmsplot_info_text.tsv"))

# --------------------------------------------------------------------------- #
# 4. Group factors behind grouped plot panels                                 #
# --------------------------------------------------------------------------- #
group_cases <- list(
  numeric = c(2, 10, 2, 10, 1, 1),
  numeric_na = c(2, NA, 1, 2, NA, 1),
  character_na = c("a", NA, "b", "a", "b", NA),
  factor_declared = factor(c("b", "a", "b", NA, "a", "b"), levels = c("c", "b", "a"))
)
lv_rows <- character()
for (name in names(group_cases)) {
  g <- group_cases[[name]]
  values_text <- paste(ifelse(is.na(g), "NA", as.character(g)), collapse = "|")
  lv_rows <- c(lv_rows, paste(
    name, "factor", values_text, paste(levels(factor(g)), collapse = "|"),
    sep = "\t"
  ))
  lv_rows <- c(lv_rows, paste(
    name, "TraMineR:::group", values_text,
    paste(levels(TraMineR:::group(g)), collapse = "|"),
    sep = "\t"
  ))
}
writeLines(c("case\tbuilder\tvalues\tlevels", lv_rows),
  file.path(out_dir, "plot_group_levels.tsv"))

cat("AUDIT20260922_S6_GOLDENS_DONE\n")
