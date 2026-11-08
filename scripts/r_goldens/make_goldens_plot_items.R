stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
# Live-R goldens for batch 13's plotting items (rulings R13-a to R13-p): the arguments TraMineR's and TraMineRextras'
# plotting code passes to the graphics functions, traced inside the packages.
#
# Each row is case, field, value. Vectors are joined by spaces (numbers in %.17g, NA as NA); text vectors by ' ;; '.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_plot_items.R
suppressMessages(library(TraMineR)); suppressMessages(library(TraMineRextras))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14", as.character(packageVersion("TraMineRextras")) == "0.6.9")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
rows <- list()
add <- function(case, field, value) rows[[length(rows) + 1L]] <<- data.frame(case = case, field = field, value = value)
num <- function(x) paste(ifelse(is.na(x), "NA", sprintf("%.17g", x)), collapse = " ")
txt <- function(x) paste(x, collapse = " ;; ")
pdf(NULL)

# --- Weighted index plots (R13-a: 048). plot.stslist.R:149 draws barplot(width = weights[idxs[sel]], horiz = TRUE);
# seqIplot passes space = 0, so the rows touch. Recorded: the widths and the bar centres barplot returns, and the y
# label positions given to axis(side = 2), for the rows as drawn (sorted when sortv is given). All-1 weights and
# weighted = FALSE draw equal rows (plot.stslist.R:103-104).
bars <- list(); ylabs <- list()
trace(graphics::barplot, print = FALSE, where = asNamespace("TraMineR"),
      tracer = quote(bars$width <<- list(...)$width), exit = quote(bars$mid <<- as.vector(returnValue())))
trace(graphics::axis, print = FALSE, where = asNamespace("TraMineR"),
      tracer = quote(if (side == 2) ylabs[[length(ylabs) + 1L]] <<- at))
index_cells <- rbind(c("A", "A", "B"), c("A", "B", "B"), c("B", "B", "A"), c("B", "A", "A"), c("A", "A", "A"), c("B", "B", "B"))
index_weights <- c(1, 3, 0.5, 2, 1, 1.5)
index_sortv <- c(3, 1, 6, 2, 5, 4)
add("index", "cells", txt(apply(index_cells, 1, paste, collapse = "")))
add("index", "weights", num(index_weights))
add("index", "sortv", num(index_sortv))
cases <- list(index_weighted = list(w = index_weights, sortv = NULL, weighted = TRUE),
              index_weighted_sorted = list(w = index_weights, sortv = index_sortv, weighted = TRUE),
              index_unweighted = list(w = index_weights, sortv = NULL, weighted = FALSE),
              index_all_one = list(w = rep(1, 6), sortv = NULL, weighted = TRUE))
for (case in names(cases)) {
  spec <- cases[[case]]
  s <- suppressMessages(seqdef(index_cells, alphabet = c("A", "B"), weights = spec$w))
  bars <- list(); ylabs <- list()
  if (is.null(spec$sortv)) seqIplot(s, weighted = spec$weighted, with.legend = FALSE)
  else seqIplot(s, weighted = spec$weighted, sortv = spec$sortv, with.legend = FALSE)
  add(case, "widths", num(bars$width))
  add(case, "centres", num(bars$mid))
  add(case, "ylabel_positions", num(ylabs[[1]]))
}
untrace(graphics::barplot, where = asNamespace("TraMineR"))
untrace(graphics::axis, where = asNamespace("TraMineR"))

# --- One frequency scale in plot.subseqelistchisq (R13-b: 049): ylim = c(min(min(freq), 0), max(freq)) over every
# group's frequency column (seqecmpgroup.R:136-138), traced per panel. 40 event sequences: in g1 most have a then b,
# in g2 most have c, and only the odd ones then b, so the two groups' largest frequencies differ; every third has d.
ids <- integer(0); times <- numeric(0); events <- character(0)
for (i in 1:40) {
  first <- if ((i <= 20 && i %% 10 != 0) || (i > 20 && i %% 10 == 0)) "a" else "c"
  ids <- c(ids, i); times <- c(times, 1); events <- c(events, first)
  if (i <= 20 || i %% 2 == 1) { ids <- c(ids, i); times <- c(times, 2); events <- c(events, "b") }
  if (i %% 3 == 0) { ids <- c(ids, i); times <- c(times, 3); events <- c(events, "d") }
}
chisq_group <- rep(c("g1", "g2"), each = 20)
add("chisq", "tse_id", num(ids))
add("chisq", "tse_time", num(times))
add("chisq", "tse_event", txt(events))
add("chisq", "group", txt(chisq_group))
chisq_fsub <- seqefsub(seqecreate(id = ids, timestamp = times, event = events), pmin.support = 0.05)
chisq_cmp <- seqecmpgroup(chisq_fsub, group = chisq_group)
panel_ylim <- list()
trace(graphics::barplot, print = FALSE, where = asNamespace("TraMineR"),
      tracer = quote(panel_ylim[[length(panel_ylim) + 1L]] <<- list(...)$ylim))
plot(chisq_cmp)
untrace(graphics::barplot, where = asNamespace("TraMineR"))
add("chisq", "labels", txt(chisq_cmp$labels))
add("chisq", "max_freq_per_group", num(apply(chisq_cmp$data[, 5:6], 2, max)))
for (k in seq_along(panel_ylim)) add("chisq", paste0("ylim_panel", k), num(panel_ylim[[k]]))

# --- seqedplot with a palette shorter than the events (R13-c: 050). R colours event k with cpal[k]
# (seqedplot.R:192-227), so with cpal = "red" and events A, B the second curve gets colour NA (drawn in no colour,
# without an error), in every group panel; the legend receives both event labels and the one colour. Recorded: the
# col each survival curve receives (survival's plot.survfit for the first event of a panel, lines.survfit for the
# others), panel after panel, and the text and colors given to TraMineR.legend.
seqed_ids <- c(1, 1, 2, 3, 3, 4, 4, 5)
seqed_times <- c(1, 3, 2, 1, 4, 1, 3, 2)
seqed_events <- c("A", "B", "A", "B", "A", "A", "A", "B")
seqed_lengths <- c(5, 5, 6, 4, 5)
seqed_group <- c("g1", "g1", "g1", "g2", "g2")
add("seqed", "tse_id", num(seqed_ids))
add("seqed", "tse_time", num(seqed_times))
add("seqed", "tse_event", txt(seqed_events))
add("seqed", "lengths", num(seqed_lengths))
add("seqed", "group", txt(seqed_group))
seqed_e <- seqecreate(id = seqed_ids, timestamp = seqed_times, event = seqed_events)
seqelength(seqed_e) <- seqed_lengths
curve_cols <- character(0); legend_args <- list()
trace(survival:::plot.survfit, print = FALSE, where = asNamespace("survival"),
      tracer = quote(curve_cols <<- c(curve_cols, as.character(col))))
trace(survival:::lines.survfit, print = FALSE, where = asNamespace("survival"),
      tracer = quote(curve_cols <<- c(curve_cols, as.character(col))))
trace("TraMineR.legend", print = FALSE, where = asNamespace("TraMineR"),
      tracer = quote(legend_args <<- list(text = text, colors = colors)))
for (case in c("seqed_short", "seqed_short_grouped")) {
  curve_cols <- character(0); legend_args <- list()
  if (case == "seqed_short") seqedplot(seqed_e, cpal = "red") else seqedplot(seqed_e, cpal = "red", group = seqed_group)
  add(case, "cpal", txt("red"))
  add(case, "curve_colours", txt(curve_cols))
  add(case, "legend_text", txt(legend_args$text))
  add(case, "legend_colours", txt(legend_args$colors))
}
untrace(survival:::plot.survfit, where = asNamespace("survival"))
untrace(survival:::lines.survfit, where = asNamespace("survival"))
untrace("TraMineR.legend", where = asNamespace("TraMineR"))

# --- The frequency plot's y axis (R13-j: 11-067 a). plot.stslist.freq.R:56-61, 102-110: two labels, 0% at the
# bottom of the first bar and the cumulative percentage of the shown sequences at the top of the last bar, titled
# "Cum. % freq. (n=N)" ("weighted n" when weighted, n rounded to 2 decimals). Recorded: the axis(side = 2) positions
# and labels, the barplot ylab, and the bars' widths and centres (barplot's value), so the label positions can be
# compared with the bars' edges.
freq_cells <- rbind(c("A", "A", "B", "B"), c("A", "B", "B", "C"), c("A", "A", "B", "B"), c("C", "C", "A", "B"),
                    c("A", "A", "B", "B"), c("A", "B", "B", "C"))
freq_weights <- c(1, 2, 0.5, 1.5, 1, 0.5)
add("freq", "cells", txt(apply(freq_cells, 1, paste, collapse = "")))
add("freq", "weights", num(freq_weights))
freq_bars <- list(); freq_axes <- list()
trace(graphics::barplot, print = FALSE, where = asNamespace("TraMineR"),
      tracer = quote(freq_bars <<- list(ylab = list(...)$ylab, width = list(...)$width)),
      exit = quote(freq_bars$mid <<- as.vector(returnValue())))
trace(graphics::axis, print = FALSE, where = asNamespace("TraMineR"),
      tracer = quote(if (side == 2) freq_axes <<- list(at = at, labels = labels)))
freq_cases <- list(freq_unweighted = list(w = NULL, idxs = 1:10), freq_weighted = list(w = freq_weights, idxs = 1:10),
                   freq_top2 = list(w = NULL, idxs = 1:2))
for (case in names(freq_cases)) {
  spec <- freq_cases[[case]]
  s <- suppressMessages(seqdef(freq_cells, alphabet = c("A", "B", "C"), weights = spec$w))
  freq_bars <- list(); freq_axes <- list()
  seqfplot(s, idxs = spec$idxs, with.legend = FALSE)
  add(case, "idxs", num(spec$idxs))
  add(case, "ylab", freq_bars$ylab)
  add(case, "ylabels", txt(freq_axes$labels))
  add(case, "ylabel_positions", num(freq_axes$at))
  add(case, "bar_widths", num(freq_bars$width))
  add(case, "bar_centres", num(freq_bars$mid))
}
untrace(graphics::barplot, where = asNamespace("TraMineR"))
untrace(graphics::axis, where = asNamespace("TraMineR"))

# --- Time ticks on R's xtstep grid (R13-k: 11-067 b). plot.stslist.R:160-161: ticks at seq(1, seql, by = xtstep),
# plus seql when tick.last, each given its column name (the device then drops overlapping labels when drawing).
# Recorded for seqIplot: the axis(side = 1) positions (tpos - 0.5) and labels. The cells do not matter; the columns
# are named c1, c2, ...
tick_axes <- list()
trace(graphics::axis, print = FALSE, where = asNamespace("TraMineR"),
      tracer = quote(if (side == 1) tick_axes <<- list(at = at, labels = labels)))
tick_cases <- list(ticks_step6_72 = list(L = 72, xtstep = 6, tick.last = FALSE),
                   ticks_default_72 = list(L = 72, xtstep = 1, tick.last = FALSE),
                   ticks_step5_last_72 = list(L = 72, xtstep = 5, tick.last = TRUE),
                   ticks_step3_last_5 = list(L = 5, xtstep = 3, tick.last = TRUE))
for (case in names(tick_cases)) {
  spec <- tick_cases[[case]]
  cells <- matrix(rep(c("A", "B"), length.out = 3 * spec$L), nrow = 3)
  s <- suppressMessages(seqdef(cells, alphabet = c("A", "B"), cnames = paste0("c", seq_len(spec$L)),
                               xtstep = spec$xtstep, tick.last = spec$tick.last))
  tick_axes <- list()
  seqIplot(s, with.legend = FALSE)
  add(case, "length", num(spec$L))
  add(case, "xtstep", num(spec$xtstep))
  add(case, "tick_last", txt(spec$tick.last))
  add(case, "positions", num(tick_axes$at))
  add(case, "labels", txt(tick_axes$labels))
}
untrace(graphics::axis, where = asNamespace("TraMineR"))

# --- Ticks at the cell centres (R13-l: 11-067 h). The frequency, modal-state, representative and
# relative-frequency plots put the tick of position t at t - 0.5, the centre of a cell drawn on [t - 1, t]
# (plot.stslist.freq.R:91-92, plot.stslist.modst.R:94, plot.stslist.rep.R:144, plot.seqrf.R through plot.stslist).
# Recorded: the first axis(side = 1) positions of each, on the frequency rows' sequences (the rf plot on those rows
# repeated to 20 sequences, so its default k = floor(20 / 10) = 2; the representative plot computes its Hamming
# distances itself, seqrep's default diss = NULL).
centre_at <- NULL
trace(graphics::axis, print = FALSE, where = asNamespace("TraMineR"),
      tracer = quote(if (side == 1 && is.null(centre_at)) centre_at <<- at))
centre_s <- suppressMessages(seqdef(freq_cells, alphabet = c("A", "B", "C")))
centre_rf_rows <- rep(seq_len(nrow(freq_cells)), length.out = 20)
centre_s20 <- suppressMessages(seqdef(freq_cells[centre_rf_rows, ], alphabet = c("A", "B", "C")))
centre_d20 <- suppressMessages(seqdist(centre_s20, method = "HAM"))
add("centres", "rf_rows", num(centre_rf_rows))
centre_calls <- list(centres_seqfplot = function() seqfplot(centre_s, with.legend = FALSE),
                     centres_seqmsplot = function() seqmsplot(centre_s, with.legend = FALSE),
                     centres_seqrplot = function() seqrplot(centre_s, method = "HAM", with.legend = FALSE),
                     centres_seqrfplot = function() seqrfplot(centre_s20, diss = centre_d20, with.legend = FALSE))
for (case in names(centre_calls)) {
  centre_at <- NULL
  suppressMessages(centre_calls[[case]]())
  add(case, "positions", num(centre_at))
}
untrace(graphics::axis, where = asNamespace("TraMineR"))

# --- seqplot's common options (R13-m: 11-067 c, d). R's seqplot takes ylab, with.legend and cpal for every type
# (seqplot.R:6-11); it drops cpal for Ht (seqplot.R:299-301) and draws no legend for Ht and pc (seqplot.R:210).
# Recorded per type, with ylab = "my y", with.legend = FALSE and cpal = red, green, blue (not for Ht): the y label R
# draws (barplot's ylab, else title()'s for the representative plot, else plot()'s for Ht), the colours barplot
# receives, and the number of state legends (TraMineR.legend calls), with with.legend = FALSE and at its default.
opt_trace <- list(); opt_legends <- 0L
trace(graphics::barplot, print = FALSE, where = asNamespace("TraMineR"),
      tracer = quote(opt_trace[[length(opt_trace) + 1L]] <<- list(ylab = list(...)$ylab, col = list(...)$col)))
trace(graphics::title, print = FALSE, where = asNamespace("TraMineR"),
      tracer = quote(if (!is.null(ylab)) opt_trace[[length(opt_trace) + 1L]] <<- list(ylab = ylab, col = NULL)))
trace(graphics::plot.default, print = FALSE,
      tracer = quote(opt_trace[[length(opt_trace) + 1L]] <<- list(ylab = ylab, col = NULL)))
trace("TraMineR.legend", print = FALSE, where = asNamespace("TraMineR"), tracer = quote(opt_legends <<- opt_legends + 1L))
opt_cpal <- c("red", "green", "blue")
add("options", "cpal", txt(opt_cpal))
for (typ in c("i", "I", "d", "dH", "f", "ms", "mt", "r", "rf", "Ht")) {
  obj <- if (typ == "rf") centre_s20 else centre_s
  extra <- if (typ == "r") list(method = "HAM") else if (typ == "rf") list(diss = centre_d20) else list()
  args <- c(list(obj, type = typ, ylab = "my y", with.legend = FALSE), extra)
  if (typ != "Ht") args$cpal <- opt_cpal
  opt_trace <- list(); opt_legends <- 0L
  suppressMessages(do.call(seqplot, args))
  ylabs <- Filter(function(v) !is.null(v) && !is.na(v), lapply(opt_trace, `[[`, "ylab"))
  cols <- Filter(Negate(is.null), lapply(opt_trace, `[[`, "col"))
  case <- paste0("options_", typ)
  add(case, "ylab", ylabs[[1]])
  add(case, "colours", if (length(cols)) txt(cols[[1]]) else "")
  add(case, "state_legends", num(opt_legends))
  opt_legends <- 0L
  suppressMessages(do.call(seqplot, c(list(obj, type = typ), extra)))
  add(case, "state_legends_default", num(opt_legends))
}
untrace(graphics::barplot, where = asNamespace("TraMineR"))
untrace(graphics::title, where = asNamespace("TraMineR"))
untrace(graphics::plot.default)
untrace("TraMineR.legend", where = asNamespace("TraMineR"))

# --- Missing cells in a solid colour (R13-o: 11-067 f). The index and distribution plots draw the missing state in
# missing.color, opaque. Recorded: the colours barplot receives (the last is the missing state's) and the alpha of
# missing.color (col2rgb, 0-255).
missing_bars <- NULL
trace(graphics::barplot, print = FALSE, where = asNamespace("TraMineR"),
      tracer = quote(if (is.null(missing_bars)) missing_bars <<- list(...)$col))
missing_cells <- rbind(c("A", "A", NA, "B"), c("A", "B", "B", "B"))
missing_s <- suppressMessages(seqdef(missing_cells, alphabet = c("A", "B"), right = NA))
add("missing", "cells", txt(apply(ifelse(is.na(missing_cells), "NA", missing_cells), 1, paste, collapse = " ")))
add("missing", "missing_color", attr(missing_s, "missing.color"))
add("missing", "missing_color_alpha", num(col2rgb(attr(missing_s, "missing.color"), alpha = TRUE)[4]))
missing_bars <- NULL
seqIplot(missing_s, with.legend = FALSE)
add("missing_index", "colours", txt(missing_bars))
missing_bars <- NULL
seqdplot(missing_s, with.legend = FALSE, with.missing = TRUE)
add("missing_distribution", "colours", txt(missing_bars))
untrace(graphics::barplot, where = asNamespace("TraMineR"))

# --- seqtrate on the transition-heatmap test's input (R13-p: TX-REP-20261019-003). Four sequences, one with a
# missing cell inside and one ending in two void cells, weights 1, 2, 1, 0.5; transition rates weighted and not, and
# the (weighted) counts, row by row over R's states.
heat_cells <- rbind(c("A", "A", "B", "B", "C"), c("A", "B", NA, "B", "C"), c("C", "C", "B", "A", "A"),
                    c("B", "B", "A", NA, NA))
heat_weights <- c(1, 2, 1, 0.5)
heat_s <- suppressMessages(seqdef(heat_cells, weights = heat_weights))
add("heatmap", "cells", txt(apply(ifelse(is.na(heat_cells), "NA", heat_cells), 1, paste, collapse = " ")))
add("heatmap", "weights", num(heat_weights))
add("heatmap", "states", txt(alphabet(heat_s)))
add("heatmap", "rates_weighted", num(t(suppressMessages(seqtrate(heat_s, weighted = TRUE)))))
add("heatmap", "rates_unweighted", num(t(suppressMessages(seqtrate(heat_s, weighted = FALSE)))))
add("heatmap", "counts_weighted", num(t(suppressMessages(seqtrate(heat_s, count = TRUE)))))

# --- The representative plot's y axis (R13-r: TX-REP-20261023-001). plot.stslist.rep.R:109-110, 132-136: the
# default y label "N representative(s) (n=N)", "weighted n=N" when weighted, N rounded to 2 decimals, drawn by title();
# the bands get no y labels (the percentage axis, 156-162, is commented out): the one side-2 axis labels the
# statistics overlay's rows B and A (232-236). Recorded on the frequency rows' sequences, unweighted and weighted:
# the title() ylab and the labels of the side-2 axis calls.
rep_ylab <- NULL; rep_side2 <- character(0)
trace(graphics::title, print = FALSE, where = asNamespace("TraMineR"),
      tracer = quote(if (!is.null(ylab)) rep_ylab <<- ylab))
trace(graphics::axis, print = FALSE, where = asNamespace("TraMineR"),
      tracer = quote(if (side == 2) rep_side2 <<- c(rep_side2, as.character(labels))))
for (case in c("rep_unweighted", "rep_weighted")) {
  w <- if (case == "rep_weighted") freq_weights else NULL
  s <- suppressMessages(seqdef(freq_cells, alphabet = c("A", "B", "C"), weights = w))
  rep_ylab <- NULL; rep_side2 <- character(0)
  suppressMessages(seqrplot(s, method = "HAM", with.legend = FALSE))
  add(case, "ylab", rep_ylab)
  add(case, "side2_labels", txt(rep_side2))
}
untrace(graphics::title, where = asNamespace("TraMineR"))
untrace(graphics::axis, where = asNamespace("TraMineR"))

invisible(dev.off())
write.csv(do.call(rbind, rows), file.path(out_dir, "plot_items.csv"), row.names = FALSE, quote = TRUE)
