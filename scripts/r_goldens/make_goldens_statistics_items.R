stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
# Live-R goldens for batch 15's sequence-statistics items (rulings R15-a to R15-k): what TraMineR and TraMineRextras
# return for the inputs the tests use, including R's errors where traminex stops as R does or corrects R.
#
# Each row is case, field, value. Vectors are joined by spaces (numbers in %.17g, NA as NA); text vectors by ' ;; '.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_statistics_items.R
suppressMessages(library(TraMineR)); suppressMessages(library(TraMineRextras))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14", as.character(packageVersion("TraMineRextras")) == "0.6.9")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
rows <- list()
add <- function(case, field, value) rows[[length(rows) + 1L]] <<- data.frame(case = case, field = field, value = value)
num <- function(x) paste(ifelse(is.na(x), "NA", sprintf("%.17g", x)), collapse = " ")
txt <- function(x) paste(ifelse(is.na(x), "NA", x), collapse = " ;; ")
err <- function(expr) tryCatch({suppressWarnings(suppressMessages(expr)); "computed"},
                               error = function(e) paste("ERROR:", conditionMessage(e)))

# --- seqdomassoc's assoc as one string (R15-a: TX-AUD-20261010-062). seqdomassoc.Rd: assoc is a character string or
# a vector with both measures. Recorded: the two domains' cells and R's df, LRT and p for assoc = "LRT".
dom1 <- rbind(c("A","B","A","B"), c("B","B","A","A"), c("A","A","B","B"), c("B","A","B","A"), c("A","B","B","B"), c("B","B","B","A"))
dom2 <- rbind(c("X","Y","Y","X"), c("Y","Y","X","X"), c("X","X","Y","Y"), c("Y","X","X","Y"), c("X","Y","Y","Y"), c("Y","X","Y","X"))
add("domassoc", "domain_1", txt(apply(dom1, 1, paste, collapse = "-")))
add("domassoc", "domain_2", txt(apply(dom2, 1, paste, collapse = "-")))
dom_r <- suppressMessages(seqdomassoc(list(suppressMessages(seqdef(dom1)), suppressMessages(seqdef(dom2))), assoc = "LRT"))
add("domassoc", "columns", txt(colnames(dom_r)))
add("domassoc", "values", num(as.vector(unclass(dom_r))))

# --- seqindic.dyn on no sequences (R15-c: TX-AUD-20261010-096). max() of no lengths is -Inf, so the window count is
# NaN and matrix() stops (seqindic-dyn.R:28-72). Recorded: the object's cells and R's error on its first 0 rows.
dyn_cells <- rbind(c("A","B","A"), c("B","B","A"), c("A","A","B"), c("B","A","A"))
add("indic_dyn_empty", "cells", txt(apply(dyn_cells, 1, paste, collapse = "-")))
add("indic_dyn_empty", "result", err(seqindic.dyn(suppressMessages(seqdef(dyn_cells))[0, ])))

# --- seqgen's alphabet type and a one-value alphabet (R15-d: TX-AUD-20261011-001). sample() keeps the alphabet's type;
# one number x >= 1 is read as 1:x, so alphabet = 5 draws 1..5 although seqgen.Rd says the sequences are generated
# from the alphabet. Recorded: R's typeof for an integer, a double and a text alphabet, and the values drawn for 5.
add("seqgen_type", "typeof_integer", typeof({set.seed(1); seqgen(2, 3, 1:3)}))
add("seqgen_type", "typeof_double", typeof({set.seed(1); seqgen(2, 3, c(1, 2, 3))}))
add("seqgen_type", "typeof_character", typeof({set.seed(1); seqgen(2, 3, c("A", "B"))}))
add("seqgen_type", "one_value_5", num(sort(unique(as.vector({set.seed(1); seqgen(50, 20, 5)})))))

# --- seqgen's p as probability weights (R15-i: TX-AUD-20261011-072). seqgen passes p to sample(), which normalizes
# weights: under one seed c(5, 3, 2) draws what c(.5, .3, .2) draws. Recorded: that identity and sample()'s
# messages for a negative, a missing, an infinite and an all-zero weight vector.
add("seqgen_weights", "same_draws", as.character(identical({set.seed(1); seqgen(50, 10, c("A", "B", "C"), p = c(5, 3, 2))},
                                                           {set.seed(1); seqgen(50, 10, c("A", "B", "C"), p = c(.5, .3, .2))})))
weight_cases <- list(negative = c(-1, 1, 1), missing = c(NA, 1, 1), infinite = c(Inf, 1, 1), zero = c(0, 0, 0),
                     negative_then_missing = c(-1, NA, 1))
for (case in names(weight_cases)) {
  add("seqgen_weights", paste0("p_", case), num(weight_cases[[case]]))
  add("seqgen_weights", paste0("result_", case), err(seqgen(2, 3, c("A", "B", "C"), p = weight_cases[[case]])))
}

# --- seqpcfilter's method through match.arg (R15-e: TX-AUD-20261011-002). match.arg(method) takes an exact name or a
# unique prefix, case-sensitive (seqpcplot.R:1100); the level check uses the unmatched name (seqpcplot.R:1101), so
# ("cum", NULL) escapes it. Recorded: R's type, value and level, or R's error, per call.
pc_cases <- list(cum = list("cum", 0.5), m = list("m", 0.5), lin = list("lin", NULL), cumfreq = list("cumfreq", 0.5),
                 empty = list("", 0.5), x = list("x", 0.5), upper = list("CUM", 0.5), cum_null = list("cum", NULL))
for (case in names(pc_cases)) {
  spec <- pc_cases[[case]]
  add(paste0("pcfilter_", case), "method", spec[[1]])
  add(paste0("pcfilter_", case), "level", if (is.null(spec[[2]])) "NULL" else num(spec[[2]]))
  r <- tryCatch(unclass(seqpcfilter(spec[[1]], spec[[2]])), error = function(e) paste("ERROR:", conditionMessage(e)))
  add(paste0("pcfilter_", case), "result",
      if (is.character(r)) r else paste(r$type, r$value, if (is.null(r$level)) "NULL" else num(r$level)))
}

# --- The number of transitions (R15-f: TX-AUD-20261011-068). R's seqtransn leaves missing cells out of the spells and
# divides by each sequence's length (missing cells counted, void not) minus 1 (seqtransn.R:11-40). Recorded: R's
# seqtransn with norm FALSE and TRUE for the gap audit's ragged object, an object with an all-missing row (right = NA)
# and one with a length-1 sequence (right = "DEL").
transn_cases <- list(
  transn_ragged = list(cells = rbind(c("A","A","B","B","C","C"), c("A","B","B",NA,"C","D"), c("B","B","C",NA,NA,NA),
                                     c("C","C","C","D","D","D"), c(NA,"A","A","B","B","B"), c("D","D","C","C","B","A"),
                                     c("A","A","A","A",NA,NA), c("B",NA,NA,"C","C","C")), right = "DEL"),
  transn_missing_row = list(cells = rbind(c("A","B","A"), c(NA,NA,NA), c("B",NA,NA)), right = NA),
  transn_length_one = list(cells = rbind(c("A","B","A"), c("B",NA,NA)), right = "DEL"))
for (case in names(transn_cases)) {
  spec <- transn_cases[[case]]
  s <- suppressMessages(seqdef(spec$cells, alphabet = sort(unique(na.omit(as.vector(spec$cells)))), right = spec$right))
  add(case, "cells", txt(apply(spec$cells, 1, function(r) paste(ifelse(is.na(r), "NA", r), collapse = " "))))
  add(case, "right", if (is.na(spec$right)) "NA" else spec$right)
  add(case, "seqconc", txt(as.vector(seqconc(s))))
  add(case, "transn", num(as.vector(seqtransn(s))))
  add(case, "transn_norm", num(as.vector(seqtransn(s, norm = TRUE))))
}

# --- seqsha's time and event (R15-g: TX-AUD-20261011-069). R stops on a missing time or event (seqsha.R:6-8) and on a
# time of 0 or below (rep() and [<- errors); seqsha.Rd: event is logical, and R keeps text events as a text column.
# Recorded: the object, times and events, and per case R's rows and event values, or R's error.
sha_cells <- rbind(c("A","A","B","B","C"), c("A","B","B","C","C"), c("B","B","B","C","C"), c("A","A","A","A","B"))
sha_seq <- suppressMessages(seqdef(sha_cells))
sha_time <- c(3, 2, 4, 5); sha_event <- c(TRUE, FALSE, TRUE, FALSE)
add("sha", "cells", txt(apply(sha_cells, 1, paste, collapse = "-")))
add("sha", "time", num(sha_time))
add("sha", "event", txt(as.character(sha_event)))
sha_cases <- list(sha_control = list(sha_time, sha_event), sha_time_zero = list(replace(sha_time, 1, 0), sha_event),
                  sha_time_negative = list(replace(sha_time, 1, -1), sha_event),
                  sha_time_missing = list(replace(sha_time, 2, NA), sha_event),
                  sha_event_missing = list(sha_time, replace(sha_event, 2, NA)),
                  sha_event_text = list(sha_time, ifelse(sha_event, "yes", "no")),
                  sha_event_numeric = list(sha_time, as.numeric(sha_event)),
                  sha_time_fraction = list(replace(sha_time, 1, 2.5), sha_event),
                  sha_time_infinite = list(replace(sha_time, 1, Inf), sha_event))
for (case in names(sha_cases)) {
  spec <- sha_cases[[case]]
  add(case, "time", num(spec[[1]]))
  add(case, "event", txt(as.character(spec[[2]])))
  r <- tryCatch(suppressMessages(seqsha(sha_seq, spec[[1]], spec[[2]])), error = function(e) paste("ERROR:", conditionMessage(e)))
  add(case, "result", if (is.character(r)) r else paste("rows", nrow(r), "| event", class(r$event), "|", txt(as.character(r$event))))
}

# --- Degenerate inputs where R stops (R15-h: TX-AUD-20261011-070). On one-state data seqidegrad, seqipos and
# seqprecarity(with.missing = TRUE) stop in their own indexing, and seqlogp(with.missing = TRUE) stops on data without
# missing cells, though R's definitions give a value. Recorded: the objects, R's four errors and seqlogp's
# with.missing = FALSE values.
one_cells <- rbind(c("A", "A"), c("A", "A"))
two_cells <- rbind(c("A", "B", "A"), c("B", "B", "A"), c("A", "A", "B"))
one_state <- suppressMessages(seqdef(one_cells, alphabet = "A"))
two_state <- suppressMessages(seqdef(two_cells, alphabet = c("A", "B")))
add("degenerate", "one_state", txt(apply(one_cells, 1, paste, collapse = "-")))
add("degenerate", "two_state", txt(apply(two_cells, 1, paste, collapse = "-")))
add("degenerate", "seqidegrad", err(seqidegrad(one_state)))
add("degenerate", "seqipos", err(seqipos(one_state, pos.states = "A")))
add("degenerate", "seqprecarity_with_missing", err(seqprecarity(one_state, with.missing = TRUE)))
add("degenerate", "seqlogp_with_missing", err(seqlogp(two_state, prob = "freq", begin = "global.freq", with.missing = TRUE)))
add("degenerate", "seqlogp_without_missing",
    num(as.vector(suppressMessages(seqlogp(two_state, prob = "freq", begin = "global.freq")))))

# --- seqimplic's group length (R15-j: TX-AUD-20261011-073). seqimplic.Rd: group gives the membership of each
# sequence. R builds its row mask from group (seqimplic.R:69-72): a shorter mask is recycled over the rows, a longer
# one indexes rows that do not exist. Recorded: the object, the groups, and per case R's error or the dimensions, the
# finite cells and the cells of R's result; with the full group, R's result for comparison.
imp_cells <- rbind(c("A","A","B","B"), c("A","B","B","B"), c("B","B","A","A"), c("A","A","A","B"), c("B","A","A","A"),
                   c("B","B","B","A"))
imp_seq <- suppressMessages(seqdef(imp_cells))
imp_group <- c("g", "h", "g", "h", "g", "h")
add("implic", "cells", txt(apply(imp_cells, 1, paste, collapse = "-")))
add("implic", "group", txt(imp_group))
imp_cases <- list(implic_full = imp_group, implic_short = imp_group[1:5], implic_long = c(imp_group, "g"))
for (case in names(imp_cases)) {
  r <- tryCatch(suppressWarnings(suppressMessages(seqimplic(imp_seq, group = imp_cases[[case]]))),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  add(case, "group_length", num(length(imp_cases[[case]])))
  add(case, "result", if (is.character(r)) r else paste("dim", num(dim(r$indices)), "| finite", num(sum(is.finite(r$indices)))))
  if (!is.character(r)) add(case, "indices", num(as.vector(r$indices)))
}

# --- seqdiff and dissassoc on degenerate inputs (batch-15 review I-1: TX-AUD-20261010-096). A window without a
# complete case gives seqdist no sequence, and its matrix() stops; with all-zero weights dissassoc divides 0 by 0 and
# returns NaN. Recorded: R's errors for no rows, an all-void row and all-missing data, and R's seqdiff and dissassoc
# values (with the LCS distances they read) for all-zero weights.
zw_cells <- rbind(c("A","B"), c("B","A"))
zw_seq <- suppressMessages(seqdef(zw_cells, weights = c(0, 0)))
add("diss_zero_weights", "cells", txt(apply(zw_cells, 1, paste, collapse = "-")))
add("diss_zero_weights", "group", txt(c("g", "h")))
zw_diff <- suppressWarnings(suppressMessages(seqdiff(zw_seq, group = c("g", "h"))))
add("diss_zero_weights", "seqdiff_stat", num(as.vector(zw_diff$stat)))
add("diss_zero_weights", "seqdiff_discrepancy", num(as.vector(zw_diff$discrepancy)))
zw_d <- suppressMessages(seqdist(zw_seq, method = "LCS", norm = "auto"))
add("diss_zero_weights", "diss", num(as.vector(zw_d)))
zw_a <- suppressWarnings(dissassoc(zw_d, c("g", "h"), R = 0, weights = c(0, 0), weight.permutation = "diss"))
add("diss_zero_weights", "groups_n", num(zw_a$groups$n))
add("diss_zero_weights", "groups_discrepancy", num(zw_a$groups$discrepancy))
add("diss_zero_weights", "stat", num(zw_a$stat[, 1]))
add("diss_zero_weights", "anova_ss", num(zw_a$anova.table$SS))
add("diss_zero_weights", "anova_df", num(zw_a$anova.table$df))
add("diss_zero_weights", "anova_mse", num(zw_a$anova.table$MSE))
nc_base <- suppressMessages(seqdef(rbind(c("A","B","A"), c("B","B","A"))))
nc_cases <- list(no_rows = list(seq = nc_base[0, ], group = character(0)),
                 all_void_row = list(seq = suppressMessages(seqdef(rbind(c("A", NA), c(NA, NA)))), group = c("g", "h")),
                 all_missing = list(seq = suppressMessages(seqdef(rbind(c(NA, NA), c(NA, NA)), right = NA, alphabet = "A")),
                                    group = c("g", "h")))
for (case in names(nc_cases)) {
  add(paste0("seqdiff_", case), "result", err(seqdiff(nc_cases[[case]]$seq, group = nc_cases[[case]]$group)))
}

# --- seqemlt without an observed state (batch-15 review I-1: TX-REP-20261023-005). With no sequence, or only missing
# cells, no situation is observed and R stops ("subscript out of bounds"). Recorded: R's error for each.
emlt_base <- suppressMessages(seqdef(rbind(c("A","B","A"), c("B","B","A"))))
add("seqemlt_no_rows", "result", err(seqemlt(emlt_base[0, ])))
add("seqemlt_all_missing", "result",
    err(seqemlt(suppressMessages(seqdef(rbind(c(NA, NA), c(NA, NA)), right = NA, alphabet = "A")))))

# --- The same one-state data with a two-state alphabet (batch-15 review M-5: TX-AUD-20261011-070). R computes
# seqidegrad and seqipos when the alphabet has a second, unobserved state; its seqprecarity stops for both
# with.missing values. Recorded: R's values or errors.
two_alph <- suppressMessages(seqdef(one_cells, alphabet = c("A", "B")))
add("degenerate_two_alphabet", "seqidegrad", num(as.vector(suppressMessages(seqidegrad(two_alph)))))
add("degenerate_two_alphabet", "seqipos", num(as.vector(suppressMessages(seqipos(two_alph, pos.states = "A")))))
add("degenerate_two_alphabet", "seqprecarity", err(seqprecarity(two_alph)))
add("degenerate_two_alphabet", "seqprecarity_with_missing", err(seqprecarity(two_alph, with.missing = TRUE)))

# --- plot.stslist's sample.meth (batch-15 review M-7: TX-REP-20261023-007). plot.stslist.R:16 reads
# match.arg(sample.meth, c("systematic", "random")): an exact name or a unique prefix, case-sensitive. Recorded: R's
# result or error per value.
for (value in c("systematic", "sys", "r", "", "SYS", "x")) {
  add(paste0("sample_meth_", if (nzchar(value)) value else "empty"), "result",
      tryCatch(match.arg(value, c("systematic", "random")), error = function(e) paste("ERROR:", conditionMessage(e))))
}

write.csv(do.call(rbind, rows), file.path(out_dir, "statistics_items.csv"), row.names = FALSE, quote = TRUE)
