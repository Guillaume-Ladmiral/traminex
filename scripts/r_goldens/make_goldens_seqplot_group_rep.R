# Live-R goldens for seqplot's grouped "r" and "rf" panels (repair of audit
# finding TX-AUD-20261010-044, batch-4 review I-1, I-2 and M-5).
#
# R's seqplot prepares types "r" and "rf" on the whole set before it splits it
# by group. For "r" it computes the distances from `method` when diss is not
# given and passes dmax = max(diss) of the whole matrix to every panel
# (seqplot.R:96-124); for "rf" it computes the default sort variable, the
# one-dimensional cmdscale (wcmdscale with weights) of the whole matrix
# (seqplot.R:127-156). Each panel then receives its subset (seqplot.R:331-339,
# :368). seqplot also stops for type "r" when a group has fewer than two cases,
# the whole set counting as one group without group (seqplot.R:97-101), and
# for type "rf" with group and which.plot = "both" (seqplot.R:53-56). The
# panel values below are traced from seqplot's own calls to seqrf and seqrep.
#
# A one-dimensional MDS axis has no sign of its own: R's cmdscale picks one and
# traminex may pick the other (TX-REP-20261017-006), and the sign changes the
# groups. "medoids" and "heights" are R's seqplot with its own axis;
# "medoids_negated" and "heights_negated" are R's seqplot given the negated axis
# as sortv. Copies of one sequence have one coordinate up to rounding, so their
# order, and which copy is a medoid, is not part of the definition. "subset_*"
# rows are what each panel gives when it computes those values on its own
# subset (R's seqrf and seqrep called per subset); they show that the fixture
# tells the two rules apart.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_seqplot_group_rep.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
icuSetCollate(locale = "root")
pdf(NULL)

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

data(biofam)
b <- biofam[1:60, ]
g <- as.character(b$sex)
levs <- levels(factor(g))
stopifnot(identical(trimws(levels(TraMineR:::group(g))), levs))
s <- suppressMessages(seqdef(b[, 10:25], alphabet = 0:7))
sw <- suppressMessages(seqdef(b[, 10:25], alphabet = 0:7, weights = b$wp00tbgs))
D <- suppressMessages(seqdist(s, method = "LCS"))

rows <- list()
add <- function(case, level, quantity, values) {
  rows[[length(rows) + 1L]] <<- data.frame(case = case, level = level, quantity = quantity,
                                           position = seq_along(values), value = as.character(values),
                                           stringsAsFactors = FALSE)
}
rep_text <- function(rep) apply(as.matrix(rep), 1, paste, collapse = "-")

# --- type "rf": seqrf's sortv, medoids and group heights per panel -----------
traced_sortv <- list()
traced_medoids <- list()
traced_heights <- list()
trace("seqrf", tracer = quote(traced_sortv <<- c(traced_sortv, list(sortv))),
      exit = quote({
        traced_medoids <<- c(traced_medoids, list(returnValue()$rf$medoids))
        traced_heights <<- c(traced_heights, list(returnValue()$rf$heights))
      }),
      print = FALSE, where = asNamespace("TraMineR"))
grouped_rf <- function(seqs, ...) {
  traced_sortv <<- list()
  traced_medoids <<- list()
  traced_heights <<- list()
  suppressMessages(seqplot(seqs, type = "rf", group = g, diss = D, k = 4, ...))
  stopifnot(length(traced_medoids) == length(levs))
  list(sortv = traced_sortv, medoids = traced_medoids, heights = traced_heights)
}
rf_cases <- list(
  list(case = "rf_k4", seqs = s, axis = cmdscale(D, k = 1)),
  list(case = "rf_weighted_k4", seqs = sw, axis = vegan::wcmdscale(D, k = 1, w = attr(sw, "weights")))
)
for (cs in rf_cases) {
  own <- grouped_rf(cs$seqs)
  negated <- grouped_rf(cs$seqs, sortv = -cs$axis)
  add(cs$case, "all", "global_axis", format(cs$axis[, 1], digits = 17))
  told_apart <- FALSE
  for (np in seq_along(levs)) {
    i <- which(g == levs[np])
    stopifnot(isTRUE(all.equal(as.numeric(own$sortv[[np]]), as.numeric(cs$axis[i, 1]))))
    subset <- suppressMessages(seqrf(cs$seqs[i, ], diss = D[i, i], k = 4))$rf$medoids
    told_apart <- told_apart || (!identical(subset, own$medoids[[np]]) && !identical(subset, negated$medoids[[np]]))
    add(cs$case, levs[np], "medoids", own$medoids[[np]])
    add(cs$case, levs[np], "medoids_negated", negated$medoids[[np]])
    add(cs$case, levs[np], "heights", format(as.numeric(own$heights[[np]]), digits = 17))
    add(cs$case, levs[np], "heights_negated", format(as.numeric(negated$heights[[np]]), digits = 17))
    add(cs$case, levs[np], "subset_medoids", subset)
  }
  stopifnot(told_apart)
}
untrace("seqrf", where = asNamespace("TraMineR"))

# --- type "r": seqrep's distances, dmax and representatives per panel ----------
traced_diss_sum <- list()
traced_dmax <- list()
traced_reps <- list()
trace_seqrep <- function() {
  trace("seqrep", tracer = quote({
          traced_diss_sum <<- c(traced_diss_sum, list(sum(diss)))
          traced_dmax <<- c(traced_dmax, list(dmax))
        }),
        exit = quote(traced_reps <<- c(traced_reps, list(returnValue()))),
        print = FALSE, where = asNamespace("TraMineR"))
}
W <- suppressMessages(seqdist(s, method = "OM", sm = "TRATE", indel = 1))
r_cases <- list(
  list(case = "r_lcs_dmax", args = list(diss = D, coverage = 0.5, pradius = 0.2),
       options = list(coverage = 0.5, pradius = 0.2), subset_diss = function(i) D[i, i]),
  list(case = "r_trate_method", args = list(method = "OM", sm = "TRATE", indel = 1),
       options = list(),
       subset_diss = function(i) suppressMessages(seqdist(s[i, ], method = "OM", sm = "TRATE", indel = 1)))
)
for (cs in r_cases) {
  traced_diss_sum <- list()
  traced_dmax <- list()
  traced_reps <- list()
  trace_seqrep()
  suppressMessages(do.call(seqplot, c(list(s, type = "r", group = g), cs$args)))
  untrace("seqrep", where = asNamespace("TraMineR"))
  stopifnot(length(traced_reps) == length(levs))
  told_apart <- FALSE
  for (np in seq_along(levs)) {
    i <- which(g == levs[np])
    Ds <- cs$subset_diss(i)
    subset <- suppressMessages(do.call(seqrep, c(list(s[i, ], diss = Ds), cs$options)))
    told_apart <- told_apart || !identical(rep_text(subset), rep_text(traced_reps[[np]])) ||
      max(Ds) != traced_dmax[[np]] || !isTRUE(all.equal(sum(Ds), traced_diss_sum[[np]]))
    add(cs$case, levs[np], "diss_sum", format(traced_diss_sum[[np]], digits = 17))
    add(cs$case, levs[np], "dmax", format(traced_dmax[[np]], digits = 17))
    add(cs$case, levs[np], "representatives", rep_text(traced_reps[[np]]))
    add(cs$case, levs[np], "na_pct", format(attr(traced_reps[[np]], "Statistics")[["na(%)"]], digits = 17))
    add(cs$case, levs[np], "subset_diss_sum", format(sum(Ds), digits = 17))
    add(cs$case, levs[np], "subset_dmax", format(max(Ds), digits = 17))
    add(cs$case, levs[np], "subset_representatives", rep_text(subset))
  }
  stopifnot(told_apart)
}
stopifnot(isTRUE(all.equal(max(W), as.numeric(traced_dmax[[1]]))))

# --- type "r" with a one-case group: R stops ---------------------------------
refused <- tryCatch({
  suppressMessages(seqplot(s, type = "r", group = c("g", rep("h", 59)), diss = D))
  "drawn"
}, error = function(e) conditionMessage(e))
stopifnot(!identical(refused, "drawn"))
add("r_one_case", "all", "error", refused)
# Without group the whole set is one group (grp <- rep(1, nrow(seqdata))).
refused <- tryCatch({
  suppressMessages(seqplot(s[1, ], type = "r", diss = D[1, 1, drop = FALSE]))
  "drawn"
}, error = function(e) conditionMessage(e))
stopifnot(!identical(refused, "drawn"))
add("r_one_sequence", "all", "error", refused)

# --- type "rf" with group and which.plot = "both": R stops (seqplot.R:53-56) ---
refused <- tryCatch({
  suppressMessages(seqplot(s, type = "rf", group = g, diss = D, k = 4, which.plot = "both"))
  "drawn"
}, error = function(e) conditionMessage(e))
stopifnot(!identical(refused, "drawn"))
add("rf_both_with_group", "all", "error", refused)

write.csv(do.call(rbind, rows), file.path(out_dir, "seqplot_group_rep.csv"), row.names = FALSE, quote = TRUE)
write.csv(data.frame(sequence = apply(as.matrix(s), 1, paste, collapse = "-"), sex = g,
                     weight = format(b$wp00tbgs, digits = 17), stringsAsFactors = FALSE),
          file.path(out_dir, "seqplot_group_rep_input.csv"), row.names = FALSE, quote = TRUE)
write.table(format(D, digits = 17), file.path(out_dir, "seqplot_group_rep_diss.csv"), sep = ",",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
