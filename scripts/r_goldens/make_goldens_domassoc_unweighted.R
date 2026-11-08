# Live-R goldens for seqdomassoc(weighted = FALSE) on weighted domains (record of
# audit finding TX-AUD-20261011-011).
#
# With weighted = FALSE, seqdomassoc puts unit weights in seqrep.args$weights
# (seqdomassoc.R:121-122), but seqrep has no weights argument and, with diss
# given, drops its ... (seqrep.R:5-24), so the representatives keep the object
# weights; for eq.group, seqrf keeps its default weighted = TRUE (seqdomassoc.R:86
# says weighted does not apply). traminex selects unweighted representatives and
# groups, a recorded correction. Its values are R's own on the same domains
# without weights ("unweighted" rows), where R's selection and tables are
# unweighted; R's values on the weighted domains with weighted = FALSE are kept
# ("r_weighted_false") to show the difference.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_domassoc_unweighted.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
icuSetCollate(locale = "root")
source(file.path("scripts", "r_goldens", "make_goldens_md_inputs.R"))

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

inp <- md_inputs()
w <- c(3, 1, 1, 1, 5, 1, 1, 1, 2, 1, 4, 1)
weighted <- lapply(inp$s, function(s) { attr(s, "weights") <- w; s })
unweighted <- inp$s
diss <- inp$d
stopifnot(length(diss) == 3, all(sapply(diss, nrow) == nrow(inp$s[[1]])))

rows <- list()
record <- function(method, kind, result) {
  stats <- unclass(result)
  for (pair in rownames(stats)) for (name in colnames(stats)) {
    rows[[length(rows) + 1L]] <<- data.frame(method = method, kind = kind, pair = pair, name = name,
                                             value = sprintf("%.17g", stats[pair, name]),
                                             stringsAsFactors = FALSE)
  }
}
# seqrf's default MDS axis has no sign of its own: R takes it from LAPACK, so the
# eq.group medoids follow the numerical library (TX-REP-20261017-006). The trace
# replaces seqrf's sortv = "mds" by R's own axis for that call (cmdscale, or
# wcmdscale with the normalised weights; seqrf.R:5-7, dissrf.R:22-25, 35-44,
# 86-91) oriented by traminex's rule: the coordinate of largest absolute value is
# positive (coordinates within a relative 1e-10 of it are tied and the first
# decides; NA ignored). Everything else is R's.
trace("seqrf", tracer = quote(if (length(sortv) == 1 && identical(as.character(sortv), "mds")) {
  w <- if (is.null(weights) && weighted) attr(seqdata, "weights") else weights
  if (!weighted) w <- NULL
  p <- 2^squared
  axis <- if (is.null(w)) cmdscale(diss^p, k = 1)[, 1] else
    vegan::wcmdscale(diss^p, k = 1, w = nrow(diss) * w / sum(w))[, 1]
  largest <- which(abs(axis) >= max(abs(axis), na.rm = TRUE) * (1 - 1e-10))[1]
  sortv <- if (axis[largest] < 0) -axis else axis
}), print = FALSE, where = asNamespace("TraMineR"))

calls <- list(
  rep = function(domains, weighted) TraMineR::seqdomassoc(domains, rep.method = "rep", diss.dom = diss,
                                                           weighted = weighted,
                                                           seqrep.args = list(coverage = .8, pradius = .1)),
  # k = 3 divides the 12 cases, so every unweighted group has the same height.
  # R's eq.group weights (heights * ncase, seqdomassoc.R:176) are then one
  # value for all medoids, and their recycling against the stacked medoid
  # sequences (:180) cannot misplace them: the unweighted twin is exact here
  # only because the heights are equal (batch-4 review, R4-l).
  eq.group = function(domains, weighted) TraMineR::seqdomassoc(domains, rep.method = "eq.group",
                                                                diss.dom = diss, weighted = weighted,
                                                                seqrf.args = list(k = 3))
)
for (method in names(calls)) {
  record(method, "r_weighted_false", suppressMessages(calls[[method]](weighted, FALSE)))
  record(method, "unweighted", suppressMessages(calls[[method]](unweighted, TRUE)))
}
write.csv(do.call(rbind, rows), file.path(out_dir, "domassoc_unweighted.csv"), row.names = FALSE, quote = TRUE)
write.csv(data.frame(weight = w), file.path(out_dir, "domassoc_unweighted_weights.csv"), row.names = FALSE)
