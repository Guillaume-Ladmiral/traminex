# Live-R goldens for seqdomassoc's seqrep.args and seqrf.args (repair of audit
# finding TX-AUD-20261011-012).
#
# seqdomassoc's defaults seqrep.args = list(coverage = .8, pradius = .1) and
# seqrf.args = list(k = 20) (seqdomassoc.R:79-80) are replaced by a list the user
# gives, so seqrep's and seqrf's own defaults then apply to the arguments the list
# leaves out (coverage = .25; k = NULL).
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_domassoc_args.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
icuSetCollate(locale = "root")
source(file.path("scripts", "r_goldens", "make_goldens_md_inputs.R"))

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

inp <- md_inputs()
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

cases <- list(
  rep_pradius = quote(TraMineR::seqdomassoc(inp$s, rep.method = "rep", diss.dom = inp$d,
                                            seqrep.args = list(pradius = .2))),
  rep_default = quote(TraMineR::seqdomassoc(inp$s, rep.method = "rep", diss.dom = inp$d)),
  eq_sortv = quote(TraMineR::seqdomassoc(inp$s, rep.method = "eq.group", diss.dom = inp$d,
                                         seqrf.args = list(sortv = "mds", k = 4))),
  eq_k_only = quote(TraMineR::seqdomassoc(inp$s, rep.method = "eq.group", diss.dom = inp$d,
                                          seqrf.args = list(grp.meth = "first")))
)
rows <- list()
for (name in names(cases)) {
  result <- tryCatch(unclass(suppressMessages(eval(cases[[name]]))), error = function(e) e)
  if (inherits(result, "error")) {
    rows[[length(rows) + 1L]] <- data.frame(case = name, pair = "", name = "error",
                                            value = conditionMessage(result), stringsAsFactors = FALSE)
    next
  }
  for (pair in rownames(result)) for (column in colnames(result)) {
    rows[[length(rows) + 1L]] <- data.frame(case = name, pair = pair, name = column,
                                            value = sprintf("%.17g", result[pair, column]),
                                            stringsAsFactors = FALSE)
  }
}
write.csv(do.call(rbind, rows), file.path(out_dir, "domassoc_args.csv"), row.names = FALSE, quote = TRUE)
