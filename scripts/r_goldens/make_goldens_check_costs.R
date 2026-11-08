# Live-R goldens for check.costs (TX-AUD-20261011-041).
#
# check.costs = TRUE (the default) runs checkcost (TraMineR-checkcost.R); FALSE skips
# it. The rows record seqdist on costs checkcost stops or flags, under both values:
# with FALSE, R computes with a negative indel, a non-zero diagonal or a matrix of
# the wrong size, and returns negative "distances" for negative costs; with TRUE an
# all-negative matrix also gives negative distances, because checkcost returns after
# its "max <= 0" message (:5-8) before the negative-cost check (:52-54). seqcost and
# seqMD rows show the same switch elsewhere.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_check_costs.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

s3 <- suppressMessages(TraMineR::seqdef(data.frame(a = c("A", "B", "A"), b = c("B", "B", "A"), c = c("A", "A", "B")),
                                        alphabet = c("A", "B")))
neg <- matrix(c(0, -1, -1, 0), 2)
cases <- list(
  negative_sm = function(cc) TraMineR::seqdist(s3, method = "OM", sm = neg, indel = 1, check.costs = cc),
  nonzero_diagonal = function(cc) TraMineR::seqdist(s3, method = "OM", sm = matrix(c(1, 2, 2, 1), 2), indel = 1, check.costs = cc),
  negative_indel = function(cc) TraMineR::seqdist(s3, method = "OM", sm = matrix(c(0, 2, 2, 0), 2), indel = -1, check.costs = cc),
  zero_indel = function(cc) TraMineR::seqdist(s3, method = "OM", sm = matrix(c(0, 2, 2, 0), 2), indel = 0, check.costs = cc),
  wrong_size = function(cc) TraMineR::seqdist(s3, method = "OM", sm = matrix(c(0, 2, 2, 2, 0, 2, 2, 2, 0), 3), indel = 1, check.costs = cc),
  triangle_violation = function(cc) TraMineR::seqdist(s3, method = "OM", sm = matrix(c(0, 5, 5, 0), 2), indel = 1, check.costs = cc),
  ham_negative_sm = function(cc) TraMineR::seqdist(s3, method = "HAM", sm = neg, check.costs = cc),
  seqcost_trate = function(cc) TraMineR::seqcost(s3, method = "TRATE", check.costs = cc)$sm,
  seqmd_negative_sm = function(cc) TraMineR::seqMD(list(s3, s3), method = "OM", sm = list(neg, neg), indel = 1,
                                                   what = "diss", check.costs = cc)
)
rows <- list()
for (name in names(cases)) {
  for (cc in c(TRUE, FALSE)) {
    result <- tryCatch(suppressWarnings(suppressMessages(cases[[name]](cc))), error = function(e) e)
    ok <- !inherits(result, "error")
    rows[[length(rows) + 1L]] <- data.frame(
      case = name, check_costs = cc, ok = ok, error = if (ok) "" else trimws(conditionMessage(result)),
      value = if (ok) paste(format(as.vector(as.matrix(result)), digits = 17), collapse = " ") else "",
      stringsAsFactors = FALSE)
  }
}
write.csv(do.call(rbind, rows), file.path(out_dir, "check_costs.csv"), row.names = FALSE, quote = TRUE)
