stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
# Live-R goldens for normalized OM-family distances under an asymmetric substitution-cost matrix
# (audit finding TX-AUD-20261019-008, an R defect corrected by traminex).
#
# R's OM-family classes take maxscost from the strict upper triangle of sm (OMdistance.cpp:26-33)
# and bound a distance by |n-m|*indel + maxscost*min(m, n) (OMdistance.cpp:124; OMPerdistance.cpp:108,
# OMPerdistanceII.cpp:109; TWEDdistance.cpp:153). seqdist.Rd defines that bound as the maximum possible
# distance, but with a larger cost below the diagonal R's normalized maxdist and gmean exceed 1. These
# rows record R's values; traminex's tests assert the corrected values (bound from the whole
# off-diagonal) and that they differ. OMloc and OMstran, whose bounds held in every case checked, are
# recorded so that their unchanged agreement with R stays tested.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_asym_bound.R

suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
quiet <- function(expr) suppressWarnings(suppressMessages(expr))
rows <- list()
record <- function(case, x, method, norm, ...) {
  d <- quiet(TraMineR::seqdist(x, method = method, norm = norm, ...))
  if (is.null(dim(d)) && !inherits(d, "dist")) {
    m <- matrix(as.numeric(d), ncol = 1)
    idx <- cbind(seq_along(d), rep(0L, length(d)))
  } else {
    m <- as.matrix(d)
    idx <- as.matrix(expand.grid(i = seq_len(nrow(m)), j = seq_len(ncol(m))))
  }
  value <- if (ncol(m) == 1L && all(idx[, 2] == 0L)) m[, 1] else m[idx]
  rows[[length(rows) + 1L]] <<- data.frame(case = case, method = method, norm = norm, i = idx[, 1],
                                          j = idx[, 2], value = sprintf("%.17g", value))
}
norms <- c("none", "maxdist", "gmean", "YujianBo", "maxlength")

# Worked by hand in the tests: two states, the larger cost below the diagonal, indel 2.
hand <- matrix(c(0, 1, 3, 0), 2, byrow = TRUE)
for (order in list(list("hand_b_first", rbind(c("B", "B"), c("A", "A"))), list("hand_a_first", rbind(c("A", "A"), c("B", "B"))))) {
  x <- quiet(TraMineR::seqdef(order[[2]], alphabet = c("A", "B")))
  for (n in norms) record(order[[1]], x, "OM", n, sm = hand, indel = 2)
}

# The audit's case (lead A1-04): four states, costs 3 below the diagonal in the first column.
A <- matrix(c(0,1,1,1, 3,0,1,1, 3,1,0,1, 3,1,1,0), 4, byrow = TRUE)
m <- rbind(strsplit("ABCD","")[[1]], strsplit("BBCD","")[[1]], strsplit("DCBA","")[[1]], strsplit("AAAA","")[[1]],
           strsplit("CCCC","")[[1]], strsplit("BADC","")[[1]])
x <- quiet(TraMineR::seqdef(m, alphabet = c("A","B","C","D")))
for (n in norms) {
  record("audit", x, "OM", n, sm = A, indel = 2)
  record("audit", x, "OMslen", n, sm = A, indel = 2)
  record("audit", x, "OMspell", n, sm = A, indel = 2)
  record("audit", x, "TWED", n, sm = A, indel = 2, h = 0.5, nu = 0.5)
  record("audit", x, "OMloc", n, sm = A)
  record("audit", x, "OMstran", n, sm = A, indel = 2, otto = 0.5)
}
record("audit_refseq", x, "OM", "maxdist", sm = A, indel = 2, refseq = 2)

# Ragged data (void at the end) and the largest TWED and OM cases of a seeded random search
# (inputs only; the values are R's).
S_twed <- matrix(c(0.0, 0.6733111123483121, 0.21390833763794703, 0.7168896906456412, 3.8856502611627395, 0.0, 0.7410285674363155, 0.43170017778345815, 2.6979008942483964, 3.1296502410905918, 0.0, 0.3394614040410233, 3.85890450310538, 3.58931482401455, 2.4626499838039315, 0.0), 4, byrow = TRUE)
x <- quiet(TraMineR::seqdef(rbind(
    c("B", "C", "D", NA, NA),
    c("C", "C", "D", "C", NA),
    c("D", "A", "D", "B", "D"),
    c("C", "B", "A", "B", NA),
    c("D", "D", "A", "A", "B"),
    c("D", "D", "A", "B", NA)), alphabet = c("A","B","C","D"), right = "DEL"))
for (n in c("none", "maxdist", "gmean")) record("search_twed", x, "TWED", n, sm = S_twed, indel = 2, h = 0.5, nu = 0.5)
S_om <- matrix(c(0.0, 0.6733111123483121, 0.21390833763794703, 0.7168896906456412, 3.8856502611627395, 0.0, 0.7410285674363155, 0.43170017778345815, 2.6979008942483964, 3.1296502410905918, 0.0, 0.3394614040410233, 3.85890450310538, 3.58931482401455, 2.4626499838039315, 0.0), 4, byrow = TRUE)
x <- quiet(TraMineR::seqdef(rbind(
    c("B", "C", "D", NA, NA),
    c("C", "C", "D", "C", NA),
    c("D", "A", "D", "B", "D"),
    c("C", "B", "A", "B", NA),
    c("D", "D", "A", "A", "B"),
    c("D", "D", "A", "B", NA)), alphabet = c("A","B","C","D"), right = "DEL"))
for (n in c("none", "maxdist", "gmean")) record("search_om", x, "OM", n, sm = S_om, indel = 2)

write.csv(do.call(rbind, rows), file.path(out_dir, "asym_bound.csv"), row.names = FALSE, quote = TRUE)
