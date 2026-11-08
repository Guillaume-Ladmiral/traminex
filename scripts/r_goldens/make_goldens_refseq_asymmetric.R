# Live-R goldens for reference distances with asymmetric costs (TX-REP-20261017-001).
#
# With an asymmetric substitution matrix, d(x, y) and d(y, x) differ. seqdist's full
# matrix holds d(x_i, x_j) for i < j over the distinct sequences in first-occurrence
# order, mirrored; a reference distance is d(x_i, reference) with the sequence first,
# and index sets give d(set-1 element, set-2 element). Rows: values (column-major,
# 17 digits). Data d5 repeats sequences before and after the reference. The
# normalized cases divide each oriented distance by R's bound (batch-6 review C-1).
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_refseq_asymmetric.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

d4 <- suppressMessages(TraMineR::seqdef(data.frame(a = c("A", "B", "C", "A"), b = c("B", "C", "A", "A"),
                                                   c = c("C", "A", "B", "B")), alphabet = c("A", "B", "C")))
d5 <- suppressMessages(TraMineR::seqdef(data.frame(a = c("A", "B", "A", "C", "B"), b = c("B", "C", "B", "A", "C"),
                                                   c = c("C", "A", "C", "B", "A")), alphabet = c("A", "B", "C")))
asym <- matrix(c(0, 1, 3, 2, 0, 1, 1, 4, 0), 3, byrow = TRUE)
om <- function(s, ...) TraMineR::seqdist(s, method = "OM", sm = asym, indel = 2, ...)
ham <- function(s, ...) TraMineR::seqdist(s, method = "HAM", sm = asym, ...)
arr <- array(asym, c(3, 3, 3))
dhd <- function(s, ...) TraMineR::seqdist(s, method = "DHD", sm = arr, ...)
sm3 <- matrix(c(0, 1, 5, 3, 0, 1, 1, 4, 0), 3, byrow = TRUE)
om3 <- function(s, indel = 10, ...) TraMineR::seqdist(s, method = "OM", sm = sm3, indel = indel, ...)
dup <- suppressMessages(TraMineR::seqdef(data.frame(a = c("A", "C", "A"), b = c("A", "C", "A"),
                                                    c = c("B", "C", "B")), alphabet = c("A", "B", "C")))
vd <- suppressMessages(TraMineR::seqdef(data.frame(a = c("A", "B", "C"), b = c("A", "C", "C"), c = c("B", "C", "C"),
                                                   e = c(NA, NA, "C")), alphabet = c("A", "B", "C"), right = "DEL"))
cases <- list(
  d4_om_full = quote(om(d4)),
  d4_om_refseq_1 = quote(om(d4, refseq = 1)), d4_om_refseq_2 = quote(om(d4, refseq = 2)),
  d4_om_refseq_3 = quote(om(d4, refseq = 3)), d4_om_refseq_4 = quote(om(d4, refseq = 4)),
  d4_om_refseq_0 = quote(om(d4, refseq = 0)),
  d4_om_sets_12_34 = quote(om(d4, refseq = list(1:2, 3:4))), d4_om_sets_34_12 = quote(om(d4, refseq = list(3:4, 1:2))),
  d4_om_external_2 = quote(om(d4, refseq = d4[2, ])),
  d4_ham_refseq_1 = quote(ham(d4, refseq = 1)), d4_ham_refseq_3 = quote(ham(d4, refseq = 3)),
  d4_ham_sets_34_12 = quote(ham(d4, refseq = list(3:4, 1:2))), d4_ham_external_2 = quote(ham(d4, refseq = d4[2, ])),
  d5_om_full = quote(om(d5)),
  d5_om_refseq_1 = quote(om(d5, refseq = 1)), d5_om_refseq_2 = quote(om(d5, refseq = 2)),
  d5_om_refseq_4 = quote(om(d5, refseq = 4)), d5_om_refseq_5 = quote(om(d5, refseq = 5)),
  d5_om_refseq_0 = quote(om(d5, refseq = 0)),
  d5_om_external_2 = quote(om(d5, refseq = d5[2, ])),
  d5_om_sets_45_12 = quote(om(d5, refseq = list(4:5, 1:2))),
  d4_om_auto_refseq_1 = quote(om(d4, norm = "auto", refseq = 1)),
  d4_om_auto_refseq_3 = quote(om(d4, norm = "auto", refseq = 3)),
  d4_om_auto_refseq_0 = quote(om(d4, norm = "auto", refseq = 0)),
  d4_om_auto_sets_34_12 = quote(om(d4, norm = "auto", refseq = list(3:4, 1:2))),
  d4_om_auto_external_2 = quote(om(d4, norm = "auto", refseq = d4[2, ])),
  d4_om_maxlength_refseq_1 = quote(om(d4, norm = "maxlength", refseq = 1)),
  d4_om_gmean_refseq_1 = quote(om(d4, norm = "gmean", refseq = 1)),
  d4_om_maxdist_refseq_3 = quote(om(d4, norm = "maxdist", refseq = 3)),
  d4_om_yujianbo_refseq_1 = quote(om(d4, norm = "YujianBo", refseq = 1)),
  d4_ham_auto_refseq_3 = quote(ham(d4, norm = "auto", refseq = 3)),
  d4_ham_auto_sets_34_12 = quote(ham(d4, norm = "auto", refseq = list(3:4, 1:2))),
  d4_dhd_refseq_3 = quote(dhd(d4, refseq = 3)),
  d4_dhd_auto_refseq_3 = quote(dhd(d4, norm = "auto", refseq = 3)),
  d4_dhd_auto_refseq_1 = quote(dhd(d4, norm = "auto", refseq = 1)),
  d5_om_auto_refseq_2 = quote(om(d5, norm = "auto", refseq = 2)),
  d5_om_auto_refseq_0 = quote(om(d5, norm = "auto", refseq = 0)),
  # The low-level pairwise(ref_sets=) and pairwise_between_sets take the same orientation
  # (TX-AUD-20261019-009): a reference set holding a duplicate of a set-1 sequence, void data
  # with a scalar and a state-vector indel, and sets that share an index.
  dup_om_sets_12_3 = quote(om3(dup, refseq = list(1:2, 3))),
  dup_om_sets_3_2 = quote(om3(dup, refseq = list(3, 2))),
  dup_om_sets_12_23 = quote(om3(dup, refseq = list(1:2, 2:3))),
  void_om_sets_3_1 = quote(om3(vd, refseq = list(3, 1))),
  void_om_ivec_sets_3_1 = quote(om3(vd, indel = c(10, 8, 12), refseq = list(3, 1))),
  void_om_sets_12_23 = quote(om3(vd, refseq = list(1:2, 2:3)))
)
rows <- list()
for (name in names(cases)) {
  result <- suppressWarnings(suppressMessages(eval(cases[[name]])))
  rows[[length(rows) + 1L]] <- data.frame(case = name, nrow = NROW(as.matrix(result)),
                                          value = paste(format(as.vector(as.matrix(result)), digits = 17), collapse = " "),
                                          stringsAsFactors = FALSE)
}
write.csv(do.call(rbind, rows), file.path(out_dir, "refseq_asymmetric.csv"), row.names = FALSE, quote = TRUE)
