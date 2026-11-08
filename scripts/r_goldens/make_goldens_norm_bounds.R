# Live-R goldens for the normalization bound of OMloc, OMspell and TWED (TX-REP-20261021-003 to -005).
#
# seqdist.Rd defines maxdist as d/m, gmean as 1 - (m - d)/(2 sqrt(l1 l2)) in R's code and YujianBo as 2d/(d + m),
# m the maximum possible distance (or an upper bound) given the two lengths; each is at most 1 exactly when d <= m.
# R's m is too small for these three methods (OMVIdistance.cpp, OMPerdistance.cpp:108, TWEDdistance.cpp:153), so
# its normalized values exceed 1. Each case records R's raw distances (norm = "none") and its maxdist, gmean and
# YujianBo values, whole matrices column by column: the raw distances are what traminex keeps, the normalized ones
# what it corrects. The last three cases (batch-10 review, TX-REP-20261022-002 and -003) record R's gmean inside
# [0, 1] where traminex's true bound makes it negative (OMloc, TWED), and R's negative OMspell distance when a spell's
# indel cost indel + expcost * (t^tpow - 1) is negative (tpow < 0), which traminex refuses.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_norm_bounds.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

build <- function(rows, alphabet) {
  width <- max(nchar(rows))
  cells <- t(sapply(rows, function(r) { v <- strsplit(r, "")[[1]]; c(v, rep(NA, width - length(v))) }))
  if (width == 1) cells <- matrix(cells, ncol = 1)
  suppressMessages(seqdef(as.data.frame(cells), alphabet = alphabet, right = "DEL"))
}
ab <- c("A", "B")
abcd <- c("A", "B", "C", "D")
data <- list(
  a_bbb = build(c("A", "BBB"), ab),
  four3 = build(c("BBB", "AAA", "BAB", "ABB"), ab),
  ragged2 = build(c("BBB", "A", "BAB", "AB"), ab),
  three = build(c("ACB", "CC", "BABC", "C", "ABCA"), c("A", "B", "C")),
  alternating = build(c("AAAAAAAAAA", "ABABABABAB"), ab),
  x1 = build(c("BDDDDC", "CBCBC", "BADBCA", "DCB", "DB", "A"), abcd),
  x2 = build(c("BDBAB", "DAA", "A", "DABDC", "CCDADDC", "ABBBDA"), abcd),
  a_aab = build(c("A", "AAB"), ab),
  a_b = build(c("A", "B"), ab),
  alt8 = build(c("AABBAABB", "BBAABBAA"), ab)
)
# An empty sequence beside ABA (R's TWED gives d = n * indel there).
data$empty_aba <- suppressMessages(seqdef(rbind(c(NA, NA, NA), c("A", "B", "A")), alphabet = ab, right = "DEL"))
asym <- matrix(c(0, 3, 0.1, 0), 2, 2)  # sm[A, B] = 0.1, sm[B, A] = 3
sym3 <- matrix(c(0, 3, 3, 0), 2, 2)
asym3 <- matrix(c(0, 2.5, 0.4, 1, 0, 3, 0.2, 1.5, 0), 3, 3)
small <- matrix(c(0, 0.1, 0.1, 0), 2, 2)
S1 <- matrix(c(0, 1.871, 1.957, 2.047, 1.871, 0, 3.639, 2.424, 1.957, 3.639, 0, 1.551, 2.047, 2.424, 1.551, 0), 4, 4, byrow = TRUE)
S2 <- matrix(c(0, 2.631, 1.99, 1.159, 2.631, 0, 3.697, 2.168, 1.99, 3.697, 0, 1.425, 1.159, 2.168, 1.425, 0), 4, 4, byrow = TRUE)
cases <- list(
  omloc_constant_e025 = list("a_bbb", list(method = "OMloc", sm = "CONSTANT", expcost = 0.25)),
  omloc_asym_e05 = list("four3", list(method = "OMloc", sm = asym, expcost = 0.5)),
  omloc_asym_e025 = list("four3", list(method = "OMloc", sm = asym, expcost = 0.25)),
  omloc_asym_e0 = list("four3", list(method = "OMloc", sm = asym, expcost = 0)),
  omloc_sym_ragged_e01 = list("ragged2", list(method = "OMloc", sm = sym3, expcost = 0.1)),
  omloc_asym3_context15 = list("three", list(method = "OMloc", sm = asym3, expcost = 0.3, context = 1.5)),
  omspell_small_e0 = list("alternating", list(method = "OMspell", sm = small, indel = 1, expcost = 0)),
  omspell_x1_tpow2 = list("x1", list(method = "OMspell", sm = S1, indel = 0.5, expcost = 1, tpow = 2)),
  omspell_x1_tpow05 = list("x1", list(method = "OMspell", sm = S1, indel = 1, expcost = 2, tpow = 0.5)),
  omspell_x1_tpowm05 = list("x1", list(method = "OMspell", sm = S1, indel = 1, expcost = 2, tpow = -0.5)),
  omspell_x1_default = list("x1", list(method = "OMspell", sm = "CONSTANT")),
  twed_x1_h2_nu2 = list("x1", list(method = "TWED", sm = S1, indel = 0.5, h = 2, nu = 2)),
  twed_x2_h0_nu05 = list("x2", list(method = "TWED", sm = S2, indel = 0.5, h = 0, nu = 0.5)),
  twed_empty = list("empty_aba", list(method = "TWED", sm = "CONSTANT", nu = 0.5, h = 0.5)),
  twed_x1_auto_indel = list("x1", list(method = "TWED", sm = S1, h = 0.5, nu = 0.5)),
  omloc_gmean_a_aab = list("a_aab", list(method = "OMloc", sm = matrix(c(0, 4, 4, 0), 2, 2), expcost = 0, context = 0.5)),
  twed_gmean_a_b = list("a_b", list(method = "TWED", sm = matrix(c(0, 2.5, 0.4, 0), 2, 2), indel = 2, nu = 1, h = 0)),
  omspell_negative_spell_indel = list("alt8", list(method = "OMspell", sm = matrix(c(0, 4, 4, 0), 2, 2), indel = 0.3, expcost = 3, tpow = -2))
)
rows <- list()
for (name in names(cases)) {
  seqs <- data[[cases[[name]][[1]]]]
  for (nm in c("none", "maxdist", "gmean", "YujianBo")) {
    d <- suppressMessages(suppressWarnings(do.call(seqdist, c(list(seqs), cases[[name]][[2]], list(norm = nm)))))
    rows[[length(rows) + 1L]] <- data.frame(case = name, norm = nm, value = paste(format(as.vector(d), digits = 17), collapse = " "))
  }
}
write.csv(do.call(rbind, rows), file.path(out_dir, "norm_bounds.csv"), row.names = FALSE, quote = TRUE)
