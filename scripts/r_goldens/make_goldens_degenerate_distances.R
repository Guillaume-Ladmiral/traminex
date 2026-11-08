# Live-R goldens for seqdist on degenerate inputs (TX-AUD-20261011-038, -039, -040).
#
# Transition-rate costs: R's seqcost(TRATE) and seqtrate stop on a one-state
# alphabet ("subscript out of bounds") and on sequences of length 1 ("undefined
# columns selected"), and so do seqdist's TRATE and DHD routes. Other degenerate
# inputs on which R's value or error departs from the method's definition: LCS
# gmean/maxdist with a one-state alphabet (R's 1x1 cost matrix), CHI2/EUCLID with
# one sequence or one state ("'x' must be an array of at least two dimensions"),
# and OMstran with a one-element transition alphabet (its loop over
# 1:(alphabet_size - 1), seqdist-OMstran.R:126). Each row is a value (column-major,
# 17 digits) or R's error.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_degenerate_distances.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ragged <- data.frame(a = c("A", "A", "A"), b = c("A", "A", "A"), c = c("A", NA, "A"), d = c("A", NA, NA))
one <- suppressMessages(TraMineR::seqdef(ragged, alphabet = "A"))
two <- suppressMessages(TraMineR::seqdef(ragged, alphabet = c("A", "B")))
one_full <- suppressMessages(TraMineR::seqdef(data.frame(a = c("A", "A"), b = c("A", "A")), alphabet = "A"))
len1 <- suppressMessages(TraMineR::seqdef(data.frame(a = c("A", "B", "A")), alphabet = c("A", "B")))
s3 <- suppressMessages(TraMineR::seqdef(data.frame(a = c("A", "B", "A"), b = c("B", "B", "A"), c = c("A", "A", "B")),
                                        alphabet = c("A", "B")))
s1row <- s3[1, ]
five <- data.frame(a = c("A", "A"), b = c("A", NA), c = c("A", NA), d = c("A", NA), e = c("A", NA))
one5 <- suppressMessages(TraMineR::seqdef(five, alphabet = "A"))
two5 <- suppressMessages(TraMineR::seqdef(five, alphabet = c("A", "B")))
cases <- list(
  # transition-rate costs (TX-AUD-20261011-039, -040)
  om_trate_one_state = quote(TraMineR::seqdist(one, method = "OM", sm = "TRATE")),
  omstran_trate_one_state = quote(TraMineR::seqdist(one, method = "OMstran", sm = "TRATE", otto = 0.5)),
  seqcost_trate_one_state = quote(TraMineR::seqcost(one, method = "TRATE")$sm),
  dhd_one_state = quote(TraMineR::seqdist(one_full, method = "DHD")),
  dhd_trate_one_state = quote(TraMineR::seqdist(one_full, method = "DHD", sm = "TRATE")),
  om_trate_len1 = quote(TraMineR::seqdist(len1, method = "OM", sm = "TRATE")),
  dhd_len1 = quote(TraMineR::seqdist(len1, method = "DHD")),
  dhd_trate_len1 = quote(TraMineR::seqdist(len1, method = "DHD", sm = "TRATE")),
  ham_trate_len1 = quote(TraMineR::seqdist(len1, method = "HAM", sm = "TRATE")),
  seqcost_trate_len1 = quote(TraMineR::seqcost(len1, method = "TRATE")$sm),
  seqtrate_len1 = quote(TraMineR::seqtrate(len1)),
  # LCS normalisations with a one-state alphabet (TX-AUD-20261011-038)
  lcs_gmean_one_state = quote(TraMineR::seqdist(one, method = "LCS", norm = "gmean")),
  lcs_maxdist_one_state = quote(TraMineR::seqdist(one, method = "LCS", norm = "maxdist")),
  lcs_yujianbo_one_state = quote(TraMineR::seqdist(one, method = "LCS", norm = "YujianBo")),
  lcs_gmean_two_states = quote(TraMineR::seqdist(two, method = "LCS", norm = "gmean")),
  lcs_maxdist_two_states = quote(TraMineR::seqdist(two, method = "LCS", norm = "maxdist")),
  lcs_yujianbo_two_states = quote(TraMineR::seqdist(two, method = "LCS", norm = "YujianBo")),
  # CHI2 / EUCLID with one sequence or one state (TX-AUD-20261011-040)
  chi2_one_row = quote(TraMineR::seqdist(s1row, method = "CHI2")),
  euclid_one_row = quote(TraMineR::seqdist(s1row, method = "EUCLID")),
  chi2_one_state = quote(TraMineR::seqdist(one_full, method = "CHI2")),
  euclid_one_state = quote(TraMineR::seqdist(one_full, method = "EUCLID")),
  # OMstran with a one-element transition alphabet (TX-AUD-20261011-039)
  omstran_constant_one_state = quote(TraMineR::seqdist(one5, method = "OMstran", sm = "CONSTANT", otto = 0.5)),
  omstran_constant_two_states = quote(TraMineR::seqdist(two5, method = "OMstran", sm = "CONSTANT", otto = 0.5))
)
rows <- list()
for (name in names(cases)) {
  result <- tryCatch(suppressWarnings(suppressMessages(eval(cases[[name]]))), error = function(e) e)
  ok <- !inherits(result, "error")
  rows[[length(rows) + 1L]] <- data.frame(
    case = name, ok = ok, error = if (ok) "" else trimws(conditionMessage(result)),
    value = if (ok) paste(format(as.vector(as.matrix(result)), digits = 17), collapse = " ") else "",
    stringsAsFactors = FALSE)
}
write.csv(do.call(rbind, rows), file.path(out_dir, "degenerate_distances.csv"), row.names = FALSE, quote = TRUE)
