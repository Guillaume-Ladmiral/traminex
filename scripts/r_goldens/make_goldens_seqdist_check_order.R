# Live-R goldens for seqdist's argument checks (batch-6 review M-8 and M-9; TX-AUD-20261019-079).
#
# R checks norm, indel, sm, prox, link and step in that order for every method (seqdist.R:127-190), so with
# two invalid arguments R names the earlier one. The "order_" rows pair every two adjacent checks and the
# four pairs of batch-6 review M-9, and a logical or complex sm after an invalid norm or indel (batch-9 review M-1,
# M-5; "sm_complex" alone), and a norm spelled with spaces (M-4). "dhd_constant" records R's refusal of sm = "CONSTANT" for DHD
# (seqdist.R:253). The "check_costs_" rows record how R reads check.costs with if(): each row's error, or
# that R computed (and the first row of the distances); the "bad_" rows give a non-zero diagonal, which R refuses
# when it reads check.costs as TRUE and computes with when it reads FALSE. The "check_costs_unread_" rows are
# routes on which R never reads check.costs (seqdist.R:370-460 reads it only for a matrix sm of the OM family and
# TWED, an array for DHD and costs seqcost generates; seqMD returns what = "MDseq" at seqMD.R:201, before its
# read at :266), so R computes whatever the value (batch-9 review I-2).
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_seqdist_check_order.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

x <- suppressMessages(TraMineR::seqdef(data.frame(a = c("A", "A", "C", "B"), b = c("A", "B", "C", "B"),
                                                  c = c("B", "B", "B", "C"), d = c("C", "C", "A", "A")),
                                       alphabet = c("A", "B", "C")))
bad <- matrix(c(1, 2, 2, 2, 1, 2, 2, 2, 1), 3, 3)  # non-zero diagonal
good <- matrix(c(0, 2, 2, 2, 0, 2, 2, 2, 0), 3, 3)
cases <- list(
  order_norm_indel = quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", norm = "bogus", indel = "bogus")),
  order_indel_sm = quote(TraMineR::seqdist(x, method = "OM", indel = "bogus", sm = "bogus")),
  order_sm_prox = quote(TraMineR::seqdist(x, method = "OM", sm = "bogus", prox = "bogus")),
  order_prox_link = quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", prox = "bogus", link = "bogus")),
  order_link_step = quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", link = "bogus", step = 0)),
  order_norm_link = quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", norm = "bogus", link = "bogus")),
  order_indel_step = quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", indel = "bogus", step = 0)),
  order_norm_step = quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", norm = "bogus", step = 0)),
  order_norm_logical_sm = quote(TraMineR::seqdist(x, method = "OM", norm = "bogus", sm = TRUE)),
  order_indel_logical_sm = quote(TraMineR::seqdist(x, method = "OM", indel = "bogus", sm = TRUE)),
  order_norm_complex_sm = quote(TraMineR::seqdist(x, method = "OM", norm = "bogus", sm = 2 + 0i)),
  sm_complex = quote(TraMineR::seqdist(x, method = "OM", sm = 2 + 0i)),
  order_norm_spaces_indel = quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", norm = " maxlength ", indel = "bogus")),
  dhd_constant = quote(TraMineR::seqdist(x, method = "DHD", sm = "CONSTANT")),
  check_costs_TRUE = quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", check.costs = TRUE)),
  check_costs_FALSE = quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", check.costs = FALSE)),
  check_costs_string_TRUE = quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", check.costs = "TRUE")),
  check_costs_string_T = quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", check.costs = "T")),
  check_costs_string_true = quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", check.costs = "true")),
  check_costs_string_FALSE = quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", check.costs = "FALSE")),
  check_costs_string_F = quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", check.costs = "F")),
  check_costs_string_false = quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", check.costs = "false")),
  check_costs_string_no = quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", check.costs = "no")),
  check_costs_number_0 = quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", check.costs = 0)),
  check_costs_number_1 = quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", check.costs = 1)),
  check_costs_number_2_5 = quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", check.costs = 2.5)),
  check_costs_NA = quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", check.costs = NA)),
  check_costs_NULL = quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", check.costs = NULL)),
  check_costs_bad_TRUE = quote(TraMineR::seqdist(x, method = "OM", sm = bad, indel = 1, check.costs = TRUE)),
  check_costs_bad_string_T = quote(TraMineR::seqdist(x, method = "OM", sm = bad, indel = 1, check.costs = "T")),
  check_costs_bad_number_1 = quote(TraMineR::seqdist(x, method = "OM", sm = bad, indel = 1, check.costs = 1)),
  check_costs_bad_FALSE = quote(TraMineR::seqdist(x, method = "OM", sm = bad, indel = 1, check.costs = FALSE)),
  check_costs_bad_string_FALSE = quote(TraMineR::seqdist(x, method = "OM", sm = bad, indel = 1, check.costs = "FALSE")),
  check_costs_bad_string_F = quote(TraMineR::seqdist(x, method = "OM", sm = bad, indel = 1, check.costs = "F")),
  check_costs_bad_number_0 = quote(TraMineR::seqdist(x, method = "OM", sm = bad, indel = 1, check.costs = 0)),
  check_costs_bad_string_TRUE = quote(TraMineR::seqdist(x, method = "OM", sm = bad, indel = 1, check.costs = "TRUE")),
  check_costs_bad_string_true = quote(TraMineR::seqdist(x, method = "OM", sm = bad, indel = 1, check.costs = "true")),
  check_costs_bad_string_false = quote(TraMineR::seqdist(x, method = "OM", sm = bad, indel = 1, check.costs = "false")),
  check_costs_bad_number_2_5 = quote(TraMineR::seqdist(x, method = "OM", sm = bad, indel = 1, check.costs = 2.5)),
  check_costs_unread_lcs_no = quote(TraMineR::seqdist(x, method = "LCS", check.costs = "no")),
  check_costs_unread_lcs_NA = quote(TraMineR::seqdist(x, method = "LCS", check.costs = NA)),
  check_costs_unread_lcs_NULL = quote(TraMineR::seqdist(x, method = "LCS", check.costs = NULL)),
  check_costs_unread_ham_matrix_no = quote(TraMineR::seqdist(x, method = "HAM", sm = good, check.costs = "no")),
  check_costs_unread_ham_default_no = quote(TraMineR::seqdist(x, method = "HAM", check.costs = "no")),
  check_costs_unread_dhd_default_no = quote(TraMineR::seqdist(x, method = "DHD", check.costs = "no")),
  check_costs_unread_chi2_no = quote(TraMineR::seqdist(x, method = "CHI2", check.costs = "no")),
  check_costs_unread_seqmd_mdseq_no = quote(as.matrix(TraMineR::seqMD(list(x, x), what = "MDseq", check.costs = "no")))
)
rows <- list()
for (name in names(cases)) {
  err <- ""
  value <- ""
  ok <- tryCatch({
    d <- suppressMessages(suppressWarnings(eval(cases[[name]])))
    value <- paste(format(as.vector(d[1, ]), digits = 15), collapse = " ")
    TRUE
  }, error = function(e) { err <<- conditionMessage(e); FALSE })
  rows[[length(rows) + 1L]] <- data.frame(case = name, ok = ok, error = err, value = value, stringsAsFactors = FALSE)
}
write.csv(do.call(rbind, rows), file.path(out_dir, "seqdist_check_order.csv"), row.names = FALSE, quote = TRUE)
