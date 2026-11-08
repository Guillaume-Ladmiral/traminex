# Live-R goldens for checkcost's messages on generated cost matrices (TX-AUD-20261011-080).
#
# seqcost runs checkcost on every matrix it generates when check.costs is TRUE
# (seqcost.R:286); seqdist resolves a cost-method string through seqcost once
# (seqdist.R:426-431; its second check, :445, is commented out). The rows record, in
# order, every "[!" message (message() or warning()) of each call; a call with none
# has one row with n = 0. The refseq cases check each message once per call on
# every reference route (batch-6 review I-1).
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_generated_cost_messages.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tri <- suppressMessages(TraMineR::seqdef(rbind(strsplit("ABABABAB", "")[[1]], strsplit("BCBCBCBC", "")[[1]]),
                                         alphabet = c("A", "B", "C")))
tri3 <- suppressMessages(TraMineR::seqdef(rbind(strsplit("ABABABAB", "")[[1]], strsplit("BCBCBCBC", "")[[1]],
                                                strsplit("ABBCABBC", "")[[1]]), alphabet = c("A", "B", "C")))
# Breaks the triangle inequality: d(A, C) = 3 > d(A, B) + d(B, C) = 2.
tri_sm <- matrix(c(0, 1, 3, 1, 0, 1, 3, 1, 0), 3, byrow = TRUE)
alt <- suppressMessages(TraMineR::seqdef(rbind(c("A", "B", "A", "B"), c("B", "A", "B", "A"), c("A", "B", "A", "B")),
                                         alphabet = c("A", "B")))
one_full <- suppressMessages(TraMineR::seqdef(data.frame(a = c("A", "A"), b = c("A", "A")), alphabet = "A"))
one <- suppressMessages(TraMineR::seqdef(data.frame(a = c("A", "A", "A"), b = c("A", "A", "A"), c = c("A", NA, "A"),
                                                    d = c("A", NA, NA)), alphabet = "A"))
cases <- list(
  seqcost_trate = quote(TraMineR::seqcost(tri, method = "TRATE")),
  seqsubm_trate = quote(TraMineR::seqsubm(tri, method = "TRATE")),
  seqdist_om_trate_indel1 = quote(TraMineR::seqdist(tri, method = "OM", sm = "TRATE", indel = 1)),
  seqdist_om_trate_auto = quote(TraMineR::seqdist(tri, method = "OM", sm = "TRATE")),
  seqcost_constant_cval0 = quote(TraMineR::seqcost(tri, method = "CONSTANT", cval = 0)),
  seqcost_constant = quote(TraMineR::seqcost(tri, method = "CONSTANT")),
  seqcost_constant_cval5 = quote(TraMineR::seqcost(tri, method = "CONSTANT", cval = 5)),
  seqdist_om_constant = quote(TraMineR::seqdist(tri, method = "OM", sm = "CONSTANT")),
  seqcost_indels = quote(TraMineR::seqcost(tri, method = "INDELS")),
  seqcost_indelslog = quote(TraMineR::seqcost(tri, method = "INDELSLOG")),
  seqdist_om_indels = quote(TraMineR::seqdist(tri, method = "OM", sm = "INDELS")),
  seqdist_om_indelslog = quote(TraMineR::seqdist(tri, method = "OM", sm = "INDELSLOG")),
  seqcost_indels_one_state = quote(TraMineR::seqcost(one, method = "INDELS")),
  seqcost_indelslog_one_state = quote(TraMineR::seqcost(one, method = "INDELSLOG")),
  seqdist_om_indels_one_state = quote(TraMineR::seqdist(one, method = "OM", sm = "INDELS")),
  seqdist_om_indelslog_one_state = quote(TraMineR::seqdist(one, method = "OM", sm = "INDELSLOG")),
  seqcost_future = quote(TraMineR::seqcost(tri, method = "FUTURE")),
  seqcost_trate_time_varying = quote(TraMineR::seqcost(tri, method = "TRATE", time.varying = TRUE)),
  seqdist_dhd = quote(TraMineR::seqdist(tri, method = "DHD")),
  seqdist_ham_trate = quote(TraMineR::seqdist(tri, method = "HAM", sm = "TRATE")),
  # LCS and HAM without sm generate constant costs through seqsubm (seqdist.R:374, 454).
  seqdist_lcs_one_state = quote(TraMineR::seqdist(one, method = "LCS")),
  seqdist_ham_one_state = quote(TraMineR::seqdist(one_full, method = "HAM")),
  seqdist_ham_default = quote(TraMineR::seqdist(tri, method = "HAM")),
  # seqMD resolves each channel's costs with suppressMessages(seqcost(..., check.costs = FALSE))
  # and checks the combined costs itself (seqMD.R:241-243, 266).
  seqmd_trate = quote(TraMineR::seqMD(list(tri, tri), method = "OM", sm = list("TRATE", "TRATE"), what = "diss")),
  # Alternating sequences: every time-varying DHD cost is 0.
  seqcost_trate_time_varying_alternating = quote(TraMineR::seqcost(alt, method = "TRATE", time.varying = TRUE)),
  seqdist_dhd_alternating = quote(TraMineR::seqdist(alt, method = "DHD")),
  seqdist_dhd_trate_alternating = quote(TraMineR::seqdist(alt, method = "DHD", sm = "TRATE")),
  seqdist_dhd_indels_alternating = quote(TraMineR::seqdist(alt, method = "DHD", sm = "INDELS")),
  seqdist_om_trate_refseq_1 = quote(TraMineR::seqdist(tri3, method = "OM", sm = "TRATE", refseq = 1)),
  seqdist_om_trate_refseq_0 = quote(TraMineR::seqdist(tri3, method = "OM", sm = "TRATE", refseq = 0)),
  seqdist_om_trate_sets = quote(TraMineR::seqdist(tri3, method = "OM", sm = "TRATE", refseq = list(1:2, 3))),
  seqdist_om_trate_external = quote(TraMineR::seqdist(tri3, method = "OM", sm = "TRATE", refseq = tri3[2, ])),
  seqdist_om_trate_norm_refseq_1 = quote(TraMineR::seqdist(tri3, method = "OM", sm = "TRATE", norm = "auto",
                                                           refseq = 1)),
  seqdist_om_constant_refseq_1 = quote(TraMineR::seqdist(tri3, method = "OM", sm = "CONSTANT", refseq = 1)),
  seqdist_ham_trate_refseq_1 = quote(TraMineR::seqdist(tri3, method = "HAM", sm = "TRATE", refseq = 1)),
  seqdist_omspell_trate_refseq_1 = quote(TraMineR::seqdist(tri3, method = "OMspell", sm = "TRATE", refseq = 1)),
  seqdist_dhd_refseq_1 = quote(TraMineR::seqdist(tri3, method = "DHD", refseq = 1)),
  seqdist_om_matrix_refseq_1 = quote(TraMineR::seqdist(tri3, method = "OM", sm = tri_sm, indel = 2, refseq = 1)),
  seqdist_om_matrix_sets = quote(TraMineR::seqdist(tri3, method = "OM", sm = tri_sm, indel = 2,
                                                   refseq = list(1:2, 3))),
  seqdist_lcs_one_state_refseq_1 = quote(TraMineR::seqdist(one, method = "LCS", refseq = 1)),
  seqdist_ham_one_state_refseq_1 = quote(TraMineR::seqdist(one_full, method = "HAM", refseq = 1))
)
rows <- list()
for (name in names(cases)) {
  msgs <- list()
  ok <- tryCatch({
    withCallingHandlers(eval(cases[[name]]),
      message = function(m) { msgs[[length(msgs) + 1L]] <<- c("message", conditionMessage(m)); invokeRestart("muffleMessage") },
      warning = function(w) { msgs[[length(msgs) + 1L]] <<- c("warning", conditionMessage(w)); invokeRestart("muffleWarning") })
    TRUE
  }, error = function(e) { msgs[[length(msgs) + 1L]] <<- c("error", conditionMessage(e)); FALSE })
  flagged <- Filter(function(m) grepl("[!", m[2], fixed = TRUE), msgs)
  if (length(flagged) == 0L) {
    rows[[length(rows) + 1L]] <- data.frame(case = name, ok = ok, n = 0L, index = 0L, kind = "", text = "",
                                            stringsAsFactors = FALSE)
  } else {
    for (i in seq_along(flagged)) {
      rows[[length(rows) + 1L]] <- data.frame(case = name, ok = ok, n = length(flagged), index = i,
                                              kind = flagged[[i]][1], text = sub("\n$", "", flagged[[i]][2]),
                                              stringsAsFactors = FALSE)
    }
  }
}
write.csv(do.call(rbind, rows), file.path(out_dir, "generated_cost_messages.csv"), row.names = FALSE, quote = TRUE)
