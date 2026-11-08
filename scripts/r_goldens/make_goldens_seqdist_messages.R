# Live-R goldens for seqdist's empty-sequence and with.missing messages (TX-AUD-20261011-042).
#
# seqdist warns "Empty sequence(s) <indices> may cause inconsistent results!" for
# sequences with no non-void element (seqdist.R:101-109; msg.warn.sempty,
# TraMineR-logging_helpers.R:89-91), except for OMloc, which stops; and it resets
# with.missing to FALSE with a message when there is nothing missing (:114-117).
# R's msg.warn is message(" [!!] ...") (logging_helpers.R:20-25). The rows record
# each call's "[!" messages in order, or its error. The refseq cases check each
# message once per call on every reference route (batch-6 review I-1).
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_seqdist_messages.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

emp <- suppressMessages(TraMineR::seqdef(data.frame(a = c("A", "B", NA), b = c("B", "B", NA)), alphabet = c("A", "B")))
emp2 <- suppressMessages(TraMineR::seqdef(data.frame(a = c(NA, "B", NA, "A"), b = c(NA, "B", NA, "A")), alphabet = c("A", "B")))
s3 <- suppressMessages(TraMineR::seqdef(data.frame(a = c("A", "B", "A"), b = c("B", "B", "A"), c = c("A", "A", "B")),
                                        alphabet = c("A", "B")))
cases <- list(
  empty_lcs = quote(TraMineR::seqdist(emp, method = "LCS")),
  empty_om = quote(TraMineR::seqdist(emp, method = "OM", sm = "CONSTANT")),
  empty_ham = quote(TraMineR::seqdist(emp, method = "HAM")),
  empty_omloc = quote(TraMineR::seqdist(emp, method = "OMloc", sm = "CONSTANT")),
  two_empty_lcs = quote(TraMineR::seqdist(emp2, method = "LCS")),
  with_missing_reset_lcs = quote(TraMineR::seqdist(s3, method = "LCS", with.missing = TRUE)),
  with_missing_reset_om = quote(TraMineR::seqdist(s3, method = "OM", sm = "CONSTANT", with.missing = TRUE)),
  empty_lcs_refseq_1 = quote(TraMineR::seqdist(emp, method = "LCS", refseq = 1)),
  empty_lcs_refseq_0 = quote(TraMineR::seqdist(emp, method = "LCS", refseq = 0)),
  empty_lcs_sets = quote(TraMineR::seqdist(emp, method = "LCS", refseq = list(1:2, 3))),
  empty_lcs_external = quote(TraMineR::seqdist(emp, method = "LCS", refseq = emp[1, ])),
  with_missing_reset_lcs_refseq_1 = quote(TraMineR::seqdist(s3, method = "LCS", with.missing = TRUE,
                                                            refseq = 1)),
  with_missing_reset_lcs_sets = quote(TraMineR::seqdist(s3, method = "LCS", with.missing = TRUE,
                                                        refseq = list(1:2, 3)))
)
rows <- list()
for (name in names(cases)) {
  msgs <- list()
  err <- ""
  ok <- tryCatch({
    withCallingHandlers(eval(cases[[name]]),
      message = function(m) { msgs[[length(msgs) + 1L]] <<- c("message", conditionMessage(m)); invokeRestart("muffleMessage") },
      warning = function(w) { msgs[[length(msgs) + 1L]] <<- c("warning", conditionMessage(w)); invokeRestart("muffleWarning") })
    TRUE
  }, error = function(e) { err <<- conditionMessage(e); FALSE })
  flagged <- Filter(function(m) grepl("[!", m[2], fixed = TRUE), msgs)
  if (length(flagged) == 0L) flagged <- list(c("", ""))
  for (i in seq_along(flagged)) {
    rows[[length(rows) + 1L]] <- data.frame(case = name, ok = ok, error = err, index = i, kind = flagged[[i]][1],
                                            text = sub("\n$", "", flagged[[i]][2]), stringsAsFactors = FALSE)
  }
}
write.csv(do.call(rbind, rows), file.path(out_dir, "seqdist_messages.csv"), row.names = FALSE, quote = TRUE)
