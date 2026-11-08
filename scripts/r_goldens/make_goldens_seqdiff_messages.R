# Live-R goldens for the conditions seqdiff lets reach its caller (TX-AUD-20261019-059).
#
# seqdiff runs seqdist in every sliding window inside suppressMessages() (seqdiff.R:62), so seqdist's
# msg.warn notices - "with.missing set as FALSE" for data without missing values, "Empty sequence(s)" - never
# reach the caller, while warning() calls would; "chi2_breaks_step" adds seqdist's msg.warn.ign2 notices
# ("'step' has been ignored because 'breaks' is specified", seqdist.R:258; batch-9 review M-2). Each row records one call: the number of messages and warnings
# that reached the caller, their texts, and the per-position statistics (stat, column by column) and
# discrepancies (column by column), so the results are pinned too.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_seqdiff_messages.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

free <- suppressMessages(TraMineR::seqdef(rbind(
  c("A", "A", "B", "B", "C"), c("A", "B", "B", "C", "C"), c("C", "C", "B", "A", "A"),
  c("B", "B", "C", "C", "A"), c("A", "A", "A", "A", "B"), c("C", "C", "C", "B", "B"))))
holes <- suppressMessages(TraMineR::seqdef(rbind(
  c("A", NA, "B", "B", "C"), c("A", "B", "B", "C", "C"), c("C", "C", NA, "A", "A"),
  c("B", "B", "C", "C", "A"), c("A", "A", "A", "A", "B"), c("C", "C", "C", "B", "B"))))
g <- c("m", "f", "m", "f", "m", "f")
cases <- list(
  missing_free = quote(TraMineR::seqdiff(free, group = g, cmprange = c(0, 1), with.missing = TRUE)),
  with_missing_cells = quote(TraMineR::seqdiff(holes, group = g, cmprange = c(0, 1), with.missing = TRUE)),
  chi2_breaks_step = quote(TraMineR::seqdiff(free, group = g, cmprange = c(0, 1),
                                             seqdist.args = list(method = "CHI2", breaks = list(c(1, 2)), step = 2)))
)
rows <- list()
for (name in names(cases)) {
  conditions <- character(0)
  r <- withCallingHandlers(eval(cases[[name]]),
    message = function(m) { conditions <<- c(conditions, paste("message:", conditionMessage(m))); invokeRestart("muffleMessage") },
    warning = function(w) { conditions <<- c(conditions, paste("warning:", conditionMessage(w))); invokeRestart("muffleWarning") })
  rows[[length(rows) + 1L]] <- data.frame(
    case = name, conditions = length(conditions), texts = paste(sub("\n$", "", conditions), collapse = " | "),
    stat = paste(format(as.vector(as.matrix(r$stat)), digits = 15), collapse = " "),
    discrepancy = paste(format(as.vector(as.matrix(r$discrepancy)), digits = 15), collapse = " "),
    stringsAsFactors = FALSE)
}
write.csv(do.call(rbind, rows), file.path(out_dir, "seqdiff_messages.csv"), row.names = FALSE, quote = TRUE)
