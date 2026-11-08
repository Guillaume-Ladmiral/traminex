# Live-R goldens for CHI2/EUCLID step and breaks (TX-AUD-20261011-048).
#
# step must be a positive integer (seqdist.R:189-190). breaks is documented as a
# "List of ordered pairs of integers" of start and end positions (seqdist.Rd:209-214);
# R does not check them, and a non-integer, unordered or out-of-range pair gives
# whatever seq(start, end) selects. The global check lets step = 0 through
# (is.a.positive.integer(0) is TRUE) and CHI2 then stops inside seq(), unless breaks
# replaces the step (seqdist.R:256-258; batch-6 review M-1). Rows: value
# (column-major, 17 digits) or R's error, for CHI2 and EUCLID.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_chi2_windows.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

s5 <- suppressMessages(TraMineR::seqdef(data.frame(a = c("A", "B", "A", "B"), b = c("B", "B", "A", "A"),
                                                   c = c("A", "A", "B", "B"), d = c("B", "A", "A", "B"),
                                                   e = c("A", "A", "A", "B")), alphabet = c("A", "B")))
windows <- list(step_2_5 = list(step = 2.5), step_2 = list(step = 2), step_2_numeric = list(step = 2.0),
                breaks_2_7_to_5 = list(breaks = list(c(2.7, 5))), breaks_2_to_4 = list(breaks = list(c(2, 4))),
                breaks_2_to_5 = list(breaks = list(c(2, 5))), breaks_3_to_5 = list(breaks = list(c(3, 5))),
                breaks_5_to_3 = list(breaks = list(c(5, 3))), breaks_0_to_5 = list(breaks = list(c(0, 5))),
                breaks_1_to_2_5_and_2_5_to_5 = list(breaks = list(c(1, 2.5), c(2.5, 5))),
                breaks_1_to_2_and_3_to_5 = list(breaks = list(c(1, 2), c(3, 5))),
                step_0 = list(step = 0), step_minus_1 = list(step = -1),
                step_0_overlap = list(step = 0, overlap = TRUE),
                step_0_with_breaks_2_to_4 = list(step = 0, breaks = list(c(2, 4))),
                step_2_5_with_breaks_2_to_4 = list(step = 2.5, breaks = list(c(2, 4))))
rows <- list()
for (method in c("CHI2", "EUCLID")) {
  for (name in names(windows)) {
    result <- tryCatch(suppressWarnings(suppressMessages(do.call(TraMineR::seqdist, c(list(seqdata = s5, method = method),
                                                                                     windows[[name]])))),
                       error = function(e) e)
    ok <- !inherits(result, "error")
    rows[[length(rows) + 1L]] <- data.frame(
      case = paste0(tolower(method), "_", name), ok = ok, error = if (ok) "" else trimws(conditionMessage(result)),
      value = if (ok) paste(format(as.vector(as.matrix(result)), digits = 17), collapse = " ") else "",
      stringsAsFactors = FALSE)
  }
}
write.csv(do.call(rbind, rows), file.path(out_dir, "chi2_windows.csv"), row.names = FALSE, quote = TRUE)
