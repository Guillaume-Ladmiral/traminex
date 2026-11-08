# Live-R goldens for OMslen's link check (TX-AUD-20261011-051).
#
# seqdist checks link before computing: every method refuses a value outside
# none/mean/gmean (seqdist.R:185-186), and OMslen also refuses "none"
# ("'link' must be specified", :229-233), whatever the number of sequences. Rows:
# value (column-major, 17 digits) or R's error.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_omslen_link.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

s3 <- suppressMessages(TraMineR::seqdef(data.frame(a = c("A", "B", "A"), b = c("B", "B", "A"), c = c("A", "A", "B")),
                                        alphabet = c("A", "B")))
same <- s3[c(1, 1), ]
rows <- list()
for (link in c("none", "bogus", "mean", "gmean")) {
  for (data_name in c("one", "identical", "distinct")) {
    data <- switch(data_name, one = s3[1, ], identical = same, distinct = s3)
    result <- tryCatch(suppressWarnings(suppressMessages(TraMineR::seqdist(data, method = "OMslen", sm = "CONSTANT",
                                                                            link = link))), error = function(e) e)
    ok <- !inherits(result, "error")
    rows[[length(rows) + 1L]] <- data.frame(
      case = paste0(link, "_", data_name), ok = ok, error = if (ok) "" else trimws(conditionMessage(result)),
      value = if (ok) paste(format(as.vector(as.matrix(result)), digits = 17), collapse = " ") else "",
      stringsAsFactors = FALSE)
  }
}
write.csv(do.call(rbind, rows), file.path(out_dir, "omslen_link.csv"), row.names = FALSE, quote = TRUE)
