# Live-R goldens for seqdist's refseq forms (TX-AUD-20261011-044).
#
# R takes refseq as one index (0 for the most frequent sequence), a list of two
# index sets, or a one-sequence state sequence object; a vector of several indices
# is none of these and stops ("'refseq' invalid type or unknown value",
# seqdist.R:67-99), as do an index above the number of sequences, a negative or
# fractional index, an index set outside the data, and a character index (batch-6
# review M-4). Rows: value (column-major, 17 digits) or R's error.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_refseq_forms.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

s3 <- suppressMessages(TraMineR::seqdef(data.frame(a = c("A", "B", "A"), b = c("B", "B", "A"), c = c("A", "A", "B")),
                                        alphabet = c("A", "B")))
cases <- list(
  vector_1_2 = quote(TraMineR::seqdist(s3, method = "LCS", refseq = c(1, 2))),
  scalar_out_of_range = quote(TraMineR::seqdist(s3, method = "LCS", refseq = 9)),
  scalar_negative = quote(TraMineR::seqdist(s3, method = "LCS", refseq = -1)),
  scalar_fraction = quote(TraMineR::seqdist(s3, method = "LCS", refseq = 1.5)),
  sets_out_of_range = quote(TraMineR::seqdist(s3, method = "LCS", refseq = list(1, 9))),
  sets_singletons = quote(TraMineR::seqdist(s3, method = "LCS", refseq = list(1, 2))),
  sets_pair_and_one = quote(TraMineR::seqdist(s3, method = "LCS", refseq = list(1:2, 3))),
  scalar_two = quote(TraMineR::seqdist(s3, method = "LCS", refseq = 2)),
  character_1 = quote(TraMineR::seqdist(s3, method = "LCS", refseq = "1"))
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
write.csv(do.call(rbind, rows), file.path(out_dir, "refseq_forms.csv"), row.names = FALSE, quote = TRUE)
