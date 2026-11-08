# Live-R goldens for when seqstatl treats the states as numbers (repair of audit
# finding TX-AUD-20261010-050).
#
# seqstatl converts the distinct values with as.numeric and, when every one
# converts, returns them as numbers in numeric order (seqstatl.R:21-30). R's
# reader rejects digit-group underscores ("1_1" stays text) and accepts hex
# ("0x10" is 16), surrounding spaces, exponents and "Inf".
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_seqstatl_r_numbers.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
icuSetCollate(locale = "root")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cases <- list(
  underscores = matrix(c("1_1", "1_2", "2_1", "1_1"), 2),
  hex = matrix(c("0x10", "2", "2", "0x10"), 2),
  spellings = matrix(c(" 2", "1e2", "0x1A", "Inf", "3", " 2"), 2)
)
rows <- list()
for (name in names(cases)) {
  result <- TraMineR::seqstatl(cases[[name]])
  rows[[length(rows) + 1L]] <- data.frame(case = name, position = seq_along(result),
                                          kind = if (is.numeric(result)) "numeric" else "character",
                                          value = as.character(result), stringsAsFactors = FALSE)
}
write.csv(do.call(rbind, rows), file.path(out_dir, "seqstatl_r_numbers.csv"), row.names = FALSE, quote = TRUE)
