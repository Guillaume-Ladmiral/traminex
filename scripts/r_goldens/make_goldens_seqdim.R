# Live-R goldens for seqdim (repair of audit finding TX-AUD-20261011-071).
#
# seqdim runs seqfcheck on as.matrix(seqdata) (seqdim.R:5-15, seqfcheck.R:3-23):
# numeric data is never decomposed; one-column character data is compressed when
# any element holds "-" (else ":") and is then split with seqdecomp; data with
# several columns is never decomposed. The dimensions are then c(1, length) for a
# vector and c(nrow, ncol) otherwise. Each row records R's seqdim for one input.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_seqdim.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cases <- list(
  plain_then_compressed = c("AAB", "A-B-C"),
  na_then_compressed = c(NA, "A-B-C"),
  single_compressed = "A-B-C",
  matrix_with_dash = matrix(c("A-B", "C", "D", "E"), 2),
  one_column_frame = data.frame(x = c("A-B", "C-D-E"), stringsAsFactors = FALSE),
  factor_column = data.frame(x = factor(c("A-B", "C-D"))),
  numeric_matrix = matrix(1:6, 2),
  all_na_vector = c(NA, NA, NA),
  colon_compressed = c("A:B", "C:D:E"),
  plain_vector = c("A", "B", "C"),
  numeric_vector_negative = c(-1, 2, 3),
  plain_scalar = "ABC",
  # c() makes these character vectors: "1", "A-B" and "-1", "A" (batch-5 review M-5)
  mixed_number_text = c(1, "A-B"),
  mixed_negative_number_text = c(-1, "A")
)
rows <- do.call(rbind, lapply(names(cases), function(case) {
  dims <- TraMineR::seqdim(cases[[case]])
  data.frame(case = case, rows = dims[["Rows"]], columns = dims[["Columns"]],
             format = TraMineR::seqfcheck(cases[[case]]), stringsAsFactors = FALSE)
}))
write.csv(rows, file.path(out_dir, "seqdim_cases.csv"), row.names = FALSE, quote = TRUE)
