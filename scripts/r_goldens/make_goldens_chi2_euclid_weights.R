# Live-R goldens for CHI2 and EUCLID with sequence weights (repair of audit
# findings TX-AUD-20261011-008 and TX-AUD-20261010-067).
#
# CHI2 weights the state distributions by the sequence weights; EUCLID never
# does: seqdist-CHI2.R:7 sets weighted <- FALSE for EUCLID. All-zero weights give
# an all-zero CHI2 matrix in R.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_chi2_euclid_weights.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

strings <- c("AABBC", "ABBCC", "CCBAA", "BBCCA")
data <- as.data.frame(do.call(rbind, lapply(strings, function(s) strsplit(s, "")[[1]])), stringsAsFactors = FALSE)
weight_sets <- list(unit = c(1, 1, 1, 1), uneven = c(1, 5, 0.2, 3), one_zero = c(0, 1, 1, 1), all_zero = c(0, 0, 0, 0))

rows <- list()
for (weights_name in names(weight_sets)) {
  seqs <- suppressMessages(TraMineR::seqdef(data, alphabet = c("A", "B", "C"), weights = weight_sets[[weights_name]]))
  for (method in c("CHI2", "EUCLID")) {
    value <- as.matrix(suppressMessages(TraMineR::seqdist(seqs, method = method, weighted = TRUE)))
    rows[[length(rows) + 1L]] <- data.frame(
      case = paste(method, weights_name, sep = "_"),
      row = rep(seq_len(nrow(value)), times = ncol(value)),
      col = rep(seq_len(ncol(value)), each = nrow(value)),
      value = sprintf("%.15g", as.vector(value)),
      stringsAsFactors = FALSE
    )
  }
}
write.csv(do.call(rbind, rows), file.path(out_dir, "chi2_euclid_weights.csv"), row.names = FALSE, quote = TRUE)
