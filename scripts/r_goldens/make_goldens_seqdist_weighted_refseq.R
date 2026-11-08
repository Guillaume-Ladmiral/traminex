# Live-R goldens for seqdist's reference forms on a weighted sequence object
# with the weight-using methods CHI2 and OMstran (repair of audit finding
# TX-AUD-20261010-003).
#
# seqdist accepts refseq for every method (seqdist.R:65-98); CHI2 and OMstran
# keep their own reference handling (seqdist.R:550-600) and use the object's
# weights unless weighted=FALSE.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_seqdist_weighted_refseq.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rows_in <- lapply(c("AABBC", "ABBCC", "CCBAA", "BBCCA"), function(s) strsplit(s, "")[[1]])
data <- as.data.frame(do.call(rbind, rows_in), stringsAsFactors = FALSE)
seqs <- suppressMessages(TraMineR::seqdef(data, alphabet = c("A", "B", "C"), weights = c(1, 2, 0.5, 1.5)))

references <- list(sets = list(1:2, 3:4), index2 = 2, mostfrequent = 0, sequence2 = seqs[2, ])
methods <- list(
  CHI2 = list(method = "CHI2"),
  OMstran = list(method = "OMstran", sm = "CONSTANT", indel = 1, otto = 0.5)
)

rows <- list()
for (method_name in names(methods)) {
  for (reference_name in names(references)) {
    for (weighted in c(TRUE, FALSE)) {
      args <- c(list(seqdata = seqs, refseq = references[[reference_name]], weighted = weighted),
                methods[[method_name]])
      value <- as.matrix(suppressMessages(do.call(TraMineR::seqdist, args)))
      if (ncol(value) == 1L || is.null(dim(value))) value <- matrix(as.vector(value), ncol = 1L)
      rows[[length(rows) + 1L]] <- data.frame(
        case = paste(method_name, reference_name, if (weighted) "weighted" else "unweighted", sep = "_"),
        row = rep(seq_len(nrow(value)), times = ncol(value)),
        col = rep(seq_len(ncol(value)), each = nrow(value)),
        value = sprintf("%.15g", as.vector(value)),
        stringsAsFactors = FALSE
      )
    }
  }
}
write.csv(do.call(rbind, rows), file.path(out_dir, "seqdist_weighted_refseq.csv"),
          row.names = FALSE, quote = TRUE)
