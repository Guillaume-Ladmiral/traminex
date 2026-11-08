# Live-R goldens for NMS and SVRspell with a proximity matrix holding negative entries (TX-REP-20261019-001,
# batch-2 review M-6): the full matrix and the distances to reference sequence 2.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_negative_prox.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

x <- suppressMessages(TraMineR::seqdef(rbind(c("A", "A", "B", "C", "C"), c("A", "B", "B", "C", "A"),
                                             c("C", "C", "B", "A", "A"), c("B", "B", "C", "A", "B"))))
prox <- matrix(c(1, -0.4, 0.2, -0.4, 1, 0.3, 0.2, 0.3, 1), 3, 3)
rows <- list()
for (m in c("NMS", "SVRspell")) {
  full <- suppressMessages(TraMineR::seqdist(x, method = m, prox = prox))
  ref <- suppressMessages(TraMineR::seqdist(x, method = m, prox = prox, refseq = 2))
  rows[[length(rows) + 1L]] <- data.frame(method = m, refseq = "none", value = paste(format(as.vector(full), digits = 15), collapse = " "))
  rows[[length(rows) + 1L]] <- data.frame(method = m, refseq = "2", value = paste(format(as.vector(ref), digits = 15), collapse = " "))
}
write.csv(do.call(rbind, rows), file.path(out_dir, "negative_prox.csv"), row.names = FALSE, quote = TRUE)
