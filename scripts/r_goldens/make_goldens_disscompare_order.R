# Live-R goldens for dissCompare's group order (repair of audit finding
# TX-AUD-20261011-022).
#
# seqCompare makes the larger group the first (seqCompare.R:143-157) because it
# draws the first group's sample from 1:max(n); dissCompare, its port to a
# dissimilarity matrix, keeps the first factor level first (dissCompare.R:
# 117-131) but draws the same way, so with the smaller group first its "first"
# sample includes cases of the second group, and the sampled statistic depends on
# the level order. The generator records R's sampled LRT for seeds 1..20 with the
# smaller group first and with the larger group first (the same cases relabelled).
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_disscompare_order.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages({library(TraMineR); library(TraMineRextras)})
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14",
          as.character(packageVersion("TraMineRextras")) == "0.6.9")
icuSetCollate(locale = "root")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

data(mvad)
cells <- mvad[1:90, 17:41]
seqs <- suppressMessages(TraMineR::seqdef(cells, alphabet = sort(unique(as.vector(as.matrix(cells))))))
D <- suppressMessages(TraMineR::seqdist(seqs, method = "LCS"))
smaller_first <- ifelse(seq_len(90) %% 3 == 0, "a", "b")  # 30 "a" (first level), 60 "b"
larger_first <- ifelse(smaller_first == "a", "z", "b")    # the same split, the larger "b" first
rows <- list()
for (order in c("smaller_first", "larger_first")) {
  group <- get(order)
  for (seed in 1:20) {
    result <- suppressMessages(TraMineRextras::dissCompare(D, group, s = 20, seed = seed, stat = "LRT"))
    rows[[length(rows) + 1L]] <- data.frame(order = order, seed = seed, LRT = sprintf("%.17g", result[["LRT"]]),
                                            stringsAsFactors = FALSE)
  }
}
write.csv(do.call(rbind, rows), file.path(out_dir, "disscompare_order.csv"), row.names = FALSE, quote = TRUE)
write.table(format(D, digits = 17), file.path(out_dir, "disscompare_order_diss.csv"), sep = ",",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
