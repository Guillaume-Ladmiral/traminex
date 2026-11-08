# Live-R goldens for the statistic TraMineRextras' seqplot.rf prints (repair of
# audit finding TX-AUD-20261010-041).
#
# seqplot.rf computes its own quality statistic from squared distances,
# R2 = 1 - sum(kmedoid.dist^2) / sum(gmedoid.dist^2) and
# F = (R2 / (k - 1)) / ((1 - R2) / (n - k)), with the global medoid
# disscenter(diss, medoids.index = "first") (seqplot-rf.R:22-24, 139-146), prints
# it with format() and writes it into the subtitle. The sequences below are
# distinct and their MDS scores have no ties, so the rank(ties.method = "random")
# that forms R's groups is deterministic. R's own scores, cmdscale(diss, k = 1),
# are passed as sortv and recorded, so both sides form the same groups whatever
# the sign of the eigenvector (TX-REP-20261017-006). seqplot.rf returns each
# sequence's medoid index, recorded to check that the groups are the same.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_seqplot_rf_stat.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages({library(TraMineR); library(TraMineRextras)})
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14",
          as.character(packageVersion("TraMineRextras")) == "0.6.9")
pdf(NULL)

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

strings <- c("AAABBC", "AABBCC", "ABBBCC", "AAAABB", "BBBCCC", "ABCABC", "AACCBB", "CCBBAA",
             "ABBBBB", "AAACCC", "BCCCCC", "AAAAAC", "BBAACC", "CABBBA")
seqs <- suppressMessages(TraMineR::seqdef(as.data.frame(do.call(rbind, strsplit(strings, ""))),
                                          alphabet = c("A", "B", "C")))
D <- suppressMessages(TraMineR::seqdist(seqs, method = "LCS"))
sv <- cmdscale(D, k = 1)[, 1]
stopifnot(!anyDuplicated(round(sv, 9)))

rows <- list()
# NA: k left at its default, floor(nrow(seqdata) / 10).
for (k in c(1L, 3L, 4L, 14L, NA)) {
  set.seed(1)
  msgs <- character(0)
  call <- if (is.na(k)) quote(TraMineRextras::seqplot.rf(seqs, diss = D, sortv = sv)) else
    quote(TraMineRextras::seqplot.rf(seqs, diss = D, k = k, sortv = sv))
  index <- withCallingHandlers(eval(call),
                               message = function(m) { msgs <<- c(msgs, conditionMessage(m))
                                                       invokeRestart("muffleMessage") })
  r2 <- sub(".*R2: ", "", trimws(grep("median-based-R2", msgs, value = TRUE)))
  fstat <- sub(".*statistic: ", "", trimws(grep("F statistic", msgs, value = TRUE)))
  rows[[length(rows) + 1L]] <- data.frame(k = if (is.na(k)) "default" else as.character(k),
                                          R2 = r2, Fstat = fstat,
                                          global_medoid = TraMineR::disscenter(D, medoids.index = "first"),
                                          kmedoid_index = paste(index, collapse = "|"),
                                          sortv = paste(sprintf("%.17g", sv), collapse = "|"),
                                          stringsAsFactors = FALSE)
}
write.csv(do.call(rbind, rows), file.path(out_dir, "seqplot_rf_stat.csv"), row.names = FALSE, quote = TRUE)
