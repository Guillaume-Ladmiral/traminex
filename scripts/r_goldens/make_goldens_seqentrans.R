# Golden generation: seqentrans (TraMineRextras) -- augments a seqefsub frequent
# subsequence table with ntrans / nevent / Avg.occ.
#
# seqentrans does NO mining: it reads fsubseq$subseq (the event-subsequence label
# strings) and adds, per row,
#   ntrans  = number of "-"-separated transitions in the label,
#   nevent  = ntrans - 1 + number of ","-separated tokens in the label,
#   Avg.occ = Count / sum(seqeweight(eseq))   (only when avg.occ=TRUE).
# We anchor traminex's augmentation to those exact formulas on TraMineR's own
# (structurally varied -- multi-event, multi-transition) subsequence labels.
#
# Requires live R + TraMineR 2.2-14 + TraMineRextras 0.6.9; see memory
# "live-traminer-r-env".
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
suppressMessages(library(TraMineRextras))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
stopifnot(as.character(packageVersion("TraMineRextras")) == "0.6.9")

data(actcal.tse, package = "TraMineR")
eseq <- seqecreate(actcal.tse)
fsub <- seqefsub(eseq, pmin.support = 0.01, max.k = 3)
ft <- seqentrans(fsub, avg.occ = TRUE)
d <- ft$data

wtot <- sum(seqeweight(fsub$eseq))

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Numeric augmentation columns in a fixed order (rows aligned with the labels file).
num <- cbind(d$Count, d$ntrans, d$nevent, d$Avg.occ)
con <- file(file.path(out_dir, "seqentrans.csv"), "w")
cat("# R: TraMineRextras::seqentrans(seqefsub(seqecreate(actcal.tse), pmin.support=0.01, max.k=3), avg.occ=TRUE)\n",
  file = con)
cat(sprintf("# columns: Count,ntrans,nevent,Avg.occ ; wtot=sum(seqeweight(eseq))=%.15g\n", wtot),
  file = con)
for (i in seq_len(nrow(num))) {
  cat(paste(sprintf("%.15g", num[i, ]), collapse = ","), "\n", file = con, sep = "")
}
close(con)

# Subsequence label strings (one per row), for the formula-faithfulness test.
writeLines(as.character(ft$subseq), file.path(out_dir, "seqentrans_subseq.txt"))
writeLines(sprintf("%.15g", wtot), file.path(out_dir, "seqentrans_wtot.txt"))

# Avg.occ after mining with weighted = FALSE (TX-AUD-20261019-016): R divides the unweighted
# Count by the total weight (seqentrans.R:11-14): 3 occurrences in 3 sequences of weight 10 give
# 0.1 where the average number of occurrences by sequence is 1. Weighted mining, where Count is
# weighted too, is recorded beside it.
ue <- seqecreate(id = c(1, 1, 2, 2, 3), timestamp = c(0, 1, 0, 2, 0), event = c("A", "B", "A", "B", "A"))
seqeweight(ue) <- c(10, 10, 10)
unweighted_rows <- list()
for (weighted in c(FALSE, TRUE)) {
  ut <- seqentrans(seqefsub(ue, pmin.support = 0.1, weighted = weighted), avg.occ = TRUE)
  unweighted_rows[[length(unweighted_rows) + 1L]] <- data.frame(
    weighted = weighted, subseq = as.character(ut$subseq), Support = sprintf("%.17g", ut$data$Support),
    Count = sprintf("%.17g", ut$data$Count), Avg.occ = sprintf("%.17g", ut$data$Avg.occ))
}
write.csv(do.call(rbind, unweighted_rows), file.path(out_dir, "seqentrans_unweighted.csv"), row.names = FALSE,
          quote = TRUE)

cat("SEQENTRANS_GOLDENS_DONE\n")
cat("nrows:", nrow(num), " wtot:", wtot, "\n")
