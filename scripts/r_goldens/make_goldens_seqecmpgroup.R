# Golden generation: seqecmpgroup -- discriminating-subsequence chi-square
# (seqecmpgroup.R). For each frequent subsequence it cross-tabulates (weighted)
# group x presence and runs chisq.test, returning per-group Freq/Resid, the
# statistic, p.value (+ Bonferroni), reordered by statistic descending.
#
# Groups are derived from the eseq POSITION (seq_along) -- identical between R and
# Python (both order sequences by id first-appearance) -- so the test reconstructs
# them without relying on id factor codes: a 2-group split (exercises Yates'
# correction for the 2x2 table) and a 3-group split (df = 2, no correction).
#
# Requires live R + TraMineR 2.2-14; see memory "live-traminer-r-env".
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
seqecmpgroup <- TraMineR:::seqecmpgroup

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

data(actcal.tse, package = "TraMineR")
ea <- seqecreate(actcal.tse)
fsub <- seqefsub(ea, pmin.support = 0.05, max.k = 2)

dump_cmp <- function(cmp, tag) {
  writeLines(as.character(cmp$subseq), file.path(out_dir, paste0("seqecmpgroup_", tag, "_subseq.txt")))
  writeLines(colnames(cmp$data), file.path(out_dir, paste0("seqecmpgroup_", tag, "_cols.txt")))
  m <- as.matrix(cmp$data)
  con <- file(file.path(out_dir, paste0("seqecmpgroup_", tag, ".csv")), "w")
  cat(paste0("# columns: ", paste(colnames(cmp$data), collapse = ","), "\n"), file = con)
  for (i in seq_len(nrow(m))) {
    cat(paste(sprintf("%.15g", m[i, ]), collapse = ","), "\n", file = con, sep = "")
  }
  close(con)
}

g2 <- seq_along(ea) %% 2
g3 <- seq_along(ea) %% 3
dump_cmp(seqecmpgroup(fsub, group = g2), "g2")
dump_cmp(seqecmpgroup(fsub, group = g3), "g3")
dump_cmp(seqecmpgroup(fsub, group = g2, method = "bonferroni"), "g2_bonf")

cat("SEQECMPGROUP_GOLDENS_DONE\n")
