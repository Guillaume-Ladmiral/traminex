# Live-R goldens for the level order of text groups in the discrepancy and tree
# functions (repair of audit finding TX-AUD-20261010-039).
#
# R orders a character group with factor(), i.e. by the session's collation
# (ICU, pinned to the root locale here): dissassoc's group table
# (dissassocweighted.R:294), disstree's split levels (disstree.R:404),
# seqrep.grp's groups and dissmergegroups' initial group numbers
# (dissmergegroups.R:105) all follow it. So do the grouped plots: seqedplot
# (seqedplot.R:53) and seqplot.tentrop (seqplot-tentrop.R:10) call factor(), and
# seqplot, behind seqHtplot, builds its panels from the internal group()
# (seqplot.R:163-184), recorded here because no exported function returns it.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_factor_level_order.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages({library(TraMineR); library(TraMineRextras)})
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14",
          as.character(packageVersion("TraMineRextras")) == "0.6.9")
icuSetCollate(locale = "root")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

strings <- c("ABCAB", "BBCCA", "CABBA", "AACCB", "BCABC", "CCABA",
             "ABBBC", "BACCA", "CBAAB", "AABBC", "BCCAA", "CABCB")
data <- as.data.frame(do.call(rbind, lapply(strings, function(s) strsplit(s, "")[[1]])), stringsAsFactors = FALSE)
seqs <- suppressMessages(TraMineR::seqdef(data, alphabet = c("A", "B", "C")))
diss <- suppressMessages(TraMineR::seqdist(seqs, method = "OM", sm = "CONSTANT", indel = 1))
group <- rep(c("a", "B", "_c", "é"), 3)
labels <- c("b", "a", "B", "_c", "é", "e", "A>%", "A>*", "10", "9", "faible", "élevé")

record <- function(case, kind, values) {
  data.frame(case = case, kind = kind, position = seq_along(values), value = as.character(values),
             stringsAsFactors = FALSE)
}
rows <- list(
  record("labels", "levels", levels(factor(labels))),
  record("labels", "seqplot_group_levels", levels(TraMineR:::group(labels))),
  record("dissassoc", "group_rows", rownames(TraMineR::dissassoc(diss, group, R = 1)$groups)),
  record("disstree", "rules",
         TraMineR::disstree.get.rules(TraMineR::disstree(diss ~ g, data = data.frame(g = group, stringsAsFactors = FALSE),
                                                          R = 1, min.size = 0.1, max.depth = 2))),
  record("seqrep_grp", "groups",
         names(suppressMessages(TraMineRextras::seqrep.grp(seqs, group = group, diss = diss, ret = "stat")))),
  record("dissmergegroups_no_merge", "group",
         suppressWarnings(TraMineR::dissmergegroups(diss, group, min.group = 4, silent = TRUE))),
  record("dissmergegroups_merged", "group", TraMineR::dissmergegroups(diss, group, min.group = 2, silent = TRUE))
)
write.csv(do.call(rbind, rows), file.path(out_dir, "factor_level_order.csv"), row.names = FALSE, quote = TRUE)
write.csv(diss, file.path(out_dir, "factor_level_order_diss.csv"), row.names = FALSE)
