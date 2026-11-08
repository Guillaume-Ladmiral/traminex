# Live-R goldens for the order of event labels and mined subsequences (repair of
# audit findings TX-AUD-20261011-015 and -017).
#
# seqecreate's dictionary is levels(factor(event)), i.e. R's collation (ICU,
# pinned to the root locale here): transitions to missing ('A>*') sort before
# transitions to void ('A>%'), and accented labels next to their base letter.
# seqefsub's tie order follows the dictionary.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_seqecreate_label_order.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
icuSetCollate(locale = "root")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cells <- data.frame(t1 = c("A", "A", "C", "B"), t2 = c("A", "A", "A", "B"), t3 = c(NA, NA, NA, "A"),
                    t4 = c("B", NA, "C", "A"), stringsAsFactors = FALSE)
seqs <- suppressMessages(TraMineR::seqdef(cells, alphabet = c("A", "B", "C"), right = "DEL"))
events <- suppressMessages(TraMineR::seqecreate(seqs))
subsequences <- suppressMessages(TraMineR::seqefsub(events, pmin.support = 0.25))
labels <- TraMineR::seqecreate(id = 1:8, timestamp = rep(0, 8),
                               event = c("e", "é", "f", "_x", "B", "a", "A>%", "A>*"))
record <- function(case, kind, values) {
  data.frame(case = case, kind = kind, position = seq_along(values), value = as.character(values),
             stringsAsFactors = FALSE)
}
rows <- list(
  record("transitions", "cells", as.vector(t(as.matrix(seqs)))),
  record("transitions", "dictionary", levels(events)),
  record("transitions", "seqefsub", as.character(subsequences$subseq)),
  record("labels", "dictionary", levels(labels))
)
write.csv(do.call(rbind, rows), file.path(out_dir, "seqecreate_label_order.csv"), row.names = FALSE, quote = TRUE)
