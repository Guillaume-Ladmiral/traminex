# Live-R goldens for subsetting a subseqelist (repair of audit finding
# TX-AUD-20261010-038).
#
# R's "[.subseqelist" (subseqelist.R:28-45) returns a subseqelist of the same
# class, subsetting subseq and data rows and keeping labels and bonferroni, so
# seqecmpgroup accepts fs[1:2].
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_subseqelist_subset.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

e2 <- TraMineR::seqecreate(id = c(1, 1, 2, 2, 3, 3, 4, 4), timestamp = c(1, 2, 1, 3, 2, 4, 1, 2),
                           event = c("A", "B", "A", "C", "B", "C", "A", "B"))
fs <- TraMineR::seqefsub(e2, pmin.support = 0.25)
record <- function(case, kind, values) {
  data.frame(case = case, kind = kind, position = seq_along(values),
             value = if (is.character(values)) values else format(values, digits = 15, trim = TRUE),
             stringsAsFactors = FALSE)
}
rows <- list(record("full", "subseq", as.character(fs$subseq)))
subsets <- list(slice_1_2 = 1:2, index_1_3 = c(1, 3), single_2 = 2,
                mask_1_3_5 = c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE))
for (name in names(subsets)) {
  sub <- fs[subsets[[name]]]
  stopifnot(inherits(sub, "subseqelist"))
  rows[[length(rows) + 1L]] <- record(name, "class", class(sub))
  rows[[length(rows) + 1L]] <- record(name, "subseq", as.character(sub$subseq))
  rows[[length(rows) + 1L]] <- record(name, "Support", sub$data$Support)
  rows[[length(rows) + 1L]] <- record(name, "Count", sub$data$Count)
}
comparisons <- list(cmp_slice_1_2 = list(index = 1:2, group = c("g1", "g1", "g2", "g2")),
                    cmp_index_1_3 = list(index = c(1, 3), group = c("g1", "g2", "g1", "g2")))
for (name in names(comparisons)) {
  cmp <- TraMineR::seqecmpgroup(fs[comparisons[[name]]$index], group = comparisons[[name]]$group)
  rows[[length(rows) + 1L]] <- record(name, "subseq", as.character(cmp$subseq))
  for (column in names(cmp$data)) rows[[length(rows) + 1L]] <- record(name, column, cmp$data[[column]])
}
whole <- TraMineR::seqecmpgroup(fs, group = c("g1", "g1", "g2", "g2"))
part <- whole[1:2]
rows[[length(rows) + 1L]] <- record("chisq_slice_1_2", "class", class(part))
rows[[length(rows) + 1L]] <- record("chisq_slice_1_2", "subseq", as.character(part$subseq))
rows[[length(rows) + 1L]] <- record("chisq_slice_1_2", "labels", part$labels)
rows[[length(rows) + 1L]] <- record("chisq_slice_1_2", "bonferroni", unlist(part$bonferroni))
write.csv(do.call(rbind, rows), file.path(out_dir, "subseqelist_subset.csv"), row.names = FALSE, quote = TRUE)
