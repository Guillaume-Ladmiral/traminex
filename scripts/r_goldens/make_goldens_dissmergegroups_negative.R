# Live-R goldens for dissmergegroups when no merge is found (record of audit
# finding TX-AUD-20261011-016).
#
# dissmergegroups starts each pass with merge.ij <- NULL and diff <- the reference
# quality and records a pair only when its loss is smaller, that is, when the
# merged quality is positive (dissmergegroups.R:127-168). When no pair qualifies,
# diff stays at the reference quality, so the stop test diff > crit * quality.ref
# fails when that quality is 0 or below (crit below 1), and the loop takes the
# merge branch with nothing to merge, printing "Merging groups and" without end.
# traminex stops when no pair qualifies. A negative reference alone does not
# loop: a merge of positive quality still qualifies (batch-4 review M-1).
# The generator runs R under setTimeLimit(elapsed = 5) and records the outcome,
# whether it printed empty "Merging groups and" lines and the reference ASW for
# three groupings: "negative" (negative reference, no merge of positive quality),
# "negative_merging" (negative reference: four hierarchical clusters, each split
# at random) and "control" (hierarchical clusters, positive quality).
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_dissmergegroups_negative.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

data(mvad)
cells <- mvad[1:60, 17:41]
seqs <- suppressMessages(TraMineR::seqdef(cells, alphabet = sort(unique(as.vector(as.matrix(cells))))))
D <- suppressMessages(TraMineR::seqdist(seqs, method = "HAM"))
set.seed(3)
negative <- sample(1:6, 60, replace = TRUE)
control <- cutree(hclust(as.dist(D), method = "ward.D"), k = 8)
set.seed(2)
negative_merging <- 2L * cutree(hclust(as.dist(D), method = "ward.D"), k = 4) - sample(0:1, 60, replace = TRUE)

run <- function(group) {
  captured <- character(0)
  con <- textConnection("captured", "w", local = TRUE)
  sink(con)
  outcome <- tryCatch({
    setTimeLimit(elapsed = 5, transient = TRUE)
    result <- TraMineR::dissmergegroups(D, group, silent = FALSE)
    setTimeLimit(elapsed = Inf)
    list(status = "returned", groups = as.integer(result))
  }, error = function(e) {
    setTimeLimit(elapsed = Inf)
    list(status = paste("error:", conditionMessage(e)), groups = integer(0))
  }, finally = {
    sink()
    close(con)
  })
  # How many lines a timed-out run prints depends on the machine, so only a
  # returned run records its count.
  outcome$merging_lines <- if (outcome$status == "returned") sum(grepl("^Merging groups", captured)) else NA
  outcome$empty_merge_seen <- any(grepl("^Merging groups +and *$", captured))
  outcome
}
rows <- list()
for (case in c("negative", "negative_merging", "control")) {
  outcome <- run(get(case))
  reference <- WeightedCluster::wcClusterQuality(D, get(case))$stats[["ASW"]]
  stopifnot((reference < 0) == (case != "control"))
  rows[[length(rows) + 1L]] <- data.frame(case = case, status = outcome$status,
                                          merging_lines = outcome$merging_lines,
                                          empty_merge_seen = outcome$empty_merge_seen,
                                          reference_asw = format(reference, digits = 17),
                                          input = paste(get(case), collapse = "|"),
                                          groups = paste(outcome$groups, collapse = "|"),
                                          stringsAsFactors = FALSE)
}
write.csv(do.call(rbind, rows), file.path(out_dir, "dissmergegroups_negative.csv"), row.names = FALSE, quote = TRUE)
write.table(format(D, digits = 17), file.path(out_dir, "dissmergegroups_negative_diss.csv"), sep = ",",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
