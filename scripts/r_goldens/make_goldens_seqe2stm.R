# Golden generation: seqe2stm (TraMineRextras) -- event -> state-machine table.
#
# seqe2stm(events, dropMatrix, dropList, firstState) builds the deterministic
# automaton whose rows are the reachable "sets of events that have occurred"
# (dot-joined names, or firstState when empty) and whose columns are the events;
# cell [s, e] is the state reached from state s when event e occurs. dropList /
# dropMatrix let an event clear previously-occurred events.
#
# Requires live R + TraMineRextras 0.6.9; see memory "live-traminer-r-env".
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
suppressMessages(library(TraMineRextras))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
stopifnot(as.character(packageVersion("TraMineRextras")) == "0.6.9")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dump_stm <- function(stm, fname, call) {
  con <- file(file.path(out_dir, fname), "w")
  cat(sprintf("# R: %s\n", call), file = con)
  cat(paste(colnames(stm), collapse = "\t"), "\n", file = con, sep = "")
  for (i in seq_len(nrow(stm))) {
    cat(rownames(stm)[i], "\t", paste(stm[i, ], collapse = "\t"), "\n",
        file = con, sep = "")
  }
  close(con)
}

# Case 1: three events, no drop rule.
dump_stm(seqe2stm(c("A", "B", "C")),
         "seqe2stm_abc.tsv", "seqe2stm(c('A','B','C'))")

# Case 2: two events with a drop rule (event B clears A).
dump_stm(seqe2stm(c("A", "B"), dropList = list(B = "A")),
         "seqe2stm_drop.tsv", "seqe2stm(c('A','B'), dropList=list(B='A'))")

# Case 3: custom firstState.
dump_stm(seqe2stm(c("X", "Y"), firstState = "start"),
         "seqe2stm_first.tsv", "seqe2stm(c('X','Y'), firstState='start')")

cat("SEQE2STM_GOLDENS_DONE\n")
