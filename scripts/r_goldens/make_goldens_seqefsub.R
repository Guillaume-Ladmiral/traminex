# Golden generation: seqefsub -- the faithful prefix-tree frequent-subsequence
# miner (TraMineR C_tmrfindsubsequences). Pins the subsequence label set/order +
# Support/Count for several (pmin.support, max.k) configs on actcal.tse, and a
# small synthetic TSE under a non-contiguous max.gap/window constraint.
#
# We also export actcal.tse itself (id,time,event) so the Python test can rebuild
# the identical eseq -- the event column is restored as a Categorical with the
# factor level order pinned in eseq_actcal_levels.txt (see make_goldens_events_factor.R).
#
# Requires live R + TraMineR 2.2-14; see memory "live-traminer-r-env".
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

data(actcal.tse, package = "TraMineR")
# Export the raw TSE for Python reconstruction (event level order in the levels golden).
write.csv(
  data.frame(
    id = actcal.tse$id, time = actcal.tse$time,
    event = as.character(actcal.tse$event), stringsAsFactors = FALSE
  ),
  file.path(out_dir, "actcal_tse.csv"), row.names = FALSE
)

e <- seqecreate(actcal.tse)

dump_fsub <- function(fsub, tag) {
  d <- fsub$data
  writeLines(as.character(fsub$subseq), file.path(out_dir, paste0("seqefsub_", tag, "_subseq.txt")))
  con <- file(file.path(out_dir, paste0("seqefsub_", tag, ".csv")), "w")
  cat("# columns: Support,Count\n", file = con)
  for (i in seq_len(nrow(d))) {
    cat(sprintf("%.15g,%.15g", d$Support[i], d$Count[i]), "\n", file = con, sep = "")
  }
  close(con)
}

dump_fsub(seqefsub(e, pmin.support = 0.05, max.k = 2), "actcal_p05_k2")
dump_fsub(seqefsub(e, pmin.support = 0.02, max.k = 3), "actcal_p02_k3")

# Small synthetic TSE: simultaneous events + a gap, mined under a max.gap window.
tse <- data.frame(
  id    = c("s1", "s1", "s1", "s1", "s2", "s2", "s2", "s3", "s3", "s3"),
  time  = c(0,    1,    1,    3,    0,    2,    4,    0,    1,    2),
  event = c("A",  "B",  "C",  "A",  "A",  "B",  "C",  "A",  "B",  "C"),
  stringsAsFactors = FALSE
)
es <- seqecreate(tse)
dump_fsub(seqefsub(es, min.support = 1, max.k = 3), "syn_k3")
dump_fsub(seqefsub(es, min.support = 1, max.k = 3,
                   constraint = seqeconstraint(max.gap = 2)), "syn_maxgap2")
dump_fsub(seqefsub(es, min.support = 1, max.k = 3,
                   constraint = seqeconstraint(count.method = "CDIST_O")), "syn_cdisto")

# str.subseq: specific (not mined) subsequences. No mining; Support/Count for the
# given labels in the given order. Default (COBJ) => Count == COBJ presence count;
# CDIST_O => Count = occurrence count, distinct from the COBJ-based Support.
dump_fsub(seqefsub(es, str.subseq = c("(A)", "(A)-(B)", "(B,C)", "(C)-(A)")), "syn_strsub")
dump_fsub(seqefsub(e, str.subseq = c("(FullTime)", "(NoActivity)-(Start)", "(Start,LowPartTime)")),
          "actcal_strsub")
dump_fsub(seqefsub(e, str.subseq = c("(NoActivity)", "(NoActivity)-(Start)"),
                   constraint = seqeconstraint(count.method = "CDIST_O")), "actcal_strsub_cdisto")

cat("SEQEFSUB_GOLDENS_DONE\n")
