# Golden generation: seqeage -- the event-age table (TraMineR:::seqeage, an
# unexported wrapper over the C routine tmreventinseq, tmrsequence.cpp:430-484).
# For one event code it returns an (n_sequences x max_occurrences) matrix whose
# cell [i,k] is the absolute time (cumulative inter-event gap) of the (k+1)-th
# occurrence of that event in sequence i, or -1 when it occurs fewer than k+1
# times; column width = max(1, max occurrences across sequences).
#
# Fixture (inline TSE) pins: repeated events (A twice in s1/s2/s3), simultaneous
# events (B,C tied at t=2 in s1), the SAME event tied at one time (A,A at t=1 in
# s3), differing per-sequence event counts (4/3/4/2), and a declared-but-absent
# level D (never occurs -> a single column of -1). It also pins the C routine's
# behaviour of using ONLY the first element of a multi-event `event.list`.
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

tse <- data.frame(
  id    = c("s1", "s1", "s1", "s1", "s2", "s2", "s2", "s3", "s3", "s3", "s3", "s4", "s4"),
  time  = c(0,    2,    2,    5,    0,    1,    3,    0,    1,    1,    4,    0,    2),
  event = factor(
    c("A", "B", "C", "A", "A", "A", "C", "B", "A", "A", "C", "A", "B"),
    levels = c("A", "B", "C", "D")
  ),
  stringsAsFactors = FALSE
)

e <- seqecreate(tse)
seqeage <- TraMineR:::seqeage

cat("levels:", paste(levels(e), collapse = ","), "\n")

dump_mat <- function(m, tag) {
  m <- as.matrix(m)
  con <- file(file.path(out_dir, paste0("seqeage_", tag, ".csv")), "w")
  cat(sprintf("# %d x %d event-age matrix\n", nrow(m), ncol(m)), file = con)
  for (i in seq_len(nrow(m))) {
    cat(paste(sprintf("%.15g", m[i, ]), collapse = ","), "\n", file = con, sep = "")
  }
  close(con)
  cat("== seqeage_", tag, " ==\n", sep = "")
  print(m)
}

dump_mat(seqeage(e, 1L), "A")   # repeated, with absent (s4 has one A)
dump_mat(seqeage(e, 2L), "B")   # single occurrence each, s2 has none
dump_mat(seqeage(e, 3L), "C")   # single occurrence each, s4 has none
dump_mat(seqeage(e, 4L), "D")   # declared level that never occurs -> all -1
dump_mat(seqeage(e, c(1L, 2L)), "list12")  # C uses event.list[1] only -> == A

cat("SEQEAGE_GOLDENS_DONE\n")
