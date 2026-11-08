# Golden generation: seqecontain -- per-sequence event-set containment.
#
# seqecontain(eseq, event.list, unknown.exclude) returns a logical per sequence:
# with unknown.exclude=FALSE (default) TRUE iff the sequence contains ANY event
# in event.list; with unknown.exclude=TRUE TRUE iff ALL of the sequence's events
# are in event.list. We dump as 1/0 for a few event lists / modes.
#
# Requires live R + TraMineR 2.2-14; see memory "live-traminer-r-env".
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

tse <- data.frame(
  id        = c("i1", "i1", "i1", "i2", "i2", "i3", "i3", "i4", "i4", "i5"),
  timestamp = c(0,    1,    2,    0,    1,    0,    1,    0,    1,    0),
  event     = c("A",  "B",  "C",  "A",  "B",  "A",  "A",  "C",  "C",  "B"),
  stringsAsFactors = FALSE
)
e <- seqecreate(tse)

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

writeLines(as.character(as.integer(seqecontain(e, c("A")))),
           file.path(out_dir, "eseq_contain_anyA.txt"))
writeLines(as.character(as.integer(seqecontain(e, c("A", "C")))),
           file.path(out_dir, "eseq_contain_anyAC.txt"))
writeLines(as.character(as.integer(seqecontain(e, c("A", "B"), unknown.exclude = TRUE))),
           file.path(out_dir, "eseq_contain_allAB.txt"))

cat("ESEQ_CONTAIN_GOLDENS_DONE\n")
cat("anyA:", as.integer(seqecontain(e, c("A"))), "\n")
cat("allAB:", as.integer(seqecontain(e, c("A", "B"), unknown.exclude = TRUE)), "\n")
