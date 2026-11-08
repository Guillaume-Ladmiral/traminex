# Golden generation: the faithful event-sequence object (eseq / seqelist).
#
# Pins as.character(seqecreate(tse)) -- the TraMineR event-label strings -- and
# levels(eseq) (the event dictionary) for a fixed TSE that exercises: a leading
# gap, simultaneous events (same timestamp -> comma group, sorted by event code),
# a repeated event, and a singleton.  These are the C-level sprint() / dictionary
# semantics that the RP5a object must reproduce.
#
# Requires live R + TraMineR 2.2-14; see memory "live-traminer-r-env".
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

tse <- data.frame(
  id        = c("i1", "i1", "i1", "i2", "i2", "i3", "i3", "i4", "i4", "i4", "i5"),
  timestamp = c(0,    2,    2,    1,    3,    0,    0,    0,    1,    2,    4),
  event     = c("A",  "B",  "C",  "A",  "B",  "B",  "A",  "C",  "C",  "A",  "B"),
  stringsAsFactors = FALSE
)
e <- seqecreate(tse)

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

writeLines(as.character(e), file.path(out_dir, "eseq_aschar.txt"))
writeLines(levels(e), file.path(out_dir, "eseq_levels.txt"))

cat("EVENTS_OBJ_GOLDENS_DONE\n")
print(as.character(e))
