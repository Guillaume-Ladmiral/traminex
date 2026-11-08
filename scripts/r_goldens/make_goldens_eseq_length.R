# Golden generation: seqelength (+ seqelength<-) on an eseq/seqelist.
#
# seqelength(eseq) returns the per-sequence *stored* length (getObsTime), which
# defaults to -1 for TSE-created sequences and is set by the replacement form
# seqelength(eseq)<-. A set length appears in as.character as a trailing gap to
# the end of observation, e.g. "(A)-2-(B)-3" for events at 0,2 with length 5.
#
# Requires live R + TraMineR 2.2-14; see memory "live-traminer-r-env".
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

tse <- data.frame(
  id        = c("i1", "i1", "i2", "i2", "i3"),
  timestamp = c(0,    2,    1,    1,    0),
  event     = c("A",  "B",  "A",  "C",  "B"),
  stringsAsFactors = FALSE
)
e <- seqecreate(tse)

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

writeLines(sprintf("%.15g", seqelength(e)), file.path(out_dir, "eseq_len_default.txt"))

seqelength(e) <- c(5, 4, 2)
writeLines(sprintf("%.15g", seqelength(e)), file.path(out_dir, "eseq_len_set.txt"))
writeLines(as.character(e), file.path(out_dir, "eseq_len_aschar.txt"))

cat("ESEQ_LENGTH_GOLDENS_DONE\n")
print(as.character(e))
