# Golden generation: seqeid / seqeweight (+ seqeweight<-) on an eseq/seqelist.
#
# seqeid(eseq) returns the per-sequence integer id (the id factor code);
# seqeweight(eseq) returns the per-sequence weight (default 1), and the
# replacement form seqeweight(eseq)<- sets them. We use character ids whose
# factor codes differ from data-appearance order so seqeid is non-trivial.
#
# Requires live R + TraMineR 2.2-14; see memory "live-traminer-r-env".
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

tse <- data.frame(
  id        = c("z", "z", "a", "a", "m"),
  timestamp = c(0,   2,   1,   3,   0),
  event     = c("A", "B", "A", "C", "B"),
  stringsAsFactors = FALSE
)
e <- seqecreate(tse)

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

writeLines(sprintf("%.15g", seqeid(e)), file.path(out_dir, "eseq_seqeid.txt"))
writeLines(sprintf("%.15g", seqeweight(e)), file.path(out_dir, "eseq_weight_default.txt"))

seqeweight(e) <- c(1.5, 2.0, 0.5)
writeLines(sprintf("%.15g", seqeweight(e)), file.path(out_dir, "eseq_weight_set.txt"))

# Integer ids: R does `id <- as.integer(id)` directly (no factor()), so seqeid
# keeps the LITERAL integer value, not a sorted-unique rank.
ei <- seqecreate(id = c(2, 10, 1, 20, 3), timestamp = c(0, 0, 0, 0, 0),
                 event = c("A", "B", "A", "B", "A"))
writeLines(sprintf("%.15g", seqeid(ei)), file.path(out_dir, "eseq_seqeid_int.txt"))

# Factor ids with a declared, non-alphabetical level order: as.integer(factor)
# follows the declared levels (z=1, a=2, m=3), NOT a sorted rank.
ef <- seqecreate(id = factor(c("z", "a", "m"), levels = c("z", "a", "m")),
                 timestamp = c(0, 1, 0), event = c("A", "B", "A"))
writeLines(sprintf("%.15g", seqeid(ef)), file.path(out_dir, "eseq_seqeid_factor.txt"))

cat("ESEQ_ACCESSORS_GOLDENS_DONE\n")
cat("seqeid:", seqeid(e), "\n")
