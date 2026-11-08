# Live-R goldens for POSIX character classes in seqpm patterns (repair of
# audit finding TX-AUD-20261010-054).
#
# seqpm greps its pattern in the seqconc strings (seqpm.R:9-11) with R's
# extended regular expressions, in which bracket classes such as [[:upper:]]
# work. Character classes follow LC_CTYPE, pinned to C.UTF-8 here.
#
# Inside a bracket a backslash is an ordinary character in R's extended regular
# expressions, so [\]A] is the bracket [\] followed by "A]" (TX-AUD-20261010-054,
# batch-3 review M1).
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_seqpm_posix.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

build <- function(strings, alphabet) {
  suppressMessages(TraMineR::seqdef(as.data.frame(do.call(rbind, strsplit(strings, ""))), alphabet = alphabet))
}
plain <- build(c("ABCCAB", "ABAACC", "CCCCCC"), c("A", "B", "C"))
accents <- build(c("AéB", "BBé", "BAB"), c("A", "B", "é"))
brackets <- build(c("A]1b", "\\Ac2", "BBBx"), c("1", "2", "A", "B", "\\", "]", "b", "c", "x"))
cases <- list(
  list(data = "plain", seqs = plain, pattern = "[[:upper:]]B"),
  list(data = "plain", seqs = plain, pattern = "[[:alpha:]]C"),
  list(data = "plain", seqs = plain, pattern = "[^[:digit:]]"),
  list(data = "plain", seqs = plain, pattern = "[[:alpha:]0-9]C"),
  list(data = "plain", seqs = plain, pattern = "A[^[:upper:]]"),
  list(data = "accents", seqs = accents, pattern = "A[[:alpha:]]B"),
  list(data = "accents", seqs = accents, pattern = "[[:lower:]]"),
  list(data = "accents", seqs = accents, pattern = "B[[:punct:]]"),
  list(data = "brackets", seqs = brackets, pattern = "[\\]A]"),
  list(data = "brackets", seqs = brackets, pattern = "[[:alpha:]\\]]"),
  list(data = "brackets", seqs = brackets, pattern = "[\\\\]"),
  list(data = "brackets", seqs = brackets, pattern = "[]A]"),
  list(data = "brackets", seqs = brackets, pattern = "[^]A]1")
)
rows <- list()
for (case in cases) {
  result <- suppressMessages(TraMineR::seqpm(case$seqs, case$pattern))
  rows[[length(rows) + 1L]] <- data.frame(data = case$data, pattern = case$pattern,
                                          nbocc = result$MTab$nbocc,
                                          index = paste(result$MIndex, collapse = "|"),
                                          stringsAsFactors = FALSE)
}
write.csv(do.call(rbind, rows), file.path(out_dir, "seqpm_posix.csv"), row.names = FALSE, quote = TRUE)
