stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
# Live-R goldens for seqdist(refseq = 0), the most frequent sequence, when frequencies tie
# (audit finding TX-AUD-20261019-001).
#
# R takes seqtab(seqnum(seqdata, with.missing), idxs = 1) and the first row seqfind matches
# (seqdist.R:514, 570-574). seqnum codes the alphabet 0..k-1 and, with with.missing = TRUE,
# the nr code as k (seqnum.R:14-24); seqtab groups the compressed SPS strings (state/duration)
# with tapply, whose factor levels follow R's collation, and sorts the weighted frequencies
# decreasingly (seqtab.R:35-45). Each case below is a tie whose candidates order differently
# under one-based codes, so it pins the zero-based rule.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_refseq_most_frequent.R

suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rows <- list()
quiet <- function(expr) suppressWarnings(suppressMessages(expr))
record <- function(case, data, alphabet, method, weights = NULL, nr = "*", ...) {
  seq <- quiet(TraMineR::seqdef(data, alphabet = alphabet, weights = weights, nr = nr))
  d <- as.numeric(quiet(TraMineR::seqdist(seq, method = method, refseq = 0, ...)))
  inputs <- apply(data, 1, function(r) paste(ifelse(is.na(r), "NA", r), collapse = "-"))
  rows[[length(rows) + 1L]] <<- data.frame(
    case = case, method = method, row = seq_along(d), sequence = inputs,
    weight = if (is.null(weights)) 1 else weights, distance = sprintf("%.17g", d),
    stringsAsFactors = FALSE
  )
}

# 11 states: K is code 10 (key 10/3, after 1/3); one-based codes made it 11/3, before 2/3.
record("eleven_states", rbind(c("K", "K", "K"), c("B", "B", "B"), c("C", "K", "K")),
       LETTERS[1:11], "LCP")
# 9 states and a missing cell. Recorded: R takes A-*-A, so its key does not write the
# missing cell as the code 9 (that would sort A-B-A, 0/1-1/1-0/1, first); the case with
# nr = "Z" tells whether the key writes the object's nr text or a fixed "*".
record("missing_tie",
       rbind(c("A", NA, "A"), c("A", "B", "A"), c("A", NA, "A"), c("A", "B", "A"), c("B", "B", "B")),
       LETTERS[1:9], "OM", with.missing = TRUE, sm = "CONSTANT", indel = 1)
record("missing_tie_nr_letter",
       rbind(c("A", NA, "A"), c("A", "B", "A"), c("A", NA, "A"), c("A", "B", "A"), c("B", "B", "B")),
       LETTERS[1:9], "OM", nr = "Z", with.missing = TRUE, sm = "CONSTANT", indel = 1)
# A leading missing cell: R's key */1-0/2 sorts before 0/3 (punctuation before digits);
# one-based codes wrote 10/1-1/2, which sorts after 1/3.
record("missing_first",
       rbind(c(NA, "A", "A"), c("A", "A", "A"), c(NA, "A", "A"), c("A", "A", "A"), c("B", "B", "B")),
       LETTERS[1:9], "OM", with.missing = TRUE, sm = "CONSTANT", indel = 1)
# Weighted frequencies tie (2 and 1 + 1) between K-K-K and B-B-B with 11 states.
record("weighted_tie", rbind(c("K", "K", "K"), c("B", "B", "B"), c("B", "B", "B"), c("C", "C", "C")),
       LETTERS[1:11], "HAM", weights = c(2, 1, 1, 1))

# seqtab drops a row only when every cell is the nr marker (seqtab.R:17: rowSums(seqdata != nr)
# != 0); a void cell counts as non-nr (batch-7 review M-1). The most frequent rows below are: void
# only (NA everywhere, void under right = "DEL"); nr then a state; nr, a state, void; nr and void
# only (kept); nr only (dropped, so A-A-A is chosen). R refuses literal marker cells beside an
# alphabet, so the last two are built with right = NA (missing cells stay nr) and, for nr and
# void, by setting the last cell to the void level. The sequence column holds the object's cells.
record_object <- function(case, seq) {
  d <- as.numeric(quiet(TraMineR::seqdist(seq, method = "OM", refseq = 0, with.missing = TRUE,
                                          sm = "CONSTANT", indel = 1)))
  rows[[length(rows) + 1L]] <<- data.frame(
    case = case, method = "OM", row = seq_along(d),
    sequence = apply(as.matrix(seq), 1, paste, collapse = "-"), weight = 1,
    distance = sprintf("%.17g", d), stringsAsFactors = FALSE
  )
}
tail_rows <- rbind(c("A", "A", "A"), c("A", "A", "A"), c("B", "A", "A"))
three <- function(r) rbind(r, r, r)
record_object("all_void_rows", quiet(TraMineR::seqdef(rbind(three(c(NA, NA, NA)), tail_rows), alphabet = c("A", "B"))))
record_object("nr_then_state", quiet(TraMineR::seqdef(rbind(three(c(NA, NA, "A")), tail_rows), alphabet = c("A", "B"))))
record_object("nr_state_void", quiet(TraMineR::seqdef(
  rbind(three(c(NA, "A", NA)), c("B", "B", "B"), c("B", "B", "B"), c("A", "A", "A")), alphabet = c("A", "B"))))
nr_only <- quiet(TraMineR::seqdef(rbind(three(c(NA, NA, NA)), tail_rows), alphabet = c("A", "B"), right = NA))
record_object("all_nr_row", nr_only)
nr_void <- nr_only
nr_void[1:3, 3] <- attr(nr_void, "void")
record_object("nr_void_row", nr_void)

write.csv(do.call(rbind, rows), file.path(out_dir, "refseq_most_frequent.csv"), row.names = FALSE, quote = TRUE)
