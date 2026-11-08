# Live-R goldens for seqnum (repair of audit finding TX-AUD-20261010-021).
#
# seqnum renumbers the states 0..n-1 in alphabet order (seqnum.R:9-37). With
# with.missing = TRUE the missing code becomes the number n (the position of nr
# in the factor levels minus one) while the alphabet stays 0..n-1; otherwise the
# nr and void symbols are kept. The rows record, per case, R's attributes and
# the cells of the returned object (as.character of each column).
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_seqnum.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rows_of <- function(strings) as.data.frame(do.call(rbind, strsplit(strings, "")), stringsAsFactors = FALSE)
letters_data <- rows_of(c("ABCA", "CCBA", "BACC"))
missing_data <- rows_of(c("A*BC", "*CBA", "BAC*", "AB**"))
numeric_data <- rows_of(c("1231", "3321", "2133"))
cases <- list(
  list(case = "unsorted_alphabet", seqs = suppressMessages(TraMineR::seqdef(letters_data, alphabet = c("C", "A", "B"),
       labels = c("lab-C", "lab-A", "lab-B"))), with.missing = FALSE),
  list(case = "missing_void", seqs = suppressMessages(TraMineR::seqdef(missing_data, alphabet = c("A", "B", "C"),
       missing = "*")), with.missing = FALSE),
  list(case = "missing_void_with_missing", seqs = suppressMessages(TraMineR::seqdef(missing_data,
       alphabet = c("A", "B", "C"), missing = "*")), with.missing = TRUE),
  list(case = "numeric_states", seqs = suppressMessages(TraMineR::seqdef(numeric_data)), with.missing = FALSE)
)
rows <- list()
add <- function(case, quantity, values) {
  rows[[length(rows) + 1L]] <<- data.frame(case = case, quantity = quantity, position = seq_along(values),
                                           value = as.character(values), stringsAsFactors = FALSE)
}
for (cs in cases) {
  # R's default (with.missing = FALSE) is taken by the cases that use it.
  result <- if (cs$with.missing) TraMineR::seqnum(cs$seqs, with.missing = TRUE) else TraMineR::seqnum(cs$seqs)
  add(cs$case, "input_alphabet", attr(cs$seqs, "alphabet"))
  add(cs$case, "alphabet", attr(result, "alphabet"))
  add(cs$case, "nr", attr(result, "nr"))
  add(cs$case, "void", attr(result, "void"))
  add(cs$case, "labels", attr(result, "labels"))
  add(cs$case, "cells", apply(sapply(result, as.character), 1, paste, collapse = "|"))
}
write.csv(do.call(rbind, rows), file.path(out_dir, "seqnum_cases.csv"), row.names = FALSE, quote = TRUE)
