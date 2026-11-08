# Live-R goldens for how seqecreate renders event values and matches end.event
# (repair of audit findings TX-AUD-20261010-002 and TX-AUD-20261010-016).
#
# seqecreate turns a non-factor event vector into factor(event), so event labels
# are R's as.character of the values (1e+05, TRUE, 0.5), and it finds end.event
# with `dictionnary[i] == end.event`, which coerces a number with as.character
# too (seqecreate.R:63-73). An end.event that matches no level is ignored: R
# starts intEvent at 0, so its is.null() stop is unreachable.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_seqecreate_rendering.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rows <- list()
record <- function(case, e) {
  add <- function(kind, values) {
    rows[[length(rows) + 1L]] <<- data.frame(
      case = case, kind = kind, position = seq_along(values),
      value = as.character(values), stringsAsFactors = FALSE
    )
  }
  add("level", levels(e))
  add("sequence", as.character(e))
  add("length", format(seqelength(e), digits = 15))
}

record("numeric_end", TraMineR::seqecreate(
  id = c(1, 1, 1, 2, 2), timestamp = c(0, 2, 5, 1, 3), event = c(1, 2, 99, 2, 99), end.event = 99
))
record("labels_1e5", TraMineR::seqecreate(
  id = c(1, 1, 1, 2, 2), timestamp = c(0, 2, 5, 1, 3), event = c(1, 100000, 7, 100000, 7), end.event = 7
))
record("logical", TraMineR::seqecreate(
  id = c(1, 1, 2), timestamp = c(0, 2, 1), event = c(TRUE, FALSE, TRUE)
))
record("end_1e5", TraMineR::seqecreate(
  id = c(1, 1, 2, 2), timestamp = c(0, 2, 1, 4), event = c(1, 100000, 2, 100000), end.event = 1e5
))
record("absent_end", TraMineR::seqecreate(
  id = c(1, 1, 2), timestamp = c(0, 2, 1), event = c(1, 2, 2), end.event = 42
))
record("logical_end", TraMineR::seqecreate(
  id = c(1, 1, 2, 2), timestamp = c(0, 3, 1, 2), event = c(TRUE, FALSE, TRUE, FALSE), end.event = FALSE
))
# R integers render without exponent (as.character(100000L) is "100000").
record("numeric_end_integer", TraMineR::seqecreate(
  id = c(1L, 1L, 1L, 2L, 2L), timestamp = c(0, 2, 5, 1, 3), event = c(1L, 2L, 99L, 2L, 99L), end.event = 99L
))
record("labels_1e5_integer", TraMineR::seqecreate(
  id = c(1L, 1L, 1L, 2L, 2L), timestamp = c(0, 2, 5, 1, 3), event = c(1L, 100000L, 7L, 100000L, 7L),
  end.event = 7L
))
record("double_end", TraMineR::seqecreate(
  id = c(1, 1, 2, 2), timestamp = c(0, 3, 1, 2), event = c(0.5, 1, 0.5, 1), end.event = 1
))

write.csv(do.call(rbind, rows), file.path(out_dir, "seqecreate_rendering.csv"),
          row.names = FALSE, quote = TRUE)
