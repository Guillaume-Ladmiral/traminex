# Public TraMineRextras 0.6.9 evidence for rowmode.
#
# Run from the repository root:
#   Rscript scripts/r_goldens/make_goldens_rowmode.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineRextras))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

if (!identical(as.character(packageVersion("TraMineRextras")), "0.6.9")) {
  stop("rowmode goldens require TraMineRextras 0.6.9")
}

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

csv_quote <- function(value) {
  paste0('"', gsub('"', '""', as.character(value), fixed = TRUE), '"')
}

rows <- list()
append_case <- function(case, call, value) {
  rows[[length(rows) + 1L]] <<- paste(
    csv_quote(case),
    csv_quote(call),
    csv_quote(typeof(value)),
    length(value),
    csv_quote(if (length(value) == 0L) "" else value[[1L]]),
    sep = ","
  )
}

append_case(
  "numeric_tie",
  "TraMineRextras::rowmode(c(2,1,2,1))",
  TraMineRextras::rowmode(c(2, 1, 2, 1))
)
append_case(
  "numeric_decimal",
  "TraMineRextras::rowmode(c(2.5,1.2,2.5))",
  TraMineRextras::rowmode(c(2.5, 1.2, 2.5))
)
append_case(
  "logical_tie",
  "TraMineRextras::rowmode(c(TRUE,FALSE,TRUE,FALSE))",
  TraMineRextras::rowmode(c(TRUE, FALSE, TRUE, FALSE))
)
append_case(
  "factor_tie",
  "TraMineRextras::rowmode(factor(c('A','B'),levels=c('B','A')))",
  TraMineRextras::rowmode(factor(c("A", "B"), levels = c("B", "A")))
)
append_case(
  "missing_character",
  "TraMineRextras::rowmode(c(NA,'B','B','A'))",
  TraMineRextras::rowmode(c(NA, "B", "B", "A"))
)
append_case(
  "scalar_except",
  "TraMineRextras::rowmode(c('AA','AA','A'),except='AA')",
  TraMineRextras::rowmode(c("AA", "AA", "A"), except = "AA")
)
append_case(
  "multiple_except",
  "TraMineRextras::rowmode(c('A','B','C','C'),except=c('C','B'))",
  TraMineRextras::rowmode(c("A", "B", "C", "C"), except = c("C", "B"))
)
append_case(
  "all_excluded",
  "TraMineRextras::rowmode(c('A','A'),except='A')",
  TraMineRextras::rowmode(c("A", "A"), except = "A")
)
append_case(
  "empty_character",
  "TraMineRextras::rowmode(character())",
  TraMineRextras::rowmode(character())
)
append_case(
  "all_na_character",
  "TraMineRextras::rowmode(c(NA_character_,NA_character_))",
  TraMineRextras::rowmode(c(NA_character_, NA_character_))
)
append_case(
  "matrix_flattened",
  "TraMineRextras::rowmode(matrix(c('B','A','A','B'),nrow=2,byrow=TRUE))",
  TraMineRextras::rowmode(matrix(c("B", "A", "A", "B"), nrow = 2, byrow = TRUE))
)
append_case(
  "character_locale_tie",
  "TraMineRextras::rowmode(c('a','B','a','B'))",
  TraMineRextras::rowmode(c("a", "B", "a", "B"))
)
append_case(
  "character_accent_tie",
  "TraMineRextras::rowmode(c('á','ä','á','ä'))",
  TraMineRextras::rowmode(c("á", "ä", "á", "ä"))
)

path <- file.path(out_dir, "rowmode_public_cases.csv")
con <- file(path, open = "wt", encoding = "UTF-8")
writeLines(
  "case,call,r_type,r_length,value",
  con = con,
  sep = "\n",
  useBytes = TRUE
)
writeLines(
  unlist(rows, use.names = FALSE),
  con = con,
  sep = "\n",
  useBytes = TRUE
)
close(con)

cat("ROWMODE_GOLDENS_DONE\n")
