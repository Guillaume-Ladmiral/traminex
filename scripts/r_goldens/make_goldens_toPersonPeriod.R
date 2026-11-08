# Golden generation for TraMineRextras 0.6.9 toPersonPeriod.
#
# Run from the repository root with:
#   /usr/bin/Rscript scripts/r_goldens/make_goldens_toPersonPeriod.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineRextras))
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineRextras")) == "0.6.9")
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

case_rows <- list()

record_case <- function(name, seqdata, time_col = "timestamp") {
  value <- toPersonPeriod(seqdata)
  stopifnot(identical(names(value), c("id", "state", "timestamp")))
  rows <- data.frame(
    case = rep(name, nrow(value)),
    row = seq_len(nrow(value)),
    id = as.character(value$id),
    state = ifelse(is.na(value$state), "<NA>", as.character(value$state)),
    timestamp = as.integer(value$timestamp),
    time_col = rep(time_col, nrow(value)),
    columns = rep(paste(c("id", "state", time_col), collapse = ";"), nrow(value)),
    contract = rep("R_EXACT", nrow(value)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  case_rows[[length(case_rows) + 1L]] <<- rows
}

custom_ids <- c("zeta", "alpha", "mu")

# Trailing missing cells become structural void under seqdef(right="DEL").
# The oracle removes those rows and ignores the sequence object's start value.
ragged <- data.frame(
  t1 = c("B", "C", "A"),
  t2 = c("A", "B", "C"),
  t3 = c(NA, "A", "B"),
  t4 = c(NA, NA, "C"),
  stringsAsFactors = FALSE
)
record_case(
  "custom_ids_ragged_void",
  suppressMessages(seqdef(ragged, id = custom_ids, start = 7, right = "DEL"))
)

# Interior missing cells are non-response rather than void. toPersonPeriod
# translates the nr marker to a genuine missing state after void filtering.
interior_missing <- data.frame(
  t1 = c("A", "B", "C"),
  t2 = c(NA, "C", "A"),
  t3 = c("B", NA, "B"),
  t4 = c("C", "A", "C"),
  stringsAsFactors = FALSE
)
record_case(
  "interior_missing",
  suppressMessages(seqdef(
    interior_missing,
    id = custom_ids,
    start = 4,
    right = "DEL"
  ))
)

# A marker-free fixture pins the ordinary conversion independently of the
# missing/void branches.
no_missing <- data.frame(
  t1 = c("C", "A", "B"),
  t2 = c("B", "C", "A"),
  t3 = c("A", "B", "C"),
  stringsAsFactors = FALSE
)
record_case(
  "no_missing",
  suppressMessages(seqdef(no_missing, id = custom_ids, start = 12))
)

# TraMineRextras names this column timestamp. The Python compatibility option
# only renames that exact R series, so retain the requested name as metadata.
custom_time <- data.frame(
  t1 = c("B", "A"),
  t2 = c("A", "C"),
  stringsAsFactors = FALSE
)
record_case(
  "custom_time_col",
  suppressMessages(seqdef(custom_time, id = c("person-b", "person-a"), start = 20)),
  time_col = "period"
)

cases <- do.call(rbind, case_rows)
write.csv(
  cases,
  file.path(out_dir, "to_person_period_cases.csv"),
  row.names = FALSE,
  quote = TRUE,
  na = ""
)
