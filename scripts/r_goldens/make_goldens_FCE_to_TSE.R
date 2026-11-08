# Golden generation for TraMineRextras 0.6.9 FCE_to_TSE.
#
# Run from the repository root with:
#   /usr/bin/Rscript scripts/r_goldens/make_goldens_FCE_to_TSE.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineRextras))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
stopifnot(as.character(packageVersion("TraMineRextras")) == "0.6.9")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

case_rows <- list()
artifact_rows <- c(
  "# TraMineRextras 0.6.9 FCE_to_TSE errors and selected Python divergences"
)

scalar_text <- function(value) {
  if (length(value) == 0L || is.na(value)) return("<NA>")
  if (is.numeric(value)) {
    return(trimws(formatC(value, format = "g", digits = 15)))
  }
  as.character(value)
}

case_frame <- function(name, value, contract) {
  stopifnot(identical(names(value), c("id", "time", "event")))
  stopifnot(nrow(value) > 0L)
  data.frame(
    case = rep(name, nrow(value)),
    row = seq_len(nrow(value)),
    id = vapply(value$id, scalar_text, character(1)),
    id_type = rep(class(value$id)[1], nrow(value)),
    time = vapply(value$time, scalar_text, character(1)),
    time_type = rep(class(value$time)[1], nrow(value)),
    event = vapply(value$event, scalar_text, character(1)),
    event_type = rep(class(value$event)[1], nrow(value)),
    columns = rep(paste(names(value), collapse = ";"), nrow(value)),
    nrow = rep(nrow(value), nrow(value)),
    contract = rep(contract, nrow(value)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

record_case <- function(name, value) {
  rows <- case_frame(name, value, "R_EXACT")
  case_rows[[length(case_rows) + 1L]] <<- rows
}

record_selected_case <- function(name, value) {
  if (nrow(value) > 0L) {
    rows <- case_frame(name, value, "SELECTED_DIVERGENCE")
  } else {
    stopifnot(identical(names(value), c("id", "time", "event")))
    rows <- data.frame(
      case = name,
      row = 0L,
      id = "<EMPTY>",
      id_type = class(value$id)[1],
      time = "<EMPTY>",
      time_type = class(value$time)[1],
      event = "<EMPTY>",
      event_type = class(value$event)[1],
      columns = paste(names(value), collapse = ";"),
      nrow = 0L,
      contract = "SELECTED_DIVERGENCE",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  case_rows[[length(case_rows) + 1L]] <<- rows
}

escaped_text <- function(message) {
  paste(strsplit(message, "\n", fixed = TRUE)[[1]], collapse = "\\n")
}

record_artifact <- function(kind, name, message) {
  stopifnot(!grepl("|", kind, fixed = TRUE))
  stopifnot(!grepl("|", name, fixed = TRUE))
  artifact_rows <<- c(
    artifact_rows,
    paste(kind, name, escaped_text(message), sep = "|")
  )
}

record_error <- function(name, expression) {
  message <- tryCatch(
    {
      suppressWarnings(suppressMessages(force(expression)))
      "<NO ERROR>"
    },
    error = function(condition) conditionMessage(condition)
  )
  record_artifact("ERROR", name, message)
}

# Row order is deliberately z, a, m. Event times tie within every nonempty
# person, and NA cells are interleaved with observed cells.
fce <- data.frame(
  pid = c("z", "a", "m"),
  eB = c(2.5, NA, 1.25),
  eA = c(2.5, 3.75, NA),
  eC = c(NA, 3.75, 1.25),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# NULL IDs become 1:nrow; character selectors and an implicit event list use
# the selected column names. Ties sort by event inside each person.
record_case(
  "null_id_character_cols_implicit",
  FCE_to_TSE(fce, cols = c("eB", "eA", "eC"))
)

# With one source row, upstream reinterprets its generated scalar ID (1) as a
# column selector. The selected Python contract keeps the documented generated
# sequence index instead of leaking the first source column into the ID.
single_null_r <- FCE_to_TSE(fce[1, , drop = FALSE], cols = c("eB", "eA"))
stopifnot(identical(single_null_r$id, c("z", "z")))
record_artifact(
  "DIVERGENCE",
  "single_row_null_id_uses_index",
  "R reinterprets generated id=1 as column 1; selected Python keeps generated id=1"
)
record_selected_case(
  "single_row_null_id_uses_index",
  data.frame(
    id = c(1, 1),
    time = c(2.5, 2.5),
    event = c("eA", "eB"),
    stringsAsFactors = FALSE
  )
)

# Numeric column selectors are one-based, and implicit event names are looked
# up from the source data frame.
record_case(
  "id_column_numeric_cols_implicit",
  FCE_to_TSE(fce, id = "pid", cols = 2:4)
)

# A length-one ID is a source-column selector, including a numeric selector.
record_case(
  "id_numeric_selector",
  FCE_to_TSE(fce, id = 1, cols = c("eB", "eA", "eC"))
)

# Longer ID vectors are used directly. Pin both character and numeric modes so
# the Python port cannot coerce identifiers through row iteration.
record_case(
  "explicit_character_ids",
  FCE_to_TSE(fce, id = c("q", "r", "s"), cols = c("eB", "eA"))
)
record_case(
  "explicit_integer_ids",
  FCE_to_TSE(fce, id = c(30L, 10L, 20L), cols = c("eB", "eA"))
)

# The birth/origin event is always at time zero. Explicit event labels also
# make the time-then-event tie ordering observable (alpha precedes beta).
record_case(
  "explicit_events_first_event",
  FCE_to_TSE(
    fce,
    id = "pid",
    cols = c("eB", "eA", "eC"),
    eventlist = c("beta", "alpha", "gamma"),
    firstEvent = "birth"
  )
)

empty_observed <- data.frame(
  pid = c("z", "a"),
  e1 = c(NA_real_, NA_real_),
  e2 = c(NA_real_, NA_real_),
  stringsAsFactors = FALSE
)

# TraMineRextras' 1:0 boundary yields one NA event row and assigns it the last
# person's ID. Preserve that raw observation in the audit artifact, but select
# a safe Python contract: no observed events means a true zero-row result with
# the canonical schema.
all_missing_r <- FCE_to_TSE(
  empty_observed,
  id = "pid",
  cols = c("e1", "e2")
)
stopifnot(
  nrow(all_missing_r) == 1L,
  identical(all_missing_r$id, "a"),
  is.na(all_missing_r$time),
  is.na(all_missing_r$event)
)
record_artifact(
  "DIVERGENCE",
  "all_missing_without_first_event",
  paste0(
    "R returns id=", scalar_text(all_missing_r$id),
    ",time=", scalar_text(all_missing_r$time),
    ",event=", scalar_text(all_missing_r$event),
    " via its 1:0 boundary; selected Python returns zero rows"
  )
)
record_selected_case(
  "all_missing_without_first_event",
  data.frame(
    id = character(),
    time = numeric(),
    event = character(),
    stringsAsFactors = FALSE
  )
)

# An eventless person in the middle triggers the same 1:0 indexing bug and
# overwrites a preceding person's identifier. The selected contract simply
# contributes no rows for that person.
eventless_mixed <- data.frame(
  pid = c("z", "a", "m"),
  e1 = c(2, NA, 1),
  stringsAsFactors = FALSE
)
eventless_r <- FCE_to_TSE(eventless_mixed, id = "pid", cols = "e1")
stopifnot(identical(eventless_r$id, c("a", "m")))
record_artifact(
  "DIVERGENCE",
  "eventless_person_does_not_corrupt_ids",
  paste0(
    "R returns ids=", paste(eventless_r$id, collapse = ","),
    " after 1:0 overwrites the preceding id; selected Python returns ids=z,m"
  )
)
record_selected_case(
  "eventless_person_does_not_corrupt_ids",
  data.frame(
    id = c("z", "m"),
    time = c(2, 1),
    event = c("e1", "e1"),
    stringsAsFactors = FALSE
  )
)

# With firstEvent, the same input has one legitimate time-zero row per person.
record_case(
  "all_missing_with_first_event",
  FCE_to_TSE(
    empty_observed,
    id = "pid",
    cols = c("e1", "e2"),
    firstEvent = "origin"
  )
)

write.csv(
  do.call(rbind, case_rows),
  file.path(out_dir, "fce_to_tse_cases.csv"),
  row.names = FALSE,
  quote = TRUE,
  na = "<NA>"
)

record_error(
  "unknown_character_cols",
  FCE_to_TSE(fce, id = "pid", cols = "not_a_column")
)
record_error(
  "numeric_col_out_of_bounds",
  FCE_to_TSE(fce, id = "pid", cols = 99)
)
record_error(
  "unknown_id_selector",
  FCE_to_TSE(fce, id = "not_an_id", cols = "eB")
)
record_error(
  "empty_cols",
  FCE_to_TSE(fce, id = "pid", cols = integer())
)
record_error(
  "zero_row_input",
  FCE_to_TSE(fce[0, , drop = FALSE], id = "pid", cols = c("eB", "eA"))
)

# R silently accepts mismatched event lists. Record the exact observations for
# auditability while selecting a strict, symmetric Python validation rule.
short_eventlist_r <- FCE_to_TSE(
  fce,
  id = "pid",
  cols = c("eB", "eA"),
  eventlist = "only"
)
stopifnot(sum(is.na(short_eventlist_r$event)) == 2L)
record_artifact(
  "DIVERGENCE",
  "eventlist_too_short",
  paste0(
    "R emits ", sum(is.na(short_eventlist_r$event)),
    " <NA> event labels; selected Python raises ValueError"
  )
)
record_artifact(
  "SELECTED_ERROR",
  "eventlist_too_short",
  "eventlist must have same length as cols"
)

long_eventlist_r <- FCE_to_TSE(
  fce,
  id = "pid",
  cols = "eB",
  eventlist = c("kept", "ignored")
)
stopifnot(identical(unique(long_eventlist_r$event), "kept"))
record_artifact(
  "DIVERGENCE",
  "eventlist_too_long",
  paste0(
    "R ignores the extra label and emits events=",
    paste(unique(long_eventlist_r$event), collapse = ","),
    "; selected Python raises ValueError"
  )
)
record_artifact(
  "SELECTED_ERROR",
  "eventlist_too_long",
  "eventlist must have same length as cols"
)

stopifnot(!any(grepl("<NO ERROR>", artifact_rows, fixed = TRUE)))
writeLines(artifact_rows, file.path(out_dir, "fce_to_tse_errors.txt"))

cat("FCE_TO_TSE_GOLDENS_DONE\n")
