# Golden generation for the complete TraMineR 2.2-14 seqecreate contract.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

case_rows <- list()
metadata_rows <- list()
case_index <- 1L
metadata_index <- 1L

record_case <- function(case_name, result) {
  for (position in seq_along(result)) {
    case_rows[[case_index]] <<- data.frame(
      case = case_name,
      position = position,
      id_code = seqeid(result)[[position]],
      as_character = as.character(result)[[position]],
      length = seqelength(result)[[position]],
      weight = seqeweight(result)[[position]],
      stringsAsFactors = FALSE
    )
    case_index <<- case_index + 1L
  }
  for (position in seq_along(levels(result))) {
    metadata_rows[[metadata_index]] <<- data.frame(
      case = case_name,
      key = "dictionary",
      position = position,
      value = levels(result)[[position]],
      stringsAsFactors = FALSE
    )
    metadata_index <<- metadata_index + 1L
  }
}

ids <- c("p1", "p1", "p1", "p2", "p2")
timestamps <- c(0, 2, 5, 0, 4)
factor_events <- factor(
  c("A", "B", "END", "B", "END"),
  levels = c("B", "A", "END", "UNUSED")
)
record_case(
  "direct_vectors_end",
  seqecreate(id = ids, timestamp = timestamps, event = factor_events, end.event = "END")
)
record_case(
  "frame_timestamp",
  seqecreate(data.frame(id = ids, timestamp = timestamps, event = factor_events))
)
record_case(
  "frame_time",
  seqecreate(data.frame(id = ids, time = timestamps, event = factor_events))
)

simultaneous <- data.frame(
  id = c("p1", "p1", "p1", "p2", "p2"),
  timestamp = c(0, 2, 2, 0, 0),
  event = factor(
    c("Beta", "Gamma", "Alpha", "Gamma", "Alpha"),
    levels = c("Gamma", "Beta", "Alpha", "UNUSED")
  )
)
record_case("factor_simultaneous", seqecreate(simultaneous))

mixed_case <- data.frame(
  id = c("a", "a", "A", "A"),
  timestamp = c(0, 1, 0, 0),
  event = c("a", "A", "A", "a"),
  stringsAsFactors = FALSE
)
record_case("character_mixed_case", seqecreate(mixed_case))

locale_labels <- c("_", "-", "!", "?", ".", "*", "%", "a", "A", "á", "ä", "b")
record_case(
  "character_locale_order",
  seqecreate(
    id = rep(1, length(locale_labels)),
    timestamp = rep(0, length(locale_labels)),
    event = locale_labels
  )
)

state_matrix <- rbind(
  c("A", "A", "B", "C"),
  c("B", "C", "C", "*"),
  c("A", "*", "B", "B")
)
state_sequences <- suppressMessages(
  seqdef(
    state_matrix,
    alphabet = c("A", "B", "C"),
    labels = c("Alpha", "Beta", "Gamma"),
    missing = "*",
    weights = c(2, 1, 3)
  )
)
record_case("sts_transition_labels", seqecreate(state_sequences, tevent = "transition"))
record_case("sts_state_labels", seqecreate(state_sequences, tevent = "state"))
record_case("sts_period_labels", seqecreate(state_sequences, tevent = "period"))
record_case(
  "sts_state_alphabet",
  seqecreate(state_sequences, tevent = "state", use.labels = FALSE)
)
record_case(
  "sts_unweighted",
  seqecreate(state_sequences, tevent = "state", use.labels = FALSE, weighted = FALSE)
)
custom_matrix <- seqetm(state_sequences, method = "state", use.labels = FALSE)
custom_matrix["A", "B"] <- "custom1,custom2"
record_case(
  "sts_explicit_matrix",
  seqecreate(state_sequences, tevent = custom_matrix, use.labels = FALSE)
)

write.csv(
  do.call(rbind, case_rows),
  file.path(out_dir, "seqecreate_r3_cases.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  do.call(rbind, metadata_rows),
  file.path(out_dir, "seqecreate_r3_metadata.csv"),
  row.names = FALSE,
  quote = TRUE
)

condition_lines <- character()
capture_conditions <- function(case_name, expression) {
  tryCatch(
    withCallingHandlers(
      {
        value <- force(expression)
        condition_lines <<- c(
          condition_lines,
          paste(case_name, "success", paste(as.character(value), collapse = ";"), sep = "|")
        )
      },
      warning = function(condition) {
        condition_lines <<- c(
          condition_lines,
          paste(case_name, "warning", conditionMessage(condition), sep = "|")
        )
        invokeRestart("muffleWarning")
      }
    ),
    error = function(condition) {
      condition_lines <<- c(
        condition_lines,
        paste(case_name, class(condition)[[1]], conditionMessage(condition), sep = "|")
      )
    }
  )
}
capture_conditions(
  "ungrouped_ids",
  seqecreate(id = c(1, 2, 1), timestamp = c(0, 0, 1), event = c("A", "B", "C"))
)
capture_conditions(
  "missing_value",
  seqecreate(id = c(1, 1), timestamp = c(0, NA), event = c("A", "B"))
)
capture_conditions(
  "length_mismatch",
  seqecreate(id = c(1, 1), timestamp = 0, event = c("A", "B"))
)
capture_conditions(
  "absent_columns",
  seqecreate(data.frame(person = c(1, 1), when = c(0, 1), what = c("A", "B")))
)
capture_conditions(
  "invalid_event_characters",
  seqecreate(id = c(1, 1), timestamp = c(0, 1), event = c("A(B", "C,D"))
)
capture_conditions(
  "unknown_end_event",
  seqecreate(id = c(1, 1), timestamp = c(0, 1), event = c("A", "B"), end.event = "NOPE")
)
writeLines(condition_lines, file.path(out_dir, "seqecreate_r3_conditions.txt"), useBytes = TRUE)

cat("SEQECREATE_R3_GOLDENS_DONE\n")
