# Golden generation for TraMineRextras 0.6.9 seqsamm.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineRextras))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
stopifnot(as.character(packageVersion("TraMineRextras")) == "0.6.9")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

raw <- data.frame(
  t1 = c("A", "A", "C", "A"),
  t2 = c("A", "B", "C", "A"),
  t3 = c("B", "B", "C", "A"),
  t4 = c("B", "C", "C", NA),
  t5 = c("C", NA, "C", "B"),
  t6 = c("C", NA, "C", "B"),
  t7 = c("C", NA, "C", NA),
  stringsAsFactors = FALSE
)
seqdata <- suppressMessages(
  seqdef(
    raw,
    alphabet = c("A", "B", "C"),
    states = c("Alpha", "Beta", "Gamma"),
    right = "DEL",
    missing = NA,
    start = 5,
    xtstep = 2,
    tick.last = TRUE,
    cpal = c("#111111", "#222222", "#333333"),
    missing.color = "#999999"
  )
)
covar <- data.frame(
  age = c(20L, 30L, 40L, 50L),
  sector = factor(c("x", "y", "x", "z"), levels = c("z", "x", "y"))
)

constant <- suppressMessages(
  seqdef(
    data.frame(
      t1 = c("A", "B"), t2 = c("A", "B"), t3 = c("A", "B"),
      t4 = c("A", "B"), t5 = c("A", "B"),
      stringsAsFactors = FALSE
    ),
    alphabet = c("A", "B"),
    start = 3
  )
)

cases <- list(
  ragged_covariates = seqsamm(seqdata, 3, covar = covar),
  no_transitions = seqsamm(constant, 2)
)

repeated <- suppressMessages(
  seqdef(
    data.frame(
      t1 = c("A", "A"), t2 = c("B", "A"), t3 = c("A", "B"),
      t4 = c("A", "B"), t5 = c("C", "A"), t6 = c("C", "C"),
      t7 = c("A", "C"), stringsAsFactors = FALSE
    ),
    alphabet = c("A", "B", "C"),
    states = c("State A", "State B", "State C"),
    start = 4,
    cpal = c("#AA0000", "#00AA00", "#0000AA")
  )
)
repeat_covar <- data.frame(
  score = c(10L, 20L),
  cohort = factor(c("old", "new"), levels = c("new", "old"))
)
repeat_samm <- seqsamm(repeated, 2, covar = repeat_covar)

frame_records <- list()
frame_index <- 1L
for (case_name in names(cases)) {
  result <- cases[[case_name]]
  for (column in names(result)) {
    frame_records[[frame_index]] <- data.frame(
      case = case_name,
      record_type = "column",
      row = "",
      column = column,
      value = as.character(match(column, names(result))),
      detail = typeof(result[[column]]),
      stringsAsFactors = FALSE
    )
    frame_index <- frame_index + 1L
    if (is.factor(result[[column]])) {
      frame_records[[frame_index]] <- data.frame(
        case = case_name,
        record_type = "levels",
        row = "",
        column = column,
        value = paste(levels(result[[column]]), collapse = "\036"),
        detail = if (is.ordered(result[[column]])) "ordered" else "factor",
        stringsAsFactors = FALSE
      )
      frame_index <- frame_index + 1L
    }
  }
  for (row in seq_len(nrow(result))) {
    for (column in names(result)) {
      frame_records[[frame_index]] <- data.frame(
        case = case_name,
        record_type = "cell",
        row = as.character(row),
        column = column,
        value = as.character(result[row, column]),
        detail = "",
        stringsAsFactors = FALSE
      )
      frame_index <- frame_index + 1L
    }
  }
}
write.csv(
  do.call(rbind, frame_records),
  file.path(out_dir, "samm_r3_frame.csv"),
  row.names = FALSE,
  quote = TRUE
)

metadata_records <- list()
metadata_index <- 1L
record_metadata <- function(case_name, section, key, value) {
  metadata_records[[metadata_index]] <<- data.frame(
    case = case_name,
    section = section,
    key = as.character(key),
    value = as.character(value),
    stringsAsFactors = FALSE
  )
  metadata_index <<- metadata_index + 1L
}

for (case_name in names(cases)) {
  result <- cases[[case_name]]
  sequences <- attr(result, "stslist")
  record_metadata(case_name, "result", "class", paste(class(result), collapse = ","))
  record_metadata(case_name, "result", "spell", paste(attr(result, "spell"), collapse = "\036"))
  record_metadata(case_name, "result", "typology", paste(attr(result, "typology"), collapse = "\036"))
  record_metadata(case_name, "result", "sname", paste(attr(result, "sname"), collapse = "\036"))
  record_metadata(case_name, "result", "sublength", attr(result, "sublength"))
  record_metadata(case_name, "sequences", "nrow", nrow(sequences))
  record_metadata(case_name, "sequences", "ncol", ncol(sequences))
  for (attribute in c(
    "start", "void", "nr", "alphabet", "labels", "cpal", "missing.color",
    "xtstep", "tick.last", "Version", "class"
  )) {
    value <- attr(sequences, attribute)
    record_metadata(
      case_name,
      "sequences_attribute",
      attribute,
      paste(as.character(value), collapse = "\036")
    )
  }
  record_metadata(
    case_name,
    "sequences_attribute",
    "names",
    paste(names(sequences), collapse = "\036")
  )
  for (row in seq_len(nrow(sequences))) {
    for (column in seq_len(ncol(sequences))) {
      record_metadata(
        case_name,
        "sequence_text",
        paste(row, column, sep = ","),
        as.character(sequences[row, column])
      )
      record_metadata(
        case_name,
        "sequence_code",
        paste(row, column, sep = ","),
        as.integer(sequences[row, column]) - 1L
      )
    }
  }
}
write.csv(
  do.call(rbind, metadata_records),
  file.path(out_dir, "samm_r3_metadata.csv"),
  row.names = FALSE,
  quote = TRUE
)

# seqsammseq oracle cases, including transition-only, repeated, and empty
# ending-spell selections.
sequence_cases <- list(
  alpha = seqsammseq(cases$ragged_covariates, "Alpha"),
  gamma_empty = seqsammseq(cases$ragged_covariates, "Gamma"),
  unknown_empty = seqsammseq(cases$ragged_covariates, "Unknown"),
  repeated_a = seqsammseq(repeat_samm, "State A")
)
sequence_records <- list()
sequence_index <- 1L
record_sequence <- function(case_name, record_type, row, column, value) {
  sequence_records[[sequence_index]] <<- data.frame(
    case = case_name,
    record_type = record_type,
    row = as.character(row),
    column = column,
    value = as.character(value),
    stringsAsFactors = FALSE
  )
  sequence_index <<- sequence_index + 1L
}
for (case_name in names(sequence_cases)) {
  result <- sequence_cases[[case_name]]
  record_sequence(case_name, "metadata", "", "nrow", nrow(result))
  record_sequence(case_name, "metadata", "", "ncol", ncol(result))
  for (attribute in c(
    "start", "void", "nr", "alphabet", "labels", "cpal", "missing.color",
    "xtstep", "tick.last", "class"
  )) {
    record_sequence(
      case_name,
      "attribute",
      "",
      attribute,
      paste(as.character(attr(result, attribute)), collapse = "\036")
    )
  }
  for (row in seq_len(nrow(result))) {
    for (column in seq_len(ncol(result))) {
      record_sequence(
        case_name, "text", row, as.character(column), as.character(result[row, column])
      )
      record_sequence(
        case_name, "code", row, as.character(column), as.integer(result[row, column]) - 1L
      )
    }
  }
}
write.csv(
  do.call(rbind, sequence_records),
  file.path(out_dir, "samm_consumers_r3_sequences.csv"),
  row.names = FALSE,
  quote = TRUE
)

# seqsammeha cases freeze character-order and factor-level indicator routes,
# person-period/spell-level frames, repeated spells, and no-transition spells.
alpha_character <- c("late", "early", "late")
alpha_factor <- factor(
  alpha_character,
  levels = c("late", "early", "unused"),
  ordered = TRUE
)
repeat_types <- c("route2", "route1", "route2", "route3")
alpha_mixed_case <- c("a", "A", "a")
alpha_locale <- c("a", "á", "b")
eha_cases <- list(
  alpha_character = seqsammeha(
    cases$ragged_covariates, "Alpha", alpha_character, persper = TRUE
  ),
  alpha_factor = seqsammeha(
    cases$ragged_covariates, "Alpha", alpha_factor, persper = TRUE
  ),
  alpha_spell = seqsammeha(
    cases$ragged_covariates, "Alpha", alpha_character, persper = FALSE
  ),
  alpha_mixed_case = seqsammeha(
    cases$ragged_covariates, "Alpha", alpha_mixed_case, persper = TRUE
  ),
  alpha_locale = seqsammeha(
    cases$ragged_covariates, "Alpha", alpha_locale, persper = TRUE
  ),
  gamma_empty = seqsammeha(
    cases$ragged_covariates,
    "Gamma",
    factor(character(), levels = c("unused", "future"), ordered = TRUE),
    persper = TRUE
  ),
  repeated_person = seqsammeha(
    repeat_samm, "State A", repeat_types, persper = TRUE
  ),
  repeated_spell = seqsammeha(
    repeat_samm, "State A", repeat_types, persper = FALSE
  )
)
eha_records <- list()
eha_index <- 1L
for (case_name in names(eha_cases)) {
  result <- eha_cases[[case_name]]
  for (column in names(result)) {
    eha_records[[eha_index]] <- data.frame(
      case = case_name,
      record_type = "column",
      row = "",
      column = column,
      value = as.character(match(column, names(result))),
      detail = typeof(result[[column]]),
      stringsAsFactors = FALSE
    )
    eha_index <- eha_index + 1L
    if (is.factor(result[[column]])) {
      eha_records[[eha_index]] <- data.frame(
        case = case_name,
        record_type = "levels",
        row = "",
        column = column,
        value = paste(levels(result[[column]]), collapse = "\036"),
        detail = if (is.ordered(result[[column]])) "ordered" else "factor",
        stringsAsFactors = FALSE
      )
      eha_index <- eha_index + 1L
    }
  }
  for (row in seq_len(nrow(result))) {
    for (column in names(result)) {
      eha_records[[eha_index]] <- data.frame(
        case = case_name,
        record_type = "cell",
        row = as.character(row),
        column = column,
        value = as.character(result[row, column]),
        detail = "",
        stringsAsFactors = FALSE
      )
      eha_index <- eha_index + 1L
    }
  }
}
write.csv(
  do.call(rbind, eha_records),
  file.path(out_dir, "samm_consumers_r3_eha.csv"),
  row.names = FALSE,
  quote = TRUE
)

condition_lines <- character()
capture_condition <- function(case_name, expression) {
  succeeded <- TRUE
  tryCatch(
    withCallingHandlers(
      force(expression),
      warning = function(condition) {
        condition_lines <<- c(
          condition_lines,
          paste(case_name, "warning", conditionMessage(condition), sep = "|")
        )
        invokeRestart("muffleWarning")
      }
    ),
    error = function(condition) {
      succeeded <<- FALSE
      condition_lines <<- c(
        condition_lines,
        paste(case_name, class(condition)[[1]], conditionMessage(condition), sep = "|")
      )
    }
  )
  if (succeeded) {
    condition_lines <<- c(condition_lines, paste(case_name, "success", sep = "|"))
  }
  invisible(NULL)
}

capture_condition("sublength_one", seqsamm(seqdata, 1))
capture_condition("sublength_width", seqsamm(seqdata, ncol(seqdata)))
capture_condition("sublength_fractional", seqsamm(seqdata, 2.5))
capture_condition("covar_short", seqsamm(seqdata, 3, covar = covar[-4, ]))
capture_condition("not_stslist", seqsamm(raw, 3))
capture_condition(
  "seqsammeha_typology_short",
  seqsammeha(cases$ragged_covariates, "Alpha", c("one", "two"))
)
capture_condition(
  "seqsammeha_typology_long",
  seqsammeha(cases$ragged_covariates, "Alpha", c("one", "two", "three", "four"))
)
capture_condition(
  "seqsammeha_unknown_spell",
  seqsammeha(cases$ragged_covariates, "Unknown", character())
)
capture_condition("seqsammseq_not_samm", seqsammseq(raw, "A"))
capture_condition("seqsammeha_not_samm", seqsammeha(raw, "A", character()))

writeLines(
  condition_lines,
  file.path(out_dir, "samm_r3_errors.txt"),
  useBytes = TRUE
)

cat("SAMM_R3_GOLDENS_DONE\n")
