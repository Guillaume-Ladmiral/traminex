# Golden generation for TraMineRextras 0.6.9 plot.SAMM selection mechanics.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineRextras))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
stopifnot(as.character(packageVersion("TraMineRextras")) == "0.6.9")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ragged_raw <- data.frame(
  t1 = c("A", "A", "C", "A"),
  t2 = c("A", "B", "C", "A"),
  t3 = c("B", "B", "C", "A"),
  t4 = c("B", "C", "C", NA),
  t5 = c("C", NA, "C", "B"),
  t6 = c("C", NA, "C", "B"),
  t7 = c("C", NA, "C", NA),
  stringsAsFactors = FALSE
)
ragged <- suppressMessages(
  seqdef(
    ragged_raw,
    alphabet = c("A", "B", "C"),
    states = c("Alpha", "Beta", "Gamma"),
    right = "DEL",
    missing = NA,
    start = 5,
    xtstep = 2,
    tick.last = TRUE,
    cpal = c("#111111", "#222222", "#333333")
  )
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

constant <- suppressMessages(
  seqdef(
    data.frame(
      t1 = c("A", "B"), t2 = c("A", "B"), t3 = c("A", "B"),
      t4 = c("A", "B"), t5 = c("A", "B"), stringsAsFactors = FALSE
    ),
    alphabet = c("A", "B"),
    start = 3
  )
)

cases <- list(
  ragged = seqsamm(ragged, 3),
  repeated = seqsamm(repeated, 2),
  no_transitions = seqsamm(constant, 2)
)

selection_records <- list()
selection_index <- 1L
metadata_records <- list()
metadata_index <- 1L
for (case_name in names(cases)) {
  result <- cases[[case_name]]
  selected_positions <- which(result$transition)
  selected_sequences <- attr(result, "stslist")[selected_positions, ]
  origins <- result[selected_positions, attr(result, "sname")[1]]
  used_origins <- factor(origins)

  metadata_records[[metadata_index]] <- data.frame(
    case = case_name,
    key = "transition_count",
    value = length(selected_positions),
    stringsAsFactors = FALSE
  )
  metadata_index <- metadata_index + 1L
  metadata_records[[metadata_index]] <- data.frame(
    case = case_name,
    key = "group_order",
    value = paste(levels(used_origins), collapse = "\036"),
    stringsAsFactors = FALSE
  )
  metadata_index <- metadata_index + 1L
  metadata_records[[metadata_index]] <- data.frame(
    case = case_name,
    key = "group_labels",
    value = paste(paste("Transition out of", levels(used_origins)), collapse = "\036"),
    stringsAsFactors = FALSE
  )
  metadata_index <- metadata_index + 1L

  if (length(selected_positions)) {
    for (selected_index in seq_along(selected_positions)) {
      for (column in seq_len(ncol(selected_sequences))) {
        selection_records[[selection_index]] <- data.frame(
          case = case_name,
          selected_index_zero = selected_index - 1L,
          source_frame_index_zero = selected_positions[selected_index] - 1L,
          origin = as.character(origins[selected_index]),
          group_label = paste("Transition out of", as.character(origins[selected_index])),
          column_index_zero = column - 1L,
          state_code_zero = as.integer(selected_sequences[selected_index, column]) - 1L,
          stringsAsFactors = FALSE
        )
        selection_index <- selection_index + 1L
      }
    }
  }
}

write.csv(
  if (length(selection_records)) do.call(rbind, selection_records) else data.frame(),
  file.path(out_dir, "plot_samm_r4_selection.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  do.call(rbind, metadata_records),
  file.path(out_dir, "plot_samm_r4_metadata.csv"),
  row.names = FALSE,
  quote = TRUE
)

conditions <- data.frame(
  case = c("supported_dispatch_types", "rf_requires_diss", "r_requires_diss"),
  value = c(
    paste(c("i", "I", "d", "Ht", "dH", "f", "mt", "ms", "r", "rf"), collapse = "\036"),
    "'diss' required for rf plots",
    "For type = 'r', you must provide a distance matrix or a method to compute it"
  ),
  stringsAsFactors = FALSE
)
write.csv(
  conditions,
  file.path(out_dir, "plot_samm_r4_conditions.csv"),
  row.names = FALSE,
  quote = TRUE
)
