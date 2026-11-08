# Golden generation for TraMineRextras 0.6.9 TSE_to_STS.
#
# Run from the repository root with:
#   /usr/bin/Rscript scripts/r_goldens/make_goldens_TSE_to_STS.R

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
alphabet_rows <- list()
error_rows <- c(
  "# TraMineRextras 0.6.9 TSE_to_STS errors and selected Python divergences"
)

scalar_text <- function(value) {
  if (length(value) == 0L || is.na(value)) return("<NA>")
  as.character(value)
}

record_case <- function(name, value, contract) {
  stopifnot(is.data.frame(value), nrow(value) > 0L, ncol(value) > 0L)
  ids <- rownames(value)
  rows <- lapply(seq_len(nrow(value)), function(i) {
    data.frame(
      case = rep(name, ncol(value)),
      row = rep(i, ncol(value)),
      id = rep(ids[i], ncol(value)),
      column_index = seq_len(ncol(value)),
      column = colnames(value),
      state = vapply(value[i, , drop = TRUE], scalar_text, character(1)),
      nrow = rep(nrow(value), ncol(value)),
      ncol = rep(ncol(value), ncol(value)),
      contract = rep(contract, ncol(value)),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  case_rows[[length(case_rows) + 1L]] <<- do.call(rbind, rows)
  # The alphabet R's seqdef gives the result: traminex's default alphabet (owner
  # decision D2 of 2026-10-19).
  states <- as.character(alphabet(suppressMessages(TraMineR::seqdef(value))))
  alphabet_rows[[length(alphabet_rows) + 1L]] <<- data.frame(
    case = name, position = seq_along(states), state = states, stringsAsFactors = FALSE
  )
}

escaped_text <- function(message) {
  paste(strsplit(message, "\n", fixed = TRUE)[[1]], collapse = "\\n")
}

record_artifact <- function(kind, name, message) {
  stopifnot(!grepl("|", kind, fixed = TRUE))
  stopifnot(!grepl("|", name, fixed = TRUE))
  error_rows <<- c(error_rows, paste(kind, name, escaped_text(message), sep = "|"))
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

# The rows are grouped by person, as required by the upstream loop, while the
# events inside each person are deliberately out of time order. B and C tie at
# 2.75 for z; A and B tie at 1 for a. Default seqe2stm event order follows the
# first occurrence order B, A, C, but event processing sorts ties by event name.
canonical_tse <- data.frame(
  id = c("z", "z", "z", "a", "a", "a"),
  time = c(2.75, 0, 2.75, 1, 3, 1),
  event = c("B", "A", "C", "A", "C", "B"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

record_case(
  "default_stm_default_selectors",
  TSE_to_STS(canonical_tse, tmax = 5),
  "R_EXACT"
)

named_tse <- data.frame(
  person = canonical_tse$id,
  moment = canonical_tse$time,
  action = canonical_tse$event,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
record_case(
  "default_stm_named_selectors",
  TSE_to_STS(
    named_tse,
    id = "person",
    timestamp = "moment",
    event = "action",
    tmax = 5
  ),
  "R_EXACT"
)

# B clears A in this explicit state machine. The input order at p's tied time
# is B then A, but R's order(time, event) applies A then B, leaving state B.
# Z is absent from the state-machine columns and is silently ignored. A custom
# firstState must name the state-machine's initial row.
drop_stm <- seqe2stm(
  c("A", "B"),
  dropList = list(B = "A"),
  firstState = "origin"
)
explicit_tse <- data.frame(
  id = c("p", "p", "p", "q", "q"),
  time = c(0, 0, 3, 1, 4),
  event = c("B", "A", "Z", "B", "A"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
record_case(
  "explicit_stm_ties_unknown_first_state",
  TSE_to_STS(
    explicit_tse,
    stm = drop_stm,
    tmin = 1,
    tmax = 5,
    firstState = "origin"
  ),
  "R_EXACT"
)

# Selectors in the R API are one-based. This case also pins final tmin:tmax
# slicing and numeric IDs becoming character row labels.
numeric_selector_tse <- data.frame(
  noise = c("x", "x", "y", "y"),
  action = c("A", "B", "B", "A"),
  pid = c(20L, 20L, 10L, 10L),
  moment = c(0, 2, 0, 3),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
record_case(
  "numeric_selectors_tmin_slice",
  TSE_to_STS(
    numeric_selector_tse,
    id = 3,
    timestamp = 4,
    event = 2,
    tmin = 2,
    tmax = 4
  ),
  "R_EXACT"
)

# Events after the requested horizon are not discarded: R clamps their
# transition into the final output cell.
late_event_tse <- data.frame(
  id = c("p", "p"),
  time = c(0, 99),
  event = c("A", "B"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
record_case(
  "late_event_clamps_to_tmax",
  TSE_to_STS(late_event_tse, tmax = 4),
  "R_EXACT"
)

# Positive infinity is also clamped to the final cell. The explicit matrix
# avoids conflating this timing contract with the one-event seqe2stm defect.
one_event_stm <- matrix(
  "A",
  nrow = 2L,
  ncol = 1L,
  dimnames = list(c("origin", "A"), "A")
)
record_case(
  "infinite_event_clamps_to_tmax",
  TSE_to_STS(
    data.frame(id = "p", time = Inf, event = "A"),
    stm = one_event_stm,
    tmax = 4,
    firstState = "origin"
  ),
  "R_EXACT"
)

# R's inclusive 1:0 boundary means events from raw time -1 through just below
# zero take effect over the complete requested grid.
record_case(
  "negative_boundary_event_applies_from_start",
  TSE_to_STS(
    data.frame(
      id = c("p", "q"),
      time = c(-0.25, -1),
      event = c("A", "A")
    ),
    stm = one_event_stm,
    tmax = 3,
    firstState = "origin"
  ),
  "R_EXACT"
)

# These are ordinary R state strings. The STSList re-expression must not
# reinterpret them as traminex's structural missing/void markers.
record_case(
  "marker_text_states_remain_real",
  TSE_to_STS(
    data.frame(
      id = c("p", "q"),
      time = c(0, 0),
      event = c("*", "%"),
      check.names = FALSE
    ),
    tmax = 2
  ),
  "R_EXACT"
)

# Upstream does not pass firstState into its automatic seqe2stm call. The
# resulting machine has a "None" initial row, so indexing it with "origin"
# fails. Derive the selected output by supplying the exact machine that the
# automatic call should have generated with firstState forwarded.
custom_first_tse <- data.frame(
  id = c("p", "p"),
  time = c(0, 2),
  event = c("A", "B"),
  stringsAsFactors = FALSE
)
custom_first_error <- tryCatch(
  {
    TSE_to_STS(custom_first_tse, tmax = 4, firstState = "origin")
    "<NO ERROR>"
  },
  error = function(condition) conditionMessage(condition)
)
stopifnot(identical(custom_first_error, "subscript out of bounds"))
record_artifact(
  "DIVERGENCE",
  "auto_stm_custom_first_state_safe",
  paste0(
    "R errors: ", custom_first_error,
    " because auto seqe2stm omits firstState; selected Python forwards first_state"
  )
)
record_case(
  "auto_stm_custom_first_state_safe",
  TSE_to_STS(
    custom_first_tse,
    stm = seqe2stm(c("A", "B"), firstState = "origin"),
    tmax = 4,
    firstState = "origin"
  ),
  "SELECTED_DIVERGENCE"
)

# seqe2stm's final one-column subset drops matrix dimensions. TSE_to_STS then
# sees no state-machine column names and silently ignores its only event. Pin
# the raw output, but select the same transition contract using the intended
# two-state, one-column matrix.
one_event_tse <- data.frame(
  id = c("p", "p"),
  time = c(0, 2),
  event = c("A", "A"),
  stringsAsFactors = FALSE
)
one_event_raw <- TSE_to_STS(one_event_tse, tmax = 4)
stopifnot(identical(as.character(one_event_raw[1, ]), rep("None", 4L)))
record_artifact(
  "DIVERGENCE",
  "one_event_auto_stm_transitions",
  paste0(
    "R returns states=", paste(one_event_raw[1, ], collapse = ","),
    " because seqe2stm drops its one-event dimensions; selected Python keeps the matrix"
  )
)
one_event_stm <- matrix(
  c("A", "A"),
  nrow = 2L,
  ncol = 1L,
  dimnames = list(c("None", "A"), "A")
)
record_case(
  "one_event_auto_stm_transitions",
  TSE_to_STS(one_event_tse, stm = one_event_stm, tmax = 4),
  "SELECTED_DIVERGENCE"
)

# The upstream while-loop only skips contiguous occurrences after processing
# all rows matching one ID. Repeated non-contiguous IDs are therefore visited
# again and overflow the preallocated nseq rows. Select stable first-appearance
# grouping, deriving the expected result from the same rows grouped once.
noncontiguous_tse <- data.frame(
  id = c("p", "q", "p", "q"),
  time = c(0, 0, 2, 2),
  event = c("A", "B", "B", "A"),
  stringsAsFactors = FALSE
)
noncontiguous_error <- tryCatch(
  {
    TSE_to_STS(noncontiguous_tse, tmax = 4)
    "<NO ERROR>"
  },
  error = function(condition) conditionMessage(condition)
)
stopifnot(identical(noncontiguous_error, "subscript out of bounds"))
record_artifact(
  "DIVERGENCE",
  "noncontiguous_ids_group_stably",
  paste0(
    "R errors: ", noncontiguous_error,
    " when IDs repeat non-contiguously; selected Python groups each ID once in stable first order"
  )
)
record_case(
  "noncontiguous_ids_group_stably",
  TSE_to_STS(noncontiguous_tse[c(1, 3, 2, 4), ], tmax = 4),
  "SELECTED_DIVERGENCE"
)

# Source errors are retained verbatim. Selector failures come from data-frame
# indexing before the tmin/tmax validation branches.
record_error(
  "tmax_required",
  TSE_to_STS(canonical_tse)
)
record_error(
  "tmin_below_one",
  TSE_to_STS(canonical_tse, tmin = 0, tmax = 5)
)
record_error(
  "unknown_named_id_selector",
  TSE_to_STS(canonical_tse, id = "missing", tmax = 5)
)
record_error(
  "unknown_numeric_event_selector",
  TSE_to_STS(canonical_tse, event = 99, tmax = 5)
)
record_error(
  "empty_auto_stm",
  TSE_to_STS(
    data.frame(id = character(), time = numeric(), event = character()),
    tmax = 3
  )
)

cases <- do.call(rbind, case_rows)
write.csv(
  cases,
  file.path(out_dir, "tse_to_sts_cases.csv"),
  row.names = FALSE,
  quote = TRUE,
  na = "<NA>"
)
write.csv(
  do.call(rbind, alphabet_rows),
  file.path(out_dir, "tse_to_sts_alphabets.csv"),
  row.names = FALSE,
  quote = TRUE
)
writeLines(error_rows, file.path(out_dir, "tse_to_sts_errors.txt"), useBytes = TRUE)

cat("TSE_TO_STS_GOLDENS_DONE\n")
