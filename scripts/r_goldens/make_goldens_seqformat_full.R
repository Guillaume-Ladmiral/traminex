# Golden generation for the complete TraMineR 2.2-14 seqformat contract.
#
# Task 12 writes the input-normalization cases and conditions below. Task 14
# extends this script with the complete output-route matrix.
# Run from the repository root with:
#   /usr/bin/Rscript scripts/r_goldens/make_goldens_seqformat_full.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

case_rows <- list()

adaptation_metadata <- function(value, nr = "*", void = "%", real_alphabet = NULL, start = 1L) {
  matrix_value <- as.matrix(value)
  if (is.null(real_alphabet)) {
    # traminex orders states the input does not declare as R's seqdef would
    # (owner decision D2 of 2026-10-19): seqstatl's sort(unique()), re-sorted by
    # value when every state is an R number (seqstatl.R:21-30), equal numbers in
    # collation order.
    flat <- as.vector(t(matrix_value))
    states <- sort(unique(as.character(flat[!is.na(flat)])))
    numbers <- suppressWarnings(as.numeric(states))
    real_alphabet <- if (length(states) && all(!is.na(numbers))) states[order(numbers)] else states
  }

  nr_used <- FALSE
  void_used <- FALSE
  for (i in seq_len(nrow(matrix_value))) {
    observed <- which(!is.na(matrix_value[i, ]))
    last_observed <- if (length(observed) == 0L) 0L else max(observed)
    missing_positions <- which(is.na(matrix_value[i, ]))
    nr_used <- nr_used || any(missing_positions <= last_observed)
    void_used <- void_used || any(missing_positions > last_observed)
  }

  alphabet <- as.character(real_alphabet)
  collision_free <- function(preferred, occupied, role) {
    if (!preferred %in% occupied) return(preferred)
    candidate <- paste0("<", role, ">")
    suffix <- 1L
    while (candidate %in% occupied) {
      candidate <- paste0("<", role, ":", suffix, ">")
      suffix <- suffix + 1L
    }
    candidate
  }
  if (nr_used) nr <- collision_free(nr, alphabet, "NR")
  if (nr_used && !nr %in% alphabet) alphabet <- c(alphabet, nr)
  if (void_used) void <- collision_free(void, alphabet, "VOID")
  if (void_used && !void %in% alphabet) alphabet <- c(alphabet, void)
  list(
    alphabet = paste(alphabet, collapse = ";"),
    nr = nr,
    void = void,
    nr_code = if (nr_used && nr %in% alphabet) match(nr, alphabet) - 1L else NA_integer_,
    void_code = if (void_used && void %in% alphabet) match(void, alphabet) - 1L else NA_integer_,
    start = as.integer(start)
  )
}

record_case <- function(
  name,
  value,
  adapt_nr = "*",
  adapt_void = "%",
  real_alphabet = NULL,
  adapt_start = 1L
) {
  issues <- attr(value, "issues")
  issue_text <- function(name) {
    if (is.null(issues)) return("<NONE>")
    paste(as.integer(issues[[name]]), collapse = ";")
  }
  value <- as.data.frame(value, stringsAsFactors = FALSE)
  metadata <- adaptation_metadata(value, adapt_nr, adapt_void, real_alphabet, adapt_start)
  ids <- rownames(value)
  if (is.null(ids)) ids <- as.character(seq_len(nrow(value)))
  rows <- vector("list", nrow(value) * ncol(value))
  k <- 1L
  for (i in seq_len(nrow(value))) {
    for (j in seq_len(ncol(value))) {
      cell <- value[[j]][i]
      rows[[k]] <- data.frame(
        case = name,
        row = i,
        id = ids[i],
        position = j,
        column = names(value)[j],
        value = if (is.na(cell)) "<NA>" else as.character(cell),
        alphabet = metadata$alphabet,
        nr = metadata$nr,
        void = metadata$void,
        nr_code = metadata$nr_code,
        void_code = metadata$void_code,
        start = metadata$start,
        issues_present = !is.null(issues),
        issues_truncated = issue_text("truncated"),
        issues_empty_seq = issue_text("empty.seq"),
        issues_start_before_birth = issue_text("start.before.birth"),
        issues_spell_before_birth = issue_text("spell.before.birth"),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  case_rows[[length(case_rows) + 1L]] <<- do.call(rbind, rows)
}

# Extended STS input plus var filtering.
sts_extended <- data.frame(
  ignored = c("x", "y"),
  t1 = c("A", "B"),
  t2 = c("A", "@"),
  t3 = c("B", "C"),
  t4 = c("@", "C"),
  stringsAsFactors = FALSE,
  row.names = c("case-z", "case-a")
)
record_case(
  "sts_extended",
  suppressMessages(seqformat(
    sts_extended,
    var = c("t1", "t2", "t3", "t4"),
    from = "STS",
    to = "STS",
    missing = "@"
  ))
)

# Compressed STS input with a non-autodetectable separator and ragged rows.
sts_compressed <- matrix(
  c("A;B;@", "C;@"),
  ncol = 1,
  dimnames = list(c("compressed-1", "compressed-2"), "sequence")
)
record_case(
  "sts_compressed_custom_right_na",
  suppressMessages(seqformat(
    sts_compressed,
    from = "STS",
    to = "STS",
    stsep = ";",
    missing = "@",
    right = NA
  ))
)
record_case(
  "sts_compressed_custom_right_del",
  suppressMessages(seqformat(
    sts_compressed,
    from = "STS",
    to = "STS",
    stsep = ";",
    missing = "@",
    right = "DEL"
  ))
)

sts_compressed_auto <- matrix(
  c("A-B-B", "C-@-A"),
  ncol = 1,
  dimnames = list(c("auto-1", "auto-2"), "sequence")
)
record_case(
  "sts_compressed_auto",
  suppressMessages(seqformat(
    sts_compressed_auto,
    from = "STS",
    to = "STS",
    missing = "@"
  ))
)

sts_compressed_empty <- matrix(
  c("A-B", ""),
  ncol = 1,
  dimnames = list(c("nonempty", "empty"), "sequence")
)
record_case(
  "sts_compressed_empty_row",
  suppressMessages(seqformat(
    sts_compressed_empty,
    from = "STS",
    to = "STS",
    missing = "@"
  ))
)

# Base R's strsplit() discards exactly one terminal empty component.  This is
# observable for compressed state strings and differs from Python's re.split().
sts_compressed_trailing <- matrix(
  c("A-", "-", "--", "A--B"),
  ncol = 1,
  dimnames = list(c("a-trailing", "one-separator", "two-separators", "interior"), "sequence")
)
record_case(
  "sts_compressed_trailing_separator",
  suppressMessages(seqformat(
    sts_compressed_trailing,
    from = "STS",
    to = "STS",
    missing = "@"
  ))
)

# Default marker spellings remain ordinary states when no raw cell is
# classified as nonresponse or structural void.
sts_marker_states <- matrix(
  c("A", "*", "%"),
  nrow = 1,
  dimnames = list("marker-states", c("t1", "t2", "t3"))
)
record_case(
  "sts_marker_states",
  suppressMessages(seqformat(
    sts_marker_states,
    from = "STS",
    to = "STS",
    missing = "@"
  ))
)

# R's %in% comparison coerces numeric states against a character missing code.
sts_numeric_missing <- matrix(
  c(1, 2, 2, 1),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(c("numeric-1", "numeric-2"), c("t1", "t2"))
)
record_case(
  "sts_numeric_missing_token",
  suppressMessages(seqformat(
    sts_numeric_missing,
    from = "STS",
    to = "STS",
    missing = "1"
  ))
)
sts_numeric_rendering <- matrix(
  c(1e5, 100001, 123456, Inf, -Inf, 1.234567890123456),
  nrow = 1,
  dimnames = list("numeric-rendering", c("t1", "t2", "t3", "t4", "t5", "t6"))
)
record_case(
  "sts_numeric_rendering",
  suppressMessages(seqformat(
    sts_numeric_rendering,
    from = "STS",
    to = "STS"
  ))
)
sts_numeric_frame_promotion <- data.frame(
  integer_state = c(100000L, 100001L),
  double_state = c(1.5, 2.5),
  row.names = c("promoted-1", "promoted-2")
)
record_case(
  "sts_numeric_frame_promotion",
  suppressMessages(seqformat(
    sts_numeric_frame_promotion,
    from = "STS",
    to = "STS"
  ))
)
sts_logical_numeric_promotion <- data.frame(
  logical_state = c(TRUE, FALSE),
  numeric_state = c(1, 2.5),
  row.names = c("logical-1", "logical-2")
)
record_case(
  "sts_logical_numeric_promotion",
  suppressMessages(seqformat(
    sts_logical_numeric_promotion,
    from = "STS",
    to = "STS"
  ))
)

# Literal default-marker states and actual missing cells must remain distinct
# in the integer-coded Python adaptation.
sts_marker_collision <- rbind(
  c("A", "*", "@", "B"),
  c("%", "A", "@", "@")
)
dimnames(sts_marker_collision) <- list(
  c("collision-1", "collision-2"),
  c("t1", "t2", "t3", "t4")
)
record_case(
  "sts_marker_collision",
  suppressMessages(seqformat(
    sts_marker_collision,
    from = "STS",
    to = "STS",
    missing = "@"
  ))
)

# A real sequence object infers from="STS" and supplies its own nr/void codes.
sts_object <- suppressMessages(seqdef(
  rbind(c("A", "M", "B"), c("B", "B", "M")),
  alphabet = c("A", "B"),
  id = c("object-b", "object-a"),
  missing = "M",
  nr = "M",
  void = "V",
  right = NA
))
record_case(
  "sts_object_inferred",
  suppressWarnings(suppressMessages(seqformat(sts_object, to = "STS", right = NA))),
  adapt_nr = "M",
  adapt_void = "V",
  real_alphabet = c("A", "B")
)

# Extended SPS with default (state,duration) couples.
sps_default <- matrix(
  c("(A,2)", "(B,1)", "(C,1)", "(@,2)"),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(c("sps-1", "sps-2"), c("spell1", "spell2"))
)
record_case(
  "sps_default",
  suppressMessages(seqformat(sps_default, from = "SPS", to = "STS", missing = "@"))
)

# Compressed SPS with custom between-spell, wrapper, and state/duration tokens.
sps_custom <- matrix(
  c("{A;2}~{B;1}", "{C;1}~{@;2}"),
  ncol = 1,
  dimnames = list(c("sps-custom-1", "sps-custom-2"), "sequence")
)
record_case(
  "sps_custom_compressed",
  suppressMessages(seqformat(
    sps_custom,
    from = "SPS",
    to = "STS",
    stsep = "~",
    SPS.in = list(xfix = "{}", sdsep = ";"),
    missing = "@"
  ))
)

# Executable SPS quirks: the final internal '-' decompression splits states,
# literal "NA" becomes missing, fractional duration truncates, and non-positive
# durations still emit once.
sps_edges <- matrix(
  c("(A-B,1)", "(NA,1)", "(F-NA,1)", "(C,0)", "(D,2.9)", "(E,-2)"),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(c("edges-1", "edges-2"), c("spell1", "spell2", "spell3"))
)
record_case(
  "sps_coercion_edges",
  suppressMessages(seqformat(sps_edges, from = "SPS", to = "STS"))
)

# The final SPS-to-STS seqdecomp() also uses base R's terminal-empty rule.
sps_trailing_state <- matrix(
  c("(A-,1)", NA_character_),
  nrow = 1,
  dimnames = list("sps-trailing", c("episode-1", "episode-2"))
)
record_case(
  "sps_trailing_state_separator",
  suppressMessages(seqformat(sps_trailing_state, from = "SPS", to = "STS"))
)

# SPELL process-time inputs exercise named selectors, overlap direction,
# fillblanks, missing statuses, and limit truncation.
spell_process <- data.frame(
  who = c("b", "b", "b", "a", "a"),
  lo = c(1L, 3L, 6L, 2L, 4L),
  hi = c(3L, 4L, 8L, 3L, 7L),
  state = c("A", "B", "C", "C", "@"),
  stringsAsFactors = FALSE
)
record_case(
  "spell_process_overwrite",
  suppressWarnings(suppressMessages(seqformat(
    spell_process,
    from = "SPELL",
    to = "STS",
    id = "who",
    begin = "lo",
    end = "hi",
    status = "state",
    process = TRUE,
    limit = 6,
    overwrite = TRUE,
    fillblanks = "_",
    missing = "@",
    right = NA
  )))
)
record_case(
  "spell_process_keep_earlier",
  suppressWarnings(suppressMessages(seqformat(
    spell_process,
    from = "SPELL",
    to = "STS",
    id = "who",
    begin = "lo",
    end = "hi",
    status = "state",
    process = TRUE,
    limit = 6,
    overwrite = FALSE,
    fillblanks = "_",
    missing = "@"
  )))
)

# Factor statuses retain the complete declared level order. When blanks are
# filled, TraMineR appends the fill state after those levels.
spell_factor <- data.frame(
  who = c("p", "p", "q"),
  begin = c(1L, 3L, 1L),
  end = c(1L, 3L, 1L),
  state = factor(c("B", "A", "A"), levels = c("A", "B", "UNUSED"))
)
record_case(
  "spell_factor_status",
  suppressMessages(seqformat(
    spell_factor,
    from = "SPELL",
    to = "STS",
    id = "who",
    begin = "begin",
    end = "end",
    status = "state",
    process = TRUE,
    limit = 3,
    fillblanks = "G"
  )),
  real_alphabet = c(levels(spell_factor$state), "G")
)

spell_numeric_status <- data.frame(
  who = c("p", "q"),
  begin = c(1L, 1L),
  end = c(1L, 1L),
  state = c(1, 2)
)
record_case(
  "spell_numeric_missing_token",
  suppressMessages(seqformat(
    spell_numeric_status,
    from = "SPELL",
    to = "STS",
    id = "who",
    begin = "begin",
    end = "end",
    status = "state",
    process = TRUE,
    limit = 2,
    missing = "1"
  ))
)

# Calendar-time SPELL input with positional selectors and explicit axis bounds.
spell_calendar <- data.frame(
  state = c("A", "B", "C"),
  who = c(10L, 10L, 20L),
  hi = c(2002L, 2004L, 2003L),
  lo = c(2001L, 2004L, 2002L),
  stringsAsFactors = FALSE
)
record_case(
  "spell_calendar",
  suppressMessages(seqformat(
    spell_calendar,
    from = "SPELL",
    to = "STS",
    id = 2,
    begin = 4,
    end = 3,
    status = 1,
    process = FALSE,
    tmin = 2000,
    tmax = 2004,
    fillblanks = "_"
  )),
  adapt_start = 2000L
)
record_case(
  "spell_calendar_var_relative",
  suppressMessages(seqformat(
    spell_calendar,
    var = c("who", "lo", "hi", "state"),
    from = "SPELL",
    to = "STS",
    id = 1,
    begin = 2,
    end = 3,
    status = 4,
    process = FALSE,
    tmin = 2000,
    tmax = 2004,
    fillblanks = "_"
  )),
  adapt_start = 2000L
)
record_case(
  "spell_calendar_inferred_one_bound",
  suppressMessages(seqformat(
    spell_calendar,
    from = "SPELL",
    to = "STS",
    id = 2,
    begin = 4,
    end = 3,
    status = 1,
    process = FALSE,
    tmin = 1900,
    tmax = NULL,
    fillblanks = "_"
  )),
  adapt_start = 2001L
)

# Calendar spells aligned to process time with external birth dates.
spell_years <- data.frame(
  who = c("p2", "p2", "p1"),
  begin = c(2002L, 2004L, 1999L),
  end = c(2003L, 2005L, 2001L),
  state = c("A", "B", "C"),
  stringsAsFactors = FALSE
)
births <- data.frame(
  person = c("p1", "p2"),
  birth = c(1998L, 2001L),
  stringsAsFactors = FALSE
)
record_case(
  "spell_year_to_age",
  suppressMessages(seqformat(
    spell_years,
    from = "SPELL",
    to = "STS",
    id = "who",
    begin = "begin",
    end = "end",
    status = "state",
    process = TRUE,
    pdata = births,
    pvar = c("person", "birth"),
    limit = 5
  ))
)
spell_fractional_birth <- data.frame(
  who = "fractional",
  begin = 2001L,
  end = 2002L,
  state = "A",
  stringsAsFactors = FALSE
)
fractional_birth <- data.frame(
  person = "fractional",
  birth = 2000.5,
  stringsAsFactors = FALSE
)
record_case(
  "spell_fractional_birth",
  suppressMessages(seqformat(
    spell_fractional_birth,
    from = "SPELL",
    to = "STS",
    id = "who",
    begin = "begin",
    end = "end",
    status = "state",
    process = TRUE,
    pdata = fractional_birth,
    pvar = c("person", "birth"),
    limit = 4
  ))
)
births_missing <- births
births_missing$birth[births_missing$person == "p2"] <- NA_integer_
record_case(
  "spell_year_to_age_missing_birth",
  suppressWarnings(suppressMessages(seqformat(
    spell_years,
    from = "SPELL",
    to = "STS",
    id = "who",
    begin = "begin",
    end = "end",
    status = "state",
    process = TRUE,
    pdata = births_missing,
    pvar = c("person", "birth"),
    limit = 5
  )))
)

spell_birth_edges <- data.frame(
  who = c("crossing", "before", "normal"),
  begin = c(1999L, 1998L, 2001L),
  end = c(2001L, 1999L, 2002L),
  state = c("X", "Y", "Z"),
  stringsAsFactors = FALSE
)
birth_edges <- data.frame(
  person = c("crossing", "before", "normal"),
  birth = c(2000L, 2000L, 2000L),
  stringsAsFactors = FALSE
)
record_case(
  "spell_birth_edges",
  suppressWarnings(suppressMessages(seqformat(
    spell_birth_edges,
    from = "SPELL",
    to = "STS",
    id = "who",
    begin = "begin",
    end = "end",
    status = "state",
    process = TRUE,
    pdata = birth_edges,
    pvar = c("person", "birth"),
    limit = 3
  )))
)

spell_no_positive_end <- data.frame(
  who = "empty",
  begin = NA_integer_,
  end = NA_integer_,
  state = "A",
  stringsAsFactors = FALSE
)
record_case(
  "spell_no_positive_end",
  suppressWarnings(suppressMessages(seqformat(
    spell_no_positive_end,
    from = "SPELL",
    to = "STS",
    id = "who",
    begin = "begin",
    end = "end",
    status = "state",
    process = TRUE,
    limit = 2
  )))
)
record_case(
  "spell_year_to_age_auto",
  suppressMessages(seqformat(
    spell_years,
    from = "SPELL",
    to = "STS",
    id = "who",
    begin = "begin",
    end = "end",
    status = "state",
    process = TRUE,
    pdata = "auto",
    limit = 5
  ))
)

spell_year_gap <- data.frame(
  who = c("gap", "gap"),
  begin = c(2000L, 2002L),
  end = c(2000L, 2002L),
  state = c("A", "B"),
  stringsAsFactors = FALSE
)
gap_birth <- data.frame(person = "gap", birth = 2000L, stringsAsFactors = FALSE)
record_case(
  "spell_year_to_age_fill_quirk",
  suppressMessages(seqformat(
    spell_year_gap,
    from = "SPELL",
    to = "STS",
    id = "who",
    begin = "begin",
    end = "end",
    status = "state",
    process = TRUE,
    pdata = gap_birth,
    pvar = c("person", "birth"),
    fillblanks = "G",
    limit = 3
  ))
)

# Age spells aligned to a calendar axis with external birth dates.
spell_ages <- data.frame(
  who = c("p2", "p2", "p1"),
  begin = c(1L, 3L, 1L),
  end = c(2L, 3L, 2L),
  state = c("A", "B", "C"),
  stringsAsFactors = FALSE
)
record_case(
  "spell_age_to_year",
  suppressMessages(seqformat(
    spell_ages,
    from = "SPELL",
    to = "STS",
    id = "who",
    begin = "begin",
    end = "end",
    status = "state",
    process = FALSE,
    pdata = births,
    pvar = c("person", "birth"),
    tmin = 1998,
    tmax = 2005
  )),
  adapt_start = 1998L
)

cases <- do.call(rbind, case_rows)
write.csv(
  cases,
  file.path(out_dir, "seqformat_input_cases.csv"),
  row.names = FALSE,
  quote = TRUE,
  na = "<NA>"
)

conditions <- character()
clean_condition <- function(x) {
  text <- gsub("[\r\n]+", "\\\\n", conditionMessage(x))
  text <- sub("\\\\n$", "", text)
  gsub("In record_(error|warnings)\\(\\)", "In seqformat()", text)
}

record_error <- function(name, expr) {
  value <- tryCatch(
    {
      suppressWarnings(suppressMessages(force(expr)))
      "<NO ERROR>"
    },
    error = function(e) clean_condition(e)
  )
  conditions <<- c(conditions, paste("ERROR", name, value, sep = "|"))
}

record_warnings <- function(name, expr) {
  found <- character()
  value <- tryCatch(
    withCallingHandlers(
      {
        force(expr)
        NULL
      },
      message = function(m) {
        text <- clean_condition(m)
        if (grepl("[!!]", text, fixed = TRUE)) found <<- c(found, text)
        invokeRestart("muffleMessage")
      }
    ),
    error = function(e) {
      found <<- c(found, paste0("ERROR: ", clean_condition(e)))
      NULL
    }
  )
  invisible(value)
  if (length(found) == 0L) found <- "<NO WARNING>"
  for (i in seq_along(found)) {
    warning_name <- if (length(found) == 1L) name else paste0(name, "_", i)
    conditions <<- c(conditions, paste("WARN", warning_name, found[[i]], sep = "|"))
  }
}

record_error(
  "missing_from_non_stslist",
  seqformat(matrix(c("A", "B"), nrow = 1), to = "STS")
)
record_error(
  "missing_to",
  seqformat(matrix(c("A", "B"), nrow = 1), from = "STS")
)
record_error(
  "unknown_from",
  seqformat(matrix(c("A", "B"), nrow = 1), from = "BAD", to = "STS")
)
record_error(
  "unknown_to",
  seqformat(matrix(c("A", "B"), nrow = 1), from = "STS", to = "BAD")
)
record_error(
  "compressed_separator_required",
  seqformat(matrix("ABC", ncol = 1), from = "STS", to = "STS")
)
record_error(
  "compressed_numeric_noncharacter",
  seqformat(matrix(c(1, 2), ncol = 1), from = "STS", to = "STS")
)
record_error(
  "compressed_numeric_negative_separator",
  seqformat(matrix(-1, ncol = 1), from = "STS", to = "STS")
)
record_error(
  "compressed_all_missing",
  seqformat(matrix(NA_character_, ncol = 1), from = "STS", to = "STS")
)
record_error(
  "compressed_partial_missing",
  seqformat(matrix(c("A", NA_character_), ncol = 1), from = "STS", to = "STS")
)
record_error(
  "compressed_partial_missing_numeric",
  seqformat(matrix(c(1, NA_real_), ncol = 1), from = "STS", to = "STS")
)
record_error(
  "compressed_empty_only",
  seqformat(matrix("", ncol = 1), from = "STS", to = "STS")
)
record_error(
  "sts_bad_data",
  seqformat(1, from = "STS", to = "STS")
)
record_error(
  "sps_rejects_stslist",
  seqformat(sts_object, from = "SPS", to = "STS")
)
record_error(
  "spell_requires_dataframe",
  seqformat(matrix(c(1, 1, 2, "A"), nrow = 1), from = "SPELL", to = "STS")
)
record_error(
  "bad_missing",
  seqformat(matrix(c("A", "B"), nrow = 1), from = "STS", to = "STS", missing = 1)
)
record_error(
  "bad_stsep",
  seqformat(matrix("A-B", ncol = 1), from = "STS", to = "STS", stsep = 1)
)
record_error(
  "multichar_stsep",
  seqformat(matrix("A--B", ncol = 1), from = "STS", to = "STS", stsep = "--")
)
record_error(
  "spell_bad_selector",
  seqformat(spell_process, from = "SPELL", to = "STS", id = c("who", "lo"))
)
spell_missing_id <- spell_process
spell_missing_id$who[2] <- NA_character_
record_error(
  "spell_missing_id",
  seqformat(
    spell_missing_id,
    from = "SPELL",
    to = "STS",
    id = "who",
    begin = "lo",
    end = "hi",
    status = "state"
  )
)
record_error(
  "spell_bad_limit",
  seqformat(spell_process, from = "SPELL", to = "STS", limit = "six")
)
record_error(
  "spell_bad_fillblanks",
  seqformat(spell_process, from = "SPELL", to = "STS", fillblanks = 0)
)
record_error(
  "spell_empty_fillblanks",
  seqformat(spell_process, from = "SPELL", to = "STS", fillblanks = "")
)
record_error(
  "spell_multichar_fillblanks",
  seqformat(spell_process, from = "SPELL", to = "STS", fillblanks = "GG")
)

spell_character_time <- transform(spell_process, lo = as.character(lo))
record_error(
  "spell_character_time",
  seqformat(
    spell_character_time,
    from = "SPELL",
    to = "STS",
    id = "who",
    begin = "lo",
    end = "hi",
    status = "state"
  )
)

spell_noninteger <- transform(spell_process, lo = as.numeric(lo))
spell_noninteger$lo[1] <- 1.5
record_error(
  "spell_noninteger_time",
  seqformat(
    spell_noninteger,
    from = "SPELL",
    to = "STS",
    id = "who",
    begin = "lo",
    end = "hi",
    status = "state"
  )
)

spell_negative <- spell_process
spell_negative$lo[1] <- 0L
record_error(
  "spell_start_before_one",
  seqformat(
    spell_negative,
    from = "SPELL",
    to = "STS",
    id = "who",
    begin = "lo",
    end = "hi",
    status = "state"
  )
)

spell_reversed <- spell_process
spell_reversed$hi[2] <- 2L
record_error(
  "spell_end_before_start",
  seqformat(
    spell_reversed,
    from = "SPELL",
    to = "STS",
    id = "who",
    begin = "lo",
    end = "hi",
    status = "state"
  )
)

record_error(
  "pdata_requires_pvar",
  seqformat(
    spell_years,
    from = "SPELL",
    to = "STS",
    id = "who",
    begin = "begin",
    end = "end",
    status = "state",
    pdata = births
  )
)

bad_births <- transform(births, birth = as.character(birth))
record_error(
  "pdata_birth_must_be_integer",
  seqformat(
    spell_years,
    from = "SPELL",
    to = "STS",
    id = "who",
    begin = "begin",
    end = "end",
    status = "state",
    pdata = bad_births,
    pvar = c("person", "birth")
  )
)

record_error(
  "pdata_missing_id",
  seqformat(
    spell_years,
    from = "SPELL",
    to = "STS",
    id = "who",
    begin = "begin",
    end = "end",
    status = "state",
    pdata = births[births$person == "p1", ],
    pvar = c("person", "birth")
  )
)
record_error(
  "pdata_invalid_string",
  seqformat(
    spell_years,
    from = "SPELL",
    to = "STS",
    id = "who",
    begin = "begin",
    end = "end",
    status = "state",
    pdata = "bogus"
  )
)
record_error(
  "spell_before_birth_missing_end",
  seqformat(
    data.frame(
      who = "missing-end",
      begin = 1999L,
      end = NA_integer_,
      state = "A",
      stringsAsFactors = FALSE
    ),
    from = "SPELL",
    to = "STS",
    id = "who",
    begin = "begin",
    end = "end",
    status = "state",
    process = TRUE,
    pdata = data.frame(
      person = "missing-end",
      birth = 2000L,
      stringsAsFactors = FALSE
    ),
    pvar = c("person", "birth"),
    limit = 3
  )
)
record_error(
  "spell_fill_after_missing_end",
  seqformat(
    data.frame(
      who = c("gap", "gap"),
      begin = c(1L, 3L),
      end = c(NA_integer_, 3L),
      state = c("A", "B"),
      stringsAsFactors = FALSE
    ),
    from = "SPELL",
    to = "STS",
    id = "who",
    begin = "begin",
    end = "end",
    status = "state",
    process = TRUE,
    fillblanks = "G",
    limit = 3
  )
)
record_error(
  "pdata_duplicate_id",
  seqformat(
    spell_years,
    from = "SPELL",
    to = "STS",
    id = "who",
    begin = "begin",
    end = "end",
    status = "state",
    process = TRUE,
    pdata = rbind(births, data.frame(person = "p2", birth = 2002L)),
    pvar = c("person", "birth"),
    limit = 5
  )
)
record_error(
  "pdata_missing_id_age2year",
  seqformat(
    spell_ages,
    from = "SPELL",
    to = "STS",
    id = "who",
    begin = "begin",
    end = "end",
    status = "state",
    process = FALSE,
    pdata = births[births$person == "p1", ],
    pvar = c("person", "birth"),
    tmin = 1998,
    tmax = 2005
  )
)

record_error(
  "sps_bad_duration",
  seqformat(matrix("(A,x)", ncol = 1), from = "SPS", to = "STS", stsep = "-")
)
record_error(
  "sps_in_missing_sdsep",
  seqformat(
    matrix(c("(A,1)", "(B,1)"), nrow = 1),
    from = "SPS",
    to = "STS",
    SPS.in = list(xfix = "()")
  )
)
record_error(
  "sps_in_empty",
  seqformat(
    matrix(c("(A,1)", "(B,1)"), nrow = 1),
    from = "SPS",
    to = "STS",
    SPS.in = list()
  )
)
record_error(
  "sps_all_missing",
  seqformat(matrix(NA_character_, nrow = 1, ncol = 2), from = "SPS", to = "STS")
)

record_warnings(
  "deprecated_compressed",
  seqformat(matrix(c("A", "B"), nrow = 1), from = "STS", to = "STS", compressed = FALSE)
)
record_warnings(
  "deprecated_nr",
  seqformat(matrix(c("A", "@"), nrow = 1), from = "STS", to = "STS", nr = "@")
)
record_warnings(
  "pvar_ignored_without_dataframe",
  seqformat(
    spell_years,
    from = "SPELL",
    to = "STS",
    id = "who",
    begin = "begin",
    end = "end",
    status = "state",
    pdata = "auto",
    pvar = c("person", "birth"),
    limit = 5
  )
)
record_warnings(
  "stslist_custom_missing_codes",
  seqformat(sts_object, to = "STS")
)
record_warnings(
  "spell_truncated",
  seqformat(
    spell_process,
    from = "SPELL",
    to = "STS",
    id = "who",
    begin = "lo",
    end = "hi",
    status = "state",
    process = TRUE,
    limit = 6,
    overwrite = TRUE,
    fillblanks = "_",
    missing = "@"
  )
)
record_warnings(
  "spell_missing_birth",
  seqformat(
    spell_years,
    from = "SPELL",
    to = "STS",
    id = "who",
    begin = "begin",
    end = "end",
    status = "state",
    process = TRUE,
    pdata = births_missing,
    pvar = c("person", "birth"),
    limit = 5
  )
)
record_warnings(
  "spell_birth_edges",
  seqformat(
    spell_birth_edges,
    from = "SPELL",
    to = "STS",
    id = "who",
    begin = "begin",
    end = "end",
    status = "state",
    process = TRUE,
    pdata = birth_edges,
    pvar = c("person", "birth"),
    limit = 3
  )
)
record_warnings(
  "spell_no_positive_end",
  seqformat(
    spell_no_positive_end,
    from = "SPELL",
    to = "STS",
    id = "who",
    begin = "begin",
    end = "end",
    status = "state",
    process = TRUE,
    limit = 2
  )
)
record_error(
  "deprecated_compressed_conflict",
  seqformat(
    matrix(c("A", "B"), nrow = 1),
    from = "STS",
    to = "STS",
    compress = TRUE,
    compressed = FALSE
  )
)
record_error(
  "deprecated_nr_conflict",
  seqformat(
    matrix(c("A", "@"), nrow = 1),
    from = "STS",
    to = "STS",
    missing = "@",
    nr = "*"
  )
)

writeLines(
  c(
    "# TraMineR 2.2-14 seqformat input errors and deprecated-alias warnings",
    conditions
  ),
  file.path(out_dir, "seqformat_input_errors.txt")
)

cat("SEQFORMAT_FULL_INPUT_GOLDENS_DONE\n")

# ---------------------------------------------------------------------------
# Task 14: the complete three-input by six-output conversion hub.
#
# The three source tables below describe the same two trajectories:
#
#   p2: A, A, missing, B
#   p1: B, C, C, missing
#
# The missing cells deliberately exercise both an interior nonresponse and a
# trailing structural void once an output converter applies seqprep().
# ---------------------------------------------------------------------------

output_sts <- data.frame(
  group = c("g2", "g1"),
  id = c("subject-2", "subject-1"),
  t1 = c("A", "B"),
  t2 = c("A", "C"),
  t3 = c("@", "C"),
  t4 = c("B", "@"),
  stringsAsFactors = FALSE,
  row.names = c("p2", "p1"),
  check.names = FALSE
)

output_sps <- matrix(
  c(
    "(A,2)", "(@,1)", "(B,1)",
    "(B,1)", "(C,2)", "(@,1)"
  ),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(c("p2", "p1"), paste0("spell", 1:3))
)

output_spell <- data.frame(
  who = c("p2", "p2", "p2", "p1", "p1", "p1"),
  begin = c(1L, 3L, 4L, 1L, 2L, 4L),
  end = c(2L, 3L, 4L, 1L, 3L, 4L),
  state = c("A", "@", "B", "B", "C", "@"),
  stringsAsFactors = FALSE
)

# Raw STS -> TSE intentionally does not replace `missing` before looking up
# transitions. Include the missing token in tevent so that this executable R
# quirk is pinned without triggering a base-R NULL lookup failure.
output_tevent <- matrix(
  NA_character_,
  nrow = 4,
  ncol = 4,
  dimnames = list(c("A", "B", "C", "@"), c("A", "B", "C", "@"))
)
# TraMineR delegates comma splitting to base R's strsplit(), which discards a
# terminal empty component. Keep trailing separators in two live cells so the
# event converter cannot accidentally emit empty event names.
diag(output_tevent) <- c("startA,", "startB", "startC", NA_character_)
output_tevent["A", "B"] <- "AtoB,"
output_tevent["B", "C"] <- "leaveB,enterC"
output_tevent["C", "A"] <- "CtoA"
output_tevent["A", "@"] <- "rawEnterMissing"
output_tevent["@", "B"] <- "rawLeaveMissing"
output_tevent["C", "@"] <- "rawEnterMissingFromC"

output_pdata <- data.frame(
  person = c("birth-p2", "birth-p1"),
  birth = c(1990L, 1980L),
  stringsAsFactors = FALSE
)

output_rows <- list()
output_schema <- character()

schema_vector <- function(value) {
  if (is.null(value)) return("<NONE>")
  paste(as.character(value), collapse = ";")
}

output_cell <- function(value) {
  if (length(value) == 0L || is.na(value)) return("<NA>")
  as.character(value)
}

record_output_case <- function(name, source, target, value, compressed = FALSE) {
  original_class <- paste(class(value), collapse = ";")
  original_rownames <- rownames(value)
  original_colnames <- colnames(value)
  table <- as.data.frame(value, stringsAsFactors = FALSE, optional = TRUE)

  factor_levels <- character()
  column_types <- character()
  for (j in seq_len(ncol(table))) {
    label <- if (is.null(original_colnames)) paste0("[", j, "]") else original_colnames[j]
    column_types <- c(column_types, paste0(label, ":", class(table[[j]])[1]))
    if (is.factor(table[[j]])) {
      factor_levels <- c(
        factor_levels,
        paste0(label, "=", paste(levels(table[[j]]), collapse = ","))
      )
    }
  }
  if (length(factor_levels) == 0L) factor_levels <- "<NONE>"

  output_schema <<- c(
    output_schema,
    paste(
      name,
      source,
      target,
      original_class,
      nrow(table),
      ncol(table),
      schema_vector(original_colnames),
      schema_vector(original_rownames),
      paste(column_types, collapse = ";"),
      paste(factor_levels, collapse = ";"),
      if (compressed) "TRUE" else "FALSE",
      sep = "|"
    )
  )

  rows <- vector("list", nrow(table) * ncol(table))
  k <- 1L
  for (i in seq_len(nrow(table))) {
    for (j in seq_len(ncol(table))) {
      rows[[k]] <- data.frame(
        case = name,
        source = source,
        target = target,
        row = i,
        row_name = if (is.null(original_rownames)) "<NONE>" else original_rownames[i],
        position = j,
        column = if (is.null(original_colnames)) "<NONE>" else original_colnames[j],
        value = output_cell(table[[j]][i]),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      k <- k + 1L
    }
  }
  output_rows[[length(output_rows) + 1L]] <<- do.call(rbind, rows)
}

output_args <- function(source, target) {
  args <- switch(
    source,
    STS = list(
      data = output_sts,
      var = c("t1", "t2", "t3", "t4"),
      from = "STS",
      to = target,
      missing = "@"
    ),
    SPS = list(
      data = output_sps,
      from = "SPS",
      to = target,
      missing = "@"
    ),
    SPELL = list(
      data = output_spell,
      from = "SPELL",
      to = target,
      id = "who",
      begin = "begin",
      end = "end",
      status = "state",
      limit = 4,
      missing = "@"
    )
  )

  if (source == "STS" && target == "STS") args$compress <- TRUE
  if (source == "SPS" && target == "DSS") args$compress <- TRUE
  # DSS intentionally ignores `right` in seqformat.R, while SPS and SPELL
  # receive it. These calls pin both sides of that routing boundary, including
  # TraMineR's named NEUTRAL policy (rendered as state "#").
  if (source == "SPELL" && target %in% c("DSS", "SPELL")) args$right <- NA
  if (source == "SPELL" && target == "SPS") args$right <- "NEUTRAL"

  if (target == "SPS") args$SPS.out <- list(xfix = "{}", sdsep = ";")
  # STS_to_SPS uses only the first and second characters of xfix. A one-char
  # wrapper therefore has a prefix but no suffix.
  if (source == "SPS" && target == "SPS") {
    args$SPS.out <- list(xfix = "[", sdsep = ":")
  }
  if (source == "SPELL" && target == "SPS") args$compress <- TRUE

  if (target == "SRS") args$nrep <- if (source == "SPS") 2 else 1
  if (source == "STS" && target == "SRS") args$covar <- c("group", "id")

  if (target == "TSE") {
    args$tevent <- output_tevent
    if (source == "STS") args$id <- c("alpha", "beta")
    if (source == "SPS") args$id <- c("sps-z", "sps-a")
  }

  if (target == "SPELL") {
    args$with.missing <- source != "STS"
    if (source == "STS") {
      args$pdata <- output_pdata
      args$pvar <- c("person", "birth")
    }
  }
  args
}

output_formats <- c("STS", "DSS", "SPS", "SRS", "TSE", "SPELL")
for (source in c("STS", "SPS", "SPELL")) {
  for (target in output_formats) {
    args <- output_args(source, target)
    value <- suppressWarnings(suppressMessages(do.call(seqformat, args)))
    record_output_case(
      paste0("route_", tolower(source), "_", tolower(target)),
      source,
      target,
      value,
      compressed = isTRUE(args$compress)
    )
  }
}

# SRS covariates are also legal for matrix and stslist sources. Record both
# source classes separately from the 18-route matrix, including the odd but
# executable base-R name assigned to a one-column matrix subset.
matrix_covar_source <- matrix(
  c("g2", "A", "B", "g1", "B", "C"),
  nrow = 2,
  byrow = TRUE
)
record_output_case(
  "matrix_covar",
  "STS_MATRIX",
  "SRS",
  suppressMessages(seqformat(
    matrix_covar_source,
    var = 2:3,
    from = "STS",
    to = "SRS",
    covar = 1
  ))
)

stslist_covar_source <- suppressMessages(seqdef(
  matrix(c("A", "B", "B", "C"), nrow = 2, byrow = TRUE)
))
record_output_case(
  "stslist_covar",
  "STSLIST",
  "SRS",
  suppressWarnings(suppressMessages(seqformat(
    stslist_covar_source,
    from = "STS",
    to = "SRS",
    covar = 1
  )))
)

stopifnot(length(output_rows) == 20L)
write.csv(
  do.call(rbind, output_rows),
  file.path(out_dir, "seqformat_output_cases.csv"),
  row.names = FALSE,
  quote = TRUE,
  na = "<NA>"
)
writeLines(
  c(
    "# TraMineR 2.2-14 seqformat output schema",
    "# case|source|target|class|nrow|ncol|columns|row_names|column_types|factor_levels|compressed",
    output_schema
  ),
  file.path(out_dir, "seqformat_output_schema.txt")
)

output_conditions <- character()
clean_output_condition <- function(x) {
  text <- gsub("[\r\n]+", "\\\\n", conditionMessage(x))
  text <- sub("\\\\n$", "", text)
  gsub("In record_output_(error|warnings)\\(\\)", "In seqformat()", text)
}

record_output_error <- function(name, expr) {
  value <- tryCatch(
    {
      suppressWarnings(suppressMessages(force(expr)))
      "<NO ERROR>"
    },
    error = function(e) clean_output_condition(e)
  )
  output_conditions <<- c(output_conditions, paste("ERROR", name, value, sep = "|"))
}

record_output_warnings <- function(name, expr) {
  found <- character()
  value <- tryCatch(
    withCallingHandlers(
      {
        force(expr)
        NULL
      },
      message = function(m) {
        text <- clean_output_condition(m)
        if (grepl("[!!]", text, fixed = TRUE)) found <<- c(found, text)
        invokeRestart("muffleMessage")
      }
    ),
    error = function(e) {
      found <<- c(found, paste0("ERROR: ", clean_output_condition(e)))
      NULL
    }
  )
  invisible(value)
  if (length(found) == 0L) found <- "<NO WARNING>"
  for (i in seq_along(found)) {
    warning_name <- if (length(found) == 1L) name else paste0(name, "_", i)
    output_conditions <<- c(
      output_conditions,
      paste("WARN", warning_name, found[[i]], sep = "|")
    )
  }
}

output_sts_base <- list(
  data = output_sts,
  var = c("t1", "t2", "t3", "t4"),
  from = "STS",
  missing = "@"
)

record_output_error(
  "missing_tevent",
  do.call(seqformat, c(output_sts_base, list(to = "TSE", id = c("a", "b"))))
)
record_output_error(
  "duplicate_tse_ids",
  do.call(seqformat, c(
    output_sts_base,
    list(to = "TSE", id = c("a", "a"), tevent = output_tevent)
  ))
)
record_output_error(
  "empty_tse_ids",
  do.call(seqformat, c(
    output_sts_base,
    list(to = "TSE", id = character(), tevent = output_tevent)
  ))
)
record_output_error(
  "incomplete_tevent",
  do.call(seqformat, c(
    output_sts_base,
    list(
      to = "TSE",
      id = c("a", "b"),
      tevent = output_tevent[c("A", "B", "C"), c("A", "B", "C")]
    )
  ))
)
record_output_error(
  "compress_srs",
  do.call(seqformat, c(output_sts_base, list(to = "SRS", compress = TRUE)))
)
record_output_error(
  "compress_tse",
  do.call(seqformat, c(
    output_sts_base,
    list(to = "TSE", id = c("a", "b"), tevent = output_tevent, compress = TRUE)
  ))
)
record_output_error(
  "compress_spell",
  do.call(seqformat, c(output_sts_base, list(to = "SPELL", compress = TRUE)))
)
record_output_error(
  "spell_bad_pdata_length",
  do.call(seqformat, c(
    output_sts_base,
    list(
      to = "SPELL",
      pdata = output_pdata[1, , drop = FALSE],
      pvar = c("person", "birth")
    )
  ))
)
record_output_error(
  "bad_covar",
  do.call(seqformat, c(output_sts_base, list(to = "SRS", covar = "notthere")))
)
record_output_error(
  "srs_short_nrep",
  seqformat(
    matrix(c("A", "B", "B", "C"), nrow = 2, byrow = TRUE),
    from = "STS",
    to = "SRS",
    nrep = 1
  )
)
record_output_error(
  "srs_negative_nrep",
  seqformat(
    matrix(c("A", "B", "C", "B", "C", "A"), nrow = 2, byrow = TRUE),
    from = "STS",
    to = "SRS",
    nrep = -1
  )
)

record_output_warnings(
  "covar_ignored_non_srs",
  do.call(seqformat, c(output_sts_base, list(to = "DSS", covar = "group")))
)
record_output_warnings(
  "deprecated_compressed_output",
  do.call(seqformat, c(output_sts_base, list(to = "DSS", compressed = TRUE)))
)
record_output_warnings(
  "deprecated_nr_output",
  seqformat(
    output_sts,
    var = c("t1", "t2", "t3", "t4"),
    from = "STS",
    to = "SPS",
    nr = "@"
  )
)
record_output_warnings(
  "tse_default_ids",
  do.call(seqformat, c(output_sts_base, list(to = "TSE", tevent = output_tevent)))
)
record_output_warnings(
  "spell_pvar_ignored_without_dataframe",
  do.call(seqformat, c(
    output_sts_base,
    list(to = "SPELL", pdata = "auto", pvar = c(1, 2))
  ))
)

writeLines(
  c(
    "# TraMineR 2.2-14 seqformat output errors and warnings",
    output_conditions
  ),
  file.path(out_dir, "seqformat_output_errors.txt")
)

cat("SEQFORMAT_FULL_OUTPUT_GOLDENS_DONE\n")
