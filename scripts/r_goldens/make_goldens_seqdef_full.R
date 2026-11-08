# Golden generation for the complete TraMineR 2.2-14 seqdef contract.
#
# Run from the repository root with:
#   /usr/bin/Rscript scripts/r_goldens/make_goldens_seqdef_full.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cell_rows <- list()
metadata_rows <- c("# TraMineR 2.2-14 seqdef metadata")
error_rows <- c("# TraMineR 2.2-14 seqdef errors")

scalar_text <- function(value) {
  if (length(value) == 0L || is.null(value)) return("<NULL>")
  if (is.na(value)) return("<NA>")
  if (is.logical(value)) return(if (value) "TRUE" else "FALSE")
  as.character(value)
}

vector_text <- function(value) {
  if (is.null(value)) return("<NULL>")
  paste(vapply(unname(value), scalar_text, character(1)), collapse = ";")
}

record_metadata <- function(case, field, value) {
  metadata_rows <<- c(metadata_rows, paste(case, field, vector_text(value), sep = "|"))
}

# traminex's stslist folds the missing (nr) and void markers into `alphabet`,
# which R keeps as separate attributes. Its real states are R's alphabet, in
# R's order: given `alphabet=` (mapped through `states=`) or R's seqstatl
# order without it (owner decision D2 of 2026-10-19; until then traminex
# kept first appearance and this function derived that order from the
# decoded matrix). The markers follow, nonresponse first.
#
# Nothing here is authored in Python; `alphabet_given` states a fact about the
# seqdef call, not an expected value.
python_alphabet <- function(value, alphabet_given) {
  nr <- attr(value, "nr")
  void <- attr(value, "void")
  r_states <- as.character(attr(value, "alphabet"))

  result <- r_states

  cells <- as.character(as.matrix(value))
  if (any(cells == nr, na.rm = TRUE) && !nr %in% result) result <- c(result, nr)
  if (any(cells == void, na.rm = TRUE) && !void %in% result) result <- c(result, void)
  result
}

record_case <- function(name, value, alphabet_given = FALSE) {
  value_matrix <- as.matrix(value)
  ids <- rownames(value)
  if (is.null(ids)) ids <- as.character(seq_len(nrow(value_matrix)))
  columns <- colnames(value)
  if (is.null(columns)) columns <- as.character(seq_len(ncol(value_matrix)))

  rows <- vector("list", nrow(value_matrix) * ncol(value_matrix))
  k <- 1L
  for (i in seq_len(nrow(value_matrix))) {
    for (j in seq_len(ncol(value_matrix))) {
      rows[[k]] <- data.frame(
        case = name,
        row = i,
        id = ids[i],
        position = j,
        column = columns[j],
        value = scalar_text(value_matrix[i, j]),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  cell_rows[[length(cell_rows) + 1L]] <<- do.call(rbind, rows)

  record_metadata(name, "r_alphabet", attr(value, "alphabet"))
  record_metadata(name, "python_alphabet", python_alphabet(value, alphabet_given))
  record_metadata(name, "labels", attr(value, "labels"))
  record_metadata(name, "cpal", attr(value, "cpal"))
  record_metadata(name, "ids", rownames(value))
  record_metadata(name, "weights", attr(value, "weights"))
  record_metadata(name, "start", attr(value, "start"))
  record_metadata(name, "missing", attr(value, "missing"))
  record_metadata(name, "void", attr(value, "void"))
  record_metadata(name, "nr", attr(value, "nr"))
  record_metadata(name, "xtlab", colnames(value))
  record_metadata(name, "xtstep", attr(value, "xtstep"))
  record_metadata(name, "tick_last", attr(value, "tick.last"))
  record_metadata(name, "missing_color", attr(value, "missing.color"))
}

record_error <- function(name, expression) {
  message <- tryCatch(
    {
      suppressWarnings(suppressMessages(force(expression)))
      "<NO ERROR>"
    },
    error = function(condition) conditionMessage(condition)
  )
  escaped <- paste(strsplit(message, "\n", fixed = TRUE)[[1]], collapse = "\\n")
  error_rows <<- c(error_rows, paste("ERROR", name, escaped, sep = "|"))
}

# Default policies: left and interior missings remain nonresponses, while
# right missings are deleted into structural voids and moved to the row end.
sts_default_input <- data.frame(
  t1 = c("B", NA, "A"),
  t2 = c(NA, "A", NA),
  t3 = c("A", NA, NA),
  stringsAsFactors = FALSE,
  row.names = c("default-b", "default-a", "default-tail")
)
sts_default_policies <- suppressMessages(seqdef(sts_default_input))
record_case("sts_default_policies", sts_default_policies)

# Applying DEL on every missing region compacts all observed states left and
# leaves only structural voids at the right edge.
sts_all_del <- suppressMessages(seqdef(
  sts_default_input,
  left = "DEL",
  gaps = "DEL",
  right = "DEL",
  void = "~",
  nr = "?"
))
record_case("sts_all_del", sts_all_del)

# The complete replacement/metadata matrix.  Explicit alphabet order and the
# positional states mapping are both observable, as are all three replacement
# policies and every display/identity attribute.
sts_replacement_input <- data.frame(
  ignored = c("x", "y", "z"),
  w1 = c("B", "MISS", "A"),
  w2 = c("MISS", "A", "MISS"),
  w3 = c("A", "MISS", "C"),
  w4 = c("MISS", "C", "MISS"),
  stringsAsFactors = FALSE,
  row.names = c("source-1", "source-2", "source-3")
)
sts_replacement_full_meta <- suppressMessages(seqdef(
  sts_replacement_input,
  var = c("w1", "w2", "w3", "w4"),
  alphabet = c("B", "A", "C"),
  states = c("Beta", "Alpha", "Gamma"),
  id = c("case-z", "case-a", "case-m"),
  weights = c(0.5, 1.25, 2.25),
  start = 1999,
  left = "LEFT",
  gaps = "GAP",
  right = "RIGHT",
  missing = "MISS",
  void = "~",
  nr = "?",
  cnames = c("Y0", "Y1", "Y2", "Y3"),
  xtstep = 2,
  tick.last = TRUE,
  cpal = c("#110000", "#001100", "#000011", "#aaaa00", "#00aaaa", "#aa00aa"),
  missing.color = "#c0ffee",
  labels = c("State B", "State A", "State C", "Before", "Between", "After")
))
record_case("sts_replacement_full_meta", sts_replacement_full_meta, alphabet_given = TRUE)

# Equal-width compressed STS with an explicitly supplied non-default separator.
sts_compressed_input <- matrix(
  c("B;A;B", "A;C;C"),
  ncol = 1,
  dimnames = list(c("compressed-b", "compressed-a"), "sequence")
)
sts_compressed <- suppressMessages(seqdef(sts_compressed_input, stsep = ";"))
record_case("sts_compressed", sts_compressed)

# Reordered explicit alphabet, an unobserved state, and positional state labels.
states_mapping_input <- matrix(
  c("b", "a", "a", "b"),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(c("mapping-1", "mapping-2"), c("first", "second"))
)
states_mapping <- suppressMessages(seqdef(
  states_mapping_input,
  alphabet = c("b", "a", "unused"),
  states = c("Bee", "Ay", "Unused")
))
record_case("states_mapping", states_mapping, alphabet_given = TRUE)

# SPS and SPELL inputs exercise seqdef's informat dispatch before seqprep.
sps_default_input <- matrix(
  c("(A,2)", "(B,1)", "(C,1)", "(@,2)"),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(c("sps-1", "sps-2"), c("spell1", "spell2"))
)
sps_default <- suppressMessages(seqdef(
  sps_default_input,
  informat = "SPS",
  missing = "@"
))
record_case("sps_default", sps_default)

spell_limit4_input <- data.frame(
  who = c("b", "b", "a"),
  begin = c(1L, 3L, 2L),
  end = c(2L, 6L, 5L),
  state = c("A", "B", "C"),
  stringsAsFactors = FALSE
)
spell_limit4 <- suppressWarnings(suppressMessages(seqdef(
  spell_limit4_input,
  informat = "SPELL",
  process = TRUE,
  limit = 4,
  overwrite = TRUE,
  fillblanks = "_"
)))
record_case("spell_limit4", spell_limit4)

# Default palettes are core seqdef metadata, independent of plotting imports.
for (size in c(1L, 2L, 3L, 8L, 9L, 12L, 13L)) {
  states <- paste0("s", seq_len(size))
  palette_case <- suppressMessages(seqdef(
    matrix(states, nrow = 1L),
    alphabet = states
  ))
  record_metadata(paste0("palette_", size), "cpal", attr(palette_case, "cpal"))
}

# Validation and base-row-name failures observable through seqdef.
validation_input <- matrix(
  c("A", "B", "B", "A"),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(c("validation-1", "validation-2"), c("t1", "t2"))
)
record_error(
  "duplicate_alphabet",
  seqdef(validation_input, alphabet = c("A", "A", "B"))
)
record_error(
  "unknown_state",
  seqdef(validation_input, alphabet = c("A"))
)
record_error(
  "weights_length",
  seqdef(validation_input, weights = 1)
)
record_error(
  "weights_na",
  seqdef(validation_input, weights = c(1, NA_real_))
)
record_error(
  "weights_negative",
  seqdef(validation_input, weights = c(1, -1))
)
record_error(
  "palette_length",
  seqdef(validation_input, cpal = "red")
)
record_error(
  "label_length",
  seqdef(validation_input, labels = "Alpha")
)
record_error(
  "void_type",
  seqdef(validation_input, void = 1)
)
record_error(
  "void_equals_left",
  seqdef(validation_input, void = "X", left = "X")
)
record_error(
  "void_equals_gaps",
  seqdef(validation_input, void = "X", gaps = "X")
)
record_error(
  "void_equals_right",
  seqdef(validation_input, void = "X", right = "X")
)
record_error(
  "id_length",
  seqdef(validation_input, id = "only-one")
)
record_error(
  "id_duplicate",
  seqdef(validation_input, id = c("same", "same"))
)
record_error(
  "id_missing",
  seqdef(validation_input, id = c("one", NA_character_))
)

cells <- do.call(rbind, cell_rows)
write.csv(
  cells,
  file.path(out_dir, "seqdef_full_cells.csv"),
  row.names = FALSE,
  na = "<NA>",
  quote = TRUE
)
writeLines(metadata_rows, file.path(out_dir, "seqdef_full_metadata.txt"), useBytes = TRUE)
writeLines(error_rows, file.path(out_dir, "seqdef_full_errors.txt"), useBytes = TRUE)
