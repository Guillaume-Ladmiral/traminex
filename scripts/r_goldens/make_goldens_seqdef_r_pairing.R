# Live-R goldens for how seqdef pairs states, labels and colours with its
# default alphabet (repair of audit finding TX-AUD-20261010-001).
#
# Without `alphabet=`, TraMineR's seqdef pairs `states=`, `labels=` and `cpal=`
# (and its default palette) with the order seqstatl returns: numeric order
# when every state converts with as.numeric, R's collation otherwise
# (seqdef.R:59-82, 110-114, 168-175, 215-232; seqstatl.R:21-30). traminex keeps
# a first-appearance alphabet order, so its tests compare each state's name,
# label and colour, not positions.
#
# The cases `defect_noncanonical_states` and `defect_hex_text_labels` record an R defect: seqdef re-renders the
# numeric seqstatl result with as.character (seqdef.R:65) and then matches
# cells by text, so a cell spelled `01` or `0x10` matches no level and becomes
# NA. Their R output is recorded so that the divergence stays visible; the
# Python test asserts the value derived from the definition instead.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_seqdef_r_pairing.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
stopifnot(capabilities("ICU"))
invisible(icuSetCollate(locale = "root"))
stopifnot(identical(icuGetCollate(), "root"))

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_rows <- list()
cell_rows <- list()
state_rows <- list()

cell_text <- function(value) ifelse(is.na(value), "<NA>", as.character(value))

record <- function(case, input, seq) {
  input <- as.matrix(input)
  input_rows[[length(input_rows) + 1L]] <<- data.frame(
    case = case,
    row = rep(seq_len(nrow(input)), times = ncol(input)),
    col = rep(seq_len(ncol(input)), each = nrow(input)),
    value = cell_text(as.vector(input)),
    stringsAsFactors = FALSE
  )
  cells <- as.matrix(seq)
  cell_rows[[length(cell_rows) + 1L]] <<- data.frame(
    case = case,
    row = rep(seq_len(nrow(cells)), times = ncol(cells)),
    col = rep(seq_len(ncol(cells)), each = nrow(cells)),
    state = cell_text(as.vector(cells)),
    stringsAsFactors = FALSE
  )
  states <- as.character(alphabet(seq))
  state_rows[[length(state_rows) + 1L]] <<- data.frame(
    case = case,
    position = seq_along(states),
    state = states,
    label = as.character(attr(seq, "labels")),
    cpal = as.character(attr(seq, "cpal")),
    stringsAsFactors = FALSE
  )
}

# biofam (the seqdef.Rd idiom): rows 1-20 plus the first row showing each of
# the eight states, so every state of the documented labels is observed.
data(biofam)
bf <- as.matrix(biofam[, 10:25])
first_rows <- vapply(0:7, function(state) which(apply(bf, 1, function(r) any(r == state)))[1], 1L)
keep <- sort(unique(c(1:20, first_rows)))
biofam_input <- bf[keep, , drop = FALSE]
stopifnot(length(unique(as.vector(biofam_input))) == 8L)
biofam_labels <- c("Parent", "Left", "Married", "Left+Marr", "Child", "Left+Child",
                   "Left+Marr+Child", "Divorced")
biofam_states <- c("P", "L", "M", "LM", "C", "LC", "LMC", "D")
record("biofam_labels_states", biofam_input, suppressMessages(TraMineR::seqdef(
  biofam_input, labels = biofam_labels, states = biofam_states
)))

char_input <- matrix(c("C", "A", "B", "A", "B", "C", "C", "A"), nrow = 2, byrow = TRUE)
record("char_states_labels_cpal", char_input, suppressMessages(TraMineR::seqdef(
  char_input, states = c("x", "y", "z"), labels = c("Lx", "Ly", "Lz"),
  cpal = c("red", "green", "blue")
)))
record("char_default", char_input, suppressMessages(TraMineR::seqdef(char_input)))

numeric_text_input <- matrix(c("10", "9", "2", "10"), nrow = 1)
record("numeric_text_states", numeric_text_input, suppressMessages(TraMineR::seqdef(
  numeric_text_input, states = c("two", "nine", "ten")
)))

underscore_input <- matrix(c("1_0", "9", "2"), nrow = 1)
record("underscore_text_labels", underscore_input, suppressMessages(TraMineR::seqdef(
  underscore_input, labels = c("L1", "L2", "L3")
)))

left_input <- matrix(c(NA, "C", "A", "B", "C", "A"), nrow = 2, byrow = TRUE)
record("left_policy_labels_cpal", left_input, suppressMessages(TraMineR::seqdef(
  left_input, left = "LL", labels = c("La", "Lb", "Lc", "Lleft"),
  cpal = c("c1", "c2", "c3", "c4")
)))

# seqstatl sees the cells spelled like `missing` before seqdef drops that code
# (seqdef.R:59-62): a text missing code makes R order numeric states as
# characters, a numeric one keeps numeric order.
missing_text_input <- matrix(c("1", "10", "x", "2", "1", "10"), nrow = 2, byrow = TRUE)
record("missing_text_code_states", missing_text_input, suppressMessages(TraMineR::seqdef(
  missing_text_input, missing = "x", states = c("s1", "s2", "s3"), labels = c("L1", "L2", "L3")
)))
missing_number_input <- matrix(c("1", "10", "99", "2", "1", "10"), nrow = 2, byrow = TRUE)
record("missing_numeric_code_states", missing_number_input, suppressMessages(TraMineR::seqdef(
  missing_number_input, missing = "99", states = c("s1", "s2", "s3"), labels = c("L1", "L2", "L3")
)))

noncanonical_input <- matrix(c("2", "01", "2"), nrow = 1)
record("defect_noncanonical_states", noncanonical_input, suppressMessages(TraMineR::seqdef(
  noncanonical_input, states = c("one", "two")
)))

hex_input <- matrix(c("0x10", "9", "2"), nrow = 1)
record("defect_hex_text_labels", hex_input, suppressMessages(TraMineR::seqdef(
  hex_input, labels = c("L1", "L2", "L3")
)))

# A states vector longer than the observed alphabet, which seqdef.Rd ("Must have
# a length equal to the size of the alphabet") excludes but R accepts: it pairs
# the first entries with the observed states and keeps the rest as unobserved
# states. Recorded to keep traminex's refusal visible as a divergence.
longer_input <- matrix(c("0", "1", "3"), nrow = 1)
record("defect_states_longer", longer_input, suppressMessages(TraMineR::seqdef(
  longer_input, states = c("P", "L", "M", "LM")
)))

write.csv(do.call(rbind, input_rows), file.path(out_dir, "seqdef_pairing_inputs.csv"),
          row.names = FALSE, quote = TRUE)
write.csv(do.call(rbind, cell_rows), file.path(out_dir, "seqdef_pairing_cells.csv"),
          row.names = FALSE, quote = TRUE)
write.csv(do.call(rbind, state_rows), file.path(out_dir, "seqdef_pairing_states.csv"),
          row.names = FALSE, quote = TRUE)
