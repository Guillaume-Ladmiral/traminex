# Golden generation: seqtab / seqtabstocc (T5) -- sequence frequency tables.
#
# seqtab pins the full public result contract: ranks, labels, weights,
# representatives, and metadata. The legacy numeric matrices remain for
# seqtabstocc and as compact compatibility evidence for fixture S.
#
# seqtabstocc pins the state co-occurrence frequencies (frequency of the set of
# distinct states each sequence visits) on the shared A/B fixtures.
#
# Run from the repository root with the oracle library (scripts/r_goldens/README.md):
#   Rscript scripts/r_goldens/make_goldens_seqtab.R
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
source(file.path("scripts", "r_goldens", "make_goldens_corestats.R"))
suppressMessages(library(TraMineRextras))

if (!identical(as.character(packageVersion("TraMineR")), "2.2.14")) {
  stop("seqtab goldens require TraMineR 2.2-14")
}
if (!identical(as.character(packageVersion("TraMineRextras")), "0.6.9")) {
  stop("seqtabstocc goldens require TraMineRextras 0.6.9")
}

# Fixture S: repeated sequences for seqtab counting/sorting/percent.
S <- rbind(
  c("A", "A", "B", "C"),
  c("A", "A", "B", "C"),
  c("A", "B", "B", "C"),
  c("A", "A", "B", "C"),
  c("A", "B", "B", "C")
)
wS <- c(1, 2, 1, 3, 1)
sdS <- seqdef(S, alphabet = c("A", "B", "C"))
sdSw <- seqdef(S, alphabet = c("A", "B", "C"), weights = wS)

emit_seqtab <- function(sd, tag) {
  tf <- TraMineR::seqtab(sd, idxs = 0)
  fr <- as.matrix(attr(tf, "freq")) # cols [Freq, Percent], sorted desc
  write_num_matrix(
    paste0("seqtab_", tag), fr,
    sprintf("# R: seqtab(%s, idxs=0) attr 'freq' = [Freq, Percent], sorted desc\n", tag)
  )
}

emit_seqtab(sdS, "S")
emit_seqtab(sdSw, "Sw")

# Asymmetric public-contract fixture: 12 unique patterns, two duplicated
# patterns, non-unit weights, an interior missing state, and trailing void.
T <- rbind(
  c("A", "A", "A", "A"),
  c("A", "A", "A", "B"),
  c("A", "A", "B", "B"),
  c("A", "B", "B", "B"),
  c("B", "B", "B", "B"),
  c("B", "B", "B", "C"),
  c("B", "B", "C", "C"),
  c("B", "C", "C", "C"),
  c("C", "C", "C", "C"),
  c("C", "C", "C", "A"),
  c("A", "*", "B", "*"),
  c("C", "A", "B", "*"),
  c("A", "A", "A", "A"),
  c("A", "A", "A", "B")
)
wT <- c(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 2)
sdT <- suppressMessages(
  TraMineR::seqdef(T, alphabet = c("A", "B", "C"), missing = "*")
)
sdTw <- suppressMessages(
  TraMineR::seqdef(
    T,
    alphabet = c("A", "B", "C"),
    missing = "*",
    weights = wT
  )
)

# Equal-frequency labels that distinguish ICU root collation from Python
# code-point ordering. Two columns avoid TraMineR's one-column subsetting drop.
V.states <- c("a", "B", "A", "b", "á", "ä", "Z")
V <- cbind(V.states, V.states)
sdV <- suppressMessages(
  TraMineR::seqdef(V, alphabet = V.states)
)

# Separate literal all-nr branch. right=NA retains terminal missing values as
# nr rather than converting them to structural void.
U <- rbind(
  c("*", "*", "*"),
  c("A", "A", "A"),
  c("A", "A", "A")
)
sdU <- suppressMessages(
  TraMineR::seqdef(
    U,
    alphabet = c("A", "B"),
    missing = "*",
    right = NA
  )
)

format_value <- function(value) {
  if (is.nan(value)) return("NaN")
  if (is.na(value)) return("NA")
  if (is.infinite(value)) return(if (value > 0) "Inf" else "-Inf")
  trimws(formatC(value, format = "g", digits = 15))
}

csv_quote <- function(value) {
  paste0('"', gsub('"', '""', as.character(value), fixed = TRUE), '"')
}

append_public_case <- function(rows, case, call, result, source_map) {
  freq <- attr(result, "freq")
  idxs <- attr(result, "idxs")
  idxf <- attr(result, "idxf")
  representative_sts <- as.vector(
    TraMineR::seqconc(result, void = attr(result, "void"))
  )
  display_sequences <- if (attr(result, "format") == "SPS") {
    suppressMessages(
      TraMineR::seqformat(
        result,
        from = "STS",
        to = "SPS",
        SPS.out = list(xfix = "", sdsep = "/"),
        compress = TRUE
      )
    )
  } else {
    representative_sts
  }
  representative_weights <- attr(result, "weights")
  stopifnot(
    nrow(freq) == nrow(result),
    length(idxs) == nrow(result),
    length(idxf) == nrow(result),
    length(representative_sts) == nrow(result),
    length(representative_weights) == nrow(result)
  )

  for (position in seq_len(nrow(result))) {
    rows[[length(rows) + 1L]] <- c(
      csv_quote(case),
      csv_quote(call),
      position,
      csv_quote(display_sequences[[position]]),
      csv_quote(rownames(freq)[[position]]),
      format_value(freq$Freq[[position]]),
      format_value(freq$Percent[[position]]),
      idxs[[position]],
      idxf[[position]],
      source_map[[idxf[[position]]]],
      csv_quote(attr(result, "format")),
      format_value(attr(result, "nbseq")),
      csv_quote(attr(result, "weighted")),
      csv_quote(representative_sts[[position]]),
      format_value(representative_weights[[position]])
    )
  }
  rows
}

capture_warning <- function(expression) {
  value <- NULL
  messages <- capture.output(
    value <- eval.parent(substitute(expression)),
    type = "message"
  )
  list(value = value, messages = messages)
}

public_rows <- list()
source_T <- seq_len(nrow(sdT))
public_rows <- append_public_case(
  public_rows,
  "default",
  "TraMineR::seqtab(sdT)",
  TraMineR::seqtab(sdT),
  source_T
)
public_rows <- append_public_case(
  public_rows,
  "all",
  "TraMineR::seqtab(sdT,idxs=0)",
  TraMineR::seqtab(sdT, idxs = 0),
  source_T
)
public_rows <- append_public_case(
  public_rows,
  "scalar_rank",
  "TraMineR::seqtab(sdT,idxs=2)",
  TraMineR::seqtab(sdT, idxs = 2),
  source_T
)
public_rows <- append_public_case(
  public_rows,
  "reversed_ranks",
  "TraMineR::seqtab(sdT,idxs=c(2,1))",
  TraMineR::seqtab(sdT, idxs = c(2, 1)),
  source_T
)
public_rows <- append_public_case(
  public_rows,
  "duplicate_ranks",
  "TraMineR::seqtab(sdT,idxs=c(2,1,2))",
  TraMineR::seqtab(sdT, idxs = c(2, 1, 2)),
  source_T
)
public_rows <- append_public_case(
  public_rows,
  "out_of_range",
  "TraMineR::seqtab(sdT,idxs=20)",
  TraMineR::seqtab(sdT, idxs = 20),
  source_T
)
public_rows <- append_public_case(
  public_rows,
  "weighted_all",
  "TraMineR::seqtab(sdTw,idxs=0,weighted=TRUE)",
  TraMineR::seqtab(sdTw, idxs = 0, weighted = TRUE),
  source_T
)
public_rows <- append_public_case(
  public_rows,
  "weighted_ignored",
  "TraMineR::seqtab(sdTw,idxs=0,weighted=FALSE)",
  TraMineR::seqtab(sdTw, idxs = 0, weighted = FALSE),
  source_T
)
public_rows <- append_public_case(
  public_rows,
  "sts_all",
  'TraMineR::seqtab(sdT,idxs=0,format="STS")',
  TraMineR::seqtab(sdT, idxs = 0, format = "STS"),
  source_T
)
public_rows <- append_public_case(
  public_rows,
  "locale_ties",
  'TraMineR::seqtab(sdV,idxs=0,weighted=FALSE,format="STS")',
  TraMineR::seqtab(sdV, idxs = 0, weighted = FALSE, format = "STS"),
  seq_len(nrow(sdV))
)

source_U <- which(rowSums(sdU != attr(sdU, "nr")) != 0)
public_rows <- append_public_case(
  public_rows,
  "all_nr_removed",
  "TraMineR::seqtab(sdU,idxs=0)",
  TraMineR::seqtab(sdU, idxs = 0),
  source_U
)

alias_only <- capture_warning(TraMineR::seqtab(sdT, tlim = c(2, 1)))
alias_default <- capture_warning(
  TraMineR::seqtab(sdT, idxs = 1:10, tlim = c(2, 1))
)
public_rows <- append_public_case(
  public_rows,
  "tlim_only",
  "TraMineR::seqtab(sdT,tlim=c(2,1))",
  alias_only$value,
  source_T
)
public_rows <- append_public_case(
  public_rows,
  "tlim_with_default",
  "TraMineR::seqtab(sdT,idxs=1:10,tlim=c(2,1))",
  alias_default$value,
  source_T
)

public_path <- file.path(out_dir, "seqtab_public_cases.csv")
public_con <- file(public_path, open = "wt", encoding = "UTF-8")
writeLines(
  paste(
    "case",
    "call",
    "position",
    "sequence",
    "r_freq_rowname",
    "Freq",
    "Percent",
    "idxs",
    "idxf",
    "source_index",
    "format",
    "nbseq",
    "weighted",
    "representative_sts",
    "representative_weight",
    sep = ","
  ),
  con = public_con,
  sep = "\n",
  useBytes = TRUE
)
writeLines(
  vapply(public_rows, paste, collapse = ",", FUN.VALUE = character(1)),
  con = public_con,
  sep = "\n",
  useBytes = TRUE
)
close(public_con)

warning_lines <- c(
  "case|message",
  paste("tlim_only", alias_only$messages, sep = "|"),
  paste("tlim_with_default", alias_default$messages, sep = "|")
)
writeLines(
  warning_lines,
  con = file.path(out_dir, "seqtab_public_warnings.txt"),
  sep = "\n",
  useBytes = TRUE
)

capture_error <- function(case, expression) {
  condition <- tryCatch(
    {
      suppressWarnings(eval.parent(substitute(expression)))
      NULL
    },
    error = function(error) error
  )
  if (is.null(condition)) stop("expected R error for case: ", case)
  paste(
    case,
    paste(class(condition), collapse = ","),
    gsub("[\r\n]+", " ", conditionMessage(condition)),
    sep = "|"
  )
}

error_lines <- c(
  "case|class|message",
  capture_error("idxs_null", TraMineR::seqtab(sdT, idxs = NULL)),
  capture_error("idxs_negative", TraMineR::seqtab(sdT, idxs = -1)),
  capture_error("idxs_mixed_zero", TraMineR::seqtab(sdT, idxs = c(0, 1))),
  capture_error(
    "idxs_tlim_conflict",
    TraMineR::seqtab(sdT, idxs = c(2, 1), tlim = 1)
  ),
  capture_error("invalid_format", TraMineR::seqtab(sdT, format = "sps"))
)
writeLines(
  error_lines,
  con = file.path(out_dir, "seqtab_public_errors.txt"),
  sep = "\n",
  useBytes = TRUE
)

emit_stocc <- function(sd, tag) {
  t <- TraMineRextras::seqtabstocc(sd, with.missing = FALSE)
  write_num_matrix(
    paste0("seqtabstocc_", tag), as.matrix(t),
    sprintf("# R: seqtabstocc(%s, with.missing=FALSE) = [Freq, Percent] of distinct-state sets\n", tag)
  )
}

fx <- corestats_fixtures()
emit_stocc(fx$A, "A")
emit_stocc(fx$B, "B")

cat("SEQTAB_GOLDENS_DONE\n")
