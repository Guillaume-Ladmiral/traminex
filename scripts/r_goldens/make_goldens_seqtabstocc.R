# Public TraMineRextras 0.6.9 evidence for seqtabstocc.
#
# Run from the repository root:
#   Rscript scripts/r_goldens/make_goldens_seqtabstocc.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
suppressMessages(library(TraMineRextras))

if (!identical(as.character(packageVersion("TraMineR")), "2.2.14")) {
  stop("seqtabstocc goldens require TraMineR 2.2-14")
}
if (!identical(as.character(packageVersion("TraMineRextras")), "0.6.9")) {
  stop("seqtabstocc goldens require TraMineRextras 0.6.9")
}

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

csv_quote <- function(value) {
  paste0('"', gsub('"', '""', as.character(value), fixed = TRUE), '"')
}

format_value <- function(value) {
  if (is.nan(value)) return("NaN")
  if (is.na(value)) return("NA")
  if (is.infinite(value)) return(if (value > 0) "Inf" else "-Inf")
  trimws(formatC(value, format = "g", digits = 15))
}

combo_rows <- list(
  c("A"),
  c("B"),
  c("C"),
  c("D"),
  c("E"),
  c("A", "B"),
  c("A", "C"),
  c("A", "D"),
  c("A", "E"),
  c("B", "C"),
  c("B", "D"),
  c("B", "E"),
  c("A", "B"),
  c("A", "B"),
  c("A", "C")
)
combo <- matrix(NA_character_, nrow = length(combo_rows), ncol = 3)
for (position in seq_along(combo_rows)) {
  combo[position, seq_along(combo_rows[[position]])] <- combo_rows[[position]]
}
sd_combo <- suppressMessages(
  TraMineR::seqdef(
    combo,
    alphabet = c("A", "B", "C", "D", "E"),
    weights = seq_along(combo_rows),
    id = paste0("combo-", seq_along(combo_rows)),
    cnames = c("t1", "t2", "t3")
  )
)

marker_data <- rbind(
  c("A", "*", "B", "*"),
  c("A", "B", "*", "*"),
  c("*", "*", "C", "*"),
  c("C", "*", "*", "*"),
  c("A", "*", "B", "*")
)
sd_markers <- suppressMessages(
  TraMineR::seqdef(
    marker_data,
    alphabet = c("A", "B", "C"),
    missing = "*",
    id = paste0("marker-", seq_len(nrow(marker_data)))
  )
)

locale_states <- c("a", "B", "A", "b", "á", "ä", "Z")
sd_locale <- suppressMessages(
  TraMineR::seqdef(
    cbind(locale_states, rep("X", length(locale_states))),
    alphabet = c(locale_states, "X")
  )
)

internal_order <- rbind(
  c("Z", "ä", "A"),
  c("é", "!", "A"),
  c("ä", "Z", "A")
)
sd_internal <- suppressMessages(
  TraMineR::seqdef(
    internal_order,
    alphabet = c("Z", "ä", "A", "é", "!")
  )
)

rows <- list()
append_case <- function(case, call, result) {
  for (position in seq_len(nrow(result))) {
    rows[[length(rows) + 1L]] <<- paste(
      csv_quote(case),
      csv_quote(call),
      position,
      csv_quote(rownames(result)[[position]]),
      format_value(result$Freq[[position]]),
      format_value(result$Percent[[position]]),
      sep = ","
    )
  }
}

append_case(
  "default_top_ten",
  "TraMineRextras::seqtabstocc(sd_combo)",
  TraMineRextras::seqtabstocc(sd_combo)
)
append_case(
  "all_patterns",
  "TraMineRextras::seqtabstocc(sd_combo,idxs=0)",
  TraMineRextras::seqtabstocc(sd_combo, idxs = 0)
)
append_case(
  "reordered_duplicate_ranks",
  "TraMineRextras::seqtabstocc(sd_combo,idxs=c(2,1,2))",
  TraMineRextras::seqtabstocc(sd_combo, idxs = c(2, 1, 2))
)
append_case(
  "missing_false",
  "TraMineRextras::seqtabstocc(sd_markers,idxs=0,with.missing=FALSE)",
  TraMineRextras::seqtabstocc(
    sd_markers,
    idxs = 0,
    with.missing = FALSE
  )
)
append_case(
  "missing_true",
  "TraMineRextras::seqtabstocc(sd_markers,idxs=0,with.missing=TRUE)",
  TraMineRextras::seqtabstocc(
    sd_markers,
    idxs = 0,
    with.missing = TRUE
  )
)
append_case(
  "locale_ties",
  "TraMineRextras::seqtabstocc(sd_locale,idxs=0)",
  TraMineRextras::seqtabstocc(sd_locale, idxs = 0)
)
append_case(
  "internal_state_order",
  "TraMineRextras::seqtabstocc(sd_internal,idxs=0)",
  TraMineRextras::seqtabstocc(sd_internal, idxs = 0)
)
append_case(
  "weighted_true_lost",
  "TraMineRextras::seqtabstocc(sd_combo,idxs=0,weighted=TRUE)",
  TraMineRextras::seqtabstocc(sd_combo, idxs = 0, weighted = TRUE)
)
append_case(
  "weighted_false_lost",
  "TraMineRextras::seqtabstocc(sd_combo,idxs=0,weighted=FALSE)",
  TraMineRextras::seqtabstocc(sd_combo, idxs = 0, weighted = FALSE)
)

captured_warnings <- capture.output(
  tlim_result <- TraMineRextras::seqtabstocc(
    sd_combo,
    tlim = c(2, 1)
  ),
  type = "message"
)
append_case(
  "deprecated_tlim",
  "TraMineRextras::seqtabstocc(sd_combo,tlim=c(2,1))",
  tlim_result
)

manifest <- file.path(out_dir, "seqtabstocc_public_cases.csv")
con <- file(manifest, open = "wt", encoding = "UTF-8")
writeLines(
  "case,call,position,r_pattern,Freq,Percent",
  con = con,
  sep = "\n",
  useBytes = TRUE
)
writeLines(
  unlist(rows, use.names = FALSE),
  con = con,
  sep = "\n",
  useBytes = TRUE
)
close(con)

writeLines(
  c(
    "case|message",
    paste("deprecated_tlim", captured_warnings, sep = "|")
  ),
  file.path(out_dir, "seqtabstocc_public_warnings.txt"),
  useBytes = TRUE
)

capture_error <- function(case, call, expression) {
  condition <- tryCatch(
    {
      suppressWarnings(force(expression))
      NULL
    },
    error = function(error) error
  )
  if (is.null(condition)) stop("expected an error for ", case)
  paste(
    case,
    call,
    class(condition)[1L],
    gsub("[\r\n]+", " ", conditionMessage(condition)),
    sep = "|"
  )
}

errors <- c(
  "case|call|condition|message",
  capture_error(
    "format_conflict",
    "TraMineRextras::seqtabstocc(sd_combo,format='STS')",
    TraMineRextras::seqtabstocc(sd_combo, format = "STS")
  ),
  capture_error(
    "unknown_forwarded",
    "TraMineRextras::seqtabstocc(sd_combo,banana=1)",
    TraMineRextras::seqtabstocc(sd_combo, banana = 1)
  ),
  capture_error(
    "idxs_tlim_conflict",
    "TraMineRextras::seqtabstocc(sd_combo,idxs=c(2,1),tlim=1)",
    TraMineRextras::seqtabstocc(sd_combo, idxs = c(2, 1), tlim = 1)
  )
)
writeLines(
  errors,
  file.path(out_dir, "seqtabstocc_public_errors.txt"),
  useBytes = TRUE
)

cat("SEQTABSTOCC_GOLDENS_DONE\n")
