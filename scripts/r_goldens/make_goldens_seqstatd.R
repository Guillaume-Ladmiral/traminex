# Golden generation: complete R2 seqstatd contract.
#
# Pins TraMineR 2.2-14 values and metadata for weighting, missing/void states,
# normalized/raw entropy, declared-but-unobserved missing, and empty positions.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
source(file.path("scripts", "r_goldens", "make_goldens_corestats.R"))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

format_num <- function(value) {
  trimws(formatC(value, format = "g", digits = 15))
}

case_rows <- list()
metadata_rows <- list()

record_metadata <- function(case, result) {
  values <- list(
    nbseq = format_num(attr(result, "nbseq")),
    weighted = as.character(attr(result, "weighted")),
    norm = as.character(attr(result, "norm")),
    xtlab = paste(attr(result, "xtlab"), collapse = "|"),
    xtstep = as.character(attr(result, "xtstep")),
    `tick.last` = as.character(attr(result, "tick.last")),
    cpal = paste(attr(result, "cpal"), collapse = "|")
  )
  metadata_rows[[length(metadata_rows) + 1L]] <<- data.frame(
    case = case,
    field = names(values),
    value = unlist(values, use.names = FALSE),
    stringsAsFactors = FALSE
  )
}

record_case <- function(case, seqdata, weighted, with_missing, norm) {
  result <- seqstatd(
    seqdata,
    weighted = weighted,
    with.missing = with_missing,
    norm = norm
  )
  states <- rownames(result$Frequencies)
  times <- colnames(result$Frequencies)
  if (length(times) > 0L) {
    rows <- vector("list", length(times) * length(states))
    index <- 1L
    for (time_index in seq_along(times)) {
      for (state_index in seq_along(states)) {
        rows[[index]] <- data.frame(
          case = case,
          time_index = time_index,
          time_label = times[[time_index]],
          state = states[[state_index]],
          value = format_num(result$Frequencies[state_index, time_index]),
          valid_states = format_num(result$ValidStates[[time_index]]),
          entropy = format_num(result$Entropy[[time_index]]),
          stringsAsFactors = FALSE
        )
        index <- index + 1L
      }
    }
    case_rows[[length(case_rows) + 1L]] <<- do.call(rbind, rows)
  }
  record_metadata(case, result)
}

invisible(capture.output(fx <- corestats_fixtures()))

# Retain the established T3 artifacts while adding the R2 long-form manifest.
emit_legacy <- function(seqdata, tag, with_missing) {
  result <- seqstatd(
    seqdata,
    weighted = TRUE,
    with.missing = with_missing,
    norm = TRUE
  )
  header <- sprintf(
    "# R: seqstatd(%s, weighted=TRUE, with.missing=%s, norm=TRUE)$Frequencies, t() -> L x K\n",
    tag,
    with_missing
  )
  write_num_matrix(
    paste0("seqstatd_", tag, "_freq"),
    t(result$Frequencies),
    header
  )
  valid_entropy <- rbind(
    as.vector(result$ValidStates),
    as.vector(result$Entropy)
  )
  header <- sprintf(
    "# R: seqstatd(%s, ..., with.missing=%s) rows=[ValidStates, Entropy]\n",
    tag,
    with_missing
  )
  write_num_matrix(
    paste0("seqstatd_", tag, "_ve"),
    valid_entropy,
    header
  )
}

emit_legacy(fx$A, "A", FALSE)
emit_legacy(fx$Aw, "Aw", FALSE)
emit_legacy(fx$B, "B_wm0", FALSE)
emit_legacy(fx$B, "B_wm1", TRUE)
emit_legacy(fx$C, "C_wm0", FALSE)
emit_legacy(fx$C, "C_wm1", TRUE)

invisible(capture.output(
  declared_nr_absent <- seqdef(
    rbind(c("A", "B"), c("B", "A")),
    alphabet = c("A", "B"),
    missing = "*",
    right = NA
  )
))

all_void <- fx$C
all_void[,] <- attr(all_void, "void")

# A zero-column stslist can be represented, but upstream seqstatd iterates
# `1:seql` and therefore refuses it. Record that oracle boundary separately;
# Python retains its safe empty-array behavior.
zero_width <- data.frame(row.names = rownames(fx$A))
for (name in setdiff(names(attributes(fx$A)), c("names", "row.names", "class"))) {
  attr(zero_width, name) <- attr(fx$A, name)
}
class(zero_width) <- class(fx$A)

record_case("A_unweighted_norm", fx$Aw, FALSE, FALSE, TRUE)
record_case("A_unweighted_raw_entropy", fx$Aw, FALSE, FALSE, FALSE)
record_case("Aw_weighted_norm", fx$Aw, TRUE, FALSE, TRUE)
record_case("Aw_weighted_raw_entropy", fx$Aw, TRUE, FALSE, FALSE)
record_case("B_wm0", fx$B, TRUE, FALSE, TRUE)
record_case("B_wm1", fx$B, TRUE, TRUE, TRUE)
record_case("C_ragged_wm0", fx$C, TRUE, FALSE, TRUE)
record_case("C_ragged_wm1", fx$C, TRUE, TRUE, TRUE)
record_case("declared_nr_absent", declared_nr_absent, TRUE, TRUE, TRUE)
record_case("all_void", all_void, TRUE, FALSE, TRUE)

zero_width_error <- tryCatch(
  {
    seqstatd(zero_width)
    ""
  },
  error = function(condition) conditionMessage(condition)
)
metadata_rows[[length(metadata_rows) + 1L]] <- data.frame(
  case = "zero_width",
  field = "error",
  value = zero_width_error,
  stringsAsFactors = FALSE
)

write.table(
  do.call(rbind, case_rows),
  file.path(out_dir, "seqstatd_r2_cases.csv"),
  sep = ",",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  na = "NA"
)
write.table(
  do.call(rbind, metadata_rows),
  file.path(out_dir, "seqstatd_r2_metadata.txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  na = "NA"
)
cat("SEQSTATD_R2_GOLDENS_DONE\n")
