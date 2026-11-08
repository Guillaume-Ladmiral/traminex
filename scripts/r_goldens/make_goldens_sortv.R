# Golden generation for TraMineRextras 0.6.9 sorti/sortv.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineRextras))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
stopifnot(as.character(packageVersion("TraMineRextras")) == "0.6.9")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

make_sts <- function(codes, states) {
  codes <- as.matrix(codes)
  levels <- c(states, "*", "%")
  if (ncol(codes) == 0) {
    frame <- data.frame(row.names = seq_len(nrow(codes)))
  } else {
    frame <- as.data.frame(
      lapply(seq_len(ncol(codes)), function(column) {
        factor(codes[, column], levels = seq_along(levels), labels = levels)
      }),
      optional = TRUE
    )
  }
  if (ncol(codes) > 0) {
    names(frame) <- paste0("T", seq_len(ncol(codes)))
  }
  class(frame) <- c("stslist", "data.frame")
  attr(frame, "start") <- 1
  attr(frame, "missing") <- NA
  attr(frame, "void") <- "%"
  attr(frame, "nr") <- "*"
  attr(frame, "alphabet") <- states
  attr(frame, "labels") <- states
  attr(frame, "cpal") <- rep("black", length(states))
  attr(frame, "xtstep") <- 1
  attr(frame, "tick.last") <- FALSE
  frame
}

cases <- list(
  declared_order = make_sts(
    matrix(c(
      2, 3, 1,
      3, 2, 1,
      2, 3, 3,
      3, 2, 3,
      2, 3, 1
    ), byrow = TRUE, ncol = 3),
    c("C", "B", "A")
  ),
  ragged_void = make_sts(
    matrix(c(
      1, 2, 4,
      1, 3, 4,
      2, 1, 2,
      4, 4, 4
    ), byrow = TRUE, ncol = 3),
    c("A", "B")
  ),
  all_void = make_sts(matrix(4, nrow = 3, ncol = 2), c("A", "B"))
)

records <- list()
record_index <- 1L
for (case_name in names(cases)) {
  seqdata <- cases[[case_name]]
  for (start in c("beg", "end")) {
    values <- list(
      sorti = sorti(seqdata, start = start, sort.index = TRUE),
      inverse = sorti(seqdata, start = start, sort.index = FALSE),
      sortv = sortv(seqdata, start = start)
    )
    for (result_type in names(values)) {
      for (position in seq_along(values[[result_type]])) {
        records[[record_index]] <- data.frame(
          case = case_name,
          start = start,
          result_type = result_type,
          position = position,
          value = values[[result_type]][[position]],
          stringsAsFactors = FALSE
        )
        record_index <- record_index + 1L
      }
    }
  }
}
write.csv(
  do.call(rbind, records),
  file.path(out_dir, "sortv_r3_cases.csv"),
  row.names = FALSE,
  quote = TRUE
)

condition_lines <- character()
capture_condition <- function(case_name, expression) {
  tryCatch(
    withCallingHandlers(
      force(expression),
      warning = function(condition) {
        condition_lines <<- c(
          condition_lines,
          paste(
            case_name,
            "warning",
            sub("[[:space:]]+$", "", conditionMessage(condition)),
            sep = "|"
          )
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

invalid <- capture_condition("invalid_start", sorti(cases$declared_order, start = "middle"))
condition_lines <- c(
  condition_lines,
  paste("invalid_start", "result", paste(invalid, collapse = ","), sep = "|")
)

zero_width <- make_sts(matrix(integer(), nrow = 3, ncol = 0), c("A", "B"))
capture_condition("zero_width_sorti", sorti(zero_width))
capture_condition("zero_width_sortv", sortv(zero_width))

writeLines(condition_lines, file.path(out_dir, "sortv_r3_conditions.txt"), useBytes = TRUE)

cat("SORTV_R3_GOLDENS_DONE\n")
