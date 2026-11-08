# Golden generation for TraMineRextras 0.6.9 group.p.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineRextras))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
stopifnot(as.character(packageVersion("TraMineRextras")) == "0.6.9")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cases <- list(
  character_unit = list(group = c("b", "a", "b", "c"), weights = NULL),
  numeric_integer = list(group = c(10, 2, 10, 3), weights = c(2, 1, 7, 2)),
  factor_unused_na = list(
    group = factor(c("b", "a", "b", NA), levels = c("z", "b", "a", "c")),
    weights = NULL
  ),
  fractional = list(
    group = factor(c("b", "a", "b", "c"), levels = c("c", "b", "a", "z")),
    weights = c(0.1, 0.2, 0.3, 0.4)
  ),
  third_rounding = list(group = c("x", "y", "z"), weights = c(1, 1, 1)),
  zero_total = list(
    group = factor(c("b", "a"), levels = c("b", "a", "z")),
    weights = c(0, 0)
  ),
  na_group = list(group = c("b", NA, "a"), weights = c(1, 2, 3)),
  na_weight = list(group = c("b", "a", "b"), weights = c(1, NA, 3)),
  mixed_case = list(group = c("a", "A"), weights = NULL),
  signed_zero = list(group = c("a", "b"), weights = c(1, -1)),
  infinite_weight = list(group = c("a", "b"), weights = c(Inf, 1)),
  locale_labels = list(
    group = c("_", "-", "!", "?", ".", "*", "%", "a", "A", "á", "ä", "b"),
    weights = NULL
  )
)

records <- list()
record_index <- 1L
for (case_name in names(cases)) {
  case <- cases[[case_name]]
  result <- group.p(case$group, weights = case$weights)
  for (position in seq_along(result)) {
    records[[record_index]] <- data.frame(
      case = case_name,
      record_type = "value",
      position = position,
      code = if (is.na(result[[position]])) "NA" else as.character(as.integer(result[[position]])),
      value = if (is.na(result[[position]])) "<NA>" else as.character(result[[position]]),
      level = "",
      stringsAsFactors = FALSE
    )
    record_index <- record_index + 1L
  }
  for (position in seq_along(levels(result))) {
    records[[record_index]] <- data.frame(
      case = case_name,
      record_type = "level",
      position = position,
      code = "",
      value = "",
      level = levels(result)[[position]],
      stringsAsFactors = FALSE
    )
    record_index <- record_index + 1L
  }
}
write.csv(
  do.call(rbind, records),
  file.path(out_dir, "group_p_r3_cases.csv"),
  row.names = FALSE,
  na = "<NA>",
  quote = TRUE
)

condition_lines <- character()
capture_error <- function(case_name, expression) {
  tryCatch(
    force(expression),
    error = function(condition) {
      condition_lines <<- c(
        condition_lines,
        paste(case_name, class(condition)[[1]], conditionMessage(condition), sep = "|")
      )
    }
  )
}
capture_error("weight_length", group.p(c("a", "b"), weights = 1))
writeLines(condition_lines, file.path(out_dir, "group_p_r3_errors.txt"), useBytes = TRUE)

cat("GROUP_P_R3_GOLDENS_DONE\n")
