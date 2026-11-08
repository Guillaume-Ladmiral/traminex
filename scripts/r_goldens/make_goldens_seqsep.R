# Golden generation for TraMineR 2.2-14 seqsep.
#
# Successful results are emitted one vector element per TSV row so empty strings
# remain unambiguous. Error cases record R's condition class and exact message.
# There are no numeric results, so the project's formatC numeric rules do not apply.
#
# Run from the repository root:
#   Rscript scripts/r_goldens/make_goldens_seqsep.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressPackageStartupMessages(library(TraMineR))

if (as.character(packageVersion("TraMineR")) != "2.2.14") {
  stop("seqsep goldens require TraMineR 2.2-14")
}

out_dir <- file.path("tests", "r_goldens")

success_cases <- list(
  list(name = "sl1_default", input = c("ABCD", "XYZ"), sl = 1, sep = "-"),
  list(name = "sl1_pipe", input = c("AA", "B"), sl = 1, sep = "|"),
  list(name = "sl2_default", input = c("AABB", "CCDD", "EE"), sl = 2, sep = "-"),
  list(name = "sl2_double_colon", input = c("AABBCC", "DDEE"), sl = 2, sep = "::"),
  list(name = "sl3_slash", input = c("AAABBB", "CCCDDDFFF", "XYZ"), sl = 3, sep = "/"),
  list(name = "empty_after_value", input = c("AB", ""), sl = 1, sep = "-")
)

error_cases <- list(
  list(
    name = "empty_vector",
    input = character(0),
    sl = 1,
    sep = "-",
    failing_index = 0
  ),
  list(name = "empty_only", input = c(""), sl = 1, sep = "-", failing_index = 1),
  list(name = "empty_first", input = c("", "AB"), sl = 1, sep = "-", failing_index = 1),
  list(name = "nondivisible_first", input = c("A", "BB"), sl = 2, sep = "-", failing_index = 1),
  list(
    name = "nondivisible_second",
    input = c("AABB", "CCC", "DDEE"),
    sl = 2,
    sep = "-",
    failing_index = 2
  ),
  list(
    name = "nondivisible_third",
    input = c("AABB", "CC", "DDD"),
    sl = 2,
    sep = "|",
    failing_index = 3
  )
)

case_rows <- do.call(
  rbind,
  lapply(success_cases, function(case) {
    result <- TraMineR::seqsep(case$input, sl = case$sl, sep = case$sep)
    data.frame(
      case = case$name,
      sl = as.character(case$sl),
      sep = case$sep,
      position = as.character(seq_along(case$input)),
      input = case$input,
      output = result,
      stringsAsFactors = FALSE
    )
  })
)

case_path <- file.path(out_dir, "seqsep_cases.txt")
case_con <- file(case_path, open = "wt")
writeLines(
  c(
    "# R oracle: TraMineR 2.2-14 seqsep",
    "# TSV columns: case, sl, sep, position, input, output"
  ),
  case_con
)
write.table(
  case_rows,
  file = case_con,
  sep = "\t",
  quote = TRUE,
  row.names = FALSE,
  col.names = FALSE,
  na = "NA"
)
close(case_con)

error_rows <- do.call(
  rbind,
  lapply(error_cases, function(case) {
    condition <- tryCatch(
      {
        TraMineR::seqsep(case$input, sl = case$sl, sep = case$sep)
        NULL
      },
      error = identity
    )
    if (is.null(condition)) {
      stop(sprintf("expected seqsep case '%s' to fail", case$name))
    }
    encoded_input <- if (length(case$input) == 0) {
      "[]"
    } else {
      paste(
        sprintf("%d=%s", seq_along(case$input), encodeString(case$input, quote = '"')),
        collapse = " | "
      )
    }
    data.frame(
      case = case$name,
      sl = as.character(case$sl),
      sep = case$sep,
      failing_index = as.character(case$failing_index),
      input = encoded_input,
      condition_class = class(condition)[1],
      message = conditionMessage(condition),
      stringsAsFactors = FALSE
    )
  })
)

error_path <- file.path(out_dir, "seqsep_errors.txt")
error_con <- file(error_path, open = "wt")
writeLines(
  c(
    "# R oracle: TraMineR 2.2-14 seqsep errors",
    "# TSV columns: case, sl, sep, failing_index, input, condition_class, message"
  ),
  error_con
)
write.table(
  error_rows,
  file = error_con,
  sep = "\t",
  quote = TRUE,
  row.names = FALSE,
  col.names = FALSE,
  na = "NA"
)
close(error_con)

cat("SEQSEP_GOLDENS_DONE\n")
