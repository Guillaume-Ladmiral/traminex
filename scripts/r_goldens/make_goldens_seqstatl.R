# Public TraMineR 2.2-14 evidence for seqstatl.
#
# Unlike the historical one-value-per-line fixtures, this manifest records
# R's output storage type.  That distinction is part of the public contract:
# wholly numeric state sets return a double vector, not character labels.
#
# Run from the repository root:
#   Rscript scripts/r_goldens/make_goldens_seqstatl.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))

if (!identical(as.character(packageVersion("TraMineR")), "2.2.14")) {
  stop("seqstatl goldens require TraMineR 2.2-14")
}

source(file.path("scripts", "r_goldens", "make_goldens_corestats.R"))

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

csv_quote <- function(value) {
  paste0('"', gsub('"', '""', as.character(value), fixed = TRUE), '"')
}

value_text <- function(value) {
  if (typeof(value) == "double") {
    return(trimws(formatC(value, format = "g", digits = 15)))
  }
  as.character(value)
}

rows <- list()
append_case <- function(case, call, value) {
  for (position in seq_along(value)) {
    rows[[length(rows) + 1L]] <<- paste(
      csv_quote(case),
      csv_quote(call),
      csv_quote(typeof(value)),
      position,
      csv_quote(value_text(value[[position]])),
      sep = ","
    )
  }
}

fx <- corestats_fixtures()
for (tag in c("A", "B", "C")) {
  append_case(
    paste0("sts_", tag),
    paste0("TraMineR::seqstatl(fx$", tag, ")"),
    TraMineR::seqstatl(fx[[tag]])
  )
}

raw_extended <- data.frame(
  t1 = c("é", "A", NA, "-"),
  t2 = c("Z", "ä", "10", "!"),
  t3 = c("2", "A", "2", "é"),
  check.names = FALSE
)
append_case(
  "raw_extended",
  "TraMineR::seqstatl(raw_extended)",
  TraMineR::seqstatl(raw_extended)
)
append_case(
  "named_var",
  "TraMineR::seqstatl(raw_extended,var=c('t3','t1'))",
  TraMineR::seqstatl(raw_extended, var = c("t3", "t1"))
)
append_case(
  "integer_var",
  "TraMineR::seqstatl(raw_extended,var=c(3,1))",
  TraMineR::seqstatl(raw_extended, var = c(3, 1))
)

append_case(
  "compressed_dash",
  "TraMineR::seqstatl(c('A-B-C','B-D'))",
  TraMineR::seqstatl(c("A-B-C", "B-D"))
)
append_case(
  "compressed_colon",
  "TraMineR::seqstatl(c('A:é:!','Z:ä:A'))",
  TraMineR::seqstatl(c("A:é:!", "Z:ä:A"))
)
append_case(
  "sps",
  "TraMineR::seqstatl(c('(A,2)-(B,1)','(B,1)-(C,2)'),format='SPS')",
  TraMineR::seqstatl(
    c("(A,2)-(B,1)", "(B,1)-(C,2)"),
    format = "SPS"
  )
)

format_data <- rbind(c("B", "A"), c("C", "B"))
for (format in c("STS", "SPELL", "DSS")) {
  append_case(
    paste0("format_", format),
    paste0("TraMineR::seqstatl(format_data,format='", format, "')"),
    TraMineR::seqstatl(format_data, format = format)
  )
}

append_case(
  "numeric_integer_text",
  "TraMineR::seqstatl(rbind(c('3','1','10'),c('2','10','1')))",
  TraMineR::seqstatl(rbind(c("3", "1", "10"), c("2", "10", "1")))
)
append_case(
  "numeric_actual",
  "TraMineR::seqstatl(matrix(c(3,1,10,2,10,1),nrow=2,byrow=TRUE))",
  TraMineR::seqstatl(matrix(c(3, 1, 10, 2, 10, 1), nrow = 2, byrow = TRUE))
)
append_case(
  "numeric_decimal_text",
  "TraMineR::seqstatl(rbind(c('2.5','1.2'),c('10','2.5')))",
  TraMineR::seqstatl(rbind(c("2.5", "1.2"), c("10", "2.5")))
)
append_case(
  "mixed_numeric_text",
  "TraMineR::seqstatl(rbind(c('10','A'),c('2','!')))",
  TraMineR::seqstatl(rbind(c("10", "A"), c("2", "!")))
)
append_case(
  "raw_missing",
  "TraMineR::seqstatl(matrix(c(NA,'B','A',NA),nrow=2))",
  TraMineR::seqstatl(matrix(c(NA, "B", "A", NA), nrow = 2))
)
append_case(
  "raw_void_literal",
  "TraMineR::seqstatl(matrix(c('%','A',NA,'%'),nrow=2))",
  TraMineR::seqstatl(matrix(c("%", "A", NA, "%"), nrow = 2))
)

manifest <- file.path(out_dir, "seqstatl_public_cases.csv")
con <- file(manifest, open = "wt", encoding = "UTF-8")
writeLines(
  "case,call,r_type,position,value",
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

errors <- list()
append_error <- function(case, call, expression) {
  condition <- tryCatch(
    {
      force(expression)
      NULL
    },
    error = function(error) error
  )
  if (is.null(condition)) {
    stop(paste("expected an error for", case))
  }
  errors[[length(errors) + 1L]] <<- paste(
    case,
    call,
    class(condition)[1L],
    conditionMessage(condition),
    sep = "|"
  )
}

append_error(
  "invalid_lowercase_format",
  "TraMineR::seqstatl(format_data,format='sts')",
  TraMineR::seqstatl(format_data, format = "sts")
)

writeLines(
  c("case|call|condition|message", unlist(errors, use.names = FALSE)),
  file.path(out_dir, "seqstatl_public_errors.txt"),
  useBytes = TRUE
)

cat("SEQSTATL_GOLDENS_DONE\n")
