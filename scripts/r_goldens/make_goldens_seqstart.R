# Golden generation for TraMineRextras::seqstart placement, correction,
# clipping, scalar/vector handling, missing fill, and error branches.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
suppressMessages(library(TraMineRextras))
stopifnot(
  as.character(packageVersion("TraMineR")) == "2.2.14",
  as.character(packageVersion("TraMineRextras")) == "0.6.9"
)

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

seqdata <- matrix(
  c("A", "B", "C", "D", "E", "F"),
  byrow = TRUE,
  nrow = 2
)

encode_vector <- function(value) {
  paste(as.character(value), collapse = ",")
}

encode_bound <- function(value) {
  if (is.null(value)) "NULL" else as.character(value)
}

encode_cells <- function(value) {
  rows <- character(nrow(value))
  for (i in seq_len(nrow(value))) {
    cells <- ifelse(is.na(value[i, ]), "<NA>", as.character(value[i, ]))
    rows[i] <- paste(cells, collapse = ",")
  }
  paste(rows, collapse = ";")
}

success_rows <- character()
emit_success <- function(
  case,
  data.start,
  new.start,
  tmin = NULL,
  tmax = NULL,
  missing = "*",
  default.missing = FALSE
) {
  args <- list(
    seqdata = seqdata,
    data.start = data.start,
    new.start = new.start,
    tmin = tmin,
    tmax = tmax
  )
  if (!default.missing) {
    args$missing <- missing
  }
  result <- do.call(TraMineRextras::seqstart, args)
  success_rows <<- c(
    success_rows,
    paste(
      case,
      encode_vector(data.start),
      encode_vector(new.start),
      encode_bound(tmin),
      encode_bound(tmax),
      if (default.missing) "<DEFAULT_NA>" else as.character(missing),
      paste(dim(result), collapse = "x"),
      encode_cells(result),
      sep = "\t"
    )
  )
}

emit_success("negative_correction", c(1998, 2000), 2000)
emit_success("scalar_scalar", 2000, 2000)
emit_success("scalar_vector_default_missing", 2000, c(2000, 1999), default.missing = TRUE)
emit_success("clip_right", c(1998, 2000), c(2000, 2000), tmin = 1, tmax = 4)
emit_success("clip_left", c(1998, 2000), c(2000, 2000), tmin = -1, tmax = 2)
emit_success("vector_vector", c(2000, 2001), c(2000, 1999))

writeLines(
  c(
    "# R 4.6.2; TraMineR 2.2.14; TraMineRextras 0.6.9",
    "# case\tdata.start\tnew.start\ttmin\ttmax\tmissing\tshape\tcells",
    success_rows
  ),
  file.path(out_dir, "seqstart_success_cases.txt"),
  useBytes = TRUE
)

error_rows <- character()
emit_error <- function(case, data.start, new.start, tmin = NULL, tmax = NULL) {
  message_text <- tryCatch(
    {
      TraMineRextras::seqstart(
        seqdata,
        data.start = data.start,
        new.start = new.start,
        tmin = tmin,
        tmax = tmax,
        missing = "*"
      )
      "<NO_ERROR>"
    },
    error = function(condition) conditionMessage(condition)
  )
  error_rows <<- c(
    error_rows,
    paste(
      case,
      encode_vector(data.start),
      encode_vector(new.start),
      encode_bound(tmin),
      encode_bound(tmax),
      gsub("[\t\r\n]+", " ", message_text),
      sep = "\t"
    )
  )
}

emit_error("length_mismatch", c(2000, 2001, 2002), 2000)
emit_error("negative_tmax", 2000, 2000, tmax = -1)

writeLines(
  c(
    "# R 4.6.2; TraMineR 2.2.14; TraMineRextras 0.6.9",
    "# case\tdata.start\tnew.start\ttmin\ttmax\terror",
    error_rows
  ),
  file.path(out_dir, "seqstart_error_cases.txt"),
  useBytes = TRUE
)

cat("SEQSTART_GOLDENS_DONE\n")
