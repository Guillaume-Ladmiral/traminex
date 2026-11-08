# Public TraMineR 2.2-14 seqsubm wrapper evidence.
#
# This generator deliberately calls TraMineR::seqsubm rather than deriving the
# expected values from Python or calling seqcost on the Python side.
#
# Run from the repository root:
#   /usr/bin/Rscript scripts/r_goldens/make_goldens_seqsubm.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))

if (!identical(as.character(packageVersion("TraMineR")), "2.2.14")) {
  stop("seqsubm goldens require TraMineR 2.2-14")
}

fixture_env <- new.env(parent = globalenv())
sys.source(
  file.path("scripts", "r_goldens", "make_goldens_corestats.R"),
  envir = fixture_env
)
fx <- fixture_env$corestats_fixtures()

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

format_value <- function(value) {
  if (is.nan(value)) return("NaN")
  if (is.na(value)) return("NA")
  if (is.infinite(value)) return(if (value > 0) "Inf" else "-Inf")
  trimws(formatC(value, format = "g", digits = 15))
}

csv_quote <- function(value) {
  paste0('"', gsub('"', '""', as.character(value), fixed = TRUE), '"')
}

append_result <- function(rows, case, options, result) {
  stopifnot(!is.null(result))
  dims <- dim(result)
  labels <- dimnames(result)
  stopifnot(length(dims) %in% c(2L, 3L))
  for (t in seq_len(if (length(dims) == 3L) dims[[3]] else 1L)) {
    time_label <- if (length(dims) == 3L) labels[[3]][[t]] else ""
    for (i in seq_len(dims[[1]])) {
      for (j in seq_len(dims[[2]])) {
        value <- if (length(dims) == 3L) result[i, j, t] else result[i, j]
        rows[[length(rows) + 1L]] <- c(
          csv_quote(case),
          csv_quote(options),
          csv_quote(labels[[1]][[i]]),
          csv_quote(labels[[2]][[j]]),
          csv_quote(time_label),
          format_value(value)
        )
      }
    }
  }
  rows
}

rows <- list()
rows <- append_result(
  rows,
  "trate_unweighted",
  "TraMineR::seqsubm(fx$A,method=TRATE,weighted=FALSE)",
  suppressMessages(TraMineR::seqsubm(fx$A, method = "TRATE", weighted = FALSE))
)
rows <- append_result(
  rows,
  "trate_weighted",
  "TraMineR::seqsubm(fx$Aw,method=TRATE,weighted=TRUE)",
  suppressMessages(TraMineR::seqsubm(fx$Aw, method = "TRATE", weighted = TRUE))
)
rows <- append_result(
  rows,
  "trate_time_missing",
  paste0(
    "TraMineR::seqsubm(fx$Bw,method=TRATE,time.varying=TRUE,",
    "with.missing=TRUE,weighted=TRUE,transition=both,lag=1,",
    "miss.cost=7|8|9,miss.cost.fixed=TRUE)"
  ),
  suppressMessages(
    TraMineR::seqsubm(
      fx$Bw,
      method = "TRATE",
      time.varying = TRUE,
      with.missing = TRUE,
      weighted = TRUE,
      transition = "both",
      lag = 1,
      miss.cost = c(7, 8, 9),
      miss.cost.fixed = TRUE
    )
  )
)

features <- data.frame(
  constant = c(5, 5, 5, 5),
  varying = c(0, 10, 0, 10)
)
rows <- append_result(
  rows,
  "features",
  paste0(
    "TraMineR::seqsubm(fx$A,method=FEATURES,",
    "state.features=constant|varying,feature.weights=2|1,",
    "feature.type=numeric:constant|varying)"
  ),
  suppressMessages(
    TraMineR::seqsubm(
      fx$A,
      method = "FEATURES",
      state.features = features,
      feature.weights = c(2, 1),
      feature.type = list(numeric = c("constant", "varying"))
    )
  )
)

path <- file.path(out_dir, "seqsubm_public_cases.csv")
con <- file(path, open = "wt", encoding = "UTF-8")
writeLines(
  "case,options,state_i,state_j,time_label,value",
  con = con,
  sep = "\n",
  useBytes = TRUE
)
writeLines(
  vapply(rows, paste, collapse = ",", FUN.VALUE = character(1)),
  con = con,
  sep = "\n",
  useBytes = TRUE
)
close(con)

prox <- suppressMessages(
  TraMineR::seqsubm(fx$A, method = "TRATE", proximities = TRUE)
)
if (!is.null(prox)) {
  stop("TraMineR::seqsubm(..., proximities=TRUE) must return NULL")
}
writeLines(
  c(
    "case|call|result",
    paste(
      "proximities_null",
      "TraMineR::seqsubm(fx$A,method=TRATE,proximities=TRUE)",
      "NULL",
      sep = "|"
    )
  ),
  con = file.path(out_dir, "seqsubm_public_outcomes.txt"),
  sep = "\n",
  useBytes = TRUE
)

cat("SEQSUBM_GOLDENS_DONE\n")
