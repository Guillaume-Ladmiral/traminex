# Public stochastic-stream evidence for the Phase 4 semantic audit.

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

rows <- list()
append_values <- function(case, field, values) {
  for (position in seq_along(values)) {
    rows[[length(rows) + 1L]] <<- data.frame(
      case = case,
      field = field,
      position = position,
      value = as.character(values[[position]]),
      stringsAsFactors = FALSE
    )
  }
}

set.seed(24680)
seqgen_before <- .Random.seed
seqgen_first <- TraMineR::seqgen(
  6L,
  7L,
  c("A", "B", "C"),
  p = c(0.55, 0.30, 0.15)
)
seqgen_after <- .Random.seed
set.seed(24680)
seqgen_second <- TraMineR::seqgen(
  6L,
  7L,
  c("A", "B", "C"),
  p = c(0.55, 0.30, 0.15)
)
stopifnot(identical(seqgen_first, seqgen_second))
append_values("seqgen", "draw", as.vector(t(seqgen_first)))
append_values("seqgen", "repeat_identical", identical(seqgen_first, seqgen_second))
append_values("seqgen", "global_rng_advanced", !identical(seqgen_before, seqgen_after))

points <- c(0.0, 0.7, 1.9, 3.6, 5.8, 8.5)
diss <- abs(outer(points, points, "-"))
set.seed(97531)
dissrep_before <- .Random.seed
dissrep_first <- suppressMessages(
  TraMineR::dissrep(
    diss,
    criterion = "random",
    nrep = 3L,
    pradius = 0
  )
)
dissrep_after <- .Random.seed
set.seed(97531)
dissrep_second <- suppressMessages(
  TraMineR::dissrep(
    diss,
    criterion = "random",
    nrep = 3L,
    pradius = 0
  )
)
stopifnot(
  identical(as.vector(dissrep_first), as.vector(dissrep_second)),
  identical(attr(dissrep_first, "Scores"), attr(dissrep_second, "Scores"))
)
append_values(
  "dissrep",
  "representative_zero_based",
  as.integer(as.vector(dissrep_first)) - 1L
)
append_values("dissrep", "score", attr(dissrep_first, "Scores"))
append_values(
  "dissrep",
  "repeat_identical",
  identical(as.vector(dissrep_first), as.vector(dissrep_second))
)
append_values(
  "dissrep",
  "global_rng_advanced",
  !identical(dissrep_before, dissrep_after)
)

output <- do.call(rbind, rows)
write.csv(
  output,
  "tests/r_goldens/stochastic_phase4.csv",
  row.names = FALSE,
  quote = TRUE,
  na = "NA"
)
