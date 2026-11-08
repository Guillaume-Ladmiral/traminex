# Golden generation for TraMineRextras 0.6.9 seqgen.missing.
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

trace_seqgen <- function(
  seqdata,
  seed,
  p.cases,
  p.left,
  p.gaps,
  p.right,
  mt.left,
  mt.gaps,
  mt.right
) {
  n <- nrow(seqdata)
  lgth <- max(seqlength(seqdata))
  nr.l <- attr(seqdata, mt.left)
  nr.g <- attr(seqdata, mt.gaps)
  nr.r <- attr(seqdata, mt.right)
  nm <- round(p.cases * n, 0)

  set.seed(seed)
  idm <- sort(sample(1:n, nm))
  rdu.l <- runif(n, min = 0, max = p.left)
  rdu.g <- runif(n, min = 0, max = p.gaps)
  rdu.r <- runif(n, min = 0, max = p.right)
  traced <- seqdata
  plans <- list()
  for (i in idm) {
    gap.count <- round(rdu.g[i] * lgth, 0)
    gaps <- sample(1:lgth, gap.count)
    traced[i, gaps] <- nr.g
    left.count <- round(rdu.l[i] * lgth, 0)
    if (left.count > 0) traced[i, 1:left.count] <- nr.l
    right.count <- round(rdu.r[i] * lgth, 0)
    if (right.count > 0) traced[i, (lgth - right.count + 1):lgth] <- nr.r
    plans[[as.character(i)]] <- list(
      gap.count = gap.count,
      gaps = gaps,
      left.count = left.count,
      right.count = right.count
    )
  }

  set.seed(seed)
  actual <- seqgen.missing(
    seqdata,
    p.cases = p.cases,
    p.left = p.left,
    p.gaps = p.gaps,
    p.right = p.right,
    mt.left = mt.left,
    mt.gaps = mt.gaps,
    mt.right = mt.right
  )
  stopifnot(identical(as.data.frame(actual), as.data.frame(traced)))
  list(result = actual, selected = idm, plans = plans, n = n, length = lgth, nm = nm)
}

records <- list()
record_index <- 1L
record <- function(case_name, record_type, index, value) {
  records[[record_index]] <<- data.frame(
    case = case_name,
    record_type = record_type,
    index = as.character(index),
    value = as.character(value),
    stringsAsFactors = FALSE
  )
  record_index <<- record_index + 1L
}

for (p.cases in c(0, 0.1, 0.3, 0.5, 1)) {
  record(
    "rounded_case_counts",
    "count",
    format(p.cases, trim = TRUE),
    round(p.cases * 5, 0)
  )
}

base <- make_sts(
  matrix(rep(c(1, 2, 1, 2, 1, 2), 6), byrow = TRUE, nrow = 6),
  c("A", "B")
)
overlap <- trace_seqgen(
  base,
  seed = 1407,
  p.cases = 2 / 3,
  p.left = 0.9,
  p.gaps = 0.9,
  p.right = 0.9,
  mt.left = "void",
  mt.gaps = "nr",
  mt.right = "nr"
)
record("seeded_overlap", "summary", "n", overlap$n)
record("seeded_overlap", "summary", "length", overlap$length)
record("seeded_overlap", "summary", "rounded_cases", overlap$nm)
record(
  "seeded_overlap",
  "summary",
  "selected_ids",
  paste(overlap$selected, collapse = ",")
)
record("seeded_overlap", "summary", "assignment_order", "gaps,left,right")
record("seeded_overlap", "summary", "overwrite_precedence", "right>left>gaps")
for (id in overlap$selected) {
  plan <- overlap$plans[[as.character(id)]]
  record("seeded_overlap", "left_count", id, plan$left.count)
  record("seeded_overlap", "gap_count", id, plan$gap.count)
  record(
    "seeded_overlap",
    "gap_positions",
    id,
    paste(sort(plan$gaps), collapse = ",")
  )
  record("seeded_overlap", "right_count", id, plan$right.count)
}
for (row in seq_len(nrow(overlap$result))) {
  for (column in seq_len(ncol(overlap$result))) {
    record(
      "seeded_overlap",
      "output",
      paste(row, column, sep = ","),
      as.character(overlap$result[row, column])
    )
  }
}

ragged <- make_sts(
  matrix(c(
    1, 2, 1, 2, 1,
    2, 1, 2, 5, 5,
    1, 2, 5, 5, 5,
    2, 1, 2, 1, 5
  ), byrow = TRUE, ncol = 5),
  c("A", "B")
)
ragged_trace <- trace_seqgen(
  ragged,
  seed = 1411,
  p.cases = 1,
  p.left = 0,
  p.gaps = 0,
  p.right = 0,
  mt.left = "nr",
  mt.gaps = "void",
  mt.right = "nr"
)
record("ragged_zero", "summary", "length", ragged_trace$length)
record("ragged_zero", "summary", "selected_ids", paste(ragged_trace$selected, collapse = ","))
record(
  "ragged_zero",
  "summary",
  "unchanged",
  identical(as.data.frame(ragged_trace$result), as.data.frame(ragged))
)

write.csv(
  do.call(rbind, records),
  file.path(out_dir, "seqgen_missing_r3_manifest.csv"),
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
          paste(case_name, "warning", conditionMessage(condition), sep = "|")
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

set.seed(12)
capture_condition(
  "invalid_missing_type",
  seqgen.missing(base, p.cases = 1, p.gaps = 1, mt.gaps = "bad")
)
set.seed(12)
capture_condition("p_cases_above_one", seqgen.missing(base, p.cases = 2))
set.seed(12)
capture_condition("p_left_below_zero", seqgen.missing(base, p.cases = 1, p.left = -0.1))
set.seed(12)
capture_condition("p_cases_nan", seqgen.missing(base, p.cases = NaN))

writeLines(
  condition_lines,
  file.path(out_dir, "seqgen_missing_r3_conditions.txt"),
  useBytes = TRUE
)

cat("SEQGEN_MISSING_R3_GOLDENS_DONE\n")
