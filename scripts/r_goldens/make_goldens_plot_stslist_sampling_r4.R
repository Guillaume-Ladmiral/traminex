# Golden generation for TraMineR 2.2-14 plot.stslist sampling mechanics.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

raw <- data.frame(
  t1 = c("B", "A", "A", "C", "A", "B"),
  t2 = c("B", "C", "A", "C", NA, "A"),
  t3 = c("A", "C", "B", "C", NA, "A"),
  t4 = c("A", "B", "B", "C", NA, NA),
  stringsAsFactors = FALSE
)
ids <- paste0("row-", seq_len(nrow(raw)))
seqdata <- suppressMessages(
  seqdef(
    raw,
    alphabet = c("A", "B", "C"),
    states = c("Alpha", "Beta", "Gamma"),
    right = "DEL",
    missing = NA,
    weights = c(1, 4, 0, 2, 3, 5),
    start = 7,
    xtstep = 2,
    tick.last = TRUE,
    cpal = c("#111111", "#555555", "#999999")
  )
)
rownames(seqdata) <- ids

# This deliberately follows the ordering and selection statements in
# SOURCES/TraMineR/R/plot.stslist.R.  The source row is retained through the
# stslist subset so the artifact can pin the pre-sampling sort contract.
ordered_rows <- function(x, sortv = NULL) {
  source_index <- seq_len(nrow(x))
  if (!is.null(sortv)) {
    if (length(sortv) == 1 && sortv %in% c("from.start", "from.end")) {
      end <- if (sortv == "from.end") max(seqlength(x)) else 1
      beg <- if (sortv == "from.end") 1 else max(seqlength(x))
      permutation <- do.call(order, unname(as.data.frame(x))[, end:beg])
    } else {
      if (is.factor(sortv)) sortv <- as.integer(sortv)
      permutation <- order(sortv)
    }
    x <- x[permutation, ]
    source_index <- source_index[permutation]
  }
  list(x = x, source_index = source_index)
}

selection_rows <- function(case_name, x, sortv = NULL, sampling = NULL,
                           sample.meth = "systematic", seed = NULL) {
  ordered <- ordered_rows(x, sortv)
  candidate_count <- nrow(ordered$x)
  if (is.null(sampling)) sampling <- candidate_count
  if (!is.null(seed)) set.seed(seed)
  if (candidate_count > sampling) {
    if (sample.meth == "systematic") {
      selected <- as.integer(seq(1, to = candidate_count, length.out = sampling))
    } else {
      selected <- sample.int(
        candidate_count,
        size = sampling,
        prob = attr(ordered$x, "weights")
      )
      selected <- selected[order(selected)]
    }
  } else {
    selected <- seq_len(candidate_count)
  }
  data.frame(
    case = case_name,
    display_index_zero = seq_along(selected) - 1L,
    candidate_index_zero = selected - 1L,
    source_index_zero = ordered$source_index[selected] - 1L,
    source_id = rownames(ordered$x)[selected],
    selected_weight = attr(ordered$x, "weights")[selected],
    stringsAsFactors = FALSE
  )
}

cases <- list(
  selection_rows("input_systematic_3", seqdata, sampling = 3),
  selection_rows("from_start_systematic_3", seqdata, sortv = "from.start", sampling = 3),
  selection_rows("from_end_systematic_4", seqdata, sortv = "from.end", sampling = 4),
  selection_rows(
    "numeric_ties_systematic_3", seqdata,
    sortv = c(2, 0, 2, 1, 0, 1), sampling = 3
  ),
  selection_rows(
    "factor_ties_systematic_4", seqdata,
    sortv = factor(c("z", "x", "z", "y", "x", "y"), levels = c("y", "z", "x")),
    sampling = 4
  ),
  selection_rows(
    "factor_missing_all", seqdata,
    sortv = factor(c("b", NA, "a", "b", NA, "a"), levels = c("a", "b")),
    sampling = nrow(seqdata)
  ),
  selection_rows("input_one", seqdata, sampling = 1),
  selection_rows("input_equal", seqdata, sampling = nrow(seqdata)),
  selection_rows("input_above", seqdata, sampling = nrow(seqdata) + 3L),
  selection_rows("input_no_sampling", seqdata, sampling = NULL),
  selection_rows("random_seed_1729", seqdata, sampling = 3, sample.meth = "random", seed = 1729)
)
write.csv(
  do.call(rbind, cases),
  file.path(out_dir, "plot_stslist_sampling_r4_selections.csv"),
  row.names = FALSE,
  quote = TRUE
)

metadata <- data.frame(
  section = c(
    rep("sequence", 8),
    rep("deliberate_divergence", 2),
    rep("source", 3)
  ),
  key = c(
    "nrow", "ncol", "alphabet", "labels", "start", "xtstep", "tick.last", "ids",
    "python_sampling_none", "python_random_stream",
    "sort_before_sampling", "systematic_expression", "random_expression"
  ),
  value = c(
    nrow(seqdata), ncol(seqdata),
    paste(alphabet(seqdata), collapse = "\036"),
    paste(attr(seqdata, "labels"), collapse = "\036"),
    attr(seqdata, "start"), attr(seqdata, "xtstep"), attr(seqdata, "tick.last"),
    paste(rownames(seqdata), collapse = "\036"),
    "all rows (Rd contract; R code defaults idxs to first ten when idxs is NULL)",
    "isolated numpy.random.Generator; no byte match to R",
    "TRUE",
    "as.integer(seq(1, length(idxs), length.out=sampling))",
    "sort(sample.int(length(idxs), size=sampling, prob=weights[idxs]))"
  ),
  stringsAsFactors = FALSE
)
write.csv(
  metadata,
  file.path(out_dir, "plot_stslist_sampling_r4_metadata.csv"),
  row.names = FALSE,
  quote = TRUE
)

condition_message <- function(expr) {
  tryCatch(
    {
      force(expr)
      "NO_ERROR"
    },
    error = function(condition) conditionMessage(condition)
  )
}
conditions <- data.frame(
  case = c("unknown_method", "zero_sampling", "negative_sampling"),
  message = c(
    condition_message(match.arg("bogus", c("systematic", "random"))),
    condition_message(selection_rows("invalid", seqdata, sampling = 0)),
    condition_message(selection_rows("invalid", seqdata, sampling = -1))
  ),
  stringsAsFactors = FALSE
)
write.csv(
  conditions,
  file.path(out_dir, "plot_stslist_sampling_r4_conditions.csv"),
  row.names = FALSE,
  quote = TRUE
)
