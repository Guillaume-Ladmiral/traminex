# Golden generation: the full TraMineR 2.2-14 seqcost result contract.
#
# Values are written in an explicit logical order rather than by flattening R
# matrices/arrays (which would silently use column-major order).  Every row is
# self-describing so the Python tests can verify values, axis labels, and the
# scalar/vector/matrix indel variants from one deterministic artifact.
#
# Run from the repository root:
#   /usr/bin/Rscript scripts/r_goldens/make_goldens_seqcost.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))

if (!identical(as.character(packageVersion("TraMineR")), "2.2.14")) {
  stop("seqcost goldens require TraMineR 2.2-14")
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
  # formatC(..., digits=15) rounds DBL_MAX to 1.79769313486232e+308,
  # which parses as Inf in IEEE-754 binary64.  Keep the finite sentinel
  # distinguishable from genuine outer-sum overflow in the same artifacts.
  if (identical(value, .Machine$double.xmax)) return("DBL_MAX")
  if (is.infinite(value)) return(if (value > 0) "Inf" else "-Inf")
  trimws(formatC(value, format = "g", digits = 15))
}

csv_quote <- function(value) {
  paste0('"', gsub('"', '""', as.character(value), fixed = TRUE), '"')
}

golden_row <- function(case, options, component, state_i = "", state_j = "",
                       time_label = "", value) {
  c(
    csv_quote(case),
    csv_quote(options),
    csv_quote(component),
    csv_quote(state_i),
    csv_quote(state_j),
    csv_quote(time_label),
    format_value(value)
  )
}

write_golden_rows <- function(name, rows) {
  path <- file.path(out_dir, name)
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(
    "case,options,component,state_i,state_j,time_label,value",
    con = con,
    sep = "\n",
    useBytes = TRUE
  )
  if (length(rows)) {
    lines <- vapply(rows, paste, collapse = ",", FUN.VALUE = character(1))
    writeLines(lines, con = con, sep = "\n", useBytes = TRUE)
  }
}

append_seqcost_result <- function(rows, case, options, result) {
  component <- if ("prox" %in% names(result)) "prox" else "sm"
  stopifnot(identical(names(result), c("indel", component)))

  values <- result[[component]]
  value_dims <- dim(values)
  value_names <- dimnames(values)
  state_labels <- value_names[[1]]
  time_labels <- if (length(value_dims) == 3L) value_names[[3]] else character()

  indel <- result$indel
  indel_dims <- dim(indel)
  if (is.null(indel_dims)) {
    if (length(indel) == 1L) {
      rows[[length(rows) + 1L]] <- golden_row(
        case, options, "indel", value = indel[[1]]
      )
    } else {
      stopifnot(length(indel) == length(state_labels))
      for (i in seq_along(state_labels)) {
        rows[[length(rows) + 1L]] <- golden_row(
          case, options, "indel", state_i = state_labels[[i]], value = indel[[i]]
        )
      }
    }
  } else {
    stopifnot(length(indel_dims) == 2L)
    indel_names <- dimnames(indel)
    for (t in seq_len(indel_dims[[2]])) {
      for (i in seq_len(indel_dims[[1]])) {
        rows[[length(rows) + 1L]] <- golden_row(
          case,
          options,
          "indel",
          state_i = indel_names[[1]][[i]],
          time_label = indel_names[[2]][[t]],
          value = indel[i, t]
        )
      }
    }
  }

  if (length(value_dims) == 2L) {
    # Explicit row-major order: source row, then destination column.
    for (i in seq_len(value_dims[[1]])) {
      for (j in seq_len(value_dims[[2]])) {
        rows[[length(rows) + 1L]] <- golden_row(
          case,
          options,
          component,
          state_i = value_names[[1]][[i]],
          state_j = value_names[[2]][[j]],
          value = values[i, j]
        )
      }
    }
  } else {
    stopifnot(length(value_dims) == 3L)
    # Explicit logical order: time slice, source row, destination column.
    for (t in seq_len(value_dims[[3]])) {
      for (i in seq_len(value_dims[[1]])) {
        for (j in seq_len(value_dims[[2]])) {
          rows[[length(rows) + 1L]] <- golden_row(
            case,
            options,
            component,
            state_i = value_names[[1]][[i]],
            state_j = value_names[[2]][[j]],
            time_label = time_labels[[t]],
            value = values[i, j, t]
          )
        }
      }
    }
  }
  rows
}

add_case <- function(rows, case, options, result) {
  append_seqcost_result(rows, case, options, suppressMessages(result))
}

# Additional fixtures for edge contracts not present in the shared A/B/C set.
absent_raw <- matrix("A", nrow = 3L, ncol = 3L)
absent <- suppressMessages(seqdef(absent_raw, alphabet = c("A", "B", "C")))

one_state <- suppressMessages(
  seqdef(matrix("A", nrow = 2L, ncol = 3L), alphabet = "A")
)

short_raw <- rbind(c("A", "B"), c("B", "A"), c("A", "A"))
short <- suppressMessages(seqdef(short_raw, alphabet = c("A", "B")))

custom_time <- fx$A
colnames(custom_time) <- paste0("age", 10:15)

features_mixed <- data.frame(
  interval = c(0, 5, 10, 10),
  category = factor(c("x", "y", "x", "z")),
  flag = c(0, 1, 0, 1),
  ordered_score = c(10, 20, 40, 30),
  ratio = c(1, 10, 100, 1000)
)
features_mixed_type <- list(
  asymm = "flag",
  ordratio = "ordered_score",
  logratio = "ratio"
)
features_constant <- data.frame(
  constant = rep(5, 4),
  varying = c(0, 10, 0, 10)
)
features_missing <- data.frame(score = c(0, 0.5, 1, 0.25))

# -------------------------------------------------------------------------
# Stationary results
# -------------------------------------------------------------------------

stationary <- list()

stationary <- add_case(
  stationary,
  "constant_default",
  "method=CONSTANT",
  seqcost(fx$A, method = "CONSTANT")
)
stationary <- add_case(
  stationary,
  "constant_custom",
  "method=CONSTANT;cval=3.5",
  seqcost(fx$A, method = "CONSTANT", cval = 3.5)
)
stationary <- add_case(
  stationary,
  "constant_no_observed_missing",
  "method=CONSTANT;with_missing=TRUE;observed_missing=FALSE",
  seqcost(fx$A, method = "CONSTANT", with.missing = TRUE)
)

stationary <- add_case(
  stationary,
  "trate_unweighted_lag1",
  "method=TRATE;weighted=FALSE;lag=1",
  seqcost(fx$A, method = "TRATE", weighted = FALSE, lag = 1)
)
stationary <- add_case(
  stationary,
  "trate_weighted_lag1",
  "method=TRATE;weighted=TRUE;lag=1",
  seqcost(fx$Aw, method = "TRATE", weighted = TRUE, lag = 1)
)
stationary <- add_case(
  stationary,
  "trate_unweighted_lag2",
  "method=TRATE;weighted=FALSE;lag=2",
  seqcost(fx$A, method = "TRATE", weighted = FALSE, lag = 2)
)
stationary <- add_case(
  stationary,
  "trate_unweighted_lag_neg1",
  "method=TRATE;weighted=FALSE;lag=-1",
  seqcost(fx$A, method = "TRATE", weighted = FALSE, lag = -1)
)
stationary <- add_case(
  stationary,
  "trate_missing_fixed",
  "method=TRATE;weighted=TRUE;with_missing=TRUE;cval=3;miss_cost=7|8|9;miss_cost_fixed=TRUE",
  seqcost(
    fx$Bw,
    method = "TRATE",
    weighted = TRUE,
    with.missing = TRUE,
    cval = 3,
    miss.cost = c(7, 8, 9),
    miss.cost.fixed = TRUE
  )
)
stationary <- add_case(
  stationary,
  "trate_missing_derived",
  "method=TRATE;weighted=FALSE;with_missing=TRUE;cval=3;miss_cost_fixed=FALSE",
  seqcost(
    fx$B,
    method = "TRATE",
    weighted = FALSE,
    with.missing = TRUE,
    cval = 3,
    miss.cost.fixed = FALSE
  )
)

stationary <- add_case(
  stationary,
  "future_unweighted_lag2",
  "method=FUTURE;weighted=FALSE;lag=2",
  seqcost(fx$A, method = "FUTURE", weighted = FALSE, lag = 2)
)
stationary <- add_case(
  stationary,
  "future_weighted_lag1",
  "method=FUTURE;weighted=TRUE;lag=1",
  seqcost(fx$Aw, method = "FUTURE", weighted = TRUE, lag = 1)
)
stationary <- add_case(
  stationary,
  "future_missing_fixed",
  "method=FUTURE;weighted=FALSE;with_missing=TRUE;miss_cost=7;miss_cost_fixed=TRUE",
  seqcost(
    fx$B,
    method = "FUTURE",
    weighted = FALSE,
    with.missing = TRUE,
    miss.cost = 7,
    miss.cost.fixed = TRUE
  )
)
stationary <- add_case(
  stationary,
  "future_missing_derived",
  "method=FUTURE;weighted=FALSE;with_missing=TRUE;miss_cost_fixed=FALSE",
  seqcost(
    fx$B,
    method = "FUTURE",
    weighted = FALSE,
    with.missing = TRUE,
    miss.cost.fixed = FALSE
  )
)

stationary <- add_case(
  stationary,
  "features_mixed",
  "method=FEATURES;feature_weights=2|1|3|1|0.5;feature_type=asymm:flag|ordratio:ordered_score|logratio:ratio",
  seqcost(
    fx$A,
    method = "FEATURES",
    state.features = features_mixed,
    feature.weights = c(2, 1, 3, 1, 0.5),
    feature.type = features_mixed_type
  )
)
stationary <- add_case(
  stationary,
  "features_constant_numeric",
  "method=FEATURES;constant_numeric=TRUE",
  seqcost(fx$A, method = "FEATURES", state.features = features_constant)
)
stationary <- add_case(
  stationary,
  "features_missing_fixed",
  "method=FEATURES;with_missing=TRUE;miss_cost=2;miss_cost_fixed=TRUE",
  seqcost(
    fx$B,
    method = "FEATURES",
    with.missing = TRUE,
    state.features = features_missing,
    miss.cost = 2,
    miss.cost.fixed = TRUE
  )
)
stationary <- add_case(
  stationary,
  "features_missing_derived",
  "method=FEATURES;with_missing=TRUE;miss_cost_fixed=FALSE",
  seqcost(
    fx$B,
    method = "FEATURES",
    with.missing = TRUE,
    state.features = features_missing,
    miss.cost.fixed = FALSE
  )
)

stationary <- add_case(
  stationary,
  "indels_unweighted",
  "method=INDELS;weighted=FALSE",
  seqcost(fx$A, method = "INDELS", weighted = FALSE)
)
stationary <- add_case(
  stationary,
  "indels_weighted",
  "method=INDELS;weighted=TRUE",
  seqcost(fx$Aw, method = "INDELS", weighted = TRUE)
)
stationary <- add_case(
  stationary,
  "indels_missing_unweighted",
  "method=INDELS;weighted=FALSE;with_missing=TRUE",
  seqcost(fx$B, method = "INDELS", weighted = FALSE, with.missing = TRUE)
)
stationary <- add_case(
  stationary,
  "indels_missing_weighted",
  "method=INDELS;weighted=TRUE;with_missing=TRUE",
  seqcost(fx$Bw, method = "INDELS", weighted = TRUE, with.missing = TRUE)
)
stationary <- add_case(
  stationary,
  "indels_ragged_missing_weighted",
  "method=INDELS;weighted=TRUE;with_missing=TRUE;void_denominator=TRUE",
  seqcost(fx$Cw, method = "INDELS", weighted = TRUE, with.missing = TRUE)
)
stationary <- add_case(
  stationary,
  "indels_missing_fixed",
  "method=INDELS;weighted=FALSE;with_missing=TRUE;miss_cost=7|8|9;miss_cost_fixed=TRUE",
  seqcost(
    fx$B,
    method = "INDELS",
    weighted = FALSE,
    with.missing = TRUE,
    miss.cost = c(7, 8, 9),
    miss.cost.fixed = TRUE
  )
)
# TraMineR 2.2-14 calls checkcost() from inside seqcost(). For every INDELS
# case below (and for time-varying INDELS on the ordinary A/Aw fixtures too)
# at least one state is absent from some slice, so their indel is .Machine$double.xmax and the
# substitution matrix overflows to Inf -- which is exactly what these fixtures
# pin. checkcost() then computes tol <- 1e-07 * max(sma) = Inf and evaluates
# any((sm - t(sm)) > tol), i.e. any(NaN > Inf) = NA, and aborts with "missing
# value where TRUE/FALSE needed". That is an upstream regression (2.2-13 had no
# such call); check.costs=FALSE is 2.2-14's own switch for it and changes only
# the diagnostics, never the returned indel/sm this fixture records.
stationary <- add_case(
  stationary,
  "indels_absent_states",
  "method=INDELS;absent_states=B|C;finite_indel=double_max;overflow_sm=Inf",
  seqcost(absent, method = "INDELS", check.costs = FALSE)
)

stationary <- add_case(
  stationary,
  "indelslog_unweighted",
  "method=INDELSLOG;weighted=FALSE",
  seqcost(fx$A, method = "INDELSLOG", weighted = FALSE)
)
stationary <- add_case(
  stationary,
  "indelslog_weighted",
  "method=INDELSLOG;weighted=TRUE",
  seqcost(fx$Aw, method = "INDELSLOG", weighted = TRUE)
)
stationary <- add_case(
  stationary,
  "indelslog_missing_weighted",
  "method=INDELSLOG;weighted=TRUE;with_missing=TRUE",
  seqcost(fx$Bw, method = "INDELSLOG", weighted = TRUE, with.missing = TRUE)
)
stationary <- add_case(
  stationary,
  "indelslog_ragged_missing_weighted",
  "method=INDELSLOG;weighted=TRUE;with_missing=TRUE;void_denominator=TRUE",
  seqcost(fx$Cw, method = "INDELSLOG", weighted = TRUE, with.missing = TRUE)
)
stationary <- add_case(
  stationary,
  "indelslog_missing_fixed",
  "method=INDELSLOG;weighted=FALSE;with_missing=TRUE;miss_cost=7|8|9;miss_cost_fixed=TRUE",
  seqcost(
    fx$B,
    method = "INDELSLOG",
    weighted = FALSE,
    with.missing = TRUE,
    miss.cost = c(7, 8, 9),
    miss.cost.fixed = TRUE
  )
)
stationary <- add_case(
  stationary,
  "indelslog_absent_states",
  "method=INDELSLOG;absent_states=B|C",
  seqcost(absent, method = "INDELSLOG")
)

write_golden_rows("seqcost_stationary_cases.csv", stationary)

# -------------------------------------------------------------------------
# Time-varying results
# -------------------------------------------------------------------------

timevarying <- list()

timevarying <- add_case(
  timevarying,
  "constant_time_default",
  "method=CONSTANT;time_varying=TRUE",
  seqcost(fx$A, method = "CONSTANT", time.varying = TRUE)
)
timevarying <- add_case(
  timevarying,
  "constant_time_missing_vector",
  "method=CONSTANT;time_varying=TRUE;with_missing=TRUE;cval=3;miss_cost=7|8|9",
  seqcost(
    fx$B,
    method = "CONSTANT",
    time.varying = TRUE,
    with.missing = TRUE,
    cval = 3,
    miss.cost = c(7, 8, 9)
  )
)

timevarying <- add_case(
  timevarying,
  "trate_time_previous",
  "method=TRATE;time_varying=TRUE;transition=previous;weighted=FALSE;lag=1",
  seqcost(
    fx$A,
    method = "TRATE",
    time.varying = TRUE,
    transition = "previous",
    weighted = FALSE,
    lag = 1
  )
)
timevarying <- add_case(
  timevarying,
  "trate_time_next",
  "method=TRATE;time_varying=TRUE;transition=next;weighted=FALSE;lag=1",
  seqcost(
    fx$A,
    method = "TRATE",
    time.varying = TRUE,
    transition = "next",
    weighted = FALSE,
    lag = 1
  )
)
timevarying <- add_case(
  timevarying,
  "trate_time_both",
  "method=TRATE;time_varying=TRUE;transition=both;weighted=FALSE;lag=1;default_cval=4",
  seqcost(
    fx$A,
    method = "TRATE",
    time.varying = TRUE,
    transition = "both",
    weighted = FALSE,
    lag = 1
  )
)
timevarying <- add_case(
  timevarying,
  "trate_time_both_weighted",
  "method=TRATE;time_varying=TRUE;transition=both;weighted=TRUE;lag=1",
  seqcost(
    fx$Aw,
    method = "TRATE",
    time.varying = TRUE,
    transition = "both",
    weighted = TRUE,
    lag = 1
  )
)
timevarying <- add_case(
  timevarying,
  "trate_time_missing_fixed",
  "method=TRATE;time_varying=TRUE;transition=both;weighted=FALSE;with_missing=TRUE;miss_cost=7|8|9;miss_cost_fixed=TRUE",
  seqcost(
    fx$B,
    method = "TRATE",
    time.varying = TRUE,
    transition = "both",
    weighted = FALSE,
    with.missing = TRUE,
    miss.cost = c(7, 8, 9),
    miss.cost.fixed = TRUE
  )
)
timevarying <- add_case(
  timevarying,
  "trate_time_missing_derived",
  "method=TRATE;time_varying=TRUE;transition=both;weighted=FALSE;with_missing=TRUE;miss_cost_fixed=FALSE",
  seqcost(
    fx$B,
    method = "TRATE",
    time.varying = TRUE,
    transition = "both",
    weighted = FALSE,
    with.missing = TRUE,
    miss.cost.fixed = FALSE
  )
)
timevarying <- add_case(
  timevarying,
  "trate_time_short_boundary",
  "method=TRATE;time_varying=TRUE;transition=both;weighted=FALSE;lag=1;width=2",
  seqcost(
    short,
    method = "TRATE",
    time.varying = TRUE,
    transition = "both",
    weighted = FALSE,
    lag = 1
  )
)
timevarying <- add_case(
  timevarying,
  "trate_time_custom_labels",
  "method=TRATE;time_varying=TRUE;transition=both;weighted=FALSE;lag=1;time_labels=age10..age15",
  seqcost(
    custom_time,
    method = "TRATE",
    time.varying = TRUE,
    transition = "both",
    weighted = FALSE,
    lag = 1
  )
)

timevarying <- add_case(
  timevarying,
  "indels_time_unweighted_object",
  "method=INDELS;time_varying=TRUE;weighted=FALSE;attached_weights=FALSE",
  seqcost(fx$A, method = "INDELS", time.varying = TRUE, weighted = FALSE,
          check.costs = FALSE)
)
timevarying <- add_case(
  timevarying,
  "indels_time_weighted_object_flag_false",
  "method=INDELS;time_varying=TRUE;weighted=FALSE;attached_weights=TRUE;weighted_flag_ignored=TRUE",
  seqcost(fx$Aw, method = "INDELS", time.varying = TRUE, weighted = FALSE,
          check.costs = FALSE)
)
timevarying <- add_case(
  timevarying,
  "indels_time_weighted_object_flag_true",
  "method=INDELS;time_varying=TRUE;weighted=TRUE;attached_weights=TRUE;weighted_flag_ignored=TRUE",
  seqcost(fx$Aw, method = "INDELS", time.varying = TRUE, weighted = TRUE,
          check.costs = FALSE)
)
timevarying <- add_case(
  timevarying,
  "indels_time_missing_weighted",
  "method=INDELS;time_varying=TRUE;weighted=TRUE;with_missing=TRUE",
  seqcost(
    fx$Bw,
    method = "INDELS",
    time.varying = TRUE,
    weighted = TRUE,
    with.missing = TRUE
  )
)
timevarying <- add_case(
  timevarying,
  "indels_time_ragged_missing_weighted",
  "method=INDELS;time_varying=TRUE;weighted=TRUE;with_missing=TRUE;void_excluded_by_position=TRUE",
  seqcost(
    fx$Cw,
    method = "INDELS",
    time.varying = TRUE,
    weighted = TRUE,
    with.missing = TRUE,
    check.costs = FALSE
  )
)
timevarying <- add_case(
  timevarying,
  "indels_time_absent_states",
  "method=INDELS;time_varying=TRUE;absent_states=B|C;finite_indel=double_max;overflow_sm=Inf",
  seqcost(absent, method = "INDELS", time.varying = TRUE, check.costs = FALSE)
)

timevarying <- add_case(
  timevarying,
  "indelslog_time_unweighted_object",
  "method=INDELSLOG;time_varying=TRUE;weighted=FALSE;attached_weights=FALSE",
  seqcost(fx$A, method = "INDELSLOG", time.varying = TRUE, weighted = FALSE)
)
timevarying <- add_case(
  timevarying,
  "indelslog_time_weighted_object_flag_false",
  "method=INDELSLOG;time_varying=TRUE;weighted=FALSE;attached_weights=TRUE;weighted_flag_ignored=TRUE",
  seqcost(fx$Aw, method = "INDELSLOG", time.varying = TRUE, weighted = FALSE)
)
timevarying <- add_case(
  timevarying,
  "indelslog_time_weighted_object_flag_true",
  "method=INDELSLOG;time_varying=TRUE;weighted=TRUE;attached_weights=TRUE;weighted_flag_ignored=TRUE",
  seqcost(fx$Aw, method = "INDELSLOG", time.varying = TRUE, weighted = TRUE)
)
timevarying <- add_case(
  timevarying,
  "indelslog_time_missing_weighted",
  "method=INDELSLOG;time_varying=TRUE;weighted=TRUE;with_missing=TRUE",
  seqcost(
    fx$Bw,
    method = "INDELSLOG",
    time.varying = TRUE,
    weighted = TRUE,
    with.missing = TRUE
  )
)
timevarying <- add_case(
  timevarying,
  "indelslog_time_ragged_missing_weighted",
  "method=INDELSLOG;time_varying=TRUE;weighted=TRUE;with_missing=TRUE;void_excluded_by_position=TRUE",
  seqcost(
    fx$Cw,
    method = "INDELSLOG",
    time.varying = TRUE,
    weighted = TRUE,
    with.missing = TRUE
  )
)
timevarying <- add_case(
  timevarying,
  "indelslog_time_absent_states",
  "method=INDELSLOG;time_varying=TRUE;absent_states=B|C",
  seqcost(absent, method = "INDELSLOG", time.varying = TRUE)
)

write_golden_rows("seqcost_timevarying_cases.csv", timevarying)

# -------------------------------------------------------------------------
# Proximities. R normalizes by one global maximum, including across time.
# -------------------------------------------------------------------------

proximity <- list()

proximity <- add_case(
  proximity,
  "prox_constant_stationary",
  "method=CONSTANT;proximities=TRUE",
  seqcost(fx$A, method = "CONSTANT", proximities = TRUE)
)
proximity <- add_case(
  proximity,
  "prox_trate_weighted_stationary",
  "method=TRATE;weighted=TRUE;proximities=TRUE",
  seqcost(fx$Aw, method = "TRATE", weighted = TRUE, proximities = TRUE)
)
proximity <- add_case(
  proximity,
  "prox_trate_missing_fixed_stationary",
  "method=TRATE;weighted=FALSE;with_missing=TRUE;miss_cost=7;proximities=TRUE",
  seqcost(
    fx$B,
    method = "TRATE",
    weighted = FALSE,
    with.missing = TRUE,
    miss.cost = 7,
    proximities = TRUE
  )
)
proximity <- add_case(
  proximity,
  "prox_features_stationary",
  "method=FEATURES;proximities=TRUE",
  seqcost(
    fx$A,
    method = "FEATURES",
    state.features = features_constant,
    proximities = TRUE
  )
)
proximity <- add_case(
  proximity,
  "prox_indelslog_stationary",
  "method=INDELSLOG;weighted=FALSE;proximities=TRUE",
  seqcost(fx$A, method = "INDELSLOG", weighted = FALSE, proximities = TRUE)
)
proximity <- add_case(
  proximity,
  "prox_trate_time_both",
  "method=TRATE;time_varying=TRUE;transition=both;weighted=FALSE;global_max=TRUE;proximities=TRUE",
  seqcost(
    fx$A,
    method = "TRATE",
    time.varying = TRUE,
    transition = "both",
    weighted = FALSE,
    proximities = TRUE
  )
)
proximity <- add_case(
  proximity,
  "prox_one_state_nan",
  "method=CONSTANT;one_state=TRUE;zero_max=TRUE;proximities=TRUE",
  suppressWarnings(seqcost(one_state, method = "CONSTANT", proximities = TRUE))
)
proximity <- add_case(
  proximity,
  "prox_indels_time_absent_nan",
  "method=INDELS;time_varying=TRUE;absent_states=B|C;Inf_over_Inf=NaN;proximities=TRUE",
  suppressWarnings(
    seqcost(absent, method = "INDELS", time.varying = TRUE, proximities = TRUE,
            check.costs = FALSE)
  )
)

write_golden_rows("seqcost_proximity_cases.csv", proximity)

# -------------------------------------------------------------------------
# Oracle refusals. These are stored separately because Python maps R's
# simpleError to TypeError/ValueError while preserving the behavioral refusal.
# -------------------------------------------------------------------------

error_lines <- "case|class|message"

capture_error <- function(case, expression) {
  condition <- tryCatch(
    {
      suppressMessages(eval.parent(substitute(expression)))
      NULL
    },
    error = function(err) err
  )
  if (is.null(condition)) stop("expected R error for case: ", case)
  message <- gsub("[\r\n]+", " ", conditionMessage(condition))
  paste(case, paste(class(condition), collapse = ","), message, sep = "|")
}

error_lines <- c(error_lines, capture_error("non_stslist", seqcost(absent_raw, "CONSTANT")))
error_lines <- c(error_lines, capture_error("missing_method", seqcost(fx$A)))
error_lines <- c(error_lines, capture_error("invalid_method", seqcost(fx$A, "constant")))
error_lines <- c(
  error_lines,
  capture_error("trate_one_state_stationary", seqcost(one_state, "TRATE"))
)
error_lines <- c(
  error_lines,
  capture_error(
    "constant_one_state_timevarying",
    seqcost(one_state, "CONSTANT", time.varying = TRUE)
  )
)
error_lines <- c(
  error_lines,
  capture_error(
    "trate_one_state_timevarying",
    seqcost(one_state, "TRATE", time.varying = TRUE)
  )
)
error_lines <- c(
  error_lines,
  capture_error("invalid_transition", seqcost(fx$A, "CONSTANT", transition = "bad"))
)
error_lines <- c(
  error_lines,
  capture_error("constant_na", seqcost(fx$A, "CONSTANT", cval = NA_real_))
)
error_lines <- c(
  error_lines,
  capture_error("features_missing", seqcost(fx$A, "FEATURES"))
)
error_lines <- c(
  error_lines,
  capture_error(
    "features_wrong_rows",
    seqcost(fx$A, "FEATURES", state.features = data.frame(x = 1:3))
  )
)
error_lines <- c(
  error_lines,
  capture_error(
    "features_zero_columns",
    seqcost(
      fx$A,
      "FEATURES",
      state.features = data.frame(row.names = seq_len(length(attr(fx$A, "alphabet"))))
    )
  )
)
error_lines <- c(
  error_lines,
  capture_error(
    "features_bad_weights",
    seqcost(fx$A, "FEATURES", state.features = features_constant, feature.weights = 1:3)
  )
)
error_lines <- c(
  error_lines,
  capture_error(
    "features_unnamed_type",
    seqcost(fx$A, "FEATURES", state.features = features_mixed, feature.type = list(3))
  )
)
error_lines <- c(
  error_lines,
  capture_error(
    "features_bad_type_column",
    seqcost(
      fx$A,
      "FEATURES",
      state.features = features_mixed,
      feature.type = list(asymm = "not_a_column")
    )
  )
)
error_lines <- c(
  error_lines,
  capture_error(
    "features_bad_binary",
    seqcost(
      fx$A,
      "FEATURES",
      state.features = features_mixed,
      feature.type = list(asymm = "category")
    )
  )
)
error_lines <- c(
  error_lines,
  capture_error("future_timevarying", seqcost(fx$A, "FUTURE", time.varying = TRUE))
)
error_lines <- c(
  error_lines,
  capture_error(
    "features_timevarying",
    seqcost(
      fx$A,
      "FEATURES",
      time.varying = TRUE,
      state.features = features_constant
    )
  )
)
for (transition in c("previous", "next", "both")) {
  error_lines <- c(
    error_lines,
    capture_error(
      paste0("trate_time_lag2_", transition),
      seqcost(
        fx$A,
        "TRATE",
        time.varying = TRUE,
        transition = transition,
        lag = 2
      )
    )
  )
}
error_lines <- c(
  error_lines,
  capture_error(
    "trate_time_negative_lag",
    seqcost(fx$A, "TRATE", time.varying = TRUE, transition = "both", lag = -1)
  )
)

writeLines(
  error_lines,
  con = file.path(out_dir, "seqcost_stationary_errors.txt"),
  sep = "\n",
  useBytes = TRUE
)

cat("SEQCOST_GOLDENS_DONE\n")
