# Golden generation: state-order canonicalization and custom-order precarity.
#
# Pins TraMineR 2.2-14 seqprecarity(seqdata, state.order=<full permutation>) on the
# void-free fixture A (alphabet A,B,C,D) for two full state-order permutations and
# both the default TRATEDSS and the RANK correction methods. The default-order
# results are already covered by make_goldens_precarity.R (seqprecarity_A*).
# Run from the repository root with the oracle library (scripts/r_goldens/README.md).
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
out_dir <- file.path("tests", "r_goldens")
if (!identical(as.character(packageVersion("TraMineR")), "2.2.14")) {
  stop("state-order goldens require TraMineR 2.2-14")
}

write_num_matrix <- function(name, m, header) {
  con <- file(file.path(out_dir, paste0(name, ".csv")), open = "wt")
  cat(header, file = con)
  txt <- apply(as.matrix(m), 1, function(r) paste(formatC(r, format = "g", digits = 15), collapse = ","))
  writeLines(txt, con)
  close(con)
}

A <- rbind(
  c("A", "A", "B", "C", "C", "C"),
  c("A", "B", "B", "B", "C", "A"),
  c("D", "D", "D", "D", "D", "D"),
  c("A", "B", "A", "B", "A", "B")
)
sdA <- seqdef(A, alphabet = c("A", "B", "C", "D"))

B <- rbind(
  c("A", "A", "*", "B", "B", "C"),
  c("*", "A", "A", "A", "*", "*"),
  c("A", "*", "A", "B", "B", "B"),
  c("C", "C", "C", "C", "C", "C")
)
sdB <- seqdef(B, alphabet = c("A", "B", "C"), missing = "*", right = NA)
weights <- c(1.5, 2, .5, 3)
sdAw <- seqdef(A, alphabet = c("A", "B", "C", "D"), weights = weights)
sdBw <- seqdef(
  B,
  alphabet = c("A", "B", "C"),
  missing = "*",
  right = NA,
  weights = weights
)
swap_bc <- function(values) {
  out <- values
  out[values == "B"] <- "C"
  out[values == "C"] <- "B"
  out
}
sdAswap <- seqdef(swap_bc(A), alphabet = c("A", "B", "C", "D"))

format_numeric <- function(value) {
  if (is.nan(value)) return("NaN")
  if (is.na(value)) return("NA")
  if (is.infinite(value)) return(if (value > 0) "Inf" else "-Inf")
  trimws(formatC(value, format = "g", digits = 15))
}

csv_quote <- function(value) {
  paste0('"', gsub('"', '""', as.character(value), fixed = TRUE), '"')
}

collapse_values <- function(values, sep = "|") {
  paste(vapply(values, format_numeric, character(1)), collapse = sep)
}

canonical_groups <- function(labels, state_order, state_equiv) {
  groups <- list()
  assigned <- character()
  equiv <- if (is.null(state_equiv)) list() else state_equiv
  for (state in state_order) {
    if (is.na(state)) next
    if (state %in% assigned) next
    containing <- Filter(function(group) state %in% group, equiv)
    group <- if (length(containing)) containing[[1]] else state
    if (length(containing) && group[[1]] %in% state_order && state != group[[1]]) next
    # Retain incomparable class members for canonical collapse; ranked membership
    # is represented separately by code_to_group=-1.
    group <- group[group %in% labels]
    groups[[length(groups) + 1L]] <- match(group, labels) - 1L
    assigned <- c(assigned, group)
  }
  groups
}

canonical_collapse <- function(seqdata, labels, state_order, state_equiv) {
  values <- as.matrix(seqdata)
  noncomp <- labels[is.na(match(labels, state_order))]
  if (ncol(values) > 1L && length(noncomp)) {
    for (column in 2:ncol(values)) {
      replace <- values[, column] %in% noncomp
      values[replace, column] <- values[replace, column - 1L]
    }
  }
  if (!is.null(state_equiv)) {
    for (group in state_equiv) {
      if (!length(group)) next
      values[values %in% group] <- group[[1]]
    }
  }
  paste(apply(values, 1L, paste, collapse = ">"), collapse = ";")
}

state_cases <- list(
  full_default = list(sd = sdA, wm = FALSE, args = list()),
  full_reverse = list(sd = sdA, wm = FALSE, args = list(state.order = c("D", "C", "B", "A"))),
  partial_omit_leading = list(sd = sdA, wm = FALSE, args = list(state.order = c("B", "C", "D"))),
  partial_omit_trailing = list(sd = sdA, wm = FALSE, args = list(state.order = c("A", "B", "C"))),
  partial_omit_middle = list(sd = sdA, wm = FALSE, args = list(state.order = c("A", "D"))),
  partial_singleton = list(sd = sdA, wm = FALSE, args = list(state.order = "A")),
  equiv_singleton = list(sd = sdA, wm = FALSE, args = list(state.equiv = list("B"))),
  equiv_multi = list(sd = sdA, wm = FALSE, args = list(state.equiv = list(c("B", "C")))),
  partial_promote_one = list(
    sd = sdA, wm = FALSE,
    args = list(state.order = c("A", "D"), state.equiv = list(c("D", "B")))
  ),
  partial_promote_multi = list(
    sd = sdA, wm = FALSE,
    args = list(state.order = c("A", "D"), state.equiv = list(c("D", "B", "C")))
  ),
  partial_equiv_noncomp = list(
    sd = sdA, wm = FALSE,
    args = list(state.order = "A", state.equiv = list(c("B", "C")))
  ),
  explicit_tied_negative = list(
    sd = sdA, wm = FALSE, args = list(stprec = c(0, .5, .5, -1))
  ),
  explicit_negative_surviving_equiv = list(
    sd = sdA, wm = FALSE,
    args = list(stprec = c(0, .2, .8, -1), state.equiv = list(c("B", "C")))
  ),
  explicit_constant = list(sd = sdA, wm = FALSE, args = list(stprec = rep(1, 4))),
  explicit_all_negative = list(sd = sdA, wm = FALSE, args = list(stprec = rep(-1, 4))),
  missing_observed = list(sd = sdB, wm = TRUE, args = list()),
  missing_observed_equiv = list(
    sd = sdB, wm = TRUE,
    args = list(state.order = c("A", "C"), state.equiv = list(c("C", "*")))
  ),
  phantom_missing = list(sd = sdA, wm = TRUE, args = list()),
  phantom_missing_equiv = list(
    sd = sdA, wm = TRUE, args = list(state.equiv = list(c("*", "A")))
  ),
  phantom_explicit_negative = list(
    sd = sdA, wm = TRUE, args = list(stprec = c(0, .2, .4, .6, -1))
  )
)

case_lines <- paste(
  "case", "fixture", "with_missing", "labels", "state_order", "state_equiv",
  "groups", "precedence", "code_to_group", "collapsed", "seqibad", "seqidegrad",
  sep = ","
)
for (case in names(state_cases)) {
  spec <- state_cases[[case]]
  call_args <- c(list(seqdata = spec$sd), spec$args, list(with.missing = spec$wm))
  precedence <- suppressMessages(do.call(seqprecstart, call_args))
  resolved_order <- attr(precedence, "state.order")
  resolved_equiv <- attr(precedence, "state.equiv")
  labels <- alphabet(spec$sd, with.missing = spec$wm)
  groups <- canonical_groups(labels, resolved_order, resolved_equiv)
  code_to_group <- rep(-1L, length(labels))
  ranked_positions <- match(resolved_order[!is.na(resolved_order)], labels) - 1L
  for (index in seq_along(groups)) {
    comparable <- groups[[index]][groups[[index]] %in% ranked_positions]
    code_to_group[comparable + 1L] <- index - 1L
  }
  group_text <- paste(vapply(groups, paste, character(1), collapse = "+"), collapse = ";")
  equiv_text <- if (is.null(resolved_equiv)) "" else {
    paste(vapply(resolved_equiv, paste, character(1), collapse = "+"), collapse = ";")
  }
  ibad <- suppressMessages(do.call(seqibad, c(list(seqdata = spec$sd), spec$args, list(
    with.missing = spec$wm
  ))))
  idegrad <- suppressMessages(do.call(seqidegrad, c(list(seqdata = spec$sd), spec$args, list(
    with.missing = spec$wm, method = "RANK", pow = FALSE
  ))))
  fields <- c(
    case,
    if (identical(spec$sd, sdA)) "A" else "B",
    if (spec$wm) "TRUE" else "FALSE",
    paste(labels, collapse = "|"),
    paste(resolved_order, collapse = "|"),
    equiv_text,
    group_text,
    collapse_values(precedence),
    paste(code_to_group, collapse = "|"),
    canonical_collapse(spec$sd, labels, resolved_order, resolved_equiv),
    collapse_values(as.vector(ibad)),
    collapse_values(as.vector(idegrad))
  )
  case_lines <- c(case_lines, paste(vapply(fields, csv_quote, character(1)), collapse = ","))
}
writeLines(case_lines, file.path(out_dir, "state_order_r2_cases.csv"), useBytes = TRUE)

capture_error <- function(case, expression) {
  condition <- tryCatch(
    {
      suppressMessages(eval.parent(substitute(expression)))
      NULL
    },
    error = function(err) err
  )
  if (is.null(condition)) stop("expected R error for case: ", case)
  message <- trimws(gsub("[\r\n]+", " ", conditionMessage(condition)))
  paste(case, message, sep = "=")
}

error_lines <- c(
  capture_error("duplicate_order", seqprecstart(sdA, state.order = c("A", "A"))),
  capture_error("unknown_order", seqprecstart(sdA, state.order = c("A", "Z"))),
  capture_error("equiv_not_list", seqprecstart(sdA, state.equiv = c("A", "B"))),
  capture_error("duplicate_equiv", seqprecstart(sdA, state.equiv = list(c("B", "B")))),
  capture_error("overlap_equiv", seqprecstart(sdA, state.equiv = list(c("A", "B"), c("B", "C")))),
  capture_error("unknown_equiv", seqprecstart(sdA, state.equiv = list(c("A", "Z")))),
  capture_error("short_stprec", seqidegrad(sdA, stprec = c(0, 1))),
  capture_error("nan_stprec", seqprecstart(sdA, stprec = c(0, 1, NaN, 2))),
  capture_error("empty_order", seqprecstart(sdA, state.order = character())),
  capture_error("null_order", seqprecstart(sdA, state.order = NULL, stprec = NULL))
)
writeLines(c("case=message", error_lines), file.path(out_dir, "state_order_r2_errors.txt"), useBytes = TRUE)

## Cross-function state-model matrix.  Every option is passed through the public
## function so this freezes the complete seqprecarity.private propagation path.
cross_cases <- list(
  prec_partial_rank = list(
    fun = seqprecarity, fixture = "A", sd = sdA, wm = FALSE,
    order = c("A", "D"), equiv = NULL,
    args = list(method = "RANK")
  ),
  prec_partial_tratedss = list(
    fun = seqprecarity, fixture = "A", sd = sdA, wm = FALSE,
    order = c("A", "D"), equiv = NULL,
    args = list(method = "TRATEDSS")
  ),
  prec_equiv_tratedss = list(
    fun = seqprecarity, fixture = "A", sd = sdA, wm = FALSE,
    order = NULL, equiv = list(c("B", "C")),
    args = list(method = "TRATEDSS")
  ),
  prec_equiv_tratedss_swapped = list(
    fun = seqprecarity, fixture = "Aswap", sd = sdAswap, wm = FALSE,
    order = NULL, equiv = list(c("B", "C")),
    args = list(method = "TRATEDSS")
  ),
  prec_partial_equiv_weighted_inv = list(
    fun = seqprecarity, fixture = "Aw", sd = sdAw, wm = FALSE,
    order = c("A", "D"), equiv = list(c("D", "B")),
    args = list(
      method = "TRATEDSS", penalized = "NEG", weight.type = "INV",
      pow = 2, otto = .35, a = .8, b = 1.1, border.effect = 5
    )
  ),
  prec_equiv_weighted_loginv = list(
    fun = seqprecarity, fixture = "Aw", sd = sdAw, wm = FALSE,
    order = NULL, equiv = list(c("B", "C")),
    args = list(
      method = "TRATEDSS", penalized = "POS", weight.type = "LOGINV",
      pow = .5, otto = .1, a = 1.3, b = .7, border.effect = 3
    )
  ),
  prec_missing_equiv_weighted = list(
    fun = seqprecarity, fixture = "Bw", sd = sdBw, wm = TRUE,
    order = c("A", "C"), equiv = list(c("C", "*")),
    args = list(
      method = "TRATEDSS", penalized = "BOTH", weight.type = "ADD",
      pow = 2, otto = .4, a = .6, b = 1.4, border.effect = 4
    )
  ),
  prec_phantom_missing_equiv = list(
    fun = seqprecarity, fixture = "A", sd = sdA, wm = TRUE,
    order = NULL, equiv = list(c("*", "A")),
    args = list(method = "RANK", penalized = "BOTH", pow = 2)
  ),
  prec_partial_scalar_correction = list(
    fun = seqprecarity, fixture = "A", sd = sdA, wm = FALSE,
    order = c("A", "D"), equiv = NULL,
    args = list(correction = .25, otto = .4, a = .2, b = 2)
  ),
  prec_equiv_vector_correction = list(
    fun = seqprecarity, fixture = "A", sd = sdA, wm = FALSE,
    order = NULL, equiv = list(c("B", "C")),
    args = list(correction = c(-.25, 0, .25, .5), otto = .3, a = .9, b = 1.1)
  ),
  insec_partial_rank = list(
    fun = seqinsecurity, fixture = "A", sd = sdA, wm = FALSE,
    order = c("A", "D"), equiv = NULL,
    args = list(method = "RANK", pow = 2, spow = 0)
  ),
  insec_partial_weighted_tratedss_bound = list(
    fun = seqinsecurity, fixture = "Aw", sd = sdAw, wm = FALSE,
    order = c("A", "D"), equiv = NULL,
    args = list(
      method = "TRATEDSS", penalized = "NEG", weight.type = "INV",
      pow = 2, spow = 1, bound = TRUE, border.effect = 3
    )
  ),
  insec_equiv_weighted_loginv = list(
    fun = seqinsecurity, fixture = "Aw", sd = sdAw, wm = FALSE,
    order = NULL, equiv = list(c("B", "C")),
    args = list(
      method = "TRATEDSS", penalized = "POS", weight.type = "LOGINV",
      pow = .5, spow = 2, border.effect = 4
    )
  ),
  insec_partial_equiv_rank_bound = list(
    fun = seqinsecurity, fixture = "A", sd = sdA, wm = FALSE,
    order = c("A", "D"), equiv = list(c("D", "B")),
    args = list(method = "RANK", penalized = "BOTH", pow = 1.5, spow = .5, bound = TRUE)
  ),
  insec_missing_equiv_weighted = list(
    fun = seqinsecurity, fixture = "Bw", sd = sdBw, wm = TRUE,
    order = c("A", "C"), equiv = list(c("C", "*")),
    args = list(
      method = "TRATEDSS", penalized = "BOTH", weight.type = "ADD",
      pow = 2, spow = 0, bound = TRUE, border.effect = 5
    )
  ),
  insec_phantom_missing_equiv = list(
    fun = seqinsecurity, fixture = "A", sd = sdA, wm = TRUE,
    order = NULL, equiv = list(c("*", "A")),
    args = list(method = "RANK", pow = 2, spow = 1, bound = TRUE)
  ),
  insec_partial_vector_correction = list(
    fun = seqinsecurity, fixture = "A", sd = sdA, wm = FALSE,
    order = c("A", "D"), equiv = NULL,
    args = list(correction = c(-.2, 0, .2, .4), pow = 2, spow = 1, bound = TRUE)
  )
)

arg_value <- function(args, name, default) {
  if (name %in% names(args)) args[[name]] else default
}
cross_lines <- paste(
  "case", "function", "fixture", "row", "weighted", "with_missing",
  "state_order", "state_equiv", "method", "penalized", "weight_type",
  "pow", "spow", "bound", "correction", "otto", "a", "b", "border_effect", "value",
  sep = ","
)
for (case in names(cross_cases)) {
  spec <- cross_cases[[case]]
  call_args <- c(
    list(seqdata = spec$sd),
    if (is.null(spec$order)) list() else list(state.order = spec$order),
    if (is.null(spec$equiv)) list() else list(state.equiv = spec$equiv),
    list(with.missing = spec$wm),
    spec$args
  )
  value <- suppressMessages(do.call(spec$fun, call_args))
  function_name <- if (identical(spec$fun, seqprecarity)) "seqprecarity" else "seqinsecurity"
  equiv_text <- if (is.null(spec$equiv)) "" else {
    paste(vapply(spec$equiv, paste, character(1), collapse = "+"), collapse = ";")
  }
  correction <- arg_value(spec$args, "correction", NULL)
  correction_text <- if (is.null(correction)) "" else collapse_values(correction)
  for (row in seq_along(value)) {
    fields <- c(
      case,
      function_name,
      spec$fixture,
      row,
      if (grepl("w$", spec$fixture)) "TRUE" else "FALSE",
      if (spec$wm) "TRUE" else "FALSE",
      if (is.null(spec$order)) "" else paste(spec$order, collapse = "|"),
      equiv_text,
      arg_value(spec$args, "method", if (function_name == "seqprecarity") "TRATEDSS" else "RANK"),
      arg_value(spec$args, "penalized", "BOTH"),
      arg_value(spec$args, "weight.type", "ADD"),
      format_numeric(arg_value(spec$args, "pow", 1)),
      format_numeric(arg_value(spec$args, "spow", arg_value(spec$args, "pow", 1))),
      if (isTRUE(arg_value(spec$args, "bound", FALSE))) "TRUE" else "FALSE",
      correction_text,
      format_numeric(arg_value(spec$args, "otto", .2)),
      format_numeric(arg_value(spec$args, "a", 1)),
      format_numeric(arg_value(spec$args, "b", 1.2)),
      format_numeric(arg_value(spec$args, "border.effect", 10)),
      format_numeric(as.vector(value)[[row]])
    )
    cross_lines <- c(cross_lines, paste(vapply(fields, csv_quote, character(1)), collapse = ","))
  }
}
writeLines(cross_lines, file.path(out_dir, "precarity_order_r2_cases.csv"), useBytes = TRUE)

ord_rev <- c("D", "C", "B", "A")
ord_mix <- c("C", "A", "D", "B")

write_num_matrix("seqprecarity_A_ord_rev", as.vector(seqprecarity(sdA, state.order = ord_rev)),
  "# R: seqprecarity(sdA, state.order=c('D','C','B','A')); TRATEDSS default\n")
write_num_matrix("seqprecarity_A_ord_mix", as.vector(seqprecarity(sdA, state.order = ord_mix)),
  "# R: seqprecarity(sdA, state.order=c('C','A','D','B')); TRATEDSS default\n")
write_num_matrix("seqprecarity_A_ord_rev_rank",
  as.vector(seqprecarity(sdA, state.order = ord_rev, method = "RANK")),
  "# R: seqprecarity(sdA, state.order=c('D','C','B','A'), method='RANK')\n")

# Also pin seqprecstart with the custom order (its state.order handling underpins the above).
write_num_matrix("seqprecstart_A_ord_rev", t(as.vector(seqprecstart(sdA, state.order = ord_rev))),
  "# R: seqprecstart(sdA, state.order=c('D','C','B','A')); states A,B,C,D order\n")

cat("PRECARITY_ORDER_GOLDENS_DONE\n")
