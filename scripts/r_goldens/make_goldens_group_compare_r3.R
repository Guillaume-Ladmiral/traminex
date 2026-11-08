# Golden generation for TraMineRextras 0.6.9 dissCompare sampling routes.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineRextras))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
stopifnot(as.character(packageVersion("TraMineRextras")) == "0.6.9")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

constant_diss <- function(n) {
  result <- matrix(1, nrow = n, ncol = n)
  diag(result) <- 0
  result
}

diss <- constant_diss(12)
group <- c(rep("a", 7), rep("b", 5))
weights <- c(1, 2, 1.5, 0.5, 3, 1, 2.5, 1, 4, 0.75, 2, 1.25)

cases <- list(
  s0_unweighted_default = dissCompare(diss, group, s = 0, weighted = FALSE),
  s0_weighted_default = dissCompare(diss, group, s = 0, weights = weights),
  s0_by_group = dissCompare(
    diss,
    group,
    s = 0,
    weighted = "by.group",
    weights = weights
  ),
  s0_squared_false = dissCompare(diss, group, s = 0, squared = FALSE),
  s0_squared_true = dissCompare(diss, group, s = 0, squared = TRUE),
  s0_lrt_inc_df = dissCompare(diss, group, s = 0, stat = "LRT", inc.df = TRUE),
  s0_bic_inc_df = dissCompare(diss, group, s = 0, stat = "BIC", inc.df = TRUE),
  sample_below_both_bf = dissCompare(diss, group, s = 3, seed = 1401, BFopt = NULL),
  sample_below_avg_bf = dissCompare(diss, group, s = 3, seed = 1401, BFopt = 1),
  sample_below_from_bic = dissCompare(diss, group, s = 3, seed = 1401, BFopt = 2),
  sample_at_min = dissCompare(diss, group, s = 5, seed = 1402, inc.df = TRUE),
  sample_at_max = dissCompare(diss, group, s = 7, seed = 1403, inc.df = TRUE),
  sample_above_max = dissCompare(diss, group, s = 9, seed = 1404, inc.df = TRUE),
  sample_weighted_unit = dissCompare(
    diss,
    group,
    s = 3,
    seed = 1405,
    weighted = TRUE,
    weights = rep(1, 12)
  ),
  sample_by_group_unit = dissCompare(
    diss,
    group,
    s = 3,
    seed = 1405,
    weighted = "by.group",
    weights = rep(1, 12)
  )
)

multi_group <- c("a", "a", "a", "a", "b", "b", "b", "a", "a", "a", "b", "b")
multi_set <- c(rep("early", 7), rep("late", 5))
cases$multiple_sets_s0 <- dissCompare(diss, multi_group, set = multi_set, s = 0)
cases$multiple_sets_sample <- dissCompare(
  diss,
  multi_group,
  set = multi_set,
  s = 2,
  seed = 1406,
  inc.df = TRUE
)

mixed_case_diss <- constant_diss(8)
mixed_case_group <- rep(c("g1", "g2"), 4)
mixed_case_set <- c(rep("a", 4), rep("A", 4))
cases$mixed_case_sets <- dissCompare(
  mixed_case_diss,
  mixed_case_group,
  set = mixed_case_set,
  s = 0,
  weighted = FALSE
)

locale_set_levels <- c("_", "-", "!", "?", ".", "*", "%", "a", "A", "á", "ä", "b")
locale_case_diss <- constant_diss(4 * length(locale_set_levels))
cases$locale_case_sets <- dissCompare(
  locale_case_diss,
  rep(c("g1", "g2"), 2 * length(locale_set_levels)),
  set = rep(locale_set_levels, each = 4),
  s = 0,
  weighted = FALSE
)

# Sequence-wrapper routes use one unique repeated state per sequence. All
# off-diagonal LCS distances are therefore exchangeable, which freezes the
# exact aggregation independently of R-versus-NumPy random-stream bytes.
state_names <- paste0("S", sprintf("%02d", 1:12))
seq_matrix <- matrix(rep(state_names, each = 2), nrow = 12, byrow = TRUE)
seqdata <- suppressMessages(
  seqdef(
    as.data.frame(seq_matrix, stringsAsFactors = FALSE),
    alphabet = state_names,
    weights = weights
  )
)
group_a <- 1:7
group_b <- 8:12
seq1 <- list(seqdata[c(1, 2, 3, 4), ], seqdata[c(5, 6, 7), ])
seq2 <- list(seqdata[c(8, 9, 10), ], seqdata[c(11, 12), ])

cases$seq_group_s0 <- seqCompare(
  seqdata, group = group, s = 0, method = "LCS"
)
cases$seq_group_s0_unweighted <- seqCompare(
  seqdata, group = group, s = 0, weighted = FALSE, method = "LCS"
)
cases$seq_group_s0_by_group <- seqCompare(
  seqdata, group = group, s = 0, weighted = "by.group", method = "LCS"
)
cases$seq_group_sample_opt1 <- seqCompare(
  seqdata, group = group, s = 3, seed = 1501, opt = 1,
  weighted = FALSE, method = "LCS"
)
cases$seq_group_sample_opt2 <- seqCompare(
  seqdata, group = group, s = 3, seed = 1501, opt = 2,
  weighted = FALSE, method = "LCS"
)
cases$seq_group_sample_auto1 <- seqCompare(
  seqdata, group = group, s = 3, seed = 1501,
  weighted = FALSE, method = "LCS"
)
cases$seq_group_sample_auto2 <- seqCompare(
  seqdata, group = group, s = 7, seed = 1502,
  weighted = FALSE, method = "LCS"
)
cases$seq_group_sets_sample <- seqCompare(
  seqdata, group = multi_group, set = multi_set, s = 2, seed = 1503,
  inc.df = TRUE, weighted = FALSE, method = "LCS"
)
cases$seq_two_object_s0 <- seqCompare(
  seqdata[group_a, ], seqdata[group_b, ], s = 0, method = "LCS"
)
cases$seq_two_object_sample <- seqCompare(
  seqdata[group_a, ], seqdata[group_b, ], s = 3, seed = 1504,
  opt = 2, BFopt = 2, weighted = FALSE, method = "LCS"
)
cases$seq_paired_lists_s0 <- seqCompare(
  seq1, seq2, s = 0, method = "LCS"
)
cases$seq_paired_lists_sample <- seqCompare(
  seq1, seq2, s = 2, seed = 1505, BFopt = 1, inc.df = TRUE,
  weighted = FALSE, method = "LCS"
)
cases$seq_lrt_wrapper <- seqLRT(
  seqdata, group = group, s = 3, seed = 1506, opt = 2,
  inc.df = TRUE, weighted = FALSE, method = "LCS"
)
cases$seq_bic_wrapper <- seqBIC(
  seqdata, group = group, s = 3, seed = 1506, opt = 2,
  BFopt = 2, inc.df = TRUE, weighted = FALSE, method = "LCS"
)

na_group <- multi_group
na_set <- multi_set
na_group[[2]] <- NA
na_set[[9]] <- NA
captured_messages <- character()
cases$na_removed <- withCallingHandlers(
  dissCompare(diss, na_group, set = na_set, s = 0),
  message = function(condition) {
    captured_messages <<- c(captured_messages, conditionMessage(condition))
    invokeRestart("muffleMessage")
  }
)

records <- list()
record_index <- 1L
for (case_name in names(cases)) {
  result <- cases[[case_name]]
  if (is.null(dim(result))) {
    result <- matrix(
      result,
      nrow = 1,
      dimnames = list("1", names(result))
    )
  }
  for (row in seq_len(nrow(result))) {
    for (column in seq_len(ncol(result))) {
      row_label <- if (is.null(rownames(result))) as.character(row) else rownames(result)[[row]]
      records[[record_index]] <- data.frame(
        case = case_name,
        record_type = "result",
        row = row_label,
        column = colnames(result)[[column]],
        value = format(result[row, column], digits = 17, scientific = FALSE),
        stringsAsFactors = FALSE
      )
      record_index <- record_index + 1L
    }
  }
}

# Freeze the upstream first-group/minority indexing quirk separately from RNG bytes.
reverse_group <- c(rep("a", 3), rep("b", 5))
m.n <- 5
n.n <- 3
s <- 3
f.n1 <- floor(s / m.n)
ff.n1 <- max(1, f.n1)
r.n1 <- ifelse(s < m.n, s - m.n %% s, s - f.n1 * m.n)
k.n <- floor((ff.n1 * m.n + r.n1) / n.n)
r.n2 <- (ff.n1 * m.n + r.n1) - k.n * n.n
set.seed(1410)
r.s1 <- c(sample(rep(1:m.n, ff.n1)), sample(1:m.n, r.n1, FALSE))
r.s2 <- c(sample(rep(1:n.n, k.n)), sample(1:n.n, r.n2, FALSE))
r.s1 <- matrix(r.s1, ncol = s)
r.s2 <- matrix(r.s2, ncol = s)
plan_facts <- c(
  first_group_size = sum(reverse_group == "a"),
  second_group_size = sum(reverse_group == "b"),
  sample_size = s,
  draw_count = nrow(r.s1),
  first_relative_max = max(r.s1),
  second_relative_max = max(r.s2)
)
for (fact_name in names(plan_facts)) {
  records[[record_index]] <- data.frame(
    case = "reversed_unbalanced_plan",
    record_type = "plan",
    row = "1",
    column = fact_name,
    value = as.character(plan_facts[[fact_name]]),
    stringsAsFactors = FALSE
  )
  record_index <- record_index + 1L
}

write.csv(
  do.call(rbind, records),
  file.path(out_dir, "group_compare_r3_cases.csv"),
  row.names = FALSE,
  quote = TRUE
)

condition_lines <- vapply(
  captured_messages,
  function(message) paste("na_removed", "message", trimws(message), sep = "|"),
  character(1)
)
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

capture_condition("set_without_group", dissCompare(diss, NULL, set = multi_set, s = 0))
capture_condition("group_length", dissCompare(diss, group[-1], s = 0))
capture_condition("one_group", dissCompare(diss, rep("a", 12), s = 0))
capture_condition("three_groups", dissCompare(diss, rep(c("a", "b", "c"), 4), s = 0))
capture_condition("weighted_value", dissCompare(diss, group, weighted = "bad", s = 0))
capture_condition("weight_length", dissCompare(diss, group, weights = weights[-1], s = 0))
capture_condition("squared_value", dissCompare(diss, group, squared = "bad", s = 0))
capture_condition("stat_value", dissCompare(diss, group, stat = "bad", s = 0))
capture_condition("negative_s", dissCompare(diss, group, s = -1))
bfopt_three <- capture_condition("bfopt_three", dissCompare(diss, group, s = 3, BFopt = 3))
if (!is.null(bfopt_three)) {
  condition_lines <- c(condition_lines, "bfopt_three|success")
}
capture_condition("seq_both_null", seqCompare(seqdata, s = 0, method = "LCS"))
capture_condition(
  "seq_set_without_group",
  seqCompare(seqdata, set = multi_set, s = 0, method = "LCS")
)
capture_condition(
  "seq_list_length",
  seqCompare(seq1, seq2[1], s = 0, method = "LCS")
)
capture_condition(
  "seq_list_element",
  seqCompare(list(seq1[[1]], data.frame(x = 1)), seq2, s = 0, method = "LCS")
)
capture_condition(
  "seq_single_bad_second",
  seqCompare(seqdata, data.frame(x = 1), s = 0, method = "LCS")
)
capture_condition(
  "seq_bad_stat",
  seqCompare(seqdata, group = group, stat = "bad", s = 0, method = "LCS")
)
bad_opt <- capture_condition(
  "seq_bad_opt",
  seqCompare(seqdata, group = group, s = 3, opt = 9, method = "LCS")
)
if (!is.null(bad_opt)) {
  condition_lines <- c(condition_lines, "seq_bad_opt|success")
}

writeLines(
  condition_lines,
  file.path(out_dir, "group_compare_r3_errors.txt"),
  useBytes = TRUE
)

cat("GROUP_COMPARE_R3_GOLDENS_DONE\n")
