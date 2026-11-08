# Golden generation: the complete public seqtrate option matrix.
#
# Pins stationary and time-varying transition probabilities/counts to
# TraMineR 2.2-14, including selected states, signed lag, weights, missing and
# void handling, empty source rows, one-transition input, and R dimension names.
#
# Run from the repository root with the oracle library (scripts/r_goldens/README.md):
#   Rscript scripts/r_goldens/make_goldens_seqtrate.R
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
source(file.path("scripts", "r_goldens", "make_goldens_corestats.R"))

if (as.character(packageVersion("TraMineR")) != "2.2.14") {
  stop("seqtrate goldens require TraMineR 2.2-14")
}

emit <- function(sd, tag, wm) {
  m <- seqtrate(sd, weighted = TRUE, lag = 1, with.missing = wm)
  write_num_matrix(
    paste0("seqtrate_", tag), matrix(as.numeric(m), nrow = nrow(m)),
    sprintf("# R: seqtrate(%s, weighted=TRUE, lag=1, with.missing=%s) KxK rates\n", tag, wm)
  )
}

fx <- corestats_fixtures()
emit(fx$A, "A", FALSE)
emit(fx$Aw, "Aw", FALSE)
emit(fx$B, "B_wm0", FALSE)
emit(fx$B, "B_wm1", TRUE)
emit(fx$C, "C_wm0", FALSE)
emit(fx$C, "C_wm1", TRUE)

raw_states <- function(sd, sel.states, with.missing) {
  states <- if (is.null(sel.states)) attr(sd, "alphabet") else sel.states
  if (with.missing) states <- c(states, attr(sd, "nr"))
  as.character(states)
}

format_values <- function(x) formatC(as.numeric(x), format = "g", digits = 15)

stationary_rows <- function(case, sd, sel.states = NULL, weighted = TRUE,
                            lag = 1, with.missing = FALSE, count = FALSE) {
  value <- suppressMessages(seqtrate(
    sd, sel.states = sel.states, weighted = weighted, lag = lag,
    with.missing = with.missing, count = count
  ))
  states <- raw_states(sd, sel.states, with.missing)
  rows <- vector("list", length(states) * length(states))
  at <- 1L
  for (i in seq_along(states)) {
    for (j in seq_along(states)) {
      rows[[at]] <- data.frame(
        case = case,
        from_state = states[[i]],
        to_state = states[[j]],
        value = format_values(value[i, j]),
        stringsAsFactors = FALSE
      )
      at <- at + 1L
    }
  }
  do.call(rbind, rows)
}

timevarying_rows <- function(case, sd, sel.states = NULL, weighted = TRUE,
                             lag = 1, with.missing = FALSE, count = FALSE) {
  value <- suppressMessages(seqtrate(
    sd, sel.states = sel.states, time.varying = TRUE, weighted = weighted,
    lag = lag, with.missing = with.missing, count = count
  ))
  states <- raw_states(sd, sel.states, with.missing)
  time.labels <- dimnames(value)[[3]]
  rows <- vector("list", length(states) * length(states) * length(time.labels))
  at <- 1L
  for (t in seq_along(time.labels)) {
    for (i in seq_along(states)) {
      for (j in seq_along(states)) {
        rows[[at]] <- data.frame(
          case = case,
          from_state = states[[i]],
          to_state = states[[j]],
          time_label = time.labels[[t]],
          value = format_values(value[i, j, t]),
          stringsAsFactors = FALSE
        )
        at <- at + 1L
      }
    }
  }
  do.call(rbind, rows)
}

absent_source <- rbind(
  c("A", "A", "B", "C"),
  c("B", "A", "B", "C")
)
sd_absent_source <- seqdef(absent_source, alphabet = c("A", "B", "C"))

one_transition <- matrix(c("A", "B"), nrow = 1L)
sd_one_transition <- seqdef(one_transition, alphabet = c("A", "B", "C"))

missing_only <- matrix(c("A", "*"), nrow = 1L)
sd_missing_only <- seqdef(
  missing_only, alphabet = c("A"), missing = "*", right = NA
)

stationary <- list(
  stationary_rows("a_unweighted_prob_lag1", fx$Aw, weighted = FALSE),
  stationary_rows("a_weighted_prob_lag1", fx$Aw, weighted = TRUE),
  stationary_rows("a_weighted_count_lag1", fx$Aw, weighted = TRUE, count = TRUE),
  stationary_rows(
    "a_subset_ba_prob_lag1", fx$A, sel.states = c("B", "A"), weighted = FALSE
  ),
  stationary_rows("a_unweighted_prob_lag2", fx$A, weighted = FALSE, lag = 2),
  stationary_rows("a_unweighted_prob_lag_neg1", fx$A, weighted = FALSE, lag = -1),
  stationary_rows("b_prob_lag1_wm0", fx$B, weighted = FALSE),
  stationary_rows(
    "b_prob_lag1_wm1", fx$B, weighted = FALSE, with.missing = TRUE
  ),
  stationary_rows("c_prob_lag1_wm0", fx$C, weighted = FALSE),
  stationary_rows(
    "c_prob_lag1_wm1", fx$C, weighted = FALSE, with.missing = TRUE
  ),
  stationary_rows(
    "a_phantom_missing_lag1", fx$A, weighted = FALSE, with.missing = TRUE
  ),
  stationary_rows("absent_source_prob_lag1", sd_absent_source, weighted = FALSE),
  stationary_rows("one_transition_prob_lag1", sd_one_transition, weighted = FALSE),
  stationary_rows(
    "one_transition_count_lag1", sd_one_transition, weighted = FALSE, count = TRUE
  ),
  stationary_rows(
    "missing_only_selected_count_wm0", sd_missing_only,
    sel.states = c("A", "*"), weighted = FALSE, count = TRUE
  )
)
write.csv(
  do.call(rbind, stationary),
  file.path(out_dir, "seqtrate_full_stationary.csv"),
  row.names = FALSE,
  quote = TRUE
)

custom <- rbind(
  c("A", "A", "B", "C", "C", "C"),
  c("A", "B", "B", "B", "C", "A"),
  c("D", "D", "D", "D", "D", "D"),
  c("A", "B", "A", "B", "A", "B")
)
colnames(custom) <- paste0("age", 10:15)
sd_custom <- seqdef(custom, alphabet = c("A", "B", "C", "D"))

timevarying <- list(
  timevarying_rows("a_tv_unweighted_prob_lag1", fx$Aw, weighted = FALSE),
  timevarying_rows("a_tv_weighted_prob_lag1", fx$Aw, weighted = TRUE),
  timevarying_rows(
    "a_tv_weighted_count_lag1", fx$Aw, weighted = TRUE, count = TRUE
  ),
  timevarying_rows(
    "a_tv_subset_ba_lag1", fx$A, sel.states = c("B", "A"), weighted = FALSE
  ),
  timevarying_rows("a_tv_unweighted_prob_lag2", fx$A, weighted = FALSE, lag = 2),
  timevarying_rows(
    "b_tv_prob_lag1_wm1", fx$B, weighted = FALSE, with.missing = TRUE
  ),
  timevarying_rows("c_tv_prob_lag1_wm0", fx$C, weighted = FALSE),
  timevarying_rows("one_tv_prob_lag1", sd_one_transition, weighted = FALSE),
  timevarying_rows("a_tv_custom_labels_lag2", sd_custom, weighted = FALSE, lag = 2)
)
write.csv(
  do.call(rbind, timevarying),
  file.path(out_dir, "seqtrate_full_timevarying.csv"),
  row.names = FALSE,
  quote = TRUE
)

negative_error <- tryCatch(
  {
    suppressMessages(seqtrate(fx$A, time.varying = TRUE, lag = -1))
    NULL
  },
  error = function(e) e
)
if (is.null(negative_error)) stop("negative time-varying lag unexpectedly succeeded")
writeLines(
  c(
    paste0("class=", class(negative_error)[[1]]),
    paste0("message=", conditionMessage(negative_error))
  ),
  file.path(out_dir, "seqtrate_full_timevarying_errors.txt")
)
cat("SEQTRATE_GOLDENS_DONE\n")
