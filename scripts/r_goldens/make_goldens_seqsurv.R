# Golden generation: RP6b seqsurv (weighted Kaplan-Meier spell survival).
# Anchors traminex's seqsurv to TraMineRextras 0.6.9 (uses survival::survfit).
# Run in the live-R env:  Rscript scripts/r_goldens/make_goldens_seqsurv.R
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
suppressMessages(library(TraMineRextras))
suppressMessages(library(survival))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
stopifnot(as.character(packageVersion("TraMineRextras")) == "0.6.9")

# `survival` is a *recommended* package that ships with the R installation, not one
# of the two retained CRAN archives, so the isolated oracle library documented in
# scripts/r_goldens/README.md structurally cannot pin it: with R_LIBS set to that
# library, `survival` still resolves from the system library. Recording its resolved
# version as fixture *data* therefore made this fixture reproducible only on a machine
# whose `survival` matched exactly, and the documented two-pass provenance proof failed
# on any other (cycle-4b F-002). The version is asserted here, where a mismatch stops
# generation, and reported to stdout; `check_r_golden_provenance.py` records it in the
# evidence JSON, where a change is informative but is not a fixture byte-diff.
#
# Compatible range: survival >= 3.8.9. The Kaplan-Meier estimator, its Greenwood
# standard error and the infinitesimal-jackknife (robust) standard error used for
# fractional weights are stable across this range; 3.8.9 produced the reviewed values
# and 3.8.11 reproduces every analytic cell byte-identically (verified: regenerating
# under 3.8.11 changed 20 of 20 differing lines, all of them this version field, and
# zero analytic values).
stopifnot(packageVersion("survival") >= "3.8.9")
cat("survival version:", as.character(packageVersion("survival")), "\n")

out_dir <- file.path("tests", "r_goldens")

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
sdAw <- seqdef(A, alphabet = c("A", "B", "C", "D"), weights = c(1.5, 2, 0.5, 3))

# Stacked summary: stratum code (1-based alphabet), time, n.risk, n.event, surv, std.err.
# Unweighted std.err is the textbook Greenwood (NaN at surv=0). For fractional weights R
# sets robust=TRUE, so the weighted std.err is the infinitesimal-jackknife SE (0 at surv=0);
# both are pinned (the weighted SE is now reproduced exactly -- audit F6).
surv_mat <- function(sd, with_se) {
  fit <- suppressMessages(seqsurv(sd))
  s <- summary(fit)
  strata_lab <- sub("spell\\$states=", "", as.character(s$strata))
  out <- cbind(match(strata_lab, alphabet(sd)), s$time, s$n.risk, s$n.event, s$surv)
  if (with_se) out <- cbind(out, s$std.err)
  out
}

write_num_matrix("seqsurv_A", surv_mat(sdA, TRUE),
  "# R: summary(seqsurv(sdA)); stratum,time,n.risk,n.event,surv,std.err (Greenwood)\n")
write_num_matrix("seqsurv_A_w", surv_mat(sdAw, TRUE),
  "# R: summary(seqsurv(sdAw weighted)); stratum,time,n.risk,n.event,surv,std.err (robust IJ)\n")

cat("SEQSURV_GOLDENS_DONE\n")

## R2: complete state/group route matrix and survfit confidence contract ------
r2_rows <- list()
r2_add <- function(case, section, curve = "", row = 0L, time = NaN,
                   n.risk = NaN, n.event = NaN, surv = NaN, std.err = NaN,
                   lower = NaN, upper = NaN, value = NaN, text = "") {
  num <- function(x) {
    if (is.nan(x)) "NaN" else if (is.na(x)) "NA" else
      formatC(as.numeric(x), digits = 15, format = "fg", flag = "#")
  }
  r2_rows[[length(r2_rows) + 1L]] <<- c(
    case = as.character(case), section = as.character(section), curve = as.character(curve),
    row = as.character(row), time = num(time), n_risk = num(n.risk),
    n_event = num(n.event), surv = num(surv), std_err = num(std.err),
    lower = num(lower), upper = num(upper), value = num(value),
    text = gsub("[\r\n,]", " ", as.character(text))
  )
}

r2_factor_levels <- function(groups, n) {
  if (is.null(groups)) levels(factor(rep("1", n))) else levels(factor(groups))
}

r2_selected_states <- function(seqdata, state) {
  obs <- seqstatl(seqdata)
  selected <- if (is.null(state)) obs else state[state %in% obs]
  alphabet(seqdata)[alphabet(seqdata) %in% selected]
}

r2_record <- function(case, fit, seqdata, groups, per.state, state) {
  selected <- r2_selected_states(seqdata, state)
  group.levels <- r2_factor_levels(groups, nrow(seqdata))
  curve.axis <- if (per.state) "group" else "state"
  curve.keys <- if (per.state) group.levels else selected
  for (i in seq_along(curve.keys)) r2_add(case, "curve", curve.keys[[i]], i)
  for (i in seq_along(selected)) r2_add(case, "state", row = i, text = selected[[i]])
  for (i in seq_along(group.levels))
    r2_add(case, "group", row = i, text = group.levels[[i]])
  r2_add(case, "curve_axis", text = curve.axis)
  # The resolved `survival` version is asserted at the top of this script and recorded
  # in the provenance evidence JSON, not written here as fixture data (F-002).
  r2_add(case, "xtstep", value = attr(fit, "xtstep"))
  r2_add(case, "tick_last", value = as.integer(isTRUE(attr(fit, "tick.last"))))
  ltext <- attr(fit, "ltext")
  for (i in seq_along(ltext)) r2_add(case, "ltext", row = i, text = ltext[[i]])
  palette <- attr(fit, "cpal")
  if (is.null(palette)) r2_add(case, "cpal", text = "NULL")
  else for (i in seq_along(palette)) r2_add(case, "cpal", row = i, text = palette[[i]])

  summary.fit <- summary(fit)
  if (length(summary.fit$time) == 0L) return(invisible(NULL))
  if (is.null(summary.fit$strata)) {
    curve.names <- rep(curve.keys[[1]], length(summary.fit$time))
  } else {
    prefix <- if (per.state) "spell$group=" else "spell$states="
    curve.names <- sub(prefix, "", as.character(summary.fit$strata), fixed = TRUE)
  }
  counters <- setNames(integer(length(curve.keys)), curve.keys)
  for (i in seq_along(summary.fit$time)) {
    key <- curve.names[[i]]
    counters[[key]] <- counters[[key]] + 1L
    r2_add(
      case, "survival", key, counters[[key]], summary.fit$time[[i]],
      summary.fit$n.risk[[i]], summary.fit$n.event[[i]], summary.fit$surv[[i]],
      summary.fit$std.err[[i]], summary.fit$lower[[i]], summary.fit$upper[[i]]
    )
  }
}

r2_messages <- c("case=message")
r2_case <- function(case, seqdata, groups = NULL, per.state = FALSE,
                    state = NULL, with.missing = FALSE) {
  messages <- character()
  fit <- withCallingHandlers(
    seqsurv(
      seqdata, groups = groups, per.state = per.state, state = state,
      with.missing = with.missing
    ),
    message = function(m) {
      text <- conditionMessage(m)
      if (grepl("isn't implemented yet", text, fixed = TRUE))
        messages <<- c(messages, text)
      invokeRestart("muffleMessage")
    }
  )
  r2_record(case, fit, seqdata, groups, per.state, state)
  for (message in messages)
    r2_messages <<- c(
      r2_messages,
      paste0(case, "=", trimws(gsub("[\r\n]", " ", message)))
    )
}

R2 <- rbind(
  c("A", "A", "B", "B", "C", "C"),
  c("A", "B", "B", "B", "C", "A"),
  c("D", "D", "D", "D", "D", "D"),
  c("A", "B", "A", "B", "A", "B"),
  c("A", "A", NA, "C", NA, NA),
  c("B", "B", "C", "C", "C", "B")
)
r2_labels <- c("Alpha", "Beta", "Gamma", "No exit")
r2_cpal <- c("#111111", "#333333", "#777777", "#AAAAAA")
r2_weights <- c(1, 2.5, 0.5, 3, 1.5, 2)
r2_sd <- suppressMessages(seqdef(
  R2, alphabet = c("A", "B", "C", "D"), labels = r2_labels, cpal = r2_cpal,
  xtstep = 2, tick.last = TRUE
))
r2_sdw <- suppressMessages(seqdef(
  R2, alphabet = c("A", "B", "C", "D"), labels = r2_labels, cpal = r2_cpal,
  weights = r2_weights, xtstep = 2, tick.last = TRUE
))
r2_groups <- factor(
  c("early", "late", "noexit", "early", "late", "early"),
  levels = c("late", "early", "noexit")
)
r2_one_group <- factor(rep("solo", nrow(r2_sd)), levels = "solo")

r2_case("state_all_unweighted", r2_sd)
r2_case("state_all_weighted", r2_sdw)
r2_case("state_one_B", r2_sd, state = "B")
r2_case("state_subset_C_A", r2_sd, state = c("C", "A"))
r2_case("state_one_group", r2_sd, groups = r2_one_group, state = c("D", "A"))
r2_case("group_all_unweighted", r2_sd, groups = r2_groups, per.state = TRUE)
r2_case(
  "group_one_state_weighted", r2_sdw, groups = r2_groups,
  per.state = TRUE, state = "B"
)
r2_case(
  "group_subset_A_C", r2_sd, groups = r2_groups,
  per.state = TRUE, state = c("C", "A")
)
r2_case("missing_false", r2_sd, state = c("A", "B", "C"), with.missing = FALSE)
r2_case("missing_true_r_route", r2_sd, state = c("A", "B", "C"), with.missing = TRUE)

r2_errors <- c("case=message")
r2_error <- function(case, expr) {
  message <- tryCatch(
    { suppressMessages(expr); "NO_ERROR" },
    error = function(e) conditionMessage(e)
  )
  r2_errors <<- c(r2_errors, paste0(case, "=", gsub("[\r\n]", " ", message)))
}
r2_error("multi_group_state_route", seqsurv(
  r2_sd, groups = r2_groups, per.state = FALSE
))
r2_error("unobserved_state", seqsurv(r2_sd, state = c("U", "V")))
r2_error("not_sequence", seqsurv(R2))

csv_quote <- function(value) paste0('"', gsub('"', '""', value, fixed = TRUE), '"')
r2_lines <- c(
  "case,section,curve,row,time,n_risk,n_event,surv,std_err,lower,upper,value,text",
  vapply(r2_rows, function(row) paste(vapply(row, csv_quote, character(1)), collapse = ","),
         character(1))
)
writeLines(r2_lines, file.path(out_dir, "seqsurv_r2_cases.csv"), useBytes = TRUE)
writeLines(r2_errors, file.path(out_dir, "seqsurv_r2_errors.txt"), useBytes = TRUE)
writeLines(r2_messages, file.path(out_dir, "seqsurv_r2_messages.txt"), useBytes = TRUE)
cat(
  "SEQSURV_R2_GOLDENS_DONE cases=", length(unique(vapply(r2_rows, `[[`, "", "case"))),
  " rows=", length(r2_rows), " errors=", length(r2_errors) - 1L,
  " messages=", length(r2_messages) - 1L, "\n", sep = ""
)

## R4: seqsplot dispatcher selection, layout, and condition contract ---------
r4_sep <- "\036"
r4_collapse <- function(values) paste(as.character(values), collapse = r4_sep)
r4_observed <- function(seqdata) {
  alphabet(seqdata)[alphabet(seqdata) %in% seqstatl(seqdata)]
}
r4_groups <- function(group, n) {
  if (is.null(group)) factor(rep("Unique group", n)) else TraMineR:::group(group)
}
r4_labels <- function(seqdata, states) {
  attr(seqdata, "labels")[match(states, alphabet(seqdata))]
}

r4_panels <- list()
r4_panel_index <- 1L
r4_record_panels <- function(case, seqdata, group = NULL, per.state = FALSE,
                             which.states = NULL, main = "auto") {
  observed <- r4_observed(seqdata)
  selected <- if (is.null(which.states)) observed else
    alphabet(seqdata)[alphabet(seqdata) %in% which.states & alphabet(seqdata) %in% observed]
  grouped <- r4_groups(group, nrow(seqdata))
  group.levels <- levels(grouped)
  if (per.state) {
    panel.keys <- selected
    panel.labels <- r4_labels(seqdata, selected)
    panel.rows <- rep(list(seq_len(nrow(seqdata))), length(selected))
    curve.keys <- rep(list(group.levels), length(selected))
    legend.labels <- group.levels
  } else {
    panel.keys <- group.levels
    panel.labels <- group.levels
    panel.rows <- lapply(group.levels, function(level) which(grouped == level))
    curve.keys <- lapply(panel.rows, function(rows) {
      local.observed <- r4_observed(seqdata[rows, ])
      selected[selected %in% local.observed]
    })
    legend.labels <- r4_labels(seqdata, selected)
  }
  titles <- if (is.null(main)) rep("", length(panel.keys)) else if (main[[1]] == "auto") {
    panel.labels
  } else if (length(main) == 1L) {
    paste(main, "-", panel.labels)
  } else {
    main
  }
  for (index in seq_along(panel.keys)) {
    r4_panels[[r4_panel_index]] <<- data.frame(
      case = case,
      per_state = per.state,
      panel_index_zero = index - 1L,
      panel_key = panel.keys[[index]],
      title = titles[[index]],
      row_indexes_one = r4_collapse(panel.rows[[index]]),
      curve_keys = r4_collapse(curve.keys[[index]]),
      legend_labels = r4_collapse(legend.labels),
      stringsAsFactors = FALSE
    )
    r4_panel_index <<- r4_panel_index + 1L
  }
}

r4_record_panels("default", r2_sd)
r4_record_panels("grouped", r2_sd, group = r2_groups)
r4_record_panels("per_state", r2_sd, group = r2_groups, per.state = TRUE)
r4_record_panels(
  "subset_declared_order", r2_sd, group = r2_groups,
  which.states = c("C", "A")
)
r4_record_panels(
  "per_state_custom_title", r2_sd, group = r2_groups, per.state = TRUE,
  which.states = c("C", "A"), main = "Survival"
)
write.csv(
  do.call(rbind, r4_panels),
  file.path(out_dir, "seqsplot_r4_panels.csv"),
  row.names = FALSE,
  quote = TRUE,
  na = "NA"
)

r4_layouts <- list()
for (nplot in c(1L, 2L, 3L, 4L)) {
  cols <- min(nplot, 2L)
  rows <- ceiling(nplot / cols)
  layout <- TraMineR:::TraMineR.setlayout(
    nplot, rows, cols, "auto", "bottom", NA
  )
  left <- unique(layout$laymat[, 1])
  left <- left[left <= nplot]
  right_small <- TraMineR:::TraMineR.setlayout(
    nplot, rows, cols, "right", "all", 0.2
  )
  right_large <- TraMineR:::TraMineR.setlayout(
    nplot, rows, cols, "right", "all", 0.4
  )
  r4_layouts[[length(r4_layouts) + 1L]] <- data.frame(
    nplot = nplot,
    rows = rows,
    cols = cols,
    bottom_panels_one = r4_collapse(layout$axisp[layout$axisp > 0]),
    left_panels_one = r4_collapse(left),
    legend_position = if (is.null(layout$legpos)) "NULL" else layout$legpos,
    right_0_2_plot_fraction = sum(
      right_small$widths[-length(right_small$widths)]
    ),
    right_0_4_plot_fraction = sum(
      right_large$widths[-length(right_large$widths)]
    ),
    stringsAsFactors = FALSE
  )
}
write.csv(
  do.call(rbind, r4_layouts),
  file.path(out_dir, "seqsplot_r4_layout.csv"),
  row.names = FALSE,
  quote = TRUE
)

r4_conditions <- list()
r4_condition_index <- 1L
r4_capture <- function(case, expression, message_if_none = "") {
  messages <- character()
  outcome <- tryCatch(
    withCallingHandlers(
      {
        suppressWarnings(force(expression))
        "OK"
      },
      message = function(condition) {
        messages <<- c(messages, trimws(conditionMessage(condition)))
        invokeRestart("muffleMessage")
      }
    ),
    error = function(condition) paste0("ERROR: ", conditionMessage(condition))
  )
  if (length(messages) == 0L) messages <- message_if_none
  for (message in messages) {
    r4_conditions[[r4_condition_index]] <<- data.frame(
      case = case,
      outcome = outcome,
      message = gsub("[\r\n]", " ", message),
      stringsAsFactors = FALSE
    )
    r4_condition_index <<- r4_condition_index + 1L
  }
}

r4_pdf <- tempfile(fileext = ".pdf")
pdf(r4_pdf)
r4_capture(
  "partial_invalid_state",
  seqsplot(r2_sd, which.states = c("C", "U"), with.legend = FALSE)
)
r4_capture(
  "all_invalid_state",
  seqsplot(r2_sd, which.states = c("U", "V"), with.legend = FALSE)
)
r4_capture(
  "ignored_state",
  seqsplot(r2_sd, state = "B", with.legend = FALSE)
)
r4_capture(
  "group_length",
  seqsplot(r2_sd, group = c("a", "b"), with.legend = FALSE)
)
r4_capture(
  "ylab_length",
  seqsplot(r2_sd, group = r2_groups, ylab = c("a", "b"), with.legend = FALSE)
)
r4_time_ticks <- function(xtstep, tick.last) {
  fit <- suppressMessages(seqsurv(r2_sd))
  maximum <- max(sort(unique(fit$time)))
  ticks <- seq(1, maximum, xtstep)
  if (tick.last && ticks[length(ticks)] < maximum) ticks <- c(ticks, maximum)
  r4_collapse(ticks)
}
r4_capture(
  "xtstep_override",
  seqsplot(r2_sd, xtstep = 3, tick.last = FALSE, with.legend = FALSE),
  r4_time_ticks(3, FALSE)
)
r4_capture(
  "tick_last_override",
  seqsplot(r2_sd, xtstep = 4, tick.last = TRUE, with.legend = FALSE),
  r4_time_ticks(4, TRUE)
)
r4_capture(
  "use_layout_true",
  seqsplot(r2_sd, use.layout = TRUE, with.legend = FALSE)
)
r4_capture(
  "use_layout_false",
  seqsplot(r2_sd, use.layout = FALSE, with.legend = FALSE)
)
dev.off()
unlink(r4_pdf)
write.csv(
  do.call(rbind, r4_conditions),
  file.path(out_dir, "seqsplot_r4_conditions.csv"),
  row.names = FALSE,
  quote = TRUE
)
cat(
  "SEQSPLOT_R4_GOLDENS_DONE panels=", length(r4_panels),
  " layouts=", length(r4_layouts), " conditions=", length(r4_conditions),
  "\n", sep = ""
)
