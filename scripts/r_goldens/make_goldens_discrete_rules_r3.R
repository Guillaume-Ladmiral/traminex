# Golden generation for TraMineRextras 0.6.9 createdatadiscrete.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineRextras))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
stopifnot(as.character(packageVersion("TraMineRextras")) == "0.6.9")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ids <- c(101, 102, 103, 104)
data <- data.frame(
  IDPERS = ids,
  A = c(20, 22, 19, 30),
  AST = c(1, 0, 1, 0),
  B = c(23, 21, 24, 20),
  BST = c(1, 1, 0, 1),
  C = c(25, 24, 22, 21),
  CST = c(0, 1, 1, 0),
  sex = c(1, 2, 1, 2),
  cohort = c(7, 7, 8, 8)
)
agemin <- data.frame(IDPERS = ids, AGE = c(18, 19, 18, 20))
agemax <- data.frame(IDPERS = ids, AGE = c(25, 24, 24, 26))

frac_ids <- c(201, 202, 203)
frac_data <- data.frame(
  IDPERS = frac_ids,
  A = c(20.5, 22, 19.75),
  AST = c(1, 0, 1),
  B = c(23.25, 21, 30),
  BST = c(1, 1, 0)
)
frac_min <- data.frame(IDPERS = frac_ids, AGE = c(18, 19.5, 18))
frac_max <- data.frame(IDPERS = frac_ids, AGE = c(25, 24, 24.6))

captured_messages <- character()
capture_messages <- function(expression) {
  withCallingHandlers(
    force(expression),
    message = function(condition) {
      captured_messages <<- c(captured_messages, trimws(conditionMessage(condition)))
      invokeRestart("muffleMessage")
    }
  )
}

cases <- list(
  no_supvar = capture_messages(
    createdatadiscrete(ids, data, c("A", "B", "C"), agemin, agemax)
  ),
  one_supvar = capture_messages(
    createdatadiscrete(ids, data, c("A", "B", "C"), agemin, agemax, "sex")
  ),
  # Fractional ages (TX-AUD-20261019-015): R counts the rows with
  # as.integer((agemax - agemin) + 1), takes the ages agemin:agemax and sets an
  # independent event from row as.integer((age - agemin) + 1) (createdatadiscrete.R:44-72).
  fractional_ages = capture_messages(
    createdatadiscrete(frac_ids, frac_data, c("A", "B"), frac_min, frac_max)
  )
)

records <- list()
record_index <- 1L
record <- function(case_name, event, record_type, row, column, value) {
  records[[record_index]] <<- data.frame(
    case = case_name,
    event = event,
    record_type = record_type,
    row = as.character(row),
    column = column,
    value = as.character(value),
    stringsAsFactors = FALSE
  )
  record_index <<- record_index + 1L
}

for (case_name in names(cases)) {
  result <- cases[[case_name]]
  record(case_name, "", "mapping", "", "keys", paste(names(result), collapse = ","))
  record(case_name, "", "mapping", "", "vars_length", length(result$vars))
  record(case_name, "", "mapping", "", "ids", paste(result$ids, collapse = ","))
  for (event in names(result$data)) {
    frame <- result$data[[event]]
    record(
      case_name,
      event,
      "columns",
      "",
      "order",
      paste(names(frame), collapse = ",")
    )
    for (column in names(frame)) {
      record(case_name, event, "typeof", "", column, typeof(frame[[column]]))
    }
    for (row in seq_len(nrow(frame))) {
      for (column in names(frame)) {
        record(case_name, event, "cell", row, column, frame[row, column])
      }
    }
  }
}

write.csv(
  do.call(rbind, records),
  file.path(out_dir, "discrete_rules_r3_cases.csv"),
  row.names = FALSE,
  quote = TRUE
)

condition_lines <- vapply(
  captured_messages,
  function(message) paste("base", "message", message, sep = "|"),
  character(1)
)
capture_condition <- function(case_name, expression) {
  succeeded <- TRUE
  tryCatch(
    withCallingHandlers(
      force(expression),
      message = function(condition) {
        condition_lines <<- c(
          condition_lines,
          paste(case_name, "message", trimws(conditionMessage(condition)), sep = "|")
        )
        invokeRestart("muffleMessage")
      },
      warning = function(condition) {
        condition_lines <<- c(
          condition_lines,
          paste(case_name, "warning", conditionMessage(condition), sep = "|")
        )
        invokeRestart("muffleWarning")
      }
    ),
    error = function(condition) {
      succeeded <<- FALSE
      condition_lines <<- c(
        condition_lines,
        paste(case_name, class(condition)[[1]], conditionMessage(condition), sep = "|")
      )
    }
  )
  if (succeeded) {
    condition_lines <<- c(condition_lines, paste(case_name, "success", sep = "|"))
  }
  invisible(NULL)
}

capture_condition(
  "multiple_supvar",
  createdatadiscrete(
    ids,
    data,
    c("A", "B", "C"),
    agemin,
    agemax,
    c("sex", "cohort")
  )
)
capture_condition(
  "single_variable",
  createdatadiscrete(ids, data, "A", agemin, agemax)
)
capture_condition(
  "empty_variables",
  createdatadiscrete(ids, data, character(), agemin, agemax)
)
capture_condition(
  "missing_window_id",
  createdatadiscrete(ids, data, c("A", "B"), agemin[-2, ], agemax)
)
capture_condition(
  "na_window",
  createdatadiscrete(
    ids,
    data,
    c("A", "B"),
    transform(agemin, AGE = replace(AGE, IDPERS == 102, NA)),
    agemax
  )
)
capture_condition(
  "missing_status_column",
  createdatadiscrete(ids, data[, names(data) != "BST"], c("A", "B"), agemin, agemax)
)

progress_ids <- 1:1000
progress_data <- data.frame(
  IDPERS = progress_ids,
  A = rep(2, 1000),
  AST = rep(0, 1000),
  B = rep(2, 1000),
  BST = rep(0, 1000)
)
progress_min <- data.frame(IDPERS = progress_ids, AGE = rep(1, 1000))
progress_max <- data.frame(IDPERS = progress_ids, AGE = rep(2, 1000))
capture_condition(
  "progress",
  createdatadiscrete(
    progress_ids,
    progress_data,
    c("A", "B"),
    progress_min,
    progress_max
  )
)

# An event before agemin (TX-AUD-20261019-015 task): every variable is the dependent one in
# turn, so R stops when its row count, as.integer((age - agemin) + 1), is negative.
capture_condition(
  "event_before_agemin",
  createdatadiscrete(c(301), data.frame(IDPERS = 301, A = 15, AST = 1, B = 21, BST = 1), c("A", "B"),
                     data.frame(IDPERS = 301, AGE = 18), data.frame(IDPERS = 301, AGE = 25))
)
# With B first, A's early event is first met as an independent variable: its row index
# as.integer((15 - 18) + 1) = -2 is negative and R stops.
capture_condition(
  "independent_event_before_agemin",
  createdatadiscrete(c(302), data.frame(IDPERS = 302, A = 15, AST = 1, B = 21, BST = 1), c("B", "A"),
                     data.frame(IDPERS = 302, AGE = 18), data.frame(IDPERS = 302, AGE = 25))
)

writeLines(
  condition_lines,
  file.path(out_dir, "discrete_rules_r3_errors.txt"),
  useBytes = TRUE
)

# seqerulesdisc oracle cases. This larger deterministic fixture avoids
# separation in the binomial models and includes one deliberately low-event
# dependent variable.
rule_ids <- 1:48
rule_data <- data.frame(
  IDPERS = rule_ids,
  A = 19 + rule_ids %% 4,
  AST = as.integer(rule_ids %% 4 != 0),
  B = 21 + rule_ids %% 4,
  BST = as.integer(rule_ids %% 5 != 0),
  C = 23 + rule_ids %% 3,
  CST = as.integer(rule_ids %% 3 != 0),
  D = rep(25, length(rule_ids)),
  DST = as.integer(rule_ids == 1),
  sex = 1 + rule_ids %% 2,
  cohort = 1 + rule_ids %% 3
)
rule_min <- data.frame(IDPERS = rule_ids, AGE = 18)
rule_max <- data.frame(IDPERS = rule_ids, AGE = 26)
rule_discrete <- suppressMessages(
  createdatadiscrete(
    rule_ids,
    rule_data,
    c("A", "B", "C", "D"),
    rule_min,
    rule_max,
    "sex"
  )
)

event_rows <- list()
event_index <- 1L
for (identifier in rule_ids) {
  for (event in c("A", "B", "C", "D")) {
    if (rule_data[identifier, paste0(event, "ST")] == 1) {
      event_rows[[event_index]] <- data.frame(
        id = identifier,
        timestamp = rule_data[identifier, event],
        event = event,
        stringsAsFactors = FALSE
      )
      event_index <- event_index + 1L
    }
  }
}
rule_tse <- do.call(rbind, event_rows)
rule_eseq <- suppressMessages(seqecreate(rule_tse))
rule_labels <- c(
  "(A)", "(B)", "(C)", "(D)", "(A)-(B)", "(B)-(C)",
  "(A)-(B)-(C)", "(A,B)-(C)", "(A)-(D)"
)
rule_fsub <- suppressMessages(seqefsub(rule_eseq, str.subseq = rule_labels))

rule_messages <- character()
rule_outputs <- character()
run_rule_case <- function(case_name, expression) {
  holder <- new.env(parent = emptyenv())
  output <- capture.output(
    holder$result <- withCallingHandlers(
      force(expression),
      message = function(condition) {
        rule_messages <<- c(
          rule_messages,
          paste(case_name, "message", trimws(conditionMessage(condition)), sep = "|")
        )
        invokeRestart("muffleMessage")
      }
    ),
    type = "output"
  )
  if (length(output)) {
    rule_outputs <<- c(
      rule_outputs,
      vapply(output, function(line) paste(case_name, "output", line, sep = "|"), character(1))
    )
  }
  holder$result
}

rule_cases <- list(
  cloglog_all = run_rule_case(
    "cloglog_all",
    seqerulesdisc(rule_fsub, rule_discrete, rule_data, pvalue = 1, adjust = FALSE)
  ),
  logit_all = run_rule_case(
    "logit_all",
    seqerulesdisc(
      rule_fsub, rule_discrete, rule_data, pvalue = 1,
      adjust = FALSE, link = "logit"
    )
  ),
  adjusted = run_rule_case(
    "adjusted",
    seqerulesdisc(rule_fsub, rule_discrete, rule_data, pvalue = 1, adjust = TRUE)
  ),
  topt = run_rule_case(
    "topt",
    seqerulesdisc(
      rule_fsub, rule_discrete, rule_data, pvalue = 1,
      adjust = FALSE, topt = TRUE
    )
  ),
  dep_c = run_rule_case(
    "dep_c",
    seqerulesdisc(
      rule_fsub, rule_discrete, rule_data, pvalue = 1,
      adjust = FALSE, dep = "C"
    )
  ),
  dep_unknown = run_rule_case(
    "dep_unknown",
    seqerulesdisc(
      rule_fsub, rule_discrete, rule_data, pvalue = 1,
      adjust = FALSE, dep = "Z"
    )
  ),
  duplicates = run_rule_case(
    "duplicates",
    seqerulesdisc(
      rule_fsub[c(5, 5, 6)], rule_discrete, rule_data,
      pvalue = 1, adjust = FALSE
    )
  ),
  low_event = run_rule_case(
    "low_event",
    seqerulesdisc(
      rule_fsub[9], rule_discrete, rule_data,
      pvalue = 1, adjust = FALSE
    )
  )
)

rule_records <- list()
rule_record_index <- 1L
record_rule <- function(case_name, record_type, row, key, value) {
  rule_records[[rule_record_index]] <<- data.frame(
    case = case_name,
    record_type = record_type,
    row = as.character(row),
    key = key,
    value = as.character(value),
    stringsAsFactors = FALSE
  )
  rule_record_index <<- rule_record_index + 1L
}

for (case_name in names(rule_cases)) {
  result <- rule_cases[[case_name]]
  graph <- attr(result, "graph")
  record_rule(case_name, "metadata", "", "columns", paste(names(result), collapse = "\036"))
  record_rule(case_name, "metadata", "", "nrow", nrow(result))
  record_rule(case_name, "metadata", "", "source_subseqs", length(rule_fsub$subseq))
  for (row in seq_len(nrow(result))) {
    record_rule(case_name, "frame", row, "Hazards", result[row, "Hazards"])
  }
  for (section in c("events", "hr", "hrindep", "pvalues")) {
    values <- graph[[section]]
    record_rule(case_name, "graph_length", "", section, length(values))
    if (length(values)) {
      for (row in seq_along(values)) {
        record_rule(
          case_name,
          paste0("graph_", section),
          row,
          section,
          paste(format(values[[row]], digits = 17), collapse = "\036")
        )
      }
    }
  }
}

# Intended safe supvars correction: ordinary adjustment covariates in the same
# discrete-time binomial GLM. Upstream seqerulesdisc's experimental Cox branch
# is captured below as a failure.
control_frame <- rule_discrete$data$B
control_frame$cohort <- rule_data$cohort[
  match(control_frame$IDPERS, rule_data$IDPERS)
]
for (control_link in c("cloglog", "logit")) {
  control_model <- glm(
    B ~ A + AGE + sex + cohort,
    data = control_frame,
    family = binomial(link = control_link)
  )
  coefficients <- summary(control_model)$coefficient
  for (term in rownames(coefficients)) {
    record_rule(
      paste0("corrected_supvars_", control_link),
      "coefficient",
      term,
      "estimate",
      format(coefficients[term, 1], digits = 17)
    )
    record_rule(
      paste0("corrected_supvars_", control_link),
      "coefficient",
      term,
      "pvalue",
      format(coefficients[term, 4], digits = 17)
    )
  }
}

write.csv(
  do.call(rbind, rule_records),
  file.path(out_dir, "discrete_seq_rules_r3_cases.csv"),
  row.names = FALSE,
  quote = TRUE
)

rule_conditions <- c(rule_messages, rule_outputs)
capture_rule_condition <- function(case_name, expression) {
  succeeded <- TRUE
  tryCatch(
    withCallingHandlers(
      force(expression),
      warning = function(condition) {
        rule_conditions <<- c(
          rule_conditions,
          paste(case_name, "warning", conditionMessage(condition), sep = "|")
        )
        invokeRestart("muffleWarning")
      }
    ),
    error = function(condition) {
      succeeded <<- FALSE
      rule_conditions <<- c(
        rule_conditions,
        paste(case_name, class(condition)[[1]], conditionMessage(condition), sep = "|")
      )
    }
  )
  if (succeeded) {
    rule_conditions <<- c(rule_conditions, paste(case_name, "success", sep = "|"))
  }
  invisible(NULL)
}

capture_rule_condition(
  "supvars",
  seqerulesdisc(
    rule_fsub, rule_discrete, rule_data,
    pvalue = 1, adjust = FALSE, supvars = "sex"
  )
)
capture_rule_condition(
  "bad_link",
  seqerulesdisc(
    rule_fsub, rule_discrete, rule_data,
    pvalue = 1, adjust = FALSE, link = "bad"
  )
)
capture_rule_condition(
  "empty_subseqs",
  seqerulesdisc(
    rule_fsub[integer(0)], rule_discrete, rule_data,
    pvalue = 1, adjust = FALSE
  )
)

writeLines(
  rule_conditions,
  file.path(out_dir, "discrete_seq_rules_r3_errors.txt"),
  useBytes = TRUE
)

cat("DISCRETE_RULES_R3_GOLDENS_DONE\n")
