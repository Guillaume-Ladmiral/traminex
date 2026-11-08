# Live-R goldens for seqdiff (RP4 T16).
#
# Sliding-window dissassoc (Pseudo F/Fbf/R2/Bartlett/Levene + per-group
# discrepancy) on a deterministic 10x6 state-sequence fixture (alphabet 1:4),
# binary group, and multiple public seqdist.args configurations.
#
# Requires R 4.6.2 + TraMineR 2.2-14.

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
source("scripts/r_goldens/make_goldens_phase2_inputs.R")  # write_phase2_csv()
suppressMessages(library(TraMineR))

if (!identical(as.character(packageVersion("TraMineR")), "2.2.14")) {
  stop("seqdiff goldens require TraMineR 2.2-14")
}

set.seed(424242)
n <- 10
L <- 6
M <- matrix(sample(1:4, n * L, replace = TRUE), nrow = n, ncol = L)
group <- rep(c(1, 2), length.out = n)
w <- 0.5 + 0.1 * (1:n)
sd <- suppressMessages(TraMineR::seqdef(M, alphabet = 1:4, weights = w))

write_phase2_csv("seqdiff_input_codes", M - 1L, "# seqdiff input: 0-based state codes (alphabet 1:4)\n")
write_phase2_csv("seqdiff_input_group", matrix(group, nrow = 1), "# seqdiff input: binary group (1,2)\n")
write_phase2_csv("seqdiff_input_weights", matrix(w, nrow = 1), "# seqdiff input: weights 0.5+0.1*i\n")

gen <- function(tag, weighted) {
  r <- suppressMessages(TraMineR::seqdiff(
    sd, group, seqdist.args = list(method = "LCS", norm = "none"), weighted = weighted
  ))
  write_phase2_csv(
    paste0("seqdiff_", tag, "_stat"), r$stat,
    "# seqdiff stat: PseudoF, PseudoFbf, PseudoR2, Bartlett, Levene (rows=positions)\n"
  )
  write_phase2_csv(
    paste0("seqdiff_", tag, "_disc"), r$discrepancy,
    "# seqdiff discrepancy: group levels + Total (rows=positions)\n"
  )
}

gen("uw", FALSE)
gen("w", TRUE)

public_rows <- list()

append_result <- function(rows, case, result) {
  append_matrix <- function(rows, kind, value) {
    for (row_index in seq_len(nrow(value))) {
      for (column_index in seq_len(ncol(value))) {
        rows[[length(rows) + 1L]] <- data.frame(
          case = case,
          kind = kind,
          row_position = row_index,
          r_row_label = rownames(value)[row_index],
          column_position = column_index,
          column = colnames(value)[column_index],
          value = as.numeric(value[row_index, column_index]),
          stringsAsFactors = FALSE
        )
      }
    }
    rows
  }
  rows <- append_matrix(rows, "stat", result$stat)
  append_matrix(rows, "discrepancy", result$discrepancy)
}

capture_messages <- function(expression) {
  messages <- character()
  value <- withCallingHandlers(
    eval.parent(substitute(expression)),
    message = function(condition) {
      messages <<- c(messages, conditionMessage(condition))
      invokeRestart("muffleMessage")
    }
  )
  list(value = value, messages = messages)
}

public_rows <- append_result(
  public_rows,
  "default_lcs_auto",
  suppressMessages(TraMineR::seqdiff(sd, group))
)
public_rows <- append_result(
  public_rows,
  "lcs_none_map",
  suppressMessages(TraMineR::seqdiff(
    sd,
    group,
    seqdist.args = list(method = "LCS", norm = "none"),
    weighted = FALSE
  ))
)

sm <- matrix(
  c(
    0, 0.5, 1.5, 2,
    0.5, 0, 2, 1,
    1.5, 2, 0, 0.75,
    2, 1, 0.75, 0
  ),
  nrow = 4,
  byrow = TRUE,
  dimnames = list(as.character(1:4), as.character(1:4))
)
public_rows <- append_result(
  public_rows,
  "om_custom",
  suppressWarnings(suppressMessages(TraMineR::seqdiff(
    sd,
    group,
    seqdist.args = list(method = "OM", norm = "none", sm = sm, indel = 0.75),
    weighted = FALSE
  )))
)
public_rows <- append_result(
  public_rows,
  "chi2_weighted_distance",
  suppressMessages(TraMineR::seqdiff(
    sd,
    group,
    seqdist.args = list(method = "CHI2", step = 2, overlap = TRUE, weighted = TRUE),
    weighted = FALSE
  ))
)
public_rows <- append_result(
  public_rows,
  "chi2_unweighted_distance",
  suppressMessages(TraMineR::seqdiff(
    sd,
    group,
    seqdist.args = list(method = "CHI2", step = 2, overlap = TRUE, weighted = FALSE),
    weighted = FALSE
  ))
)

group_missing <- group
group_missing[c(3, 8)] <- NA
public_rows <- append_result(
  public_rows,
  "ordinary_group_missing",
  suppressMessages(TraMineR::seqdiff(sd, group_missing))
)
group_locale <- rep(c("a", "B"), length.out = n)
public_rows <- append_result(
  public_rows,
  "ordinary_group_locale",
  suppressMessages(TraMineR::seqdiff(sd, group_locale))
)

group_matrix <- rbind(
  c("z", "*", "m", "z", "a", "m"),
  c("a", "a", "m", "*", "*", "*"),
  c("m", "z", "z", "a", "a", "m"),
  c("z", "m", "a", "a", "m", "z"),
  c("a", "z", "m", "z", "a", "m"),
  c("m", "a", "z", "m", "z", "a"),
  c("z", "m", "a", "z", "m", "a"),
  c("a", "z", "m", "a", "z", "m"),
  c("m", "a", "z", "m", "a", "z"),
  c("z", "a", "m", "z", "m", "a")
)
group_sts <- suppressMessages(TraMineR::seqdef(
  group_matrix,
  alphabet = c("z", "a", "m"),
  missing = "*"
))
public_rows <- append_result(
  public_rows,
  "sts_group",
  suppressMessages(TraMineR::seqdiff(
    sd,
    group_sts,
    seqdist.args = list(method = "LCS", norm = "none")
  ))
)

ragged_matrix <- rbind(
  c("1", "2", "3", "4", "1", "2"),
  c("2", "1", "*", "3", "4", "1"),
  c("3", "3", "1", "2", "2", "3"),
  c("4", "1", "2", "*", "2", "1"),
  c("1", "1", "2", "3", "4", "2"),
  c("2", "3", "4", "1", "1", "3"),
  c("3", "2", "1", "4", "3", "*"),
  c("4", "4", "1", "1", "1", "*"),
  c("1", "3", "2", "4", "3", "*"),
  c("2", "2", "3", "3", "4", "*"),
  c("3", "1", "2", "2", "1", "*"),
  c("4", "2", "1", "3", "2", "*")
)
ragged_group <- rep(c(1, 2), length.out = nrow(ragged_matrix))
ragged_weights <- 0.5 + 0.1 * (1:nrow(ragged_matrix))
ragged_sd <- suppressMessages(TraMineR::seqdef(
  ragged_matrix,
  alphabet = 1:4,
  missing = "*",
  weights = ragged_weights
))
public_rows <- append_result(
  public_rows,
  "missing_override",
  suppressMessages(TraMineR::seqdiff(
    ragged_sd,
    ragged_group,
    seqdist.args = list(
      method = "LCS",
      norm = "none",
      with.missing = FALSE
    ),
    with.missing = TRUE
  ))
)

deprecated <- capture_messages(TraMineR::seqdiff(
  sd,
  group,
  seqdist_arg = list(method = "LCS", norm = "none"),
  weighted = FALSE
))
public_rows <- append_result(public_rows, "deprecated_alias", deprecated$value)

public_table <- do.call(rbind, public_rows)
write.csv(
  public_table,
  file = "tests/r_goldens/seqdiff_public_cases.csv",
  row.names = FALSE,
  na = "NA"
)
writeLines(
  c(
    "case|message",
    paste(
      "deprecated_alias",
      trimws(gsub("[\r\n]+", " ", deprecated$messages)),
      sep = "|"
    )
  ),
  con = "tests/r_goldens/seqdiff_public_warnings.txt",
  sep = "\n"
)

capture_error <- function(case, expression) {
  condition <- tryCatch(
    {
      suppressMessages(eval.parent(substitute(expression)))
      NULL
    },
    error = function(condition) condition
  )
  if (is.null(condition)) stop("expected R error for case: ", case)
  paste(
    case,
    paste(class(condition), collapse = ","),
    gsub("[\r\n]+", " ", conditionMessage(condition)),
    sep = "|"
  )
}
both_map_error <- tryCatch(
  {
    suppressMessages(TraMineR::seqdiff(
      sd,
      group,
      seqdist.args = list(method = "LCS"),
      seqdist_arg = list(method = "LCS")
    ))
    NULL
  },
  error = function(condition) condition
)
if (is.null(both_map_error)) stop("expected both-map conflict")
writeLines(
  c(
    "case|class|message",
    paste(
      "both_maps",
      paste(class(both_map_error), collapse = ","),
      gsub("[\r\n]+", " ", conditionMessage(both_map_error)),
      sep = "|"
    ),
    capture_error(
      "empty_map",
      TraMineR::seqdiff(sd, group, seqdist.args = list())
    ),
    capture_error(
      "cmprange_reversed",
      TraMineR::seqdiff(
        sd,
        group,
        cmprange = c(1, 0),
        seqdist.args = list(method = "LCS", norm = "none")
      )
    ),
    capture_error(
      "cmprange_outside",
      TraMineR::seqdiff(
        sd,
        group,
        cmprange = c(0, ncol(sd) + 1),
        seqdist.args = list(method = "LCS", norm = "none")
      )
    )
  ),
  con = "tests/r_goldens/seqdiff_public_errors.txt",
  sep = "\n"
)

cat("SEQDIFF_GOLDENS_DONE\n")
