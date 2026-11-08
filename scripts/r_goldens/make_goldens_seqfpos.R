# Golden generation: R2 seqfpos and seqfposend single-state contracts.
#
# Pins raw rectangular first positions plus first-spell end positions to
# TraMineR 2.2-14, including missing, void-padded, empty, absent, and invalid
# state routes.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

raw <- rbind(
  c("A", "A", "B", "C", "C"),
  c("B", "A", "A", "C", "B"),
  c("B", "B", "C", "C", "C"),
  c("*", "A", "*", "B", "B"),
  c("A", "B", "*", "*", "*"),
  c("*", "*", "*", "*", "*")
)
sd <- suppressMessages(seqdef(
  raw,
  alphabet = c("A", "B", "C", "D"),
  missing = "*"
))

rows <- list()
append_values <- function(case, function_name, state, value) {
  rows[[length(rows) + 1L]] <<- data.frame(
    case = case,
    `function` = function_name,
    state = state,
    row = seq_along(value),
    value = as.vector(value),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

states <- c("A", "B", "C", "*", "D", "Z")
for (state in states) {
  append_values(
    paste0("raw_", ifelse(state == "*", "missing", state)),
    "seqfpos",
    state,
    seqfpos(sd, state)
  )
}

for (wm in c(FALSE, TRUE)) {
  for (lead in c(-1L, 0L, 1L, 2L)) {
    for (origin in c(FALSE, TRUE)) {
      for (state in states) {
        state_tag <- ifelse(state == "*", "missing", state)
        case <- paste0(
          "spells_", state_tag,
          "_wm", wm,
          "_lead", lead,
          "_origin", origin
        )
        append_values(
          case,
          "seqfposend",
          state,
          seqfposend(
            sd,
            state,
            with.missing = wm,
            lead = lead,
            from.seq.start = origin
          )
        )
      }
    }
  }
}

utils::write.csv(
  do.call(rbind, rows),
  file.path(out_dir, "seqfpos_r2_cases.csv"),
  row.names = FALSE,
  na = "NA",
  quote = FALSE
)

normalize_error <- function(expr) {
  tryCatch(
    {
      force(expr)
      "VALUE"
    },
    error = function(e) paste0("ERROR:", conditionMessage(e))
  )
}
writeLines(
  c(
    paste0("seqfpos_non_stslist=", normalize_error(seqfpos(raw, "A"))),
    paste0(
      "seqfposend_non_stslist=",
      normalize_error(seqfposend(raw, "A"))
    ),
    paste0("seqfpos_invalid_state=", paste(seqfpos(sd, "Z"), collapse = ",")),
    paste0(
      "seqfposend_invalid_state=",
      paste(seqfposend(sd, "Z"), collapse = ",")
    )
  ),
  file.path(out_dir, "seqfpos_r2_errors.txt")
)

cat("SEQFPOS_R2_GOLDENS_DONE\n")
