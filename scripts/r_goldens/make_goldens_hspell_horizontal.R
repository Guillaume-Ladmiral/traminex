# Golden generation for TraMineRextras::HSPELL_to_STS horizontal overlap
# ownership and declared spell-column order.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
suppressMessages(library(TraMineRextras))
stopifnot(
  as.character(packageVersion("TraMineR")) == "2.2.14",
  as.character(packageVersion("TraMineRextras")) == "0.6.9"
)

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir, "hspell_horizontal_cases.txt")

rows <- character()
emit_case <- function(case, data, begin, end, status, overwrite) {
  result <- TraMineRextras::HSPELL_to_STS(
    data,
    begin = begin,
    end = end,
    status = status,
    overwrite = overwrite,
    id = "id"
  )
  spell_order <- paste(sprintf("%s:%s:%s", begin, end, status), collapse = ">")
  rows <<- c(
    rows,
    paste(
      case,
      row.names(result)[1],
      spell_order,
      toupper(as.character(overwrite)),
      paste(names(result), collapse = ","),
      paste(as.character(unlist(result[1, ], use.names = FALSE)), collapse = ","),
      sep = "\t"
    )
  )
}

emit_case(
  "right_tail",
  data.frame(id = "right", b1 = 1, e1 = 3, s1 = "A", b2 = 2, e2 = 4, s2 = "B"),
  c("b1", "b2"), c("e1", "e2"), c("s1", "s2"), TRUE
)
emit_case(
  "right_tail_no_overwrite",
  data.frame(id = "first", b1 = 1, e1 = 3, s1 = "A", b2 = 2, e2 = 4, s2 = "B"),
  c("b1", "b2"), c("e1", "e2"), c("s1", "s2"), FALSE
)
emit_case(
  "out_of_chronological_order",
  data.frame(id = "ordered", b1 = 3, e1 = 4, s1 = "B", b2 = 1, e2 = 3, s2 = "A"),
  c("b1", "b2"), c("e1", "e2"), c("s1", "s2"), TRUE
)
emit_case(
  "nested",
  data.frame(id = "nested", b1 = 1, e1 = 6, s1 = "A", b2 = 2, e2 = 4, s2 = "B"),
  c("b1", "b2"), c("e1", "e2"), c("s1", "s2"), TRUE
)
emit_case(
  "chained",
  data.frame(
    id = "chain",
    b1 = 1, e1 = 5, s1 = "A",
    b2 = 2, e2 = 3, s2 = "B",
    b3 = 3, e3 = 6, s3 = "C"
  ),
  c("b1", "b2", "b3"), c("e1", "e2", "e3"), c("s1", "s2", "s3"), TRUE
)
emit_case(
  "same_bounds",
  data.frame(id = "tie", b1 = 1, e1 = 4, s1 = "A", b2 = 1, e2 = 4, s2 = "B"),
  c("b1", "b2"), c("e1", "e2"), c("s1", "s2"), TRUE
)
emit_case(
  "adjacent",
  data.frame(id = "adjacent", b1 = 1, e1 = 2, s1 = "A", b2 = 3, e2 = 4, s2 = "B"),
  c("b1", "b2"), c("e1", "e2"), c("s1", "s2"), TRUE
)

writeLines(
  c(
    "# R 4.6.2; TraMineR 2.2.14; TraMineRextras 0.6.9",
    "# case\tid\tspell_order\toverwrite\tR_columns\tcells",
    rows
  ),
  out_file,
  useBytes = TRUE
)

pvar_data <- data.frame(
  id = c("case-z", "case-a"),
  begin1 = c(2000, 2000),
  end1 = c(2002, 2002),
  status1 = c("A", "B"),
  origin = c(1998, 2000),
  stringsAsFactors = FALSE
)
pvar_result <- TraMineRextras::HSPELL_to_STS(
  pvar_data,
  begin = "begin1",
  end = "end1",
  status = "status1",
  pvar = "origin",
  fillblanks = "_",
  id = "id"
)
pvar_rows <- character()
for (i in seq_len(nrow(pvar_result))) {
  input_row <- match(row.names(pvar_result)[i], pvar_data$id)
  pvar_rows <- c(
    pvar_rows,
    paste(
      row.names(pvar_result)[i],
      pvar_data$origin[input_row],
      pvar_data$begin1[input_row],
      pvar_data$end1[input_row],
      pvar_data$status1[input_row],
      paste(names(pvar_result), collapse = ","),
      paste(as.character(unlist(pvar_result[i, ], use.names = FALSE)), collapse = ","),
      sep = "\t"
    )
  )
}
writeLines(
  c(
    "# R 4.6.2; TraMineR 2.2.14; TraMineRextras 0.6.9",
    "# id\tpvar\tbegin\tend\tstatus\tR_columns\tcells",
    pvar_rows
  ),
  file.path(out_dir, "hspell_pvar_case.txt"),
  useBytes = TRUE
)

# Whole-matrix cases (2026-09-22 audit S2-01..S2-04): row order, missing ends,
# the inferred time grid, blank cells, endObs, fixed.status, numeric statuses,
# overwrite = FALSE and explicit bounds.  Every row of every result is written
# in R's own row order, with R's column names; NA cells are written "NA".
matrix_rows <- character()
emit_matrix <- function(case, result) {
  for (i in seq_len(nrow(result))) {
    matrix_rows <<- c(
      matrix_rows,
      paste(
        case,
        i,
        row.names(result)[i],
        paste(names(result), collapse = ","),
        paste(ifelse(is.na(unlist(result[i, ], use.names = FALSE)), "NA",
          as.character(unlist(result[i, ], use.names = FALSE))), collapse = ","),
        sep = "\t"
      )
    )
  }
}
two_spells <- data.frame(
  id = c("p1", "p2"), begin1 = c(1, 1), end1 = c(3, 2), status1 = c("A", "A"),
  begin2 = c(2, 3), end2 = c(4, 4), status2 = c("B", "B"), stringsAsFactors = FALSE
)
bs <- c("begin1", "begin2"); es <- c("end1", "end2"); ss <- c("status1", "status2")
emit_matrix("fillblanks", suppressMessages(TraMineRextras::HSPELL_to_STS(
  two_spells, begin = bs, end = es, status = ss, fillblanks = "_", id = "id")))
emit_matrix("no_id", suppressMessages(TraMineRextras::HSPELL_to_STS(
  two_spells, begin = bs, end = es, status = ss)))
emit_matrix("overwrite_false", suppressMessages(TraMineRextras::HSPELL_to_STS(
  two_spells, begin = bs, end = es, status = ss, overwrite = FALSE, id = "id")))
emit_matrix("explicit_bounds", suppressMessages(TraMineRextras::HSPELL_to_STS(
  two_spells, begin = bs, end = es, status = ss, tmin = 0, tmax = 6, id = "id")))
emit_matrix("endobs_scalar", suppressMessages(TraMineRextras::HSPELL_to_STS(
  two_spells, begin = bs, end = es, status = ss, endObs = 3, id = "id")))
endobs_frame <- two_spells
endobs_frame$obs <- c(3, 4)
emit_matrix("endobs_column", suppressMessages(TraMineRextras::HSPELL_to_STS(
  endobs_frame, begin = bs, end = es, status = ss, fillblanks = "_", endObs = "obs", id = "id")))
twelve <- data.frame(id = 1:12, b = 1, e = 3, s = LETTERS[1:12], stringsAsFactors = FALSE)
emit_matrix("row_order", suppressMessages(TraMineRextras::HSPELL_to_STS(
  twelve, begin = "b", end = "e", status = "s", id = "id")))
na_ends <- data.frame(
  id = 1:3, b1 = c(1, 2, 1), e1 = c(2, NA, 3), s1 = c("A", "B", "A"),
  b2 = c(3, NA, 4), e2 = c(NA, NA, 5), s2 = c("C", "D", "C"), stringsAsFactors = FALSE
)
emit_matrix("missing_ends", suppressMessages(TraMineRextras::HSPELL_to_STS(
  na_ends, begin = c("b1", "b2"), end = c("e1", "e2"), status = c("s1", "s2"), id = "id")))
emit_matrix("fixed_status_vector", suppressMessages(TraMineRextras::HSPELL_to_STS(
  na_ends, begin = c("b1", "b2"), end = c("e1", "e2"), fixed.status = c("X", "Y"), id = "id")))
emit_matrix("fixed_status_scalar", suppressMessages(TraMineRextras::HSPELL_to_STS(
  na_ends, begin = c("b1", "b2"), end = c("e1", "e2"), fixed.status = "X", id = "id")))
numeric_status <- data.frame(
  id = 1:2, b1 = c(0, 1), e1 = c(2, 2), s1 = c(1, 2),
  b2 = c(3, 3), e2 = c(5, 5), s2 = c(2.5, NA)
)
emit_matrix("numeric_status", suppressMessages(TraMineRextras::HSPELL_to_STS(
  numeric_status, begin = c("b1", "b2"), end = c("e1", "e2"), status = c("s1", "s2"), id = "id")))
# R sets a missing end to the *unshifted* tmax (HSPELL_to_STS.R:37) after
# shifting every time by 1 - tmin, so with tmin < 1 an open spell stops
# 1 - tmin periods early.  traminex deliberately ends it at tmax; this case
# records what R actually does.
open_from_zero <- data.frame(id = 1:2, b1 = c(0, 1), e1 = c(3, NA), s1 = c("A", "B"))
emit_matrix("open_end_tmin_zero", suppressMessages(TraMineRextras::HSPELL_to_STS(
  open_from_zero, begin = "b1", end = "e1", status = "s1", id = "id")))
help_example <- data.frame(
  begin1 = rep(1, 5), end1 = c(2:5, NA), status1 = 1:5,
  begin2 = c(3:6, NA), end2 = rep(NA, 5), status2 = 5:1
)
emit_matrix("help_example", suppressMessages(TraMineRextras::HSPELL_to_STS(
  help_example, begin = c("begin1", "begin2"), end = c("end1", "end2"),
  status = c("status1", "status2"))))
writeLines(
  c(
    "# R 4.6.2; TraMineR 2.2.14; TraMineRextras 0.6.9",
    "# case\trow\trow_name\tR_columns\tcells",
    matrix_rows
  ),
  file.path(out_dir, "hspell_horizontal_matrix_cases.txt"),
  useBytes = TRUE
)

cat("HSPELL_HORIZONTAL_GOLDENS_DONE\n")
