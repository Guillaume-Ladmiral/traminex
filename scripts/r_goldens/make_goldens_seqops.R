# Golden generation: sequence-operation indicators (RP6a BUILD) vs TraMineR 2.2-14.
#
#   seqconc   per-sequence string of non-void states joined by sep
#   seqdecomp character strings -> padded state matrix (miss -> NA)
#   seqdim    c(Rows, Columns) of a sequence object / decompressed strings
#   seqmpos   # matching positions between two sequences (void/nr aware)
#   seqcomp   TRUE iff two sequences have equal length and all positions equal
#   seqsubsn  # distinct subsequences per sequence (DSS by default)
#   seqpm     regex pattern match over seqconc(sep) strings -> nbocc + indices
#
# Reuses the corestats A/B/C fixtures. Text goldens (.txt) carry strings/NA;
# numeric goldens (.csv) use %.15g via write_num_matrix.
#
# Run from the repository root with the oracle library (scripts/r_goldens/README.md):
#   Rscript scripts/r_goldens/make_goldens_seqops.R
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
source(file.path("scripts", "r_goldens", "make_goldens_corestats.R"))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")

write_lines_golden <- function(name, lines, header) {
  con <- file(file.path(out_dir, paste0(name, ".txt")), open = "wt")
  cat(header, file = con)
  writeLines(lines, con)
  close(con)
}

fx <- corestats_fixtures()

## ---- seqconc (default sep "-") ----
emit_conc <- function(sd, tag) {
  s <- as.character(seqconc(sd, sep = "-"))
  write_lines_golden(paste0("seqconc_", tag), s, sprintf("# R: seqconc(%s, sep='-')\n", tag))
}
emit_conc(fx$A, "A")
emit_conc(fx$B, "B")
emit_conc(fx$C, "C")

## Complete public seqconc input, selector, void, and result-metadata surface.
csv_quote <- function(value) {
  paste0('"', gsub('"', '""', as.character(value), fixed = TRUE), '"')
}

seqconc_rows <- list()
append_seqconc_case <- function(rows, case, call, value) {
  for (position in seq_len(nrow(value))) {
    rows[[length(rows) + 1L]] <- paste(
      csv_quote(case),
      csv_quote(call),
      position,
      csv_quote(value[position, 1L]),
      csv_quote(rownames(value)[position]),
      csv_quote(colnames(value)[1L]),
      nrow(value),
      ncol(value),
      sep = ","
    )
  }
  rows
}

for (tag in c("A", "B", "C")) {
  value <- fx[[tag]]
  seqconc_rows <- append_seqconc_case(
    seqconc_rows,
    paste0("sts_", tag),
    paste0("TraMineR::seqconc(fx$", tag, ",void='ignored')"),
    TraMineR::seqconc(value, void = "ignored")
  )
}

raw_vector <- c("A", NA, "B")
seqconc_rows <- append_seqconc_case(
  seqconc_rows,
  "raw_vector",
  "TraMineR::seqconc(c('A',NA,'B'))",
  TraMineR::seqconc(raw_vector)
)

raw_matrix <- matrix(
  c("A", NA, "VOID", "B", "C", "VOID"),
  nrow = 2,
  byrow = TRUE
)
seqconc_rows <- append_seqconc_case(
  seqconc_rows,
  "raw_matrix_void",
  "TraMineR::seqconc(raw_matrix,sep=':',vname='Joined',void='VOID')",
  TraMineR::seqconc(raw_matrix, sep = ":", vname = "Joined", void = "VOID")
)

raw_frame <- data.frame(
  first = c("A", "B"),
  second = c("C", "D"),
  third = c("E", "F"),
  row.names = c("id-1", "id-2"),
  check.names = FALSE
)
seqconc_rows <- append_seqconc_case(
  seqconc_rows,
  "dataframe_named_var",
  "TraMineR::seqconc(raw_frame,var=c('third','first'),sep='|',vname='Selected')",
  TraMineR::seqconc(
    raw_frame,
    var = c("third", "first"),
    sep = "|",
    vname = "Selected"
  )
)
seqconc_rows <- append_seqconc_case(
  seqconc_rows,
  "dataframe_integer_var",
  "TraMineR::seqconc(raw_frame,var=c(3,1),sep='|')",
  TraMineR::seqconc(raw_frame, var = c(3, 1), sep = "|")
)
seqconc_rows <- append_seqconc_case(
  seqconc_rows,
  "missing_only",
  "TraMineR::seqconc(matrix(NA_character_,nrow=1))",
  TraMineR::seqconc(matrix(NA_character_, nrow = 1))
)
seqconc_rows <- append_seqconc_case(
  seqconc_rows,
  "one_row",
  "TraMineR::seqconc(matrix(c('A',NA,'B'),nrow=1))",
  TraMineR::seqconc(
    matrix(
      c("A", NA, "B"),
      nrow = 1,
      dimnames = list("only", c("x", "y", "z"))
    )
  )
)
seqconc_rows <- append_seqconc_case(
  seqconc_rows,
  "zero_columns",
  "TraMineR::seqconc(matrix(character(0),nrow=2,ncol=0))",
  TraMineR::seqconc(matrix(character(0), nrow = 2, ncol = 0))
)
seqconc_rows <- append_seqconc_case(
  seqconc_rows,
  "zero_selector",
  "TraMineR::seqconc(raw_frame,var=0)",
  TraMineR::seqconc(raw_frame, var = 0)
)

seqconc_path <- file.path(out_dir, "seqconc_public_cases.csv")
seqconc_con <- file(seqconc_path, open = "wt", encoding = "UTF-8")
writeLines(
  "case,call,position,value,row_name,column_name,nrow,ncol",
  con = seqconc_con,
  sep = "\n",
  useBytes = TRUE
)
writeLines(
  unlist(seqconc_rows, use.names = FALSE),
  con = seqconc_con,
  sep = "\n",
  useBytes = TRUE
)
close(seqconc_con)

capture_seqconc_error <- function(case, expression) {
  condition <- tryCatch(
    {
      eval.parent(substitute(expression))
      NULL
    },
    error = function(error) error
  )
  if (is.null(condition)) stop("expected seqconc error for case: ", case)
  paste(
    case,
    paste(class(condition), collapse = ","),
    gsub("[\r\n]+", " ", conditionMessage(condition)),
    sep = "|"
  )
}

writeLines(
  c(
    "case|class|message",
    capture_seqconc_error(
      "explicit_null_void",
      TraMineR::seqconc(matrix(c("A", "B"), nrow = 1), void = NULL)
    ),
    capture_seqconc_error(
      "unknown_column",
      TraMineR::seqconc(raw_frame, var = "missing")
    ),
    capture_seqconc_error(
      "vname_length",
      TraMineR::seqconc(raw_matrix, vname = c("A", "B"))
    )
  ),
  con = file.path(out_dir, "seqconc_public_errors.txt"),
  sep = "\n",
  useBytes = TRUE
)

## ---- seqdecomp (on the "-"-joined seqconc strings; ragged -> NA pad) ----
emit_decomp <- function(sd, tag) {
  s <- as.character(seqconc(sd, sep = "-"))
  m <- seqdecomp(as.matrix(s), sep = "-")
  lines <- apply(m, 1, function(r) paste(ifelse(is.na(r), "NA", as.character(r)), collapse = ","))
  write_lines_golden(
    paste0("seqdecomp_", tag), lines,
    sprintf("# R: seqdecomp(seqconc(%s, sep='-'), sep='-'); NA verbatim\n", tag)
  )
}
emit_decomp(fx$A, "A")
emit_decomp(fx$C, "C")

## Complete public seqdecomp tokenization, ordering, and metadata surface.
seqdecomp_rows <- list()
append_seqdecomp_case <- function(rows, case, call, value) {
  for (row_position in seq_len(nrow(value))) {
    for (column_position in seq_len(ncol(value))) {
      cell <- value[row_position, column_position]
      rows[[length(rows) + 1L]] <- paste(
        csv_quote(case),
        csv_quote(call),
        row_position,
        column_position,
        csv_quote(if (is.na(cell)) "" else cell),
        is.na(cell),
        csv_quote(rownames(value)[row_position]),
        csv_quote(colnames(value)[column_position]),
        nrow(value),
        ncol(value),
        sep = ","
      )
    }
  }
  rows
}

append_decomp <- function(case, call, value) {
  seqdecomp_rows <<- append_seqdecomp_case(
    seqdecomp_rows,
    case,
    call,
    value
  )
}

append_decomp(
  "standard",
  "TraMineR::seqdecomp(c('A-B-C','D-E'),sep='-')",
  TraMineR::seqdecomp(c("A-B-C", "D-E"), sep = "-")
)
append_decomp(
  "colon",
  "TraMineR::seqdecomp(c('A:B','C:D:E'),sep=':')",
  TraMineR::seqdecomp(c("A:B", "C:D:E"), sep = ":")
)
append_decomp(
  "empty_separator",
  "TraMineR::seqdecomp(c('ABC','DE'),sep='')",
  TraMineR::seqdecomp(c("ABC", "DE"), sep = "")
)
append_decomp(
  "regex_pipe",
  "TraMineR::seqdecomp(c('|A||B|','C|D'),sep='\\\\|')",
  TraMineR::seqdecomp(c("|A||B|", "C|D"), sep = "\\|")
)
append_decomp(
  "regex_dot",
  "TraMineR::seqdecomp(c('ABC','D'),sep='.')",
  TraMineR::seqdecomp(c("ABC", "D"), sep = ".")
)
append_decomp(
  "multiple_missing",
  "TraMineR::seqdecomp(c('A-X-NA','B-C'),miss=c('X','NA'))",
  TraMineR::seqdecomp(c("A-X-NA", "B-C"), miss = c("X", "NA"))
)
append_decomp(
  "leading_separator",
  "TraMineR::seqdecomp(c('-A-B','C-D'))",
  TraMineR::seqdecomp(c("-A-B", "C-D"))
)
append_decomp(
  "trailing_separator",
  "TraMineR::seqdecomp(c('A-B-','C-D'))",
  TraMineR::seqdecomp(c("A-B-", "C-D"))
)
append_decomp(
  "repeated_separator",
  "TraMineR::seqdecomp(c('A--B','C-D'))",
  TraMineR::seqdecomp(c("A--B", "C-D"))
)
append_decomp(
  "empty_string_mixed",
  "TraMineR::seqdecomp(c('','A-B'))",
  TraMineR::seqdecomp(c("", "A-B"))
)

decomp_frame <- data.frame(
  first = c("A-B", "C-D"),
  second = c("E-F", "G-H"),
  row.names = c("r1", "r2"),
  check.names = FALSE
)
append_decomp(
  "dataframe_selected",
  "TraMineR::seqdecomp(decomp_frame,var='second',vnames=c('left','right'))",
  TraMineR::seqdecomp(
    decomp_frame,
    var = "second",
    vnames = c("left", "right")
  )
)

decomp_matrix <- matrix(
  c("A-B", "C-D", "E-F", "G-H"),
  nrow = 2,
  byrow = TRUE
)
append_decomp(
  "matrix_column_major",
  "TraMineR::seqdecomp(decomp_matrix)",
  TraMineR::seqdecomp(decomp_matrix)
)
append_decomp(
  "matrix_var_order",
  "TraMineR::seqdecomp(decomp_matrix,var=c(2,1))",
  TraMineR::seqdecomp(decomp_matrix, var = c(2, 1))
)

seqdecomp_path <- file.path(out_dir, "seqdecomp_public_cases.csv")
seqdecomp_con <- file(seqdecomp_path, open = "wt", encoding = "UTF-8")
writeLines(
  paste(
    "case",
    "call",
    "row_position",
    "column_position",
    "value",
    "is_na",
    "row_name",
    "column_name",
    "nrow",
    "ncol",
    sep = ","
  ),
  con = seqdecomp_con,
  sep = "\n",
  useBytes = TRUE
)
writeLines(
  unlist(seqdecomp_rows, use.names = FALSE),
  con = seqdecomp_con,
  sep = "\n",
  useBytes = TRUE
)
close(seqdecomp_con)

capture_seqdecomp_error <- function(case, expression) {
  condition <- tryCatch(
    {
      eval.parent(substitute(expression))
      NULL
    },
    error = function(error) error
  )
  if (is.null(condition)) stop("expected seqdecomp error for case: ", case)
  paste(
    case,
    paste(class(condition), collapse = ","),
    gsub("[\r\n]+", " ", conditionMessage(condition)),
    sep = "|"
  )
}

writeLines(
  c(
    "case|class|message",
    capture_seqdecomp_error(
      "empty_string_only",
      TraMineR::seqdecomp("", sep = "-")
    ),
    capture_seqdecomp_error(
      "dataframe_multicolumn_rows",
      TraMineR::seqdecomp(decomp_frame)
    ),
    capture_seqdecomp_error(
      "vnames_length",
      TraMineR::seqdecomp(c("A-B", "C-D"), vnames = "only")
    ),
    capture_seqdecomp_error(
      "unknown_column",
      TraMineR::seqdecomp(decomp_frame, var = "missing")
    )
  ),
  con = file.path(out_dir, "seqdecomp_public_errors.txt"),
  sep = "\n",
  useBytes = TRUE
)

## ---- seqdim (on the "-"-joined strings -> decompressed dims) ----
emit_dim <- function(sd, tag) {
  s <- as.character(seqconc(sd, sep = "-"))
  d <- seqdim(as.matrix(s))
  write_num_matrix(paste0("seqdim_", tag), matrix(d, nrow = 1), sprintf("# R: seqdim(seqconc(%s,'-'))\n", tag))
}
emit_dim(fx$A, "A")
emit_dim(fx$C, "C")

## ---- seqmpos (pairwise; +with.missing for C) ----
emit_mpos <- function(sd, tag, wm) {
  n <- nrow(sd)
  m <- matrix(0L, n, n)
  for (i in 1:n) for (j in 1:n) m[i, j] <- seqmpos(sd[i, ], sd[j, ], with.missing = wm)
  write_num_matrix(paste0("seqmpos_", tag), m, sprintf("# R: seqmpos(%s[i,],%s[j,],with.missing=%s)\n", tag, tag, wm))
}
emit_mpos(fx$A, "A", FALSE)
emit_mpos(fx$B, "B_wm0", FALSE)
emit_mpos(fx$B, "B_wm1", TRUE)
emit_mpos(fx$C, "C_wm0", FALSE)
emit_mpos(fx$C, "C_wm1", TRUE)

## ---- seqcomp (pairwise bool -> 1/0) ----
emit_comp <- function(sd, tag) {
  n <- nrow(sd)
  m <- matrix(0L, n, n)
  for (i in 1:n) for (j in 1:n) m[i, j] <- as.integer(seqcomp(sd[i, ], sd[j, ]))
  write_num_matrix(paste0("seqcomp_", tag), m, sprintf("# R: as.integer(seqcomp(%s[i,],%s[j,]))\n", tag, tag))
}
emit_comp(fx$A, "A")
emit_comp(fx$C, "C")

## ---- seqsubsn (DSS=TRUE/FALSE, +with.missing) ----
emit_subsn <- function(sd, tag, dss, wm) {
  m <- seqsubsn(sd, DSS = dss, with.missing = wm)
  write_num_matrix(paste0("seqsubsn_", tag), m, sprintf("# R: seqsubsn(%s, DSS=%s, with.missing=%s)\n", tag, dss, wm))
}
emit_subsn(fx$A, "A_dss", TRUE, FALSE)
emit_subsn(fx$A, "A_full", FALSE, FALSE)
emit_subsn(fx$C, "C_dss_wm0", TRUE, FALSE)
emit_subsn(fx$C, "C_dss_wm1", TRUE, TRUE)

## ---- seqpm (nbocc + 1-based indices) ----
emit_pm <- function(sd, tag, pattern, sep = "") {
  res <- seqpm(sd, pattern, sep = sep)
  lines <- c(as.character(res$MTab$nbocc), paste(res$MIndex, collapse = " "))
  write_lines_golden(
    paste0("seqpm_", tag), lines,
    sprintf("# R: seqpm(%s, '%s', sep='%s'); line1=nbocc line2=MIndex(1-based)\n", tag, pattern, sep)
  )
}
emit_pm(fx$A, "A_BC", "BC", "")
emit_pm(fx$A, "A_DD", "DD", "")
emit_pm(fx$A, "A_A.C", "A.C", "")

cat("SEQOPS_GOLDENS_DONE\n")
