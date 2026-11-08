# Golden generation: seqistatd (T8) -- per-sequence state durations (counts).
#
# Pins traminex.seqistatd to TraMineR 2.2-14 seqistatd(with.missing): n x nbstat
# counts of each state in statl (alphabet, +nr if with.missing). void is never a
# column. Fixture A is tested with.missing=FALSE only (R fabricates an all-zero nr
# column under with.missing=TRUE for missing-free data; traminex does not).
#
# Run from the repository root with the oracle library (scripts/r_goldens/README.md):
#   Rscript scripts/r_goldens/make_goldens_seqistatd.R
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
source(file.path("scripts", "r_goldens", "make_goldens_corestats.R"))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

emit <- function(sd, tag, wm) {
  m <- seqistatd(sd, with.missing = wm)
  hdr <- sprintf("# R: seqistatd(%s, with.missing=%s); n x nbstat counts\n", tag, wm)
  write_num_matrix(paste0("seqistatd_", tag), m, hdr)
}

fx <- corestats_fixtures()
emit(fx$A, "A", FALSE)
emit(fx$B, "B_wm0", FALSE)
emit(fx$B, "B_wm1", TRUE)
emit(fx$C, "C_wm0", FALSE)
emit(fx$C, "C_wm1", TRUE)

## ---- R2: prop= branches, progress message, and empty-shape disposition ----
## Long form keeps the R state labels explicit and lets the Python test compare
## only TraMineR's statl columns (Python additionally carries a void column).
r2_rows <- list()
append_r2 <- function(case, sd, wm) {
  counts <- suppressMessages(seqistatd(sd, with.missing = wm, prop = FALSE))
  props <- suppressMessages(seqistatd(sd, with.missing = wm, prop = TRUE))
  for (row in seq_len(nrow(counts))) {
    for (col in seq_len(ncol(counts))) {
      r2_rows[[length(r2_rows) + 1L]] <<- data.frame(
        case = case,
        row = row,
        state = colnames(counts)[col],
        count = counts[row, col],
        proportion = props[row, col],
        stringsAsFactors = FALSE
      )
    }
  }
}

# R represents a nonzero-width all-void sequence as an stslist whose cells are
# the declared void code.  Preserve the fixture metadata after slicing.
all_void <- fx$A[1:2, 1:3, drop = FALSE]
all_void[,] <- attr(fx$A, "void")
for (name in c("alphabet", "states", "labels", "nr", "void")) {
  attr(all_void, name) <- attr(fx$A, name)
}
class(all_void) <- class(fx$A)

append_r2("dense_wm0", fx$A, FALSE)
append_r2("dense_wm1", fx$A, TRUE)
append_r2("missing_wm0", fx$B, FALSE)
append_r2("missing_wm1", fx$B, TRUE)
append_r2("ragged_wm0", fx$C, FALSE)
append_r2("ragged_wm1", fx$C, TRUE)
append_r2("all_void_wm0", all_void, FALSE)

r2 <- do.call(rbind, r2_rows)
utils::write.csv(
  r2,
  file.path(out_dir, "seqistatd_r2_prop.csv"),
  row.names = FALSE,
  na = "NaN",
  quote = FALSE
)

message_lines <- capture.output(
  invisible(seqistatd(fx$A, with.missing = FALSE, prop = FALSE)),
  type = "message"
)
zero_width <- fx$A[, FALSE, drop = FALSE]
zero_disposition <- function(prop) {
  tryCatch(
    {
      suppressMessages(seqistatd(zero_width, prop = prop))
      "VALUE"
    },
    error = function(e) paste0("ERROR:", conditionMessage(e))
  )
}
writeLines(
  c(
    paste0("progress=", message_lines),
    paste0("zero_width_count=", zero_disposition(FALSE)),
    paste0("zero_width_prop=", zero_disposition(TRUE))
  ),
  file.path(out_dir, "seqistatd_r2_messages.txt")
)
cat("SEQISTATD_GOLDENS_DONE\n")
