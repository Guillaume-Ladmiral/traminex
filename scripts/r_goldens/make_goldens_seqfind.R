# Public TraMineR 2.2-14 seqfind reference-major matching evidence.
#
# Run from the repository root:
#   /usr/bin/Rscript scripts/r_goldens/make_goldens_seqfind.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))

if (!identical(as.character(packageVersion("TraMineR")), "2.2.14")) {
  stop("seqfind goldens require TraMineR 2.2-14")
}

reference_raw <- rbind(
  c("A", "B", "*", "*"),
  c("A", "B", "*", "*"),
  c("C", "*", "A", "*"),
  c("B", "A", "*", "*"),
  c("A", "C", "*", "*")
)
target_raw <- rbind(
  c("A", "B", "*"),
  c("C", "*", "A"),
  c("A", "B", "*"),
  c("B", "A", "*"),
  c("C", "C", "*")
)

reference <- suppressMessages(
  TraMineR::seqdef(
    reference_raw,
    alphabet = c("A", "B", "C"),
    missing = "*"
  )
)
target <- suppressMessages(
  TraMineR::seqdef(
    target_raw,
    alphabet = c("C", "B", "A"),
    missing = "*"
  )
)

reference_labels <- as.vector(
  TraMineR::seqconc(reference, void = attr(reference, "void"))
)
target_labels <- as.vector(
  TraMineR::seqconc(target, void = attr(target, "void"))
)
indexes <- TraMineR::seqfind(reference, target)
no_match <- TraMineR::seqfind(reference[5, ], target)

stopifnot(
  identical(indexes, c(1L, 3L, 1L, 3L, 2L, 4L)),
  identical(no_match, integer(0))
)

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

lines <- "object|row|width|conc|alphabet|void"
for (i in seq_along(reference_labels)) {
  lines <- c(
    lines,
    paste(
      "reference",
      i,
      ncol(reference),
      reference_labels[[i]],
      paste(attr(reference, "alphabet"), collapse = ","),
      attr(reference, "void"),
      sep = "|"
    )
  )
}
for (i in seq_along(target_labels)) {
  lines <- c(
    lines,
    paste(
      "target",
      i,
      ncol(target),
      target_labels[[i]],
      paste(attr(target, "alphabet"), collapse = ","),
      attr(target, "void"),
      sep = "|"
    )
  )
}
writeLines(
  lines,
  con = file.path(out_dir, "seqfind_public_inputs.txt"),
  sep = "\n",
  useBytes = TRUE
)

writeLines(
  c(
    "case|call|r_indexes",
    paste(
      "reference_major",
      "TraMineR::seqfind(reference,target)",
      paste(indexes, collapse = ","),
      sep = "|"
    ),
    paste(
      "no_match",
      "TraMineR::seqfind(reference[5,],target)",
      paste(no_match, collapse = ","),
      sep = "|"
    )
  ),
  con = file.path(out_dir, "seqfind_public_results.txt"),
  sep = "\n",
  useBytes = TRUE
)

cat("SEQFIND_GOLDENS_DONE\n")
