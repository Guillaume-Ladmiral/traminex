# Golden generation: seqdist's with.missing reset (seqdist.R:111-118).
#
# R resets with.missing to FALSE when neither seqdata nor a sequence-type
# refseq contains a missing value, and stops when missing values are present
# without with.missing.  traminex's plain-integer-array input has no R
# counterpart, so its R analogue is seqdef(X, alphabet = 0:2): the same codes
# with no missing cell.  These fixtures pin, for every method that resolves a
# cost-method string through seqcost, the full distance matrix with
# with.missing = TRUE and FALSE (identical in R, by the reset), and the
# distances to an external reference that does carry a missing cell.
#
# Consumed by tests/test_seqdist_plain_array_with_missing.py.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

X <- matrix(
  c(0, 1, 2, 1,
    1, 1, 0, 2,
    2, 0, 1, 1,
    0, 2, 2, 0,
    1, 0, 0, 1),
  ncol = 4, byrow = TRUE
)
s <- suppressMessages(TraMineR::seqdef(X, alphabet = c(0, 1, 2)))

# Method-specific required arguments; every other argument is R's default.
method_args <- list(
  OM = list(),
  OMslen = list(),
  OMspell = list(),
  OMloc = list(),
  OMstran = list(otto = 0.5),
  TWED = list(nu = 0.5, h = 0.5),
  HAM = list()
)
cost_methods <- c("CONSTANT", "TRATE", "INDELS", "INDELSLOG")

fmt <- function(x) sprintf("%.17g", x)
rows <- character()
add_rows <- function(m, sm, wm, indel_label, d) {
  for (i in seq_len(nrow(d))) {
    for (j in seq_len(ncol(d))) {
      rows <<- c(rows, paste(m, sm, if (wm) "TRUE" else "FALSE", indel_label, i - 1, j - 1, fmt(d[i, j]), sep = ","))
    }
  }
}
for (m in names(method_args)) {
  for (sm in cost_methods) {
    for (wm in c(TRUE, FALSE)) {
      args <- c(list(seqdata = s, method = m, sm = sm, with.missing = wm), method_args[[m]])
      add_rows(m, sm, wm, "auto", as.matrix(suppressMessages(do.call(TraMineR::seqdist, args))))
    }
  }
}
# An explicit per-state indel vector over the three real states.  Its length
# must equal nstates *after* the reset (seqdist.R:131-135), so it is accepted
# with with.missing = TRUE only because R resets the flag.  TWED warns and uses
# max(indel); the warning is not part of this fixture.
indel_vector <- c(1, 1.5, 2)
for (m in c("OM", "OMslen", "OMspell", "OMstran", "TWED")) {
  for (wm in c(TRUE, FALSE)) {
    args <- c(
      list(seqdata = s, method = m, sm = "CONSTANT", indel = indel_vector, with.missing = wm),
      method_args[[m]]
    )
    d <- as.matrix(suppressWarnings(suppressMessages(do.call(TraMineR::seqdist, args))))
    add_rows(m, "CONSTANT", wm, "1;1.5;2", d)
  }
}

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
writeLines(
  c("method,sm,with_missing,indel,row,col,value", rows),
  file.path(out_dir, "seqdist_with_missing_reset.csv")
)

# External reference carrying one missing cell ('*').  R forbids
# with.missing = FALSE here and widens the cost matrix when it is TRUE.
# INDELS is omitted: the missing state has zero frequency in seqdata, so R's
# result saturates at .Machine$double.xmax, a separate contract.
ref <- suppressMessages(TraMineR::seqdef(matrix(c(0, 1, NA, 2), nrow = 1), alphabet = c(0, 1, 2)))
stopifnot(identical(attr(ref, "nr"), attr(s, "nr")))
ref_rows <- character()
for (sm in c("CONSTANT", "TRATE", "INDELSLOG")) {
  for (wm in c(TRUE, FALSE)) {
    result <- tryCatch(
      suppressMessages(TraMineR::seqdist(s, method = "OM", sm = sm, refseq = ref, with.missing = wm)),
      error = function(e) conditionMessage(e)
    )
    if (is.character(result)) {
      ref_rows <- c(ref_rows, paste(sm, if (wm) "TRUE" else "FALSE", "error", NA, sep = ","))
    } else {
      for (i in seq_along(result)) {
        ref_rows <- c(ref_rows, paste(sm, if (wm) "TRUE" else "FALSE", i - 1, fmt(result[i]), sep = ","))
      }
    }
  }
}
writeLines(
  c("sm,with_missing,row,value", ref_rows),
  file.path(out_dir, "seqdist_with_missing_reset_refseq.csv")
)
# The refusal message itself, recorded once (identical for every sm).
refusal <- tryCatch(
  suppressMessages(TraMineR::seqdist(s, method = "OM", sm = "CONSTANT", refseq = ref, with.missing = FALSE)),
  error = function(e) conditionMessage(e)
)
writeLines(refusal, file.path(out_dir, "seqdist_with_missing_reset_refusal.txt"))

# A reference whose alphabet lacks a state seqdata uses: the reference's missing
# code would coincide with that state's code.  R refuses the pair outright.
s3 <- suppressMessages(TraMineR::seqdef(
  matrix(c(0, 1, 2, 0, 1, 2, 0, 0, 0, 0, 1, 2), ncol = 4, byrow = TRUE)
))
narrow_ref <- suppressMessages(TraMineR::seqdef(matrix(c(0, 1, NA, 1), nrow = 1), alphabet = c(0, 1)))
alphabet_refusal <- tryCatch(
  suppressMessages(TraMineR::seqdist(s3, method = "OM", sm = "CONSTANT", refseq = narrow_ref, with.missing = TRUE)),
  error = function(e) conditionMessage(e)
)
stopifnot(is.character(alphabet_refusal))
writeLines(alphabet_refusal, file.path(out_dir, "seqdist_with_missing_reset_alphabet_refusal.txt"))

cat("SEQDIST_WITH_MISSING_RESET_GOLDENS_DONE\n")
