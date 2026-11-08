# Golden generation: argument forms seqdist accepts, maps or refuses.
#
# seqdist.R maps the deprecated methods OMopt/LCSopt to OM/LCS and a logical
# norm to "auto"/"none" (each with a warning), uses the first two elements of
# a longer refseq list (with a warning), refuses refseq = TRUE, and -- because
# is.a.positive.integer counts 0 and R indexing drops it -- silently ignores a
# 0 inside refseq index sets.  traminex refused the first three and accepted
# TRUE as row 1 (2026-09-22 audit S1a-12).
#
# Consumed by tests/test_seqdist_argument_forms.py.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

X <- matrix(c(0, 1, 2, 1, 1, 1, 0, 2, 2, 0, 1, 1, 0, 2, 2, 0), ncol = 4, byrow = TRUE)
s <- suppressMessages(seqdef(X, alphabet = 0:2))
lines <- character()
record <- function(case, expr) {
  seen <- character()
  result <- tryCatch(
    withCallingHandlers(suppressMessages(expr), warning = function(w) {
      seen <<- c(seen, trimws(conditionMessage(w)))
      invokeRestart("muffleWarning")
    }),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    lines <<- c(lines, paste(case, "error", trimws(conditionMessage(result)), sep = "\t"))
  } else {
    lines <<- c(lines, paste(
      case, "value", paste(sprintf("%.17g", as.vector(as.matrix(result))), collapse = " "),
      if (length(seen)) paste(unique(seen), collapse = " | ") else "no warning",
      sep = "\t"
    ))
  }
}
record("om", seqdist(s, method = "OM", sm = "CONSTANT"))
record("omopt", seqdist(s, method = "OMopt", sm = "CONSTANT"))
record("lcs", seqdist(s, method = "LCS"))
record("lcsopt", seqdist(s, method = "LCSopt"))
record("norm_auto", seqdist(s, method = "OM", sm = "CONSTANT", norm = "auto"))
record("norm_true", seqdist(s, method = "OM", sm = "CONSTANT", norm = TRUE))
record("norm_false", seqdist(s, method = "OM", sm = "CONSTANT", norm = FALSE))
record("refseq_two_sets", seqdist(s, method = "OM", sm = "CONSTANT", refseq = list(c(1, 2), 3)))
record("refseq_three_sets", seqdist(s, method = "OM", sm = "CONSTANT", refseq = list(c(1, 2), 3, 4)))
record("refseq_true", seqdist(s, method = "OM", sm = "CONSTANT", refseq = TRUE))
record("refseq_set_with_zero", seqdist(s, method = "OM", sm = "CONSTANT", refseq = list(c(0, 1), 3)))

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
writeLines(c("# case\tkind\tR result (column-major) or error\tR warnings", lines),
  file.path(out_dir, "seqdist_argument_forms.txt"))
cat("SEQDIST_ARGUMENT_FORMS_GOLDENS_DONE\n")
