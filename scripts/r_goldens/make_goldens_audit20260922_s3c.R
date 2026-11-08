# Golden generation: follow-up statistics findings of the 2026-09-22 audit
# (remediation package s3c: S3-20, S3-21, S3-22).
#
# Each record is one R call on a small hand-built object: the case name, then
# "text" and the character result joined by "|", "value" and the numeric
# result (column-major, %.17g, so NaN/NA are spelled out), or "error" and R's
# message; the last field lists R's warnings ("none" when there were none).
#
# S3-21: seqtab.R:12 refuses NULL, a negative value, and a vector longer than
#   one containing a value below 1; a single value is a scalar, so c(0) is 0
#   and means every rank (seqtab.R:50), as does any rank beyond the number of
#   distinct sequences. Records are the row names of attr(, "freq").
# S3-20: seqimplic's levels are levels(factor(group)), i.e. as.character of
#   the sorted unique values: 15 significant digits for doubles ("2", "1e+05",
#   and 0.1 + 0.2 and 0.3 share the level "0.3"), "FALSE"/"TRUE" for logicals.
# S3-22: seqsurv's groups are factor(groups) (seqsurv.R:38): character levels
#   in R's collation (a before A before b), numbers in numeric order with
#   as.character labels, logicals as FALSE/TRUE. Records are the curve labels
#   (ltext), the strata in fit order, and each curve's summary rows
#   (row-major: time, n.risk, n.event, surv).
#
# Consumed by tests/test_seqtab_single_rank_vector.py (S3-21),
# tests/test_seqimplic_number_labels.py (S3-20) and
# tests/test_seqsurv_group_levels.py (S3-22).
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

lines <- character()
record_text <- function(case, expr) {
  seen <- character()
  result <- tryCatch(
    withCallingHandlers(suppressMessages(expr), warning = function(w) {
      seen <<- c(seen, trimws(conditionMessage(w)))
      invokeRestart("muffleWarning")
    }),
    error = function(e) e
  )
  warned <- if (length(seen)) paste(unique(seen), collapse = " | ") else "none"
  if (inherits(result, "error")) {
    lines <<- c(lines, paste(case, "error", trimws(conditionMessage(result)), warned, sep = "\t"))
  } else {
    values <- as.character(result)
    payload <- if (length(values)) paste(values, collapse = "|") else "none"
    lines <<- c(lines, paste(case, "text", payload, warned, sep = "\t"))
  }
}

record_values <- function(case, values) {
  payload <- if (length(values)) paste(sprintf("%.17g", as.vector(values)), collapse = " ") else "none"
  lines <<- c(lines, paste(case, "value", payload, "none", sep = "\t"))
}

# ---- S3-21: seqtab idxs forms ------------------------------------------ #
tab_s <- suppressMessages(seqdef(
  rbind(c("A", "A", "B", "C"), c("B", "B", "C", "C"), c("A", "B", "A", "B"), c("A", "A", "B", "C")),
  alphabet = c("A", "B", "C")
))
tab_ranks <- function(idxs) rownames(attr(seqtab(tab_s, idxs = idxs), "freq"))
record_text("s321_idxs_scalar_0", tab_ranks(0))
record_text("s321_idxs_vector_0", tab_ranks(c(0)))
record_text("s321_idxs_vector_2", tab_ranks(c(2)))
record_text("s321_idxs_vector_7", tab_ranks(c(7)))
record_text("s321_idxs_vector_2_1", tab_ranks(c(2, 1)))
record_text("s321_idxs_vector_2_9", tab_ranks(c(2, 9)))
record_text("s321_idxs_vector_minus_1", tab_ranks(c(-1)))
record_text("s321_idxs_vector_0_1", tab_ranks(c(0, 1)))

# ---- S3-20: seqimplic labels of numeric and logical groups ------------ #
suppressMessages(library(TraMineRextras))
stopifnot(as.character(packageVersion("TraMineRextras")) == "0.6.9")
implic_s <- suppressMessages(seqdef(
  rbind(c("A", "B", "A"), c("B", "A", "A"), c("A", "A", "B"), c("B", "B", "B")),
  alphabet = c("A", "B")
))
implic_groups <- list(
  fractional = c(1.5, 2, 10, 2),
  whole_doubles = c(2, 10, 2, 10),
  large = c(1e5, 2, 1e5, 2),
  float_neighbours = c(0.1 + 0.2, 0.3, 1, 1),
  logical = c(TRUE, FALSE, TRUE, FALSE)
)
for (name in names(implic_groups)) {
  res <- suppressWarnings(seqimplic(implic_s, group = implic_groups[[name]]))
  record_text(paste0("s320_", name, "_levels"), res$levels)
  record_values(paste0("s320_", name, "_indices"), res$indices)
}

# ---- S3-22: seqsurv group levels --------------------------------------- #
# `survival` is a recommended package the isolated oracle library cannot pin
# (see make_goldens_seqsurv.R); assert the documented floor and report the
# resolved version, which check_r_golden_provenance.py records as provenance.
suppressMessages(library(survival))
stopifnot(packageVersion("survival") >= "3.8.9")
cat("survival version:", as.character(packageVersion("survival")), "\n")
surv_s <- suppressMessages(seqdef(
  rbind(
    c("A", "A", "B", "B"), c("B", "B", "A", "A"), c("A", "B", "B", "B"),
    c("B", "A", "A", "A"), c("A", "A", "A", "B"), c("B", "B", "B", "A")
  ),
  alphabet = c("A", "B")
))
surv_groups <- list(
  letter_case = c("b", "A", "a", "b", "A", "a"),
  mixed_case_words = c("ab", "Ab", "aB", "ab", "Ab", "aB"),
  fractional = c(1.5, 2, 10, 2, 1.5, 10),
  logical = c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)
)
for (name in names(surv_groups)) {
  fit <- suppressMessages(seqsurv(surv_s, groups = surv_groups[[name]], per.state = TRUE))
  prefix <- paste0("s322_", name)
  record_text(paste0(prefix, "_ltext"), attr(fit, "ltext"))
  strata <- sub("spell$group=", "", names(fit$strata), fixed = TRUE)
  record_text(paste0(prefix, "_strata"), strata)
  s <- summary(fit)
  curve <- sub("spell$group=", "", as.character(s$strata), fixed = TRUE)
  for (i in seq_along(strata)) {
    take <- curve == strata[[i]]
    record_values(
      paste0(prefix, "_curve", i),
      t(cbind(s$time, s$n.risk, s$n.event, s$surv)[take, , drop = FALSE])
    )
  }
}

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
writeLines(c("# case\tkind (text, value or error)\tpayload\tR warnings", lines),
  file.path(out_dir, "audit20260922_s3c.txt"))
cat("AUDIT20260922_S3C_GOLDENS_DONE\n")
