# Golden generation: degenerate and boundary cases of the sequence statistics
# (2026-09-22 whole-codebase audit, remediation package s3b).
#
# Each record is one R call on a small hand-built object: the case name, then
# "value" and the result (column-major, %.17g, so NaN/NA/Inf are spelled out),
# "error" and R's message, or "text" and character values joined by "|"; the
# last field lists R's warnings ("none" when there were none).
#
# S3-11: a sequence with no valid state (every cell missing, with.missing =
#   FALSE). seqivardur type 1 is NaN (the variance of an empty vector); type 2
#   counts the K unvisited states as zero-length spells, so the value and mean
#   are 0 and vmax is 0/(K-2): 0 for K >= 3 and NaN for K <= 2. seqST inherits
#   these through log2(phi * (vmax + 1) / (var + 1)).
# S3-10: degenerate normalisers. With a one-state alphabet seqient divides by
#   log(1) = 0 and seqST(norm = TRUE) by maxT - 1 = 0; when every sequence has
#   length 1, maxT is again 1. R returns the resulting NaN (0/0); seqici and
#   seqindic inherit it.
# S3-13: seqauto's table(srs[, nc], srs[, nc - i]) keeps every level either
#   column takes, including one seen only beside an NA in the other column; an
#   all-zero row or column then makes chisq.test's statistic 0/0 = NaN. A table
#   with one row or column becomes a goodness-of-fit test against equal
#   probabilities (Cramer's V is then Inf), and a 1 x 1 table is refused.
# S3-14: seqsurv with zero case weights. summary.survfit keeps only the times
#   whose (weighted) number of events is positive, so a time at which only
#   zero-weight spells end is not a row. Each record is one curve's summary
#   rows, row-major: time, n.risk, n.event, surv, std.err, lower, upper.
# S3-15: seqimplic orders group levels as factor(group) does (numeric order for
#   numbers, the locale collation for text, the declared order for a factor),
#   and with.missing = TRUE always appends the missing state to the alphabet
#   (label "Missing", the missing colour, an all-NA slice) even when no cell
#   is missing. Text records join values with "|"; the indices are
#   column-major (level fastest), as R's as.vector.
# S3-16: the "weighted" attribute of seqmodst and seqmeant is FALSE when the
#   object has no weights, when weighted = FALSE, or when every weight is 1.
# S3-17: argument checks. seqST passes type to seqivardur, which stops unless
#   it is 1 or 2; seqivolatility stops on a non-logical adjust; seqintegr
#   stops on more than one state or a state outside alphabet(with.missing)
#   (the missing state is accepted only with with.missing = TRUE, even when
#   no cell is missing); seqipos maps the obsolete index "integration" to
#   "integr" (with a message, not a warning).
# S3-12 (recorded divergence): seqipos(index = "volatility") recodes the states
#   outside pos.states/neg.states (the missing state when with.missing =
#   FALSE) to void, which can then sit inside a sequence; seqdss reads only the
#   first seqlength cells, so the spells after the gap are lost. traminex
#   counts every observed cell; its expected value is R's seqivolatility on the
#   same binary sequence without the gap.
#
# Consumed by tests/test_degenerate_statistics_nan.py (S3-10, S3-11),
# tests/test_seqauto_table_levels.py (S3-13) and
# tests/test_seqsurv_zero_weights.py (S3-14) and
# tests/test_seqimplic_levels_and_missing.py (S3-15) and
# tests/test_weighted_metadata_flag.py (S3-16) and
# tests/test_statistics_argument_checks.py (S3-17) and
# tests/test_seqipos_volatility_gaps.py (S3-12).
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
suppressMessages(library(TraMineRextras))
stopifnot(as.character(packageVersion("TraMineRextras")) == "0.6.9")
# `survival` is a recommended package that the isolated oracle library cannot
# pin (see make_goldens_seqsurv.R); assert the documented floor and report the
# resolved version, which check_r_golden_provenance.py records as provenance.
suppressMessages(library(survival))
stopifnot(packageVersion("survival") >= "3.8.9")
cat("survival version:", as.character(packageVersion("survival")), "\n")

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
  warned <- if (length(seen)) paste(unique(seen), collapse = " | ") else "none"
  if (inherits(result, "error")) {
    lines <<- c(lines, paste(case, "error", trimws(conditionMessage(result)), warned, sep = "\t"))
  } else {
    values <- as.vector(as.matrix(result))
    payload <- if (length(values)) paste(sprintf("%.17g", values), collapse = " ") else "none"
    lines <<- c(lines, paste(case, "value", payload, warned, sep = "\t"))
  }
}
record_text <- function(case, values) {
  payload <- if (length(values)) paste(as.character(values), collapse = "|") else "none"
  lines <<- c(lines, paste(case, "text", payload, "none", sep = "\t"))
}
sts <- function(rows, alphabet) {
  suppressMessages(seqdef(
    as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE),
    alphabet = alphabet, right = NA
  ))
}

# ---- S3-11: rows with no valid state ------------------------------------ #
empty_rows <- list(
  k1 = sts(list(c("A", "A", "A"), c(NA, NA, NA)), "A"),
  k2 = sts(list(c("A", "B", "B"), c(NA, NA, NA)), c("A", "B")),
  k3 = sts(list(c("A", NA, "B", "B"), c(NA, NA, NA, NA), c("B", NA, "C", "C")), c("A", "B", "C")),
  k4 = sts(list(c("A", "B", "C", "D"), c(NA, NA, NA, NA)), c("A", "B", "C", "D"))
)
for (name in names(empty_rows)) {
  s <- empty_rows[[name]]
  for (type in 1:2) {
    prefix <- paste0("s311_", name, "_type", type)
    v <- seqivardur(s, type = type, with.missing = FALSE)
    record(paste0(prefix, "_seqivardur_value"), as.vector(v))
    record(paste0(prefix, "_seqivardur_vmax"), attr(v, "vmax"))
    record(paste0(prefix, "_seqivardur_meand"), attr(v, "meand"))
    record(paste0(prefix, "_seqST"), seqST(s, norm = FALSE, with.missing = FALSE, type = type))
    record(paste0(prefix, "_seqST_norm"), seqST(s, norm = TRUE, with.missing = FALSE, type = type))
  }
  indic <- seqindic(s, indic = c("turb", "turbn", "turb2", "turb2n", "dustd", "dustd2"))
  for (column in colnames(indic)) {
    record(paste0("s311_", name, "_seqindic_", column), indic[[column]])
  }
}

# ---- S3-10: degenerate normalisers -------------------------------------- #
degenerate <- list(
  k1 = sts(list(c("A", "A", "A"), c("A", "A", "A")), "A"),
  k1miss = sts(list(c("A", "A", "A"), c("A", NA, "A")), "A"),
  len1 = sts(list("A", "B", "A"), c("A", "B")),
  len1k1 = sts(list("A", "A"), "A"),
  control = sts(list(c("A", "B", "B", "C"), c("C", "C", "A", "A"), c("A", "A", "A", "A")), c("A", "B", "C"))
)
for (name in names(degenerate)) {
  s <- degenerate[[name]]
  for (wm in c(FALSE, TRUE)) {
    prefix <- paste0("s310_", name, "_wm", as.integer(wm))
    record(paste0(prefix, "_seqient_norm"), seqient(s, norm = TRUE, with.missing = wm))
    record(paste0(prefix, "_seqient_raw"), seqient(s, norm = FALSE, with.missing = wm))
    record(paste0(prefix, "_seqici"), seqici(s, with.missing = wm))
    for (type in 1:2) {
      record(paste0(prefix, "_type", type, "_seqST"), seqST(s, norm = FALSE, with.missing = wm, type = type))
      record(paste0(prefix, "_type", type, "_seqST_norm"), seqST(s, norm = TRUE, with.missing = wm, type = type))
    }
    indic <- seqindic(s, indic = c("entr", "cplx", "turbn", "turb2n"), with.missing = wm)
    for (column in colnames(indic)) {
      record(paste0(prefix, "_seqindic_", column), indic[[column]])
    }
  }
}

# ---- S3-13: seqauto contingency-table levels ---------------------------- #
auto_cases <- list(
  level_seen_only_beside_na = list(rbind(c("X", "A", "A", "B"), c("X", "B", "B", "A"), c("A", "A", "B", "B"), c("B", "A", "A", "A")), c("A", "B", "X"), 1),
  zero_row_after_first_position = list(rbind(c("B", "A", "A"), c("A", "A", "A")), c("A", "B"), 1),
  level_only_in_current_column = list(rbind(c("A", "B", "A", "B"), c("B", "A", "B", "A"), c("A", "A", "B", "B"), c("B", "B", "A", "C")), c("A", "B", "C"), 1),
  one_column = list(rbind(c("A", "A", "B"), c("A", "A", "A")), c("A", "B"), 1),
  one_by_one = list(rbind(c("A", "A", "A"), c("A", "A", "A")), c("A", "B"), 1),
  control_order2 = list(rbind(c("A", "B", "A", "B"), c("B", "A", "B", "A"), c("A", "A", "B", "B"), c("B", "B", "A", "A")), c("A", "B"), 2),
  control_lag2_uneven = list(rbind(c("B", "A", "A", "B"), c("A", "A", "A", "A"), c("A", "B", "A", "A")), c("A", "B"), 2),
  unobserved_alphabet_state = list(rbind(c("A", "B", "A", "B"), c("B", "A", "B", "A"), c("A", "A", "B", "B"), c("B", "B", "A", "A")), c("A", "B", "Z"), 1)
)
for (name in names(auto_cases)) {
  case <- auto_cases[[name]]
  s <- suppressMessages(seqdef(case[[1]], alphabet = case[[2]]))
  # Row-major (lag by lag): auto-association then p-value for each order.
  record(paste0("s313_", name), t(seqauto(s, order = case[[3]])))
}
missing_case <- sts(
  list(c("A", NA, "A", "B"), c("B", "A", "B", "A"), c("A", "A", "B", "B"), c("B", "B", "A", "A")),
  c("A", "B")
)
record("s313_missing_cells", t(seqauto(missing_case, order = 1)))

# ---- S3-14: seqsurv with zero case weights ------------------------------ #
surv_rows <- rbind(
  c("A", "A", "B", "B", "C", "C"), c("A", "B", "B", "C", "C", "C"),
  c("B", "B", "B", "A", "A", "C"), c("C", "C", "A", "A", "B", "B"),
  c("A", "A", "A", "A", "A", "A"), c("B", "C", "B", "C", "B", "C")
)
surv_case <- function(name, fit, keys, prefix) {
  s <- summary(fit)
  curve <- if (is.null(s$strata)) rep(keys[[1]], length(s$time)) else
    sub(prefix, "", as.character(s$strata), fixed = TRUE)
  for (key in keys) {
    take <- curve == key
    record(
      paste0("s314_", name, "_", key),
      t(cbind(s$time, s$n.risk, s$n.event, s$surv, s$std.err, s$lower, s$upper)[take, , drop = FALSE])
    )
  }
}
integer_w <- suppressMessages(seqdef(surv_rows, alphabet = c("A", "B", "C", "D"), weights = c(0, 0, 1, 1, 0, 2)))
surv_case("integer", seqsurv(integer_w), c("A", "B", "C"), "spell$states=")
fractional_w <- suppressMessages(seqdef(surv_rows, alphabet = c("A", "B", "C", "D"), weights = c(0, 0, 1.5, 2, 0, 0.5)))
surv_case("fractional", seqsurv(fractional_w), c("A", "B", "C"), "spell$states=")
groups <- c("g1", "g1", "g2", "g2", "g1", "g2")
surv_case(
  "per_state", seqsurv(integer_w, groups = groups, per.state = TRUE, state = c("A", "B")),
  c("g1", "g2"), "spell$group="
)
after_end <- suppressMessages(seqdef(
  rbind(c("A", "B", "B", "B"), c("A", "A", "A", "B"), c("B", "A", "B", "B")),
  alphabet = c("A", "B"), weights = c(1, 0, 1)
))
surv_case("after_weighted_end", seqsurv(after_end), c("A", "B"), "spell$states=")

# ---- S3-15: seqimplic group levels and the missing state --------------- #
implic_rows <- rbind(c("A", "B", "A"), c("B", "A", "A"), c("A", "A", "B"), c("B", "B", "B"))
implic_s <- suppressMessages(seqdef(implic_rows, alphabet = c("A", "B")))
implic_cases <- list(
  numeric_groups = list(implic_s, c(1, 2, 10, 2), FALSE),
  character_case = list(implic_s, c("b", "A", "a", "b"), FALSE),
  declared_factor = list(implic_s, factor(c("x", "y", "x", "y"), levels = c("y", "x")), FALSE),
  missing_state_absent = list(implic_s, c(1, 2, 1, 2), TRUE),
  missing_state_present = list(
    sts(list(c("A", NA, "A"), c("B", "A", "A"), c("A", "A", NA), c("B", "B", "B")), c("A", "B")),
    c(1, 2, 1, 2), TRUE
  )
)
for (name in names(implic_cases)) {
  case <- implic_cases[[name]]
  res <- suppressWarnings(seqimplic(case[[1]], group = case[[2]], with.missing = case[[3]]))
  record_text(paste0("s315_", name, "_levels"), res$levels)
  record_text(paste0("s315_", name, "_alphabet"), res$alphabet)
  record_text(paste0("s315_", name, "_labels"), res$labels)
  record_text(paste0("s315_", name, "_cpal"), res$cpal)
  record(paste0("s315_", name, "_indices"), as.vector(res$indices))
}
# na.rm = FALSE keeps a row whose group is NA: NA is not a level, and every sum
# over that row is NA, so each position it enters is NA throughout.
na_group <- suppressWarnings(seqimplic(implic_s, group = c(1, NA, 2, 2), na.rm = FALSE))
record_text("s315_na_group_kept_levels", na_group$levels)
record("s315_na_group_kept_indices", as.vector(na_group$indices))

# ---- S3-16: the weighted attribute ------------------------------------- #
flag_rows <- rbind(c("A", "A", "B"), c("B", "B", "A"), c("A", "B", "B"))
flag_objects <- list(
  no_weights = suppressMessages(seqdef(flag_rows, alphabet = c("A", "B"))),
  unit_weights = suppressMessages(seqdef(flag_rows, alphabet = c("A", "B"), weights = c(1, 1, 1))),
  other_weights = suppressMessages(seqdef(flag_rows, alphabet = c("A", "B"), weights = c(1, 2, 0.5)))
)
for (name in names(flag_objects)) {
  s <- flag_objects[[name]]
  for (weighted in c(TRUE, FALSE)) {
    prefix <- paste0("s316_", name, "_weighted", as.integer(weighted))
    modst <- seqmodst(s, weighted = weighted)
    record_text(paste0(prefix, "_seqmodst_weighted"), attr(modst, "weighted"))
    record(paste0(prefix, "_seqmodst_nbseq"), attr(modst, "nbseq"))
    meant <- seqmeant(s, weighted = weighted)
    record_text(paste0(prefix, "_seqmeant_weighted"), attr(meant, "weighted"))
    record(paste0(prefix, "_seqmeant_nbseq"), attr(meant, "nbseq"))
  }
}

# ---- S3-17: argument checks -------------------------------------------- #
check_s <- suppressMessages(seqdef(rbind(c("A", "A", "B", "C"), c("B", "B", "C", "C"), c("A", "B", "A", "B")), alphabet = c("A", "B", "C")))
check_m <- sts(list(c("A", NA, "B", "C"), c("B", "B", NA, "C"), c("A", "B", "A", "B")), c("A", "B", "C"))
check_r <- suppressMessages(seqdef(
  as.data.frame(rbind(c("A", "B", "C", "C"), c("B", "B", NA, NA)), stringsAsFactors = FALSE),
  alphabet = c("A", "B", "C")
))
record("s317_seqST_type3", seqST(check_s, norm = FALSE, with.missing = FALSE, type = 3))
record("s317_seqST_type0_norm", seqST(check_s, norm = TRUE, with.missing = FALSE, type = 0))
record("s317_seqivolatility_adjust_text", seqivolatility(check_s, adjust = "yes"))
record("s317_seqivolatility_adjust_number", seqivolatility(check_s, adjust = 1))
record("s317_seqintegr_nr_without_missing_flag", seqintegr(check_m, state = "*", with.missing = FALSE))
record("s317_seqintegr_nr_with_missing_flag", seqintegr(check_m, state = "*", with.missing = TRUE))
record("s317_seqintegr_nr_absent_from_data", seqintegr(check_s, state = "*", with.missing = TRUE))
record("s317_seqintegr_void", seqintegr(check_r, state = "%"))
record("s317_seqintegr_unknown", seqintegr(check_s, state = "Z"))
record("s317_seqintegr_two_states", seqintegr(check_s, state = c("A", "B")))
record("s317_seqintegr_valid", seqintegr(check_m, state = "A"))
record("s317_seqipos_integration", seqipos(check_s, pos.states = "A", index = "integration"))
record("s317_seqipos_integr", seqipos(check_s, pos.states = "A", index = "integr"))

# ---- TX-AUD-20261019-012: void positions in the integration weights ---- #
# seqintegr.Rd numbers the positions 1..L of the sequence and treats only non-void missing
# values as a state with with.missing = TRUE, but seqintegr.R divides by sum((1:ncol)^pow) then,
# so the void padding of a-b counts: a-b gives 0.1 / 0.2 in data of width 4 and 1/6 / 2/6 in
# data of width 3. R's values recorded; traminex's corrected ones are asserted in
# tests/test_integr_void.py.
integr_x <- suppressMessages(seqdef(rbind(c("a", "b", "b", "a"), c("a", "b", NA, NA), c("b", NA, "a", "a"))))
integr_w3 <- suppressMessages(seqdef(rbind(c("a", "b", NA), c("a", "b", "b"))))
record("s8_012_integr_with_missing", seqintegr(integr_x, with.missing = TRUE))
record("s8_012_integr_without_missing", seqintegr(integr_x, with.missing = FALSE))
record("s8_012_integr_width3_with_missing", seqintegr(integr_w3, with.missing = TRUE))
record("s8_012_ibad_with_missing", seqibad(integr_x, with.missing = TRUE))
record("s8_012_ipos_integr_with_missing", seqipos(integr_x, pos.states = "a", index = "integr", with.missing = TRUE))
# TX-REP-20261025-006: an all-missing row is void throughout (seqdef's right = "DEL"), an empty sequence. With
# with.missing = FALSE R gives NaN (no positions); with with.missing = TRUE its every-column denominator gives 0.
integr_e <- suppressMessages(suppressWarnings(seqdef(rbind(c("a", "b", NA, NA), c(NA, NA, NA, NA), c("b", "b", "a", "b")),
  alphabet = c("a", "b"))))
record("s8_rep006_integr_empty_without_missing", seqintegr(integr_e, state = "a", with.missing = FALSE))
record("s8_rep006_integr_empty_with_missing", seqintegr(integr_e, state = "a", with.missing = TRUE))
record("s8_rep006_ibad_empty_without_missing", seqibad(integr_e, with.missing = FALSE))
record("s8_rep006_ibad_empty_with_missing", seqibad(integr_e, with.missing = TRUE))
record("s8_rep006_ipos_integr_empty_with_missing", seqipos(integr_e, pos.states = "a", index = "integr", with.missing = TRUE))
# seqindic's Integr and Bad columns are seqintegr and seqibad (batch-19 review M-1); Integr of the state "a".
record("s8_rep006_seqindic_empty_without_missing",
  seqindic(integr_e, indic = c("integr", "bad"), ipos.args = list(pos.states = "a"), with.missing = FALSE))
record("s8_rep006_seqindic_empty_with_missing",
  seqindic(integr_e, indic = c("integr", "bad"), ipos.args = list(pos.states = "a"), with.missing = TRUE))

# ---- S3-12: seqipos volatility across an internal gap ----------------- #
gap_missing <- suppressMessages(seqdef(
  as.data.frame(rbind(c("A", "B", NA, NA, "C", "A")), stringsAsFactors = FALSE),
  alphabet = c("A", "B", "C")
))
gap_unlisted <- suppressMessages(seqdef(rbind(c("A", "B", "C", "C", "B", "A")), alphabet = c("A", "B", "C")))
record("s312_r_seqipos_missing_gap", seqipos(gap_missing, pos.states = "A", index = "volatility"))
record("s312_r_seqipos_unlisted_gap", seqipos(gap_unlisted, pos.states = "A", neg.states = "B", index = "volatility"))
record(
  "s312_r_seqindic_vpos_missing_gap",
  seqindic(gap_missing, indic = "vpos", ipos.args = list(pos.states = "A"))$Vpos
)
gap_free <- suppressMessages(seqdef(rbind(c("p", "n", "n", "p")), alphabet = c("n", "p")))
record("s312_seqivolatility_without_gap", seqivolatility(gap_free))

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
writeLines(c("# case\tkind (value, error or text)\tpayload\tR warnings", lines),
  file.path(out_dir, "audit20260922_s3b.txt"))
cat("AUDIT20260922_S3B_GOLDENS_DONE\n")
