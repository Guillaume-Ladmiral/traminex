# Golden generation: statistics edge cases from the 2026-09-22 audit (package s3a).
#
# Each record pins one TraMineR 2.2-14 result that traminex got wrong:
#
# * S3-01 seqcomp/seqmpos compare the two objects' states as factor labels
#   (Ops.factor), so objects with different alphabets or alphabet orders are
#   compared by state, not by integer code.
# * S3-02 seqidegrad (and seqprecarity/seqinsecurity/seqindic through it)
#   accepts a logical penalized, validates its options, and implements the
#   "ONE" and "+" methods.
# * S3-03 seqidegrad(method = "TRATE") does not bridge missing cells.
# * S3-04 seqstatd(with.missing = TRUE) counts the missing state in the
#   entropy maximum and the palette even when no cell is missing.
# * S3-05 seqalign(with.missing = TRUE) sizes sm with the nr state even when
#   no cell is missing.
# * S3-06 seqtransn(pweight = TRUE, with.missing = TRUE) is NA for a sequence
#   with a transition to or from the missing state.
# * S3-07 seqmodst has no modal state (NA) where no cell is valid.
# * S3-08 seqmeant(prop = TRUE) is NaN when a sequence has no valid cell.
# * S3-09 seqstatf(with.missing = TRUE) has an nr row even when no cell is
#   missing, and its Percent is NaN when the total weight is zero.
#
# Consumed by the tests named in each section.  Records are
# "case<TAB>kind<TAB>values<TAB>warnings": kind is "value" (values in R's
# column-major order, numbers as %.17g, logicals as TRUE/FALSE) or "error"
# (the R error message).
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

lines <- character()
render <- function(result) {
  values <- as.vector(as.matrix(result))
  if (is.logical(values)) {
    return(paste(ifelse(is.na(values), "NA", ifelse(values, "TRUE", "FALSE")), collapse = " "))
  }
  if (is.numeric(values)) {
    return(paste(sprintf("%.17g", values), collapse = " "))
  }
  paste(ifelse(is.na(values), "NA", as.character(values)), collapse = " ")
}
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
    lines <<- c(lines, paste(case, "value", render(result), warned, sep = "\t"))
  }
}
sts <- function(rows, ...) {
  data <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  suppressMessages(seqdef(data, ...))
}

# ---- S3-01: seqcomp / seqmpos across objects with different codebooks ----
# Consumed by tests/test_seqcomp_seqmpos_state_labels.py.
x <- sts(list(c("A", "B", "B")))
y <- sts(list(c("B", "A", "A")))
z <- sts(list(c("A", "B", "B")), alphabet = c("B", "A"))
w <- sts(list(c("A", "C", "C")))
xs <- sts(list(c("A", "B", "B")), alphabet = c("A", "B", "C"))
m <- sts(list(c("A", NA, "B")))
n <- sts(list(c("A", "B", "B")))
a <- sts(list(c("A", "B", "C"), c("C", "C", "A")), alphabet = c("A", "B", "C"))
b <- sts(list(c("A", "C", "C"), c("C", "B", "A")), alphabet = c("C", "B", "A"))
v1 <- sts(list(c("A", "B", "B"), c("A", "B", NA)))
v2 <- sts(list(c("B", "B", "B"), c("A", "B", NA)), alphabet = c("B", "A"))
record("s3_01_seqcomp_x_y", seqcomp(x, y))
record("s3_01_seqcomp_x_z", seqcomp(x, z))
record("s3_01_seqcomp_x_w", seqcomp(x, w))
record("s3_01_seqcomp_x_xs", seqcomp(x, xs))
record("s3_01_seqcomp_m_n", seqcomp(m, n))
record("s3_01_seqcomp_m_m", seqcomp(m, m))
record("s3_01_seqcomp_v1r2_v2r2", seqcomp(v1[2, ], v2[2, ]))
record("s3_01_seqmpos_x_y", seqmpos(x, y))
record("s3_01_seqmpos_x_z", seqmpos(x, z))
record("s3_01_seqmpos_x_w", seqmpos(x, w))
record("s3_01_seqmpos_x_xs", seqmpos(x, xs))
record("s3_01_seqmpos_m_n", seqmpos(m, n))
record("s3_01_seqmpos_m_n_wm", seqmpos(m, n, with.missing = TRUE))
record("s3_01_seqmpos_m_m_wm", seqmpos(m, m, with.missing = TRUE))
record("s3_01_seqmpos_a_b", seqmpos(a, b))
record("s3_01_seqmpos_v1_v2", seqmpos(v1, v2))
# A real state spelled like the other object's nr text: each object's missing
# cells are identified by its own nr attribute only.
st <- sts(list(c("A", "*", "B")), alphabet = c("A", "*", "B"), nr = "?")
record("s3_01_seqmpos_st_m", seqmpos(st, m))
record("s3_01_seqmpos_m_st", seqmpos(m, st))
record("s3_01_seqmpos_st_m_wm", seqmpos(st, m, with.missing = TRUE))
# The same for a real state spelled like the other object's void text.
sv <- sts(list(c("A", "B", "%")), alphabet = c("A", "B", "%"), void = "#")
vr <- sts(list(c("A", "B", NA)))
record("s3_01_seqmpos_sv_vr", seqmpos(sv, vr))
record("s3_01_seqmpos_vr_sv", seqmpos(vr, sv))

# ---- S3-02: the degradation family's method / penalized options ----
# seqidegrad.R:74-119 maps a logical penalized, validates method, weight.type,
# border.effect and penalized, gives every transition weight 1 for "ONE" and
# scales by the mean transition weight per spell for the "+" methods.
# Consumed by tests/test_seqidegrad_method_options.py.
dg <- sts(list(
  c("A", "A", "B", "B", "C", "C"), c("A", "B", "B", "C", "C", "C"),
  c("B", "B", "B", "A", "A", "C"), c("C", "C", "A", "A", "B", "B"),
  c("A", "A", "A", "A", "A", "A"), c("B", "C", "B", "C", "B", "C")
), alphabet = c("A", "B", "C"), weights = c(1, 2, 0.5, 1.5, 3, 1))
for (meth in c("RANK+", "TRATE+", "TRATEDSS+", "FREQ+", "ONE")) {
  record(paste0("s3_02_degrad_", meth), seqidegrad(dg, method = meth))
  record(paste0("s3_02_degrad_", meth, "_nointegr"), seqidegrad(dg, method = meth, pow = FALSE))
}
record("s3_02_degrad_penalized_true", seqidegrad(dg, penalized = TRUE))
record("s3_02_degrad_one_neg", seqidegrad(dg, method = "ONE", penalized = "NEG"))
record("s3_02_degrad_one_partial", seqidegrad(dg, method = "ONE", state.order = c("A", "B")))
record("s3_02_degrad_one_equiv", seqidegrad(dg, method = "ONE", state.equiv = list(c("B", "C"))))
record("s3_02_degrad_rank_bogus_weight", seqidegrad(dg, method = "RANK", weight.type = "BOGUS"))
record("s3_02_degrad_rank_border_half", seqidegrad(dg, method = "RANK", border.effect = 0.5))
record("s3_02_degrad_one_bogus_weight", seqidegrad(dg, method = "ONE", weight.type = "BOGUS"))
record("s3_02_degrad_bogus_method", seqidegrad(dg, method = "BOGUS"))
record("s3_02_degrad_lowercase_method", seqidegrad(dg, method = "rank"))
record("s3_02_degrad_bogus_penalized", seqidegrad(dg, penalized = "BOGUS"))
record("s3_02_degrad_trateplus_bogus_weight", seqidegrad(dg, method = "TRATE+", weight.type = "BOGUS"))
record("s3_02_degrad_freq_border_one", seqidegrad(dg, method = "FREQ", border.effect = 1))
record("s3_02_prec_one", seqprecarity(dg, method = "ONE"))
record("s3_02_prec_tratedssplus", seqprecarity(dg, method = "TRATEDSS+"))
record("s3_02_prec_penalized_true", seqprecarity(dg, penalized = TRUE))
record("s3_02_prec_penalized_false", seqprecarity(dg, penalized = FALSE))
record("s3_02_prec_no_bogus_method", seqprecarity(dg, penalized = "NO", method = "BOGUS"))
record("s3_02_insec_one", seqinsecurity(dg, method = "ONE"))
record("s3_02_insec_rankplus", seqinsecurity(dg, method = "RANK+"))
record("s3_02_indic_one", seqindic(dg, indic = c("degrad", "prec", "insec"),
  prec.args = list(method = "ONE")))
record("s3_02_indic_trateplus_logical", seqindic(dg, indic = c("degrad", "prec", "insec"),
  prec.args = list(method = "TRATE+", penalized = TRUE)))
dgm <- sts(list(c("A", "B", NA, "C"), c(NA, NA, NA, NA), c("C", "C", "A", "B")),
  alphabet = c("A", "B", "C"), right = NA)
record("s3_02_degrad_missing_RANK+", seqidegrad(dgm, method = "RANK+"))
record("s3_02_degrad_missing_RANK+_wm", seqidegrad(dgm, method = "RANK+", with.missing = TRUE))

# ---- S3-03: TRATE counts only adjacent pairs, as seqtrate does ----
# seqidegrad.R:233-234 calls seqtrate(seqdata), which pairs column t with
# column t + 1 only, so a missing cell breaks a transition instead of being
# bridged.  Consumed by tests/test_seqidegrad_trate_missing_cells.py.
tb <- sts(list(c("A", NA, "A", "B", "C"), c("C", "B", "B", "A", "A"),
  c("B", "C", "C", "A", "B")), alphabet = c("A", "B", "C"))
tbw <- sts(list(c("A", NA, "A", "B", "C"), c("C", "B", "B", "A", "A"),
  c("B", "C", "C", "A", "B")), alphabet = c("A", "B", "C"), weights = c(1, 2, 0.5))
record("s3_03_degrad_trate", seqidegrad(tb, method = "TRATE"))
record("s3_03_degrad_trate_nointegr", seqidegrad(tb, method = "TRATE", pow = FALSE))
record("s3_03_degrad_trate_weighted", seqidegrad(tbw, method = "TRATE"))
record("s3_03_degrad_trateplus", seqidegrad(tb, method = "TRATE+"))
record("s3_03_degrad_trate_inv", seqidegrad(tb, method = "TRATE", weight.type = "INV"))
record("s3_03_degrad_trate_wm", seqidegrad(tb, method = "TRATE", with.missing = TRUE))
record("s3_03_prec_trate", seqprecarity(tb, method = "TRATE"))
record("s3_03_insec_trate", seqinsecurity(tb, method = "TRATE"))

# ---- S3-04: seqstatd(with.missing = TRUE) on data without missing cells ----
# seqstatd.R:12-17 always appends nr (and missing.color) to statl under
# with.missing, so the entropy maximum is log(K + 1) even when no cell is
# missing.  Consumed by tests/test_seqstatd_missing_state_entropy.py.
sf <- sts(list(c("A", "A", "B"), c("B", "C", "C")), alphabet = c("A", "B", "C"))
sfw <- sts(list(c("A", "A", "B"), c("B", "C", "C"), c("C", "A", "A")),
  alphabet = c("A", "B", "C"), weights = c(1, 2.5, 0.5))
sk1 <- sts(list(c("A", "A", "A"), c("A", "A", "A")), alphabet = "A")
record("s3_04_statd_entropy_wm", seqstatd(sf, with.missing = TRUE)$Entropy)
record("s3_04_statd_entropy_wm_raw", seqstatd(sf, with.missing = TRUE, norm = FALSE)$Entropy)
record("s3_04_statd_entropy_wm_weighted", seqstatd(sfw, with.missing = TRUE)$Entropy)
record("s3_04_statd_entropy_k1_wm", seqstatd(sk1, with.missing = TRUE)$Entropy)
record("s3_04_statd_freq_wm", seqstatd(sf, with.missing = TRUE)$Frequencies)
record("s3_04_statd_valid_wm", seqstatd(sf, with.missing = TRUE)$ValidStates)
record("s3_04_statd_cpal_wm", attr(seqstatd(sf, with.missing = TRUE), "cpal"))
record("s3_04_statd_cpal_k1_wm", attr(seqstatd(sk1, with.missing = TRUE), "cpal"))

# ---- S3-05: seqalign(with.missing = TRUE) on data without missing cells ----
# seqalign.R:9-14 names sm by alphabet(seqdata) plus nr under with.missing,
# so it takes a (K + 1) x (K + 1) matrix (the size seqsubm(with.missing =
# TRUE) builds) and refuses a K x K one even when no cell is missing.
# Consumed by tests/test_seqalign_missing_state_costs.py.
sa <- sts(list(c("A", "B", "C"), c("A", "C", "C"), c("C", "A", "B")),
  alphabet = c("A", "B", "C"))
sa_const <- suppressMessages(seqsubm(sa, method = "CONSTANT", cval = 2, with.missing = TRUE))
sa_var <- matrix(c(0, 1, 3, 5, 1, 0, 2, 5, 3, 2, 0, 5, 5, 5, 5, 0), nrow = 4)
record("s3_05_align_const_sm", sa_const)
for (case in list(list("const_12", c(1, 2), 1, sa_const), list("var_13", c(1, 3), 4, sa_var))) {
  al <- seqalign(sa, case[[2]], indel = case[[3]], sm = case[[4]], with.missing = TRUE)
  record(paste0("s3_05_align_", case[[1]], "_operation"), al$operation)
  record(paste0("s3_05_align_", case[[1]], "_seq1"), al$seq1)
  record(paste0("s3_05_align_", case[[1]], "_seq2"), al$seq2)
  record(paste0("s3_05_align_", case[[1]], "_cost"), al$cost)
  record(paste0("s3_05_align_", case[[1]], "_costmatrix"), al$costmatrix)
}
record("s3_05_align_k_by_k_sm", seqalign(sa, c(1, 2), indel = 1,
  sm = seqsubm(sa, method = "CONSTANT", cval = 2), with.missing = TRUE))

# ---- S3-06: seqtransn(pweight = TRUE, with.missing = TRUE) ----
# seqtransn.R:15-27 weighs DSS transitions with seqtrate(seqdata), which has
# no missing row or column, and indexes it with seqasnum(dss) + 1, which is NA
# for a missing cell; a transition touching the missing state therefore makes
# the sequence's count NA.  Consumed by tests/test_seqtransn_pweight_missing.py.
tm <- sts(list(c("A", NA, "B", "B"), c(NA, NA, NA, NA), c("B", NA, "C", "C"),
  c("A", "B", "B", "C"), c("C", "C", NA, NA)), alphabet = c("A", "B", "C"), right = NA)
record("s3_06_transn_pweight_wm", seqtransn(tm, with.missing = TRUE, pweight = TRUE))
record("s3_06_transn_pweight_wm_norm", seqtransn(tm, with.missing = TRUE, pweight = TRUE, norm = TRUE))
record("s3_06_transn_pweight", seqtransn(tm, pweight = TRUE))
record("s3_06_transn_wm", seqtransn(tm, with.missing = TRUE))

# ---- S3-07: seqmodst at a position without any valid state ----
# seqmodst.R:28-32 takes which(freq == max(freq))[1], which is NA on an
# all-NaN frequency column (every cell excluded, or zero total weight), so the
# modal cell and its frequency are NA and the modal sequence matches nothing.
# Consumed by tests/test_seqmodst_no_valid_state.py.
mo <- sts(list(c("A", NA, "B", "B"), c(NA, NA, NA, NA), c("B", NA, "C", "C")),
  alphabet = c("A", "B", "C"), right = NA)
mz <- sts(list(c("A", "B", "A", "B"), c("B", "B", "A", "A")), alphabet = c("A", "B"),
  weights = c(0, 0))
mr <- sts(list(c("A", "B", NA), c("A", NA, NA)), alphabet = c("A", "B"))
for (case in list(list("missing", mo, TRUE, FALSE), list("missing_wm", mo, TRUE, TRUE),
  list("zero_weights", mz, TRUE, FALSE), list("zero_weights_unweighted", mz, FALSE, FALSE),
  list("ragged", mr, TRUE, FALSE))) {
  modal <- seqmodst(case[[2]], weighted = case[[3]], with.missing = case[[4]])
  record(paste0("s3_07_modst_", case[[1]], "_states"), as.matrix(modal))
  record(paste0("s3_07_modst_", case[[1]], "_freq"), attr(modal, "Frequencies"))
  record(paste0("s3_07_modst_", case[[1]], "_occurrences"), attr(modal, "Occurrences"))
}

# ---- S3-08: seqmeant(prop = TRUE) with a sequence without valid cells ----
# seqmeant.R:10-26 averages seqistatd(prop = TRUE), whose row for a sequence
# without any valid cell is 0/0 = NaN in every state, so every mean is NaN.
# Consumed by tests/test_seqmeant_prop_empty_sequence.py.
mm <- sts(list(c("A", NA, "B", "B"), c(NA, NA, NA, NA), c("B", NA, "C", "C")),
  alphabet = c("A", "B", "C"), right = NA)
mmw <- sts(list(c("A", NA, "B", "B"), c(NA, NA, NA, NA), c("B", NA, "C", "C")),
  alphabet = c("A", "B", "C"), right = NA, weights = c(2, 0, 1))
record("s3_08_meant_prop", seqmeant(mm, prop = TRUE))
record("s3_08_meant_prop_serr", seqmeant(mm, prop = TRUE, serr = TRUE))
record("s3_08_meant_prop_zero_weight", seqmeant(mmw, prop = TRUE))
record("s3_08_meant_prop_wm", seqmeant(mm, prop = TRUE, with.missing = TRUE))
record("s3_08_meant_counts", seqmeant(mm))

# ---- S3-09: seqstatf with the missing state absent, and zero total weight ----
# seqstatf.R:10-22 sums seqistatd(with.missing), which has an nr column
# whenever with.missing = TRUE, and computes Percent as Freq / sum(Freq) * 100,
# which is 0/0 = NaN when the total weight is zero.
# Consumed by tests/test_seqstatf_missing_row_and_zero_weight.py.
tf <- sts(list(c("A", "A", "B"), c("B", "C", "C")), alphabet = c("A", "B", "C"))
tz <- sts(list(c("A", "B"), c("B", "B")), alphabet = c("A", "B"), weights = c(0, 0))
for (case in list(list("wm", tf, TRUE, TRUE), list("zero_weight", tz, TRUE, FALSE),
  list("zero_weight_wm", tz, TRUE, TRUE), list("zero_weight_unweighted", tz, FALSE, FALSE))) {
  table <- seqstatf(case[[2]], weighted = case[[3]], with.missing = case[[4]])
  record(paste0("s3_09_statf_", case[[1]]), as.matrix(table))
  record(paste0("s3_09_statf_", case[[1]], "_rownames"), rownames(table))
}

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
writeLines(c("# case\tkind\tR result (column-major) or error\tR warnings", lines),
  file.path(out_dir, "audit20260922_s3a.txt"))
cat("AUDIT20260922_S3A_GOLDENS_DONE\n")
