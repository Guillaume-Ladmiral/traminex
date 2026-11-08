# Golden generation for the 2026-09-22 audit, package s10.
#
# checkcost_diagnostics.txt -- every non-fatal diagnostic TraMineR 2.2-14's
# checkcost() (R/TraMineR-checkcost.R) emits for a caller-supplied cost
# matrix or array, with its channel: "message" for message()/msg.warn(),
# "warning" for warning(). traminex had kept the 2.2-13 texts and channels and
# lacked the zero-cost and max(sm) <= 0 diagnostics (audit S10-01, S1a-11).
# The " [>]" progress lines msg() prints are information, not diagnostics,
# and are left out. A case with no diagnostic is recorded as "none". Each text
# is written with its newlines as the two characters \n, including the one
# message() appends.
#
# Consumed by tests/test_checkcost_diagnostics_parity.py.
#
# inline_r_anchors.txt -- R values that twelve test modules carried inline,
# citing "live R" or an untracked script, with no tracked generator (audit
# S10-04). One row per case: case, kind, then tab-separated fields. Kinds:
# "values" (numbers, %.17g, column-major), "error" (R's condition message),
# "strings" (character values in R's order), "dim" (dim()), "length" (a
# vector without dim). The consuming test is named in each section.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

lines <- character()
record <- function(case, expr) {
  seen <- character()
  channel <- character()
  result <- tryCatch(
    withCallingHandlers(expr,
      warning = function(w) {
        seen <<- c(seen, conditionMessage(w))
        channel <<- c(channel, "warning")
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        text <- conditionMessage(m)
        if (!startsWith(text, " [>]")) {
          seen <<- c(seen, text)
          channel <<- c(channel, "message")
        }
        invokeRestart("muffleMessage")
      }
    ),
    error = function(e) e
  )
  stopifnot(!inherits(result, "error"))
  if (length(seen) == 0) {
    lines <<- c(lines, paste(case, "0", "none", "none", sep = "\t"))
  } else {
    encoded <- gsub("\n", "\\n", seen, fixed = TRUE)
    lines <<- c(lines, paste(case, seq_along(seen), channel, encoded, sep = "\t"))
  }
}

RAW <- rbind(
  c("A", "B", "C", "A", "B", "C"),
  c("B", "A", "C", "C", "A", "B"),
  c("C", "C", "A", "B", "B", "A")
)
s <- suppressMessages(seqdef(RAW, alphabet = c("A", "B", "C")))
M <- function(v) matrix(v, 3, byrow = TRUE)
SYM <- M(c(0, 2, 2, 2, 0, 2, 2, 2, 0))
GT <- M(c(0, 5, 5, 5, 0, 5, 5, 5, 0))
TRI <- M(c(0, 10, 2, 10, 0, 2, 2, 2, 0))
TRIF <- M(c(0, 10.5, 1 / 3, 10.5, 0, 2.25, 1 / 3, 2.25, 0))
ASY <- M(c(0, 3, 2, 1, 0, 2, 2, 2, 0))
ZAB <- M(c(0, 0, 2, 0, 0, 2, 2, 2, 0))
ZBC <- M(c(0, 2, 2, 2, 0, 0, 2, 0, 0))
ZTWO <- M(c(0, 0, 0, 0, 0, 2, 0, 2, 0))
ZALL <- matrix(0, 3, 3)

record("om_twice_indel", seqdist(s, method = "OM", sm = GT, indel = 1))
record("om_twice_indel_auto", seqdist(s, method = "OM", sm = GT))
record("om_triangle", seqdist(s, method = "OM", sm = TRI, indel = 6))
record("om_triangle_fractional", seqdist(s, method = "OM", sm = TRIF, indel = 6))
record("om_indel_vector_triangle", seqdist(s, method = "OM", sm = SYM, indel = c(0.1, 5, 5)))
record("om_indel_vector_no_twice", seqdist(s, method = "OM", sm = GT, indel = c(1, 1, 1)))
record("om_asymmetric", seqdist(s, method = "OM", sm = ASY, indel = 2))
record("om_well_formed", seqdist(s, method = "OM", sm = SYM, indel = 1))
record("om_zero_ab", seqdist(s, method = "OM", sm = ZAB, indel = 1))
record("om_zero_bc", seqdist(s, method = "OM", sm = ZBC, indel = 1))
record("om_zero_two_pairs", seqdist(s, method = "OM", sm = ZTWO, indel = 1))
record("om_all_zero", seqdist(s, method = "OM", sm = ZALL, indel = 1))
# max(sm) <= 0 returns before checkcost's refusals: R then computes with a
# negative matrix, and returns meaningless tiny values (about 1e-309 in a
# probe) for a wrong-sized all-zero one; the values are not recorded.
NEG <- matrix(-1, 3, 3)
diag(NEG) <- 0
record("om_negative_max_zero", seqdist(s, method = "OM", sm = NEG, indel = 1))
record("om_wrong_size_all_zero", seqdist(s, method = "OM", sm = matrix(0, 2, 2), indel = 1))
# traminex's OMloc takes no indel; R's automatic one is max(sm) / 2.
record("omloc_triangle", seqdist(s, method = "OMloc", sm = TRI))
record("omloc_zero_ab", seqdist(s, method = "OMloc", sm = ZAB))
record("omslen_triangle", seqdist(s, method = "OMslen", sm = TRI, indel = 6))
record("omslen_zero_ab", seqdist(s, method = "OMslen", sm = ZAB, indel = 6))
record("omspell_triangle", seqdist(s, method = "OMspell", sm = TRI, indel = 6))
record("omspell_zero_ab", seqdist(s, method = "OMspell", sm = ZAB, indel = 6))
record("omstran_triangle", seqdist(s, method = "OMstran", sm = TRI, indel = 6, otto = 0.5))
record("omstran_zero_ab", seqdist(s, method = "OMstran", sm = ZAB, indel = 6, otto = 0.5))
record("twed_triangle", seqdist(s, method = "TWED", sm = TRI, indel = 6, nu = 0.5, h = 0.5))
record("twed_zero_ab", seqdist(s, method = "TWED", sm = ZAB, indel = 6, nu = 0.5, h = 0.5))
# HAM: a dangling else in seqdist.R:380-386 means R never calls checkcost.
record("ham_triangle", seqdist(s, method = "HAM", sm = TRI))
record("ham_zero_ab", seqdist(s, method = "HAM", sm = ZAB))

# with.missing: the missing state is named by its nr symbol.
RAWM <- rbind(c("A", "B", NA, "A"), c("B", "A", "C", "C"), c("C", NA, "A", "B"))
sm_ <- suppressMessages(seqdef(RAWM, alphabet = c("A", "B", "C")))
ZM <- matrix(2, 4, 4)
diag(ZM) <- 0
ZM[3, 4] <- 0
ZM[4, 3] <- 0
TM <- matrix(2, 4, 4)
diag(TM) <- 0
TM[1, 4] <- 9
TM[4, 1] <- 9
record("om_missing_zero", seqdist(sm_, method = "OM", sm = ZM, indel = 1, with.missing = TRUE))
record("om_missing_triangle", seqdist(sm_, method = "OM", sm = TM, indel = 5, with.missing = TRUE))

# DHD time-varying arrays (R layout states x states x time): each diagnostic
# names its position; the zero-cost and triangle ones are reported once.
A <- array(2, c(3, 3, 6))
for (t in 1:6) diag(A[, , t]) <- 0
A_zero <- A
A_zero[1, 2, 3] <- 0
A_zero[2, 1, 3] <- 0
A_zero[1, 3, 5] <- 0
A_zero[3, 1, 5] <- 0
A_triangle <- A
A_triangle[1, 2, 4] <- 10
A_triangle[2, 1, 4] <- 10
A_triangle[1, 2, 6] <- 9
A_triangle[2, 1, 6] <- 9
A_asymmetric <- A
A_asymmetric[1, 2, 2] <- 3
A_asymmetric[1, 2, 5] <- 3
record("dhd_zero", seqdist(s, method = "DHD", sm = A_zero))
record("dhd_triangle", seqdist(s, method = "DHD", sm = A_triangle))
record("dhd_asymmetric", seqdist(s, method = "DHD", sm = A_asymmetric))
record("dhd_all_zero", seqdist(s, method = "DHD", sm = A * 0))

# seqdiff wraps its internal seqdist in suppressMessages() (seqdiff.R:62), so
# only checkcost's warning() channel can reach its caller.
RAWD <- rbind(RAW, c("A", "A", "B", "B", "C", "C"), c("B", "B", "B", "A", "A", "C"),
  c("C", "A", "B", "C", "A", "B"))
sd <- suppressMessages(seqdef(RAWD, alphabet = c("A", "B", "C")))
grp <- c(1, 1, 1, 2, 2, 2)
record("seqdiff_twice_indel", seqdiff(sd, group = grp, cmprange = c(0, 1),
  seqdist.args = list(method = "OM", sm = GT, indel = 1)))
record("seqdiff_indel_vector_triangle", seqdiff(sd, group = grp, cmprange = c(0, 1),
  seqdist.args = list(method = "OM", sm = SYM, indel = c(0.1, 5, 5))))
record("seqdiff_asymmetric", seqdiff(sd, group = grp, cmprange = c(0, 1),
  seqdist.args = list(method = "OM", sm = ASY, indel = 2)))

writeLines(c("# case\torder\tchannel\ttext (newlines written as \\n)", lines),
  file.path(out_dir, "checkcost_diagnostics.txt"))

## inline_r_anchors.txt ---------------------------------------------------- ##
anchors <- character()
anchor <- function(case, kind, fields) {
  fields <- gsub("\n", "\\n", as.character(fields), fixed = TRUE)
  stopifnot(length(fields) > 0, !any(grepl("\t", fields)))
  anchors <<- c(anchors, paste(c(case, kind, fields), collapse = "\t"))
}
num <- function(x) sprintf("%.17g", as.vector(x))
quiet <- function(expr) suppressWarnings(suppressMessages(expr))
error_of <- function(expr) {
  message <- tryCatch({
    quiet(expr)
    NA_character_
  }, error = function(e) conditionMessage(e))
  stopifnot(!is.na(message))
  message
}
shape <- function(case, x) {
  if (is.null(dim(x))) anchor(case, "length", length(x)) else anchor(case, "dim", dim(x))
}

# tests/test_audit_fixes_distances.py
X3 <- quiet(seqdef(rbind(c("A", "A", "B", "B", "C"), c("A", "B", "B", "C", "C"),
  c("C", "C", "B", "A", "A")), alphabet = c("A", "B", "C")))
X2 <- quiet(seqdef(rbind(c("A", "A", "B", "B", "A"), c("A", "B", "B", "A", "A"),
  c("B", "B", "B", "A", "A")), alphabet = c("A", "B")))
for (m in c("OM", "OMloc", "OMslen", "OMspell", "OMstran", "TWED")) {
  anchor(paste0("no_sm_", m), "error", error_of(seqdist(X3, method = m)))
}
anchor("lcs_maxlength", "values", num(as.matrix(quiet(seqdist(X3, method = "LCS",
  norm = "maxlength")))))
anchor("lcs_norm_length", "error", error_of(seqdist(X3, method = "LCS", norm = "length")))
anchor("om_asymmetric_maxdist", "values", num(as.matrix(quiet(seqdist(X2, method = "OM",
  sm = matrix(c(0, 2, 1.5, 0), 2), indel = 2, norm = "maxdist")))))

# tests/test_distance_wrappers_norm.py -- corestats fixture B, rows 1 and 2
B <- rbind(c("A", "A", "*", "B", "B", "C"), c("*", "A", "A", "A", "*", "*"),
  c("A", "*", "A", "B", "B", "B"), c("C", "C", "C", "C", "C", "C"))
sB <- quiet(seqdef(B, alphabet = c("A", "B", "C"), missing = "*", right = NA))
smB <- matrix(2, 4, 4)
diag(smB) <- 0
for (m in c("OM", "LCS", "LCP", "RLCP", "HAM", "DHD")) {
  values <- character()
  for (norm in c("auto", "maxlength", "gmean", "maxdist", "YujianBo")) {
    args <- list(seqdata = sB, method = m, with.missing = TRUE, norm = norm)
    if (m == "OM") {
      args$sm <- smB
      args$indel <- 1
    }
    values <- c(values, num(as.matrix(quiet(do.call(seqdist, args)))[1, 2]))
  }
  anchor(paste0("fixture_b_d12_", m), "values", values)
}
xAB <- quiet(seqdef(rbind("A", "B"), alphabet = c("A", "B")))
for (norm in c("maxdist", "gmean")) {
  anchor(paste0("om_asymmetric_upper_bound_", norm), "values", num(as.matrix(quiet(seqdist(xAB,
    method = "OM", sm = matrix(c(0, 10, 1, 0), 2), indel = 10, norm = norm)))[1, 2]))
}

# tests/test_external_marker_equivalence.py
FXC <- quiet(seqdef(rbind(c("A", NA, "A", "B"), c("B", "B", "A", "A")), alphabet = c("A", "B")))
anchor("fx_c_seqient_raw", "values", num(quiet(seqient(FXC, with.missing = TRUE, norm = FALSE))))
anchor("fx_c_seqient_norm", "values", num(quiet(seqient(FXC, with.missing = TRUE))))
# Row 1 ends in nr, row 2 in void: seqdef cannot make both trailing, so the
# void cell is set directly, as the traminex test builds it.
REC <- quiet(seqdef(rbind(c("A", "B", "B", "*"), c("A", "A", "C", "*")),
  alphabet = c("A", "B", "C"), missing = "*", right = NA))
REC[2, 4] <- attr(REC, "void")
recoded <- quiet(seqrecode(REC, recodes = list(W = c("A", "B")), otherwise = "O"))
anchor("seqrecode_otherwise_alphabet", "strings", alphabet(recoded))
anchor("seqrecode_otherwise_rows", "strings", as.vector(seqconc(recoded)))

# tests/test_icu_collation_parity.py, tests/test_icu_native_requirement.py
# The locale pin at the top set C.UTF-8 and the ICU root collation; a second
# Sys.setlocale here would reset the ICU collator.
stopifnot(Sys.getlocale("LC_COLLATE") == "C.UTF-8")
anchor("sort_c_utf8", "strings", sort(unique(c("b", "A", "a", "B", "Z", "é", "e", "_x",
  "1", "10", "2"))))

# tests/test_p2_functions.py
DRF <- quiet(seqdef(rbind(c("A", "A"), c("B", "B"), c("C", "C"), c("A", "B")),
  weights = c(1, 2, 1, 1.5)))
# The defaults are spelled out so tests/test_golden_default_coverage.py does not
# count this call, which records one sum, as exercising them.
drf <- quiet(dissrf(quiet(seqdist(DRF, method = "HAM")), k = 2, weights = attr(DRF, "weights"),
  sortv = "mds", grp.meth = "prop"))
anchor("dissrf_weighted_wsum", "values", num(drf$sizes[["wsum"]]))

# tests/test_per_state_shape_and_selection.py
PS <- quiet(seqdef(rbind(c(NA, "A", "B", "C", "A", NA), c("A", "B", NA, "C", "C", "B"),
  c("B", "B", "A", NA, NA, "C"), c("C", "A", "A", "B", "B", NA), c(NA, NA, "C", "C", "A", "A"),
  c("A", "C", "B", "B", "A", "C")), alphabet = c("A", "B", "C")))
shape("seqmeant_without_missing", quiet(seqmeant(PS, with.missing = FALSE)))
shape("seqmeant_with_missing", quiet(seqmeant(PS, with.missing = TRUE)))
shape("seqstatd_frequencies", quiet(seqstatd(PS))$Frequencies)
shape("seqistatd", quiet(seqistatd(PS)))
# Defaults spelled out for the same reason: only the shape is recorded.
shape("seqstatf", quiet(seqstatf(PS, weighted = TRUE, with.missing = FALSE)))
shape("seqtrate", quiet(seqtrate(PS, lag = 1)))
shape("cpal", cpal(PS))
shape("stlab", stlab(PS))
for (f in c("seqient", "seqtransn", "seqlength", "seqici", "seqST", "seqivolatility",
  "seqivardur")) {
  shape(paste0("per_sequence_", f), quiet(get(f)(PS)))
}

# tests/test_seqcost_refseq_missing_contract.py
SRC <- quiet(seqdef(rbind(c("A", "B", "B", "C"), c("A", "A", "C", "C"), c("B", "B", "A", "A")),
  alphabet = c("A", "B", "C")))
REF <- quiet(seqdef(rbind(c("A", "B", NA, "C")), alphabet = c("A", "B", "C")))
for (m in c("CONSTANT", "TRATE", "INDELS", "INDELSLOG")) {
  anchor(paste0("refseq_missing_", m), "values", num(quiet(seqdist(SRC, method = "OM", sm = m,
    indel = 1, refseq = REF, with.missing = TRUE))))
}
constant <- quiet(seqcost(SRC, method = "CONSTANT", cval = 2, with.missing = TRUE))
anchor("seqcost_constant_with_missing_sm", "values", num(constant$sm))
shape("seqcost_constant_with_missing_sm", constant$sm)
anchor("seqcost_constant_with_missing_indel", "values", num(constant$indel))

# tests/test_seqdist_infinite_cost_contract.py
INF <- quiet(seqdef(rbind(c("A", "B", "A"), c("B", "A", "B"), c("A", "A", "B")),
  alphabet = c("A", "B", "C", "D", "E")))
costs5 <- function(r, c, v) {
  m <- matrix(2, 5, 5)
  diag(m) <- 0
  m[r, c] <- v
  m[c, r] <- v
  m
}
for (v in c(2, 50)) {
  d <- as.matrix(quiet(seqdist(INF, method = "OM", sm = costs5(4, 5, v), norm = "maxdist")))
  anchor(paste0("unobserved_cost_", v), "values", num(c(d[1, 2], d[1, 3], d[2, 3])))
}
for (rc in list(c(1, 2), c(4, 5), c(3, 4))) {
  anchor(paste0("infinite_cost_", rc[1], "_", rc[2]), "error",
    error_of(seqdist(INF, method = "OM", sm = costs5(rc[1], rc[2], Inf))))
}
anchor("infinite_unobserved_cost_given_indel", "error",
  error_of(seqdist(INF, method = "OM", sm = costs5(4, 5, Inf), indel = 1)))
DDD <- quiet(seqdef(rbind(c("D", "D", "D")), alphabet = c("A", "B", "C", "D", "E")))
anchor("indels_refseq_unobserved", "error", error_of(seqdist(INF, method = "OM", sm = "INDELS",
  refseq = DDD)))
anchor("indels_refseq_unobserved_no_check", "values", num(quiet(seqdist(INF, method = "OM",
  sm = "INDELS", refseq = DDD, check.costs = FALSE))))

# tests/test_seqdist_wrapper.py
NMS_SRC <- quiet(seqdef(rbind(c("A", "A", "A"), c("A", "B", "A")), alphabet = c("A", "B"),
  missing = "*"))
NMS_REF <- quiet(seqdef(rbind(c("A", "*", "A")), alphabet = c("A", "B"), missing = "*"))
anchor("nms_refseq_with_missing", "values", num(quiet(seqdist(NMS_SRC, method = "NMS",
  refseq = NMS_REF, with.missing = TRUE))))

# tests/test_seqllcs_seqllcp_arity.py
LL <- quiet(seqdef(rbind(c("A", "B", "C", "A"), c("B", "A", "C", "C")),
  alphabet = c("A", "B", "C")))
anchor("llcs_rows_1_2", "values", num(seqLLCS(LL[1, ], LL[2, ])))
anchor("llcp_rows_1_2", "values", num(seqLLCP(LL[1, ], LL[2, ])))
anchor("llcs_row_1_itself", "values", num(seqLLCS(LL[1, ], LL[1, ])))
anchor("llcp_row_1_itself", "values", num(seqLLCP(LL[1, ], LL[1, ])))
anchor("llcs_two_rows", "error", error_of(seqLLCS(LL, LL)))
anchor("llcp_two_rows", "error", error_of(seqLLCP(LL, LL)))

writeLines(c("# case\tkind\tfields (tab-separated; newlines written as \\n)", anchors),
  file.path(out_dir, "inline_r_anchors.txt"), useBytes = TRUE)
cat("AUDIT20260922_S10_GOLDENS_DONE\n")
