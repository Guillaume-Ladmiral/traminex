# Golden generation for the 2026-09-22 audit, package s10b.
#
# seqMD_checkcost_diagnostics.txt -- the cost diagnostics TraMineR 2.2-14's
# seqMD emits for what = "cost" and "diss". seqMD runs checkcost() on every
# channel's costs before weighting them (seqMD.R:266-267): an explicit matrix
# or array, and equally the costs it generates with seqcost(..., check.costs =
# FALSE) "because costs are checked below" (seqMD.R:241-243). With what =
# "diss" the combined costs then pass through seqdist(), which runs its own
# checkcost (none for HAM). traminex's seqMD emitted none of the channel
# diagnostics (audit S10-11).
#
# One row per diagnostic, in the order R emits them: case, order, channel
# ("message" for message()/msg.warn(), "warning" for warning()), text with
# newlines written as the two characters \n (including the one message()
# appends). Only diagnostics are kept: warnings and messages containing
# " [!!]"; seqMD's progress messages are left out. A case with no diagnostic
# is recorded as "none".
#
# Consumed by tests/test_seqmd_cost_diagnostics.py.
#
# seqMD_cat_costs.txt -- two CAT-cost contracts of R's seqMD. (1) With
# what = "cost", seqMD ignores `method`: whether the costs are time-varying
# depends only on the dimensions of sm (seqMD.R:37; the DHD switch at :93
# applies to what = "diss" only), where traminex made method = "DHD"
# time-varying for what = "cost" too (audit S10-12; cases trate_method_*).
# (2) The combined cost of states i < j is computed once and mirrored
# (seqMD.R:307-316 and the time-varying loop), so an asymmetric channel
# matrix still gives symmetric CAT costs, where traminex kept the transposed
# channel costs below the diagonal (audit S10-13; cases asym_*). Rows: case,
# kind ("dim", "values" at %.17g column-major, "indel", "alphabet"), then
# tab-separated fields.
# Consumed by tests/test_seqmd_cost_method_ignored.py and
# tests/test_seqmd_asymmetric_channel_costs.py.
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
        if (grepl(" [!!]", text, fixed = TRUE)) {
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

q <- function(expr) suppressMessages(expr)
M <- function(v) matrix(v, 3, byrow = TRUE)
SYM <- M(c(0, 2, 2, 2, 0, 2, 2, 2, 0))
GT <- M(c(0, 5, 5, 5, 0, 5, 5, 5, 0))
TRI <- M(c(0, 10, 2, 10, 0, 2, 2, 2, 0))
ASY <- M(c(0, 3, 2, 1, 0, 2, 2, 2, 0))
ZAB <- M(c(0, 0, 2, 0, 0, 2, 2, 2, 0))
ZALL <- matrix(0, 3, 3)

# Three states each, one declared but unobserved (as make_goldens_seqMD.R).
chA <- q(seqdef(rbind(c("A", "B", "A", "B"), c("B", "A", "B", "A"), c("A", "A", "B", "B")),
  alphabet = c("A", "B", "U")))
chX <- q(seqdef(rbind(c("X", "Y", "X", "Y"), c("Y", "X", "Y", "X"), c("X", "X", "Y", "Y")),
  alphabet = c("X", "Y", "V")))
# chA with a missing value: its costs include the missing state "*".
chM <- q(seqdef(rbind(c("A", NA, "B", "A"), c("B", "A", "B", "A"), c("A", "B", "A", "B")),
  alphabet = c("A", "B", "U")))
# Corpora whose generated TRATE costs break the triangle inequality.
chT <- q(seqdef(rbind(c("A", "B", "A", "B"), c("B", "A", "B", "A"), c("C", "C", "A", "B"),
  c("A", "A", "C", "C")), alphabet = c("A", "B", "C")))
chS <- q(seqdef(rbind(c("X", "Y", "X", "Y"), c("Y", "X", "Y", "X"), c("X", "X", "Z", "Z"),
  c("Z", "Y", "X", "X")), alphabet = c("X", "Y", "Z")))

record("cost_explicit_zero_then_triangle",
  seqMD(list(chA, chX), what = "cost", sm = list(ZAB, TRI), indel = c(1, 6)))
record("cost_explicit_twice_indel_then_asymmetric",
  seqMD(list(chA, chX), what = "cost", sm = list(GT, ASY), indel = c(1, 2)))
record("cost_explicit_indel_vectors",
  seqMD(list(chA, chX), what = "cost", sm = list(SYM, SYM),
    indel = list(c(0.1, 5, 5), c(1, 1, 1))))
record("cost_explicit_all_zero_channel",
  seqMD(list(chA, chX), what = "cost", sm = list(ZALL, SYM), indel = c(1, 1)))
record("cost_explicit_auto_indel",
  seqMD(list(chA, chX), what = "cost", sm = list(GT, TRI)))
ZM <- matrix(2, 4, 4)
diag(ZM) <- 0
ZM[3, 4] <- 0
ZM[4, 3] <- 0
record("cost_explicit_missing_state_label",
  seqMD(list(chM, chX), what = "cost", sm = list(ZM, SYM), indel = c(1, 1),
    with.missing = c(TRUE, FALSE)))
record("cost_generated_trate",
  seqMD(list(chT, chS), what = "cost", sm = "TRATE"))
# what = "cost" ignores method = "DHD": these costs are stationary.
record("cost_method_dhd_generated_trate",
  seqMD(list(chM, chX), method = "DHD", what = "cost", sm = "TRATE",
    with.missing = c(TRUE, FALSE)))
A <- array(2, c(3, 3, 4))
for (t in 1:4) diag(A[, , t]) <- 0
A_zero <- A
A_zero[1, 2, 2] <- 0
A_zero[2, 1, 2] <- 0
A_asym <- A
A_asym[1, 2, 1] <- 3
A_asym[1, 2, 3] <- 3
record("cost_explicit_arrays",
  seqMD(list(chA, chX), method = "DHD", what = "cost", sm = list(A_zero, A_asym),
    indel = c(1, 1)))
record("diss_om_explicit_zero",
  seqMD(list(chA, chX), method = "OM", what = "diss", sm = list(ZAB, SYM), indel = c(1, 1)))
record("diss_ham_explicit_triangle",
  seqMD(list(chA, chX), method = "HAM", what = "diss", sm = list(TRI, SYM)))
record("diss_dhd_generated_trate",
  seqMD(list(chM, chX), method = "DHD", what = "diss", sm = "TRATE",
    with.missing = c(TRUE, FALSE)))
record("mdseq_checks_no_costs",
  seqMD(list(chA, chX), what = "MDseq", sm = list(ZAB, TRI)))

writeLines(c("# case\torder\tchannel\ttext (newlines written as \\n)", lines),
  file.path(out_dir, "seqMD_checkcost_diagnostics.txt"))

## seqMD_cat_costs.txt ----------------------------------------------------- ##
rows <- character()
add <- function(case, kind, fields) {
  stopifnot(length(fields) > 0)
  rows <<- c(rows, paste(c(case, kind, as.character(fields)), collapse = "\t"))
}
add_cost <- function(case, cost) {
  add(case, "dim", if (is.null(dim(cost))) length(cost) else dim(cost))
  add(case, "values", sprintf("%.17g", as.vector(cost)))
  add(case, "indel", if (is.null(attr(cost, "indel"))) "NULL" else
    sprintf("%.17g", as.vector(attr(cost, "indel"))))
  add(case, "alphabet", attr(cost, "alphabet"))
}
for (method in c("none", "DHD", "OM", "HAM")) {
  args <- list(channels = list(chM, chX), what = "cost", sm = "TRATE",
    with.missing = c(TRUE, FALSE))
  if (method != "none") args$method <- method
  add_cost(paste0("trate_method_", method), q(do.call(seqMD, args)))
}
quiet_all <- function(expr) suppressWarnings(suppressMessages(expr))
add_cost("asym_stationary", quiet_all(seqMD(list(chA, chX), what = "cost",
  sm = list(ASY, SYM), indel = c(1, 1))))
add_cost("asym_time_varying", quiet_all(seqMD(list(chA, chX), method = "DHD", what = "cost",
  sm = list(A_asym, A), indel = c(1, 1))))
d <- quiet_all(seqMD(list(chA, chX), method = "OM", what = "diss", sm = list(ASY, SYM),
  indel = c(1, 1)))
add("asym_diss_om", "dim", dim(as.matrix(d)))
add("asym_diss_om", "values", sprintf("%.17g", as.vector(as.matrix(d))))
writeLines(c("# case\tkind\tfields", rows), file.path(out_dir, "seqMD_cat_costs.txt"))
cat("AUDIT20260922_S10B_GOLDENS_DONE\n")
