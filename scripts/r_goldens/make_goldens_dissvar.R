# Golden generation: dissvar (Phase 2.1 / DUPLICATES A5).
#
# Anchors traminex's R-faithful weighted discrepancy `dissvar` to
# TraMineR 2.2-14 `dissvar(diss, weights=, squared=)`:
#   - unweighted        == sum(diss) / (2*n^2)
#   - weighted          == tmrWeightedInertiaDist(..., var=TRUE)
#   - weighted, squared == same on diss^2
#
# Inputs come from make_goldens_phase2_inputs.R (shared Phase 2 fixture).
# Requires R + TraMineR 2.2-14; see scripts/r_goldens/README.md.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
source(file.path("scripts", "r_goldens", "make_goldens_phase2_inputs.R"))

inp <- phase2_inputs()
D <- inp$D

v_unw <- dissvar(D)
v_w <- dissvar(D, weights = inp$wreal)
v_w_sq <- dissvar(D, weights = inp$wreal, squared = TRUE)
v_unw_sq <- dissvar(D, squared = TRUE)
v_wint <- dissvar(D, weights = inp$wint)

write_phase2_csv(
  "dissvar_phase2", matrix(c(v_unw, v_w, v_w_sq, v_unw_sq, v_wint), ncol = 1),
  paste0(
    "# R: TraMineR::dissvar on phase2_input_D; rows = [dissvar(D),",
    " dissvar(D,w=wreal), dissvar(D,w=wreal,squared=TRUE), dissvar(D,squared=TRUE),",
    " dissvar(D,w=wint)]\n"
  )
)
## ---- cycle-6 F-003: non-finite and non-positive total case weights ----
# `tmrWeightedInertiaDist` guards BOTH of its divisions with `if (totweights>0)`
# and otherwise returns the *undivided* sum
# (SOURCES/TraMineR/src/weightedinertia.cpp:34-38), so a NaN, Inf, all-negative
# or zero-sum weight vector is never mapped to a finite 0. R does not validate
# weights here and neither does traminex; only the value must agree.
#
# Rows are the weight vectors in WEIGHT_EDGE_CASES order (see
# tests/test_diss_weight_edge_cases.py); columns are
#   [dissvar(D,w), dissvar(D,w,squared=TRUE), tmrWeightedInertiaDist(var=FALSE)].
# The third column calls the registered C symbol directly: traminex's
# `weighted_inertia` is a direct port of that kernel at var=FALSE, and no
# exported R function surfaces it on its own.
Dw <- matrix(c(0, 1, 2, 1, 0, 3, 2, 3, 0), nrow = 3, byrow = TRUE)
weight_edge_cases <- list(
  ones = c(1, 1, 1),
  varied = c(2, 1, 3),
  negative_middle = c(1, -1, 1),
  all_negative = c(-1, -1, -1),
  zeros = c(0, 0, 0),
  one_zero = c(0, 1, 1),
  sum_zero = c(1, -1, 0),
  inf_first = c(Inf, 1, 1),
  nan_first = c(NaN, 1, 1),
  nan_two = c(NaN, NaN, 1),
  nan_all = c(NaN, NaN, NaN),
  # A negative infinity is not symmetric with `inf_first`. The C kernel divides
  # by the total weight only under `if (totweights > 0)`
  # (SOURCES/TraMineR/src/weightedinertia.cpp:34), which is false for -Inf, so
  # the undivided sum is returned and keeps its sign -- where +Inf passes the
  # guard and collapses to Inf/Inf = NaN. The kernel's inner loop visits i < j;
  # column-major matrix storage selects the logical lower triangle. It never
  # evaluates a diagonal `w_i * d_ii * w_i` and
  # never forms the +-Inf * 0 = NaN that a whole-matrix product would.
  neg_inf_middle = c(1, -Inf, 1),
  neg_inf_all = c(-Inf, -Inf, -Inf),
  mixed_inf = c(Inf, -Inf, 1)
)
weight_edge_rows <- t(vapply(
  weight_edge_cases,
  function(w) {
    c(
      dissvar(Dw, weights = w),
      dissvar(Dw, weights = w, squared = TRUE),
      .Call(TraMineR:::C_tmrWeightedInertiaDist, Dw, 3L, 0L, 1:3, as.double(w), 0L)
    )
  },
  numeric(3)
))
write_phase2_csv(
  "dissvar_weight_edge_cases", weight_edge_rows,
  paste0(
    "# R: TraMineR::dissvar on D=matrix(c(0,1,2,1,0,3,2,3,0),3,byrow=TRUE);",
    " rows = ", paste(names(weight_edge_cases), collapse = ","),
    "; cols = [dissvar, dissvar(squared=TRUE), tmrWeightedInertiaDist(var=FALSE)]\n"
  )
)

cat("DISSVAR_GOLDENS_DONE\n")
