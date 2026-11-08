# Golden generation: disscenter (Phase 2.3 / DUPLICATES A7).
#
# Anchors traminex's R-faithful `disscenter` (distance of each case to its
# group centre — NOT matrix double-centering, which traminex now calls
# `double_center_matrix`) to TraMineR 2.2-14
# `disscenter(diss, group=, medoids.index=, allcenter=, weights=, squared=)`.
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
grp <- inp$grp

write_phase2_csv(
  "disscenter_nogrp_unw", matrix(disscenter(D), ncol = 1),
  "# R: disscenter(D) -> 15x1 distances to the overall centre\n"
)
write_phase2_csv(
  "disscenter_grp_unw", matrix(disscenter(D, group = grp), ncol = 1),
  "# R: disscenter(D, group=grp) -> 15x1 distances to own group centre\n"
)
write_phase2_csv(
  "disscenter_grp_wreal", matrix(disscenter(D, group = grp, weights = inp$wreal), ncol = 1),
  "# R: disscenter(D, group=grp, weights=wreal) -> 15x1 weighted distances to own group centre\n"
)
write_phase2_csv(
  "disscenter_grp_wreal_sq",
  matrix(disscenter(D, group = grp, weights = inp$wreal, squared = TRUE), ncol = 1),
  "# R: disscenter(D, group=grp, weights=wreal, squared=TRUE) -> 15x1\n"
)
write_phase2_csv(
  "disscenter_grp_wreal_allcenter",
  as.matrix(disscenter(D, group = grp, weights = inp$wreal, allcenter = TRUE)),
  "# R: disscenter(D, group=grp, weights=wreal, allcenter=TRUE) -> 15x3 distances of every case to each group centre\n"
)
med <- disscenter(D, group = grp, weights = inp$wreal, medoids.index = "first")
write_phase2_csv(
  "disscenter_grp_wreal_medoids_first", matrix(as.numeric(med), ncol = 1),
  "# R: disscenter(D, group=grp, weights=wreal, medoids.index='first') -> 3x1 medoid indices (1-based)\n"
)

## ---- cycle-6 F-003: non-finite and non-positive total case weights ----
# `disscenter` subtracts `weighted.mean(dc, weights[cond])/2` *unconditionally*
# (SOURCES/TraMineR/R/disscenter.R:87), so a zero, NaN or zero-sum weight
# vector propagates NaN/Inf into every member of the group even though the C
# kernel's own division is guarded by `if (totweights>0)`. The allcenter route
# additionally divides the external contribution by `iweight*totweights`
# (SOURCES/TraMineR/src/weightedinertia.cpp:150), un-guarded, so a zero-weight
# outsider yields 0/0 = NaN rather than the algebraically cancelled value.
#
# Rows are the weight vectors in WEIGHT_EDGE_CASES order (see
# tests/test_diss_weight_edge_cases.py). Groups are c("a","a","b").
Dw <- matrix(c(0, 1, 2, 1, 0, 3, 2, 3, 0), nrow = 3, byrow = TRUE)
grp_edge <- c("a", "a", "b")
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
  # Kept in step with make_goldens_dissvar.R's list; see the note there for why
  # -Inf is not symmetric with +Inf in the weighted-inertia kernel.
  neg_inf_middle = c(1, -Inf, 1),
  neg_inf_all = c(-Inf, -Inf, -Inf),
  mixed_inf = c(Inf, -Inf, 1)
)
write_phase2_csv(
  "disscenter_weight_edge_cases_nogrp",
  t(vapply(weight_edge_cases, function(w) as.numeric(disscenter(Dw, weights = w)), numeric(3))),
  paste0(
    "# R: TraMineR::disscenter(D, weights=w) on D=matrix(c(0,1,2,1,0,3,2,3,0),3,byrow=TRUE);",
    " rows = ", paste(names(weight_edge_cases), collapse = ","), "; cols = the 3 cases\n"
  )
)
write_phase2_csv(
  "disscenter_weight_edge_cases_grp",
  t(vapply(
    weight_edge_cases,
    function(w) as.numeric(disscenter(Dw, group = grp_edge, weights = w)),
    numeric(3)
  )),
  paste0(
    "# R: TraMineR::disscenter(D, group=c('a','a','b'), weights=w);",
    " rows = ", paste(names(weight_edge_cases), collapse = ","), "; cols = the 3 cases\n"
  )
)
write_phase2_csv(
  "disscenter_weight_edge_cases_allcenter",
  t(vapply(
    weight_edge_cases,
    function(w) {
      as.numeric(unlist(disscenter(Dw, group = grp_edge, weights = w, allcenter = TRUE)))
    },
    numeric(6)
  )),
  paste0(
    "# R: TraMineR::disscenter(D, group=c('a','a','b'), weights=w, allcenter=TRUE);",
    " rows = ", paste(names(weight_edge_cases), collapse = ","),
    "; cols = column-major 3x2 (group a then group b)\n"
  )
)

cat("DISSCENTER_GOLDENS_DONE\n")
