# Live-R goldens for dissvar_grp (RP4, TraMineRextras dissvar.grp).
#
# Per-group discrepancy (pseudo-variance) on the shared phase2 15x15 distance
# matrix (3 latent groups of 5). The unweighted goldens come straight from
# TraMineRextras dissvar.grp; the weighted golden is the per-group dissvar with
# subset weights (R's dissvar.grp does not subset case weights).
#
# Requires R 4.6.2 + TraMineR 2.2-14 + TraMineRextras 0.6.9.

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
source("scripts/r_goldens/make_goldens_phase2_inputs.R")  # phase2_inputs(), write_phase2_csv()
suppressMessages(library(TraMineR))
suppressMessages(library(TraMineRextras))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
stopifnot(as.character(packageVersion("TraMineRextras")) == "0.6.9")

inp <- phase2_inputs()
D <- inp$D
grp <- inp$grp
w <- inp$wreal

write_phase2_csv(
  "dissvar_grp_uw", matrix(as.vector(dissvar.grp(D, group = grp)), nrow = 1),
  "# dissvar.grp(D, group) unweighted; group order = sorted levels\n"
)
write_phase2_csv(
  "dissvar_grp_uw_sq", matrix(as.vector(dissvar.grp(D, group = grp, squared = TRUE)), nrow = 1),
  "# dissvar.grp(D, group, squared=TRUE)\n"
)

levs <- sort(unique(grp))
vw <- sapply(levs, function(l) {
  ig <- which(grp == l)
  dissvar(D[ig, ig], weights = w[ig])
})
write_phase2_csv(
  "dissvar_grp_w", matrix(vw, nrow = 1),
  "# per-group dissvar with subset weights w[ig]\n"
)

cat("DISSVAR_GRP_GOLDENS_DONE\n")
