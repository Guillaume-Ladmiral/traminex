# Live-R goldens for dissindic (RP4, TraMineRextras dissindic).
#
# Individual marginality / gain indicators on the shared phase2 15x15 distance
# matrix (3 latent groups of 5).
#
# AUTHENTICITY (RP3/RP4 audit, F3): TraMineRextras 0.6.9 dissindic() works with a
# FACTOR grouping variable -- only the *numeric* path errors (its internal
# dissmfacwresid does `data[seqok, ]` on a one-column data.frame, which R collapses
# to a vector). So these goldens come from the DIRECT `dissindic(D, factor(group))`
# call (not a replicated formula). traminex's dissindic always treats `group`
# categorically, which is the only working R semantics (R's numeric->continuous
# model.matrix path is non-functional), so there is no numeric-vs-factor fork.
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

gen <- function(tag, weights, squared) {
  r <- suppressMessages(
    dissindic(D, factor(grp), gower = FALSE, squared = squared, weights = weights)
  )
  write_phase2_csv(
    paste0("dissindic_", tag, "_marg"), matrix(r$marginality, nrow = 1),
    sprintf("# R: dissindic(D, factor(grp), squared=%s, weights=%s)$marginality\n",
            squared, if (is.null(weights)) "NULL" else "wreal")
  )
  write_phase2_csv(
    paste0("dissindic_", tag, "_gain"), matrix(r$gain, nrow = 1),
    sprintf("# R: dissindic(D, factor(grp), squared=%s, weights=%s)$gain\n",
            squared, if (is.null(weights)) "NULL" else "wreal")
  )
}

gen("uw", NULL, FALSE)
gen("w", inp$wreal, FALSE)
gen("w_sq", inp$wreal, TRUE)

cat("DISSINDIC_GOLDENS_DONE\n")
