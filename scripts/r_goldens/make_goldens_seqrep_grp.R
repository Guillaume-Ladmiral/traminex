# Golden generation: seqrep.grp -> seqrep_grp (RP6a BUILD) vs TraMineRextras 0.6.9.
#
# seqrep.grp(seqdata, group, diss, ret) runs seqrep() per group with a SHARED
# dmax = max(diss) and returns (for ret="stat") the per-group quality Statistics.
# We reuse the exact seqrep fixture of make_goldens_dissrep.R (set.seed(99001),
# 10 seqs x 6, alphabet 1:4, weights 0.5+0.1*i, LCS norm=none diss) and split it
# into two groups of five. The fixture (codes/diss/weights/group) is dumped so the
# Python test reconstructs the identical inputs.
#
# Run from the repository root with the oracle library (scripts/r_goldens/README.md):
#   Rscript scripts/r_goldens/make_goldens_seqrep_grp.R
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
source("scripts/r_goldens/make_goldens_phase2_inputs.R")  # write_phase2_csv()
suppressMessages(library(TraMineR))
suppressMessages(library(TraMineRextras))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
stopifnot(as.character(packageVersion("TraMineRextras")) == "0.6.9")

set.seed(99001)
nr <- 10
Lr <- 6
Mr <- matrix(sample(1:4, nr * Lr, replace = TRUE), nrow = nr, ncol = Lr)
wr <- 0.5 + 0.1 * (1:nr)
sdr <- suppressMessages(seqdef(Mr, alphabet = 1:4, weights = wr))
Dr <- suppressMessages(seqdist(sdr, method = "LCS", norm = "none"))
grp <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2)

write_phase2_csv("seqrep_grp_input_group", matrix(grp, nrow = 1), "# seqrep_grp group vector\n")

stat_w <- seqrep.grp(sdr, group = grp, diss = Dr, ret = "stat", criterion = "density",
  coverage = 0.25, weighted = TRUE)
stat_uw <- seqrep.grp(sdr, group = grp, diss = Dr, ret = "stat", criterion = "density",
  coverage = 0.25, weighted = FALSE)

emit <- function(stat, suff) {
  levs <- names(stat)
  for (k in seq_along(levs)) {
    write_phase2_csv(
      paste0("seqrep_grp_", suff, "_g", k, "_stats"), as.matrix(stat[[k]]),
      sprintf("# seqrep.grp group %s Statistics: na,na(%%),nb,nb(%%),SD,MD,DC,V,Q\n", trimws(levs[k]))
    )
  }
}
emit(stat_w, "w")
emit(stat_uw, "uw")

cat("SEQREP_GRP_GOLDENS_DONE\n")
