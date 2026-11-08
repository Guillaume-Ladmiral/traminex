# Golden generation: seqmeant (T6) -- mean time spent in each state.
#
# Pins traminex.seqmeant to TraMineR 2.2-14 seqmeant(weighted, with.missing,
# prop, serr): index = statl (real alphabet, +nr if with.missing), void never a
# row.
#   *_ms  : prop=FALSE, serr=TRUE  -> columns [Mean, Var, Stdev, SE]
#   *_prop: prop=TRUE,  serr=FALSE -> column  [Mean] (per-seq proportions averaged)
#
# Run from the repository root with the oracle library (scripts/r_goldens/README.md):
#   Rscript scripts/r_goldens/make_goldens_seqmeant.R
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
source(file.path("scripts", "r_goldens", "make_goldens_corestats.R"))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

emit <- function(sd, tag, wm) {
  ms <- seqmeant(sd, weighted = TRUE, with.missing = wm, prop = FALSE, serr = TRUE)
  write_num_matrix(
    paste0("seqmeant_", tag, "_ms"), as.matrix(ms),
    sprintf("# R: seqmeant(%s, with.missing=%s, serr=TRUE) cols=[Mean,Var,Stdev,SE]\n", tag, wm)
  )
  pr <- seqmeant(sd, weighted = TRUE, with.missing = wm, prop = TRUE, serr = FALSE)
  write_num_matrix(
    paste0("seqmeant_", tag, "_prop"), as.matrix(pr),
    sprintf("# R: seqmeant(%s, with.missing=%s, prop=TRUE) col=[Mean]\n", tag, wm)
  )
}

fx <- corestats_fixtures()
emit(fx$A, "A", FALSE)
emit(fx$Aw, "Aw", FALSE)
emit(fx$B, "B_wm0", FALSE)
emit(fx$B, "B_wm1", TRUE)
emit(fx$C, "C_wm0", FALSE)
emit(fx$C, "C_wm1", TRUE)
cat("SEQMEANT_GOLDENS_DONE\n")
