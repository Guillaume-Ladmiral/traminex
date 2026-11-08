# Golden generation: seqstatf (T10) -- overall state frequency table.
#
# Pins traminex.seqstatf (re-targeted from first-position to overall frequency)
# to TraMineR 2.2-14 seqstatf(weighted, with.missing): Freq = colSums(seqistatd *
# weights), Percent = Freq/sum(Freq)*100, index = statl (alphabet, +nr if
# with.missing; never void).
#
# Run from the repository root with the oracle library (scripts/r_goldens/README.md):
#   Rscript scripts/r_goldens/make_goldens_seqstatf.R
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
source(file.path("scripts", "r_goldens", "make_goldens_corestats.R"))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

emit <- function(sd, tag, wm) {
  f <- seqstatf(sd, weighted = TRUE, with.missing = wm)
  write_num_matrix(
    paste0("seqstatf_", tag), as.matrix(f),
    sprintf("# R: seqstatf(%s, weighted=TRUE, with.missing=%s) cols=[Freq, Percent]\n", tag, wm)
  )
}

fx <- corestats_fixtures()
emit(fx$A, "A", FALSE)
emit(fx$Aw, "Aw", FALSE)
emit(fx$B, "B_wm0", FALSE)
emit(fx$B, "B_wm1", TRUE)
emit(fx$C, "C_wm0", FALSE)
emit(fx$C, "C_wm1", TRUE)
cat("SEQSTATF_GOLDENS_DONE\n")
