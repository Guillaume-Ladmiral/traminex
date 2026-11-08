# Golden generation: seqmodst (T7) -- modal state sequence.
#
# Pins traminex.seqmodst to TraMineR 2.2-14 seqmodst(weighted, with.missing):
#   *_cf  : row0 = modal-state codes seqasnum(m, with.missing=TRUE) (void -> -1),
#           row1 = attr(m,"Frequencies") (max freq at each position)
#   *_occ : attr(m,"Occurrences") (# sequences exactly matching the modal seq)
#
# Run from the repository root with the oracle library (scripts/r_goldens/README.md):
#   Rscript scripts/r_goldens/make_goldens_seqmodst.R
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
source(file.path("scripts", "r_goldens", "make_goldens_corestats.R"))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

# seqasnum is an internal TraMineR function (not exported) -> reach with :::.
asnum <- TraMineR:::seqasnum

emit <- function(sd, tag, wm) {
  m <- seqmodst(sd, weighted = TRUE, with.missing = wm)
  codes <- as.vector(asnum(m, with.missing = TRUE))
  codes[is.na(codes)] <- -1
  freq <- as.vector(attr(m, "Frequencies"))
  write_num_matrix(
    paste0("seqmodst_", tag, "_cf"), rbind(codes, freq),
    sprintf("# R: seqmodst(%s, with.missing=%s) rows=[seqasnum codes (void=-1), Frequencies]\n", tag, wm)
  )
  write_num_matrix(
    paste0("seqmodst_", tag, "_occ"), matrix(attr(m, "Occurrences"), nrow = 1),
    sprintf("# R: seqmodst(%s, with.missing=%s) Occurrences\n", tag, wm)
  )
}

fx <- corestats_fixtures()
emit(fx$A, "A", FALSE)
emit(fx$Aw, "Aw", FALSE)
emit(fx$B, "B_wm1", TRUE)
emit(fx$C, "C_wm1", TRUE)
cat("SEQMODST_GOLDENS_DONE\n")
