# Golden generation: seqtransn (T4) -- per-sequence transition counts.
#
# Pins traminex.seqtransn (re-pointed from the K x K matrix to the per-sequence
# count) to TraMineR 2.2-14 seqtransn(with.missing, norm, pweight):
#   *_<tag> : n x 2 = [norm=FALSE (counts), norm=TRUE (counts/(seqlength-1))]
#   *_pweight_<tag> : n x 1 = seqtransn(pweight=TRUE) (1 - transition-rate weighted)
# norm divides by the non-void length - 1.
#
# Run from the repository root with the oracle library (scripts/r_goldens/README.md):
#   Rscript scripts/r_goldens/make_goldens_seqtransn.R
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
source(file.path("scripts", "r_goldens", "make_goldens_corestats.R"))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

emit <- function(sd, tag, wm) {
  c0 <- as.vector(seqtransn(sd, with.missing = wm, norm = FALSE))
  c1 <- as.vector(seqtransn(sd, with.missing = wm, norm = TRUE))
  write_num_matrix(
    paste0("seqtransn_", tag), cbind(c0, c1),
    sprintf("# R: seqtransn(%s, with.missing=%s) cols=[norm=FALSE, norm=TRUE]\n", tag, wm)
  )
}

fx <- corestats_fixtures()
emit(fx$A, "A", FALSE)
emit(fx$B, "B_wm0", FALSE)
emit(fx$B, "B_wm1", TRUE)
emit(fx$C, "C_wm0", FALSE)
emit(fx$C, "C_wm1", TRUE)

# pweight (with.missing=FALSE; uses default seqtrate weighted=TRUE)
write_num_matrix(
  "seqtransn_pweight_A", matrix(as.vector(seqtransn(fx$A, pweight = TRUE)), ncol = 1),
  "# R: seqtransn(A, pweight=TRUE)\n"
)
write_num_matrix(
  "seqtransn_pweight_C", matrix(as.vector(seqtransn(fx$C, pweight = TRUE)), ncol = 1),
  "# R: seqtransn(C, pweight=TRUE) ragged/void\n"
)
cat("SEQTRANSN_GOLDENS_DONE\n")
