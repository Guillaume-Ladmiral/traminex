# Golden generation: seqalign (T17 VERIFY) vs TraMineR 2.2-14.
#
# seqalign(seqdata, indices, indel, sm, with.missing) returns the optimal-matching
# alignment of two sequences: the operation list (E/S/I/D), the aligned seq1/seq2
# (with "-" for gaps), and the per-operation cost (sum = OM distance). R aligns over
# seqlength() positions, i.e. trailing void is excluded; nr "*" is kept as a state
# (so the missing path needs with.missing=TRUE + a matching sm). We emit one text
# golden per case: 4 lines = operations / seq1c / seq2c / costs (%.15g).
#
# Run from the repository root with the oracle library (scripts/r_goldens/README.md):
#   Rscript scripts/r_goldens/make_goldens_seqalign.R
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
source(file.path("scripts", "r_goldens", "make_goldens_corestats.R"))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")

emit_align <- function(name, sd, indices, indel, sm, with.missing = FALSE) {
  sa <- TraMineR:::seqalign(sd, indices, indel = indel, sm = sm, with.missing = with.missing)
  con <- file(file.path(out_dir, paste0(name, ".txt")), open = "wt")
  writeLines(paste(sa$operation, collapse = " "), con)
  writeLines(paste(sa$seq1, collapse = " "), con)
  writeLines(paste(sa$seq2, collapse = " "), con)
  writeLines(paste(formatC(sa$cost, format = "g", digits = 15), collapse = " "), con)
  close(con)
  # cycle-6 section 9: the two DP matrices are part of R's public return value and
  # were unpinned, which hid an origin-cell difference. R fills the first column
  # with "D" and *then* the first row with "I" (seqalign.R:22-23), so [1,1] ends
  # up "I". Empty cells are written as "." so the row stays token-aligned.
  om <- sa$opmatrix
  om[om == ""] <- "."
  con <- file(file.path(out_dir, paste0(name, "_matrices.txt")), open = "wt")
  writeLines(paste(dim(om), collapse = " "), con)
  for (r in seq_len(nrow(om))) writeLines(paste(om[r, ], collapse = " "), con)
  for (r in seq_len(nrow(sa$costmatrix))) {
    writeLines(paste(formatC(sa$costmatrix[r, ], format = "g", digits = 15), collapse = " "), con)
  }
  close(con)
}

fx <- corestats_fixtures()

# Clean (void-free, no missing) -- the canonical seqalign use.
smA_const <- seqsubm(fx$A, method = "CONSTANT", cval = 2)
smA_trate <- suppressMessages(seqsubm(fx$A, method = "TRATE"))
# Dump the TRATE sm so the Python test feeds the same matrix (decoupled from seqsubm).
write_num_matrix("seqalign_smA_trate", smA_trate, "# R: seqsubm(A, method='TRATE')\n")
emit_align("seqalign_A12_const", fx$A, c(1, 2), 1, smA_const)
emit_align("seqalign_A12_trate", fx$A, c(1, 2), 1, smA_trate)
emit_align("seqalign_A24_const", fx$A, c(2, 4), 1, smA_const)

# Ragged (trailing void) + interior missing -- with.missing=TRUE keeps nr as a state,
# trailing void is dropped. Pins the void-stripping path.
smC_const <- seqsubm(fx$C, method = "CONSTANT", cval = 2, with.missing = TRUE)
emit_align("seqalign_C13_const_wm", fx$C, c(1, 3), 1, smC_const, with.missing = TRUE)
emit_align("seqalign_C12_const_wm", fx$C, c(1, 2), 1, smC_const, with.missing = TRUE)

cat("SEQALIGN_GOLDENS_DONE\n")
