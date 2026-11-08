# Golden generation: seqLLCS / seqLLCP (T12) -- scalar common-subsequence / prefix length.
#
# Pins the BREAKING re-point of traminex.seqLLCS / seqLLCP onto the live R contract
# (TraMineR 2.2-14): each is a 2-arg function returning the *scalar* length of the
# longest common subsequence (seqLLCS) / longest common prefix (seqLLCP) between two
# single sequences. R compares over seqlength() positions (trailing void excluded;
# nr "*" kept as a state). We emit the full pairwise matrix over the A/B/C fixtures so
# the void/missing paths are all covered.
#
# Run from the repository root with the oracle library (scripts/r_goldens/README.md):
#   Rscript scripts/r_goldens/make_goldens_seqLLCS.R
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
source(file.path("scripts", "r_goldens", "make_goldens_corestats.R"))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

emit_pairs <- function(sd, tag) {
  n <- nrow(sd)
  llcs <- matrix(0L, n, n)
  llcp <- matrix(0L, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      llcs[i, j] <- seqLLCS(sd[i, ], sd[j, ])
      llcp[i, j] <- seqLLCP(sd[i, ], sd[j, ])
    }
  }
  write_num_matrix(
    paste0("seqLLCS_", tag), llcs,
    sprintf("# R: seqLLCS(%s[i,], %s[j,]) pairwise length-of-LCS matrix\n", tag, tag)
  )
  write_num_matrix(
    paste0("seqLLCP_", tag), llcp,
    sprintf("# R: seqLLCP(%s[i,], %s[j,]) pairwise length-of-LCP matrix\n", tag, tag)
  )
}

fx <- corestats_fixtures()
emit_pairs(fx$A, "A")
emit_pairs(fx$B, "B")
emit_pairs(fx$C, "C")
cat("SEQLLCS_GOLDENS_DONE\n")
