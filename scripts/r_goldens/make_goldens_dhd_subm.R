# Golden generation: DHD / HAM with a caller-supplied substitution matrix.
#
# Anchors traminex's pairwise(method="DHD", subm=...) explicit-substitution path
# (added in Phase 1 / A1 so traminex can faithfully serve weightedcluster's DHD)
# to TraMineR 2.2-14:
#   - a static (K,K) sm   == R seqdist(method="HAM", sm=sm, norm="none")
#   - a time-varying sm   == R seqdist(method="DHD", sm=sm, norm="none")
#
# Requires R + TraMineR 2.2-14; see scripts/r_goldens/README.md for the run
# command. Writes one CSV per golden into tests/r_goldens/.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

seqs <- rbind(
  c("A", "A", "A", "B", "B", "C", "C", "C"),
  c("A", "B", "A", "B", "C", "B", "C", "A"),
  c("C", "C", "C", "C", "A", "A", "A", "A"),
  c("A", "A", "B", "B", "B", "B", "C", "C"),
  c("B", "B", "B", "B", "B", "B", "B", "B"),
  c("C", "A", "C", "A", "B", "C", "B", "A")
)
sts <- suppressMessages(seqdef(seqs, alphabet = c("A", "B", "C")))
K <- 3
L <- ncol(seqs)

write_csv <- function(name, D, header) {
  D <- as.matrix(D)
  path <- file.path(out_dir, paste0(name, ".csv"))
  con <- file(path, "w")
  cat(header, file = con)
  for (i in seq_len(nrow(D))) {
    cat(paste(sprintf("%.15g", D[i, ]), collapse = ","), "\n", file = con, sep = "")
  }
  close(con)
}

# Static (K,K) substitution matrix: off-diagonal 2, zero diagonal.
sm2 <- matrix(2, K, K)
diag(sm2) <- 0
rownames(sm2) <- colnames(sm2) <- c("A", "B", "C")
D_ham <- suppressMessages(seqdist(sts, method = "HAM", sm = sm2, norm = "none"))
write_csv(
  "DHD_subm_static_sm2_none", D_ham,
  "# R: seqdist(method='HAM', sm=matrix(2,3,3) diag 0, norm='none') on the 6x8 ABC sequences\n"
)

# Time-varying substitution matrix, R layout dim=c(K,K,L): base 2 (diag 0),
# with the mismatch at position 5 made expensive (cost 5).
smtv <- array(2, dim = c(K, K, L))
for (t in seq_len(L)) diag(smtv[, , t]) <- 0
smtv[1, 2, 5] <- smtv[2, 1, 5] <- 5
smtv[1, 3, 5] <- smtv[3, 1, 5] <- 5
smtv[2, 3, 5] <- smtv[3, 2, 5] <- 5
D_dhd <- suppressMessages(seqdist(sts, method = "DHD", sm = smtv, norm = "none"))
write_csv(
  "DHD_subm_timevarying_none", D_dhd,
  "# R: seqdist(method='DHD', sm=array(K,K,L) base 2 diag 0, pos5 offdiag 5, norm='none')\n"
)

cat("DHD_SUBM_GOLDENS_DONE\n")
