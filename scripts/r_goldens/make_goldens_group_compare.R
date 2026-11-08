# Live-R goldens for dissCompare / seqLRT / seqBIC / seqCompare (RP4, s=0).
#
# dissCompare on the phase2 15x15 distance matrix with a binary group; the
# seqCompare/seqLRT/seqBIC goldens (sequence inputs) are appended in the
# seqCompare commit. Only the deterministic s=0 path is pinned.
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
g2 <- ifelse(grp <= 1, 1L, 2L)  # binary group: 5 vs 10

r_uw <- suppressMessages(dissCompare(D, group = g2, s = 0, stat = "all", weighted = FALSE))
r_w <- suppressMessages(
  dissCompare(D, group = g2, s = 0, stat = "all", weighted = TRUE, weights = w)
)

write_phase2_csv(
  "dissCompare_uw", matrix(as.vector(r_uw), nrow = 1),
  "# dissCompare(D, g2, s=0, stat=all, weighted=FALSE): LRT, p, dBIC, BF\n"
)
write_phase2_csv(
  "dissCompare_w", matrix(as.vector(r_w), nrow = 1),
  "# dissCompare(D, g2, s=0, stat=all, weighted=TRUE, weights=w): LRT, p, dBIC, BF\n"
)

cat("GROUP_COMPARE_DISS_GOLDENS_DONE\n")

## seqCompare / seqLRT / seqBIC goldens (sequence inputs, s=0) ------------------
set.seed(7777)
nS <- 12
LS <- 6
MS <- matrix(sample(1:4, nS * LS, replace = TRUE), nrow = nS, ncol = LS)
gS <- rep(c(1, 2), length.out = nS)
wS <- 0.5 + 0.1 * (1:nS)
sdS <- suppressMessages(seqdef(MS, alphabet = 1:4, weights = wS))

write_phase2_csv("seqcompare_input_codes", MS - 1L, "# seqCompare input: 0-based codes (alphabet 1:4)\n")
write_phase2_csv("seqcompare_input_group", matrix(gS, nrow = 1), "# seqCompare input: binary group\n")
write_phase2_csv("seqcompare_input_weights", matrix(wS, nrow = 1), "# seqCompare input: weights\n")

genS <- function(tag, weighted) {
  lrt <- suppressMessages(seqLRT(sdS, group = gS, s = 0, method = "LCS", weighted = weighted))
  bic <- suppressMessages(seqBIC(sdS, group = gS, s = 0, method = "LCS", weighted = weighted))
  cmp <- suppressMessages(
    seqCompare(sdS, group = gS, s = 0, stat = "all", method = "LCS", weighted = weighted)
  )
  write_phase2_csv(paste0("seqLRT_", tag), matrix(as.vector(lrt), nrow = 1), "# seqLRT: LRT, p\n")
  write_phase2_csv(paste0("seqBIC_", tag), matrix(as.vector(bic), nrow = 1), "# seqBIC: dBIC, BF\n")
  write_phase2_csv(
    paste0("seqCompare_", tag), matrix(as.vector(cmp), nrow = 1), "# seqCompare: LRT,p,dBIC,BF\n"
  )
}

genS("uw", FALSE)
genS("w", TRUE)

cat("GROUP_COMPARE_SEQ_GOLDENS_DONE\n")
