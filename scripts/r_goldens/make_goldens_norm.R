# Golden generation: OM-family distance NORMALIZATION (TWED/OMslen/OMspell/OMstran).
# Produces the matrices pinned by tests/test_om_normalization.py.
#
# Requires R + TraMineR 2.2-14. On the dev machine this was run against TraMineR
# installed in the Ubuntu R 4.6.2; see scripts/r_goldens/README.md for the run
# command.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

seqs <- rbind(
  c("A","A","A","B","B","C","C","C"),
  c("A","B","A","B","C","B","C","A"),
  c("C","C","C","C","A","A","A","A"),
  c("A","A","B","B","B","B","C","C"),
  c("B","B","B","B","B","B","B","B"),
  c("C","A","C","A","B","C","B","A")
)
sts <- suppressMessages(seqdef(seqs, alphabet = c("A","B","C")))
cmat <- matrix(2, 3, 3); diag(cmat) <- 0
rownames(cmat) <- colnames(cmat) <- c("A","B","C")

dump_mat <- function(name, D) {
  D <- as.matrix(D)
  cat("BEGIN", name, "\n")
  for (i in 1:nrow(D)) cat(paste(sprintf("%.17g", D[i, ]), collapse = ","), "\n")
  cat("END", name, "\n")
}

run <- function(name, ...) dump_mat(name, suppressMessages(seqdist(sts, ...)))

norms <- c("maxlength","gmean","maxdist","YujianBo","auto")
for (nm in norms) run(sprintf("TWED_nu05_h1_%s", nm), method="TWED", nu=0.5, h=1.0, sm=cmat, norm=nm)
for (nm in norms) run(sprintf("OMslen_mean_%s", nm), method="OMslen", sm=cmat, indel=1, link="mean", norm=nm)
for (nm in norms) run(sprintf("OMslen_gmean_%s", nm), method="OMslen", sm=cmat, indel=1, link="gmean", norm=nm)
for (nm in norms) run(sprintf("OMspell_exp05_%s", nm), method="OMspell", sm=cmat, indel=1, expcost=0.5, norm=nm)
for (nm in norms) run(sprintf("OMstran_otto05_%s", nm), method="OMstran", sm=cmat, indel=1, otto=0.5, norm=nm)

# Refusal parity: norm must be rejected for NMS/NMSMST/SVRspell (seqdist.R line 205).
for (mth in c("NMS","NMSMST","SVRspell")) {
  res <- tryCatch({ suppressMessages(seqdist(sts, method=mth, kweights=rep(1,8), norm="maxdist")); "NO_ERROR" },
                  error=function(e) paste0("ERROR: ", conditionMessage(e)))
  cat(sprintf("REFUSE %s: %s\n", mth, res))
}
