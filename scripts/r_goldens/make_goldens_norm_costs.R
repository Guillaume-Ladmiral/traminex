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
  for (i in 1:nrow(D)) cat(paste(sprintf("%.15g", D[i, ]), collapse = ","), "\n")
  cat("END", name, "\n")
}

dump_mat("OM_norm_maxdist", seqdist(sts, method = "OM", sm = cmat, indel = 1, norm = "maxdist"))
dump_mat("OM_norm_maxlength", seqdist(sts, method = "OM", sm = cmat, indel = 1, norm = "maxlength"))

sc_log <- seqcost(sts, method = "INDELSLOG")
cat("INDELSLOG_indel\n")
cat(paste(sprintf("%.15g", sc_log$indel), collapse = ","), "\n")
dump_mat("INDELSLOG_sm", sc_log$sm)

sc_fut <- seqcost(sts, method = "FUTURE")
dump_mat("FUTURE_sm", sc_fut$sm)

feat <- data.frame(x = c(0, 0.5, 1))
sc_feat <- seqcost(sts, method = "FEATURES", state.features = feat)
dump_mat("FEATURES_sm", sc_feat$sm)

cat("ALL_GOLDENS2B_DONE\n")
