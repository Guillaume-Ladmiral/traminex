# Golden generation round 2: TWED / OMslen / OMspell / OMstran / OM norm / seqcost.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(packageVersion("TraMineR") == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

write_omspell_fixture <- function(D) {
  header <- c(
    "# R: TraMineR 2.2-14 (R 4.6.2, Ubuntu 24.04), values dumped with %.15g",
    "# sts <- seqdef(rbind(",
    "#   c(\"A\",\"A\",\"A\",\"B\",\"B\",\"C\",\"C\",\"C\"),",
    "#   c(\"A\",\"B\",\"A\",\"B\",\"C\",\"B\",\"C\",\"A\"),",
    "#   c(\"C\",\"C\",\"C\",\"C\",\"A\",\"A\",\"A\",\"A\"),",
    "#   c(\"A\",\"A\",\"B\",\"B\",\"B\",\"B\",\"C\",\"C\"),",
    "#   c(\"B\",\"B\",\"B\",\"B\",\"B\",\"B\",\"B\",\"B\"),",
    "#   c(\"C\",\"A\",\"C\",\"A\",\"B\",\"C\",\"B\",\"A\")), alphabet=c(\"A\",\"B\",\"C\"))",
    "# cmat <- matrix(2, 3, 3); diag(cmat) <- 0",
    "# seqdist(sts, method=\"OMspell\", sm=cmat, indel=1, expcost=0.5)"
  )
  rows <- apply(
    as.matrix(D),
    1,
    function(row) {
      paste(trimws(formatC(row, format = "g", digits = 15)), collapse = ",")
    }
  )
  writeLines(
    c(header, rows),
    file.path(out_dir, "OMspell_exp05.csv"),
    useBytes = TRUE
  )
}

# --- TWED ---
dump_mat("TWED_nu05_h1", seqdist(sts, method = "TWED", nu = 0.5, h = 1.0, sm = cmat))
dump_mat("TWED_large_nu50_h25", seqdist(sts, method = "TWED", nu = 50, h = 25, sm = cmat))

# --- OMslen (R default h = 0.5) ---
dump_mat("OMslen_mean", seqdist(sts, method = "OMslen", sm = cmat, indel = 1, link = "mean"))
dump_mat("OMslen_gmean", seqdist(sts, method = "OMslen", sm = cmat, indel = 1, link = "gmean"))
dump_mat("OMslen_mean_h1", seqdist(sts, method = "OMslen", sm = cmat, indel = 1, link = "mean", h = 1.0))

# --- OMspell ---
omspell_exp05 <- TraMineR::seqdist(
  sts,
  method = "OMspell",
  sm = cmat,
  indel = 1,
  expcost = 0.5
)
dump_mat("OMspell_exp05", omspell_exp05)
write_omspell_fixture(omspell_exp05)
dump_mat("OMspell_exp15", seqdist(sts, method = "OMspell", sm = cmat, indel = 1, expcost = 1.5))

# --- OMstran ---
dump_mat("OMstran_otto05", seqdist(sts, method = "OMstran", sm = cmat, indel = 1, otto = 0.5))
dump_mat("OMstran_otto1", seqdist(sts, method = "OMstran", sm = cmat, indel = 1, otto = 1.0))
dump_mat("OMstran_otto05_prob", seqdist(sts, method = "OMstran", sm = cmat, indel = 1, otto = 0.5, transindel = "prob"))

# --- OM with maximum-distance normalization ---
dump_mat("OM_norm_max", seqdist(sts, method = "OM", sm = cmat, indel = 1, norm = "maxdist"))

# --- seqcost goldens ---
sc_log <- seqcost(sts, method = "INDELSLOG")
cat("INDELSLOG_indel\n")
cat(paste(sprintf("%.15g", sc_log$indel), collapse = ","), "\n")
dump_mat("INDELSLOG_sm", sc_log$sm)

sc_fut <- seqcost(sts, method = "FUTURE")
dump_mat("FUTURE_sm", sc_fut$sm)

feat <- data.frame(x = c(0, 0.5, 1))
sc_feat <- seqcost(sts, method = "FEATURES", state.features = feat)
dump_mat("FEATURES_sm", sc_feat$sm)

cat("ALL_GOLDENS2_DONE\n")
