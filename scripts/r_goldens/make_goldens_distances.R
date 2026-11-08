# Golden generation for traminex verification.
# TraMineR 2.2-14 compiled from SOURCES/TraMineR_2.2-14.tar.gz (vendored ground truth).
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

dump_mat <- function(name, D) {
  D <- as.matrix(D)
  cat("BEGIN", name, "\n")
  for (i in 1:nrow(D)) cat(paste(sprintf("%.15g", D[i, ]), collapse = ","), "\n")
  cat("END", name, "\n")
}

# --- CHI2: R defaults (step=1, overlap=FALSE, norm='none') ---
dump_mat("CHI2_default", seqdist(sts, method = "CHI2"))
# CHI2 normalized ('auto' divides by sqrt(#intervals))
dump_mat("CHI2_norm_auto", seqdist(sts, method = "CHI2", norm = "auto"))
# CHI2 step=2 (interval aggregation), and overlap variant
dump_mat("CHI2_step2", seqdist(sts, method = "CHI2", step = 2))
dump_mat("CHI2_step2_overlap", seqdist(sts, method = "CHI2", step = 2, overlap = TRUE))
# Full-length single interval (what the implementer's property tests used)
dump_mat("CHI2_stepL", seqdist(sts, method = "CHI2", step = 8))

# --- EUCLID: same family, default + stepL ---
dump_mat("EUCLID_default", seqdist(sts, method = "EUCLID"))
dump_mat("EUCLID_stepL", seqdist(sts, method = "EUCLID", step = 8))

# --- DHD: TraMineR's transition-rate Dynamic Hamming (time-varying TRATE, cval=4) ---
dump_mat("DHD_default", seqdist(sts, method = "DHD"))

# --- HAM: R-faithful raw count (norm='none' default, unit cost), the /L proportion
# (norm='maxlength'), and a custom substitution matrix (off-diagonal = 3) ---
dump_mat("HAM_default", seqdist(sts, method = "HAM"))
dump_mat("HAM_norm_maxlength", seqdist(sts, method = "HAM", norm = "maxlength"))
cmat3 <- matrix(3, 3, 3); diag(cmat3) <- 0
dump_mat("HAM_sm_off3", seqdist(sts, method = "HAM", sm = cmat3))

# --- NMS / NMSMST / SVRspell ---
kw <- rep(1, 8)
dump_mat("NMS_default", seqdist(sts, method = "NMS", kweights = kw))
dump_mat("NMSMST_default", seqdist(sts, method = "NMSMST", kweights = kw, tpow = 1))
dump_mat("NMSMST_tpow2", seqdist(sts, method = "NMSMST", kweights = kw, tpow = 2))
dump_mat("SVRspell_default", seqdist(sts, method = "SVRspell", kweights = kw, tpow = 1))
dump_mat("SVRspell_tpow2", seqdist(sts, method = "SVRspell", kweights = kw, tpow = 2))

cat("ALL_GOLDENS_DONE\n")
