# Golden generation: seqemlt (TraMineRextras) -- EMLT situations / profiles /
# coordinates / correlations.
#
# Goldens the sign-stable quantities exactly (the situation correlation matrix
# sit.cor, the profile squared-distance matrix c, and the princomp standard
# deviations sdev). The sequence coordinates (coord) come from princomp on the
# square, rank-deficient Benzecri matrix, so their last component is ~0 with an
# arbitrary sign; they are compared up to a per-component sign on the
# non-degenerate components (plan decision 0.2.3).
#
# Small synthetic dataset (every situation occurs, so no NA), weighted=FALSE.
# Requires live R + TraMineR 2.2-14 + TraMineRextras 0.6.9; see memory
# "live-traminer-r-env".
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
suppressMessages(library(TraMineRextras))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
stopifnot(as.character(packageVersion("TraMineRextras")) == "0.6.9")

write_mat <- function(name, M, header) {
  out_dir <- file.path("tests", "r_goldens")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  M <- as.matrix(M)
  dimnames(M) <- NULL
  con <- file(file.path(out_dir, paste0(name, ".csv")), "w")
  cat(header, file = con)
  for (r in seq_len(nrow(M))) {
    cat(paste(sprintf("%.15g", M[r, ]), collapse = ","), "\n", file = con, sep = "")
  }
  close(con)
}

seqmat <- rbind(
  c("A", "A", "A", "A", "A"),
  c("B", "B", "B", "B", "B"),
  c("C", "C", "C", "C", "C"),
  c("A", "B", "C", "A", "B"),
  c("C", "A", "B", "C", "A"),
  c("B", "C", "A", "B", "C"),
  c("A", "C", "B", "A", "C"),
  c("B", "A", "C", "B", "A")
)
sts <- suppressMessages(seqdef(seqmat, alphabet = c("A", "B", "C")))
r <- seqemlt(sts, a = 1, b = 1, weighted = FALSE)

write_mat("seqemlt_corr", r$sit.cor,
  "# R: TraMineRextras::seqemlt(sts, weighted=FALSE)$sit.cor (situation correlations); see make_goldens_seqemlt.R\n")
write_mat("seqemlt_distc", r$c,
  "# R: seqemlt(sts, weighted=FALSE)$c (profile squared-distance matrix); see make_goldens_seqemlt.R\n")
write_mat("seqemlt_coords", r$coord,
  "# R: seqemlt(sts, weighted=FALSE)$coord (sequence coordinates on PCA axes); compare up to sign\n")
write_mat("seqemlt_sdev", matrix(r$pca$sdev, nrow = 1),
  "# R: seqemlt(sts, weighted=FALSE)$pca$sdev (princomp standard deviations)\n")

cat("SEQEMLT_GOLDENS_DONE\n")
cat("dim sit.cor:", dim(r$sit.cor), " dim coord:", dim(r$coord), " sdev:", round(r$pca$sdev, 6), "\n")
