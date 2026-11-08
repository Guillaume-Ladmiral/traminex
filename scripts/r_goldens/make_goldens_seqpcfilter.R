# Golden generation: RP6b seqpcfilter (parallel-coordinate line filter functions).
# Anchors traminex's seqpcfilter line rules to TraMineR 2.2-14 minfreq/cumfreq/linear.
# Run in the live-R env:  Rscript scripts/r_goldens/make_goldens_seqpcfilter.R
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
out_dir <- file.path("tests", "r_goldens")

write_num_matrix <- function(name, m, header) {
  con <- file(file.path(out_dir, paste0(name, ".csv")), open = "wt")
  cat(header, file = con)
  txt <- apply(as.matrix(m), 1, function(r) paste(formatC(r, format = "g", digits = 15), collapse = ","))
  writeLines(txt, con)
  close(con)
}

minfreq <- TraMineR:::minfreq
cumfreq <- TraMineR:::cumfreq
linear <- TraMineR:::linear

x <- c(0.05, 0.15, 0.3, 0.1, 0.4)  # distinct line frequencies (sum = 1)

write_num_matrix("seqpcfilter_minfreq_01", t(minfreq(x, 0.1)),  "# R: minfreq(c(.05,.15,.3,.1,.4), 0.1)\n")
write_num_matrix("seqpcfilter_minfreq_02", t(minfreq(x, 0.2)),  "# R: minfreq(..., 0.2)\n")
write_num_matrix("seqpcfilter_cumfreq_05", t(cumfreq(x, 0.5)),  "# R: cumfreq(..., 0.5)\n")
write_num_matrix("seqpcfilter_cumfreq_075", t(cumfreq(x, 0.75)), "# R: cumfreq(..., 0.75)\n")
write_num_matrix("seqpcfilter_linear", t(linear(x)),            "# R: linear(...)\n")

cat("SEQPCFILTER_GOLDENS_DONE\n")
