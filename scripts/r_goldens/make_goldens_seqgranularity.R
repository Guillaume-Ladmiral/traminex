# Golden generation: RP6b seqgranularity (time-granularity coarsening).
# Anchors traminex's seqgranularity to TraMineRextras 0.6.9. Run in the live-R env:
#   Rscript scripts/r_goldens/make_goldens_seqgranularity.R
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
suppressMessages(library(TraMineRextras))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
stopifnot(as.character(packageVersion("TraMineRextras")) == "0.6.9")
out_dir <- file.path("tests", "r_goldens")

write_num_matrix <- function(name, m, header) {
  con <- file(file.path(out_dir, paste0(name, ".csv")), open = "wt")
  cat(header, file = con)
  txt <- apply(as.matrix(m), 1, function(r) paste(formatC(r, format = "g", digits = 15), collapse = ","))
  writeLines(txt, con)
  close(con)
}

A <- rbind(
  c("A", "A", "B", "C", "C", "C"),
  c("A", "B", "B", "B", "C", "A"),
  c("D", "D", "D", "D", "D", "D"),
  c("A", "B", "A", "B", "A", "B")
)
sdA <- seqdef(A, alphabet = c("A", "B", "C", "D"))
C <- rbind(
  c("A", "A", "B", "C", "C", "C"),
  c("A", "B", "B", "*", "*", "*"),
  c("A", "*", "A", "B", "*", "*"),
  c("C", "C", "*", "*", "*", "*")
)
sdC <- seqdef(C, alphabet = c("A", "B", "C"), missing = "*")

# 1-based code matrix of the coarsened object (void/missing -> NA via match against alphabet).
gran_mat <- function(sd, tspan, method) {
  res <- seqgranularity(sd, tspan = tspan, method = method)
  m <- as.matrix(res)
  matrix(match(m, alphabet(sd)), nrow = nrow(m))
}

for (meth in c("first", "last", "first.valid", "last.valid", "mostfreq")) {
  write_num_matrix(paste0("seqgran_A_", gsub("\\.", "", meth)), gran_mat(sdA, 3, meth),
    paste0("# R: seqgranularity(sdA, tspan=3, '", meth, "')\n"))
}
write_num_matrix("seqgran_C_last", gran_mat(sdC, 2, "last"),
  "# R: seqgranularity(sdC ragged, tspan=2, 'last'); void handling\n")
write_num_matrix("seqgran_C_mostfreq", gran_mat(sdC, 2, "mostfreq"),
  "# R: seqgranularity(sdC ragged, tspan=2, 'mostfreq')\n")

cat("SEQGRAN_GOLDENS_DONE\n")
