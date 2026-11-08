# Golden generation: RP6b seqindic.dyn (dynamic/windowed sequence indicator).
# Anchors traminex's seqindic_dyn to TraMineRextras 0.6.9. Run in the live-R env:
#   Rscript scripts/r_goldens/make_goldens_seqindic_dyn.R
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

dump <- function(name, dyn, header) {
  m <- matrix(as.numeric(dyn), nrow = nrow(dyn), ncol = ncol(dyn))
  write_num_matrix(name, m, header)
}

dump("seqindic_dyn_A_cplx_slide", seqindic.dyn(sdA, indic = "cplx", window.size = 0.5),
  "# R: seqindic.dyn(sdA, cplx, ws=0.5); sliding step=3; 4 windows\n")
dump("seqindic_dyn_A_cplx_incr", seqindic.dyn(sdA, indic = "cplx", window.size = 0.5, sliding = FALSE),
  "# R: seqindic.dyn(sdA, cplx, ws=0.5, sliding=F); incremental\n")
dump("seqindic_dyn_A_trans_slide", seqindic.dyn(sdA, indic = "trans", window.size = 2),
  "# R: seqindic.dyn(sdA, trans, ws=2); sliding step=2\n")
dump("seqindic_dyn_C_cplx_slide", seqindic.dyn(sdC, indic = "cplx", window.size = 0.5),
  "# R: seqindic.dyn(sdC ragged, cplx, ws=0.5); NA where window-end is void\n")

cat("SEQINDIC_DYN_GOLDENS_DONE\n")
