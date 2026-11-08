# Golden generation: seqST (Elzinga turbulence) void/missing faithfulness.
#
# Pins traminex.seqST to TraMineR 2.2-14 across rectangular void-free (A, B) and
# ragged void-padded (C) fixtures, for norm in {FALSE, TRUE} x with.missing in
# {FALSE, TRUE}, type=1. This anchors the fix for the void-counting defect under
# with.missing=TRUE (seqST counted structural void as a state) and the norm
# reference alphabet (it must exclude void and exclude nr unless with.missing).
#
# Reuses the Phase 1 fixtures A/B/C (see make_goldens_seqdss_seqdur.R). Fixture B
# is built with right=NA (void-free, trailing '*' kept as nr); fixture C uses the
# default right="DEL" so trailing '*' become void (ragged).
#
# Run from the repository root with the oracle library (scripts/r_goldens/README.md):
#   Rscript scripts/r_goldens/make_goldens_seqst.R
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
  txt <- apply(m, 1, function(r) paste(formatC(r, format = "g", digits = 15), collapse = ","))
  writeLines(txt, con)
  close(con)
}

st_rows <- function(sd, type = 1) {
  rbind(
    as.vector(seqST(sd, norm = FALSE, with.missing = FALSE, type = type)),
    as.vector(seqST(sd, norm = FALSE, with.missing = TRUE, type = type)),
    as.vector(seqST(sd, norm = TRUE, with.missing = FALSE, type = type)),
    as.vector(seqST(sd, norm = TRUE, with.missing = TRUE, type = type))
  )
}

HEADER <- "# R: seqST rows=[norm=F/wm=F, norm=F/wm=T, norm=T/wm=F, norm=T/wm=T], type=1\n"
HEADER2 <- "# R: seqST rows=[norm=F/wm=F, norm=F/wm=T, norm=T/wm=F, norm=T/wm=T], type=2\n"

## ---- Fixture A (no missing) ----
A <- rbind(
  c("A", "A", "B", "C", "C", "C"),
  c("A", "B", "B", "B", "C", "A"),
  c("D", "D", "D", "D", "D", "D"),
  c("A", "B", "A", "B", "A", "B")
)
sdA <- seqdef(A, alphabet = c("A", "B", "C", "D"))
write_num_matrix("seqST_A", st_rows(sdA), HEADER)
write_num_matrix("seqST_A_t2", st_rows(sdA, type = 2), HEADER2)

## ---- Fixture B (interior missing '*', void-free via right=NA) ----
B <- rbind(
  c("A", "A", "*", "B", "B", "C"),
  c("*", "A", "A", "A", "*", "*"),
  c("A", "*", "A", "B", "B", "B"),
  c("C", "C", "C", "C", "C", "C")
)
sdB <- seqdef(B, alphabet = c("A", "B", "C"), missing = "*", right = NA)
write_num_matrix("seqST_B", st_rows(sdB), HEADER)
write_num_matrix("seqST_B_t2", st_rows(sdB, type = 2), HEADER2)

## ---- Fixture C (ragged: trailing '*' -> void via default right="DEL") ----
C <- rbind(
  c("A", "A", "B", "C", "C", "C"),
  c("A", "B", "B", "*", "*", "*"),
  c("A", "*", "A", "B", "*", "*"),
  c("C", "C", "*", "*", "*", "*")
)
sdC <- seqdef(C, alphabet = c("A", "B", "C"), missing = "*")
write_num_matrix("seqST_C", st_rows(sdC), HEADER)
write_num_matrix("seqST_C_t2", st_rows(sdC, type = 2), HEADER2)

cat("SEQST_GOLDENS_DONE\n")
