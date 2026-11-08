# Golden generation: seqrecode (T17 VERIFY) vs TraMineR 2.2-14.
#
# seqrecode(seqdata, recodes, otherwise) collapses the alphabet into a new one.
# New alphabet order = names(recodes) then, if otherwise is NULL, the unlisted
# original states in their original order (else `otherwise`). We pin the recoded
# sequence codes (seqasnum, 0-based over the new alphabet) and the new alphabet.
#
# Run from the repository root with the oracle library (scripts/r_goldens/README.md):
#   Rscript scripts/r_goldens/make_goldens_seqrecode.R
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
source(file.path("scripts", "r_goldens", "make_goldens_corestats.R"))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")

emit_recode <- function(name, res) {
  con <- file(file.path(out_dir, paste0(name, "_alph.txt")), open = "wt")
  writeLines(paste(alphabet(res), collapse = " "), con)
  close(con)
  codes <- TraMineR:::seqasnum(res)
  codes[is.na(codes)] <- -1
  write_num_matrix(paste0(name, "_codes"), codes, sprintf("# R: seqasnum(%s)\n", name))
}

fx <- corestats_fixtures()

# Scheme 1: collapse A,B -> Work; C -> Leisure; D unlisted (otherwise=NULL, kept).
r1 <- seqrecode(fx$A, recodes = list(Work = c("A", "B"), Leisure = "C"))
emit_recode("seqrecode_keep", r1)

# Scheme 2: A,B,C -> Active; D -> Inactive (otherwise).
r2 <- seqrecode(fx$A, recodes = list(Active = c("A", "B", "C")), otherwise = "Inactive")
emit_recode("seqrecode_otherwise", r2)

cat("SEQRECODE_GOLDENS_DONE\n")
