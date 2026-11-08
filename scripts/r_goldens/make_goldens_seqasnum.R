# Golden generation: seqasnum / seqnum (T11) -- numeric recoding.
#
# Pins traminex.seqasnum to TraMineR 2.2-14 seqasnum(with.missing): 0-based state
# codes (i-1 for the i-th alphabet state), NA for non-states (void always; nr
# unless with.missing, in which case nr takes the next code). NA is dumped as the
# sentinel -1 (matching the seqdss goldens). seqasnum is internal in TraMineR.
#
# Run from the repository root with the oracle library (scripts/r_goldens/README.md):
#   Rscript scripts/r_goldens/make_goldens_seqasnum.R
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
source(file.path("scripts", "r_goldens", "make_goldens_corestats.R"))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

asnum <- TraMineR:::seqasnum

emit <- function(sd, tag, wm) {
  m <- asnum(sd, with.missing = wm)
  m[is.na(m)] <- -1
  write_num_matrix(
    paste0("seqasnum_", tag), m,
    sprintf("# R: seqasnum(%s, with.missing=%s); NA (void / excluded nr) -> -1\n", tag, wm)
  )
}

fx <- corestats_fixtures()
emit(fx$A, "A", FALSE)
emit(fx$B, "B_wm0", FALSE)
emit(fx$B, "B_wm1", TRUE)
emit(fx$C, "C_wm0", FALSE)
emit(fx$C, "C_wm1", TRUE)
cat("SEQASNUM_GOLDENS_DONE\n")
