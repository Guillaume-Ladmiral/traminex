# Golden generation: RP6b seqidegrad (degradation index) + seqinsecurity (type=2).
# Anchors traminex's seqidegrad / seqinsecurity to TraMineR 2.2-14. Reuses the
# precarity A/B fixtures. Run in the live-R env (memory "live-traminer-r-env"):
#   Rscript scripts/r_goldens/make_goldens_idegrad.R
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

A <- rbind(
  c("A", "A", "B", "C", "C", "C"),
  c("A", "B", "B", "B", "C", "A"),
  c("D", "D", "D", "D", "D", "D"),
  c("A", "B", "A", "B", "A", "B")
)
sdA <- seqdef(A, alphabet = c("A", "B", "C", "D"))
B <- rbind(
  c("A", "A", "*", "B", "B", "C"),
  c("*", "A", "A", "A", "*", "*"),
  c("A", "*", "A", "B", "B", "B"),
  c("C", "C", "C", "C", "C", "C")
)
sdB <- seqdef(B, alphabet = c("A", "B", "C"), missing = "*", right = NA)

## ---- seqidegrad (Degrad column) ----
write_num_matrix("seqidegrad_A",          as.vector(seqidegrad(sdA)),                       "# R: seqidegrad(sdA); RANK, spell.integr=T\n")
write_num_matrix("seqidegrad_A_tratedss", as.vector(seqidegrad(sdA, method = "TRATEDSS")),  "# R: seqidegrad(sdA, method=TRATEDSS)\n")
write_num_matrix("seqidegrad_A_nospell",  as.vector(seqidegrad(sdA, pow = FALSE)),          "# R: seqidegrad(sdA, pow=FALSE); RANK, spell.integr=F\n")
write_num_matrix("seqidegrad_B_wm",       as.vector(seqidegrad(sdB, with.missing = TRUE)),  "# R: seqidegrad(sdB, wm=T)\n")

## ---- seqinsecurity (type=2; commit 3) ----
write_num_matrix("seqinsecurity_A",        as.vector(seqinsecurity(sdA)),                       "# R: seqinsecurity(sdA); type2 RANK, start/spell.integr=T\n")
write_num_matrix("seqinsecurity_A_bound",  as.vector(seqinsecurity(sdA, bound = TRUE)),         "# R: seqinsecurity(sdA, bound=T)\n")
write_num_matrix("seqinsecurity_B_wm",     as.vector(seqinsecurity(sdB, with.missing = TRUE)),  "# R: seqinsecurity(sdB, wm=T)\n")

## ---- WEIGHTED rate methods (audit F4): seqtrate() inside FREQ/TRATE/TRATEDSS is
## weighted=TRUE, so the transition counts must use the sequence weights. The
## unweighted goldens above cannot see this; use a weighted seqdef of fixture A.
wA <- c(1.5, 2, 0.5, 3)
sdAw <- seqdef(A, alphabet = c("A", "B", "C", "D"), weights = wA)
write_num_matrix("seqidegrad_Aw_tratedss",    as.vector(seqidegrad(sdAw, method = "TRATEDSS")), "# R: seqidegrad(sdAw weighted, method=TRATEDSS)\n")
write_num_matrix("seqidegrad_Aw_trate",       as.vector(seqidegrad(sdAw, method = "TRATE")),    "# R: seqidegrad(sdAw weighted, method=TRATE)\n")
write_num_matrix("seqidegrad_Aw_freq",        as.vector(seqidegrad(sdAw, method = "FREQ")),     "# R: seqidegrad(sdAw weighted, method=FREQ)\n")
write_num_matrix("seqinsecurity_Aw_tratedss", as.vector(seqinsecurity(sdAw, method = "TRATEDSS")), "# R: seqinsecurity(sdAw weighted, method=TRATEDSS)\n")

cat("IDEGRAD_GOLDENS_DONE\n")
