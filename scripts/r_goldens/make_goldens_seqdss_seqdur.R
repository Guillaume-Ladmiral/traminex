# Golden generation: seqdss / seqdur / seqivolatility (roadmap Phase 1).
#
# Anchors traminex's seqdss/seqdur/seqivolatility to TraMineR 2.2-14.
# Fixtures are void-free and rectangular so seqlength == ncol. NB: R seqdef()
# defaults to right="DEL", which converts a row's *trailing* missing to void and
# would break void-freeness (and disagree with traminex's seqdef, which keeps
# trailing missing as nr). We therefore pass right=NA on the missing fixture.
# Run from the repository root with the oracle library
# (scripts/r_goldens/README.md): Rscript scripts/r_goldens/make_goldens_seqdss_seqdur.R
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

# seqasnum is an internal TraMineR function (NOT in the package NAMESPACE), so it
# must be reached with the ::: operator; a bare seqasnum() call errors with
# "could not find function 'seqasnum'".
asnum <- TraMineR:::seqasnum

out_dir <- file.path("tests", "r_goldens")

write_int_matrix <- function(name, m, header) {
  con <- file(file.path(out_dir, paste0(name, ".csv")), open = "wt")
  cat(header, file = con)
  utils::write.table(m, con, sep = ",", row.names = FALSE, col.names = FALSE)
  close(con)
}

write_num_matrix <- function(name, m, header) {
  con <- file(file.path(out_dir, paste0(name, ".csv")), open = "wt")
  cat(header, file = con)
  # %.15g formatting to match the repo convention
  txt <- apply(m, 1, function(r) paste(formatC(r, format = "g", digits = 15), collapse = ","))
  writeLines(txt, con)
  close(con)
}

## ---- Fixture A (no missing) ----
A <- rbind(
  c("A","A","B","C","C","C"),
  c("A","B","B","B","C","A"),
  c("D","D","D","D","D","D"),
  c("A","B","A","B","A","B")
)
sdA <- seqdef(A, alphabet = c("A","B","C","D"))

dssA <- asnum(seqdss(sdA))               # 0-based codes, NA for right-pad
dssA[is.na(dssA)] <- -1
write_int_matrix("seqdss_A", dssA, "# R: seqasnum(seqdss(sdA)); NA->-1; alphabet A,B,C,D\n")

durA <- seqdur(sdA)
durA[is.na(durA)] <- -1
write_int_matrix("seqdur_A", durA, "# R: seqdur(sdA); NA->-1\n")

volA <- rbind(
  seqivolatility(sdA, w = 0.5, adjust = TRUE),
  seqivolatility(sdA, w = 0.5, adjust = FALSE),
  seqivolatility(sdA, w = 0.0, adjust = TRUE),   # isolates transp
  seqivolatility(sdA, w = 1.0, adjust = TRUE)    # isolates visit term
)
write_num_matrix("seqivol_A", volA,
  "# R: seqivolatility(sdA) rows=[w.5/adj, w.5/noadj, w0/adj, w1/adj]\n")

## ---- Fixture B (interior missing '*') ----
B <- rbind(
  c("A","A","*","B","B","C"),
  c("*","A","A","A","*","*"),
  c("A","*","A","B","B","B"),
  c("C","C","C","C","C","C")
)
sdB <- seqdef(B, alphabet = c("A","B","C"), missing = "*", right = NA)

dssB_nm <- asnum(seqdss(sdB, with.missing = FALSE), with.missing = FALSE)
dssB_nm[is.na(dssB_nm)] <- -1
write_int_matrix("seqdss_B_nomiss", dssB_nm,
  "# R: seqasnum(seqdss(sdB, with.missing=FALSE)); NA->-1; alphabet A,B,C\n")

dssB_m <- asnum(seqdss(sdB, with.missing = TRUE), with.missing = TRUE)
dssB_m[is.na(dssB_m)] <- -1
write_int_matrix("seqdss_B_miss", dssB_m,
  "# R: seqasnum(seqdss(sdB, with.missing=TRUE), with.missing=TRUE); NA->-1; missing code=3\n")

durB_nm <- seqdur(sdB, with.missing = FALSE); durB_nm[is.na(durB_nm)] <- -1
write_int_matrix("seqdur_B_nomiss", durB_nm, "# R: seqdur(sdB, with.missing=FALSE); NA->-1\n")

durB_m <- seqdur(sdB, with.missing = TRUE); durB_m[is.na(durB_m)] <- -1
write_int_matrix("seqdur_B_miss", durB_m, "# R: seqdur(sdB, with.missing=TRUE); NA->-1\n")

volB <- rbind(
  seqivolatility(sdB, w = 0.5, with.missing = FALSE, adjust = TRUE),
  seqivolatility(sdB, w = 0.5, with.missing = TRUE,  adjust = TRUE),
  seqivolatility(sdB, w = 0.5, with.missing = FALSE, adjust = FALSE)
)
write_num_matrix("seqivol_B", volB,
  "# R: seqivolatility(sdB) rows=[wm=F/adj, wm=T/adj, wm=F/noadj]\n")

## ---- Fixture C (ragged: trailing '*' -> void via default right="DEL") ----
## Short rows are padded with '*'; the DEFAULT right="DEL" turns the *trailing*
## '*' into void, while the interior '*' of C3 stays nr. (No right=NA here.)
C <- rbind(
  c("A","A","B","C","C","C"),
  c("A","B","B","*","*","*"),
  c("A","*","A","B","*","*"),
  c("C","C","*","*","*","*")
)
sdC <- seqdef(C, alphabet = c("A","B","C"), missing = "*")

dssC_nm <- asnum(seqdss(sdC, with.missing = FALSE), with.missing = FALSE)
dssC_nm[is.na(dssC_nm)] <- -1
write_int_matrix("seqdss_C_nomiss", dssC_nm,
  "# R: seqasnum(seqdss(sdC, with.missing=FALSE)); ragged/void; NA->-1; alphabet A,B,C\n")

dssC_m <- asnum(seqdss(sdC, with.missing = TRUE), with.missing = TRUE)
dssC_m[is.na(dssC_m)] <- -1
write_int_matrix("seqdss_C_miss", dssC_m,
  "# R: seqasnum(seqdss(sdC, with.missing=TRUE), with.missing=TRUE); ragged/void; NA->-1; missing code=3\n")

durC_nm <- seqdur(sdC, with.missing = FALSE); durC_nm[is.na(durC_nm)] <- -1
write_int_matrix("seqdur_C_nomiss", durC_nm, "# R: seqdur(sdC, with.missing=FALSE); ragged/void; NA->-1\n")

durC_m <- seqdur(sdC, with.missing = TRUE); durC_m[is.na(durC_m)] <- -1
write_int_matrix("seqdur_C_miss", durC_m, "# R: seqdur(sdC, with.missing=TRUE); ragged/void; NA->-1\n")

volC <- rbind(
  seqivolatility(sdC, w = 0.5, with.missing = FALSE, adjust = TRUE),
  seqivolatility(sdC, w = 0.5, with.missing = TRUE,  adjust = TRUE),
  seqivolatility(sdC, w = 0.5, with.missing = FALSE, adjust = FALSE)
)
write_num_matrix("seqivol_C", volC,
  "# R: seqivolatility(sdC) ragged/void rows=[wm=F/adj, wm=T/adj, wm=F/noadj]\n")

cat("SEQDSS_SEQDUR_GOLDENS_DONE\n")
