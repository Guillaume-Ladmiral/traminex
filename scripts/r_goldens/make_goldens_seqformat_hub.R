# Golden generation: RP6b faithful seqformat STS->SPELL / STS->SRS hub (closes T18).
# Anchors traminex's STS_to_SPELL / STS_to_SRS to TraMineR 2.2-14 seqformat().
# Run from the repository root with the oracle library (scripts/r_goldens/README.md):
#   Rscript scripts/r_goldens/make_goldens_seqformat_hub.R
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

## ---- Fixtures (A void/nr-free; B interior nr, full length; C ragged void) ----
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
C <- rbind(
  c("A", "A", "B", "C", "C", "C"),
  c("A", "B", "B", "*", "*", "*"),
  c("A", "*", "A", "B", "*", "*"),
  c("C", "C", "*", "*", "*", "*")
)
sdC <- seqdef(C, alphabet = c("A", "B", "C"), missing = "*")

## ---- STS -> SPELL : columns id, begin, end, state(1-based factor code) ----
spell_mat <- function(sd, wm) {
  sp <- suppressMessages(seqformat(sd, from = "STS", to = "SPELL", with.missing = wm))
  cbind(
    as.numeric(as.character(sp$id)),
    as.numeric(sp$begin),
    as.numeric(sp$end),
    as.integer(sp$states)
  )
}
write_num_matrix("seqformat_spell_A",    spell_mat(sdA, TRUE),  "# R: seqformat(sdA,STS,SPELL,wm=T); cols id,begin,end,state(1-based)\n")
write_num_matrix("seqformat_spell_B_wm", spell_mat(sdB, TRUE),  "# R: seqformat(sdB,STS,SPELL,wm=T); nr kept as state; levels A,B,C,*\n")
write_num_matrix("seqformat_spell_B_nm", spell_mat(sdB, FALSE), "# R: seqformat(sdB,STS,SPELL,wm=F); nr->void; levels A,B,C\n")
write_num_matrix("seqformat_spell_C_wm", spell_mat(sdC, TRUE),  "# R: seqformat(sdC,STS,SPELL,wm=T); ragged void; levels A,B,C,*\n")

## Interior void (audit F5): R's seqdef rejects a literal void '%' in data, but
## seqformat internally unifies nr+void -> NA and re-derives interior missing as
## nr (kept) / trailing missing as void (dropped). Build that reference via an
## interior nr '*' + a trailing '*' (-> void by right="DEL"). traminex's
## interior-VOID object (A A % B C *) must produce the SAME spell output.
E <- rbind(c("A", "A", "*", "B", "C", "*"))
sdE <- seqdef(E, alphabet = c("A", "B", "C"), missing = "*")
write_num_matrix("seqformat_spell_E_wm", spell_mat(sdE, TRUE),  "# R: interior missing kept as nr '*' (pos3); trailing ->void dropped; levels A,B,C,*\n")
write_num_matrix("seqformat_spell_E_nm", spell_mat(sdE, FALSE), "# R: interior nr->void dropped; levels A,B,C\n")

## ---- STS -> SRS : id, idx, then T-(L-1)..T-1,T as 1-based state codes (NA kept) ----
srs_mat <- function(sd) {
  srs <- suppressMessages(seqformat(sd, from = "STS", to = "SRS"))
  alph <- alphabet(sd)
  nc <- ncol(srs)
  out <- matrix(NA_real_, nrow(srs), nc)
  out[, 1] <- as.numeric(srs$id)
  out[, 2] <- as.numeric(srs$idx)
  for (j in 3:nc) out[, j] <- match(as.character(srs[[j]]), alph)
  out
}
write_num_matrix("seqformat_srs_A", srs_mat(sdA), "# R: seqformat(sdA,STS,SRS); id,idx,T-5..T-1,T as 1-based codes (NA=leading)\n")

cat("SEQFORMAT_HUB_GOLDENS_DONE\n")
