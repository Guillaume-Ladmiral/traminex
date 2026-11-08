# Golden generation: RP6b seqcta (competing-trajectory-analysis converter).
# Anchors traminex's seqcta to TraMineRextras 0.6.9. Run in the live-R env:
#   Rscript scripts/r_goldens/make_goldens_seqcta.R
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

# Numeric block: id, time, event(0/1), lastobs(0/1), then T1..Tsubseq as 1-based codes.
cta_mat <- function(res, sd, subseq) {
  alph <- alphabet(sd)
  tcols <- paste0("T", 1:subseq)
  pm <- sapply(res[tcols], function(col) match(as.character(col), alph))
  cbind(res$id, res$time, as.integer(res$event), as.integer(res$lastobs), pm)
}

r1 <- seqcta(sdA, subseq = 2, time = c(4, 3, 5, 2), event = c(TRUE, FALSE, TRUE, TRUE))
write_num_matrix("seqcta_A_given", cta_mat(r1, sdA, 2),
  "# R: seqcta(sdA, subseq=2, time=c(4,3,5,2), event); id,time,event,lastobs,T1,T2\n")

r2 <- seqcta(sdA, subseq = 2, initial.state = "A")
write_num_matrix("seqcta_A_init", cta_mat(r2, sdA, 2),
  "# R: seqcta(sdA, subseq=2, initial.state='A')\n")

cat("SEQCTA_GOLDENS_DONE\n")
