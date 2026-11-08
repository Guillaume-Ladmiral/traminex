# Golden generation: RP6b seqsha (sequence-history person-period converter).
# Anchors traminex's seqsha to TraMineRextras 0.6.9. Run in the live-R env:
#   Rscript scripts/r_goldens/make_goldens_seqsha.R
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
ev <- c(TRUE, FALSE, TRUE, TRUE)

# Numeric block: id, time, event(0/1), then past-state columns as 1-based codes (NA kept).
sha_mat <- function(sd, time, event, ...) {
  res <- seqsha(sd, time = time, event = event, ...)
  alph <- alphabet(sd)
  past_cols <- setdiff(colnames(res), c("id", "time", "event"))
  pm <- sapply(res[past_cols], function(col) match(as.character(col), alph))
  cbind(res$id, res$time, as.integer(res$event), pm)
}

write_num_matrix("seqsha_A_start", sha_mat(sdA, c(3, 6, 2, 4), ev),
  "# R: seqsha(sdA, time=c(3,6,2,4), event); start-align; id,time,event,T1..T6 codes\n")
write_num_matrix("seqsha_A_end", sha_mat(sdA, c(3, 4, 2, 5), ev, align.end = TRUE),
  "# R: seqsha(sdA, time=c(3,4,2,5), align.end=T); Tm6..Tm1 codes\n")
write_num_matrix("seqsha_A_end_present", sha_mat(sdA, c(3, 4, 2, 5), ev, align.end = TRUE, include.present = TRUE),
  "# R: seqsha(sdA, align.end=T, include.present=T)\n")

cat("SEQSHA_GOLDENS_DONE\n")
