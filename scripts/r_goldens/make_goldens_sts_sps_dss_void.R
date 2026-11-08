# Live-R goldens for void cells in the SPS and DSS forms, and for the ids
# R's STS_to_TSE gives by default (repair of audit finding TX-AUD-20261010-033).
#
# seqformat(to = "SPS") and seqformat(to = "DSS") turn void cells into NA, so a
# shorter sequence ends early; a missing cell stays a state (seqformat.R:240-253).
# STS_to_TSE, defined but not exported (called through ::: with its tevent),
# numbers the sequences 1..n when id is NULL (seqformat-STS_to_TSE.R:14).
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_sts_sps_dss_void.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rag <- suppressMessages(TraMineR::seqdef(data.frame(a = c("A", "A", "B"), b = c("B", NA, "B"), c = c(NA, "B", NA),
                                                    d = c(NA, "B", NA)), alphabet = c("A", "B")))
sps <- suppressMessages(TraMineR::seqformat(rag, from = "STS", to = "SPS"))
dss <- suppressMessages(TraMineR::seqformat(rag, from = "STS", to = "DSS"))
rows <- list()
for (i in seq_len(nrow(rag))) {
  rows[[length(rows) + 1L]] <- data.frame(form = "SPS", sequence = i, position = seq_len(ncol(sps)),
                                          token = ifelse(is.na(as.vector(sps[i, ])), "<NA>", as.character(sps[i, ])))
  rows[[length(rows) + 1L]] <- data.frame(form = "DSS", sequence = i, position = seq_len(ncol(dss)),
                                          token = ifelse(is.na(as.vector(as.matrix(dss)[i, ])), "<NA>", as.character(as.matrix(dss)[i, ])))
}
ids <- suppressMessages(TraMineR::seqdef(data.frame(a = c("A", "A"), b = c("A", "B"), c = c("B", "B")), id = c("p1", "p2")))
tev <- suppressMessages(TraMineR::seqetm(ids, method = "state"))
tse <- TraMineR:::STS_to_TSE(ids, NULL, tev)
rows[[length(rows) + 1L]] <- data.frame(form = "TSE_id_default", sequence = 0L, position = seq_len(nrow(tse)),
                                        token = as.character(tse$id))
write.csv(do.call(rbind, rows), file.path(out_dir, "sts_sps_dss_void.csv"), row.names = FALSE, quote = TRUE)
