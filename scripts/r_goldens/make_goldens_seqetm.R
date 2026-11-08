# Golden generation: seqetm (TraMineR) -- build the n x n transition/event matrix
# from a state-sequence object's alphabet (+ optional state labels).
#
# seqetm(seqdata, method, use.labels, sep, bp, ep) returns a character matrix
# with rownames/colnames = the alphabet; the diagonal is the (labelled) state and
# the off-diagonal cells encode method="transition" ("A>B"), "period"
# ("endA,beginB") or "state" ("B"). We pin all three methods, with and without
# state labels, as exact string matrices.
#
# Requires live R + TraMineR 2.2-14; see memory "live-traminer-r-env".
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

# Small fixed stslist: 3 states A/B/C with long labels, 4 short sequences.
mat <- rbind(
  c("A", "A", "B", "C"),
  c("B", "B", "B", "A"),
  c("C", "A", "A", "B"),
  c("A", "B", "C", "C")
)
sd <- seqdef(mat, alphabet = c("A", "B", "C"),
             labels = c("Apple", "Banana", "Cherry"))

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dump_tm <- function(tm, fname, call) {
  con <- file(file.path(out_dir, fname), "w")
  cat(sprintf("# R: %s\n", call), file = con)
  # header: tab-joined colnames
  cat(paste(colnames(tm), collapse = "\t"), "\n", file = con, sep = "")
  for (i in seq_len(nrow(tm))) {
    cat(rownames(tm)[i], "\t", paste(tm[i, ], collapse = "\t"), "\n",
        file = con, sep = "")
  }
  close(con)
}

dump_tm(seqetm(sd, method = "transition", use.labels = TRUE),
        "seqetm_transition.tsv", "seqetm(sd, method='transition', use.labels=TRUE)")
dump_tm(seqetm(sd, method = "period", use.labels = TRUE),
        "seqetm_period.tsv", "seqetm(sd, method='period', use.labels=TRUE)")
dump_tm(seqetm(sd, method = "state", use.labels = TRUE),
        "seqetm_state.tsv", "seqetm(sd, method='state', use.labels=TRUE)")
dump_tm(seqetm(sd, method = "transition", use.labels = FALSE),
        "seqetm_transition_nolab.tsv", "seqetm(sd, method='transition', use.labels=FALSE)")

# Data carrying a missing (nr '*') AND a void ('%') state: R's seqetm appends
# each ONCE -> statl = A,B,C,*,% (three real states plus two markers). Interior '*' stays nr; the
# trailing '*' in row 2 is converted to void '%' by seqdef(right="DEL").
matm <- rbind(
  c("A", "A", "*", "B"),
  c("B", "B", "A", "*"),
  c("C", "A", "A", "B")
)
sdm <- seqdef(matm, alphabet = c("A", "B", "C"), missing = "*")
dump_tm(seqetm(sdm, method = "transition", use.labels = FALSE),
        "seqetm_transition_miss.tsv", "seqetm(sdm, method='transition', use.labels=FALSE)")
dump_tm(seqetm(sdm, method = "state", use.labels = FALSE),
        "seqetm_state_miss.tsv", "seqetm(sdm, method='state', use.labels=FALSE)")

cat("SEQETM_GOLDENS_DONE\n")
