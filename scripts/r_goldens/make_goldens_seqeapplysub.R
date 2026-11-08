# Golden generation: seqeapplysub -- the (n_seq x n_subseq) subsequence-count
# matrix (TraMineR C_tmrmatrixsubseqinseq). This is where the SIX count methods
# (COBJ/CDIST_O/CWIN/CMINWIN/CDIST + age) are exercised end to end.
#
# We mine a small synthetic TSE (reconstructible in Python) with seqefsub, then
# apply each method and dump the raw matrix (rows = sequences in seqecreate order,
# cols = subsequences in seqefsub Count-desc order). CWIN/CMINWIN are applied via
# a window.size constraint.
#
# Requires live R + TraMineR 2.2-14; see memory "live-traminer-r-env".
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

tse <- data.frame(
  id    = c("s1", "s1", "s1", "s1", "s2", "s2", "s2", "s3", "s3", "s3"),
  time  = c(0,    1,    1,    3,    0,    2,    4,    0,    1,    2),
  event = c("A",  "B",  "C",  "A",  "A",  "B",  "C",  "A",  "B",  "C"),
  stringsAsFactors = FALSE
)
es <- seqecreate(tse)
fsub <- seqefsub(es, min.support = 1, max.k = 3)

dump_mat <- function(m, tag) {
  con <- file(file.path(out_dir, paste0("seqeapplysub_", tag, ".csv")), "w")
  cat(sprintf("# %d x %d matrix (rows=sequences, cols=subsequences)\n", nrow(m), ncol(m)),
    file = con)
  for (i in seq_len(nrow(m))) {
    cat(paste(sprintf("%.15g", m[i, ]), collapse = ","), "\n", file = con, sep = "")
  }
  close(con)
}

dump_mat(seqeapplysub(fsub, method = "presence"), "presence")
dump_mat(seqeapplysub(fsub, method = "count"), "count")
dump_mat(seqeapplysub(fsub, method = "age"), "age")
dump_mat(seqeapplysub(fsub, method = "CDIST"), "cdist")
dump_mat(seqeapplysub(fsub, constraint = seqeconstraint(window.size = 3, count.method = "CWIN")),
  "cwin_w3")
dump_mat(seqeapplysub(fsub, constraint = seqeconstraint(window.size = 3, count.method = "CMINWIN")),
  "cminwin_w3")
dump_mat(seqeapplysub(fsub, method = "presence", rules = TRUE), "rules_presence")

# Fractional timestamps expose the age-method integer truncation: C++
# Sequence::count accumulates the first-occurrence age into an `int c`, so e.g.
# an age of 1.9999 is truncated (not rounded) to 1. The integer-time fixture
# above cannot see this.
tsef <- data.frame(
  id    = c("p1", "p1", "p2", "p2", "p3", "p3"),
  time  = c(0.7,  1.5,  2.9,  5.0,  1.9999, 2.0),
  event = c("A",  "B",  "A",  "B",  "A",    "B"),
  stringsAsFactors = FALSE
)
ef <- seqecreate(tsef)
fsubf <- seqefsub(ef, min.support = 1, max.k = 2)
dump_mat(seqeapplysub(fsubf, method = "age"), "age_frac")
writeLines(as.character(fsubf$subseq), file.path(out_dir, "seqeapplysub_age_frac_cols.txt"))

cat("SEQEAPPLYSUB_GOLDENS_DONE\n")
cat("nsub:", length(fsub$subseq), " nseq:", length(es), "\n")
