# Golden generation: seqerules -- association rules from frequent event
# subsequences with the FULL TraMineR statistics (seqerules.R + implicativestat.R):
# Support / Conf / Lift / Standardlift / JMeasure / ImplicStat / p-value /
# p.valueB1 / p.valueB2. Pinned on actcal.tse and the small synthetic TSE.
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

# Note: data.frame() name-mangles "p-value" -> "p.value".
NUMCOLS <- c("Support", "Conf", "Lift", "Standardlift", "JMeasure", "ImplicStat",
  "p.value", "p.valueB1", "p.valueB2")

dump_rules <- function(rules, tag) {
  writeLines(as.character(rules$Rules), file.path(out_dir, paste0("seqerules_", tag, "_rules.txt")))
  con <- file(file.path(out_dir, paste0("seqerules_", tag, ".csv")), "w")
  cat(paste0("# columns: ", paste(NUMCOLS, collapse = ","), "\n"), file = con)
  m <- as.matrix(rules[, NUMCOLS])
  for (i in seq_len(nrow(m))) {
    cat(paste(sprintf("%.15g", m[i, ]), collapse = ","), "\n", file = con, sep = "")
  }
  close(con)
}

seqerules <- TraMineR:::seqerules

data(actcal.tse, package = "TraMineR")
ea <- seqecreate(actcal.tse)
dump_rules(seqerules(seqefsub(ea, pmin.support = 0.02, max.k = 3)), "actcal")

tse <- data.frame(
  id    = c("s1", "s1", "s1", "s1", "s2", "s2", "s2", "s3", "s3", "s3"),
  time  = c(0,    1,    1,    3,    0,    2,    4,    0,    1,    2),
  event = c("A",  "B",  "C",  "A",  "A",  "B",  "C",  "A",  "B",  "C"),
  stringsAsFactors = FALSE
)
es <- seqecreate(tse)
dump_rules(seqerules(seqefsub(es, min.support = 1, max.k = 3)), "syn")

cat("SEQERULES_GOLDENS_DONE\n")
