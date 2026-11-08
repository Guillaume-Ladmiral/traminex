# Golden generation: seqecreate must honor a *factor* event column's LEVEL ORDER
# as the event dictionary (R keeps factor levels; it does NOT sort them).
#
# This is the RP5b precursor fix: the event-mining goldens (seqefsub etc.) are
# generated on actcal.tse, whose event column is a factor with NON-alphabetical
# levels (PartTime, NoActivity, Start, ...). The factor codes drive (a) the
# comma-group ordering of simultaneous events and (b) Count-tie ordering in the
# miner, so traminex's seqecreate must use the same dictionary order.
#
# We pin two things:
#   1. a small synthetic factor TSE with non-alphabetical levels and simultaneous
#      events (sorted by factor code, NOT alphabetically) -- self-contained;
#   2. levels(seqecreate(actcal.tse)) -- the real 8-level factor order.
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

# (1) Synthetic factor TSE: levels Gamma<Beta<Alpha (NOT alphabetical).
# p1 has Gamma & Alpha at the same time -> comma group by factor code -> (Gamma,Alpha).
tse <- data.frame(
  id    = c("p1", "p1", "p1", "p2", "p2", "p3", "p3"),
  time  = c(0,    2,    2,    1,    3,    0,    0),
  event = factor(
    c("Beta", "Gamma", "Alpha", "Beta", "Alpha", "Gamma", "Alpha"),
    levels = c("Gamma", "Beta", "Alpha")
  ),
  stringsAsFactors = FALSE
)
e <- seqecreate(tse)
writeLines(as.character(e), file.path(out_dir, "eseq_factor_aschar.txt"))
writeLines(levels(e), file.path(out_dir, "eseq_factor_levels.txt"))

# (2) actcal.tse real factor level order.
data(actcal.tse, package = "TraMineR")
ea <- seqecreate(actcal.tse)
writeLines(levels(ea), file.path(out_dir, "eseq_actcal_levels.txt"))

cat("EVENTS_FACTOR_GOLDENS_DONE\n")
print(as.character(e))
cat("levels:", paste(levels(e), collapse = ","), "\n")
cat("actcal levels:", paste(levels(ea), collapse = ","), "\n")
