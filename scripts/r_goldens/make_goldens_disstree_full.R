# Golden generation: disstree leaf helpers (RP2 / B-5).
#
# Anchors traminex's R-faithful leaf-membership / classification-rule helpers to
# TraMineR 2.2-14 `disstreeleaf(label=TRUE)`, `disstree.get.rules`,
# `disstree.assign`. Reuses the shared Phase-2 fixture and the same deterministic
# tree as make_goldens_disstree.R (R=1 -> no permutation gate -> deterministic).
#
# Outputs:
#   disstree_full_leaflabels.txt  one human-readable leaf label per case (15 rows)
#   disstree_full_rules.txt       the classification rules (R syntax), one per line
#   disstree_full_assign.csv      disstree.assign(rules, pred): 1-based leaf index per case
#
# Requires R + TraMineR 2.2-14; see memory "live-traminer-r-env".
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
source(file.path("scripts", "r_goldens", "make_goldens_phase2_inputs.R"))

inp <- phase2_inputs()
D <- inp$D
i <- 1:inp$n
pred <- data.frame(
  A = factor(c("a", "a", "b", "b", "a", "b", "b", "b", "c", "c", "c", "c", "c", "a", "c")),
  X = ((i * 7) %% 15) + 0.1 * i
)

tr <- disstree(D ~ A + X, data = pred, R = 1, min.size = 2, max.depth = 4)

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

leaf_lab <- as.character(disstreeleaf(tr, label = TRUE, collapse = "; "))
rules <- disstree.get.rules(tr, collapse = "; ")
assign_id <- disstree.assign(rules, pred)

writeLines(leaf_lab, file.path(out_dir, "disstree_full_leaflabels.txt"))
writeLines(as.character(rules), file.path(out_dir, "disstree_full_rules.txt"))
write_phase2_csv(
  "disstree_full_assign", matrix(as.numeric(assign_id), ncol = 1),
  "# R: disstree.assign(disstree.get.rules(tr), pred) -> 1-based leaf rule index per case\n"
)

cat("=== leaf labels (per case) ===\n")
print(leaf_lab)
cat("=== rules (levels order) ===\n")
print(rules)
cat("=== assign (per case) ===\n")
print(assign_id)
cat("DISSTREE_FULL_GOLDENS_DONE\n")
