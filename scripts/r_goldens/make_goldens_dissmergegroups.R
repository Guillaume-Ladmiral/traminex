# Live-R goldens for dissmergegroups (RP4, TraMineR + WeightedCluster).
#
# Greedy ASW-quality group merging on the shared phase2 15x15 distance matrix,
# starting from a 6-cluster Ward cut.
#
# Oracle verified on Linux with R 4.6.2, TraMineR 2.2-14, and
# WeightedCluster 2.0. The generator must run against those exact package
# authorities in the isolated R library used by the golden-generation loop.

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
source("scripts/r_goldens/make_goldens_phase2_inputs.R")  # phase2_inputs(), write_phase2_csv()
suppressMessages(library(TraMineR))
suppressMessages(library(WeightedCluster))

if (getRversion() != "4.6.2" ||
    packageVersion("TraMineR") != "2.2.14" ||
    packageVersion("WeightedCluster") != "2.0") {
  stop("dissmergegroups goldens require R 4.6.2, TraMineR 2.2-14, and WeightedCluster 2.0")
}

inp <- phase2_inputs()
D <- inp$D
w <- inp$wreal

hc <- hclust(as.dist(D), method = "ward.D2")
g0 <- cutree(hc, 6)
write_phase2_csv("dissmerge_input_group", matrix(g0, nrow = 1), "# initial grouping (ward.D2 cutree 6)\n")

# Direct wcClusterQuality ASW pins for the initial grouping.
write_phase2_csv(
  "dissmerge_asw_uw", matrix(wcClusterQuality(as.dist(D), g0)$stats["ASW"], nrow = 1),
  "# wcClusterQuality ASW, unweighted\n"
)
write_phase2_csv(
  "dissmerge_asw_w", matrix(wcClusterQuality(as.dist(D), g0, weights = w)$stats["ASW"], nrow = 1),
  "# wcClusterQuality ASW, weighted\n"
)

gen <- function(tag, weights) {
  mg <- suppressMessages(
    dissmergegroups(D, group = g0, weights = weights, min.group = 2, small = .05, silent = TRUE)
  )
  write_phase2_csv(
    paste0("dissmerge_", tag, "_group"), matrix(as.vector(mg), nrow = 1), "# merged labels (1..K)\n"
  )
  write_phase2_csv(
    paste0("dissmerge_", tag, "_qual"), matrix(attr(mg, "qual"), nrow = 1), "# final ASW quality\n"
  )
}

gen("uw", NULL)
gen("w", w)

cat("DISSMERGE_GOLDENS_DONE\n")
