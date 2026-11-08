# Live-R goldens for gower_matrix (RP4 T15).
#
# TraMineR's gower_matrix (internal, dissmfacw.R:11) double-centres a distance
# matrix with optional squared/weights, on the shared phase2 15x15 matrix.
#
# Requires R 4.6.2 + TraMineR 2.2-14; see scripts/r_goldens/README.md.

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
source("scripts/r_goldens/make_goldens_phase2_inputs.R")  # phase2_inputs(), write_phase2_csv()
suppressMessages(library(TraMineR))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

inp <- phase2_inputs()
D <- inp$D
w <- inp$wreal

gen <- function(tag, squared, weights) {
  m <- TraMineR:::gower_matrix(D, squared = squared, weights = weights)
  wlab <- if (is.null(weights)) "NULL" else "wreal"
  write_phase2_csv(
    paste0("gower_", tag), m,
    paste0("# TraMineR:::gower_matrix(D, squared=", squared, ", weights=", wlab, ")\n")
  )
}

gen("sq_uw", TRUE, NULL)
gen("raw_uw", FALSE, NULL)
gen("sq_w", TRUE, w)
gen("raw_w", FALSE, w)

cat("GOWER_GOLDENS_DONE\n")
