# Live-R goldens for seqdist(method = "TWED") with the automatic indel (repair of
# audit finding TX-AUD-20261010-011).
#
# seqdist's default indel is "auto"; for TWED it resolves to 2 * max(sm) + nu + h
# over the whole substitution matrix, for an explicit matrix and for a generated
# cost method alike (seqdist.R:359-364, 433-439).
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_seqdist_twed_auto_indel.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

as_rows <- function(strings) {
  cells <- lapply(strings, function(s) strsplit(s, "")[[1]])
  width <- max(lengths(cells))
  as.data.frame(do.call(rbind, lapply(cells, function(x) c(x, rep(NA, width - length(x))))),
                stringsAsFactors = FALSE)
}
s4 <- suppressMessages(TraMineR::seqdef(as_rows(c("AABBC", "ABBCC", "CCBAA", "BBCCA")),
                                        alphabet = c("A", "B", "C")))
s5 <- suppressMessages(TraMineR::seqdef(as_rows(c("ABBCDA", "BBCCDD", "DACBBA", "AADDCC", "CBADCB")),
                                        alphabet = c("A", "B", "C", "D")))
s5r <- suppressMessages(TraMineR::seqdef(as_rows(c("ABBCDA", "BBCC", "DACBBA", "AAD")),
                                         alphabet = c("A", "B", "C", "D"), right = "DEL"))
sm5 <- abs(outer(1:4, 1:4, "-")) * 1.0  # double: checkcost's C helper refuses an integer matrix

cases <- list(
  s4_constant_auto = list(seqdata = s4, sm = "CONSTANT", nu = 0.5, h = 0.5),
  s5_matrix_auto = list(seqdata = s5, sm = sm5, nu = 0.3, h = 0.2),
  s5_trate_auto = list(seqdata = s5, sm = "TRATE", nu = 0.3, h = 0.2),
  s5_matrix_auto_yujianbo = list(seqdata = s5, sm = sm5, nu = 0.3, h = 0.2, norm = "YujianBo"),
  s5r_matrix_auto = list(seqdata = s5r, sm = sm5, nu = 0.3, h = 0.2)
)

rows <- list()
for (case_name in names(cases)) {
  args <- c(cases[[case_name]], list(method = "TWED", indel = "auto"))
  value <- as.matrix(suppressMessages(do.call(TraMineR::seqdist, args)))
  rows[[length(rows) + 1L]] <- data.frame(
    case = case_name,
    row = rep(seq_len(nrow(value)), times = ncol(value)),
    col = rep(seq_len(ncol(value)), each = nrow(value)),
    value = sprintf("%.15g", as.vector(value)),
    stringsAsFactors = FALSE
  )
}
write.csv(do.call(rbind, rows), file.path(out_dir, "seqdist_twed_auto_indel.csv"),
          row.names = FALSE, quote = TRUE)
