# Golden generation: the state space CHI2/EUCLID and the NMS family are sized on.
#
# R sizes global.pdotj and prox by the alphabet (plus the missing state when
# with.missing), never by the codes a corpus happens to use.  These fixtures
# pin "obs" and explicit vectors/matrices on an object whose last state is
# unobserved, on a ragged (void-bearing) object, and on the analogue of a plain
# integer matrix (seqdef with alphabet 0:max); plus R's refusals of vectors and
# matrices sized to the observed codes only (2026-09-22 audit S1a-01, S1a-04,
# S1b-04, S1b-05).
#
# Consumed by tests/test_seqdist_state_space.py.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

X <- matrix(c(0, 1, 2, 1, 1, 1, 0, 2, 2, 0, 1, 1, 0, 2, 2, 0, 1, 0, 0, 1), ncol = 4, byrow = TRUE)
plain_analogue <- suppressMessages(seqdef(X, alphabet = 0:2))
unobserved_last <- suppressMessages(seqdef(X, alphabet = 0:3))
Xr <- X
Xr[1, 4] <- NA
Xr[3, 3:4] <- NA
ragged <- suppressMessages(seqdef(Xr, alphabet = 0:2, right = "DEL"))
Xm <- X
Xm[2, 3] <- NA
with_missing <- suppressMessages(seqdef(Xm, alphabet = 0:2))
letters_dense <- suppressMessages(seqdef(
  data.frame(a = c("A", "A", "C"), b = c("A", "B", "C"), c = c("B", "B", "A"), d = c("C", "C", "B")),
  alphabet = c("A", "B", "C", "D")
))
P4 <- matrix(c(1, .5, 0, .2, .5, 1, .3, .2, 0, .3, 1, .2, .2, .2, .2, 1), 4, byrow = TRUE)

rows <- character()
refusals <- character()
record <- function(case, expr) {
  result <- tryCatch(suppressWarnings(suppressMessages(expr)), error = function(e) e)
  if (inherits(result, "error")) {
    refusals <<- c(refusals, paste(case, conditionMessage(result), sep = "\t"))
    return(invisible())
  }
  d <- as.matrix(result)
  for (i in seq_len(nrow(d))) {
    for (j in seq_len(ncol(d))) {
      rows <<- c(rows, paste(case, i - 1, j - 1, sprintf("%.17g", d[i, j]), sep = ","))
    }
  }
}

record("chi2_obs_plain", seqdist(plain_analogue, method = "CHI2", global.pdotj = "obs"))
record("chi2_obs_plain_auto", seqdist(plain_analogue, method = "CHI2", global.pdotj = "obs", norm = "auto"))
record("euclid_obs_plain", seqdist(plain_analogue, method = "EUCLID", global.pdotj = "obs"))
record("chi2_obs_unobserved_last", seqdist(unobserved_last, method = "CHI2", global.pdotj = "obs"))
record("chi2_vec4_unobserved_last", seqdist(unobserved_last, method = "CHI2", global.pdotj = c(0.2, 0.3, 0.4, 0.1)))
record("chi2_vec3_unobserved_last", seqdist(unobserved_last, method = "CHI2", global.pdotj = c(0.2, 0.3, 0.5)))
record("euclid_obs_unobserved_last", seqdist(unobserved_last, method = "EUCLID", global.pdotj = "obs"))
record("chi2_obs_ragged", seqdist(ragged, method = "CHI2", global.pdotj = "obs"))
record("euclid_obs_ragged", seqdist(ragged, method = "EUCLID", global.pdotj = "obs"))
record("chi2_obs_with_missing", seqdist(with_missing, method = "CHI2", global.pdotj = "obs", with.missing = TRUE))
record("nms_prox4_unobserved_last", seqdist(letters_dense, method = "NMS", prox = P4))
record("nms_prox3_unobserved_last", seqdist(letters_dense, method = "NMS", prox = P4[1:3, 1:3]))
record("svrspell_prox4_unobserved_last", seqdist(letters_dense, method = "SVRspell", prox = P4))
record("svrspell_prox3_unobserved_last", seqdist(letters_dense, method = "SVRspell", prox = P4[1:3, 1:3]))

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
writeLines(c("case,row,col,value", rows), file.path(out_dir, "seqdist_state_space_cases.csv"))
writeLines(c("# case\tR error message", refusals), file.path(out_dir, "seqdist_state_space_refusals.txt"))
cat("SEQDIST_STATE_SPACE_GOLDENS_DONE\n")
