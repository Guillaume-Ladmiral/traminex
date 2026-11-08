# Live-R goldens for dissrep / seqrep (RP4 T13).
#
# dissrep: representative-set selection + quality Statistics/Quality, on the
# shared phase2 15x15 distance matrix (deterministic, no ties in D). seqrep
# goldens (Task 5) are appended to this file when that task lands.
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

gen <- function(tag, crit, weights) {
  r <- suppressMessages(
    dissrep(D, criterion = crit, coverage = 0.5, weights = weights)
  )
  # representative object indices, R is 1-based -> dump 0-based
  write_phase2_csv(
    paste0("dissrep_", tag, "_idx"), matrix(as.vector(r) - 1L, nrow = 1),
    paste0("# dissrep(criterion=", crit, ", coverage=0.5) representative indices (0-based)\n")
  )
  st <- attr(r, "Statistics")
  write_phase2_csv(
    paste0("dissrep_", tag, "_stats"), as.matrix(st),
    "# dissrep Statistics: cols na,na(%),nb,nb(%),SD,MD,DC,V,Q ; rows r1..rk,Total\n"
  )
  write_phase2_csv(
    paste0("dissrep_", tag, "_quality"), matrix(attr(r, "Quality"), nrow = 1),
    "# dissrep overall Quality\n"
  )
}

gen("density_uw", "density", NULL)
gen("density_w", "density", w)
gen("dist_uw", "dist", NULL)
gen("dist_w", "dist", w)

cat("DISSREP_GOLDENS_DONE\n")

## seqrep goldens (RP4 T13c): seqrep wraps dissrep on a sequence diss matrix ----
set.seed(99001)
nr <- 10
Lr <- 6
Mr <- matrix(sample(1:4, nr * Lr, replace = TRUE), nrow = nr, ncol = Lr)
wr <- 0.5 + 0.1 * (1:nr)
sdr <- suppressMessages(seqdef(Mr, alphabet = 1:4, weights = wr))
Dr <- suppressMessages(seqdist(sdr, method = "LCS", norm = "none"))

write_phase2_csv("seqrep_input_codes", Mr - 1L, "# seqrep input: 0-based codes (alphabet 1:4)\n")
write_phase2_csv("seqrep_input_diss", Dr, "# seqrep input: LCS norm=none distance matrix\n")
write_phase2_csv("seqrep_input_weights", matrix(wr, nrow = 1), "# seqrep input weights 0.5+0.1*i\n")

genrep <- function(tag, crit, weighted) {
  r <- suppressMessages(
    seqrep(sdr, diss = Dr, criterion = crit, coverage = 0.25, weighted = weighted)
  )
  write_phase2_csv(
    paste0("seqrep_", tag, "_idx"), matrix(attr(r, "Index") - 1L, nrow = 1),
    paste0("# seqrep(criterion=", crit, ") Index (0-based)\n")
  )
  write_phase2_csv(
    paste0("seqrep_", tag, "_stats"), as.matrix(attr(r, "Statistics")),
    "# seqrep Statistics: na,na(%),nb,nb(%),SD,MD,DC,V,Q\n"
  )
  write_phase2_csv(
    paste0("seqrep_", tag, "_quality"), matrix(attr(r, "Quality"), nrow = 1),
    "# seqrep Quality\n"
  )
}

genrep("density_uw", "density", FALSE)
genrep("density_w", "density", TRUE)
genrep("dist_w", "dist", TRUE)

cat("SEQREP_GOLDENS_DONE\n")

## R2 seqrep criteria/options -------------------------------------------------
r2_rows <- list()
append_r2 <- function(case, section, row, col, value = NA_real_, text = "") {
  r2_rows[[length(r2_rows) + 1L]] <<- data.frame(
    case = case,
    section = section,
    row = row,
    col = col,
    value = value,
    text = text,
    stringsAsFactors = FALSE
  )
}

r2_messages <- character()
run_r2 <- function(
  case,
  criterion,
  weighted,
  score = NULL,
  decreasing = TRUE,
  coverage = .25,
  nrep = NULL,
  pradius = .10,
  dmax = NULL
) {
  messages <- character()
  warnings <- character()
  r <- withCallingHandlers(
    seqrep(
      sdr,
      diss = Dr,
      criterion = criterion,
      score = score,
      decreasing = decreasing,
      coverage = coverage,
      nrep = nrep,
      pradius = pradius,
      dmax = dmax,
      weighted = weighted
    ),
    message = function(m) {
      messages <<- c(
        messages,
        trimws(gsub("[\r\n]+", " ", conditionMessage(m)))
      )
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      warnings <<- c(
        warnings,
        trimws(gsub("[\r\n]+", " ", conditionMessage(w)))
      )
      invokeRestart("muffleWarning")
    }
  )

  index <- attr(r, "Index")
  for (position in seq_along(index)) {
    append_r2(case, "rep_index", position, 1, index[position] - 1L)
    for (time in seq_len(ncol(r))) {
      append_r2(
        case,
        "representative",
        position,
        time,
        TraMineR:::seqasnum(r[position, , drop = FALSE])[1, time]
      )
    }
  }
  scores <- attr(r, "Scores")
  for (i in seq_along(scores)) append_r2(case, "score", i, 1, scores[i])
  statistics <- as.matrix(attr(r, "Statistics"))
  for (i in seq_len(nrow(statistics))) {
    for (j in seq_len(ncol(statistics))) {
      append_r2(case, "statistic", i, j, statistics[i, j], colnames(statistics)[j])
    }
  }
  append_r2(case, "quality", 1, 1, attr(r, "Quality"))
  rep_group <- attr(r, "Rep.group")
  for (i in seq_along(rep_group)) {
    append_r2(case, "rep_group", i, 1, rep_group[i] - 1L)
  }
  distances <- as.matrix(attr(r, "Distances"))
  for (i in seq_len(nrow(distances))) {
    for (j in seq_len(ncol(distances))) {
      append_r2(case, "distance", i, j, distances[i, j])
    }
  }
  idx_rep <- attr(r, "idx.rep")
  for (i in seq_along(idx_rep)) {
    for (j in seq_along(idx_rep[[i]])) {
      append_r2(case, "idx_rep", i, j, idx_rep[[i]][j] - 1L)
    }
  }
  append_r2(case, "attribute", 1, 1, attr(r, "nbseq"), "nbseq")
  append_r2(case, "attribute", 2, 1, attr(r, "dmax"), "dmax")
  append_r2(case, "attribute", 3, 1, NA_real_, paste0("criterion=", criterion))
  append_r2(
    case,
    "attribute",
    4,
    1,
    NA_real_,
    paste0("weighted=", as.character(attr(r, "weighted")))
  )
  r2_messages <<- c(
    r2_messages,
    paste0(case, ":messages=", paste(messages, collapse = " | ")),
    paste0(case, ":warnings=", paste(warnings, collapse = " | "))
  )
}

explicit_ties <- c(4, 4, 3, 3, 2, 2, 1, 1, 0, 0)
run_r2("density_uw", "density", FALSE)
run_r2("density_w", "density", TRUE)
run_r2("dist_uw", "dist", FALSE)
run_r2("dist_w", "dist", TRUE)
run_r2("prob_uw", "prob", FALSE)
run_r2("prob_w", "prob", TRUE)
run_r2("mscore_uw", "mscore", FALSE)
run_r2("mscore_w", "mscore", TRUE)
run_r2(
  "explicit_mscore_desc_ties",
  "mscore",
  TRUE,
  score = explicit_ties,
  decreasing = TRUE,
  coverage = .4
)
run_r2(
  "explicit_prob_asc_ties_nrep2",
  "prob",
  FALSE,
  score = explicit_ties,
  decreasing = FALSE,
  nrep = 2L
)
run_r2(
  "density_nrep2_radius0_dmax12",
  "density",
  TRUE,
  nrep = 2L,
  pradius = 0,
  dmax = 12
)
run_r2(
  "dist_nrep3_radius20_dmax12",
  "dist",
  FALSE,
  nrep = 3L,
  pradius = .20,
  dmax = 12
)

utils::write.csv(
  do.call(rbind, r2_rows),
  file.path("tests", "r_goldens", "seqrep_r2_cases.csv"),
  row.names = FALSE,
  na = "NaN",
  quote = FALSE
)

normalize_seqrep_error <- function(expr) {
  tryCatch(
    {
      suppressMessages(force(expr))
      "VALUE"
    },
    error = function(e) paste0("ERROR:", conditionMessage(e))
  )
}
writeLines(
  c(
    r2_messages,
    paste0(
      "invalid_criterion=",
      normalize_seqrep_error(seqrep(sdr, diss = Dr, criterion = "bogus"))
    ),
    paste0(
      "invalid_score_length=",
      normalize_seqrep_error(seqrep(
        sdr,
        diss = Dr,
        criterion = "prob",
        score = c(1, 2)
      ))
    ),
    paste0(
      "invalid_pradius=",
      normalize_seqrep_error(seqrep(sdr, diss = Dr, pradius = 1.1))
    ),
    paste0(
      "invalid_coverage_zero_nrep_null=",
      normalize_seqrep_error(seqrep(sdr, diss = Dr, coverage = 0, nrep = NULL))
    )
  ),
  file.path("tests", "r_goldens", "seqrep_r2_errors.txt")
)

cat("SEQREP_R2_GOLDENS_DONE\n")
