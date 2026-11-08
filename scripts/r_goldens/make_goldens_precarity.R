# Golden generation: Phase 2 R-faithful precarity/complexity indicators.
# Anchors traminex's seqintegr/seqprecstart/seqibad/seqivardur/seqmaintokens/
# seqipos/seqlogp/seqprecarity to TraMineR 2.2-14. Reuses Phase 1 fixtures A/B/C
# (see make_goldens_seqdss_seqdur.R) plus a small seqlogp fixture L.
# Run from the repository root with the oracle library (scripts/r_goldens/README.md):
#   Rscript scripts/r_goldens/make_goldens_precarity.R
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
asnum <- TraMineR:::seqasnum
out_dir <- file.path("tests", "r_goldens")
if (!identical(as.character(packageVersion("TraMineR")), "2.2.14")) {
  stop("precarity goldens require TraMineR 2.2-14")
}

write_num_matrix <- function(name, m, header) {
  con <- file(file.path(out_dir, paste0(name, ".csv")), open = "wt")
  cat(header, file = con)
  txt <- apply(as.matrix(m), 1, function(r) paste(formatC(r, format = "g", digits = 15), collapse = ","))
  writeLines(txt, con)
  close(con)
}
write_int_matrix <- function(name, m, header) {
  con <- file(file.path(out_dir, paste0(name, ".csv")), open = "wt")
  cat(header, file = con)
  utils::write.table(as.matrix(m), con, sep = ",", row.names = FALSE, col.names = FALSE)
  close(con)
}

## ---- Fixtures (A, B void-free; C ragged; L small for seqlogp) ----
A <- rbind(
  c("A", "A", "B", "C", "C", "C"),
  c("A", "B", "B", "B", "C", "A"),
  c("D", "D", "D", "D", "D", "D"),
  c("A", "B", "A", "B", "A", "B")
)
sdA <- seqdef(A, alphabet = c("A", "B", "C", "D"))
B <- rbind(
  c("A", "A", "*", "B", "B", "C"),
  c("*", "A", "A", "A", "*", "*"),
  c("A", "*", "A", "B", "B", "B"),
  c("C", "C", "C", "C", "C", "C")
)
sdB <- seqdef(B, alphabet = c("A", "B", "C"), missing = "*", right = NA)
C <- rbind(
  c("A", "A", "B", "C", "C", "C"),
  c("A", "B", "B", "*", "*", "*"),
  c("A", "*", "A", "B", "*", "*"),
  c("C", "C", "*", "*", "*", "*")
)
sdC <- seqdef(C, alphabet = c("A", "B", "C"), missing = "*")
L <- rbind(c("A", "B"), c("A", "A"))
sdL <- seqdef(L, alphabet = c("A", "B"))
# Fixture M: most-frequent states are OUT of alphabet order (C >> A > B), so a
# faithful seqmaintokens (sorted-ascending indices) diverges from the proxy
# (descending-by-frequency). meant: C=3.33, A=2.0, B=0.67 -> top-2 = {A,C}.
M <- rbind(
  c("C", "C", "C", "A", "A", "B"),
  c("C", "C", "C", "A", "A", "A"),
  c("C", "C", "C", "C", "A", "B")
)
sdM <- seqdef(M, alphabet = c("A", "B", "C"))

# Fixture Kr (F1): variable-length sequences, void padding (right="DEL"), NO interior
# missing -> seqlogp default path (with.missing=FALSE). Pins void-faithful surprisal.
Kr <- rbind(
  c("A", "B", "A", "B"),
  c("A", "A", "B", NA),
  c("B", "A", NA, NA),
  c("A", "B", "B", "A")
)
sdKr <- seqdef(Kr, alphabet = c("A", "B"))
# Fixture J (F2): meant out of order so the seqmaintokens mint+(#>=mint > k) branch
# (R returns thresholded-subset positions) is exercised. meant: A=.333,B=1.667,C=1.333,D=2.667.
J <- rbind(
  c("A", "B", "B", "D", "D", "D"),
  c("B", "B", "D", "D", "C", "C"),
  c("B", "D", "D", "D", "C", "C")
)
sdJ <- seqdef(J, alphabet = c("A", "B", "C", "D"))
# Fixture V (F4): an all-void row (row 3) and a ragged row, plus full rows. R divides
# by zero on the empty row -> NaN; pins the NaN (not 0) behaviour.
V <- rbind(
  c("A", "A", "B", "B"),
  c("B", "A", "B", "A"),
  c(NA, NA, NA, NA),
  c("A", "B", NA, NA)
)
sdV <- seqdef(V, alphabet = c("A", "B"))
# Fixture W (F5): void-padded (ragged) rows incl. a single-state+void row, no nr.
# Pins void-faithful seqici / seqient (and thus seqprecarity on ragged data).
W <- rbind(
  c("A", "A", "B", NA, NA),
  c("A", "B", "A", "B", "B"),
  c("A", "A", "A", NA, NA)
)
sdW <- seqdef(W, alphabet = c("A", "B"))

## ---- seqintegr (cols = alphabet(with.missing)) ----
write_num_matrix("seqintegr_A",    seqintegr(sdA, pow = 1, with.missing = FALSE), "# R: seqintegr(sdA, pow=1, wm=F); cols A,B,C,D\n")
write_num_matrix("seqintegr_A_p2", seqintegr(sdA, pow = 2, with.missing = FALSE), "# R: seqintegr(sdA, pow=2, wm=F)\n")
write_num_matrix("seqintegr_B_nm", seqintegr(sdB, pow = 1, with.missing = FALSE), "# R: seqintegr(sdB, pow=1, wm=F); cols A,B,C\n")
write_num_matrix("seqintegr_B_m",  seqintegr(sdB, pow = 1, with.missing = TRUE),  "# R: seqintegr(sdB, pow=1, wm=T); cols A,B,C,nr\n")

## ---- seqprecstart (length-K_eff vector; write as a 1-row matrix) ----
write_num_matrix("seqprecstart_A",   t(as.vector(seqprecstart(sdA))),                      "# R: seqprecstart(sdA); default; states A,B,C,D\n")
write_num_matrix("seqprecstart_B_m", t(as.vector(seqprecstart(sdB, with.missing = TRUE))), "# R: seqprecstart(sdB, wm=T); states A,B,C,nr\n")

## ---- seqibad ----
write_num_matrix("seqibad_A",    as.vector(seqibad(sdA, pow = 1, with.missing = FALSE)), "# R: seqibad(sdA, pow=1, wm=F)\n")
write_num_matrix("seqibad_A_p2", as.vector(seqibad(sdA, pow = 2, with.missing = FALSE)), "# R: seqibad(sdA, pow=2, wm=F)\n")
write_num_matrix("seqibad_B_m",  as.vector(seqibad(sdB, pow = 1, with.missing = TRUE)),  "# R: seqibad(sdB, pow=1, wm=T)\n")

## ---- seqivardur (per-sequence variance) ----
write_num_matrix("seqivardur_A_t1",  as.vector(seqivardur(sdA, type = 1, with.missing = FALSE)), "# R: seqivardur(sdA, type=1, wm=F)\n")
write_num_matrix("seqivardur_A_t2",  as.vector(seqivardur(sdA, type = 2, with.missing = FALSE)), "# R: seqivardur(sdA, type=2, wm=F)\n")
write_num_matrix("seqivardur_B_t1m", as.vector(seqivardur(sdB, type = 1, with.missing = TRUE)),  "# R: seqivardur(sdB, type=1, wm=T)\n")
write_num_matrix("seqivardur_C_t1m", as.vector(seqivardur(sdC, type = 1, with.missing = TRUE)),  "# R: seqivardur(sdC, type=1, wm=T)\n")

## ---- seqmaintokens (1-based R indices; subtract 1 in Python) ----
write_int_matrix("seqmaintokens_A_k2",   t(as.integer(seqmaintokens(sdA, k = 2))),      "# R: seqmaintokens(sdA, k=2); 1-based\n")
write_int_matrix("seqmaintokens_A_mint", t(as.integer(seqmaintokens(sdA, mint = 1.5))), "# R: seqmaintokens(sdA, mint=1.5); 1-based\n")
write_int_matrix("seqmaintokens_M_k2",   t(as.integer(seqmaintokens(sdM, k = 2))),      "# R: seqmaintokens(sdM, k=2); 1-based; top states out of alpha order\n")

## ---- seqipos (index=share default; also integr, volatility) ----
write_num_matrix("seqipos_A_share",  as.vector(seqipos(sdA, pos.states = c("A", "B"), neg.states = c("C", "D"), index = "share")),      "# R: seqipos(sdA, pos=A,B, neg=C,D, share)\n")
write_num_matrix("seqipos_A_integr", as.vector(seqipos(sdA, pos.states = c("A", "B"), neg.states = c("C", "D"), index = "integr")),     "# R: seqipos(sdA, integr)\n")
write_num_matrix("seqipos_A_vol",    as.vector(seqipos(sdA, pos.states = c("A", "B"), neg.states = c("C", "D"), index = "volatility")), "# R: seqipos(sdA, volatility)\n")

## ---- seqlogp (Markov surprisal) ----
write_num_matrix("seqlogp_L",   as.vector(seqlogp(sdL)),                        "# R: seqlogp(sdL)\n")
write_num_matrix("seqlogp_A",   as.vector(seqlogp(sdA)),                        "# R: seqlogp(sdA)\n")
write_num_matrix("seqlogp_B_m", as.vector(seqlogp(sdB, with.missing = TRUE)),   "# R: seqlogp(sdB, wm=T)\n")

## ---- seqprecarity (type=1 default path; plus method=RANK) ----
write_num_matrix("seqprecarity_A",      as.vector(seqprecarity(sdA)),                      "# R: seqprecarity(sdA); TRATEDSS default\n")
write_num_matrix("seqprecarity_A_rank", as.vector(seqprecarity(sdA, method = "RANK")),     "# R: seqprecarity(sdA, method=RANK)\n")
write_num_matrix("seqprecarity_B_m",    as.vector(seqprecarity(sdB, with.missing = TRUE)), "# R: seqprecarity(sdB, wm=T)\n")

## ---- F1: seqlogp on ragged (void-padded) data (default with.missing=FALSE) ----
write_num_matrix("seqlogp_Kr", as.vector(seqlogp(sdKr)), "# R: seqlogp(sdKr); ragged void padding, default wm=F\n")

## ---- F2: seqmaintokens mint + (#>=mint > k) branch (R subset positions) ----
write_int_matrix("seqmaintokens_J_mint_k2", t(as.integer(seqmaintokens(sdJ, mint = 1.0, k = 2))), "# R: seqmaintokens(sdJ, mint=1.0, k=2); 1-based; mint+(#>k) subset-position branch\n")

## ---- F6 (cycle-6 F-001): the DEFAULT k=8 over statl, on missing-bearing data ----
# R indexes `statl` = alphabet(sd), or c(alphabet(sd), nr) when with.missing=TRUE --
# never the internal codebook, so the nr/void markers are not selectable states.
# With k=8 > length(statl), R's `[1:k]` yields NA beyond the end and sort() drops
# them, so the result is the whole of statl. sdB carries nr; sdC carries nr AND
# void; sdA carries neither, so `with.missing=TRUE` there pins the phantom-nr row
# (mean 0) that R still counts as a selectable token.
write_int_matrix("seqmaintokens_B_default",   t(as.integer(seqmaintokens(sdB))),                       "# R: seqmaintokens(sdB); DEFAULT k=8; 1-based into alphabet(sdB)=A,B,C\n")
write_int_matrix("seqmaintokens_B_default_m", t(as.integer(seqmaintokens(sdB, with.missing = TRUE))),  "# R: seqmaintokens(sdB, wm=T); DEFAULT k=8; 1-based into c(A,B,C,nr)\n")
write_int_matrix("seqmaintokens_C_default",   t(as.integer(seqmaintokens(sdC))),                       "# R: seqmaintokens(sdC); DEFAULT k=8; ragged+nr; 1-based into alphabet(sdC)=A,B,C\n")
write_int_matrix("seqmaintokens_C_default_m", t(as.integer(seqmaintokens(sdC, with.missing = TRUE))),  "# R: seqmaintokens(sdC, wm=T); DEFAULT k=8; ragged+nr; 1-based into c(A,B,C,nr)\n")
write_int_matrix("seqmaintokens_A_default_m", t(as.integer(seqmaintokens(sdA, with.missing = TRUE))),  "# R: seqmaintokens(sdA, wm=T); DEFAULT k=8; phantom nr (mean 0) still selected\n")

## ---- F3: include_missing=TRUE on missing-free data (phantom nr; sdA has no missing) ----
write_num_matrix("seqintegr_A_m",    seqintegr(sdA, pow = 1, with.missing = TRUE),        "# R: seqintegr(sdA, wm=T); phantom nr col; cols A,B,C,D,nr\n")
write_num_matrix("seqprecstart_A_m", t(as.vector(seqprecstart(sdA, with.missing = TRUE))), "# R: seqprecstart(sdA, wm=T); phantom nr; states A,B,C,D,nr\n")
write_num_matrix("seqibad_A_m",      as.vector(seqibad(sdA, pow = 1, with.missing = TRUE)), "# R: seqibad(sdA, wm=T); phantom nr\n")
write_num_matrix("seqprecarity_A_m", as.vector(seqprecarity(sdA, with.missing = TRUE)),     "# R: seqprecarity(sdA, wm=T); phantom nr\n")

## ---- F4: all-void / empty rows -> R NaN ----
write_num_matrix("seqintegr_V",   seqintegr(sdV, pow = 1, with.missing = FALSE), "# R: seqintegr(sdV, wm=F); row3 all-void -> NaN\n")
write_num_matrix("seqibad_V",     as.vector(seqibad(sdV, pow = 1)),              "# R: seqibad(sdV); row3 NaN\n")
write_num_matrix("seqivardur_V",  as.vector(seqivardur(sdV, type = 1)),          "# R: seqivardur(sdV, type=1); row3 NaN\n")
write_num_matrix("seqprecarity_V", as.vector(seqprecarity(sdV)),                 "# R: seqprecarity(sdV); row3 NaN\n")

## ---- F5: void-faithful seqici / seqient (and seqprecarity) on ragged data ----
write_num_matrix("seqient_W", as.vector(seqient(sdW, norm = TRUE)), "# R: seqient(sdW, norm=T); void-padded\n")
write_num_matrix("seqici_W",  as.vector(seqici(sdW)),               "# R: seqici(sdW); void-padded\n")
write_num_matrix("seqprecarity_W", as.vector(seqprecarity(sdW)),    "# R: seqprecarity(sdW); void-padded\n")

## ---- R2: seqient logarithm bases, natural-log norm, and messages ----
# sdD declares nr but contains no missing observations (phantom nr route).
sdD <- seqdef(A[, 1:4], alphabet = c("A", "B", "C", "D"), missing = "*")
ient_rows <- list()
append_ient <- function(case, sd, base, norm, wm) {
  value <- suppressMessages(seqient(
    sd,
    base = base,
    norm = norm,
    with.missing = wm,
    silent = TRUE
  ))
  ient_rows[[length(ient_rows) + 1L]] <<- data.frame(
    case = case,
    row = seq_along(value),
    base = base,
    norm = norm,
    with_missing = wm,
    value = as.vector(value),
    stringsAsFactors = FALSE
  )
}

for (base in c(exp(1), 2, 10)) {
  for (norm in c(FALSE, TRUE)) {
    append_ient(
      paste0(
        "dense_base",
        trimws(formatC(base, format = "g", digits = 15)),
        "_norm",
        norm
      ),
      sdA,
      base,
      norm,
      FALSE
    )
  }
}
append_ient("declared_unobserved_missing_wm0", sdD, 2, TRUE, FALSE)
append_ient("declared_unobserved_missing_wm1", sdD, 2, TRUE, TRUE)
ient <- do.call(rbind, ient_rows)
utils::write.csv(
  ient,
  file.path(out_dir, "seqient_r2_base.csv"),
  row.names = FALSE,
  na = "NaN",
  quote = FALSE
)

ient_message <- capture.output(
  invisible(seqient(sdA, silent = FALSE)),
  type = "message"
)
ient_silent <- capture.output(
  invisible(seqient(sdA, silent = TRUE)),
  type = "message"
)
invalid_disposition <- function(base) {
  warnings <- character()
  value <- tryCatch(
    withCallingHandlers(
      seqient(sdL, base = base, silent = TRUE),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  outcome <- if (inherits(value, "error")) {
    paste0("ERROR:", conditionMessage(value))
  } else {
    paste(format(value, digits = 17), collapse = ",")
  }
  paste0(outcome, ";warnings=", paste(warnings, collapse = " | "))
}
invalid_bases <- list(zero = 0, negative = -1, one = 1, infinite = Inf, nan = NaN)
writeLines(
  c(
    paste0("progress=", ient_message),
    paste0("silent=", paste(ient_silent, collapse = "")),
    vapply(
      names(invalid_bases),
      function(name) paste0(
        "invalid_", name, "=", invalid_disposition(invalid_bases[[name]])
      ),
      character(1)
    )
  ),
  file.path(out_dir, "seqient_r2_messages.txt")
)

## ---- R2: seqlogp probability and starting-distribution routes ----
r2_weights <- c(1.5, 2.0, 0.5, 3.0)
sdAw <- seqdef(A, alphabet = c("A", "B", "C", "D"), weights = r2_weights)
sdBw <- seqdef(
  B,
  alphabet = c("A", "B", "C"),
  missing = "*",
  right = NA,
  weights = r2_weights
)

custom_stationary <- matrix(
  c(
    .55, .25, .15, .05,
    .10, .60, .20, .10,
    .15, .15, .60, .10,
    .20, .10, .10, .60
  ),
  nrow = 4,
  byrow = TRUE,
  dimnames = list(
    paste0("from_", rev(c("A", "B", "C", "D"))),
    paste0("to_", rev(c("A", "B", "C", "D")))
  )
)
custom_timevarying <- array(0, dim = c(4, 4, 5))
for (time in seq_len(5)) {
  for (from in seq_len(4)) {
    values <- ((from + 2 * seq_len(4) + time) %% 7) + 1
    custom_timevarying[from, , time] <- values / sum(values)
  }
}
dimnames(custom_timevarying) <- list(
  paste0("from_", rev(c("A", "B", "C", "D"))),
  paste0("to_", rev(c("A", "B", "C", "D"))),
  paste0("age_", 5:1)
)
custom_begin <- c(.1, .2, .3, .4)

logp_rows <- list()
append_logp <- function(
  case,
  sd,
  prob_tag,
  prob_value,
  time_varying,
  begin_tag,
  begin_value,
  weighted,
  wm
) {
  value <- suppressMessages(seqlogp(
    sd,
    prob = prob_value,
    time.varying = time_varying,
    begin = begin_value,
    weighted = weighted,
    with.missing = wm
  ))
  state_order <- if (wm) {
    c(alphabet(sd), attr(sd, "nr"))
  } else {
    alphabet(sd)
  }
  nstate <- max(TraMineR:::seqasnum(sd, with.missing = wm) + 1, na.rm = TRUE)
  ntransition <- max(seqlength(sd, with.missing = wm)) - 1
  logp_rows[[length(logp_rows) + 1L]] <<- data.frame(
    case = case,
    row = seq_along(value),
    prob = prob_tag,
    time_varying = time_varying,
    begin = begin_tag,
    weighted = weighted,
    include_missing = wm,
    nstate = nstate,
    ntransition = ntransition,
    state_order = paste(state_order, collapse = "|"),
    value = as.vector(value),
    stringsAsFactors = FALSE
  )
}

model_specs <- list(
  list(tag = "trate", value = "trate", tv = TRUE),
  list(tag = "trate", value = "trate", tv = FALSE),
  list(tag = "freq", value = "freq", tv = TRUE),
  list(tag = "freq", value = "freq", tv = FALSE),
  list(tag = "matrix", value = custom_stationary, tv = FALSE),
  list(tag = "array", value = custom_timevarying, tv = TRUE)
)
begin_specs <- list(
  list(tag = "freq", value = "freq"),
  list(tag = "global.freq", value = "global.freq"),
  list(tag = "array", value = custom_begin)
)

for (wm in c(FALSE, TRUE)) {
  sd <- if (wm) sdBw else sdAw
  for (weighted in c(FALSE, TRUE)) {
    for (model in model_specs) {
      for (starting in begin_specs) {
        case <- paste0(
          "model_", model$tag,
          "_tv", model$tv,
          "_begin_", starting$tag,
          "_weighted", weighted,
          "_wm", wm
        )
        append_logp(
          case,
          sd,
          model$tag,
          model$value,
          model$tv,
          starting$tag,
          starting$value,
          weighted,
          wm
        )
      }
    }
  }
}

# R bases custom-model width on max observed numeric code, not on a phantom
# declared missing state.  This valid route freezes that source-level quirk.
append_logp(
  "phantom_missing_matrix4",
  sdAw,
  "matrix",
  custom_stationary,
  FALSE,
  "array",
  custom_begin,
  TRUE,
  TRUE
)

# No smoothing: non-positive and NaN probabilities are skipped, while Inf
# contributes -Inf.  Use one row per starting state to expose each route.
edge <- rbind(
  c("A", "B", "C"),
  c("B", "C", "D"),
  c("C", "D", "A"),
  c("D", "A", "B")
)
sdEdge <- seqdef(edge, alphabet = c("A", "B", "C", "D"))
edge_prob <- matrix(.25, nrow = 4, ncol = 4)
edge_prob[1, 2] <- 0
edge_prob[2, 3] <- NaN
edge_prob[3, 4] <- Inf
edge_prob[4, 1] <- .5
append_logp(
  "edge_zero_nan_inf",
  sdEdge,
  "matrix_nonfinite",
  edge_prob,
  FALSE,
  "array_nonfinite",
  c(0, NaN, Inf, .5),
  FALSE,
  FALSE
)

utils::write.csv(
  do.call(rbind, logp_rows),
  file.path(out_dir, "seqlogp_r2_cases.csv"),
  row.names = FALSE,
  na = "NaN",
  quote = FALSE
)

normalize_logp <- function(expr) {
  tryCatch(
    {
      value <- suppressMessages(force(expr))
      paste0("VALUE:", paste(format(value, digits = 17), collapse = ","))
    },
    error = function(e) paste0("ERROR:", conditionMessage(e))
  )
}
writeLines(
  c(
    "nomissing_model_shape=4x4x5;states=A|B|C|D",
    "missing_model_shape=4x4x5;states=A|B|C|*",
    "phantom_missing_alphabet_width=5;observed_model_width=4",
    paste0(
      "invalid_prob_global.trate=",
      normalize_logp(seqlogp(sdAw, prob = "global.trate"))
    ),
    paste0(
      "invalid_prob_global.freq=",
      normalize_logp(seqlogp(sdAw, prob = "global.freq"))
    ),
    paste0(
      "malformed_prob_vector=",
      normalize_logp(seqlogp(sdAw, prob = rep(.25, 4)))
    ),
    paste0(
      "malformed_prob_matrix=",
      normalize_logp(seqlogp(sdAw, prob = diag(3)))
    ),
    paste0(
      "malformed_prob_array=",
      normalize_logp(seqlogp(sdAw, prob = array(0, dim = c(4, 4, 4))))
    ),
    paste0(
      "malformed_prob_rank4=",
      normalize_logp(seqlogp(sdAw, prob = array(0, dim = c(4, 4, 5, 1))))
    ),
    paste0(
      "malformed_begin_length=",
      normalize_logp(seqlogp(sdAw, begin = c(.2, .3, .5)))
    ),
    paste0(
      "malformed_begin_string=",
      normalize_logp(seqlogp(sdAw, begin = "trate"))
    ),
    paste0(
      "phantom_missing_begin5=",
      normalize_logp(seqlogp(
        sdAw,
        prob = diag(5),
        begin = rep(.2, 5),
        with.missing = TRUE
      ))
    )
  ),
  file.path(out_dir, "seqlogp_r2_errors.txt")
)

cat("PRECARITY_GOLDENS_DONE\n")
