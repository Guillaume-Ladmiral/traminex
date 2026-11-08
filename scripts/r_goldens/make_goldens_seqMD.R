# Golden generation: seqMD / seqdistmc (D2 -- CAT combined-OM multidomain).
#
# Anchors the faithful CAT engine to TraMineR 2.2-14 seqdistmc / seqMD(what="diss").
# CAT = "additive trick": build MD sequences of combined domain states on the
# expanded alphabet, derive MD substitution/indel costs by summing the per-domain
# (channel-weighted) costs, then run one OM on the combined sequences.
#
# Requires live R + TraMineR 2.2-14; see memory "live-traminer-r-env".
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

write_csv <- function(name, M, header) {
  out_dir <- file.path("tests", "r_goldens")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  M <- as.matrix(M)
  dimnames(M) <- NULL
  path <- file.path(out_dir, paste0(name, ".csv"))
  con <- file(path, "w")
  cat(header, file = con)
  for (r in seq_len(nrow(M))) {
    cat(paste(sprintf("%.15g", M[r, ]), collapse = ","), "\n", file = con, sep = "")
  }
  close(con)
}

ch1 <- suppressMessages(seqdef(rbind(
  c("A", "A", "B", "B", "C", "C"), c("A", "B", "B", "C", "C", "A"),
  c("C", "C", "B", "B", "A", "A"), c("B", "B", "B", "C", "C", "C")),
  alphabet = c("A", "B", "C")))
ch2 <- suppressMessages(seqdef(rbind(
  c("X", "X", "Y", "Y", "Y", "Z"), c("Y", "Y", "Z", "Z", "X", "X"),
  c("Z", "Z", "Z", "X", "X", "Y"), c("X", "Y", "Y", "Y", "Z", "Z")),
  alphabet = c("X", "Y", "Z")))

# what="diss", default link="sum", equal channel weights, CONSTANT sm, indel 1
write_csv(
  "seqMD_OM_sum",
  suppressMessages(seqdistmc(list(ch1, ch2), method = "OM",
    sm = "CONSTANT", indel = 1, cweight = c(1, 1))),
  "# R: TraMineR::seqdistmc(list(ch1,ch2), method='OM', sm='CONSTANT', indel=1, cweight=c(1,1)); see make_goldens_seqMD.R\n"
)
# cweighted channels c(2,1)
write_csv(
  "seqMD_OM_cw",
  suppressMessages(seqdistmc(list(ch1, ch2), method = "OM",
    sm = "CONSTANT", indel = 1, cweight = c(2, 1))),
  "# R: TraMineR::seqdistmc(list(ch1,ch2), method='OM', sm='CONSTANT', indel=1, cweight=c(2,1)); see make_goldens_seqMD.R\n"
)
# LCS family (handled by R as OM with CONSTANT sm, indel 1, cval 2)
write_csv(
  "seqMD_LCS",
  suppressMessages(seqdistmc(list(ch1, ch2), method = "LCS", cweight = c(1, 1))),
  "# R: TraMineR::seqdistmc(list(ch1,ch2), method='LCS', cweight=c(1,1)); see make_goldens_seqMD.R\n"
)
# seqMD(what="diss") must equal seqdistmc (alias identity)
write_csv(
  "seqMD_alias",
  suppressMessages(seqMD(list(ch1, ch2), method = "OM",
    sm = "CONSTANT", indel = 1, what = "diss", cweight = c(1, 1))),
  "# R: TraMineR::seqMD(list(ch1,ch2), method='OM', sm='CONSTANT', indel=1, what='diss', cweight=c(1,1)); see make_goldens_seqMD.R\n"
)
cat("SEQMD_GOLDENS_DONE\n")

## R2: declared channel states, ragged channels, and DHD CAT arrays -----------
##
## The long manifest deliberately records both alphabets involved in CAT:
## channel costs use every declared channel state, while the combined alphabet
## contains only joined states observed in the MD sequence.  Array storage is
## the public R order [combined state, combined state, time].
r2_rows <- list()
r2_add <- function(case, section, row = 0L, col = 0L, time = 0L,
                   value = NaN, text = "") {
  r2_rows[[length(r2_rows) + 1L]] <<- c(
    case = as.character(case), section = as.character(section),
    row = as.character(row), col = as.character(col), time = as.character(time),
    value = if (is.nan(value)) "NaN" else if (is.na(value)) "NA" else
      formatC(as.numeric(value), digits = 15, format = "fg", flag = "#"),
    text = gsub("[\r\n,]", " ", as.character(text))
  )
}

r2_record_channel <- function(case, channels, sm, time.varying) {
  for (chan in seq_along(channels)) {
    channel <- channels[[chan]]
    alph <- attr(channel, "alphabet")
    for (pos in seq_along(alph))
      r2_add(case, "channel_alphabet", chan, pos, text = alph[[pos]])
    has.miss <- any(channel == attr(channel, "nr"))
    if (is.character(sm[[chan]])) {
      cost <- suppressMessages(seqcost(
        channel, sm[[chan]], with.missing = has.miss,
        time.varying = time.varying
      ))$sm
    } else {
      cost <- sm[[chan]]
    }
    dims <- dim(cost)
    for (axis in seq_along(dims))
      r2_add(case, "channel_cost_dim", chan, axis, value = dims[[axis]])
  }
}

r2_record_mdseq <- function(case, out) {
  alph <- attr(out, "alphabet")
  for (pos in seq_along(alph))
    r2_add(case, "md_alphabet", 1L, pos, text = alph[[pos]])
  for (i in seq_len(nrow(out)))
    for (t in seq_len(ncol(out)))
      r2_add(case, "md_state", i, t, text = as.character(out[i, t]))
}

r2_record_cost <- function(case, out) {
  dims <- dim(out)
  for (axis in seq_along(dims))
    r2_add(case, "cat_dim", 1L, axis, value = dims[[axis]])
  alph <- attr(out, "alphabet")
  for (pos in seq_along(alph))
    r2_add(case, "md_alphabet", 1L, pos, text = alph[[pos]])
  if (length(dims) == 2L) {
    for (i in seq_len(dims[[1]]))
      for (j in seq_len(dims[[2]]))
        r2_add(case, "cat_sm", i, j, 1L, out[i, j])
  } else {
    for (t in seq_len(dims[[3]]))
      for (i in seq_len(dims[[1]]))
        for (j in seq_len(dims[[2]]))
          r2_add(case, "cat_sm", i, j, t, out[i, j, t])
  }
  cat.indel <- attr(out, "indel")
  if (is.null(cat.indel)) {
    r2_add(case, "cat_indel", text = "NULL")
  } else {
    for (pos in seq_along(cat.indel))
      r2_add(case, "cat_indel", 1L, pos, value = cat.indel[[pos]])
  }
  cweight <- attr(out, "cweight")
  for (pos in seq_along(cweight))
    r2_add(case, "cweight", 1L, pos, value = cweight[[pos]])
}

r2_record_distance <- function(case, out) {
  if (is.matrix(out)) {
    for (i in seq_len(nrow(out)))
      for (j in seq_len(ncol(out)))
        r2_add(case, "distance", i, j, value = out[i, j])
  } else {
    values <- as.numeric(out)
    for (pos in seq_along(values))
      r2_add(case, "distance", 1L, pos, value = values[[pos]])
  }
}

r2_capture <- function(case, expr, recorder) {
  warnings <- character()
  out <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      text <- conditionMessage(m)
      if (grepl("[!!]", text, fixed = TRUE))
        warnings <<- c(warnings, text)
      invokeRestart("muffleMessage")
    }
  )
  recorder(case, out)
  for (pos in seq_along(warnings))
    r2_add(case, "warning", 1L, pos, text = warnings[[pos]])
  invisible(out)
}

## Fixed-length channels.  U/V/R are declared but deliberately unobserved.
r2_a <- suppressMessages(seqdef(rbind(
  c("A", "B", "A", "B"), c("B", "A", "B", "A"), c("A", "A", "B", "B")
), alphabet = c("A", "B", "U")))
r2_b <- suppressMessages(seqdef(rbind(
  c("X", "Y", "X", "Y"), c("Y", "X", "Y", "X"), c("X", "X", "Y", "Y")
), alphabet = c("X", "Y", "V")))
r2_c <- suppressMessages(seqdef(rbind(
  c("P", "Q", "Q", "P"), c("Q", "P", "P", "Q"), c("P", "P", "Q", "Q")
), alphabet = c("P", "Q", "R")))

## One channel has an interior missing value; both have structural right voids.
r2_m_ragged <- suppressMessages(seqdef(rbind(
  c("A", NA, "B", NA), c("B", "A", "B", "A"), c("A", "B", NA, NA)
), alphabet = c("A", "B", "U")))
r2_b_ragged <- suppressMessages(seqdef(rbind(
  c("X", "Y", NA, NA), c("Y", "X", "Y", "X"), c("X", "X", "Y", NA)
), alphabet = c("X", "Y", "V")))

## Fixed-width missing fixture for DHD.  Its active cost width is the three
## declared states plus the observed missing state.
r2_m_fixed <- suppressMessages(seqdef(rbind(
  c("A", NA, "B", "A"), c("B", "A", "B", "A"), c("A", "B", "A", "B")
), alphabet = c("A", "B", "U")))

case <- "mdseq_ragged_2"
r2_capture(case, seqMD(list(r2_m_ragged, r2_b_ragged), what = "MDseq"), r2_record_mdseq)

case <- "cost_trate_declared_2"
sm <- rep("TRATE", 2L)
r2_record_channel(case, list(r2_a, r2_b), sm, FALSE)
r2_capture(case, seqMD(
  list(r2_a, r2_b), what = "cost", sm = sm, indel = c(1, 2),
  cweight = c(2, 1), link = "sum"
), r2_record_cost)

case <- "cost_indels_vector_mean"
sm <- list(
  suppressMessages(seqcost(r2_a, "INDELS"))$sm,
  suppressMessages(seqcost(r2_b, "INDELS"))$sm
)
r2_record_channel(case, list(r2_a, r2_b), sm, FALSE)
r2_capture(case, seqMD(
  list(r2_a, r2_b), what = "cost", sm = sm,
  indel = list(c(1, 2, 3), c(4, 5, 6)), cweight = c(2, 1), link = "mean"
), r2_record_cost)

case <- "cost_indelslog_auto"
sm <- rep("INDELSLOG", 2L)
r2_record_channel(case, list(r2_a, r2_b), sm, FALSE)
r2_capture(case, seqMD(
  list(r2_a, r2_b), what = "cost", sm = sm, indel = "auto"
), r2_record_cost)

case <- "cost_trate_declared_3"
sm <- rep("TRATE", 3L)
r2_record_channel(case, list(r2_a, r2_b, r2_c), sm, FALSE)
r2_capture(case, seqMD(
  list(r2_a, r2_b, r2_c), what = "cost", sm = sm,
  indel = 1.5, cweight = c(1, 2, 3)
), r2_record_cost)

case <- "dhd_trate_cost_missing"
sm <- list(
  suppressMessages(seqcost(
    r2_m_fixed, "TRATE", with.missing = TRUE, time.varying = TRUE
  ))$sm,
  suppressMessages(seqcost(r2_b, "TRATE", time.varying = TRUE))$sm
)
r2_record_channel(case, list(r2_m_fixed, r2_b), sm, TRUE)
r2_capture(case, seqMD(
  list(r2_m_fixed, r2_b), method = "DHD", what = "cost", sm = sm,
  indel = c(1, 2), with.missing = c(TRUE, FALSE), cweight = c(1, 2)
), r2_record_cost)

case <- "dhd_indels_cost_3"
sm <- lapply(
  list(r2_a, r2_b, r2_c),
  function(channel) suppressMessages(seqcost(
    channel, "INDELS", time.varying = TRUE
  ))$sm
)
r2_record_channel(case, list(r2_a, r2_b, r2_c), sm, TRUE)
r2_capture(case, seqMD(
  list(r2_a, r2_b, r2_c), method = "DHD", what = "cost", sm = sm,
  indel = c(1, 1, 1), cweight = c(1, 0.5, 2), link = "sum"
), r2_record_cost)

case <- "dhd_indelslog_cost_mean"
sm <- lapply(
  list(r2_a, r2_b),
  function(channel) suppressMessages(seqcost(
    channel, "INDELSLOG", time.varying = TRUE
  ))$sm
)
r2_record_channel(case, list(r2_a, r2_b), sm, TRUE)
r2_capture(case, seqMD(
  list(r2_a, r2_b), method = "DHD", what = "cost", sm = sm,
  indel = c(1, 1), cweight = c(2, 1), link = "mean"
), r2_record_cost)

make_dhd_sm <- function(k, width, scale) {
  out <- array(0, dim = c(k, k, width))
  for (t in seq_len(width))
    for (i in seq_len(k))
      for (j in seq_len(k))
        out[i, j, t] <- scale * abs(i - j) * t
  out
}
explicit_missing <- list(make_dhd_sm(4L, 4L, 0.5), make_dhd_sm(3L, 4L, 0.25))
case <- "dhd_explicit_cost_missing"
r2_record_channel(case, list(r2_m_fixed, r2_b), explicit_missing, TRUE)
r2_capture(case, seqMD(
  list(r2_m_fixed, r2_b), method = "DHD", what = "cost", sm = explicit_missing,
  indel = c(1, 2), with.missing = c(TRUE, FALSE), cweight = c(2, 1)
), r2_record_cost)

case <- "dhd_trate_diss_full"
r2_record_channel(case, list(r2_m_fixed, r2_b), rep("TRATE", 2L), TRUE)
r2_capture(case, seqMD(
  list(r2_m_fixed, r2_b), method = "DHD", what = "diss", sm = "TRATE",
  with.missing = c(TRUE, FALSE), full.matrix = TRUE
), r2_record_distance)

case <- "dhd_indels_diss_full_3"
r2_record_channel(case, list(r2_a, r2_b, r2_c), rep("INDELS", 3L), TRUE)
r2_capture(case, seqMD(
  list(r2_a, r2_b, r2_c), method = "DHD", what = "diss", sm = "INDELS",
  indel = 1, cweight = c(1, 0.5, 2), full.matrix = TRUE
), r2_record_distance)

case <- "dhd_indelslog_diss_condensed"
r2_record_channel(case, list(r2_a, r2_b), rep("INDELSLOG", 2L), TRUE)
r2_capture(case, seqMD(
  list(r2_a, r2_b), method = "DHD", what = "diss", sm = "INDELSLOG",
  indel = 1, cweight = c(2, 1), link = "mean", full.matrix = FALSE
), r2_record_distance)

explicit_three <- list(
  make_dhd_sm(3L, 4L, 0.5), make_dhd_sm(3L, 4L, 0.25),
  make_dhd_sm(3L, 4L, 0.75)
)
case <- "dhd_explicit_diss_condensed_3"
r2_record_channel(case, list(r2_a, r2_b, r2_c), explicit_three, TRUE)
r2_capture(case, seqMD(
  list(r2_a, r2_b, r2_c), method = "DHD", what = "diss", sm = explicit_three,
  indel = 1, cweight = c(1, 2, 1), full.matrix = FALSE
), r2_record_distance)

## Channel-specific declarations are checked before what="MDseq" returns, but
## TraMineR merely warns when they disagree with observed missingness.
case <- "warning_with_missing"
r2_capture(case, seqMD(
  list(r2_m_fixed, r2_b), what = "MDseq", with.missing = c(FALSE, TRUE)
), r2_record_mdseq)

r2_error_rows <- c("case=message")
r2_error <- function(case, expr) {
  message <- tryCatch(
    { suppressMessages(expr); "NO_ERROR" },
    error = function(e) conditionMessage(e)
  )
  message <- gsub("[\r\n]", " ", message)
  r2_error_rows <<- c(r2_error_rows, paste0(case, "=", message))
}
r2_error("dhd_state_width", seqMD(
  list(r2_m_fixed, r2_b), method = "DHD", what = "cost",
  sm = list(make_dhd_sm(3L, 4L, 1), make_dhd_sm(3L, 4L, 1)),
  with.missing = c(TRUE, FALSE)
))
r2_error("dhd_time_width", seqMD(
  list(r2_a, r2_b), method = "DHD", what = "cost",
  sm = list(make_dhd_sm(3L, 3L, 1), make_dhd_sm(3L, 3L, 1))
))
r2_error("stationary_state_width", seqMD(
  list(r2_a, r2_b), what = "cost",
  sm = list(matrix(1, 2, 2), matrix(1, 3, 3)), indel = 1
))
r2_error("indel_state_width", seqMD(
  list(r2_a, r2_b), what = "cost", sm = "TRATE",
  indel = list(c(1, 2), c(1, 2, 3))
))

csv_quote <- function(value) paste0('"', gsub('"', '""', value, fixed = TRUE), '"')
r2_lines <- c(
  "case,section,row,col,time,value,text",
  vapply(r2_rows, function(row) paste(vapply(row, csv_quote, character(1)), collapse = ","),
         character(1))
)
writeLines(r2_lines, file.path("tests", "r_goldens", "seqMD_r2_cases.csv"), useBytes = TRUE)
writeLines(r2_error_rows, file.path("tests", "r_goldens", "seqMD_r2_errors.txt"), useBytes = TRUE)
cat("SEQMD_R2_GOLDENS_DONE cases=", length(unique(vapply(r2_rows, `[[`, "", "case"))),
    " rows=", length(r2_rows), " errors=", length(r2_error_rows) - 1L, "\n", sep = "")
