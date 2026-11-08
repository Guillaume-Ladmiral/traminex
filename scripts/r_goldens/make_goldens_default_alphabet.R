stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
# Live-R goldens for the default alphabet order (owner decision D2 of 2026-10-19; audit findings
# TX-AUD-20261019-071, -038 and -013).
#
# Without `alphabet=`, TraMineR's seqdef takes `plevels <- as.character(seqstatl(seqdata))`
# (seqdef.R:59-66): seqstatl sorts the distinct values and re-sorts them numerically when every one
# converts with as.numeric (seqstatl.R:21-30); a `missing=` code takes part in that decision
# (seqdef.R:59-62); replacement states of left/gaps/right follow (seqdef.R:116-132). Every traminex
# constructor that derives an alphabet now uses that order, so the tests compare alphabets by
# position. The cases cover numeric, text (case, punctuation, accents), mixed and missing-code data,
# the SPS, SPELL (text and factor status), long (base-R reshape) and TSE event routes, seqrecode,
# seqdss and seqMD's combined alphabet (finding 038), the tie rules of seqmodst and seqgranularity
# (finding 013), and order-dependent indicators of a default object.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_default_alphabet.R

suppressMessages({
  library(TraMineR)
  library(TraMineRextras)
})
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
stopifnot(as.character(packageVersion("TraMineRextras")) == "0.6.9")
stopifnot(identical(icuGetCollate(), "root"))

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_rows <- list()
result_rows <- list()

cell_text <- function(value) ifelse(is.na(value), "<NA>", as.character(value))
num_text <- function(value) ifelse(is.na(value), "NA", sprintf("%.17g", value))

record_input <- function(case, table, frame) {
  frame <- as.data.frame(frame, stringsAsFactors = FALSE)
  for (column in names(frame)) {
    input_rows[[length(input_rows) + 1L]] <<- data.frame(
      case = case, table = table, row = seq_len(nrow(frame)), column = column,
      value = cell_text(as.character(frame[[column]])), stringsAsFactors = FALSE
    )
  }
}

record_result <- function(case, route, kind, values) {
  values <- if (is.numeric(values)) num_text(values) else cell_text(values)
  result_rows[[length(result_rows) + 1L]] <<- data.frame(
    case = case, route = route, kind = kind, position = seq_along(values), value = values,
    stringsAsFactors = FALSE
  )
}

quiet <- function(expr) suppressWarnings(suppressMessages(expr))

sts_case <- function(case, rows, ...) {
  frame <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(frame) <- paste0("t", seq_len(ncol(frame)))
  record_input(case, "sts", frame)
  seq <- quiet(TraMineR::seqdef(frame, ...))
  record_result(case, "seqdef", "alphabet", alphabet(seq))
  invisible(seq)
}

# --- STS data -----------------------------------------------------------------------------------
sts_case("numeric", list(c("10", "2", "1"), c("2", "9", "10")))
sts_case("text_case", list(c("b", "B", "a"), c("A", "_x", "b")))
sts_case("text_accents", list(c("é", "e", "f"), c("ä", "a", "z")))
sts_case("text_punct", list(c("a.b", "a_b", "aB"), c("a b", "a.b", "aB")))
sts_case("mixed_number_text", list(c("1", "a", "10"), c("2", "1", "a")))
sts_case("missing_text_code", list(c("1", "10", "x"), c("2", "1", "10")), missing = "x")
sts_case("missing_numeric_code", list(c("1", "10", "99"), c("2", "1", "10")), missing = "99")
sts_case("left_policy", list(c(NA, "c", "a"), c("b", "c", "a")), left = "LL")
sts_case("right_na", list(c("B", NA, NA), c("A", "B", NA)), right = NA)
# A logical column beside an integer one: R's seqdef reads the data frame as numbers (TRUE as 1).
promotion <- data.frame(logical = c(TRUE, FALSE), integer = c(1L, 2L))
record_input("logical_promotion", "frame", promotion)
record_result("logical_promotion", "seqdef", "alphabet", alphabet(quiet(TraMineR::seqdef(promotion))))
# A numeric factor column: seqdef takes the factor's levels, which R prints as "1" and "2" for the
# doubles 1 and 2; the codes are the level positions (TX-AUD-20261019-088: a test claimed this from R
# without a golden).
numeric_factor <- data.frame(state = factor(c(1, 2, 1)))
record_input("numeric_factor", "frame", numeric_factor)
numeric_factor_seq <- quiet(TraMineR::seqdef(numeric_factor))
record_result("numeric_factor", "seqdef", "alphabet", alphabet(numeric_factor_seq))
record_result("numeric_factor", "seqdef", "codes", as.integer(numeric_factor_seq[[1]]))

# --- the census frame of lead A4-05: first appearance b, a, c -------------------------------------
census <- sts_case("census_wide", list(c("b", "a", "c"), c("a", "c", "c")))

sps <- data.frame(x = c("(b,1)-(a,1)-(c,1)", "(a,1)-(c,2)"), stringsAsFactors = FALSE)
record_input("sps", "sps", sps)
record_result("sps", "seqdef_SPS", "alphabet", alphabet(quiet(TraMineR::seqdef(sps, informat = "SPS"))))
record_result("sps", "seqformat_SPS_STS", "alphabet",
              alphabet(quiet(TraMineR::seqdef(TraMineR::seqformat(sps, from = "SPS", to = "STS")))))

spell <- data.frame(id = c(1, 1, 1, 2, 2), begin = c(1, 2, 3, 1, 2), end = c(1, 2, 3, 1, 3),
                    status = c("b", "a", "c", "a", "c"), stringsAsFactors = FALSE)
record_input("spell_text", "spell", spell)
record_result("spell_text", "seqdef_SPELL", "alphabet",
              alphabet(quiet(TraMineR::seqdef(spell, informat = "SPELL", process = FALSE))))
spell_factor <- spell
spell_factor$status <- factor(spell$status, levels = c("c", "b", "a"))
record_input("spell_factor", "spell", spell)
record_input("spell_factor", "levels", data.frame(level = levels(spell_factor$status)))
record_result("spell_factor", "seqdef_SPELL", "alphabet",
              alphabet(quiet(TraMineR::seqdef(spell_factor, informat = "SPELL", process = FALSE))))
formatted <- quiet(TraMineR::seqformat(spell_factor, from = "SPELL", to = "STS", process = FALSE))
record_result("spell_factor", "seqformat_SPELL_STS", "levels", levels(formatted[, 1]))

long <- data.frame(id = c(1, 1, 1, 2, 2, 2), time = c(1, 2, 3, 1, 2, 3),
                   state = c("b", "a", "c", "a", "c", "c"), stringsAsFactors = FALSE)
record_input("long", "long", long)
reshaped <- reshape(long, idvar = "id", timevar = "time", direction = "wide")
record_result("long", "seqdef_reshape", "alphabet",
              alphabet(quiet(TraMineR::seqdef(reshaped[, -1, drop = FALSE]))))
# A long table holding a state spelled like the nr marker (batch-7 review M-3). R's own route
# lists * as a state yet counts its cells as missing; seqdef_from_long and TSE_to_STS keep * a
# real state (TX-AUD-20261019-005), so their lengths differ from these.
for (star_case in list(list("long_star", c("A", "*", "B", "A", "A", "B")),
                       list("long_star_gap", c("A", "*", "B", "A", NA, "B")))) {
  star <- data.frame(id = c(1, 1, 1, 2, 2, 2), time = c(1, 2, 3, 1, 2, 3), state = star_case[[2]],
                     stringsAsFactors = FALSE)
  record_input(star_case[[1]], "long", star)
  wide <- quiet(TraMineR::seqdef(reshape(star, idvar = "id", timevar = "time", direction = "wide")[, -1, drop = FALSE]))
  record_result(star_case[[1]], "seqdef_reshape", "alphabet", alphabet(wide))
  record_result(star_case[[1]], "seqdef_reshape", "seqlength_without_missing",
                as.numeric(TraMineR::seqlength(wide, with.missing = FALSE)))
}

events <- data.frame(id = c(1, 1, 2), time = c(0, 1, 0), event = c("eb", "ea", "ea"),
                     stringsAsFactors = FALSE)
record_input("tse_event", "events", events)
tse <- quiet(TraMineRextras::TSE_to_STS(events, id = 1, timestamp = 2, event = 3, tmax = 3))
record_result("tse_event", "TSE_to_STS_seqdef", "alphabet", alphabet(quiet(TraMineR::seqdef(tse))))

# The state grid traminex's events_to_state_sequences builds for p1: Zeta@0, alpha@1, Zeta@2 and
# p2: alpha@1, with fill "_" (code points give Zeta, _, alpha; R's collation _, alpha, Zeta).
grid <- data.frame(t0 = c("Zeta", "_"), t1 = c("alpha", "alpha"), t2 = c("Zeta", "alpha"),
                   stringsAsFactors = FALSE)
record_input("events_grid", "sts", grid)
record_result("events_grid", "seqdef", "alphabet", alphabet(quiet(TraMineR::seqdef(grid))))

record_result("census_wide", "seqrecode_Z_c", "alphabet",
              alphabet(quiet(TraMineR::seqrecode(census, recodes = list(Z = "c")))))
record_result("census_wide", "seqdss", "alphabet", alphabet(quiet(TraMineR::seqdss(census))))

# --- seqMD's combined alphabet (finding 038) -------------------------------------------------------
md_a <- quiet(TraMineR::seqdef(rbind(c("10", "2", "1"), c("2", "10", "1"))))
md_b <- quiet(TraMineR::seqdef(rbind(c("1", "1", "1"), c("1", "1", "1"))))
record_input("seqmd_numeric", "channel1", rbind(c("10", "2", "1"), c("2", "10", "1")))
record_input("seqmd_numeric", "channel2", rbind(c("1", "1", "1"), c("1", "1", "1")))
record_result("seqmd_numeric", "seqMD_MDseq", "alphabet",
              alphabet(quiet(TraMineR::seqMD(list(md_a, md_b), what = "MDseq", ch.sep = "."))))
record_result("seqmd_numeric", "seqMD_cost", "rownames",
              rownames(quiet(TraMineR::seqMD(list(md_a, md_b), what = "cost", sm = "CONSTANT",
                                             ch.sep = "."))))
md_c <- quiet(TraMineR::seqdef(rbind(c("B", "a", "c"), c("a", "c", "c"))))
md_d <- quiet(TraMineR::seqdef(rbind(c("x", "Y", "x"), c("Y", "x", "x"))))
record_input("seqmd_text", "channel1", rbind(c("B", "a", "c"), c("a", "c", "c")))
record_input("seqmd_text", "channel2", rbind(c("x", "Y", "x"), c("Y", "x", "x")))
record_result("seqmd_text", "seqMD_MDseq", "alphabet",
              alphabet(quiet(TraMineR::seqMD(list(md_c, md_d), what = "MDseq"))))

# --- tie rules (finding 013) ------------------------------------------------------------------------
ties <- sts_case("ties", list(c("b", "b", "a"), c("a", "b", "b")))
row_strings <- function(seq) {
  m <- as.matrix(seq)
  apply(m, 1, function(r) paste(r[!(r %in% c(attr(seq, "void")))], collapse = "-"))
}
record_result("ties", "seqmodst", "string", row_strings(quiet(TraMineR::seqmodst(ties))))
record_result("ties", "seqgranularity_mostfreq_2", "string",
              row_strings(quiet(TraMineRextras::seqgranularity(ties, tspan = 2, method = "mostfreq"))))

# --- order-dependent results of a default object ---------------------------------------------------
indicators <- sts_case("indicators", list(
  c("d", "d", "b", "b", "a", "a", "c"),
  c("b", "a", "a", "c", "c", "d", "d"),
  c("a", "a", "a", "b", "b", "b", "b"),
  c("c", "c", "d", "d", "a", "b", "b"),
  c("d", "b", "b", "a", "c", "c", "c"),
  c("b", "b", "d", "d", "d", "a", "a")
))
record_result("indicators", "seqibad", "value", as.numeric(quiet(TraMineR::seqibad(indicators))))
record_result("indicators", "seqprecstart", "value", as.numeric(quiet(TraMineR::seqprecstart(indicators))))
record_result("indicators", "seqidegrad", "value", as.numeric(quiet(TraMineR::seqidegrad(indicators))))
record_result("indicators", "seqprecarity", "value", as.numeric(quiet(TraMineR::seqprecarity(indicators))))
record_result("indicators", "seqinsecurity", "value", as.numeric(quiet(TraMineR::seqinsecurity(indicators))))
record_result("indicators", "seqmodst", "string", row_strings(quiet(TraMineR::seqmodst(indicators))))
record_result("indicators", "seqgranularity_mostfreq_2", "string",
              row_strings(quiet(TraMineRextras::seqgranularity(indicators, tspan = 2, method = "mostfreq"))))
record_result("indicators", "seqtrate", "rownames", rownames(quiet(TraMineR::seqtrate(indicators))))

write.csv(do.call(rbind, input_rows), file.path(out_dir, "default_alphabet_inputs.csv"),
          row.names = FALSE, quote = TRUE, fileEncoding = "UTF-8")
write.csv(do.call(rbind, result_rows), file.path(out_dir, "default_alphabet.csv"),
          row.names = FALSE, quote = TRUE, fileEncoding = "UTF-8")
