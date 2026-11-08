stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
# Live-R goldens for batch 14's sequence-construction, format and multidomain items (rulings R14-a to R14-m): what
# TraMineR and TraMineRextras return for the inputs the tests use, including R's errors where traminex corrects R.
#
# Each row is case, field, value. Vectors are joined by spaces (numbers in %.17g, NA as NA); text vectors by ' ;; '.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_construction_items.R
suppressMessages(library(TraMineR)); suppressMessages(library(TraMineRextras))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14", as.character(packageVersion("TraMineRextras")) == "0.6.9")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
rows <- list()
add <- function(case, field, value) rows[[length(rows) + 1L]] <<- data.frame(case = case, field = field, value = value)
num <- function(x) paste(ifelse(is.na(x), "NA", sprintf("%.17g", x)), collapse = " ")
txt <- function(x) paste(ifelse(is.na(x), "NA", x), collapse = " ;; ")

# --- seqrecode's otherwise (R14-a: TX-AUD-20261019-033). R recodes every cell recodes does not list, the missing and
# void markers included (seqrecode.R:69-77), so the void padding becomes the otherwise state and the shorter
# sequence lengthens. Recorded: R's cells (seqconc), lengths and alphabet, with and without the missing marker named.
recode_cells <- rbind(c("a", "b", NA, "b"), c("a", "c", NA, NA))
recode_s <- suppressMessages(seqdef(recode_cells, right = "DEL"))
add("recode", "cells", txt(apply(recode_cells, 1, paste, collapse = " ")))
add("recode", "input_lengths", num(seqlength(recode_s)))
recode_cases <- list(recode_otherwise = list("A" = "a"), recode_missing_named = list("A" = "a", "M" = "*"))
for (case in names(recode_cases)) {
  r <- suppressMessages(seqrecode(recode_s, recodes = recode_cases[[case]], otherwise = "O"))
  add(case, "sequences", txt(seqconc(r)))
  add(case, "lengths", num(seqlength(r)))
  add(case, "alphabet", txt(alphabet(r)))
}

# --- Labels and colours of long data (R14-b: TX-AUD-20261019-035, TX-REP-20261019-013). R's route for long data is
# reshape to wide, then seqdef, which pairs labels and cpal with the real states (one entry each) and gives the
# markers their own (missing.color; void is left undrawn). Recorded: the alphabet, stlab, the cpal attribute and
# missing.color, with seqdef's defaults and with labels given for the real states only.
long_cases <- list(
  long_defaults = list(long = data.frame(id = c(1, 1, 1, 2, 2), time = c(1, 2, 3, 1, 2), state = c("b", "a", "c", "a", "c")),
                       labels = NULL),
  long_real_labels = list(long = data.frame(id = c(1, 1, 1, 2, 2, 2), time = c(1, 2, 3, 1, 2, 3),
                                            state = c("A", NA, "B", "B", "C", "A")),
                          labels = c("la", "lb", "lc")))
for (case in names(long_cases)) {
  spec <- long_cases[[case]]
  add(case, "id", num(spec$long$id))
  add(case, "time", num(spec$long$time))
  add(case, "state", txt(spec$long$state))
  wide <- reshape(spec$long, direction = "wide", idvar = "id", timevar = "time")
  s <- if (is.null(spec$labels)) suppressMessages(seqdef(wide[, -1])) else suppressMessages(seqdef(wide[, -1], labels = spec$labels))
  if (!is.null(spec$labels)) add(case, "labels_given", txt(spec$labels))
  add(case, "alphabet", txt(alphabet(s)))
  add(case, "stlab", txt(stlab(s)))
  add(case, "cpal", txt(attr(s, "cpal")))
  add(case, "missing_color", attr(s, "missing.color"))
  add(case, "cells", txt(seqconc(s)))
}

# --- The persons' order in long data (R14-d: TX-AUD-20261019-039). R's route (reshape to wide, then seqdef) and
# TSE_to_STS's event mode keep the persons in order of first appearance. Recorded: the ids of each, for text ids and
# for numeric ids.
order_cases <- list(order_text = c("p10", "p10", "p2", "p2", "p1", "p1"), order_numeric = c(10, 10, 2, 2, 1, 1))
for (case in names(order_cases)) {
  ids <- order_cases[[case]]
  long <- data.frame(id = ids, time = rep(c(1, 2), 3), state = c("A", "B", "B", "B", "A", "A"))
  add(case, "id", if (is.numeric(ids)) num(ids) else txt(ids))
  add(case, "time", num(long$time))
  add(case, "state", txt(long$state))
  wide <- reshape(long, direction = "wide", idvar = "id", timevar = "time")
  s <- suppressMessages(seqdef(wide[, -1], id = wide$id))
  add(case, "route_ids", txt(rownames(s)))
  ev <- data.frame(id = unique(ids), time = c(0, 0, 0), event = c("A", "B", "A"))
  sts <- suppressMessages(TSE_to_STS(ev, id = "id", timestamp = "time", event = "event", tmax = 2))
  add(case, "event_ids", txt(rownames(sts)))
}

# --- The nr marker as a replacement state (R14-e: TX-AUD-20261019-040). R's seqdef accepts left, gaps or right equal
# to the nr marker "*" and lists it as a state of the alphabet, twice when two of them name it; traminex refuses.
# Recorded: R's alphabet and cells for the audit's data.
policy_d <- data.frame(a = c(NA, "A"), b = c("A", "B"), c = c("B", NA), stringsAsFactors = FALSE)
add("policy", "cells", txt(apply(as.matrix(policy_d), 1, function(r) paste(ifelse(is.na(r), "NA", r), collapse = " "))))
policy_cases <- list(policy_left = list(left = "*"), policy_gaps_left = list(gaps = "*", left = "*"),
                     policy_right = list(right = "*"))
for (case in names(policy_cases)) {
  s <- suppressMessages(do.call(seqdef, c(list(policy_d), policy_cases[[case]])))
  add(case, "options", txt(paste(names(policy_cases[[case]]), "=", unlist(policy_cases[[case]]))))
  add(case, "alphabet", txt(alphabet(s)))
  add(case, "cells", txt(seqconc(s)))
}

# --- Mixed numeric and text columns (R14-h: TX-AUD-20261019-062). R's seqdef turns a data frame into a matrix with
# as.matrix, which pads the numbers of a numeric column (" 1", " 2" beside "10") or makes NA cells, so 1 and "1" are
# different states or a cell is lost; traminex reads each cell with R's as.character. Recorded: R's alphabet and cells.
mixed_cases <- list(
  mixed_int_chr = data.frame(a = c(1L, 2L, 10L), b = c("1", "2", "x")),
  mixed_int_chr2 = data.frame(a = c(1L, 10L, 2L), b = c("10", "1", "2")),
  mixed_dbl_chr = data.frame(a = c(1, 2.5, 3), b = c("1", "2.5", "3")))
for (case in names(mixed_cases)) {
  d <- mixed_cases[[case]]
  add(case, "column_a", num(d$a))
  add(case, "column_b", txt(d$b))
  s <- suppressMessages(seqdef(d))
  add(case, "alphabet", txt(alphabet(s)))
  add(case, "cells", txt(seqconc(s)))
}

# --- Numeric states of long data (R14-j: TX-AUD-20261010-073). R's route (reshape, then seqdef) prints numeric
# states with as.character: 1 and 2, not 1.0 and 2.0. Recorded: the alphabet for whole and fractional floats.
float_cases <- list(float_whole = c(1.0, 2.0, 2.0, 1.0), float_fraction = c(1.0, 2.5, 2.5, 1.0))
for (case in names(float_cases)) {
  long <- data.frame(id = c(1, 1, 2, 2), time = c(1, 2, 1, 2), state = float_cases[[case]])
  add(case, "state", num(long$state))
  wide <- reshape(long, direction = "wide", idvar = "id", timevar = "time")
  add(case, "alphabet", txt(alphabet(suppressMessages(seqdef(wide[, -1])))))
}

# --- seqMD with cost-method strings and per-state list indels (R14-k: TX-AUD-20261011-060). With sm given as
# strings, R's seqMD never copies a list indel into indel_list (seqMD.R:119-124, 226-264), so checkcost indexes an
# empty list and R stops; with the same costs given as matrices it copies them (seqMD.R:261-264) and computes.
# Recorded: R's error on the string route and R's distances on the matrix route, for the gap audit's channels.
md_c1 <- suppressMessages(seqdef(data.frame(a = c("A", "B", "A"), b = c("B", "B", "A"), c = c("A", "A", "B")), alphabet = c("A", "B")))
md_c2 <- suppressMessages(seqdef(data.frame(a = c("X", "Y", "Y"), b = c("X", "X", "Y"), c = c("Y", "X", "X")), alphabet = c("X", "Y")))
add("md_indel", "channel_1", txt(seqconc(md_c1)))
add("md_indel", "channel_2", txt(seqconc(md_c2)))
md_indel <- list(c(2, 1), c(1, 1))
add("md_indel", "indel_1", num(md_indel[[1]]))
add("md_indel", "indel_2", num(md_indel[[2]]))
md_string <- tryCatch(suppressMessages(seqMD(list(md_c1, md_c2), method = "OM", sm = list("CONSTANT", "CONSTANT"),
                                             indel = md_indel, what = "diss")),
                      error = function(e) paste("ERROR:", conditionMessage(e)))
add("md_indel", "string_route", if (is.character(md_string)) md_string else "computed")
md_sm <- list(suppressMessages(seqcost(md_c1, "CONSTANT"))$sm, suppressMessages(seqcost(md_c2, "CONSTANT"))$sm)
md_matrix <- suppressMessages(seqMD(list(md_c1, md_c2), method = "OM", sm = md_sm, indel = md_indel, what = "diss"))
add("md_indel", "matrix_route_diss", num(as.vector(as.matrix(md_matrix))))
# The same route with one indel per channel (batch-14 review M-3): the test needs per-state costs whose distances
# differ from these, or it cannot tell that the per-state costs are used.
md_scalar <- suppressMessages(seqMD(list(md_c1, md_c2), method = "OM", sm = md_sm, indel = c(1, 1), what = "diss"))
add("md_indel", "matrix_route_scalar_diss", num(as.vector(as.matrix(md_scalar))))

# --- seqstart's empty and inverted windows (R14-m: TX-AUD-20261011-062). With no rows R's min and max of the empty
# index are Inf and -Inf (seqstart.R:18-25), so the default tmax stops with "no data" and a default tmin leaves a frame
# matrix() refuses; a frame from tmin to tmax of negative length is refused by matrix() too (seqstart.R:31), while one
# of length 0 (tmin = tmax + 1) is returned. Recorded: R's error or result dimensions, with data.start = new.start = 1.
start_d <- data.frame(a = c("A", "B"), b = c("B", "A"))
add("seqstart", "cells", txt(apply(as.matrix(start_d), 1, paste, collapse = " ")))
start_cases <- list(seqstart_empty = list(rows = 0, tmin = NULL, tmax = NULL),
                    seqstart_empty_tmax = list(rows = 0, tmin = NULL, tmax = 2),
                    seqstart_inverted = list(rows = 2, tmin = 3, tmax = 1),
                    seqstart_zero_width = list(rows = 2, tmin = 2, tmax = 1))
for (case in names(start_cases)) {
  spec <- start_cases[[case]]
  add(case, "rows", num(spec$rows))
  add(case, "tmin", if (is.null(spec$tmin)) "NULL" else num(spec$tmin))
  add(case, "tmax", if (is.null(spec$tmax)) "NULL" else num(spec$tmax))
  r <- tryCatch(suppressWarnings(suppressMessages(seqstart(start_d[seq_len(spec$rows), ], data.start = 1, new.start = 1,
                                                          tmin = spec$tmin, tmax = spec$tmax))),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  add(case, "result", if (is.character(r) && length(r) == 1 && startsWith(r, "ERROR")) r else paste("dim", num(dim(r))))
}

# --- A missing id in long data (batch-14 review I-2, TX-AUD-20261019-039). R's route reshapes the person with no id
# like the others, then seqdef with the ids as row names stops. Recorded: the ids and R's error.
missing_id_long <- data.frame(id = c(3, 3, NA, NA, 1, 1), time = rep(1:2, 3), state = c("A", "B", "B", "B", "A", "A"))
add("long_missing_id", "id", num(missing_id_long$id))
add("long_missing_id", "time", num(missing_id_long$time))
add("long_missing_id", "state", txt(missing_id_long$state))
missing_id_wide <- reshape(missing_id_long, direction = "wide", idvar = "id", timevar = "time")
missing_id_s <- tryCatch(suppressMessages(seqdef(missing_id_wide[, -1], id = missing_id_wide$id)),
                         error = function(e) paste("ERROR:", conditionMessage(e)))
add("long_missing_id", "result", if (is.character(missing_id_s)) missing_id_s else "computed")

write.csv(do.call(rbind, rows), file.path(out_dir, "construction_items.csv"), row.names = FALSE, quote = TRUE)
