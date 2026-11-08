# Golden generation: RP6b seqpolyads (linked-polyad analysis) -- DETERMINISTIC parts only.
# Pins observed.dist (per-polyad mean pairwise distance) and mean.dist["Obs"]. The
# Monte-Carlo random.dist / U / V part is stochastic (R RNG) and only structure-tested.
# Run in the live-R env:  Rscript scripts/r_goldens/make_goldens_seqpolyads.R
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
suppressMessages(library(TraMineRextras))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
stopifnot(as.character(packageVersion("TraMineRextras")) == "0.6.9")
out_dir <- file.path("tests", "r_goldens")

write_num_matrix <- function(name, m, header) {
  con <- file(file.path(out_dir, paste0(name, ".csv")), open = "wt")
  cat(header, file = con)
  txt <- apply(as.matrix(m), 1, function(r) paste(formatC(r, format = "g", digits = 15), collapse = ","))
  writeLines(txt, con)
  close(con)
}

alph <- c("A", "B", "C", "D")
S1 <- rbind(
  c("A", "A", "B", "C", "C", "C"),
  c("A", "B", "B", "B", "C", "A"),
  c("D", "D", "D", "D", "D", "D"),
  c("A", "B", "A", "B", "A", "B")
)
S2 <- rbind(
  c("B", "B", "A", "A", "D", "D"),
  c("C", "C", "C", "A", "A", "A"),
  c("A", "A", "B", "B", "C", "C"),
  c("D", "C", "B", "A", "D", "C")
)
S3 <- rbind(
  c("C", "A", "C", "A", "C", "A"),
  c("D", "D", "A", "A", "B", "B"),
  c("B", "B", "B", "C", "C", "D"),
  c("A", "B", "C", "D", "A", "B")
)
sd1 <- seqdef(S1, alphabet = alph)
sd2 <- seqdef(S2, alphabet = alph)
sd3 <- seqdef(S3, alphabet = alph)

# method="OM" (default indel=1, CONSTANT sm=2): traminex's seqdist OM matches R's raw OM.
res <- suppressMessages(seqpolyads(list(sd1, sd2, sd3),
  a = 1, method = "OM", sm = "CONSTANT", indel = 1))

write_num_matrix("seqpolyads_obs", t(as.vector(res$observed.dist)),
  "# R: seqpolyads(list(sd1,sd2,sd3), a=1, OM)$observed.dist\n")
write_num_matrix("seqpolyads_meanobs", matrix(res$mean.dist["Obs"], 1, 1),
  "# R: ...$mean.dist['Obs']\n")

# method="HAM" is R's default; traminex's seqdist HAM is now R-faithful (raw count), so
# seqpolyads(method="HAM") matches R exactly too (previously diverged by 1/L -- audit P1).
resH <- suppressMessages(seqpolyads(list(sd1, sd2, sd3), a = 1, method = "HAM"))
write_num_matrix("seqpolyads_obs_ham", t(as.vector(resH$observed.dist)),
  "# R: seqpolyads(list(sd1,sd2,sd3), a=1, HAM)$observed.dist\n")
write_num_matrix("seqpolyads_meanobs_ham", matrix(resH$mean.dist["Obs"], 1, 1),
  "# R: ...HAM...$mean.dist['Obs']\n")

cat("SEQPOLYADS_GOLDENS_DONE\n")

## R2: complete controls, a=2, seed defect, and normalized timing transcript --
## Keep core="auto" deterministic and lightweight in the oracle process: the
## fallback reports two logical cores, so upstream resolves auto to core=1.
options(
  parallelly.availableCores.methods = "fallback",
  parallelly.availableCores.fallback = 2
)

r2_rows <- list()
r2_add <- function(case, section, index = 0L, value = NaN, text = "") {
  r2_rows[[length(r2_rows) + 1L]] <<- c(
    case = as.character(case), section = as.character(section), index = as.character(index),
    value = if (is.nan(value)) "NaN" else if (is.na(value)) "NA" else
      formatC(as.numeric(value), digits = 15, format = "fg", flag = "#"),
    text = gsub("[\r\n,]", " ", as.character(text))
  )
}

r2_messages <- c("case=message")
r2_record <- function(case, out, config) {
  for (name in names(config)) r2_add(case, "config", text = paste0(name, "=", config[[name]]))
  components <- list(
    mean_dist = out$mean.dist, U = out$U, U_tp = out$U.tp, V = out$V,
    V_95 = out$V.95, observed = out$observed.dist, random = out$random.dist
  )
  for (name in names(components))
    for (i in seq_along(components[[name]]))
      r2_add(case, name, i, components[[name]][[i]])
}

r2_run <- function(case, channels, ..., show.time = FALSE) {
  args <- list(...)
  args$seqlist <- channels
  args$show.time <- show.time
  messages <- character()
  stdout <- capture.output(
    out <- withCallingHandlers(
      do.call(seqpolyads, args),
      message = function(m) {
        messages <<- c(messages, trimws(conditionMessage(m)))
        invokeRestart("muffleMessage")
      }
    )
  )
  config <- c(
    a = if (is.null(args$a)) 1 else args$a,
    method = if (is.null(args$method)) "HAM" else args$method,
    roles = length(channels), width = ncol(channels[[1]]),
    seed_arg = if (is.null(args$s)) "NULL" else args$s,
    T = if (is.null(args$T)) 1000 else args$T,
    core = if (is.null(args$core)) 1 else args$core,
    replace = if (is.null(args$replace)) TRUE else args$replace,
    weighted = if (is.null(args$weighted)) TRUE else args$weighted,
    with_missing = if (is.null(args$with.missing)) FALSE else args$with.missing,
    rand_weight_type = if (is.null(args$rand.weight.type)) 1 else args$rand.weight.type
  )
  r2_record(case, out, config)
  for (message in messages)
    r2_messages <<- c(r2_messages, paste0(case, "=", message))
  if (show.time) {
    timing <- trimws(stdout[nzchar(trimws(stdout))])
    if (length(timing) >= 1L)
      r2_messages <<- c(r2_messages, paste0(case, "=timing_header:", timing[[1]]))
    if (length(timing) >= 2L)
      r2_messages <<- c(r2_messages, paste0(case, "=timing_values:<elapsed>"))
  }
  invisible(out)
}

rw <- c(1, 2, 3, 4)
r2_sd1 <- suppressMessages(seqdef(S1, alphabet = alph, weights = rw))
r2_sd2 <- suppressMessages(seqdef(S2, alphabet = alph, weights = rw))
r2_sd3 <- suppressMessages(seqdef(S3, alphabet = alph, weights = rw))
r2_two <- list(r2_sd1, r2_sd2)
r2_three <- list(r2_sd1, r2_sd2, r2_sd3)

r2_run("a1_ham_two", r2_two, a = 1, method = "HAM", T = 12, s = 11)
r2_run(
  "a1_om_three_pair_coreauto", r2_three, a = 1, method = "OM",
  w = c(1, 2, 0.5), indel = 1, sm = "CONSTANT", weighted = FALSE,
  T = 12, s = 22, core = "auto"
)
r2_run(
  "a1_weighted_type1", r2_three, a = 1, method = "HAM",
  w = c(2, 1, 3), role.weights = c(0.2, 0.3, 0.5),
  weighted = TRUE, rand.weight.type = 1, T = 12, s = 33
)
r2_run(
  "a1_weighted_type2", r2_three, a = 1, method = "HAM",
  w = c(2, 1, 3), role.weights = c(0.2, 0.3, 0.5),
  weighted = TRUE, rand.weight.type = 2, T = 12, s = 33
)
r2_run(
  "a2_ham_replace_true", r2_two, a = 2, method = "HAM",
  replace = TRUE, T = 12, s = 44
)
r2_run(
  "a2_ham_replace_false", r2_three, a = 2, method = "HAM",
  w = c(1, 2, 1), replace = FALSE, weighted = FALSE, T = 12, s = 55
)
r2_run(
  "a2_om_weighted", r2_three, a = 2, method = "OM",
  indel = 1, sm = "CONSTANT", role.weights = c(0.5, 0.3, 0.2),
  rand.weight.type = 2, replace = TRUE, weighted = TRUE, T = 12, s = 66
)

S1m <- S1
S2m <- S2
S1m[1, 3] <- NA
S1m[4, 6] <- NA
S2m[2, 2] <- NA
S2m[4, 6] <- NA
r2_m1 <- suppressMessages(seqdef(S1m, alphabet = alph, weights = rw))
r2_m2 <- suppressMessages(seqdef(S2m, alphabet = alph, weights = rw))
r2_run(
  "a2_missing_void", list(r2_m1, r2_m2), a = 2, method = "OM",
  indel = 1, sm = "CONSTANT", replace = TRUE, with.missing = TRUE, T = 8, s = 77
)

## Upstream overwrites s with sequence width before set.seed(s).
r2_run("seed_defect_s1", r2_two, a = 1, method = "HAM", T = 10, s = 1)
r2_run("seed_defect_s999", r2_two, a = 1, method = "HAM", T = 10, s = 999)
r2_sd1_w5 <- suppressMessages(seqdef(S1[, 1:5], alphabet = alph, weights = rw))
r2_sd2_w5 <- suppressMessages(seqdef(S2[, 1:5], alphabet = alph, weights = rw))
r2_run(
  "seed_width5", list(r2_sd1_w5, r2_sd2_w5), a = 1, method = "HAM",
  T = 10, s = 1
)
r2_run(
  "show_time", r2_two, a = 1, method = "HAM", T = 4, s = 88,
  show.time = TRUE
)
# Four roles with pair weights 1..6 (TX-AUD-20261019-011): R pairs them with upper.tri's
# column-wise order (1,2), (1,3), (2,3), (1,4), (2,4), (3,4) (seqpolyads.R:120, 154, 176), where the
# Rd (item w) states the row-wise order; owner decision D3 keeps the Rd's order, so these
# observed distances (6 and 6) are R's, recorded to keep the divergence visible.
four_role <- function(rows) suppressMessages(seqdef(do.call(rbind, lapply(rows, function(r) strsplit(r, "")[[1]])),
                                                    alphabet = c("A", "B")))
r2_four <- list(four_role(c("AAA", "AAA")), four_role(c("AAB", "AAA")), four_role(c("ABB", "BBB")),
                four_role(c("BBB", "ABA")))
r2_run("four_roles_pair_weights", r2_four, a = 1, method = "HAM", w = 1:6, weighted = FALSE, T = 4, s = 99)

r2_errors <- c("case=message")
r2_error <- function(case, expr) {
  message <- tryCatch(
    { suppressMessages(expr); "NO_ERROR" },
    error = function(e) conditionMessage(e)
  )
  r2_errors <<- c(r2_errors, paste0(case, "=", trimws(gsub("[\r\n]", " ", message))))
}
r2_error("bad_a", seqpolyads(r2_two, a = 3, T = 2))
r2_error("bad_core", seqpolyads(r2_two, core = "many", T = 2))
r2_error("bad_role_weights", seqpolyads(r2_three, role.weights = c(1, 2), T = 2))
r2_error("not_list", seqpolyads(r2_sd1, T = 2))
r2_error("not_sts", seqpolyads(list(r2_sd1, S2), T = 2))
r2_error("mismatched_width", seqpolyads(list(r2_sd1, r2_sd2_w5), T = 2))

csv_quote <- function(value) paste0('"', gsub('"', '""', value, fixed = TRUE), '"')
r2_lines <- c(
  "case,section,index,value,text",
  vapply(r2_rows, function(row) paste(vapply(row, csv_quote, character(1)), collapse = ","),
         character(1))
)
writeLines(r2_lines, file.path(out_dir, "seqpolyads_r2_cases.csv"), useBytes = TRUE)
writeLines(r2_errors, file.path(out_dir, "seqpolyads_r2_errors.txt"), useBytes = TRUE)
writeLines(r2_messages, file.path(out_dir, "seqpolyads_r2_messages.txt"), useBytes = TRUE)
cat(
  "SEQPOLYADS_R2_GOLDENS_DONE cases=", length(unique(vapply(r2_rows, `[[`, "", "case"))),
  " rows=", length(r2_rows), " errors=", length(r2_errors) - 1L,
  " messages=", length(r2_messages) - 1L, "\n", sep = ""
)
