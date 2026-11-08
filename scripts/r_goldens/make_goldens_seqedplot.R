# Golden generation: seqedplot producing SERIES (TraMineRextras seqedplot.R).
# seqedplot is a thin plot over two computed series; this pins those series to
# live R, exactly as seqedplot computes them (the renderer is structure-tested
# against these, not pixel baselines -- the RP7a rule).
#
#   type="survival" (Ph2): per non-ignored event, the Kaplan-Meier survival curve
#     of the time to the event's FIRST occurrence (seqeage column 1); a sequence
#     in which the event never occurs is right-censored at seqelength. R call:
#     survfit(Surv(time, status) ~ 1) -- UNWEIGHTED (seqedplot.R:201-235). The
#     weighted fixture proves the survival branch ignores seqeweight.
#   type="hazard"  (Ph3): added in the hazard commit.
#
# Requires live R + TraMineR 2.2-14 / TraMineRextras 0.6.9 / survival 3.8.6;
# see memory "live-traminer-r-env".
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
suppressMessages(library(TraMineRextras))
suppressMessages(library(survival))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
stopifnot(as.character(packageVersion("TraMineRextras")) == "0.6.9")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

seqeage <- TraMineR:::seqeage

# --- survival fixture: events A (code 1) and B (code 2), with sequences in which
# one event never occurs (-> right-censored at the observation length). -------
tse <- data.frame(
  id    = c("s1", "s1", "s2", "s3", "s3", "s4", "s4", "s5"),
  time  = c(1,    3,    2,    1,    4,    1,    3,    2),
  event = c("A",  "B",  "A",  "B",  "A",  "A",  "A",  "B"),
  stringsAsFactors = FALSE
)
obs_len <- c(5, 5, 6, 4, 5)   # seqelength per sequence (s1..s5)
wts <- c(1.5, 2, 0.5, 3, 1)   # weights -- survival must IGNORE these

dump_mat <- function(m, tag) {
  m <- as.matrix(m)
  con <- file(file.path(out_dir, paste0("seqedplot_", tag, ".csv")), "w")
  cat(sprintf("# %d x %d : %s\n", nrow(m), ncol(m), tag), file = con)
  for (i in seq_len(nrow(m))) {
    cat(paste(sprintf("%.15g", m[i, ]), collapse = ","), "\n", file = con, sep = "")
  }
  close(con)
  cat("== seqedplot_", tag, " ==\n", sep = "")
  print(m)
}

# Replicate seqedplot's per-event survival series exactly (seqedplot.R:204-209,
# no ages clipping -> ages=NULL default): first-occurrence age, censor=seqelength.
surv_series <- function(e, ev) {
  am <- seqeage(e, as.integer(ev))
  am[am == -1] <- NA
  time <- am[, 1]
  endobs <- seqelength(e)
  status <- !is.na(time)
  time[!status] <- endobs[!status]
  fit <- survfit(Surv(time, as.integer(status)) ~ 1)
  s <- summary(fit)
  cbind(s$time, s$n.risk, s$n.event, s$surv, s$std.err)
}

e <- seqecreate(tse)
seqelength(e) <- obs_len
cat("levels:", paste(levels(e), collapse = ","), "\n")

dump_mat(surv_series(e, 1), "surv_A")  # time,n.risk,n.event,surv,std.err
dump_mat(surv_series(e, 2), "surv_B")

# Weighted fixture: identical events + lengths, but with weights set. The
# survival series must be IDENTICAL (survfit is called unweighted).
ew <- seqecreate(tse)
seqelength(ew) <- obs_len
seqeweight(ew) <- wts
dump_mat(surv_series(ew, 1), "surv_A_wtd")  # == surv_A
dump_mat(surv_series(ew, 2), "surv_B_wtd")  # == surv_B

cat("SEQEDPLOT_SURVIVAL_GOLDENS_DONE\n")

# ---------------------------------------------------------------------------
# Hazard series (Ph3): the weighted mean number of occurrences per equal-width
# age period. The delicate part is reproducing R's cut(ages, breaks=N) interval
# boundaries (range extended by dx/1000 at each end, then formatC dig.lab
# rounding), which become the bin breakpoints. This replicates seqedplot.R
# (:119-172) verbatim, including the regex parse of the cut() labels.
# ---------------------------------------------------------------------------
htse <- data.frame(
  id    = c("s1", "s1", "s1", "s2", "s2", "s3", "s3", "s4", "s4", "s4", "s5", "s5", "s6", "s6"),
  time  = c(1,    3,    5,    2,    4,    1,    6,    1,    2,    4,    3,    5,    5,    2),
  event = c("A",  "A",  "B",  "A",  "B",  "B",  "A",  "A",  "A",  "A",  "B",  "B",  "A",  "B"),
  stringsAsFactors = FALSE
)
hgroup <- c("X", "X", "X", "Y", "Y", "Y")          # s1..s6
hwts <- c(1.0, 2.0, 0.5, 1.5, 2.5, 1.0)
he <- seqecreate(htse)
seqeweight(he) <- hwts
cat("hazard levels:", paste(levels(he), collapse = ","), "\n")

# seqedplot.R:119-172 replicated: global ages range -> cutpoints -> per-(group,
# event) xtabs(rep(w, ncol) ~ cut(...))/group size. Returns onx + the ony rows.
hazard_all <- function(e, group, breaks, ignore = NULL) {
  allevents <- which(!levels(e) %in% ignore)
  nevent <- length(allevents)
  agesmatrices <- list()
  minage <- NA
  maxage <- NA
  for (event in 1:nevent) {
    am <- seqeage(e, as.integer(allevents[event]))
    am[am == -1] <- NA
    suppressWarnings(minage <- min(c(minage, am), na.rm = TRUE))
    suppressWarnings(maxage <- max(c(maxage, am), na.rm = TRUE))
    agesmatrices[[event]] <- am
  }
  ages <- c(minage, maxage)
  labs <- levels(cut(ages, breaks))
  cutpoints <- cbind(
    lower = as.numeric(sub("\\((.+),.*", "\\1", labs)),
    upper = as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", labs))
  )
  cutpoints <- c(cutpoints[, 1], max(cutpoints[, 2]))
  onx <- cutpoints[-1]
  glev <- levels(factor(group))
  gindex <- lapply(glev, function(L) which(group == L))
  rows <- list(onx)
  names_rows <- "onx"
  for (np in seq_along(gindex)) {
    for (event in 1:nevent) {
      ony <- as.numeric(xtabs(rep(seqeweight(e)[gindex[[np]]], ncol(agesmatrices[[event]])) ~
        cut(agesmatrices[[event]][gindex[[np]], ],
          breaks = cutpoints, ordered_result = TRUE
        )) / length(gindex[[np]]))
      rows[[length(rows) + 1]] <- ony
      names_rows <- c(names_rows, paste0("g", np, "_", levels(e)[allevents[event]]))
    }
  }
  m <- do.call(rbind, rows)
  rownames(m) <- names_rows
  m
}

# Grouped (X/Y), 2 break counts.
dump_mat(hazard_all(he, hgroup, 5), "haz_b5")  # rows: onx, g1_A, g1_B, g2_A, g2_B
dump_mat(hazard_all(he, hgroup, 3), "haz_b3")  # rows: onx, g1_A, g1_B, g2_A, g2_B
# Ungrouped (one group = all sequences), breaks=5.
dump_mat(hazard_all(he, rep("all", 6), 5), "haz_all_b5")  # rows: onx, g1_A, g1_B

cat("SEQEDPLOT_HAZARD_GOLDENS_DONE\n")

# ---------------------------------------------------------------------------
# Corner cases (Ph6 mini-audit): ages clipping for both series. seqedplot clips
# the age matrix to [ages[1], ages[2]] when `ages` is given (seqedplot.R:131-134):
# out-of-window first occurrences become censored (survival) / dropped (hazard),
# and the hazard cut-points come from the provided `ages`, not the data range.
# ---------------------------------------------------------------------------
surv_series_ages <- function(e, ev, ages) {
  am <- seqeage(e, as.integer(ev))
  am[am == -1] <- NA
  am[am < ages[1]] <- NA
  am[am > ages[2]] <- NA
  time <- am[, 1]
  endobs <- seqelength(e)
  status <- !is.na(time)
  time[!status] <- endobs[!status]
  s <- summary(survfit(Surv(time, as.integer(status)) ~ 1))
  cbind(s$time, s$n.risk, s$n.event, s$surv, s$std.err)
}

# Survival with ages=c(1.5, 4.5): clips first-occurrence ages of 1 (-> censored).
dump_mat(surv_series_ages(e, 1, c(1.5, 4.5)), "surv_clip_A")
dump_mat(surv_series_ages(e, 2, c(1.5, 4.5)), "surv_clip_B")

hazard_ages <- function(e, breaks, ages) {
  allevents <- seq_along(levels(e))
  nevent <- length(allevents)
  agesmatrices <- list()
  for (event in 1:nevent) {
    am <- seqeage(e, as.integer(allevents[event]))
    am[am == -1] <- NA
    am[am < ages[1]] <- NA
    am[am > ages[2]] <- NA
    agesmatrices[[event]] <- am
  }
  labs <- levels(cut(ages, breaks))
  cutpoints <- cbind(
    lower = as.numeric(sub("\\((.+),.*", "\\1", labs)),
    upper = as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", labs))
  )
  cutpoints <- c(cutpoints[, 1], max(cutpoints[, 2]))
  onx <- cutpoints[-1]
  rows <- list(onx)
  names_rows <- "onx"
  for (event in 1:nevent) {
    ony <- as.numeric(xtabs(rep(seqeweight(e), ncol(agesmatrices[[event]])) ~
      cut(agesmatrices[[event]], breaks = cutpoints, ordered_result = TRUE)) / length(e))
    rows[[length(rows) + 1]] <- ony
    names_rows <- c(names_rows, levels(e)[allevents[event]])
  }
  m <- do.call(rbind, rows)
  rownames(m) <- names_rows
  m
}

# Hazard with ages=c(1.5, 5.5): drops the occurrences at ages 1 and 6.
dump_mat(hazard_ages(he, 5, c(1.5, 5.5)), "haz_clip_b5")  # rows: onx, A, B

cat("SEQEDPLOT_CORNER_GOLDENS_DONE\n")
