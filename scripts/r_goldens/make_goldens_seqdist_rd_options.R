# Golden generation: every seqdist option R documents, at a NON-DEFAULT value.
#
# Cycle-7 W-05. traminex's seqdist takes its per-method options through
# **kwargs, so inspect.signature cannot see them and the cycle-6 argument-name
# sweep skipped them entirely. This generator walks man/seqdist.Rd option by
# option, and for every (method, option) pair the Rd documents it records R's
# answer with that option set away from its default -- so a Python side that
# accepted the name but ignored the value would be caught.
#
# The per-method applicability below is taken from the \arguments section of
# SOURCES/TraMineR/man/seqdist.Rd:
#
#   norm          OM OMloc OMslen OMspell OMstran TWED HAM DHD LCS LCP RLCP
#                 CHI2 EUCLID                                    (Rd:48-63)
#   indel         OM OMslen OMspell OMstran                      (Rd:65-81)
#   sm            OM OMloc OMslen OMspell OMstran HAM TWED (matrix),
#                 DHD (array)                                    (Rd:83-115)
#   kweights      NMS NMSMST SVRspell                            (Rd:131-141)
#   tpow          OMspell NMSMST SVRspell                        (Rd:143-148)
#   expcost       OMloc OMspell                                  (Rd:150-156)
#   context       OMloc                                          (Rd:158-162)
#   link          OMslen                                         (Rd:164-170)
#   h             OMslen TWED                                    (Rd:172-181)
#   nu            TWED                                           (Rd:183-188)
#   transindel    OMstran                                        (Rd:190-196)
#   otto          OMstran                                        (Rd:198-201)
#   previous      OMstran                                        (Rd:203-208)
#   add.column    OMstran                                        (Rd:210-216)
#   breaks        CHI2 EUCLID                                    (Rd:218-222)
#   step          CHI2 EUCLID                                    (Rd:224-230)
#   overlap       CHI2 EUCLID                                    (Rd:232-238)
#   weighted      CHI2, and whenever sm is a cost-method string  (Rd:240-246)
#   global.pdotj  CHI2 only                                      (Rd:248-253)
#   prox          NMS SVRspell                                   (Rd:255-260)
#
# Values are written with %.15g. Requires live R + TraMineR 2.2-14.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(packageVersion("TraMineR") == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## Corpus ------------------------------------------------------------------ ##
set.seed(20260830)
n_seq <- 10L
seq_len_ <- 8L
alphabet <- c("A", "B", "C", "D")
codes <- matrix(sample(seq_along(alphabet), n_seq * seq_len_, replace = TRUE),
                nrow = n_seq, ncol = seq_len_)
states <- matrix(alphabet[codes], nrow = n_seq, ncol = seq_len_)
weights <- round(seq(0.5, 2.3, length.out = n_seq), 3)
sdata <- suppressMessages(
  seqdef(as.data.frame(states, stringsAsFactors = FALSE),
         alphabet = alphabet, weights = weights)
)
write.csv(as.data.frame(codes - 1L),
          file.path(out_dir, "seqdist_rd_codes.csv"), row.names = FALSE)
writeLines(alphabet, file.path(out_dir, "seqdist_rd_alphabet.txt"))
writeLines(sprintf("%.15g", weights), file.path(out_dir, "seqdist_rd_weights.txt"))

K <- length(alphabet)
# A fixed, asymmetric-free substitution matrix so `sm` is unambiguous.
sm_mat <- matrix(0, K, K)
for (i in 1:K) for (j in 1:K) if (i != j) sm_mat[i, j] <- 1 + abs(i - j) / 2
write.csv(as.data.frame(sm_mat), file.path(out_dir, "seqdist_rd_subm.csv"),
          row.names = FALSE)
# A DHD cost ARRAY: one K x K matrix per POSITION, so the third dimension is
# ncol(seqdata) -- R rejects K x K x (L-1) with "size of substitution cost
# matrix must be 4x4x8".
sm_arr <- array(0, dim = c(K, K, seq_len_))
for (t in seq_len(seq_len_)) sm_arr[, , t] <- sm_mat * (1 + t / 10)
write.csv(as.data.frame(matrix(as.vector(sm_arr), nrow = K * K)),
          file.path(out_dir, "seqdist_rd_subm_array.csv"), row.names = FALSE)
prox_mat <- diag(K)
for (i in 1:K) for (j in 1:K) if (i != j) prox_mat[i, j] <- 1 / (1 + abs(i - j))
write.csv(as.data.frame(prox_mat), file.path(out_dir, "seqdist_rd_prox.csv"),
          row.names = FALSE)
kweights_vec <- c(1.0, 0.5, 0.25, 0.125, 0.1, 0.1, 0.1, 0.1)
writeLines(sprintf("%.15g", kweights_vec),
           file.path(out_dir, "seqdist_rd_kweights.txt"))

## The (method, option, non-default value) cases --------------------------- ##
# Each entry: label, method, and the extra arguments beyond seqdata/method.
cases <- list()
add <- function(label, method, args) {
  cases[[label]] <<- list(method = method, args = args)
}

# --- norm, for each of the 13 methods the Rd lists ------------------------ #
add("norm|OM",      "OM",      list(sm = sm_mat, indel = 1.0, norm = "maxlength"))
add("norm|OMloc",   "OMloc",   list(sm = sm_mat, expcost = 0.5, norm = "YujianBo"))
add("norm|OMslen",  "OMslen",  list(sm = sm_mat, indel = 1.0, norm = "YujianBo"))
add("norm|OMspell", "OMspell", list(sm = sm_mat, indel = 1.0, norm = "YujianBo"))
add("norm|OMstran", "OMstran", list(sm = sm_mat, indel = 1.0, otto = 0.5, norm = "YujianBo"))
add("norm|TWED",    "TWED",    list(sm = sm_mat, nu = 0.5, h = 0.5, norm = "YujianBo"))
add("norm|HAM",     "HAM",     list(sm = sm_mat, norm = "maxlength"))
add("norm|DHD",     "DHD",     list(norm = "maxlength"))
add("norm|LCS",     "LCS",     list(norm = "gmean"))
add("norm|LCP",     "LCP",     list(norm = "gmean"))
add("norm|RLCP",    "RLCP",    list(norm = "gmean"))
add("norm|CHI2",    "CHI2",    list(step = 1, norm = "auto"))
add("norm|EUCLID",  "EUCLID",  list(step = 1, norm = "auto"))

# --- indel ---------------------------------------------------------------- #
add("indel|OM",      "OM",      list(sm = sm_mat, indel = 1.75))
add("indel|OMslen",  "OMslen",  list(sm = sm_mat, indel = 1.75))
add("indel|OMspell", "OMspell", list(sm = sm_mat, indel = 1.75))
add("indel|OMstran", "OMstran", list(sm = sm_mat, indel = 1.75, otto = 0.5))
# vector-valued indel (one per state) -- the Rd's second form
add("indel_vec|OM",  "OM",      list(sm = sm_mat, indel = c(0.5, 1.0, 1.5, 2.0)))

# --- sm ------------------------------------------------------------------- #
add("sm|OM",      "OM",      list(sm = sm_mat, indel = 1.0))
add("sm|OMloc",   "OMloc",   list(sm = sm_mat, expcost = 0.5))
add("sm|OMslen",  "OMslen",  list(sm = sm_mat, indel = 1.0))
add("sm|OMspell", "OMspell", list(sm = sm_mat, indel = 1.0))
add("sm|OMstran", "OMstran", list(sm = sm_mat, indel = 1.0, otto = 0.5))
add("sm|HAM",     "HAM",     list(sm = sm_mat))
add("sm|TWED",    "TWED",    list(sm = sm_mat, nu = 0.5, h = 0.5))
add("sm_array|DHD", "DHD",   list(sm = sm_arr))
# sm as a cost-method string
add("sm_str_TRATE|OM",      "OM", list(sm = "TRATE", indel = 1.0))
add("sm_str_CONSTANT|OM",   "OM", list(sm = "CONSTANT", indel = 1.0))
add("sm_str_INDELS|OM",     "OM", list(sm = "INDELS", indel = 1.0))
add("sm_str_INDELSLOG|OM",  "OM", list(sm = "INDELSLOG", indel = 1.0))

# --- kweights ------------------------------------------------------------- #
add("kweights|NMS",      "NMS",      list(kweights = kweights_vec, prox = prox_mat))
add("kweights|NMSMST",   "NMSMST",   list(kweights = kweights_vec))
add("kweights|SVRspell", "SVRspell", list(kweights = kweights_vec, prox = prox_mat))

# --- tpow ----------------------------------------------------------------- #
add("tpow|OMspell",  "OMspell",  list(sm = sm_mat, indel = 1.0, tpow = 1.5))
add("tpow|NMSMST",   "NMSMST",   list(tpow = 1.5))
add("tpow|SVRspell", "SVRspell", list(prox = prox_mat, tpow = 1.5))

# --- expcost / context ---------------------------------------------------- #
add("expcost|OMloc",   "OMloc",   list(sm = sm_mat, expcost = 0.3))
add("expcost|OMspell", "OMspell", list(sm = sm_mat, indel = 1.0, expcost = 0.3))
add("context|OMloc",   "OMloc",   list(sm = sm_mat, expcost = 0.3, context = 0.2))

# --- link / h / nu -------------------------------------------------------- #
add("link|OMslen", "OMslen", list(sm = sm_mat, indel = 1.0, link = "gmean"))
add("h|OMslen",    "OMslen", list(sm = sm_mat, indel = 1.0, h = 0.75))
add("h|TWED",      "TWED",   list(sm = sm_mat, nu = 0.5, h = 0.75))
add("nu|TWED",     "TWED",   list(sm = sm_mat, nu = 0.25, h = 0.5))

# --- OMstran-only options ------------------------------------------------- #
add("transindel|OMstran", "OMstran",
    list(sm = sm_mat, indel = 1.0, otto = 0.5, transindel = "subcost"))
add("transindel_prob|OMstran", "OMstran",
    list(sm = sm_mat, indel = 1.0, otto = 0.5, transindel = "prob"))
add("otto|OMstran",       "OMstran", list(sm = sm_mat, indel = 1.0, otto = 0.25))
add("previous|OMstran",   "OMstran",
    list(sm = sm_mat, indel = 1.0, otto = 0.5, previous = TRUE))
add("add_column|OMstran", "OMstran",
    list(sm = sm_mat, indel = 1.0, otto = 0.5, add.column = FALSE))

# --- CHI2 / EUCLID interval options --------------------------------------- #
add("step|CHI2",     "CHI2",   list(step = 2))
add("step|EUCLID",   "EUCLID", list(step = 2))
add("overlap|CHI2",  "CHI2",   list(step = 2, overlap = TRUE))
add("overlap|EUCLID","EUCLID", list(step = 2, overlap = TRUE))
add("breaks|CHI2",   "CHI2",   list(breaks = list(c(1, 3), c(4, 8))))
add("breaks|EUCLID", "EUCLID", list(breaks = list(c(1, 3), c(4, 8))))
add("breaks_overlapping|CHI2", "CHI2",
    list(breaks = list(c(1, 4), c(3, 6), c(5, 8))))
add("global_pdotj_obs|CHI2", "CHI2", list(step = 2, global.pdotj = "obs"))
add("global_pdotj_vec|CHI2", "CHI2",
    list(step = 2, global.pdotj = c(0.4, 0.3, 0.2, 0.1)))
# R renormalises a supplied vector by its sum (seqdist-CHI2.R:64), so this
# unnormalised vector must give exactly the same answer as the one above.
add("global_pdotj_unnormalised|CHI2", "CHI2",
    list(step = 2, global.pdotj = c(4, 3, 2, 1)))
# breaks wins over step/overlap, with a warning (seqdist.R:256-258).
add("breaks_beats_step|CHI2", "CHI2",
    list(breaks = list(c(1, 3), c(4, 8)), step = 5, overlap = TRUE))
add("weighted_false|CHI2", "CHI2", list(step = 2, weighted = FALSE))
add("weighted_false|OM_TRATE", "OM", list(sm = "TRATE", indel = 1.0, weighted = FALSE))

# --- prox ----------------------------------------------------------------- #
add("prox|NMS",      "NMS",      list(prox = prox_mat))
add("prox|SVRspell", "SVRspell", list(prox = prox_mat))

## Run --------------------------------------------------------------------- ##
rows <- character(0)
errors <- character(0)
for (label in names(cases)) {
  entry <- cases[[label]]
  value <- tryCatch(
    {
      d <- suppressMessages(do.call(
        seqdist, c(list(seqdata = sdata, method = entry$method), entry$args)
      ))
      paste(sprintf("%.15g", as.vector(as.matrix(d))), collapse = ";")
    },
    error = function(e) paste0("ERROR: ", conditionMessage(e))
  )
  if (startsWith(value, "ERROR:")) {
    errors <- c(errors, sprintf("%s,%s", label, value))
  } else {
    rows <- c(rows, sprintf("%s,%s,%s", label, entry$method, value))
  }
}

writeLines(c("# label,method,column-major distance matrix", rows),
           file.path(out_dir, "seqdist_rd_option_values.csv"))
writeLines(c("# label,error", errors),
           file.path(out_dir, "seqdist_rd_option_errors.csv"))

cat("SEQDIST_RD_OPTIONS_GOLDENS_DONE\n")
cat("cases:", length(cases), " values:", length(rows), " errors:", length(errors), "\n")
for (e in errors) cat("  ", e, "\n")
