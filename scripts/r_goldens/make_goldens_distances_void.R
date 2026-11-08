# Complete R1 seqdist contract matrix from TraMineR 2.2-14.
#
# The manifest deliberately contains both successful numeric results and R
# refusals.  Numeric values are written in logical row-major order with
# formatC(..., format="g", digits=15); errors are reduced to stable categories
# so the Python parity test can compare supported/refused dispositions without
# depending on TraMineR's message prefix formatting.
#
# Run from the repository root:
#   /usr/bin/Rscript scripts/r_goldens/make_goldens_distances_void.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))

if (!identical(as.character(packageVersion("TraMineR")), "2.2.14")) {
  stop("distance R1 goldens require TraMineR 2.2.14")
}

fixture_env <- new.env(parent = globalenv())
suppressMessages(sys.source(
  file.path("scripts", "r_goldens", "make_goldens_corestats.R"),
  envir = fixture_env
))
fx <- suppressMessages(fixture_env$corestats_fixtures())

# A ragged fixture with an unambiguous most-frequent sequence (row 2, repeated
# as row 5).  It also retains the interior missing state in row 3.
freq_raw <- rbind(
  c("A", "A", "B", "C", "C", "C"),
  c("A", "B", "B", "*", "*", "*"),
  c("A", "*", "A", "B", "*", "*"),
  c("C", "C", "*", "*", "*", "*"),
  c("A", "B", "B", "*", "*", "*")
)
freq_seq <- suppressMessages(seqdef(
  freq_raw,
  alphabet = c("A", "B", "C"),
  missing = "*"
))

# A separately constructed, shorter reference sequence.  The first "*" is a
# retained missing state; the trailing "*" values are structural void.
external_ref <- suppressMessages(seqdef(
  matrix(c("A", "*", "C", "C", "*", "*"), nrow = 1),
  alphabet = c("A", "B", "C"),
  missing = "*"
))

# A complete source plus a separately observed interior-missing reference.
# Their logical alphabets agree even though only the reference materializes nr.
asym_seq <- suppressMessages(seqdef(
  rbind(c("A", "A", "A"), c("A", "B", "A")),
  alphabet = c("A", "B"),
  missing = "*"
))
asym_ref <- suppressMessages(seqdef(
  matrix(c("A", "*", "A"), nrow = 1),
  alphabet = c("A", "B"),
  missing = "*"
))

# Equally frequent sequences deliberately start out of lexical numeric-SPS
# order, so refseq=0 must select AA rather than the first row BB.
tie_seq <- suppressMessages(seqdef(
  rbind(c("B", "B"), c("A", "A"), c("C", "C")),
  alphabet = c("A", "B", "C")
))

# Empty/short rows exercise edit-distance borders and transition tokenization.
empty_seq <- suppressMessages(seqdef(
  rbind(c("*", "*", "*"), c("A", "B", "A"), c("A", "B", "B")),
  alphabet = c("A", "B"),
  missing = "*"
))

asymmetric_seq <- suppressMessages(seqdef(
  rbind("A", "B"),
  alphabet = c("A", "B")
))
asymmetric_sm <- matrix(
  c(0, 10, 1, 0),
  nrow = 2,
  dimnames = list(c("A", "B"), c("A", "B"))
)
symmetric4_sm <- matrix(
  c(0, 4, 4, 0),
  nrow = 2,
  dimnames = list(c("A", "B"), c("A", "B"))
)
missing_indel_seq <- suppressMessages(seqdef(
  rbind(c("A", "*"), c("*", "A")),
  alphabet = "A",
  missing = "*",
  right = NA
))
empty_frequency_seq <- suppressMessages(seqdef(
  rbind(
    c("*", "*", "*"),
    c("*", "*", "*"),
    c("A", "*", "A"),
    c("A", "B", "*")
  ),
  alphabet = c("A", "B"),
  missing = "*"
))
vector_indel_seq <- suppressMessages(seqdef(
  rbind(c("A", "A", "B"), c("A", "B", "B")),
  alphabet = c("A", "B")
))
declared_state_dhd_seq <- suppressMessages(seqdef(
  rbind(c("A", "A", "A"), c("A", "B", "A")),
  alphabet = c("A", "B", "C")
))
missing_ragged_dhd_seq <- suppressMessages(seqdef(
  rbind(c("A", "*", "A", "*"), c("B", "A", "B", "*")),
  alphabet = c("A", "B"),
  missing = "*"
))
short_ref <- suppressMessages(seqdef(
  matrix(c("A", "B"), nrow = 1),
  alphabet = c("A", "B"),
  missing = "*"
))

fixtures <- list(
  B = fx$B, C = fx$C, F = freq_seq, G = asym_seq,
  H = tie_seq, J = asymmetric_seq, K = empty_seq, L = missing_indel_seq,
  M = empty_frequency_seq, N = vector_indel_seq,
  O = declared_state_dhd_seq, P = missing_ragged_dhd_seq,
  W = fx$Bw, CW = fx$Cw
)
sms <- list(
  B = suppressMessages(seqsubm(
    fx$B, "CONSTANT", with.missing = TRUE, cval = 2, miss.cost = 2
  )),
  C = suppressMessages(seqsubm(
    fx$C, "CONSTANT", with.missing = TRUE, cval = 2, miss.cost = 2
  ))
)
empty_sm <- suppressMessages(seqsubm(empty_seq, "CONSTANT", cval = 2))

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

methods <- c(
  "OM", "OMloc", "OMslen", "OMspell", "OMstran", "HAM", "DHD",
  "CHI2", "EUCLID", "LCS", "LCP", "RLCP", "NMS", "NMSMST",
  "SVRspell", "TWED"
)
supported_auto <- setdiff(methods, c("NMS", "NMSMST", "SVRspell"))
all_norms <- c("none", "auto", "maxlength", "gmean", "maxdist", "YujianBo")

format_value <- function(value) {
  if (is.nan(value)) return("NaN")
  if (is.na(value)) return("NA")
  if (is.infinite(value)) return(if (value > 0) "Inf" else "-Inf")
  trimws(formatC(value, format = "g", digits = 15))
}

csv_quote <- function(value) {
  paste0('"', gsub('"', '""', as.character(value), fixed = TRUE), '"')
}

case_specs <- list()

add_case <- function(
  case, method, fixture, options, python_call, args,
  allow_warnings = FALSE
) {
  if (!is.null(case_specs[[case]])) stop("duplicate distance case: ", case)
  case_specs[[case]] <<- list(
    case = case,
    method = method,
    fixture = fixture,
    options = options,
    python_call = python_call,
    args = args,
    allow_warnings = allow_warnings
  )
}

method_args <- function(fixture, method, norm = "none") {
  x <- fixtures[[fixture]]
  args <- list(
    seqdata = x,
    method = method,
    with.missing = TRUE,
    norm = norm,
    full.matrix = TRUE
  )
  if (method == "OM") {
    args <- c(args, list(sm = sms[[fixture]], indel = 1))
  } else if (method == "OMloc") {
    args <- c(args, list(sm = sms[[fixture]], expcost = 0.25, context = 0.5))
  } else if (method == "OMslen") {
    args <- c(args, list(sm = sms[[fixture]], indel = 1, link = "mean", h = 0.5))
  } else if (method == "OMspell") {
    args <- c(args, list(sm = sms[[fixture]], indel = 1, expcost = 0.5, tpow = 1))
  } else if (method == "OMstran") {
    args <- c(args, list(
      sm = sms[[fixture]], indel = 1, otto = 0.5,
      transindel = "constant", previous = FALSE, add.column = TRUE
    ))
  } else if (method == "TWED") {
    args <- c(args, list(sm = sms[[fixture]], nu = 0.5, h = 0.5))
  }
  args
}

method_option_suffix <- function(method) {
  switch(
    method,
    OM = ";sm=constant2;indel=1",
    OMloc = ";sm=constant2;expcost=0.25;context=0.5",
    OMslen = ";sm=constant2;indel=1;link=mean;h=0.5",
    OMspell = ";sm=constant2;indel=1;expcost=0.5;tpow=1",
    OMstran = paste0(
      ";sm=constant2;indel=1;otto=0.5;transindel=constant;",
      "previous=FALSE;add.column=TRUE"
    ),
    TWED = ";sm=constant2;nu=0.5;h=0.5",
    ""
  )
}

python_kw_suffix <- function(method) {
  switch(
    method,
    OM = ", subm=2.0, indel=1.0",
    OMloc = ", subm=2.0, expcost=0.25, context=0.5",
    OMslen = ", subm=2.0, indel=1.0, link=\"mean\", h=0.5",
    OMspell = ", subm=2.0, indel=1.0, expcost=0.5, tpow=1.0",
    OMstran = paste0(
      ", subm=2.0, indel=1.0, otto=0.5, transindel=\"constant\"",
      ", previous=False, add_column=True"
    ),
    TWED = ", subm=2.0, nu=0.5, h=0.5",
    ""
  )
}

add_standard_case <- function(case, fixture, method, norm) {
  fixture_call <- if (fixture == "B") "fixture_b()" else "fixture_c()"
  options <- paste0(
    "method=", method, ";with.missing=TRUE;norm=", norm,
    ";full.matrix=TRUE", method_option_suffix(method)
  )
  python_call <- paste0(
    "seqdist(", fixture_call, ", method=\"", method,
    "\", with_missing=True, norm=\"", norm, "\"",
    python_kw_suffix(method), ")"
  )
  add_case(
    case, method, fixture, options, python_call,
    method_args(fixture, method, norm)
  )
}

# Required matrix: every official method on one interior-missing fixture and
# one ragged fixture.  HAM and DHD intentionally become R refusal rows on C.
for (method in methods) {
  slug <- tolower(method)
  add_standard_case(paste0("interior_", slug), "B", method, "none")
  add_standard_case(paste0("ragged_", slug), "C", method, "none")
}

# R's auto selection for every method where normalization is implemented.
# HAM/DHD use the equal-length B fixture; all other methods exercise ragged C.
for (method in supported_auto) {
  fixture <- if (method %in% c("HAM", "DHD")) "B" else "C"
  add_standard_case(
    paste0("auto_", tolower(method), "_", if (fixture == "B") "interior" else "ragged"),
    fixture,
    method,
    "auto"
  )
}

# All remaining normalization formula names on the same ragged OM inputs.
for (norm in setdiff(all_norms, c("none", "auto"))) {
  add_standard_case(paste0("norm_om_ragged_", tolower(norm)), "C", "OM", norm)
}

# Exhaust the unsupported normalization boundary: five non-none norms for the
# three subsequence-vector methods, and the four non-auto formula names for
# CHI2/EUCLID.
for (method in c("NMS", "NMSMST", "SVRspell")) {
  for (norm in setdiff(all_norms, "none")) {
    add_standard_case(
      paste0("refuse_norm_", tolower(method), "_", tolower(norm)),
      "B",
      method,
      norm
    )
  }
}
for (method in c("CHI2", "EUCLID")) {
  for (norm in c("maxlength", "gmean", "maxdist", "YujianBo")) {
    add_standard_case(
      paste0("refuse_norm_", tolower(method), "_", tolower(norm)),
      "B",
      method,
      norm
    )
  }
}

# Output and reference contracts.  Two Python interfaces (R-compatible refseq
# and existing zero-based ref_sets) intentionally share one R oracle result.
add_case(
  "ragged_lcs_full_matrix_false",
  "LCS",
  "C",
  "method=LCS;with.missing=TRUE;norm=none;full.matrix=FALSE",
  "seqdist(fixture_c(), method=\"LCS\", with_missing=True, norm=\"none\", full_matrix=False)",
  list(
    seqdata = fx$C, method = "LCS", with.missing = TRUE,
    norm = "none", full.matrix = FALSE
  )
)
add_case(
  "ragged_lcs_condensed_alias",
  "LCS",
  "C",
  "method=LCS;with.missing=TRUE;norm=none;full.matrix=FALSE",
  "seqdist(fixture_c(), method=\"LCS\", with_missing=True, norm=\"none\", condensed=True)",
  list(
    seqdata = fx$C, method = "LCS", with.missing = TRUE,
    norm = "none", full.matrix = FALSE
  )
)
for (case in c("ragged_lcs_refseq_sets", "ragged_lcs_ref_sets")) {
  python_ref <- if (case == "ragged_lcs_refseq_sets") {
    "refseq=[[1, 3], [2, 4]]"
  } else {
    "ref_sets=([0, 2], [1, 3])"
  }
  add_case(
    case,
    "LCS",
    "C",
    "method=LCS;with.missing=TRUE;norm=none;refseq=list(1|3,2|4)",
    paste0(
      "seqdist(fixture_c(), method=\"LCS\", with_missing=True, norm=\"none\", ",
      python_ref, ")"
    ),
    list(
      seqdata = fx$C, method = "LCS", with.missing = TRUE,
      norm = "none", refseq = list(c(1L, 3L), c(2L, 4L))
    )
  )
}
add_case(
  "ragged_om_refseq_single",
  "OM",
  "C",
  paste0(
    "method=OM;with.missing=TRUE;norm=none;refseq=2;",
    "sm=constant2;indel=1"
  ),
  paste0(
    "seqdist(fixture_c(), method=\"OM\", with_missing=True, norm=\"none\", ",
    "subm=2.0, indel=1.0, refseq=2)"
  ),
  list(
    seqdata = fx$C, method = "OM", with.missing = TRUE, norm = "none",
    sm = sms$C, indel = 1, refseq = 2L
  )
)
add_case(
  "ragged_lcp_refseq_most_frequent",
  "LCP",
  "F",
  "method=LCP;with.missing=TRUE;norm=none;refseq=0",
  "seqdist(fixture_f(), method=\"LCP\", with_missing=True, norm=\"none\", refseq=0)",
  list(
    seqdata = freq_seq, method = "LCP", with.missing = TRUE,
    norm = "none", refseq = 0L
  )
)
add_case(
  "ragged_lcs_refseq_external",
  "LCS",
  "C",
  "method=LCS;with.missing=TRUE;norm=none;refseq=external_missing_ragged",
  paste0(
    "seqdist(fixture_c(), method=\"LCS\", with_missing=True, norm=\"none\", ",
    "refseq=external_ref())"
  ),
  list(
    seqdata = fx$C, method = "LCS", with.missing = TRUE,
    norm = "none", refseq = external_ref
  )
)

# External references are a generic seqdist contract, not an LCS-only route.
# Cover every additional method that accepts the shorter reference (HAM/DHD
# retain their unequal-length refusal and are covered separately).
for (method in c(
  "OMloc", "OMslen", "OMspell", "OMstran", "CHI2", "EUCLID",
  "NMS", "NMSMST", "SVRspell", "TWED"
)) {
  slug <- tolower(method)
  args <- method_args("C", method, "none")
  args$refseq <- external_ref
  add_case(
    paste0("ragged_", slug, "_refseq_external"),
    method,
    "C",
    paste0(
      "method=", method,
      ";with.missing=TRUE;norm=none;full.matrix=TRUE;",
      "refseq=external_missing_ragged",
      method_option_suffix(method)
    ),
    paste0(
      "seqdist(fixture_c(), method=\"", method,
      "\", with_missing=True, norm=\"none\"",
      python_kw_suffix(method), ", refseq=external_ref())"
    ),
    args
  )
}

asym_args <- method_args("G", "NMS", "none")
asym_args$refseq <- asym_ref
add_case(
  "complete_nms_refseq_external_missing",
  "NMS",
  "G",
  paste0(
    "method=NMS;with.missing=TRUE;norm=none;full.matrix=TRUE;",
    "refseq=external_interior_missing"
  ),
  paste0(
    "seqdist(fixture_g(), method=\"NMS\", with_missing=True, ",
    "norm=\"none\", refseq=asymmetric_external_ref())"
  ),
  asym_args
)

soft_nms_args <- method_args("G", "NMS", "none")
soft_nms_args$prox <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, byrow = TRUE)
add_case(
  "complete_nms_nonidentity_prox",
  "NMS",
  "G",
  "method=NMS;with.missing=TRUE;norm=none;prox=nonidentity",
  paste0(
    "seqdist(fixture_g(), method=\"NMS\", with_missing=True, ",
    "norm=\"none\", prox=soft_prox())"
  ),
  soft_nms_args
)

# Source-review edge contracts discovered while closing the complete matrix.
# These freeze asymmetric OM bounds, TWED's empty border, TRATE preprocessing,
# attached weights, all OMstran token forms, and refseq=0 tie-breaking.
for (norm in c("maxdist", "gmean")) {
  add_case(
    paste0("asymmetric_om_", norm),
    "OM",
    "J",
    paste0("method=OM;norm=", norm, ";sm=asymmetric;indel=10"),
    paste0(
      "seqdist(fixture_j(), method=\"OM\", norm=\"", norm,
      "\", subm=asymmetric_subm(), indel=10.0)"
    ),
    list(
      seqdata = asymmetric_seq, method = "OM", norm = norm,
      sm = asymmetric_sm, indel = 10
    ),
    allow_warnings = TRUE
  )
}

for (method in c("OM", "OMslen", "OMspell", "OMstran")) {
  method_args_auto <- switch(
    method,
    OM = list(),
    OMslen = list(link = "mean", h = 0.5),
    OMspell = list(expcost = 0.5, tpow = 1),
    OMstran = list(
      otto = 0.5, transindel = "constant",
      previous = FALSE, add.column = TRUE
    )
  )
  python_suffix <- switch(
    method,
    OM = "",
    OMslen = ", link=\"mean\", h=0.5",
    OMspell = ", expcost=0.5, tpow=1.0",
    OMstran = paste0(
      ", otto=0.5, transindel=\"constant\", previous=False, ",
      "add_column=True"
    )
  )
  add_case(
    paste0("explicit_sm_auto_indel_", tolower(method)),
    method,
    "J",
    paste0("method=", method, ";norm=none;sm=constant4;indel=auto"),
    paste0(
      "seqdist(fixture_j(), method=\"", method,
      "\", norm=\"none\", subm=4.0", python_suffix, ")"
    ),
    c(
      list(
        seqdata = asymmetric_seq, method = method, norm = "none",
        sm = symmetric4_sm
      ),
      method_args_auto
    )
  )
}

add_case(
  "empty_twed_auto_indel_border",
  "TWED",
  "K",
  "method=TWED;norm=none;sm=constant2;nu=0.5;h=0.5",
  paste0(
    "seqdist(fixture_k(), method=\"TWED\", norm=\"none\", ",
    "subm=2.0, nu=0.5, h=0.5)"
  ),
  list(
    seqdata = empty_seq, method = "TWED", norm = "none",
    sm = empty_sm, nu = 0.5, h = 0.5
  )
)

trate_specs <- list(
  list(
    case = "ragged_om_trate", fixture = "C", seqdata = fx$C,
    suffix = "", refseq = NULL, weighted = NULL
  ),
  list(
    case = "ragged_om_trate_refseq_single", fixture = "C", seqdata = fx$C,
    suffix = ", refseq=2", refseq = 2L, weighted = NULL
  ),
  list(
    case = "ragged_om_trate_refseq_external", fixture = "C", seqdata = fx$C,
    suffix = ", refseq=external_ref()", refseq = external_ref, weighted = NULL
  ),
  list(
    case = "ragged_om_trate_weighted", fixture = "CW", seqdata = fx$Cw,
    suffix = "", refseq = NULL, weighted = NULL
  ),
  list(
    case = "ragged_om_trate_attached_unweighted", fixture = "CW", seqdata = fx$Cw,
    suffix = ", weighted=False", refseq = NULL, weighted = FALSE
  )
)
for (spec in trate_specs) {
  args <- list(
    seqdata = spec$seqdata, method = "OM", sm = "TRATE",
    with.missing = TRUE, norm = "none"
  )
  if (!is.null(spec$refseq)) args$refseq <- spec$refseq
  if (!is.null(spec$weighted)) args$weighted <- spec$weighted
  fixture_call <- if (spec$fixture == "CW") "fixture_cw()" else "fixture_c()"
  add_case(
    spec$case,
    "OM",
    spec$fixture,
    paste0(
      "method=OM;with.missing=TRUE;norm=none;sm=TRATE",
      if (is.null(spec$weighted)) "" else ";weighted=FALSE",
      if (is.null(spec$refseq)) "" else ";refseq"
    ),
    paste0(
      "seqdist(", fixture_call,
      ", method=\"OM\", subm=\"TRATE\", with_missing=True, norm=\"none\"",
      spec$suffix, ")"
    ),
    args
  )
}

add_case(
  "interior_om_trate_zero_pre_missing_indel",
  "OM",
  "L",
  "method=OM;with.missing=TRUE;norm=none;sm=TRATE;auto_indel_before_missing",
  paste0(
    "seqdist(fixture_l(), method=\"OM\", subm=\"TRATE\", ",
    "with_missing=True, norm=\"none\")"
  ),
  list(
    seqdata = missing_indel_seq, method = "OM", sm = "TRATE",
    with.missing = TRUE, norm = "none"
  )
)

for (method in c("DHD", "CHI2")) {
  for (use_weights in c(TRUE, FALSE)) {
    add_case(
      paste0(
        "interior_", tolower(method), "_attached_",
        if (use_weights) "weighted" else "unweighted"
      ),
      method,
      "W",
      paste0(
        "method=", method,
        ";with.missing=TRUE;norm=none;weighted=", use_weights
      ),
      paste0(
        "seqdist(fixture_bw(), method=\"", method,
        "\", with_missing=True, norm=\"none\", weighted=",
        if (use_weights) "True" else "False", ")"
      ),
      list(
        seqdata = fx$Bw, method = method, with.missing = TRUE,
        norm = "none", weighted = use_weights
      )
    )
  }
}

for (transindel in c("constant", "prob", "subcost")) {
  for (previous in c(FALSE, TRUE)) {
    for (add_column in c(FALSE, TRUE)) {
      if (transindel == "constant" && !previous && add_column) next
      add_case(
        paste0(
          "ragged_omstran_", transindel,
          "_previous_", tolower(as.character(previous)),
          "_add_column_", tolower(as.character(add_column))
        ),
        "OMstran",
        "C",
        paste0(
          "method=OMstran;with.missing=TRUE;norm=none;sm=constant2;",
          "indel=1;otto=0.5;transindel=", transindel,
          ";previous=", previous, ";add.column=", add_column
        ),
        paste0(
          "seqdist(fixture_c(), method=\"OMstran\", with_missing=True, ",
          "norm=\"none\", subm=2.0, indel=1.0, otto=0.5, transindel=\"",
          transindel, "\", previous=", if (previous) "True" else "False",
          ", add_column=", if (add_column) "True" else "False", ")"
        ),
        list(
          seqdata = fx$C, method = "OMstran", with.missing = TRUE,
          norm = "none", sm = sms$C, indel = 1, otto = 0.5,
          transindel = transindel, previous = previous,
          add.column = add_column
        )
      )
    }
  }
}

for (use_weights in c(TRUE, FALSE)) {
  add_case(
    paste0(
      "ragged_omstran_prob_attached_",
      if (use_weights) "weighted" else "unweighted"
    ),
    "OMstran",
    "CW",
    paste0(
      "method=OMstran;with.missing=TRUE;norm=none;sm=constant2;",
      "indel=1;otto=0.5;transindel=prob;weighted=", use_weights
    ),
    paste0(
      "seqdist(fixture_cw(), method=\"OMstran\", ",
      "with_missing=True, norm=\"none\", subm=2.0, indel=1.0, ",
      "otto=0.5, transindel=\"prob\", weighted=",
      if (use_weights) "True" else "False", ")"
    ),
    list(
      seqdata = fx$Cw, method = "OMstran", with.missing = TRUE,
      norm = "none", sm = sms$C, indel = 1, otto = 0.5,
      transindel = "prob", weighted = use_weights
    )
  )
}

prob_external_args <- method_args("C", "OMstran", "none")
prob_external_args$transindel <- "prob"
prob_external_args$refseq <- external_ref
add_case(
  "ragged_omstran_prob_refseq_external",
  "OMstran",
  "C",
  paste0(
    "method=OMstran;with.missing=TRUE;norm=none;sm=constant2;",
    "indel=1;otto=0.5;transindel=prob;refseq=external"
  ),
  paste0(
    "seqdist(fixture_c(), method=\"OMstran\", with_missing=True, ",
    "norm=\"none\", subm=2.0, indel=1.0, otto=0.5, ",
    "transindel=\"prob\", refseq=external_ref())"
  ),
  prob_external_args
)

for (transindel in c("constant", "prob")) {
  add_case(
    paste0("empty_omstran_", transindel),
    "OMstran",
    "K",
    paste0(
      "method=OMstran;norm=none;sm=constant2;indel=1;otto=0.5;",
      "transindel=", transindel
    ),
    paste0(
      "seqdist(fixture_k(), method=\"OMstran\", norm=\"none\", ",
      "subm=2.0, indel=1.0, otto=0.5, transindel=\"",
      transindel, "\")"
    ),
    list(
      seqdata = empty_seq, method = "OMstran", norm = "none",
      sm = empty_sm, indel = 1, otto = 0.5,
      transindel = transindel
    )
  )
}

add_case(
  "refseq_most_frequent_tie_numeric_sps",
  "LCS",
  "H",
  "method=LCS;norm=none;refseq=0;frequency_tie=numeric_SPS",
  "seqdist(fixture_h(), method=\"LCS\", norm=\"none\", refseq=0)",
  list(seqdata = tie_seq, method = "LCS", norm = "none", refseq = 0L)
)
add_case(
  "refseq_most_frequent_structural_empty",
  "LCS",
  "M",
  "method=LCS;with.missing=TRUE;norm=none;refseq=0;empty_frequency=2",
  paste0(
    "seqdist(fixture_m(), method=\"LCS\", with_missing=True, ",
    "norm=\"none\", refseq=0)"
  ),
  list(
    seqdata = empty_frequency_seq, method = "LCS",
    with.missing = TRUE, norm = "none", refseq = 0L
  )
)

# DHD fits its time-varying TRATE costs over the complete declared alphabet,
# including states that happen not to occur.  Missing is then overwritten with
# miss.cost, and a shorter external HAM/DHD reference is compared only through
# its own length (the equal-length check applies to source rows before append).
add_case(
  "dhd_declared_unobserved_state_auto",
  "DHD",
  "O",
  "method=DHD;norm=auto;declared_unobserved=C",
  "seqdist(fixture_o(), method=\"DHD\", norm=\"auto\")",
  list(seqdata = declared_state_dhd_seq, method = "DHD", norm = "auto")
)
for (norm in c("none", "auto")) {
  add_case(
    paste0("dhd_ragged_missing_", norm),
    "DHD",
    "P",
    paste0("method=DHD;with.missing=TRUE;norm=", norm, ";trailing_void"),
    paste0(
      "seqdist(fixture_p(), method=\"DHD\", with_missing=True, norm=\"",
      norm, "\")"
    ),
    list(
      seqdata = missing_ragged_dhd_seq, method = "DHD",
      with.missing = TRUE, norm = norm
    )
  )
}
for (method in c("CHI2", "EUCLID")) {
  add_case(
    paste0(tolower(method), "_ragged_empty_interval_auto"),
    method,
    "P",
    paste0(
      "method=", method,
      ";with.missing=TRUE;norm=auto;trailing_void_interval"
    ),
    paste0(
      "seqdist(fixture_p(), method=\"", method,
      "\", with_missing=True, norm=\"auto\")"
    ),
    list(
      seqdata = missing_ragged_dhd_seq, method = method,
      with.missing = TRUE, norm = "auto"
    )
  )
}
for (method in c("HAM", "DHD")) {
  add_case(
    paste0("complete_", tolower(method), "_refseq_shorter"),
    method,
    "G",
    paste0("method=", method, ";norm=none;refseq=shorter"),
    paste0(
      "seqdist(fixture_g(), method=\"", method,
      "\", norm=\"none\", refseq=short_ref())"
    ),
    list(
      seqdata = asym_seq, method = method, norm = "none", refseq = short_ref
    )
  )
}

# State-dependent indels and the only explicit position-dependent substitution
# cost shape accepted by seqdist: a DHD (states x states x positions) array.
add_case(
  "ragged_om_vector_indel",
  "OM",
  "C",
  "method=OM;with.missing=TRUE;norm=none;sm=constant2;indel=1|1.25|1.5|1.75",
  paste0(
    "seqdist(fixture_c(), method=\"OM\", with_missing=True, norm=\"none\", ",
    "subm=2.0, indel=[1.0, 1.25, 1.5, 1.75])"
  ),
  list(
    seqdata = fx$C, method = "OM", with.missing = TRUE, norm = "none",
    sm = sms$C, indel = c(1, 1.25, 1.5, 1.75)
  )
)
for (method in c("OMslen", "OMspell", "TWED")) {
  method_options <- switch(
    method,
    OMslen = list(sm = matrix(c(0, 2, 2, 0), nrow = 2), link = "mean", h = 0.5),
    OMspell = list(sm = matrix(c(0, 2, 2, 0), nrow = 2), expcost = 0.5, tpow = 1),
    TWED = list(sm = matrix(c(0, 2, 2, 0), nrow = 2), nu = 0.5, h = 0.5)
  )
  python_options <- switch(
    method,
    OMslen = ", subm=2.0, indel=[1.0, 2.0], link=\"mean\", h=0.5",
    OMspell = ", subm=2.0, indel=[1.0, 2.0], expcost=0.5, tpow=1.0",
    TWED = ", subm=2.0, indel=[1.0, 2.0], nu=0.5, h=0.5"
  )
  add_case(
    paste0("vector_indel_", tolower(method)),
    method,
    "N",
    paste0("method=", method, ";norm=none;sm=constant2;indel=1|2"),
    paste0(
      "seqdist(fixture_n(), method=\"", method, "\", norm=\"none\"",
      python_options, ")"
    ),
    c(
      list(
        seqdata = vector_indel_seq, method = method, norm = "none",
        indel = c(1, 2)
      ),
      method_options
    ),
    allow_warnings = method == "TWED"
  )
}
dhd_position_sm <- array(
  0,
  dim = c(4, 4, 6),
  dimnames = list(
    c("A", "B", "C", "*"),
    c("A", "B", "C", "*"),
    as.character(seq_len(6))
  )
)
for (time in seq_len(6)) {
  dhd_position_sm[, , time] <- time
  diag(dhd_position_sm[, , time]) <- 0
}
add_case(
  "interior_dhd_position_subm",
  "DHD",
  "B",
  "method=DHD;with.missing=TRUE;norm=none;sm=position_offdiag_1_to_6",
  paste0(
    "seqdist(fixture_b(), method=\"DHD\", with_missing=True, norm=\"none\", ",
    "subm=position_subm())"
  ),
  list(
    seqdata = fx$B, method = "DHD", with.missing = TRUE,
    norm = "none", sm = dhd_position_sm
  )
)

# DHD accepts generated time-varying costs by method name, but never a static
# matrix or an array whose time axis differs from the original stored width.
# INDELS/INDELSLOG deliberately inherit seqstatd's attached-weight default,
# even when seqdist itself is asked to be unweighted.
for (cost_method in c("TRATE", "INDELS", "INDELSLOG")) {
  add_case(
    paste0("interior_dhd_", tolower(cost_method), "_method"),
    "DHD",
    "B",
    paste0(
      "method=DHD;with.missing=TRUE;norm=none;sm=", cost_method
    ),
    paste0(
      "seqdist(fixture_b(), method=\"DHD\", with_missing=True, ",
      "norm=\"none\", subm=\"", cost_method, "\")"
    ),
    list(
      seqdata = fx$B, method = "DHD", with.missing = TRUE,
      norm = "none", sm = cost_method
    )
  )
}
add_case(
  "interior_dhd_indels_attached_unweighted",
  "DHD",
  "W",
  paste0(
    "method=DHD;with.missing=TRUE;norm=none;sm=INDELS;",
    "weighted=FALSE;attached_weights_still_used"
  ),
  paste0(
    "seqdist(fixture_bw(), method=\"DHD\", with_missing=True, ",
    "norm=\"none\", subm=\"INDELS\", weighted=False)"
  ),
  list(
    seqdata = fx$Bw, method = "DHD", with.missing = TRUE,
    norm = "none", sm = "INDELS", weighted = FALSE
  )
)
add_case(
  "complete_dhd_indels_refseq_shorter",
  "DHD",
  "G",
  "method=DHD;norm=none;sm=INDELS;refseq=shorter",
  paste0(
    "seqdist(fixture_g(), method=\"DHD\", norm=\"none\", ",
    "subm=\"INDELS\", refseq=short_ref())"
  ),
  list(
    seqdata = asym_seq, method = "DHD", norm = "none",
    sm = "INDELS", refseq = short_ref
  )
)

dhd_ragged_position_sm <- array(
  0,
  dim = c(3, 3, 4),
  dimnames = list(
    c("A", "B", "*"), c("A", "B", "*"), as.character(seq_len(4))
  )
)
for (time in seq_len(4)) {
  dhd_ragged_position_sm[, , time] <- time
  diag(dhd_ragged_position_sm[, , time]) <- 0
}
add_case(
  "dhd_ragged_position_subm_stored_width",
  "DHD",
  "P",
  paste0(
    "method=DHD;with.missing=TRUE;norm=none;sm=3x3x4;",
    "logical_width=3;stored_width=4"
  ),
  paste0(
    "seqdist(fixture_p(), method=\"DHD\", with_missing=True, ",
    "norm=\"none\", subm=ragged_position_subm())"
  ),
  list(
    seqdata = missing_ragged_dhd_seq, method = "DHD",
    with.missing = TRUE, norm = "none", sm = dhd_ragged_position_sm
  )
)
add_case(
  "refuse_dhd_canonical_width_subm_on_ragged",
  "DHD",
  "P",
  paste0(
    "method=DHD;with.missing=TRUE;norm=none;sm=3x3x3;",
    "logical_width=3;stored_width=4"
  ),
  paste0(
    "seqdist(fixture_p(), method=\"DHD\", with_missing=True, ",
    "norm=\"none\", subm=ragged_canonical_subm())"
  ),
  list(
    seqdata = missing_ragged_dhd_seq, method = "DHD",
    with.missing = TRUE, norm = "none",
    sm = dhd_ragged_position_sm[, , seq_len(3), drop = FALSE]
  )
)

dhd_wrong_time_sm <- array(
  0,
  dim = c(2, 2, 4),
  dimnames = list(c("A", "B"), c("A", "B"), as.character(seq_len(4)))
)
for (time in seq_len(4)) {
  dhd_wrong_time_sm[, , time] <- 2
  diag(dhd_wrong_time_sm[, , time]) <- 0
}
add_case(
  "refuse_dhd_static_subm",
  "DHD",
  "N",
  "method=DHD;norm=none;sm=static_2x2",
  paste0(
    "seqdist(fixture_n(), method=\"DHD\", norm=\"none\", ",
    "subm=static_subm())"
  ),
  list(
    seqdata = vector_indel_seq, method = "DHD", norm = "none",
    sm = matrix(c(0, 2, 2, 0), nrow = 2)
  )
)
for (external in c(FALSE, TRUE)) {
  add_case(
    paste0(
      "refuse_dhd_wrong_time_subm",
      if (external) "_refseq_shorter" else ""
    ),
    "DHD",
    if (external) "G" else "N",
    paste0(
      "method=DHD;norm=none;sm=2x2x4;stored_width=3",
      if (external) ";refseq=shorter" else ""
    ),
    paste0(
      "seqdist(", if (external) "fixture_g()" else "fixture_n()",
      ", method=\"DHD\", norm=\"none\", subm=wrong_time_subm()",
      if (external) ", refseq=short_ref()" else "", ")"
    ),
    c(
      list(
        seqdata = if (external) asym_seq else vector_indel_seq,
        method = "DHD", norm = "none", sm = dhd_wrong_time_sm
      ),
      if (external) list(refseq = short_ref) else list()
    )
  )
}

dhd_negative_sm <- array(0, dim = c(2, 2, 3))
dhd_negative_sm[1, 2, ] <- 2
dhd_negative_sm[2, 1, ] <- 2
dhd_negative_sm[1, 2, 1] <- -1
dhd_negative_sm[2, 1, 1] <- -1
dhd_diagonal_sm <- array(0, dim = c(2, 2, 3))
dhd_diagonal_sm[1, 2, ] <- 2
dhd_diagonal_sm[2, 1, ] <- 2
dhd_diagonal_sm[1, 1, 1] <- 1
for (invalid_kind in c("negative", "diagonal")) {
  add_case(
    paste0("refuse_dhd_", invalid_kind, "_subm"),
    "DHD",
    "N",
    paste0("method=DHD;norm=none;sm=2x2x3;", invalid_kind),
    paste0(
      "seqdist(fixture_n(), method=\"DHD\", norm=\"none\", subm=",
      invalid_kind, "_dhd_subm())"
    ),
    list(
      seqdata = vector_indel_seq, method = "DHD", norm = "none",
      sm = if (invalid_kind == "negative") dhd_negative_sm else dhd_diagonal_sm
    )
  )
}
for (nonfinite_kind in c("nan", "infinite")) {
  nonfinite_sm <- array(0, dim = c(2, 2, 3))
  nonfinite_sm[1, 2, ] <- 2
  nonfinite_sm[2, 1, ] <- 2
  nonfinite_value <- if (nonfinite_kind == "nan") NaN else Inf
  nonfinite_sm[1, 2, 1] <- nonfinite_value
  nonfinite_sm[2, 1, 1] <- nonfinite_value
  add_case(
    paste0("refuse_dhd_", nonfinite_kind, "_subm"),
    "DHD",
    "N",
    paste0("method=DHD;norm=none;sm=2x2x3;", nonfinite_kind),
    paste0(
      "seqdist(fixture_n(), method=\"DHD\", norm=\"none\", subm=",
      nonfinite_kind, "_dhd_subm())"
    ),
    list(
      seqdata = vector_indel_seq, method = "DHD", norm = "none",
      sm = nonfinite_sm
    )
  )
}

for (method in c("OMslen", "OMspell", "OMstran")) {
  extra_args <- switch(
    method,
    OMslen = list(link = "mean", h = 0.5),
    OMspell = list(expcost = 0.5, tpow = 1),
    OMstran = list(otto = 0.5, transindel = "constant")
  )
  python_extra <- switch(
    method,
    OMslen = ", link=\"mean\", h=0.5",
    OMspell = ", expcost=0.5, tpow=1.0",
    OMstran = ", otto=0.5, transindel=\"constant\""
  )
  add_case(
    paste0("refuse_nonpositive_indel_", tolower(method)),
    method,
    "N",
    paste0("method=", method, ";norm=none;sm=constant2;indel=[1,-1]"),
    paste0(
      "seqdist(fixture_n(), method=\"", method,
      "\", norm=\"none\", subm=2.0, indel=[1.0, -1.0]",
      python_extra, ")"
    ),
    c(
      list(
        seqdata = vector_indel_seq, method = method, norm = "none",
        sm = matrix(c(0, 2, 2, 0), nrow = 2), indel = c(1, -1)
      ),
      extra_args
    )
  )
}

# Positive infinity passes checkcost for indels but each compiled method then
# has its own established disposition. Freeze scalar and state-vector forms.
for (indel_kind in c("scalar", "vector")) {
  infinite_methods <- if (indel_kind == "scalar") {
    c("OMslen", "OMspell", "OMstran", "TWED")
  } else {
    c("OM", "OMslen", "OMspell", "OMstran", "TWED")
  }
  for (method in infinite_methods) {
    method_extra <- switch(
      method,
      OM = list(),
      OMslen = list(link = "mean", h = 0.5),
      OMspell = list(expcost = 0.5, tpow = 1),
      OMstran = list(
        otto = 0.5, transindel = "constant",
        previous = FALSE, add.column = TRUE
      ),
      TWED = list(nu = 0.5, h = 0.5)
    )
    python_extra <- switch(
      method,
      OM = "",
      OMslen = ", link=\"mean\", h=0.5",
      OMspell = ", expcost=0.5, tpow=1.0",
      OMstran = paste0(
        ", otto=0.5, transindel=\"constant\", previous=False, ",
        "add_column=True"
      ),
      TWED = ", nu=0.5, h=0.5"
    )
    python_indel <- if (indel_kind == "scalar") {
      "infinite_value()"
    } else {
      "[infinite_value(), 1.0]"
    }
    r_indel <- if (indel_kind == "scalar") Inf else c(Inf, 1)
    add_case(
      paste0(
        if (method == "OMstran") "refuse_" else "",
        "infinite_indel_", indel_kind, "_", tolower(method)
      ),
      method,
      "N",
      paste0(
        "method=", method, ";norm=none;sm=constant2;indel=",
        if (indel_kind == "scalar") "Inf" else "[Inf,1]"
      ),
      paste0(
        "seqdist(fixture_n(), method=\"", method,
        "\", norm=\"none\", subm=2.0, indel=", python_indel,
        python_extra, ")"
      ),
      c(
        list(
          seqdata = vector_indel_seq, method = method, norm = "none",
          sm = matrix(c(0, 2, 2, 0), nrow = 2), indel = r_indel
        ),
        method_extra
      ),
      allow_warnings = indel_kind == "vector"
    )
  }
}
add_case(
  "refuse_nonpositive_indel_twed",
  "TWED",
  "N",
  "method=TWED;norm=none;sm=constant2;indel=[-1,-2];max_coercion",
  paste0(
    "seqdist(fixture_n(), method=\"TWED\", norm=\"none\", ",
    "subm=2.0, indel=[-1.0, -2.0], nu=0.5, h=0.5)"
  ),
  list(
    seqdata = vector_indel_seq, method = "TWED", norm = "none",
    sm = matrix(c(0, 2, 2, 0), nrow = 2), indel = c(-1, -2),
    nu = 0.5, h = 0.5
  ),
  allow_warnings = TRUE
)
add_case(
  "refuse_zero_auto_indel_om",
  "OM",
  "N",
  "method=OM;norm=none;sm=zeros;indel=auto",
  "seqdist(fixture_n(), method=\"OM\", norm=\"none\", subm=0.0)",
  list(
    seqdata = vector_indel_seq, method = "OM", norm = "none",
    sm = matrix(0, nrow = 2, ncol = 2)
  )
)
add_case(
  "refuse_nan_indel_om",
  "OM",
  "N",
  "method=OM;norm=none;sm=constant2;indel=NaN",
  paste0(
    "seqdist(fixture_n(), method=\"OM\", norm=\"none\", ",
    "subm=2.0, indel=nan_value())"
  ),
  list(
    seqdata = vector_indel_seq, method = "OM", norm = "none",
    sm = matrix(c(0, 2, 2, 0), nrow = 2), indel = NaN
  )
)
add_case(
  "accept_infinite_indel_om",
  "OM",
  "N",
  "method=OM;norm=none;sm=constant2;indel=Inf",
  paste0(
    "seqdist(fixture_n(), method=\"OM\", norm=\"none\", ",
    "subm=2.0, indel=infinite_value())"
  ),
  list(
    seqdata = vector_indel_seq, method = "OM", norm = "none",
    sm = matrix(c(0, 2, 2, 0), nrow = 2), indel = Inf
  )
)
for (nonfinite_kind in c("nan", "infinite")) {
  nonfinite_value <- if (nonfinite_kind == "nan") NaN else Inf
  nonfinite_sm <- matrix(
    c(0, nonfinite_value, nonfinite_value, 0),
    nrow = 2
  )
  add_case(
    paste0("refuse_om_", nonfinite_kind, "_subm"),
    "OM",
    "N",
    paste0("method=OM;norm=none;sm=", nonfinite_kind, ";indel=1"),
    paste0(
      "seqdist(fixture_n(), method=\"OM\", norm=\"none\", subm=",
      nonfinite_kind, "_value(), indel=1.0)"
    ),
    list(
      seqdata = vector_indel_seq, method = "OM", norm = "none",
      sm = nonfinite_sm, indel = 1
    ),
    allow_warnings = nonfinite_kind == "infinite"
  )
}
add_case(
  "refuse_nan_indel_twed",
  "TWED",
  "N",
  "method=TWED;norm=none;sm=constant2;indel=[NaN,1];max_coercion",
  paste0(
    "seqdist(fixture_n(), method=\"TWED\", norm=\"none\", ",
    "subm=2.0, indel=[nan_value(), 1.0], nu=0.5, h=0.5)"
  ),
  list(
    seqdata = vector_indel_seq, method = "TWED", norm = "none",
    sm = matrix(c(0, 2, 2, 0), nrow = 2), indel = c(NaN, 1),
    nu = 0.5, h = 0.5
  ),
  allow_warnings = TRUE
)

# Explicit missing-state refusal (all successful missing/ragged rows above use
# with.missing=TRUE).
add_case(
  "interior_lcs_missing_refused",
  "LCS",
  "B",
  "method=LCS;with.missing=FALSE;norm=none",
  "seqdist(fixture_b(), method=\"LCS\", with_missing=False, norm=\"none\")",
  list(seqdata = fx$B, method = "LCS", with.missing = FALSE, norm = "none")
)

normalize_error <- function(message, method, case) {
  clean <- trimws(gsub("[[:space:]]+", " ", message))
  if (grepl("missing value where TRUE/FALSE needed", clean, fixed = TRUE)) {
    if (method == "DHD") {
      return(c(
        "invalid_dhd_costs",
        "DHD substitution costs must be finite, nonnegative, and zero-diagonal"
      ))
    }
    if (grepl("_subm", case, fixed = TRUE)) {
      return(c(
        "invalid_subm",
        "substitution costs must be finite"
      ))
    }
    if (method == "OMstran") {
      return(c(
        "invalid_indel",
        "indel cost must be finite and strictly positive"
      ))
    }
    return(c(
      "invalid_indel",
      "indel cost must be strictly positive"
    ))
  }
  # TraMineR 2.2-14 rewrote checkcost()'s diagnostics. The contracts are the
  # same; only the wording changed, so both spellings map to one class:
  #   2.2-13 "indel cost should be positive"      -> 2.2-14 "negative indel value(s)!"
  #   2.2-13 "All element on the diagonal of sm"  -> 2.2-14 "Non-zero diagonal elements ..."
  #   2.2-13 "You should specify either one ..."  -> 2.2-14 "You must specify either a single ..."
  if (grepl("indel cost should be positive", clean, ignore.case = TRUE) ||
      grepl("negative indel value", clean, ignore.case = TRUE)) {
    return(c(
      "invalid_indel",
      "indel cost must be strictly positive"
    ))
  }
  if (grepl("must specify either a single indel", clean, ignore.case = TRUE) ||
      grepl("should specify either one global indel", clean, ignore.case = TRUE)) {
    return(c(
      "invalid_indel",
      "indel cost must be strictly positive"
    ))
  }
  if (method == "DHD" &&
      (grepl("'sm' invalid type", clean, fixed = TRUE) ||
       grepl("size of substitution cost matrix", clean, fixed = TRUE) ||
       # TraMineR 2.2-14 reports a 2-D sm handed to DHD as "incorrect number of
       # dimensions"; 2.2-13 reached the same refusal through the opaque
       # "missing value where TRUE/FALSE needed". Same contract, clearer text.
       grepl("incorrect number of dimensions", clean, fixed = TRUE))) {
    return(c(
      "invalid_dhd_subm",
      "DHD subm must match the active states and original stored width"
    ))
  }
  if (method == "DHD" &&
      (grepl("Negative substitution costs", clean, fixed = TRUE) ||
       grepl("diagonal", clean, ignore.case = TRUE))) {
    return(c(
      "invalid_dhd_costs",
      "DHD substitution costs must be finite, nonnegative, and zero-diagonal"
    ))
  }
  if (grepl("with.missing", clean, fixed = TRUE) &&
      grepl("must be TRUE", clean, fixed = TRUE)) {
    return(c(
      "missing_requires_opt_in",
      "with_missing must be true when sequences contain missing values"
    ))
  }
  if (grepl("different length", clean, fixed = TRUE)) {
    return(c(
      "unequal_length",
      paste0(method, " is not defined for sequences of different length")
    ))
  }
  if (grepl("norm", clean, ignore.case = TRUE) &&
      (grepl("implemented", clean, ignore.case = TRUE) ||
       grepl("only be one", clean, ignore.case = TRUE))) {
    return(c(
      "unsupported_norm",
      paste0("norm is not supported for method ", method)
    ))
  }
  stop("unclassified R error for ", method, ": ", clean)
}

numeric_rows <- list()
error_rows <- character()

append_numeric <- function(spec, result) {
  kind <- if (inherits(result, "dist")) {
    "condensed"
  } else if (is.matrix(result)) {
    "matrix"
  } else {
    "vector"
  }

  common <- c(
    csv_quote(spec$case),
    csv_quote(spec$method),
    csv_quote(spec$fixture),
    csv_quote(spec$options),
    csv_quote(spec$python_call),
    csv_quote(kind)
  )

  if (kind == "matrix") {
    nr <- nrow(result)
    nc <- ncol(result)
    for (i in seq_len(nr)) {
      for (j in seq_len(nc)) {
        numeric_rows[[length(numeric_rows) + 1L]] <<- c(
          common,
          csv_quote(nr),
          csv_quote(nc),
          csv_quote((i - 1L) * nc + (j - 1L)),
          csv_quote(i - 1L),
          csv_quote(j - 1L),
          format_value(result[i, j])
        )
      }
    }
  } else {
    values <- as.numeric(result)
    for (i in seq_along(values)) {
      numeric_rows[[length(numeric_rows) + 1L]] <<- c(
        common,
        csv_quote(length(values)),
        csv_quote(""),
        csv_quote(i - 1L),
        csv_quote(""),
        csv_quote(""),
        format_value(values[[i]])
      )
    }
  }
}

for (spec in case_specs) {
  captured_warnings <- character()
  captured_error <- NULL
  result <- tryCatch(
    suppressMessages(withCallingHandlers(
      do.call(seqdist, spec$args),
      warning = function(warning) {
        captured_warnings <<- c(captured_warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }
    )),
    error = function(error) {
      captured_error <<- error
      NULL
    }
  )

  if (length(captured_warnings) && !isTRUE(spec$allow_warnings)) {
    stop(
      "unexpected warning in distance case ", spec$case, ": ",
      paste(captured_warnings, collapse = " || ")
    )
  }

  if (is.null(captured_error)) {
    append_numeric(spec, result)
  } else {
    normalized <- normalize_error(
      conditionMessage(captured_error), spec$method, spec$case
    )
    error_rows <- c(error_rows, paste(
      spec$case,
      spec$method,
      spec$fixture,
      spec$options,
      spec$python_call,
      normalized[[1]],
      normalized[[2]],
      sep = "|"
    ))
  }
}

numeric_path <- file.path(out_dir, "distances_r1_cases.csv")
numeric_con <- file(numeric_path, open = "wt", encoding = "UTF-8")
on.exit(close(numeric_con), add = TRUE)
writeLines(
  paste(
    c(
      "case", "method", "fixture", "options", "python_call", "output_kind",
      "shape_0", "shape_1", "index", "row", "col", "value"
    ),
    collapse = ","
  ),
  con = numeric_con,
  useBytes = TRUE
)
if (length(numeric_rows)) {
  writeLines(
    vapply(numeric_rows, paste, collapse = ",", FUN.VALUE = character(1)),
    con = numeric_con,
    useBytes = TRUE
  )
}
close(numeric_con)
on.exit(NULL, add = FALSE)

writeLines(
  c(
    "case|method|fixture|options|python_call|category|message",
    error_rows
  ),
  con = file.path(out_dir, "distances_r1_errors.txt"),
  useBytes = TRUE
)

numeric_cases <- unique(vapply(numeric_rows, `[[`, character(1), 1L))
numeric_cases <- gsub('^"|"$', "", numeric_cases)
stopifnot(
  length(unique(c(numeric_cases, sub("\\|.*$", "", error_rows)))) ==
    length(case_specs)
)
cat(
  "DISTANCES_R1_GOLDENS_DONE cases=", length(case_specs),
  " numeric=", length(unique(numeric_cases)),
  " errors=", length(error_rows),
  " rows=", length(numeric_rows), "\n",
  sep = ""
)
