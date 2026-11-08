# Live-R goldens for seqdist's argument surface (TX-AUD-20261010-071).
#
# seqdist accepts every one of its arguments for every method (seqdist.R:4-11). It
# validates some for every method (norm :130-132, indel :137-147, sm :154-168, prox
# :172-180, link :185-186, step :189-190), some only for the method that uses them
# (:206-281), and otherwise ignores an argument the method does not use. For each
# of R's 16 methods the rows record: the call with only the method's required
# arguments; the call with every defaulted argument given explicitly; the call with
# each argument the method does not use at a non-default value (some of them values
# R would refuse for a method that uses them); and R's checks that stop any method.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_seqdist_r_arguments.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sx <- suppressMessages(TraMineR::seqdef(
  data.frame(a = c("A", "B", "A", "C"), b = c("B", "B", "A", "C"), c = c("A", "C", "B", "C"), d = c("C", "A", "B", "A")),
  alphabet = c("A", "B", "C")))
methods <- c("OM", "OMloc", "OMslen", "OMspell", "OMstran", "HAM", "DHD", "CHI2", "EUCLID", "LCS", "LCP", "RLCP",
             "NMS", "NMSMST", "SVRspell", "TWED")
required <- list(OM = list(sm = "CONSTANT"), OMloc = list(sm = "CONSTANT"), OMslen = list(sm = "CONSTANT"),
                 OMspell = list(sm = "CONSTANT"), OMstran = list(sm = "CONSTANT", otto = 0.5),
                 TWED = list(sm = "CONSTANT", nu = 0.5))
# The arguments each method uses (seqdist.R:206-281, 604-640, 700-815); all others are ignored by R.
uses <- list(indel = c("OM", "OMslen", "OMspell", "OMstran", "TWED"),
             # LCS replaces sm by constant costs (seqdist.R:370-374); NMS and SVRspell stop on sm (:293-294).
             sm = c("OM", "OMloc", "OMslen", "OMspell", "OMstran", "HAM", "DHD", "TWED", "NMS", "SVRspell"),
             kweights = c("NMS", "NMSMST", "SVRspell"), tpow = c("OMspell", "NMSMST", "SVRspell"),
             expcost = c("OMloc", "OMspell"), context = "OMloc", link = "OMslen", h = c("OMslen", "TWED"),
             nu = "TWED", transindel = "OMstran", otto = "OMstran", previous = "OMstran", add.column = "OMstran",
             breaks = c("CHI2", "EUCLID"), step = c("CHI2", "EUCLID"), overlap = c("CHI2", "EUCLID"),
             global.pdotj = c("CHI2", "EUCLID"), prox = c("NMS", "SVRspell"), opt.args = "OMspell")
defaults <- list(refseq = NULL, norm = "none", indel = "auto", sm = NULL, with.missing = FALSE, full.matrix = TRUE,
                 kweights = rep(1.0, ncol(sx)), tpow = 1.0, expcost = 0.5, link = "mean", h = 0.5,
                 transindel = "constant", previous = FALSE, add.column = TRUE, breaks = NULL, step = 1,
                 overlap = FALSE, weighted = TRUE, global.pdotj = NULL, prox = NULL, check.max.size = TRUE,
                 check.costs = TRUE, opt.args = list())
# A non-default value for an argument the method does not use; the second set holds values R refuses for a
# method that uses them (seqdist.R:213-214, 234-236, 243-245, 266-268, 279-283), so R's ignoring is visible.
ignored_values <- list(indel = 3, sm = "INDELS", kweights = rep(2.0, ncol(sx)), tpow = 2, expcost = 0.25,
                       context = 0.3, link = "gmean", h = 1, nu = 0.7, transindel = "subcost", otto = 0.3,
                       previous = TRUE, add.column = FALSE, breaks = list(c(1, 2)), step = 2, overlap = TRUE,
                       global.pdotj = c(0.5, 0.25, 0.25), prox = diag(3), opt.args = list(tokdep.coeff = c(1, 1, 1)),
                       check.max.size = FALSE)
invalid_ignored <- list(expcost = -1, h = -1, nu = 0, otto = 2, kweights = -1, tpow = -1)
global_stops <- list(
  list(case = "global_link_bogus", method = "OM", args = list(sm = "CONSTANT", link = "bogus")),
  list(case = "global_step_fraction", method = "LCS", args = list(step = 2.5)),
  list(case = "global_prox_not_matrix", method = "OM", args = list(sm = "CONSTANT", prox = 1)),
  list(case = "global_sm_bogus", method = "CHI2", args = list(sm = "bogus")),
  list(case = "global_indel_bogus", method = "HAM", args = list(indel = "bogus")),
  list(case = "global_indel_vector_length", method = "HAM", args = list(indel = c(1, 2))),
  list(case = "nms_sm_given", method = "NMS", args = list(sm = "CONSTANT")),
  list(case = "svrspell_sm_given", method = "SVRspell", args = list(sm = "CONSTANT")),
  list(case = "nms_norm_auto", method = "NMS", args = list(norm = "auto")),
  list(case = "nmsmst_norm_auto", method = "NMSMST", args = list(norm = "auto")))
# R's global sm check stops on anything but NULL, a matrix, an array or a method name,
# for every method (seqdist.R:153-168); opt.args is read only when it is a list, and
# only its tokdep.coeff entry (:690-691) (batch-6 review M-2, M-3).
argument_forms <- list(
  list(case = "global_sm_number", method = "CHI2", args = list(sm = 2)),
  list(case = "global_sm_number_lcs", method = "LCS", args = list(sm = 2)),
  list(case = "global_sm_vector", method = "LCP", args = list(sm = c(1, 2))),
  list(case = "omspell_opt_args_number", method = "OMspell", args = list(sm = "CONSTANT", opt.args = 5)),
  list(case = "omspell_opt_args_unnamed_list", method = "OMspell",
       args = list(sm = "CONSTANT", opt.args = list("tokdep.coeff"))),
  list(case = "omspell_opt_args_null_coeff", method = "OMspell",
       args = list(sm = "CONSTANT", opt.args = list(tokdep.coeff = NULL))))

rows <- list()
run <- function(case, method, argument, args) {
  result <- tryCatch(suppressWarnings(suppressMessages(do.call(TraMineR::seqdist, c(list(seqdata = sx, method = method), args)))),
                     error = function(e) e)
  ok <- !inherits(result, "error")
  value <- if (ok) paste(format(as.vector(as.matrix(result)), digits = 17), collapse = " ") else ""
  rows[[length(rows) + 1L]] <<- data.frame(case = case, method = method, argument = argument, ok = ok,
                                           error = if (ok) "" else conditionMessage(result), value = value,
                                           stringsAsFactors = FALSE)
}
for (m in methods) {
  base <- if (is.null(required[[m]])) list() else required[[m]]
  run(paste0(m, "_base"), m, "", base)
  all_args <- base
  for (a in names(defaults)) if (!(a == "sm" && !is.null(base$sm))) all_args[a] <- list(defaults[[a]])
  run(paste0(m, "_all_defaults"), m, "", all_args)
  for (a in names(ignored_values)) {
    if (a %in% names(uses) && m %in% uses[[a]]) next
    if (a %in% names(base)) next
    args <- base; args[a] <- list(ignored_values[[a]])
    run(paste0(m, "_ignored_", a), m, a, args)
  }
  for (a in names(invalid_ignored)) {
    if (a %in% names(uses) && m %in% uses[[a]]) next
    if (a %in% names(base)) next
    args <- base; args[a] <- list(invalid_ignored[[a]])
    run(paste0(m, "_ignored_invalid_", a), m, a, args)
  }
}
for (s in global_stops) run(s$case, s$method, names(s$args)[length(s$args)], s$args)
for (s in argument_forms) run(s$case, s$method, names(s$args)[length(s$args)], s$args)

write.csv(do.call(rbind, rows), file.path(out_dir, "seqdist_r_arguments.csv"), row.names = FALSE, quote = TRUE)
