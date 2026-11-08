"""Which Python defaults do the R goldens never exercise, and are they right? (W-07)

A golden that passes an argument explicitly cannot test that argument's
default. That is the pattern that hid the ``seqmaintokens`` bug in cycle 6, and
a second instance was found in the same cycle; nobody had swept for the rest.

This module does the sweep mechanically and keeps doing it. It reads every
generator in ``scripts/r_goldens/``, extracts each R call with its named
arguments, maps the R function to its traminex counterpart, and classifies
every parameter that has a default:

* **never exercised** -- every call passes it, so the default is never taken;
* **only the default exercised** -- no call passes it;
* **both** -- some calls pass it and some do not.

For the first group there is no golden left to appeal to, so the check becomes
whether traminex's default agrees with R's own formal default, read from
``tests/r_goldens/r_formals.csv`` (generator:
``scripts/r_goldens/make_goldens_r_formals.R``).

The denominator and the exact never-exercised set are asserted, so adding a
generator that pins a new argument, or adding a parameter with a default, fails
here until it is looked at.
"""

from __future__ import annotations

import csv
import inspect
import re
from collections import defaultdict
from pathlib import Path

import pytest

import traminex as tx
from traminex.api_manifest import PUBLIC_API

ROOT = Path(__file__).resolve().parents[1]
GENERATOR_DIR = ROOT / "scripts" / "r_goldens"
FORMALS = ROOT / "tests" / "r_goldens" / "r_formals.csv"

_CALL = re.compile(r"(?<![\w.$@])([A-Za-z][A-Za-z0-9._]*)\s*\(")
_NAMED = re.compile(r"^\s*([A-Za-z.][A-Za-z0-9._]*)\s*=(?!=)")


# --------------------------------------------------------------------------- #
# A small R reader: enough to find calls and their named arguments.
# --------------------------------------------------------------------------- #
def _strip_comments(text: str) -> str:
    out: list[str] = []
    quote = None
    index = 0
    while index < len(text):
        char = text[index]
        if quote:
            out.append(char)
            if char == "\\" and index + 1 < len(text):
                out.append(text[index + 1])
                index += 2
                continue
            if char == quote:
                quote = None
        elif char in "\"'":
            quote = char
            out.append(char)
        elif char == "#":
            while index < len(text) and text[index] != "\n":
                index += 1
            out.append("\n")
            continue
        else:
            out.append(char)
        index += 1
    return "".join(out)


def _split_top_level(text: str) -> list[str]:
    parts: list[str] = []
    current: list[str] = []
    depth = 0
    quote = None
    index = 0
    while index < len(text):
        char = text[index]
        if quote:
            current.append(char)
            if char == "\\" and index + 1 < len(text):
                current.append(text[index + 1])
                index += 2
                continue
            if char == quote:
                quote = None
        elif char in "\"'":
            quote = char
            current.append(char)
        elif char in "([{":
            depth += 1
            current.append(char)
        elif char in ")]}":
            depth -= 1
            current.append(char)
        elif char == "," and depth == 0:
            parts.append("".join(current))
            current = []
        else:
            current.append(char)
        index += 1
    if current:
        parts.append("".join(current))
    return parts


def _calls(text: str) -> list[tuple[str, set[str]]]:
    found: list[tuple[str, set[str]]] = []
    for match in _CALL.finditer(text):
        opening = match.end() - 1
        depth = 0
        index = opening
        quote = None
        while index < len(text):
            char = text[index]
            if quote:
                if char == "\\":
                    index += 2
                    continue
                if char == quote:
                    quote = None
            elif char in "\"'":
                quote = char
            elif char == "(":
                depth += 1
            elif char == ")":
                depth -= 1
                if depth == 0:
                    break
            index += 1
        if index >= len(text):
            continue
        named = set()
        for part in _split_top_level(text[opening + 1 : index]):
            hit = _NAMED.match(part)
            if hit:
                named.add(hit.group(1))
        found.append((match.group(1), named))
    return found


# --------------------------------------------------------------------------- #
# Name mapping (docs/API_CONSISTENCY.md section 1b).
# --------------------------------------------------------------------------- #
def _spellings(parameter: str) -> set[str]:
    out = {parameter, parameter.replace("_", "."), parameter.rstrip("_")}
    camel = re.sub(r"_([a-z])", lambda m: m.group(1).upper(), parameter)
    out.add(camel)
    if parameter == "subm":
        out.add("sm")
    if parameter in {"with_missing", "include_missing"}:
        out.update({"with.missing", "withmissing"})
    return out


def _python_callable(r_name: str) -> str | None:
    for candidate in (r_name, r_name.replace(".", "_")):
        if candidate in PUBLIC_API and callable(getattr(tx, candidate, None)):
            return candidate
    return None


def _defaults(func) -> dict[str, object]:
    try:
        signature = inspect.signature(func)
    except (TypeError, ValueError):  # pragma: no cover - builtins
        return {}
    return {
        parameter.name: parameter.default
        for parameter in signature.parameters.values()
        if parameter.default is not inspect.Parameter.empty
        and parameter.kind not in (parameter.VAR_POSITIONAL, parameter.VAR_KEYWORD)
    }


# --------------------------------------------------------------------------- #
# The sweep, computed once.
# --------------------------------------------------------------------------- #
def _sweep():
    generators = sorted(GENERATOR_DIR.glob("*.R"))
    per_function: dict[str, list[set[str]]] = defaultdict(list)
    for path in generators:
        for name, named in _calls(_strip_comments(path.read_text(errors="replace"))):
            per_function[name].append(named)

    never: dict[tuple[str, str], object] = {}
    r_names: dict[str, str] = {}
    only_default = 0
    both = 0
    functions = 0
    total = 0
    for r_name in sorted(per_function):
        python_name = _python_callable(r_name)
        if python_name is None:
            continue
        defaults = _defaults(getattr(tx, python_name))
        if not defaults:
            continue
        functions += 1
        r_names.setdefault(python_name, r_name)
        call_list = per_function[r_name]
        for parameter, default in sorted(defaults.items()):
            total += 1
            aliases = _spellings(parameter)
            passing = sum(1 for named in call_list if named & aliases)
            if passing == len(call_list):
                never[(python_name, parameter)] = default
            elif passing == 0:
                only_default += 1
            else:
                both += 1
    return {
        "generators": len(generators),
        "functions": functions,
        "total": total,
        "never": never,
        "r_names": r_names,
        "only_default": only_default,
        "both": both,
    }


SWEEP = _sweep()


def _r_formals() -> dict[tuple[str, str], str]:
    with FORMALS.open() as handle:
        return {(row["fn"], row["arg"]): row["default"] for row in csv.DictReader(handle)}


R_FORMALS = _r_formals()


# --------------------------------------------------------------------------- #
# 1. The denominator, stated rather than assumed.
# --------------------------------------------------------------------------- #
def test_the_sweep_denominator_is_what_the_report_claims():
    # History of these counts, never / only default / both (all 2026-09-22):
    # - 79 / 258 / 324 over 114 functions (661 defaults) once the HSPELL
    #   whole-matrix goldens call HSPELL_to_STS without status; the four seqdist
    #   generators added earlier that day changed nothing here.
    # - Each audit fix branch added one generator; measured on master after each
    #   merge: s3a + s3b 76 / 249 / 336 (seqmpos, seqidegrad, seqprecarity,
    #   seqici, seqimplic, seqintegr, seqmodst and seqmeant options);
    #   + s5 76 / 248 / 337 (seqecmpgroup's pvalue.limit); + s2 72 / 248 / 341
    #   (seqMD without what, seqstart without tmin/tmax/missing); + s4
    #   72 / 247 / 342 (dissrep's dmax); + s6 72 / 252 / 343 over 116 functions
    #   (667 defaults), since make_goldens_audit20260922_s6.R adds seqpcplot and
    #   seqpcfilter (+2 functions, +6 parameters). The s10, s3c, s2b, s5b, s4b
    #   and s10b generators left the sweep unchanged. Each branch's own deltas
    #   are in its commit messages.
    # 2026-10-08: restored group=None adds one R-NULL default. Generators
    # still pass group explicitly; independent default regressions cover it.
    # 2026-10-12 repair batch 1: make_goldens_seqdef_r_pairing.R,
    # make_goldens_seqecreate_rendering.R and
    # make_goldens_seqdist_weighted_refseq.R and
    # make_goldens_seqidegrad_logical_pow.R (+4 generators) leave the function
    # and default counts unchanged. make_goldens_disstree_assign_names.R (+1)
    # calls disstree without max.depth and disstreeleaf without collapse (never
    # -> both, 73 -> 71), and passes pval to disstree and covar to
    # disstree.assign (only-default -> both, 252 -> 250; both 343 -> 347).
    # make_goldens_seqplotmd_sort_order.R (+1) adds seqplotMD (+1 function, +4
    # parameters: type, cpal, dnames and figsize, all only-default, 250 -> 254).
    # 2026-10-13 repair batch 2: make_goldens_seqdist_twed_auto_indel.R (+1
    # generator) leaves the function and default counts unchanged.
    # make_goldens_dhd_explicit_norms.R (+1 generator) leaves the function and
    # default counts unchanged.
    # make_goldens_seqdist_nan_parameters.R (+1 generator) leaves the function and default counts
    # unchanged.
    # make_goldens_chi2_euclid_weights.R (+1 generator) leaves the function and default counts
    # unchanged.
    # make_goldens_seqedist_empty.R (+1 generator) leaves the function and default counts
    # unchanged.
    # make_goldens_subseqelist_subset.R (+1 generator) leaves the function and default counts
    # unchanged.
    # make_goldens_factor_level_order.R (+1 generator) calls dissmergegroups without small and
    # weights and seqrep.grp without coverage, criterion and weighted, so those five defaults are
    # now taken (never -> both, 71 -> 66; both 347 -> 352).
    # make_goldens_sortv_character_order.R (+1 generator) leaves the function and default counts
    # unchanged.
    # make_goldens_sort_marker_layouts.R (+1 generator) leaves the function and default counts
    # unchanged.
    # make_goldens_seqecreate_label_order.R (+1 generator) leaves the function and default counts
    # unchanged.
    # make_goldens_sts_spell_ids.R (+1 generator) is the first to call STS_to_SPELL (via :::),
    # adding the 4 defaulted parameters of the Python STS_to_SPELL, each seen only at its
    # default: functions 117 -> 118, total 672 -> 676, only_default 254 -> 258.
    # make_goldens_sts_sps_dss_void.R (+1 generator) is the first to call STS_to_TSE (via :::),
    # adding the 2 defaulted parameters of the Python STS_to_TSE, each seen only at its default,
    # and calls seqetm with use.labels at its default: functions 118 -> 119, total 676 -> 678,
    # never 66 -> 65 (seqetm use_labels), only_default 258 -> 260, both 352 -> 353.
    # make_goldens_hspell_alphabet_order.R (+1 generator) leaves the function and default counts
    # unchanged.
    # make_goldens_seqdef_nr_void_neutral.R (+1 generator) leaves the function and default counts
    # unchanged.
    # make_goldens_group_label_rendering.R (+1 generator) calls group.p with weights at its
    # default: never 65 -> 64 (group_p weights), both 353 -> 354.
    # make_goldens_seqstatl_r_numbers.R (+1 generator) leaves the function and default counts
    # unchanged.
    # make_goldens_seqindic_dyn_inputs.R (+1 generator) leaves the function and default counts
    # unchanged.
    # make_goldens_seqpolyads_negative_u.R (+1 generator) calls seqpolyads with T and s set, so
    # one parameter seen only at its default before is now also seen set: only_default 260 -> 259,
    # both 354 -> 355.
    # make_goldens_seqpm_posix.R (+1 generator) calls seqpm with sep at its default: never
    # 64 -> 63 (seqpm sep), both 355 -> 356.
    # make_goldens_internal_void_indicators.R (+1 generator) calls seqST with include_missing, norm
    # and type at their defaults: never 63 -> 60, both 356 -> 359. Its batch-3 review fix also calls
    # seqsubsn with with_missing at its default: never 60 -> 59, both 359 -> 360.
    # make_goldens_seqrecode_palette.R (+1 generator) leaves the function and default counts
    # unchanged.
    # make_goldens_datetime_labels.R (+1 generator) leaves the function and default counts
    # unchanged.
    # make_goldens_disstree_node_labels.R (+1 generator) calls disstree2dot with use_title at its
    # default and image_format set: never 59 -> 58 (disstree2dot use_title), only_default
    # 259 -> 258 (disstree2dot image_format), both 360 -> 362.
    # seqrep gains random_state (TX-AUD-20261010-015), which no generator passes: total
    # 678 -> 679, only_default 258 -> 259.
    # seqplot gains group and main (TX-AUD-20261010-044), so it now has defaulted parameters and
    # enters the sweep; make_goldens_seqplot_group_titles.R (+1 generator) passes type, group and
    # main in every call: functions 119 -> 120, total 679 -> 682, never 58 -> 61 (seqplot type,
    # group, main; all equal R's defaults).
    # make_goldens_dissmfacw_rank.R (+1 generator) calls dissmfacw with squared at its default:
    # never 61 -> 60 (dissmfacw squared), both 362 -> 363.
    # make_goldens_seqplot_rf_stat.R (+1 generator) is the first to call seqplot.rf, so seqplot_rf
    # enters the sweep with its eight defaulted parameters (one call takes k's default): functions
    # 120 -> 121, total 682 -> 690, never 60 -> 61 (seqplot_rf sortv), only_default 259 -> 265,
    # both 363 -> 364.
    # make_goldens_survival_extent.R (+1 generator) is the first to call seqedplot (type always
    # passed) and calls seqsplot with with_legend at its default: functions 121 -> 122, total
    # 690 -> 698, never 61 -> 61 (+ seqedplot type, - seqsplot with_legend), only_default
    # 265 -> 272, both 364 -> 365.
    # make_goldens_seqmeant_rows.R (+1 generator) leaves the function and default counts
    # unchanged.
    # make_goldens_ecmpgroup_colours.R (+1 generator) leaves the function and default counts
    # unchanged.
    # make_goldens_seqpcplot_positions.R (+1 generator) leaves the function and default counts
    # unchanged.
    # make_goldens_index_orientation.R (+1 generator) leaves the function and default counts
    # unchanged.
    # make_goldens_domassoc_window.R (+1 generator) calls seqdomassoc with cross_table, with_missing
    # and dnames set: only_default 272 -> 269, both 365 -> 368.
    # make_goldens_domassoc_unweighted.R (+1 generator) calls seqdomassoc with seqrf_args set:
    # only_default 269 -> 268, both 368 -> 369.
    # make_goldens_domassoc_args.R (+1 generator) leaves the function and default counts
    # unchanged.
    # make_goldens_domassoc_missing_level.R (+1 generator) calls seqdomassoc with struct_zero set:
    # only_default 268 -> 267, both 369 -> 370.
    # make_goldens_erulesdisc_aliased.R (+1 generator) leaves the function and default counts
    # unchanged.
    # make_goldens_dissmergegroups_negative.R (+1 generator) calls dissmergegroups with min_group at
    # its default: never 61 -> 60 (dissmergegroups min_group), both 370 -> 371.
    # make_goldens_disscompare_order.R (+1 generator) leaves the function and default counts
    # unchanged.
    # make_goldens_seqplot_group_rep.R (+1 generator) calls seqplot without main and once without
    # group, and seqrf without sortv, so those three defaults are now taken: never 60 -> 57
    # (seqplot main and group, seqrf sortv), both 371 -> 374.
    # make_goldens_seqnum.R (+1 generator) is the first to call seqnum, with with.missing set and
    # at its default: functions 122 -> 123, total 698 -> 699, both 374 -> 375.
    # make_goldens_seqdim.R (+1 generator) leaves the function and default counts unchanged.
    # make_goldens_mds_orientation.R (+1 generator) leaves the function and default counts
    # unchanged.
    # + make_goldens_seqdist_r_arguments.R (seqdist with every R argument for every R method; TX-
    # AUD-20261010-071).
    # + make_goldens_check_costs.R (seqdist, seqcost and seqMD under check.costs TRUE and FALSE;
    # TX-AUD-20261011-041).
    # + make_goldens_generated_cost_messages.R (checkcost messages on generated costs; TX-
    # AUD-20261011-080).
    # + make_goldens_seqdist_messages.R (seqdist's empty-sequence and with.missing messages; TX-
    # AUD-20261011-042).
    # + make_goldens_degenerate_distances.R (seqdist and seqcost on degenerate inputs; TX-
    # AUD-20261011-038, -039, -040).
    # + make_goldens_refseq_forms.R (seqdist's refseq forms and refusals; TX-AUD-20261011-044).
    # + make_goldens_refseq_asymmetric.R (reference distances with asymmetric costs; TX-
    # REP-20261017-001).
    # + make_goldens_chi2_windows.R (CHI2/EUCLID step and breaks; TX-AUD-20261011-048).
    # + make_goldens_omslen_link.R (OMslen's link check; TX-AUD-20261011-051).
    # + make_goldens_default_alphabet.R (the default alphabet; owner decision D2): the first
    # generator to pass seqMD's ch.sep, so seqMD ch_sep moves from only_default to both:
    # only_default 267 -> 266, both 378 -> 379.
    # + make_goldens_seqpolyads_roles.R (TX-AUD-20261019-004): the first generator to pass
    # seqpolyads' with.missing, so it moves from only_default to both: 266 -> 265, 379 -> 380.
    # + make_goldens_refseq_most_frequent.R (TX-AUD-20261019-001) leaves the function and
    # default counts unchanged.
    # + make_goldens_indels_unobserved.R (TX-AUD-20261019-002) leaves the function and
    # default counts unchanged.
    # + make_goldens_asym_list_matrix.R (TX-AUD-20261019-003) leaves them unchanged too.
    # + make_goldens_asym_bound.R (TX-AUD-20261019-008) leaves them unchanged too.
    # + make_goldens_negative_prox.R (TX-REP-20261019-001) leaves them unchanged too.
    # + make_goldens_spell_powers.R (TX-AUD-20261019-030) leaves them unchanged too.
    # + make_goldens_subseq_overflow.R (TX-AUD-20261019-029) leaves them unchanged too.
    # + make_goldens_norm_bounds.R (TX-REP-20261021-003 to -005) leaves them unchanged too.
    # + make_goldens_spelling_convention.R (TX-AUD-20261019-061) leaves them unchanged too.
    # + make_goldens_event_arguments.R (batch 11).
    # + make_goldens_rank_and_convergence.R (TX-REP-20261019-002) leaves them unchanged too.
    # + make_goldens_diss_trees_reps.R (batch 12).
    # + make_goldens_plot_items.R (batch 13); its index and chi-squared rows move no parameter
    # between classes, its seqedplot rows (TX-AUD-20261019-050) do (below), and its seqfplot rows
    # (TX-AUD-20261011-067) are the first to call seqfplot: functions 123 -> 124; its tick-centre
    # rows the first to call seqmsplot, seqrplot and seqrfplot: 124 -> 127.
    # + make_goldens_construction_items.R (batch 14).
    # + make_goldens_statistics_items.R (batch 15).
    # + make_goldens_ledger_items.R (batch 18).
    # + make_goldens_disstree_weighted_gate.R and make_goldens_consumed_items.R (2026-11-02).
    # + make_goldens_study_findings.R (2026-11-07).
    assert SWEEP["generators"] == 209
    assert SWEEP["functions"] == 127
    # make_goldens_check_costs.R exercises check.costs = FALSE of seqcost and seqMD:
    # total 699 -> 701, both 375 -> 377.
    # dissvar_grp's group gained R's default NULL (None; TX-AUD-20261011-004), a new pair, which 4 of
    # its 7 calls pass: total 701 -> 702, both +1.
    # seqfplot's 8 parameters (batch 13, make_goldens_plot_items.R): total 702 -> 710; seqmsplot's 8,
    # seqrfplot's 14 and seqrplot's 11: 710 -> 743.
    assert SWEEP["total"] == 743
    # make_goldens_degenerate_distances.R calls seqtrate: never 57 -> 56, both 377 -> 378.
    # make_goldens_diss_trees_reps.R calls seqrep.grp with mdis and without diss or ret: seqrep_grp's
    # diss and ret leave never (56 -> 54) for both, and its mdis leaves only_default for both.
    # Its dissrf calls (TX-AUD-20261011-023) are the first to take grp.meth's default: never 54 -> 53;
    # its dissrf scale rows (TX-AUD-20261011-052) the first to take sortv's: never 53 -> 52; its
    # seqtree rows (TX-AUD-20261010-080) the first to take diss, max.depth, min.size and
    # weight.permutation's: never 52 -> 48; its dissindic rows (TX-AUD-20261011-066) the first to take
    # squared and weights': never 48 -> 46.
    # make_goldens_plot_items.R's seqedplot rows (TX-AUD-20261019-050) are the first to take
    # seqedplot's type default: never 46 -> 45.
    # make_goldens_statistics_items.R's seqgen rows (TX-AUD-20261011-001) are the first to take seqgen's p
    # default: never 45 -> 44. Its seqemlt rows (TX-REP-20261023-005) are the first to take seqemlt's weighted
    # default: never 44 -> 43.
    # make_goldens_ledger_items.R's gower_matrix rows (TX-AUD-20261011-005) are the first to take
    # gower_matrix's squared default: never 43 -> 42.
    # 2026-11-07: make_goldens_study_findings.R's dissrep row is the first to take dissrep's criterion
    # default: never 42 -> 41.
    assert len(SWEEP["never"]) == 41
    # make_goldens_audit20260922_s3b.R is the first generator to pass seqipos' with.missing
    # (TX-AUD-20261019-012): only_default 265 -> 264, both 380 -> 381.
    # make_goldens_sts_spell_ids.R is the first generator to pass STS_to_SPELL's with.missing
    # (TX-AUD-20261019-014): only_default 264 -> 263, both 381 -> 382.
    # make_goldens_spelling_convention.R's seqMD norm rows are the first to pass seqMD's norm
    # (0 of 58 seqMD calls before, 8 of 66 after; follow-up review I-1): only_default 263 -> 262,
    # both 382 -> 383.
    # Batch 12 (make_goldens_diss_trees_reps.R and dissvar_grp's default): only_default 262 -> 261,
    # both 383 -> 387 (seqrep_grp diss, ret, mdis; dissvar_grp group), then 387 -> 389 (dissrf grp_meth
    # and sortv), then 389 -> 394 (seqtree's four, and seqdist_args, which its seqdist.args = list()
    # call is the first to pass: only_default 261 -> 260), then 394 -> 396 (dissindic squared and weights).
    # Batch 13: make_goldens_plot_items.R's seqedplot rows are the first to pass seqedplot's cpal and
    # group (only_default 260 -> 258) and to take its type default: both 396 -> 399. Its seqfplot
    # calls pass none of seqfplot's parameters by their Python spellings: only_default 258 -> 266;
    # its seqmsplot, seqrfplot and seqrplot calls take all 33 of their parameters' defaults (seqrplot's
    # distances computed from method, seqrfplot's default k): 266 -> 299.
    # 2026-11-02: make_goldens_disstree_weighted_gate.R is the first to pass disstree's first and squared
    # and seqtree's first, pval and squared: 299 -> 294.
    assert SWEEP["only_default"] == 294
    # Batch 15: make_goldens_statistics_items.R's seqgen rows take p's default, which the other seqgen calls
    # pass: (seqgen, p) leaves never for both, 399 -> 400; its seqemlt rows take weighted's default, which other
    # seqemlt calls pass: (seqemlt, weighted) leaves never for both, 400 -> 401.
    # Batch 18: make_goldens_ledger_items.R's gower_matrix rows take squared's default, which
    # make_goldens_gower.R passes: (gower_matrix, squared) leaves never for both, 401 -> 402.
    # 2026-11-02: those five of disstree and seqtree, which other calls leave at their defaults: 402 -> 407.
    # 2026-11-07: make_goldens_study_findings.R's dissrep row takes criterion's default, which the other
    # dissrep calls pass: (dissrep, criterion) leaves never for both, 407 -> 408.
    assert SWEEP["both"] == 408
    assert SWEEP["total"] == len(SWEEP["never"]) + SWEEP["only_default"] + SWEEP["both"]


#: Every parameter whose default no generator ever takes, frozen so that a new
#: one has to be looked at rather than accumulating silently.
NEVER_EXERCISED = {
    ("FCE_to_TSE", "cols"),
    ("HSPELL_to_STS", "begin"),
    ("HSPELL_to_STS", "end"),
    ("dissCompare", "s"),
    ("dissassoc", "R"),
    ("dissmergegroups", "silent"),
    ("dissrf", "k"),
    ("disstree", "R"),
    ("disstree", "min_size"),
    ("disstree2dot", "image_fun"),
    ("disstreeleaf", "label"),
    ("gower_matrix", "weights"),
    ("seqBIC", "group"),
    ("seqBIC", "method"),
    ("seqBIC", "s"),
    ("seqBIC", "weighted"),
    ("seqCompare", "method"),
    ("seqCompare", "s"),
    ("seqLRT", "group"),
    ("seqLRT", "method"),
    ("seqLRT", "s"),
    ("seqLRT", "weighted"),
    ("seqalign", "include_missing"),
    ("seqalign", "indel"),
    ("seqalign", "sm"),
    ("seqauto", "order"),
    ("seqcta", "subseq"),
    ("seqerulesdisc", "adjust"),
    ("seqerulesdisc", "pvalue"),
    ("seqetm", "method"),
    ("seqgen_missing", "p_cases"),
    ("seqgranularity", "tspan"),
    ("seqindic", "indic"),
    ("seqplot", "type"),
    ("seqplot_rf", "sortv"),
    ("seqrep_grp", "group"),
    ("seqsep", "sep"),
    ("seqsep", "sl"),
    ("seqstatf", "include_missing"),
    ("seqstatf", "weighted"),
    ("seqtree", "R"),
}


def test_the_never_exercised_set_has_not_drifted():
    assert set(SWEEP["never"]) == NEVER_EXERCISED


# --------------------------------------------------------------------------- #
# 2. For those, the only remaining check: does traminex's default match R's?
# --------------------------------------------------------------------------- #
def _normalise(value: str) -> str:
    text = str(value).strip()
    if text.startswith("c(") and text.endswith(")"):
        items = [part.strip().strip('"').strip("'") for part in text[2:-1].split(",")]
        return "|".join(items)
    if text.startswith("(") and text.endswith(")"):
        items = [
            part.strip().strip('"').strip("'") for part in text[1:-1].split(",") if part.strip()
        ]
        return "|".join(items)
    text = text.strip('"').strip("'")
    lowered = text.lower()
    mapping = {"true": "True", "false": "False", "none": "None", "null": "None", "na": "None"}
    if lowered in mapping:
        return mapping[lowered]
    try:
        number = float(text)
    except ValueError:
        return text
    return repr(int(number)) if number == int(number) else repr(number)


#: The never-exercised defaults for which R has no comparable formal: either R
#: requires the argument and traminex supplies a convenience default, or the
#: parameter is a traminex extension. Each is covered by a test instead.
NO_R_COUNTERPART = {
    ("FCE_to_TSE", "cols"),
    ("HSPELL_to_STS", "begin"),
    ("HSPELL_to_STS", "end"),
    ("seqBIC", "method"),
    ("seqCompare", "method"),
    ("seqLRT", "method"),
    ("seqalign", "sm"),
}

#: Never-exercised defaults where traminex deliberately differs from R's.
KNOWN_DEFAULT_DIVERGENCES = {
    # R prints each merge step (silent=FALSE); a library should be quiet by
    # default. Values are identical either way -- proven in
    # test_dissmergegroups_default_is_quiet_but_otherwise_identical.
    ("dissmergegroups", "silent"): ("True", "FALSE"),
}


def _r_argument(r_fn: str, parameter: str) -> tuple[str, str] | None:
    for candidate in sorted(_spellings(parameter)):
        if (r_fn, candidate) in R_FORMALS:
            return candidate, R_FORMALS[(r_fn, candidate)]
    return None


@pytest.mark.parametrize("python_name,parameter", sorted(NEVER_EXERCISED), ids=lambda v: str(v))
def test_a_never_exercised_default_matches_r_or_is_a_recorded_divergence(python_name, parameter):
    """No golden can vouch for these, so R's own formals must."""
    if (python_name, parameter) in NO_R_COUNTERPART:
        pytest.skip("R has no comparable formal; covered by a dedicated test")

    r_name = SWEEP["r_names"][python_name]
    found = _r_argument(r_name, parameter)
    if found is None:
        pytest.skip(f"{r_name} is not an exported R name with this argument")
    r_argument, r_default = found
    if r_default == "<no default>":
        pytest.skip("R requires this argument")

    python_default = _normalise(repr(SWEEP["never"][(python_name, parameter)]))
    recorded = KNOWN_DEFAULT_DIVERGENCES.get((python_name, parameter))
    if recorded is not None:
        assert (python_default, r_default) == recorded
        return

    assert python_default == _normalise(r_default), (
        f"{python_name}.{parameter} defaults to {python_default!r} but R's "
        f"{r_argument} defaults to {r_default!r}, and no generator exercises it"
    )


def test_every_never_exercised_default_is_looked_up_under_the_r_name_the_sweep_read():
    """R's formals are keyed by R's name, e.g. ``seqrep.grp``, not ``seqrep_grp``.

    The lookup used the Python name, so every function whose R name contains a
    dot was skipped as "not an exported R name" and its defaults were never
    compared with R's (2026-09-22 audit S10-02).
    """
    unresolved = sorted(
        (python_name, parameter)
        for python_name, parameter in NEVER_EXERCISED - NO_R_COUNTERPART
        if _r_argument(SWEEP["r_names"][python_name], parameter) is None
    )
    assert unresolved == []


def test_exactly_one_never_exercised_default_diverges_from_r():
    """The headline number, so a second divergence cannot slip in unnoticed."""
    assert len(KNOWN_DEFAULT_DIVERGENCES) == 1


# --------------------------------------------------------------------------- #
# 3. Coverage for the never-exercised defaults that matter: the one that
#    diverges from R, and the ones traminex invented where R requires the
#    argument. Nothing else can vouch for these.
# --------------------------------------------------------------------------- #
def _merging_corpus():
    """Six labelled groups over three real clusters, so merges actually occur.

    ``min_group`` defaults to 4, so the merge loop only runs with more than
    four groups -- an earlier three-group attempt merged nothing and would have
    proved nothing about ``silent``.
    """
    import numpy as np

    rng = np.random.default_rng(11)
    centres = np.array([0.0, 0.1, 5.0, 5.1, 10.0, 10.1])
    labels = np.repeat(np.arange(6), 5)
    points = centres[labels] + rng.normal(scale=0.15, size=labels.size)
    return np.abs(points[:, None] - points[None, :]), labels


def test_dissmergegroups_is_quiet_by_default_but_otherwise_identical():
    """The one never-exercised default that differs from R's.

    ``dissmergegroups.R:77`` defaults ``silent`` to FALSE and prints a line per
    merge; traminex defaults it to True. Only the printed progress differs.
    """
    import contextlib
    import io

    import numpy as np

    matrix, labels = _merging_corpus()

    quiet_output = io.StringIO()
    with contextlib.redirect_stdout(quiet_output):
        quiet = tx.dissmergegroups(matrix, labels.copy())
    loud_output = io.StringIO()
    with contextlib.redirect_stdout(loud_output):
        loud = tx.dissmergegroups(matrix, labels.copy(), silent=False)

    assert quiet_output.getvalue() == ""
    assert "Merging groups" in loud_output.getvalue()
    # ...and the merge really happened, so the comparison is not vacuous
    assert len(np.unique(quiet["group"])) < 6

    np.testing.assert_array_equal(quiet["group"], loud["group"])
    assert quiet["quality"] == loud["quality"]


@pytest.mark.parametrize("entry", ["seqCompare", "seqLRT", "seqBIC"])
def test_the_invented_lcs_method_default_equals_passing_lcs(entry):
    """R gives ``method`` no default; traminex defaults it to "LCS".

    Every generator passes it, so nothing pinned the default until now.
    """
    import numpy as np
    import pandas as pd

    sequences = tx.seqdef(
        pd.DataFrame(
            [
                list("ABAB"),
                list("BABA"),
                list("AABB"),
                list("BBAA"),
                list("ABBA"),
                list("BAAB"),
            ]
        )
    )
    group = ["g1", "g2", "g1", "g2", "g1", "g2"]
    caller = {"seqCompare": tx.seqCompare, "seqLRT": tx.seqLRT, "seqBIC": tx.seqBIC}[entry]

    implicit = caller(sequences, group=group, s=0)
    explicit = caller(sequences, group=group, s=0, method="LCS")

    np.testing.assert_allclose(
        implicit.to_numpy(dtype=float), explicit.to_numpy(dtype=float), rtol=0, atol=0
    )
    assert list(implicit.columns) == list(explicit.columns)


def test_the_invented_empty_cols_default_refuses_with_rs_own_wording():
    """R makes ``cols`` mandatory; traminex defaults it to ``()``.

    An omitted ``cols`` therefore behaves like R's ``cols = integer(0)``, whose
    message is reproduced exactly, rather than R's "argument "cols" is missing,
    with no default", which a defaulted parameter cannot produce.
    """
    import pandas as pd

    frame = pd.DataFrame({"id": [1, 2], "a": [1.0, 2.0], "b": [2.0, 3.0]})

    with pytest.raises(ValueError, match="argument is of length zero"):
        tx.FCE_to_TSE(frame)


def test_the_invented_hspell_column_defaults_refuse_naming_the_columns():
    """R requires the spell columns; traminex defaults them to None."""
    import pandas as pd

    frame = pd.DataFrame({"id": [1], "b": [1], "e": [2], "s": ["A"]})

    with pytest.raises(ValueError, match="start_col.*end_col.*state_col"):
        tx.HSPELL_to_STS(frame)


def test_the_invented_seqalign_cost_defaults_produce_a_unit_cost_alignment():
    """R makes ``sm`` mandatory; traminex defaults it to None and ``indel`` to 1.0."""
    import numpy as np
    import pandas as pd

    sequences = tx.seqdef(pd.DataFrame([list("ABAB"), list("BABA")]))

    implicit = tx.seqalign(sequences, (0, 1))
    explicit = tx.seqalign(sequences, (0, 1), indel=1.0, sm=None)

    assert implicit.total_cost == explicit.total_cost
    assert np.isfinite(implicit.total_cost)
    np.testing.assert_array_equal(implicit.cost_matrix, explicit.cost_matrix)
