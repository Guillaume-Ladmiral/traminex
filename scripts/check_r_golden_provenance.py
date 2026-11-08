#!/usr/bin/env python3
"""Validate and, optionally, execute the complete R-golden provenance contract."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

# 2026-10-12 repair batch 1: + make_goldens_seqdef_r_pairing.R (TX-AUD-20261010-001),
# make_goldens_seqecreate_rendering.R (-002/-016), make_goldens_seqdist_weighted_refseq.R (-003),
# make_goldens_seqidegrad_logical_pow.R (-004), make_goldens_disstree_assign_names.R (-006),
# make_goldens_seqplotmd_sort_order.R (-008).
# 2026-10-13 repair batch 2: + make_goldens_seqdist_twed_auto_indel.R (TX-AUD-20261010-011).
# + make_goldens_dhd_explicit_norms.R.
# + make_goldens_seqdist_nan_parameters.R.
# + make_goldens_chi2_euclid_weights.R.
# + make_goldens_seqedist_empty.R.
# + make_goldens_subseqelist_subset.R.
# + make_goldens_factor_level_order.R.
# + make_goldens_sortv_character_order.R.
# + make_goldens_sort_marker_layouts.R.
# + make_goldens_seqecreate_label_order.R.
# + make_goldens_sts_spell_ids.R.
# + make_goldens_sts_sps_dss_void.R.
# + make_goldens_hspell_alphabet_order.R.
# + make_goldens_seqdef_nr_void_neutral.R.
# + make_goldens_group_label_rendering.R.
# + make_goldens_seqstatl_r_numbers.R.
# + make_goldens_seqindic_dyn_inputs.R.
# + make_goldens_seqpolyads_negative_u.R.
# + make_goldens_seqpm_posix.R.
# + make_goldens_internal_void_indicators.R.
# + make_goldens_seqrecode_palette.R.
# + scripts/r_goldens/make_goldens_datetime_labels.R.
# + scripts/r_goldens/make_goldens_disstree_node_labels.R.
# + scripts/r_goldens/make_goldens_seqplot_group_titles.R.
# + scripts/r_goldens/make_goldens_dissmfacw_rank.R.
# + scripts/r_goldens/make_goldens_seqplot_rf_stat.R.
# + scripts/r_goldens/make_goldens_survival_extent.R.
# + scripts/r_goldens/make_goldens_seqmeant_rows.R.
# + scripts/r_goldens/make_goldens_ecmpgroup_colours.R.
# + scripts/r_goldens/make_goldens_seqpcplot_positions.R.
# + scripts/r_goldens/make_goldens_index_orientation.R.
# + scripts/r_goldens/make_goldens_domassoc_window.R.
# + scripts/r_goldens/make_goldens_domassoc_unweighted.R.
# + scripts/r_goldens/make_goldens_domassoc_args.R.
# + scripts/r_goldens/make_goldens_domassoc_missing_level.R.
# + scripts/r_goldens/make_goldens_erulesdisc_aliased.R.
# + scripts/r_goldens/make_goldens_dissmergegroups_negative.R.
# + scripts/r_goldens/make_goldens_disscompare_order.R.
# + scripts/r_goldens/make_goldens_seqplot_group_rep.R.
# + scripts/r_goldens/make_goldens_seqnum.R.
# + scripts/r_goldens/make_goldens_seqdim.R.
# + scripts/r_goldens/make_goldens_mds_orientation.R.
# + make_goldens_seqdist_r_arguments.R.
# + make_goldens_check_costs.R.
# + make_goldens_generated_cost_messages.R.
# + make_goldens_seqdist_messages.R.
# + make_goldens_degenerate_distances.R.
# + make_goldens_refseq_forms.R.
# + make_goldens_refseq_asymmetric.R.
# + make_goldens_chi2_windows.R.
# + make_goldens_omslen_link.R.
# + make_goldens_default_alphabet.R (owner decision D2).
# + make_goldens_seqpolyads_roles.R (TX-AUD-20261019-004).
# + make_goldens_refseq_most_frequent.R (TX-AUD-20261019-001).
# + make_goldens_indels_unobserved.R (TX-AUD-20261019-002).
# + make_goldens_asym_list_matrix.R (TX-AUD-20261019-003).
# + make_goldens_asym_bound.R (TX-AUD-20261019-008).
# + make_goldens_seqdist_check_order.R (TX-REP-20261019-010).
# + make_goldens_seqdiff_messages.R (TX-AUD-20261019-059).
# + make_goldens_seqmd_details.R (TX-AUD-20261010-074).
# + make_goldens_degenerate_omstran.R (TX-AUD-20261019-080).
# + make_goldens_negative_prox.R (TX-REP-20261019-001).
# + make_goldens_spell_powers.R (TX-AUD-20261019-030).
# + make_goldens_subseq_overflow.R (TX-AUD-20261019-029, TX-AUD-20261011-047).
# + make_goldens_norm_bounds.R (TX-REP-20261021-003 to -005).
# + make_goldens_spelling_convention.R (TX-AUD-20261019-061, owner decision of 2026-10-22).
# + make_goldens_event_arguments.R (batch 11: TX-AUD-20261019-043 and the event items).
# + make_goldens_rank_and_convergence.R (batch 11: TX-REP-20261019-002).
# + make_goldens_diss_trees_reps.R (batch 12: TX-AUD-20261019-046 and the dissimilarity-analysis items).
# + make_goldens_plot_items.R (batch 13: TX-AUD-20261019-048 and the plotting items).
# + make_goldens_construction_items.R (batch 14: sequence construction, formats and multidomain).
# + make_goldens_statistics_items.R (batch 15: sequence statistics).
# + make_goldens_ledger_items.R (batch 18: TX-AUD-20261011-005, the ledger items).
# + make_goldens_disstree_weighted_gate.R and make_goldens_consumed_items.R (2026-11-02: what the
#   other libraries use and no test here ran).
# + make_goldens_study_findings.R (2026-11-07: the 2026-10-25 study's findings).
EXPECTED_SCRIPT_COUNT = 209
# 612 + 5 cycle-6 F-001 seqmaintokens default-k fixtures (make_goldens_precarity.R)
# + 4 cycle-6 F-003 case-weight edge-case fixtures (make_goldens_dissvar.R,
# make_goldens_disscenter.R) + 5 cycle-6 section-9 seqalign DP-matrix fixtures
# (make_goldens_seqalign.R) + 5 cycle-7 W-01 deep/degenerate event-engine
# fixtures (make_goldens_event_engine_edges.R) + 12 cycle-7 W-02 deterministic
# group-comparison fixtures (make_goldens_group_compare_deterministic.R)
# + 9 cycle-7 W-05 seqdist Rd-option fixtures
# (make_goldens_seqdist_rd_options.R) + 1 cycle-7 W-07 R-formals fixture
# (make_goldens_r_formals.R) + 15 cycle-7 W-04 wide event and
# group-comparison fixtures (make_goldens_events_wide_corpus.R) + 4 seqdist
# with.missing-reset fixtures (make_goldens_seqdist_with_missing_reset.R)
# + 1 HSPELL whole-matrix fixture (make_goldens_hspell_horizontal.R)
# + 2 seqdist state-space fixtures (make_goldens_seqdist_state_space.R)
# + 2 seqdist cost-check fixtures (make_goldens_seqdist_cost_checks.R)
# + 2 asymmetric-duplicate fixtures (make_goldens_seqdist_asymmetric_duplicates.R)
# + 2 kweights fixtures (make_goldens_seqdist_kweights.R)
# + 1 numeric-edge fixture (make_goldens_seqdist_numeric_edges.R)
# + 1 argument-form fixture (make_goldens_seqdist_argument_forms.R)
# + 1 2026-09-22 audit s3a statistics fixture (make_goldens_audit20260922_s3a.R).
# + 1 audit 2026-09-22 s3b degenerate and boundary statistics fixture
# (make_goldens_audit20260922_s3b.R).
# + 4 2026-09-22 audit S5 event and SAMM fixtures (make_goldens_audit20260922_s5.R).
# + 1 audit-20260922 S2 construction/format fixture (make_goldens_audit20260922_s2.R).
# + 1 checkcost diagnostic fixture (make_goldens_audit20260922_s10.R)
# + 1 inline-R-anchor fixture (make_goldens_audit20260922_s10.R).
# + 1 audit 2026-09-22 s3c follow-up statistics fixture
# (make_goldens_audit20260922_s3c.R).
# + 11 2026-09-22 audit s4 discrepancy-analysis fixtures (make_goldens_audit20260922_s4.R).
# + 1 audit-20260922 S2-19 seqMD void-component cost fixture (make_goldens_audit20260922_s2b.R).
# + 1 2026-09-22 audit S5-16 chi-square arithmetic fixture (make_goldens_audit20260922_s5b.R).
# + 4 2026-09-22 audit S6 plotting-semantics fixtures (make_goldens_audit20260922_s6.R)
# + 2 2026-09-22 audit s4b disstree/dissrep follow-up fixtures (make_goldens_audit20260922_s4b.R).
# + 2 2026-09-22 audit S10-11/S10-12 seqMD cost fixtures (make_goldens_audit20260922_s10b.R).
# + 3 seqdef state/label/colour pairing fixtures (make_goldens_seqdef_r_pairing.R).
# + 1 seqecreate event-rendering fixture (make_goldens_seqecreate_rendering.R).
# + 1 weighted CHI2/OMstran refseq fixture (make_goldens_seqdist_weighted_refseq.R).
# + 1 seqidegrad logical-pow fixture (make_goldens_seqidegrad_logical_pow.R).
# + 2 disstree.assign covariate-name fixtures (make_goldens_disstree_assign_names.R).
# + 1 seqplotMD drawn-order fixture (make_goldens_seqplotmd_sort_order.R).
# + 1 TWED automatic-indel fixture (make_goldens_seqdist_twed_auto_indel.R).
# + 2 DHD normalisation and cost-refusal fixtures (make_goldens_dhd_explicit_norms.R).
# + 1 seqdist NaN-parameter fixture (make_goldens_seqdist_nan_parameters.R).
# + 1 CHI2/EUCLID weights fixture (make_goldens_chi2_euclid_weights.R).
# + 1 seqedist empty-sequence fixture (make_goldens_seqedist_empty.R).
# + 1 subseqelist subset fixture (make_goldens_subseqelist_subset.R).
# + 2 text group level-order fixtures (make_goldens_factor_level_order.R).
# + 1 sortv order() fixture (make_goldens_sortv_character_order.R).
# + 4 sort-with-markers and seqrf fixtures (make_goldens_sort_marker_layouts.R).
# + 1 event label order fixture (make_goldens_seqecreate_label_order.R).
# + 1 STS_to_SPELL id fixture (make_goldens_sts_spell_ids.R).
# + 1 SPS/DSS void and TSE id fixture (make_goldens_sts_sps_dss_void.R).
# + 1 HSPELL_to_STS alphabet fixture (make_goldens_hspell_alphabet_order.R).
# + 1 seqdef nr/void and NEUTRAL fixture (make_goldens_seqdef_nr_void_neutral.R).
# + 1 group label rendering fixture (make_goldens_group_label_rendering.R).
# + 1 seqstatl R numbers fixture (make_goldens_seqstatl_r_numbers.R).
# + 1 seqindic.dyn indic fixture (make_goldens_seqindic_dyn_inputs.R).
# + 1 seqpolyads negative U fixture (make_goldens_seqpolyads_negative_u.R).
# + 1 seqpm POSIX class fixture (make_goldens_seqpm_posix.R).
# + 1 internal void indicator fixture (make_goldens_internal_void_indicators.R).
# + 1 seqrecode palette fixture (make_goldens_seqrecode_palette.R).
# + 1 seqrecode warning fixture (make_goldens_seqrecode_palette.R).
# + 1 Live-R goldens for date and date-time time labels fixture (scripts/r_goldens/make_goldens_datetime_labels.R).
# + 1 Live-R goldens for discrepancy-tree node text fixture (scripts/r_goldens/make_goldens_disstree_node_labels.R).
# + 1 Live-R goldens for seqplot's group panel titles fixture (scripts/r_goldens/make_goldens_seqplot_group_titles.R).
# + 3 Live-R goldens for dissmfacw on rank-deficient designs fixtures (scripts/r_goldens/make_goldens_dissmfacw_rank.R).
# + 1 Live-R goldens for seqplot.rf's printed statistic fixture (scripts/r_goldens/make_goldens_seqplot_rf_stat.R).
# + 2 Live-R goldens for how far survival curves are drawn fixtures (scripts/r_goldens/make_goldens_survival_extent.R).
# + 1 Live-R goldens for seqmtplot's bars fixture (scripts/r_goldens/make_goldens_seqmeant_rows.R).
# + 2 Live-R goldens for plot.subseqelistchisq's bar colours fixtures (scripts/r_goldens/make_goldens_ecmpgroup_colours.R).
# + 2 Live-R goldens for seqpcplot's event positions fixtures (scripts/r_goldens/make_goldens_seqpcplot_positions.R).
# + 1 Live-R goldens for the index plot's row order fixture (scripts/r_goldens/make_goldens_index_orientation.R).
# + 1 Live-R goldens for seqdomassoc's windowed tables fixture (scripts/r_goldens/make_goldens_domassoc_window.R).
# + 2 Live-R goldens for seqdomassoc weighted=FALSE fixtures (scripts/r_goldens/make_goldens_domassoc_unweighted.R).
# + 1 Live-R goldens for seqdomassoc's argument lists fixture (scripts/r_goldens/make_goldens_domassoc_args.R).
# + 1 Live-R goldens for seqdomassoc's missing level fixture (scripts/r_goldens/make_goldens_domassoc_missing_level.R).
# + 1 Live-R goldens for seqerulesdisc with an aliased predictor fixture (scripts/r_goldens/make_goldens_erulesdisc_aliased.R).
# + 2 Live-R goldens for dissmergegroups without a merge fixtures (scripts/r_goldens/make_goldens_dissmergegroups_negative.R).
# + 2 Live-R goldens for dissCompare's group order fixtures (scripts/r_goldens/make_goldens_disscompare_order.R).
# + 3 Live-R goldens for seqplot's grouped r and rf panels fixtures (scripts/r_goldens/make_goldens_seqplot_group_rep.R).
# + 1 Live-R goldens for seqnum fixture (scripts/r_goldens/make_goldens_seqnum.R).
# + 1 Live-R goldens for seqdim fixture (scripts/r_goldens/make_goldens_seqdim.R).
# + 1 Live-R goldens for the MDS axis orientation fixture (scripts/r_goldens/make_goldens_mds_orientation.R).
# + 1 seqdist R-argument surface fixture (make_goldens_seqdist_r_arguments.R).
# + 1 check.costs fixture (make_goldens_check_costs.R).
# + 1 generated-cost message fixture (make_goldens_generated_cost_messages.R).
# + 1 seqdist message fixture (make_goldens_seqdist_messages.R).
# + 1 degenerate-input distance fixture (make_goldens_degenerate_distances.R).
# + 1 refseq form fixture (make_goldens_refseq_forms.R).
# + 1 asymmetric refseq fixture (make_goldens_refseq_asymmetric.R).
# + 1 CHI2 window fixture (make_goldens_chi2_windows.R).
# + 1 OMslen link fixture (make_goldens_omslen_link.R).
# + 2 default-alphabet fixtures (make_goldens_default_alphabet.R) and 1 TSE_to_STS
#   alphabet fixture (make_goldens_TSE_to_STS.R) (owner decision D2).
# + 1 seqpolyads role fixture (make_goldens_seqpolyads_roles.R).
# + 1 refseq=0 fixture (make_goldens_refseq_most_frequent.R).
# + 1 INDELS unobserved-state fixture (make_goldens_indels_unobserved.R).
# + 1 nested-list matrix fixture (make_goldens_asym_list_matrix.R).
# + 1 seqindic row-name fixture (make_goldens_seqindic.R; TX-AUD-20261010-093).
# + 1 asymmetric-bound fixture (make_goldens_asym_bound.R; TX-AUD-20261019-008).
# + 1 seqentrans fixture (make_goldens_seqentrans.R; TX-AUD-20261019-016).
# + 1 argument-check fixture (make_goldens_seqdist_check_order.R; TX-REP-20261019-010).
# + 1 seqdiff-conditions fixture (make_goldens_seqdiff_messages.R; TX-AUD-20261019-059).
# + 1 seqMD-details fixture (make_goldens_seqmd_details.R; TX-AUD-20261010-074).
# + 1 degenerate-OMstran fixture (make_goldens_degenerate_omstran.R; TX-AUD-20261019-080).
# + 1 negative-proximity fixture (make_goldens_negative_prox.R; TX-REP-20261019-001).
# + 1 spell-power fixture (make_goldens_spell_powers.R; TX-AUD-20261019-030).
# + 1 subsequence-overflow fixture (make_goldens_subseq_overflow.R; TX-AUD-20261019-029).
# + 1 normalization-bound fixture (make_goldens_norm_bounds.R; TX-REP-20261021-003 to -005).
# + 1 spelling-convention fixture (make_goldens_spelling_convention.R; TX-AUD-20261019-061).
# + 1 event-arguments fixture (make_goldens_event_arguments.R; batch 11).
# + 1 rank-and-convergence fixture (make_goldens_rank_and_convergence.R; TX-REP-20261019-002).
# + 1 dissimilarity trees and representatives fixture (make_goldens_diss_trees_reps.R; batch 12).
# + 1 plotting-items fixture (make_goldens_plot_items.R; batch 13).
# + 1 construction-items fixture (make_goldens_construction_items.R; batch 14).
# + 1 statistics-items fixture (make_goldens_statistics_items.R; batch 15).
# + 1 ledger-items fixture (make_goldens_ledger_items.R; batch 18).
# + 1 weighted-gate fixture (make_goldens_disstree_weighted_gate.R) and 1 consumed-items fixture
#   (make_goldens_consumed_items.R; 2026-11-02).
# + 1 study-findings fixture (make_goldens_study_findings.R; 2026-11-07).
EXPECTED_FIXTURE_COUNT = 819
HELPER_SCRIPTS = frozenset({"make_goldens_corestats.R"})
DIAGNOSTIC_SCRIPTS = frozenset(
    {
        "make_goldens_distances.R",
        "make_goldens_norm.R",
        "make_goldens_norm_costs.R",
    }
)

# Every generator that source()s make_goldens_md_inputs.R or make_goldens_phase2_inputs.R
# rewrites their CSVs; tests/test_r_golden_provenance.py derives both sets from those calls
# (TX-AUD-20261019-023 added the domassoc_* and disstree_node_labels generators of batch 4).
_MD_SHARED_WRITERS = frozenset(
    {
        "make_goldens_dissdomassoc.R",
        "make_goldens_domassoc_args.R",
        "make_goldens_domassoc_unweighted.R",
        "make_goldens_domassoc_window.R",
        "make_goldens_md_inputs.R",
        "make_goldens_seqdomassoc.R",
    }
)
_PHASE2_SHARED_WRITERS = frozenset(
    {
        "make_goldens_diss_mfa.R",
        "make_goldens_diss_rf.R",
        "make_goldens_dissassoc.R",
        "make_goldens_disscenter.R",
        "make_goldens_dissdomassoc.R",
        "make_goldens_dissindic.R",
        "make_goldens_dissmergegroups.R",
        "make_goldens_dissrep.R",
        "make_goldens_disstree.R",
        "make_goldens_disstree_dot.R",
        "make_goldens_disstree_full.R",
        "make_goldens_disstree_node_labels.R",
        "make_goldens_dissvar.R",
        "make_goldens_dissvar_grp.R",
        "make_goldens_domassoc_args.R",
        "make_goldens_domassoc_unweighted.R",
        "make_goldens_domassoc_window.R",
        "make_goldens_gower.R",
        "make_goldens_group_compare.R",
        "make_goldens_md_inputs.R",
        "make_goldens_phase2_inputs.R",
        "make_goldens_seqdiff.R",
        "make_goldens_seqdiff_void.R",
        "make_goldens_seqdomassoc.R",
        "make_goldens_seqrep_grp.R",
    }
)
EXPECTED_SHARED_OWNERS: dict[str, frozenset[str]] = {
    **{
        fixture: _MD_SHARED_WRITERS
        for fixture in (
            "md_input_d1.csv",
            "md_input_d2.csv",
            "md_input_d3.csv",
            "md_input_s1.csv",
            "md_input_s2.csv",
            "md_input_s3.csv",
            "md_input_weights.csv",
        )
    },
    **{
        fixture: _PHASE2_SHARED_WRITERS
        for fixture in (
            "phase2_input_D.csv",
            "phase2_input_weights_int.csv",
            "phase2_input_weights_real.csv",
        )
    },
}

_QUOTED_PATH_RE = re.compile(r'"((?:[^"\\]|\\.)*)"')
_WRITE_SYSCALL_MARKERS = ("O_WRONLY", "O_RDWR", "O_CREAT", "O_TRUNC", "rename(")


@dataclass(frozen=True)
class Inventory:
    """The complete classified script and fixture inventory."""

    scripts: tuple[Path, ...]
    fixtures: tuple[Path, ...]
    writers: tuple[Path, ...]
    helpers: tuple[Path, ...]
    diagnostics: tuple[Path, ...]


class ProvenanceError(RuntimeError):
    """Raised when the R-golden provenance contract is violated."""


def _repo_root(candidate: Path | None = None) -> Path:
    if candidate is not None:
        return candidate.resolve()
    return Path(__file__).resolve().parents[1]


def static_inventory(repo: Path | None = None) -> Inventory:
    """Return the exact classified inventory, rejecting denominator drift."""
    root = _repo_root(repo)
    script_dir = root / "scripts" / "r_goldens"
    fixture_dir = root / "tests" / "r_goldens"
    scripts = tuple(sorted(script_dir.glob("*.R"), key=lambda path: path.name.encode()))
    fixtures = tuple(
        sorted(
            (path for path in fixture_dir.rglob("*") if path.is_file()),
            key=lambda path: path.relative_to(fixture_dir).as_posix().encode(),
        )
    )

    if len(scripts) != EXPECTED_SCRIPT_COUNT:
        raise ProvenanceError(
            f"R-generator denominator drift: expected {EXPECTED_SCRIPT_COUNT}, "
            f"found {len(scripts)}"
        )
    if len(fixtures) != EXPECTED_FIXTURE_COUNT:
        raise ProvenanceError(
            f"R-golden denominator drift: expected {EXPECTED_FIXTURE_COUNT}, "
            f"found {len(fixtures)}"
        )

    by_name = {path.name: path for path in scripts}
    special_names = HELPER_SCRIPTS | DIAGNOSTIC_SCRIPTS
    missing_special = special_names - by_name.keys()
    if missing_special:
        raise ProvenanceError(
            "declared non-writer scripts are missing: " + ", ".join(sorted(missing_special))
        )

    helpers = tuple(by_name[name] for name in sorted(HELPER_SCRIPTS))
    diagnostics = tuple(by_name[name] for name in sorted(DIAGNOSTIC_SCRIPTS))
    writers = tuple(path for path in scripts if path.name not in special_names)
    return Inventory(
        scripts=scripts,
        fixtures=fixtures,
        writers=writers,
        helpers=helpers,
        diagnostics=diagnostics,
    )


def _r_string_literals(text: str) -> set[str]:
    """Every string literal in R source, ignoring ``#`` comments."""
    literals: set[str] = set()
    index = 0
    while index < len(text):
        char = text[index]
        if char == "#":
            newline = text.find("\n", index)
            index = len(text) if newline < 0 else newline
            continue
        if char in "\"'":
            end = index + 1
            while end < len(text) and text[end] != char:
                end += 2 if text[end] == "\\" else 1
            literals.add(text[index + 1 : end])
            index = end + 1
            continue
        index += 1
    return literals


#: A literal that spells a fixture's file name in full (a suffix such as
#: "_matrices.txt" starts with "_" and is not one).
_FULL_FIXTURE_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.+-]*\.(?:csv|txt|tsv|dot)")
#: Full names writers use for files that are not fixtures: disstree2dot's
#: temporary DOT file under tempdir().
_NON_FIXTURE_FILE_NAMES = frozenset({"tree.dot"})


def _name_prefix(item: str, stem: str) -> bool:
    """Whether literal ``item`` can begin a built fixture name ``stem``.

    A built name's literal part contains ``_`` or ``-`` and the name goes on at a
    boundary: the literal ends in ``_``/``-``, or the next character of the stem
    is ``_``, ``-``, ``.`` or a digit (``paste0("md_input_d", k)``). An ordinary
    word such as ``"diss"`` from ``what = "diss"`` does not qualify
    (TX-AUD-20261010-056).
    """
    if len(item) < 3 or len(stem) <= len(item) or not stem.startswith(item):
        return False
    if "_" not in item and "-" not in item:
        return False
    return item[-1] in "_-" or stem[len(item)] in "_-.0123456789"


def static_attribution(repo: Path, inventory: Inventory) -> dict[str, set[str]]:
    """Name, for every fixture, the writer scripts whose source names it.

    A writer names a fixture when one of its string literals (or the last path
    component of one) is the fixture's relative path, file name or stem, or --
    for names built with ``paste0`` -- a name prefix of its stem (see
    :func:`_name_prefix`). This proves attribution by name only; that a writer
    actually writes a fixture is proved by ``--execute``. A fixture no writer
    names fails here, so a fixture written by anything else cannot pass on a
    bumped count alone (2026-09-22 audit S10-06); and a fixture name a writer
    spells in full must exist, so replacing such a fixture by an orphan of the
    same count fails too (TX-AUD-20261010-056). A built name is named by its
    prefix only: an orphan with that prefix (``md_input_d9.csv`` for
    ``md_input_d1.csv``) passes here, and ``--execute`` rejects it when the
    writer writes the real fixture again (batch-5 review M-2).
    """
    root = _repo_root(repo)
    fixture_dir = root / "tests" / "r_goldens"
    literals: dict[str, set[str]] = {}
    for path in inventory.writers:
        found = _r_string_literals(path.read_text(encoding="utf-8", errors="replace"))
        literals[path.name] = found | {item.rsplit("/", 1)[-1] for item in found}
    relatives = {fixture.relative_to(fixture_dir).as_posix() for fixture in inventory.fixtures}
    names = {fixture.name for fixture in inventory.fixtures}
    missing = sorted(
        f"{writer} names {item}"
        for writer, found in literals.items()
        for item in found
        if _FULL_FIXTURE_NAME.fullmatch(item.rsplit("/", 1)[-1])
        and item.rsplit("/", 1)[-1] not in _NON_FIXTURE_FILE_NAMES
        and item.removeprefix("tests/r_goldens/") not in relatives
        and item.rsplit("/", 1)[-1] not in names
        and "/" not in item.removeprefix("tests/r_goldens/")
    )
    attribution: dict[str, set[str]] = {}
    for fixture in inventory.fixtures:
        relative = fixture.relative_to(fixture_dir).as_posix()
        exact = {relative, fixture.name, fixture.stem, f"tests/r_goldens/{relative}"}
        owners = {name for name, found in literals.items() if found & exact}
        if not owners:
            owners = {
                name
                for name, found in literals.items()
                if any(_name_prefix(item, fixture.stem) for item in found)
            }
        attribution[relative] = owners
    unattributed = sorted(name for name, owners in attribution.items() if not owners)
    problems = []
    if missing:
        problems.append(
            "fixtures that a writer names in full but that do not exist: " + ", ".join(missing)
        )
    if unattributed:
        problems.append(
            "fixtures that no R writer names (by path, file name, stem or name prefix): "
            + ", ".join(unattributed)
        )
    if problems:
        raise ProvenanceError("; ".join(problems))
    return attribution


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _fixture_snapshot(repo: Path, inventory: Inventory) -> dict[str, str]:
    fixture_dir = repo / "tests" / "r_goldens"
    return {
        path.relative_to(fixture_dir).as_posix(): _sha256(path.read_bytes())
        for path in inventory.fixtures
    }


def _tree_hash(snapshot: Mapping[str, str]) -> str:
    digest = hashlib.sha256()
    for relative_path, file_hash in sorted(snapshot.items()):
        digest.update(relative_path.encode("utf-8"))
        digest.update(b"\0")
        digest.update(file_hash.encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def _script_role(path: Path) -> str:
    if path.name in HELPER_SCRIPTS:
        return "helper"
    if path.name in DIAGNOSTIC_SCRIPTS:
        return "diagnostic"
    return "writer"


def _decode_strace_path(value: str) -> str:
    return bytes(value, "utf-8").decode("unicode_escape")


def _direct_fixture_writes(trace: str, repo: Path, fixture_names: set[str]) -> set[str]:
    """Extract write-opened fixture paths from Linux strace output."""
    fixture_dir = repo / "tests" / "r_goldens"
    written: set[str] = set()
    for line in trace.splitlines():
        if not any(marker in line for marker in _WRITE_SYSCALL_MARKERS):
            continue
        for encoded_path in _QUOTED_PATH_RE.findall(line):
            raw_path = Path(_decode_strace_path(encoded_path))
            candidate = raw_path if raw_path.is_absolute() else repo / raw_path
            try:
                relative = candidate.resolve(strict=False).relative_to(fixture_dir)
            except ValueError:
                continue
            fixture_name = relative.as_posix()
            if fixture_name in fixture_names:
                written.add(fixture_name)
    return written


def validate_ownership(
    ownership: Mapping[str, set[str] | frozenset[str]],
    inventory: Inventory,
) -> None:
    """Require complete, singular ownership except for declared shared inputs."""
    fixture_dir = inventory.fixtures[0].parent
    expected_fixtures = {path.relative_to(fixture_dir).as_posix() for path in inventory.fixtures}
    actual_fixtures = set(ownership)
    missing = expected_fixtures - actual_fixtures
    unexpected = actual_fixtures - expected_fixtures
    if missing or unexpected:
        raise ProvenanceError(
            f"fixture ownership mismatch: missing={sorted(missing)!r}, "
            f"unexpected={sorted(unexpected)!r}"
        )

    expected_writer_names = {path.name for path in inventory.writers}
    writers_with_outputs = {owner for owners in ownership.values() for owner in owners}
    missing_writers = expected_writer_names - writers_with_outputs
    unexpected_writers = writers_with_outputs - expected_writer_names
    if missing_writers or unexpected_writers:
        raise ProvenanceError(
            f"canonical writer ownership mismatch: missing={sorted(missing_writers)!r}, "
            f"unexpected={sorted(unexpected_writers)!r}"
        )

    actual_shared = {
        fixture: frozenset(owners) for fixture, owners in ownership.items() if len(owners) > 1
    }
    declared_shared = {
        fixture: owners
        for fixture, owners in EXPECTED_SHARED_OWNERS.items()
        if fixture in expected_fixtures
    }
    if actual_shared != declared_shared:
        raise ProvenanceError(
            "shared-writer contract mismatch: "
            f"expected={declared_shared!r}, actual={actual_shared!r}"
        )


#: Third-party R packages that a generator depends on but that the retained CRAN
#: archives cannot pin. ``survival`` is a *recommended* package supplied by the R
#: installation, so it resolves from the system library even under ``R_LIBS``
#: pointing at the isolated oracle library. Its version is recorded here as
#: provenance -- informative when it changes, but never a fixture byte-diff, which
#: is what made the documented two-pass proof unrunnable (cycle-4b F-002).
UNPINNED_DEPENDENCY_FLOORS = {"survival": "3.8.9"}

#: R's base packages (``installed.packages(priority = "base")`` under R 4.6.2), which
#: ship with R itself and are recorded through R's version.
R_BASE_PACKAGES = frozenset(
    {
        "base",
        "compiler",
        "datasets",
        "graphics",
        "grDevices",
        "grid",
        "methods",
        "parallel",
        "splines",
        "stats",
        "stats4",
        "tcltk",
        "tools",
        "utils",
    }
)
_R_PACKAGE_LOAD = re.compile(
    r"\b(?:library|require|requireNamespace|loadNamespace)\(\s*[\"']?([A-Za-z][A-Za-z0-9.]*)"
)
_R_QUALIFIED_CALL = re.compile(r"(?<![A-Za-z0-9._])([A-Za-z][A-Za-z0-9.]*):::?(?=[A-Za-z.`])")


def generator_r_packages(scripts: Iterable[Path]) -> set[str]:
    """The non-base R packages the scripts load or call with ``::`` (comments ignored)."""
    packages: set[str] = set()
    for script in scripts:
        for line in script.read_text(encoding="utf-8").splitlines():
            code = line.split("#", 1)[0]
            for pattern in (_R_PACKAGE_LOAD, _R_QUALIFIED_CALL):
                packages.update(match.group(1) for match in pattern.finditer(code))
    return packages - R_BASE_PACKAGES


#: The R the oracle runs on (``scripts/r_goldens/README.md``: "Ubuntu 24.04 with R 4.6.2").
#: Asserted with the package versions, so a regeneration under another R stops
#: (TX-AUD-20261010-090 (d)).
ORACLE_R_VERSION = "4.6.2"


def _oracle_versions(
    rscript: Path, env: Mapping[str, str], repo: Path, packages: Iterable[str] = ()
) -> dict[str, Any]:
    """Assert the oracle's R and package versions, and record the R package closure.

    ``dependencies`` maps every package of the recursive Depends/Imports closure of
    TraMineR, TraMineRextras, the unpinned dependencies and ``packages`` (the ones the
    generators load), computed by R at run time, to the version R loads; a package of
    that closure that R cannot find stops the proof (TX-AUD-20261011-063).
    """
    unpinned = sorted(UNPINNED_DEPENDENCY_FLOORS)
    direct = sorted({"TraMineR", "TraMineRextras", *unpinned, *packages} - R_BASE_PACKAGES)
    reported = ["R.version.string", 'as.character(packageVersion("TraMineR"))']
    reported.append('as.character(packageVersion("TraMineRextras"))')
    reported.extend(f'as.character(packageVersion("{name}"))' for name in unpinned)
    floors = ", ".join(
        f'packageVersion("{name}") >= "{floor}"'
        for name, floor in sorted(UNPINNED_DEPENDENCY_FLOORS.items())
    )
    closure = (
        f'direct <- c({", ".join(json.dumps(name) for name in direct)}); '
        "ip <- installed.packages(); "
        'base <- rownames(ip)[ip[, "Priority"] %in% "base"]; '
        "reached <- unique(c(direct, unlist(tools::package_dependencies("
        'direct, db = ip, which = c("Depends", "Imports"), recursive = TRUE)))); '
        'reached <- setdiff(reached, c(base, "R")); '
        "found <- vapply(reached, function(p) tryCatch("
        'as.character(packageVersion(p)), error = function(e) "MISSING"), ""); '
        'cat(sprintf("dependency\\t%s\\t%s", reached, found), sep="\\n")'
    )
    expression = (
        f'cat({", ".join(reported)}, sep="\\n"); '
        f'stopifnot(getRversion() == "{ORACLE_R_VERSION}", '
        'packageVersion("TraMineR") == "2.2.14", '
        f'packageVersion("TraMineRextras") == "0.6.9", {floors}); '
        f"{closure}"
    )
    result = subprocess.run(
        [os.fspath(rscript), "-e", expression],
        cwd=repo,
        env=env,
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise ProvenanceError(f"exact R oracle validation failed:\n{result.stderr}")
    lines = result.stdout.splitlines()
    expected = 3 + len(unpinned)
    if len(lines) < expected:
        raise ProvenanceError(f"unexpected R oracle version output: {result.stdout!r}")
    versions: dict[str, Any] = {
        "r": lines[0],
        "TraMineR": lines[1],
        "TraMineRextras": lines[2],
    }
    for offset, name in enumerate(unpinned):
        versions[name] = lines[3 + offset]
    dependencies = {
        name: version
        for _, name, version in (
            line.split("\t") for line in lines[expected:] if line.startswith("dependency\t")
        )
    }
    missing = sorted(name for name, version in dependencies.items() if version == "MISSING")
    if missing:
        raise ProvenanceError(f"R packages of the oracle closure are not installed: {missing!r}")
    versions["dependencies"] = dict(sorted(dependencies.items(), key=lambda item: item[0].lower()))
    return versions


def _run_script(
    *,
    repo: Path,
    script: Path,
    rscript: Path,
    env: Mapping[str, str],
    fixture_names: set[str],
    trace_dir: Path,
) -> tuple[set[str], dict[str, Any]]:
    trace_path = trace_dir / f"{script.stem}.strace"
    command = [
        "strace",
        "-f",
        "-qq",
        "-e",
        "trace=%file",
        "-o",
        os.fspath(trace_path),
        os.fspath(rscript),
        os.fspath(script.relative_to(repo)),
    ]
    result = subprocess.run(
        command,
        cwd=repo,
        env=env,
        check=False,
        capture_output=True,
    )
    trace = trace_path.read_text(encoding="utf-8", errors="replace")
    direct_writes = _direct_fixture_writes(trace, repo, fixture_names)
    record = {
        "script": script.relative_to(repo).as_posix(),
        "role": _script_role(script),
        "returncode": result.returncode,
        "direct_writes": sorted(direct_writes),
        "stdout_sha256": _sha256(result.stdout),
        "stderr_sha256": _sha256(result.stderr),
    }
    if result.returncode != 0:
        stderr = result.stderr.decode("utf-8", errors="replace")
        raise ProvenanceError(f"{script.name} exited {result.returncode}:\n{stderr}")
    if record["role"] == "writer" and not direct_writes:
        raise ProvenanceError(f"canonical writer {script.name} wrote no fixture")
    if record["role"] != "writer" and direct_writes:
        raise ProvenanceError(
            f"{record['role']} {script.name} unexpectedly wrote {sorted(direct_writes)!r}"
        )
    if record["role"] == "diagnostic" and not result.stdout:
        raise ProvenanceError(f"stdout-only diagnostic {script.name} produced no stdout")
    if script.name in HELPER_SCRIPTS and b"CORESTATS_FIXTURES_OK" not in result.stdout:
        raise ProvenanceError(f"sourced helper health check failed for {script.name}")
    return direct_writes, record


def _refuse_in_place(repo: Path) -> None:
    """Refuse an ``--execute`` run that could overwrite fixtures no copy protects.

    ``--execute`` runs every generator with ``repo`` as its working directory,
    so a byte change is written into ``repo``'s fixtures before the check
    fails. ``repo`` must therefore not be the checkout holding this script
    (also not through a symlink), nor a git work tree with uncommitted changes
    under tests/r_goldens (TX-AUD-20261010-020). A plain copy, or a clean
    clone whose fixtures git can restore, is accepted.
    """
    resolved = repo.resolve()
    checkout = Path(__file__).resolve().parents[1]
    if resolved == checkout:
        raise ProvenanceError(
            f"--execute would rewrite the fixtures of the checkout {checkout}; run it on a "
            "disposable copy (--repo <copy>, scripts/r_goldens/README.md) or pass "
            "--allow-in-place"
        )
    inside = subprocess.run(
        ["git", "-C", os.fspath(resolved), "rev-parse", "--is-inside-work-tree"],
        capture_output=True,
        text=True,
        check=False,
    )
    if inside.returncode != 0 or inside.stdout.strip() != "true":
        return
    status = subprocess.run(
        ["git", "-C", os.fspath(resolved), "status", "--porcelain", "--", "tests/r_goldens"],
        capture_output=True,
        text=True,
        check=False,
    )
    if status.stdout.strip():
        raise ProvenanceError(
            f"{resolved} has uncommitted changes under tests/r_goldens that --execute would "
            "overwrite; commit or copy them first, or pass --allow-in-place"
        )


#: The locale and message language every script runs under. The scripts also
#: set both themselves (their first code line), so a script run by hand from
#: another shell writes the same bytes (TX-AUD-20261010-019).
_PINNED_LOCALE = {"LC_ALL": "C.UTF-8", "LANGUAGE": "en"}


def _child_environment(r_lib: Path | None) -> dict[str, str]:
    """The environment of the R child processes: the caller's, with the locale pinned."""
    env = os.environ.copy()
    env.update(_PINNED_LOCALE)
    if r_lib is not None:
        env["R_LIBS"] = os.fspath(r_lib.resolve())
    return env


def execute_provenance_check(
    repo: Path,
    *,
    rscript: Path,
    r_lib: Path | None,
    passes: int,
) -> dict[str, Any]:
    """Execute every declared role and prove deterministic writer coverage."""
    if passes < 2:
        raise ProvenanceError("at least two regeneration passes are required")
    if shutil.which("strace") is None:
        raise ProvenanceError("strace is required to prove direct fixture writes")

    inventory = static_inventory(repo)
    fixture_dir = repo / "tests" / "r_goldens"
    fixture_names = {path.relative_to(fixture_dir).as_posix() for path in inventory.fixtures}
    initial_snapshot = _fixture_snapshot(repo, inventory)
    initial_tree_hash = _tree_hash(initial_snapshot)
    env = _child_environment(r_lib)
    oracle = _oracle_versions(rscript, env, repo, generator_r_packages(inventory.scripts))

    pass_records: list[dict[str, Any]] = []
    expected_ownership: dict[str, set[str]] | None = None
    with tempfile.TemporaryDirectory(prefix="traminex-r-golden-traces-") as temp_dir:
        trace_root = Path(temp_dir)
        for pass_number in range(1, passes + 1):
            ownership: dict[str, set[str]] = {name: set() for name in fixture_names}
            script_records: list[dict[str, Any]] = []
            for script in inventory.scripts:
                before = _fixture_snapshot(repo, inventory)
                direct_writes, record = _run_script(
                    repo=repo,
                    script=script,
                    rscript=rscript,
                    env=env,
                    fixture_names=fixture_names,
                    trace_dir=trace_root,
                )
                current_inventory = static_inventory(repo)
                after = _fixture_snapshot(repo, current_inventory)
                changed = {name for name in fixture_names if before[name] != after[name]}
                untraced_changes = changed - direct_writes
                if untraced_changes:
                    raise ProvenanceError(
                        f"{script.name} changed fixtures without traced direct writes: "
                        f"{sorted(untraced_changes)!r}"
                    )
                for fixture in direct_writes:
                    ownership[fixture].add(script.name)
                script_records.append(record)

            validate_ownership(ownership, inventory)
            final_snapshot = _fixture_snapshot(repo, inventory)
            final_tree_hash = _tree_hash(final_snapshot)
            if final_snapshot != initial_snapshot:
                changed = sorted(
                    name
                    for name in fixture_names
                    if final_snapshot.get(name) != initial_snapshot.get(name)
                )
                raise ProvenanceError(
                    f"pass {pass_number} differs from reviewed fixture tree: {changed!r}"
                )
            if expected_ownership is not None and ownership != expected_ownership:
                raise ProvenanceError(
                    f"pass {pass_number} direct ownership differs from the first pass"
                )
            expected_ownership = ownership
            pass_records.append(
                {
                    "pass": pass_number,
                    "fixture_tree_sha256": final_tree_hash,
                    "scripts": script_records,
                }
            )

    assert expected_ownership is not None
    return {
        "schema_version": 1,
        "repo": os.fspath(repo),
        "oracle": oracle,
        "locale": {name: env[name] for name in _PINNED_LOCALE},
        "script_count": len(inventory.scripts),
        "writer_count": len(inventory.writers),
        "helper_count": len(inventory.helpers),
        "diagnostic_count": len(inventory.diagnostics),
        "fixture_count": len(inventory.fixtures),
        "reviewed_fixture_tree_sha256": initial_tree_hash,
        "ownership": {
            fixture: sorted(owners) for fixture, owners in sorted(expected_ownership.items())
        },
        "passes": pass_records,
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=None, help="repository root")
    parser.add_argument(
        "--execute",
        action="store_true",
        help="execute all roles and prove direct ownership and byte stability",
    )
    parser.add_argument(
        "--allow-in-place",
        action="store_true",
        help="allow --execute to rewrite the fixtures of the checkout or of a modified work tree",
    )
    parser.add_argument("--passes", type=int, default=2, help="regeneration passes")
    parser.add_argument("--rscript", type=Path, default=Path("/usr/bin/Rscript"))
    parser.add_argument("--r-lib", type=Path, default=None, help="isolated R library")
    parser.add_argument("--evidence", type=Path, default=None, help="write JSON evidence")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    repo = _repo_root(args.repo)
    try:
        inventory = static_inventory(repo)
        attribution = static_attribution(repo, inventory)
        if args.execute:
            if not args.allow_in_place:
                _refuse_in_place(repo)
            evidence = execute_provenance_check(
                repo,
                rscript=args.rscript,
                r_lib=args.r_lib,
                passes=args.passes,
            )
        else:
            evidence = {
                "schema_version": 1,
                "repo": os.fspath(repo),
                "script_count": len(inventory.scripts),
                "writer_count": len(inventory.writers),
                "helper_count": len(inventory.helpers),
                "diagnostic_count": len(inventory.diagnostics),
                "fixture_count": len(inventory.fixtures),
                "statically_attributed_fixture_count": len(attribution),
            }
        rendered = json.dumps(evidence, indent=2, sort_keys=True) + "\n"
        if args.evidence is not None:
            args.evidence.parent.mkdir(parents=True, exist_ok=True)
            args.evidence.write_text(rendered, encoding="utf-8")
        print(rendered, end="")
    except (OSError, ProvenanceError, subprocess.SubprocessError) as exc:
        print(f"R-golden provenance check failed: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
