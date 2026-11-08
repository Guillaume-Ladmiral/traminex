"""Regression tests for the complete R-golden writer provenance contract."""

import ast
import re
import shutil
from pathlib import Path

import pytest

from scripts import check_r_golden_provenance as provenance

REPO = Path(__file__).resolve().parents[1]


def test_static_inventory_classifies_every_script_and_fixture():
    inventory = provenance.static_inventory(REPO)

    # 2026-10-12 repair batch 1: + seqdef_r_pairing, seqecreate_rendering,
    # seqdist_weighted_refseq, seqidegrad_logical_pow, disstree_assign_names,
    # seqplotmd_sort_order generators.
    # 2026-10-13 repair batch 2: + seqdist_twed_auto_indel generator.
    # + dhd_explicit_norms generator.
    # + seqdist_nan_parameters generator.
    # + chi2_euclid_weights generator.
    # + seqedist_empty generator.
    # + subseqelist_subset generator.
    # + factor_level_order generator.
    # + sortv_character_order generator.
    # + sort_marker_layouts generator.
    # + seqecreate_label_order generator.
    # + sts_spell_ids generator.
    # + sts_sps_dss_void generator.
    # + hspell_alphabet_order generator.
    # + seqdef_nr_void_neutral generator.
    # + group_label_rendering generator.
    # + seqstatl_r_numbers generator.
    # + seqindic_dyn_inputs generator.
    # + seqpolyads_negative_u generator.
    # + seqpm_posix generator.
    # + internal_void_indicators generator.
    # + seqrecode_palette generator.
    # + scripts/r_goldens/make_goldens_datetime_labels generator.
    # + scripts/r_goldens/make_goldens_disstree_node_labels generator.
    # + scripts/r_goldens/make_goldens_seqplot_group_titles generator.
    # + scripts/r_goldens/make_goldens_dissmfacw_rank generator.
    # + scripts/r_goldens/make_goldens_seqplot_rf_stat generator.
    # + scripts/r_goldens/make_goldens_survival_extent generator.
    # + scripts/r_goldens/make_goldens_seqmeant_rows generator.
    # + scripts/r_goldens/make_goldens_ecmpgroup_colours generator.
    # + scripts/r_goldens/make_goldens_seqpcplot_positions generator.
    # + scripts/r_goldens/make_goldens_index_orientation generator.
    # + scripts/r_goldens/make_goldens_domassoc_window generator.
    # + scripts/r_goldens/make_goldens_domassoc_unweighted generator.
    # + scripts/r_goldens/make_goldens_domassoc_args generator.
    # + scripts/r_goldens/make_goldens_domassoc_missing_level generator.
    # + scripts/r_goldens/make_goldens_erulesdisc_aliased generator.
    # + scripts/r_goldens/make_goldens_dissmergegroups_negative generator.
    # + scripts/r_goldens/make_goldens_disscompare_order generator.
    # + scripts/r_goldens/make_goldens_seqplot_group_rep generator.
    # + scripts/r_goldens/make_goldens_seqnum generator.
    # + scripts/r_goldens/make_goldens_seqdim generator.
    # + scripts/r_goldens/make_goldens_mds_orientation generator.
    # + seqdist_r_arguments generator.
    # + check_costs generator.
    # + generated_cost_messages generator.
    # + seqdist_messages generator.
    # + degenerate_distances generator.
    # + refseq_forms generator.
    # + refseq_asymmetric generator.
    # + chi2_windows generator.
    # + omslen_link generator.
    # + default_alphabet generator (owner decision D2).
    # + seqpolyads_roles generator (TX-AUD-20261019-004).
    # + refseq_most_frequent generator (TX-AUD-20261019-001).
    # + indels_unobserved generator (TX-AUD-20261019-002).
    # + asym_list_matrix generator (TX-AUD-20261019-003).
    # + asym_bound generator (TX-AUD-20261019-008).
    # + seqdist_check_order generator (TX-REP-20261019-010).
    # + seqdiff_messages generator (TX-AUD-20261019-059).
    # + seqmd_details generator (TX-AUD-20261010-074).
    # + degenerate_omstran generator (TX-AUD-20261019-080).
    # + negative_prox generator (TX-REP-20261019-001).
    # + spell_powers generator (TX-AUD-20261019-030).
    # + subseq_overflow generator (TX-AUD-20261019-029).
    # + norm_bounds generator (TX-REP-20261021-003 to -005).
    # + spelling_convention generator (TX-AUD-20261019-061).
    # + event_arguments generator (batch 11).
    # + rank_and_convergence generator (TX-REP-20261019-002).
    # + diss_trees_reps generator (batch 12).
    # + plot_items generator (batch 13).
    # + construction_items generator (batch 14).
    # + statistics_items generator (batch 15).
    # + ledger_items generator (batch 18).
    # + disstree_weighted_gate and consumed_items generators (2026-11-02).
    # + study_findings generator (2026-11-07).
    assert len(inventory.scripts) == 209
    # 612 + 5 cycle-6 F-001 seqmaintokens default-k fixtures
    # + 4 cycle-6 F-003 case-weight edge-case fixtures
    # + 5 cycle-6 section-9 seqalign DP-matrix fixtures
    # + 5 cycle-7 W-01 deep/degenerate event-engine fixtures
    # + 12 cycle-7 W-02 deterministic group-comparison fixtures
    # + 9 cycle-7 W-05 seqdist Rd-option fixtures
    # + 1 cycle-7 W-07 R-formals fixture
    # + 15 cycle-7 W-04 wide event and group-comparison fixtures
    # + 4 seqdist with.missing-reset fixtures
    # + 1 HSPELL whole-matrix fixture
    # + 2 seqdist state-space fixtures
    # + 2 seqdist cost-check fixtures
    # + 2 asymmetric-duplicate fixtures
    # + 2 kweights fixtures
    # + 1 numeric-edge fixture
    # + 1 argument-form fixture
    # + 1 2026-09-22 audit s3a statistics fixture.
    # + 1 audit 2026-09-22 s3b degenerate and boundary statistics records.
    # + 4 2026-09-22 audit S5 event and SAMM fixtures.
    # + 1 audit-20260922 S2 construction/format fixture.
    # + 1 checkcost diagnostic fixture
    # + 1 inline-R-anchor fixture.
    # + 1 audit 2026-09-22 s3c follow-up statistics fixture.
    # + 11 2026-09-22 audit s4 discrepancy-analysis fixtures.
    # + 1 audit-20260922 S2-19 seqMD void-component cost fixture.
    # + 1 2026-09-22 audit S5-16 chi-square arithmetic fixture.
    # + 4 2026-09-22 audit S6 plotting-semantics fixtures
    # + 2 2026-09-22 audit s4b disstree/dissrep follow-up fixtures.
    # + 2 2026-09-22 audit S10-11/S10-12 seqMD cost fixtures.
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
    #   alphabet fixture (make_goldens_TSE_to_STS.R); +1 writer (owner decision D2).
    # + 1 seqpolyads role fixture (make_goldens_seqpolyads_roles.R); +1 writer.
    # + 1 refseq=0 fixture (make_goldens_refseq_most_frequent.R); +1 writer.
    # + 1 INDELS unobserved-state fixture (make_goldens_indels_unobserved.R); +1 writer.
    # + 1 nested-list matrix fixture (make_goldens_asym_list_matrix.R); +1 writer.
    # + 1 seqindic row-name fixture (make_goldens_seqindic.R; TX-AUD-20261010-093).
    # + 1 asymmetric-bound fixture (make_goldens_asym_bound.R; TX-AUD-20261019-008); +1 writer.
    # + 1 seqentrans fixture (make_goldens_seqentrans.R; TX-AUD-20261019-016).
    # + 1 argument-check fixture (make_goldens_seqdist_check_order.R); +1 writer.
    # + 1 seqdiff-conditions fixture (make_goldens_seqdiff_messages.R); +1 writer.
    # + 1 seqMD-details fixture (make_goldens_seqmd_details.R); +1 writer.
    # + 1 degenerate-OMstran fixture (make_goldens_degenerate_omstran.R); +1 writer.
    # + 1 negative-proximity fixture (make_goldens_negative_prox.R); +1 writer.
    # + 1 spell-power fixture (make_goldens_spell_powers.R); +1 writer.
    # + 1 subsequence-overflow fixture (make_goldens_subseq_overflow.R); +1 writer.
    # + 1 normalization-bound fixture (make_goldens_norm_bounds.R); +1 writer.
    # + 1 spelling-convention fixture (make_goldens_spelling_convention.R); +1 writer.
    # + 1 event-arguments fixture (make_goldens_event_arguments.R); +1 writer.
    # + 1 rank-and-convergence fixture (make_goldens_rank_and_convergence.R); +1 writer.
    # + 1 diss-trees-reps fixture (make_goldens_diss_trees_reps.R); +1 writer.
    # + disstree_weighted_gate.txt and consumed_items.csv (2026-11-02).
    # + study_findings.csv (2026-11-07).
    assert len(inventory.fixtures) == 819
    assert len(inventory.writers) == 205  # + their two generators, + make_goldens_study_findings.R
    assert [path.name for path in inventory.helpers] == ["make_goldens_corestats.R"]
    assert [path.name for path in inventory.diagnostics] == [
        "make_goldens_distances.R",
        "make_goldens_norm.R",
        "make_goldens_norm_costs.R",
    ]
    assert {path.name for path in inventory.writers}.isdisjoint(
        provenance.HELPER_SCRIPTS | provenance.DIAGNOSTIC_SCRIPTS
    )


def _tiny_repo(tmp_path: Path, fixtures: dict[str, str], writer_text: str) -> Path:
    """A repository with the four declared non-writers and one writer."""
    scripts = tmp_path / "scripts" / "r_goldens"
    scripts.mkdir(parents=True)
    for name in provenance.HELPER_SCRIPTS | provenance.DIAGNOSTIC_SCRIPTS:
        (scripts / name).write_text('cat("no fixture")\n', encoding="utf-8")
    (scripts / "make_goldens_example.R").write_text(writer_text, encoding="utf-8")
    goldens = tmp_path / "tests" / "r_goldens"
    goldens.mkdir(parents=True)
    for name, text in fixtures.items():
        (goldens / name).write_text(text, encoding="utf-8")
    return tmp_path


def test_static_mode_rejects_a_fixture_no_writer_names(tmp_path, monkeypatch, capsys):
    """Matching counts must not let an unattributed fixture through.

    Static mode checked only the script and fixture counts, so a fixture
    written by anything else plus a bumped count passed every automated gate
    (2026-09-22 audit S10-06).
    """
    repo = _tiny_repo(
        tmp_path,
        {"example_values.csv": "1\n", "orphan.csv": "2\n"},
        'writeLines("1", file.path("tests", "r_goldens", "example_values.csv"))\n'
        "# orphan.csv is only named in a comment\n",
    )
    monkeypatch.setattr(provenance, "EXPECTED_SCRIPT_COUNT", 5)
    monkeypatch.setattr(provenance, "EXPECTED_FIXTURE_COUNT", 2)

    assert provenance.main(["--repo", str(repo)]) == 1
    assert "orphan.csv" in capsys.readouterr().err


def test_static_attribution_accepts_literal_and_prefix_names(tmp_path, monkeypatch):
    repo = _tiny_repo(
        tmp_path,
        {"example_values.csv": "1\n", "example_case_a.csv": "2\n"},
        'writeLines("1", file.path("tests", "r_goldens", "example_values.csv"))\n'
        'for (case in c("a")) write.csv(1, paste0("tests/r_goldens/example_case_", case, ".csv"))\n',
    )
    monkeypatch.setattr(provenance, "EXPECTED_SCRIPT_COUNT", 5)
    monkeypatch.setattr(provenance, "EXPECTED_FIXTURE_COUNT", 2)

    attribution = provenance.static_attribution(repo, provenance.static_inventory(repo))

    assert attribution == {
        "example_case_a.csv": {"make_goldens_example.R"},
        "example_values.csv": {"make_goldens_example.R"},
    }


def test_every_tracked_fixture_is_named_by_a_writer():
    inventory = provenance.static_inventory(REPO)
    attribution = provenance.static_attribution(REPO, inventory)

    assert len(attribution) == len(inventory.fixtures)
    writers = {path.name for path in inventory.writers}
    assert all(owners and owners <= writers for owners in attribution.values())


def test_direct_writer_trace_ignores_reads_and_detects_relative_and_absolute_writes(
    tmp_path,
):
    fixture_dir = tmp_path / "tests" / "r_goldens"
    fixture_dir.mkdir(parents=True)
    absolute_fixture = fixture_dir / "absolute.csv"
    trace = "\n".join(
        [
            'openat(AT_FDCWD, "tests/r_goldens/read-only.csv", O_RDONLY) = 3',
            'openat(AT_FDCWD, "tests/r_goldens/relative.csv", '
            "O_WRONLY|O_CREAT|O_TRUNC, 0666) = 4",
            f'openat(AT_FDCWD, "{absolute_fixture}", ' "O_WRONLY|O_CREAT|O_TRUNC, 0666) = 5",
        ]
    )

    actual = provenance._direct_fixture_writes(
        trace,
        tmp_path,
        {"read-only.csv", "relative.csv", "absolute.csv"},
    )

    assert actual == {"relative.csv", "absolute.csv"}


def _small_inventory(tmp_path: Path) -> provenance.Inventory:
    fixture_dir = tmp_path / "tests" / "r_goldens"
    fixture_dir.mkdir(parents=True)
    fixtures = (fixture_dir / "one.csv", fixture_dir / "two.csv")
    for fixture in fixtures:
        fixture.touch()
    scripts = (tmp_path / "writer_one.R", tmp_path / "writer_two.R")
    return provenance.Inventory(
        scripts=scripts,
        fixtures=fixtures,
        writers=scripts,
        helpers=(),
        diagnostics=(),
    )


def test_validate_ownership_requires_every_fixture_and_writer(tmp_path):
    inventory = _small_inventory(tmp_path)

    with pytest.raises(provenance.ProvenanceError, match="fixture ownership mismatch"):
        provenance.validate_ownership({"one.csv": {"writer_one.R"}}, inventory)


def test_validate_ownership_rejects_undeclared_shared_writers(tmp_path):
    inventory = _small_inventory(tmp_path)
    ownership = {
        "one.csv": {"writer_one.R", "writer_two.R"},
        "two.csv": {"writer_two.R"},
    }

    with pytest.raises(provenance.ProvenanceError, match="shared-writer contract mismatch"):
        provenance.validate_ownership(ownership, inventory)


def test_validate_ownership_accepts_complete_single_writer_contract(tmp_path):
    inventory = _small_inventory(tmp_path)
    ownership = {
        "one.csv": {"writer_one.R"},
        "two.csv": {"writer_two.R"},
    }

    provenance.validate_ownership(ownership, inventory)


#: The oracle every generator must be pinned to, mirroring ``SOURCES/README.md``
#: and the ``--execute`` guard in ``check_r_golden_provenance.py``.
ORACLE_VERSION = "2.2.14"


def _r_code(text: str) -> str:
    """R source with its comments removed: a ``#`` outside a string starts one."""
    kept = []
    for line in text.splitlines():
        quote = None
        escaped = False
        for index, char in enumerate(line):
            if escaped:
                escaped = False
            elif char == "\\" and quote:
                escaped = True
            elif quote:
                quote = None if char == quote else quote
            elif char in "\"'":
                quote = char
            elif char == "#":
                line = line[:index]
                break
        kept.append(line)
    return "\n".join(kept)


#: ``packageVersion("TraMineR")``, possibly inside ``as.character()``, compared with the
#: oracle version by ``==``, ``!=`` or ``identical(..., version)``.
_ORACLE_PIN = re.compile(
    r'packageVersion\("TraMineR"\)\)?\s*(?P<op>==|!=|,)\s*"' + re.escape("2.2.14") + '"'
)


def _enclosing_calls(code: str, position: int) -> list[tuple[str, int]]:
    """The calls (and ``if``) around ``position``, innermost first: name and the index of
    its ``(``; a leading ``!`` is kept with the name."""
    calls = []
    depth = 0
    for index in range(position - 1, -1, -1):
        if code[index] == ")":
            depth += 1
        elif code[index] == "(":
            if depth:
                depth -= 1
                continue
            match = re.search(r"(!?)\s*([A-Za-z_.][A-Za-z0-9_.]*)\s*$", code[:index])
            calls.append(((match.group(1) + match.group(2)) if match else "", index))
    return calls


def _stops_after(code: str, open_index: int) -> bool:
    """Whether the body of the ``if`` whose condition opens at ``open_index`` calls stop()."""
    depth = 0
    for index in range(open_index, len(code)):
        depth += {"(": 1, ")": -1}.get(code[index], 0)
        if depth == 0:
            break
    rest = code[index + 1 :].lstrip()
    if rest.startswith("{"):
        depth = 0
        for end, char in enumerate(rest):
            depth += {"{": 1, "}": -1}.get(char, 0)
            if depth == 0:
                break
        body = rest[: end + 1]
    else:
        body = rest.split("\n", 1)[0]
    return re.search(r"\bstop\s*\(", body) is not None


def _pins_the_oracle(text: str) -> bool:
    """Whether a generator's code asserts the oracle's TraMineR version.

    A pin is ``stopifnot()`` of an equality with "2.2.14" (``==``, or an ``identical()``
    that is not negated), or an ``if`` on an inequality (``!=`` or ``!identical()``) whose
    body calls ``stop()``. The text ``packageVersion("TraMineR")`` alone is not one: in a
    comment, printed, stored, or compared with another version, it lets a regeneration
    run against the wrong TraMineR (TX-AUD-20261010-090 (c); batch-16 review M-3).
    """
    code = _r_code(text)
    for match in _ORACLE_PIN.finditer(code):
        calls = [
            call for call in _enclosing_calls(code, match.start()) if call[0] != "as.character"
        ]
        operator = match.group("op")
        if operator == ",":
            if not calls or calls[0][0] not in ("identical", "!identical"):
                continue
            negated = calls[0][0].startswith("!")
            calls = calls[1:]
        else:
            negated = operator == "!="
        if not calls:
            continue
        name, open_index = calls[0]
        if name == "stopifnot" and not negated:
            return True
        if name == "if" and negated and _stops_after(code, open_index):
            return True
    return False


def test_every_generator_pins_the_oracle_version():
    """A generator that does not assert its oracle can silently use the wrong one.

    ``scripts/r_goldens/README.md`` requires that "Every new or changed
    generator must also assert the package version inside the script", but 73
    of the 108 generators never mentioned ``packageVersion`` at all. That is
    harmless only while one R library is installed; during an oracle
    substitution it lets a regeneration run against the wrong TraMineR and
    succeed, writing fixtures that claim a provenance they do not have
    (2026-09-20 G-01).
    """
    generators = sorted((REPO / "scripts" / "r_goldens").glob("*.R"))
    assert generators, "no golden generators found"

    unpinned = [
        path.name for path in generators if not _pins_the_oracle(path.read_text(encoding="utf-8"))
    ]
    assert not unpinned, (
        f"{len(unpinned)} of {len(generators)} generators do not assert the TraMineR "
        f"oracle version:\n  " + "\n  ".join(unpinned)
    )


def test_no_generator_pins_a_superseded_oracle_version():
    """Every asserted version must be the current oracle, not a stale one."""
    import re

    stale: list[str] = []
    for path in sorted((REPO / "scripts" / "r_goldens").glob("*.R")):
        text = path.read_text(encoding="utf-8")
        for found in re.findall(r'packageVersion\("TraMineR"\)\)?\s*[!=]=\s*"([0-9.]+)"', text):
            if found != ORACLE_VERSION:
                stale.append(f"{path.name}: {found}")
    assert not stale, "generators pinned to a superseded oracle:\n  " + "\n  ".join(stale)


_SOURCE_CALL = re.compile(r'source\(\s*(?:file\.path\(([^)]*)\)|"([^"]+)")')


def _sourced_generators(path: Path) -> set[str]:
    """The generator files a script ``source()``s, read from its uncommented lines."""
    names = set()
    code = "\n".join(
        line
        for line in path.read_text(encoding="utf-8").splitlines()
        if not line.lstrip().startswith("#")
    )
    for parts, literal in _SOURCE_CALL.findall(code):
        target = literal or "/".join(re.findall(r'"([^"]+)"', parts))
        names.add(target.rsplit("/", 1)[-1])
    return names


def test_shared_writer_sets_follow_source_calls():
    # A generator that source()s make_goldens_md_inputs.R or make_goldens_phase2_inputs.R
    # runs it, and so rewrites the shared input CSVs (TX-AUD-20261019-023: four batch-4
    # generators were missing from the sets, which only the two-pass proof showed).
    generators = sorted((REPO / "scripts" / "r_goldens").glob("make_goldens_*.R"))
    direct = {path.name: _sourced_generators(path) for path in generators}

    def reached(name: str) -> set[str]:
        seen: set[str] = set()
        stack = [name]
        while stack:
            for target in direct.get(stack.pop(), set()):
                if target not in seen:
                    seen.add(target)
                    stack.append(target)
        return seen

    def writers_of(producer: str) -> set[str]:
        return {producer} | {name for name in direct if producer in reached(name)}

    assert set(provenance._MD_SHARED_WRITERS) == writers_of("make_goldens_md_inputs.R")
    assert set(provenance._PHASE2_SHARED_WRITERS) == writers_of("make_goldens_phase2_inputs.R")


@pytest.mark.parametrize(
    ("text", "pinned"),
    [
        ('stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")\n', True),
        (
            'if (!identical(as.character(packageVersion("TraMineR")), "2.2.14")) {\n  stop()\n}\n',
            True,
        ),
        ('stopifnot(packageVersion("TraMineR") == "2.2.14",\n          TRUE)\n', True),
        ('if (packageVersion("TraMineR") != "2.2.14") stop("wrong oracle")\n', True),
        (
            'stopifnot(\n  as.character(packageVersion("TraMineR")) == "2.2.14",\n  TRUE\n)\n',
            True,
        ),
        (
            'if (getRversion() != "4.6.2" ||\n    packageVersion("TraMineR") != "2.2.14" ||\n'
            '    FALSE) {\n  stop("wrong oracle")\n}\n',
            True,
        ),
        ('# stopifnot(packageVersion("TraMineR") == "2.2.14")\nlibrary(TraMineR)\n', False),
        ('x <- 1  # packageVersion("TraMineR") == "2.2.14"\n', False),
        ('message(packageVersion("TraMineR"))\n', False),
        ('stopifnot(packageVersion("TraMineR") == "2.2.13")\n', False),
        # Code that compares or prints but asserts nothing (batch-16 review M-3).
        ('cat("# text", as.character(packageVersion("TraMineR")) == "2.2.14")\n', False),
        ('message(packageVersion("TraMineR"), "2.2.14")\n', False),
        ('ok <- packageVersion("TraMineR") == "2.2.14"\n', False),
        ('stopifnot(packageVersion("TraMineR") != "2.2.14")\n', False),
        ('if (packageVersion("TraMineR") != "2.2.14") message("old oracle")\n', False),
    ],
)
def test_the_pin_check_reads_an_uncommented_comparison_with_the_oracle(text, pinned):
    """A pin is an assertion: ``stopifnot`` of an equality, or ``if`` of an inequality that
    stops. A comment, a printed or stored comparison, or another version is not
    (TX-AUD-20261010-090 (c); batch-16 review M-3)."""
    assert _pins_the_oracle(text) is pinned


_RSCRIPT = shutil.which("Rscript")


@pytest.mark.skipif(_RSCRIPT is None, reason="Rscript is not installed")
def test_the_oracle_check_refuses_another_r_version(monkeypatch, tmp_path):
    """``--execute`` checks R itself, not only the packages (TX-AUD-20261010-090 (d)).

    scripts/r_goldens/README.md names R 4.6.2 as the oracle; the check printed
    R.version.string but asserted the package versions only.
    """
    monkeypatch.setattr(provenance, "ORACLE_R_VERSION", "9.9.9", raising=False)

    with pytest.raises(provenance.ProvenanceError, match=r"getRversion\(\) == \"9\.9\.9\""):
        provenance._oracle_versions(Path(_RSCRIPT), provenance._child_environment(None), tmp_path)


#: Test names that claim agreement with R (the A8 agent's scan of the 2026-10-19 audit), and
#: the names that show a module reads a recorded R value or runs R. The audit's scan also
#: accepted any ".csv" or ".txt" in the module, which an unrelated file name satisfies.
_R_CLAIM = re.compile(
    r"(_|^)(as_r|as_in_r|match(es)?_r|like_r|r_value|r_result|matches_r|agrees_with_r"
    r"|equals_r|r_parity|same_as_r|with_r)(_|$)"
)
_GOLDEN_SOURCE = re.compile(
    r"load_golden|r_anchor|R_GOLDENS_DIR|r_goldens|assert_matches_golden|Rscript"
)

#: Modules whose R-named tests rest on R's source or platform, not on an R output value.
_R_CLAIM_EXEMPT = {
    "test_r_sum_shared.py": "checks that the 80-bit fixture is skipped without R's 63-bit mantissa",
    "test_seqrfplot_box_geometry.py": "R's box geometry from its source formula (dissrf.R:249; "
    "plot.seqrf.R:79-90), a definition",
}


def _modules_claiming_r_without_a_golden(paths: list[Path]) -> list[str]:
    """Modules with a test named as agreeing with R that never read a golden or run R."""
    flagged = []
    for path in paths:
        source = path.read_text(encoding="utf-8")
        claims = [
            node.name
            for node in ast.walk(ast.parse(source))
            if isinstance(node, ast.FunctionDef)
            and node.name.startswith("test")
            and _R_CLAIM.search(node.name)
        ]
        if claims and not _GOLDEN_SOURCE.search(source):
            flagged.append(f"{path.name}: {', '.join(claims)}")
    return flagged


def test_modules_whose_tests_claim_r_read_a_golden():
    """A value a test calls R's comes from a generator's golden (TX-AUD-20261019-088).

    Values written into a test by hand sit outside the provenance gate: a regenerated
    golden does not reach them. The check is by module and by name, as the audit's scan
    was; it does not prove that every R claim is backed.
    """
    modules = sorted(Path(__file__).resolve().parent.glob("test_*.py"))
    flagged = {entry.split(":")[0] for entry in _modules_claiming_r_without_a_golden(modules)}

    assert flagged == set(_R_CLAIM_EXEMPT)


def test_the_r_claim_scan_flags_a_module_without_a_golden(tmp_path):
    module = tmp_path / "test_example.py"
    module.write_text("def test_distance_matches_r():\n    assert 1.5 == 1.5\n", encoding="utf-8")
    assert _modules_claiming_r_without_a_golden([module]) == [
        "test_example.py: test_distance_matches_r"
    ]

    module.write_text(
        "from r_goldens_utils import load_golden\n\n\ndef test_distance_matches_r():\n"
        "    assert load_golden('x')\n",
        encoding="utf-8",
    )
    assert _modules_claiming_r_without_a_golden([module]) == []


# TX-AUD-20261011-063: the --execute proof recorded R, TraMineR, TraMineRextras and
# survival only. The generators also load WeightedCluster, vegan, colorspace and
# RColorBrewer, which load further packages (cluster among them); the proof now records
# the version of every package in that Depends/Imports closure.


def test_generator_r_packages_reads_loads_and_qualified_calls(tmp_path):
    script = tmp_path / "make_goldens_example.R"
    script.write_text(
        'library(TraMineR)\nrequire("vegan")\nsuppressMessages(requireNamespace("colorspace"))\n'
        "fit <- WeightedCluster::wcKMedoids\nutils::head(fit)\n# library(NotLoaded)\n",
        encoding="utf-8",
    )

    assert provenance.generator_r_packages([script]) == {
        "TraMineR",
        "WeightedCluster",
        "colorspace",
        "vegan",
    }


def test_the_generators_load_these_r_packages_directly():
    scripts = provenance.static_inventory(REPO).scripts

    assert provenance.generator_r_packages(scripts) == {
        "RColorBrewer",
        "TraMineR",
        "TraMineRextras",
        "WeightedCluster",
        "colorspace",
        "survival",
        "vegan",
    }


@pytest.mark.skipif(_RSCRIPT is None, reason="Rscript is not installed")
def test_the_oracle_record_holds_the_r_package_closure(tmp_path):
    packages = provenance.generator_r_packages(provenance.static_inventory(REPO).scripts)

    versions = provenance._oracle_versions(
        Path(_RSCRIPT), provenance._child_environment(None), tmp_path, packages
    )

    dependencies = versions["dependencies"]
    reached = {"WeightedCluster", "vegan", "cluster", "colorspace", "RColorBrewer", "survival"}
    assert reached <= set(dependencies)
    assert dependencies["TraMineR"] == versions["TraMineR"]
    assert all(dependencies.values())


@pytest.mark.skipif(_RSCRIPT is None, reason="Rscript is not installed")
def test_the_oracle_record_names_a_package_r_cannot_find(tmp_path):
    with pytest.raises(provenance.ProvenanceError, match="NoSuchPackageForTraminex"):
        provenance._oracle_versions(
            Path(_RSCRIPT),
            provenance._child_environment(None),
            tmp_path,
            ["NoSuchPackageForTraminex"],
        )


def test_the_ci_workflow_stops_on_an_r_other_than_the_oracle():
    # TX-AUD-20261011-079, option (b): the apt install stays, and the oracle step fails when R's version
    # differs from the oracle's, so a run against another R cannot pass silently.
    workflow = (Path(__file__).resolve().parents[1] / ".github" / "workflows" / "ci.yml").read_text(
        encoding="utf-8"
    )
    assert f'getRversion() == "{provenance.ORACLE_R_VERSION}"' in workflow
