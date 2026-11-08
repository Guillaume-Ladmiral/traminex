"""Generator-authored fixture records are labelled (TX-AUD-20261011-020).

scripts/r_goldens/README.md says fixtures hold R output, but some records are
written by the generator itself: R statements re-implemented because the R
function draws and returns nothing (plot.stslist sampling, seqedplot series,
seqsplot panels), hand-typed messages, traminex-side metadata, a per-group
computation that is not dissvar.grp's output, R errors mapped to traminex
codes, and the MDS axis oriented by traminex's rule (TX-REP-20261017-006).
scripts/r_goldens/AUTHORED_RECORDS.csv names each one; this test keeps
the register true to the generators and fixtures it names.
"""

from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REGISTER = ROOT / "scripts" / "r_goldens" / "AUTHORED_RECORDS.csv"
COLUMNS = ["generator", "construct", "fixture", "records", "authored", "reason"]
CITED = {
    "make_goldens_plot_stslist_sampling_r4.R",
    "make_goldens_plot_samm_r4.R",
    "make_goldens_seqformat_full.R",
    "make_goldens_precarity_order.R",
    "make_goldens_seqedplot.R",
    "make_goldens_seqsurv.R",
    "make_goldens_dissvar_grp.R",
    "make_goldens_distances_void.R",
    # The MDS axis oriented by traminex's rule (TX-REP-20261017-006; batch-5 review M-7).
    "make_goldens_mds_orientation.R",
    "make_goldens_domassoc_args.R",
    "make_goldens_domassoc_unweighted.R",
}


def _rows() -> list[dict[str, str]]:
    with REGISTER.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames == COLUMNS
        return list(reader)


def test_every_row_names_a_real_construct_and_fixture() -> None:
    for row in _rows():
        generator = ROOT / "scripts" / "r_goldens" / row["generator"]
        assert row["construct"] in generator.read_text(encoding="utf-8"), row
        assert (ROOT / "tests" / "r_goldens" / row["fixture"]).is_file(), row
        assert row["records"] and row["authored"] and row["reason"], row


def test_the_cited_generators_are_covered() -> None:
    assert {row["generator"] for row in _rows()} >= CITED


def test_the_readme_rule_names_the_register() -> None:
    readme = (ROOT / "scripts" / "r_goldens" / "README.md").read_text(encoding="utf-8")
    assert "AUTHORED_RECORDS.csv" in readme
