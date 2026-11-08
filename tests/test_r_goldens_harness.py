"""Phase 0.1 smoke tests: the shared R-golden harness works end to end.

The golden used here (``tests/r_goldens/OMspell_exp05.csv``) carries the same
TraMineR 2.2-14 values as the inline ``OMspell_exp05`` case in
``test_r_goldens_distances.py``; this module validates the file-based loading
path that the Phase 2 goldens (dissvar/dissassoc/disscenter/disstree) will use.
"""

from pathlib import Path

import numpy as np
import pytest
from r_goldens_utils import (
    R_GOLDENS_DIR,
    assert_matches_golden,
    load_golden,
    parse_r_matrix,
)

import traminex as tx

REPO = Path(__file__).resolve().parents[1]


def test_parse_r_matrix_skips_comments_and_blanks():
    M = parse_r_matrix("# generating R call\n\n0,1.5\n2.5,0\n")
    np.testing.assert_array_equal(M, np.array([[0.0, 1.5], [2.5, 0.0]]))


def test_parse_r_matrix_rejects_empty_input():
    with pytest.raises(ValueError):
        parse_r_matrix("# only a comment\n")


def test_load_golden_missing_name_is_actionable():
    with pytest.raises(FileNotFoundError, match="scripts/r_goldens"):
        load_golden("no_such_golden")


def test_goldens_dir_exists():
    assert R_GOLDENS_DIR.is_dir()


def test_load_golden_omspell_is_a_distance_matrix():
    G = load_golden("OMspell_exp05")
    assert G.shape == (6, 6)
    np.testing.assert_array_equal(np.diag(G), np.zeros(6))
    np.testing.assert_array_equal(G, G.T)


def test_traminex_omspell_matches_file_golden():
    seqs = np.array(
        [
            ["A", "A", "A", "B", "B", "C", "C", "C"],
            ["A", "B", "A", "B", "C", "B", "C", "A"],
            ["C", "C", "C", "C", "A", "A", "A", "A"],
            ["A", "A", "B", "B", "B", "B", "C", "C"],
            ["B", "B", "B", "B", "B", "B", "B", "B"],
            ["C", "A", "C", "A", "B", "C", "B", "A"],
        ],
        dtype=object,
    )
    sts = tx.seqdef(seqs, alphabet=["A", "B", "C"])
    D = tx.seqdist(sts, method="OMspell", indel=1.0, subm=2.0, expcost=0.5)
    assert_matches_golden(D, "OMspell_exp05")


def test_omspell_golden_has_an_authoritative_direct_r_writer():
    source = (REPO / "scripts" / "r_goldens" / "make_goldens_spell_methods.R").read_text(
        encoding="utf-8"
    )

    assert 'packageVersion("TraMineR") == "2.2.14"' in source
    assert 'file.path(out_dir, "OMspell_exp05.csv")' in source
    assert "formatC" in source


def test_r_golden_operations_define_canonical_script_roles_and_coverage_check():
    operations = (REPO / "scripts" / "r_goldens" / "README.md").read_text(encoding="utf-8")

    assert "canonical fixture writers" in operations
    assert "sourced helper" in operations
    assert "stdout-only diagnostics" in operations
    assert "check_r_golden_provenance.py" in operations
