"""Regressions for the 2026-09-02 audit's distance-family findings.

R anchors come from live TraMineR 2.2-14, recorded by
``scripts/r_goldens/make_goldens_audit20260922_s10.R`` in
``tests/r_goldens/inline_r_anchors.txt``, and from the vendored ``seqdist.R``.
"""

from __future__ import annotations

import subprocess
import sys

import numpy as np
import pytest
from r_goldens_utils import r_anchor, r_anchor_values

import traminex as tx
from traminex.distances import pairwise
from traminex.sequences import STSList

MAT = np.array([[0, 0, 1, 1, 2], [0, 1, 1, 2, 2], [2, 2, 1, 0, 0]])
X3 = STSList(MAT.copy(), alphabet=["A", "B", "C"])
MAT2 = np.array([[0, 0, 1, 1, 0], [0, 1, 1, 0, 0], [1, 1, 1, 0, 0]])
X2 = STSList(MAT2.copy(), alphabet=["A", "B"])


def test_omspell_rejects_negative_expcost():
    # R (seqdist.R:212-213): "'expcost' must be positive" for OMloc and OMspell.
    with pytest.raises(ValueError, match="expcost"):
        tx.seqdist(X3, method="OMspell", subm=2.0, indel=1.0, expcost=-5.0)
    with pytest.raises(ValueError, match="expcost"):
        tx.seqdist(X3, method="OMloc", subm=2.0, expcost=-5.0)


@pytest.mark.parametrize("method", ["OM", "OMloc", "OMslen", "OMspell", "OMstran", "TWED"])
def test_seqdist_requires_explicit_subm_like_r(method):
    # R refuses all six calls: "'sm' must be specified", except that OMstran
    # first asks for 'otto' and TWED for 'nu'.
    assert r_anchor(f"no_sm_{method}", "error")[0].endswith("must be specified")
    with pytest.raises(ValueError, match="subm"):
        tx.seqdist(X3, method=method)


def test_seqdist_autogenerates_subm_for_ham_dhd_lcs_like_r():
    # R auto-creates sm for HAM (CONSTANT 1), DHD (TRATE cval=4), LCS (CONSTANT 2).
    for method in ("HAM", "DHD", "LCS"):
        d = tx.seqdist(X3, method=method)
        assert d.shape == (3, 3)


def test_norm_length_is_a_true_maxlength_alias_for_the_classic_family():
    # R rejects the 'length' spelling, which traminex documents as an alias of
    # maxlength; the value is R's maxlength one.
    assert r_anchor("lcs_norm_length", "error")[0].startswith("'norm' must be one of")
    for method, kwargs in (
        ("LCS", {}),
        ("LCP", {}),
        ("RLCP", {}),
        ("OM", {"subm": 2.0, "indel": 1.0}),
    ):
        np.testing.assert_allclose(
            tx.seqdist(X3, method=method, norm="length", **kwargs),
            tx.seqdist(X3, method=method, norm="maxlength", **kwargs),
            atol=1e-15,
        )
    d = tx.seqdist(X3, method="LCS", norm="length")
    r_maxlength = r_anchor_values("lcs_maxlength").reshape(3, 3, order="F")
    assert d[0, 1] == pytest.approx(r_maxlength[0, 1], abs=1e-12)


def test_norm_max_is_a_true_maxdist_alias_even_for_asymmetric_subm():
    subm = np.array([[0.0, 1.5], [2.0, 0.0]])
    with pytest.warns(UserWarning):
        d_max = tx.seqdist(X2, method="OM", subm=subm, indel=2.0, norm="max")
    with pytest.warns(UserWarning):
        d_maxdist = tx.seqdist(X2, method="OM", subm=subm, indel=2.0, norm="maxdist")
    np.testing.assert_allclose(d_max, d_maxdist, atol=1e-15)
    # The bound takes the whole off-diagonal maximum, min(2.0, 2*indel) -> 2.0, so the value
    # is 3.5 / (5 * 2.0); R scans the strict upper triangle (1.5) and gives 3.5 / (5 * 1.5)
    # (TX-AUD-20261019-008, an R defect corrected).
    r_maxdist = r_anchor_values("om_asymmetric_maxdist").reshape(3, 3, order="F")
    assert d_maxdist[0, 1] == pytest.approx(3.5 / (5 * 2.0), abs=1e-12)
    assert r_maxdist[0, 1] == pytest.approx(3.5 / (5 * 1.5), abs=1e-12)


def test_norm_max_is_accepted_for_lcs_as_maxdist():
    np.testing.assert_allclose(
        tx.seqdist(X3, method="LCS", norm="max"),
        tx.seqdist(X3, method="LCS", norm="maxdist"),
        atol=1e-15,
    )


def test_pairwise_norm_max_uses_the_whole_off_diagonal_maxscost():
    # R's upper-triangle maxscost (1.5) would give 3.5 / (5 * 1.5) (TX-AUD-20261019-008).
    subm = np.array([[0.0, 1.5], [2.0, 0.0]])
    with pytest.warns(UserWarning):
        d = pairwise(MAT2, method="OM", subm=subm, indel=2.0, norm="max")
    assert d[0, 1] == pytest.approx(3.5 / (5 * 2.0), abs=1e-12)


def test_trate_smoothing_defaults_are_r_faithful_zero_everywhere():
    from traminex.distances import om_trate_pairwise
    from traminex.distances.omloc import omloc

    d_seqdist = tx.seqdist(X2, method="OM", subm="TRATE", indel=10.0)
    d_pairwise = pairwise(MAT2, method="OM", subm="TRATE", indel=10.0)
    d_om_trate = om_trate_pairwise(MAT2, indel=10.0, norm=None)
    np.testing.assert_array_equal(d_seqdist, d_pairwise)
    np.testing.assert_array_equal(d_seqdist, d_om_trate)
    d_default = omloc(MAT2)
    d_zero = omloc(MAT2, trate_smoothing=0.0)
    np.testing.assert_array_equal(d_default, d_zero)


def test_chi2_euclid_reject_odd_step_with_overlap_like_r():
    # R (seqdist.R:259-260): "'step' must be even when 'overlap = TRUE'".
    for method in ("CHI2", "EUCLID"):
        with pytest.raises(ValueError, match="step"):
            tx.seqdist(X3, method=method, step=3, overlap=True)
    d = tx.seqdist(X3, method="CHI2", step=2, overlap=True)
    assert d.shape == (3, 3)


def test_chi2_euclid_norm_flag_matches_r_vocabulary():
    # R accepts only "auto"/"none", and maps a logical norm onto them with a
    # deprecation notice (seqdist.R:34-37; 2026-09-22 audit S1a-12).
    for norm in ("true", "yes", "false", "no"):
        with pytest.raises(ValueError, match="norm"):
            tx.seqdist(X3, method="EUCLID", norm=norm)
    for norm in ("auto", "none", None):
        d = tx.seqdist(X3, method="EUCLID", step=1, norm=norm)
        assert d.shape == (3, 3)
    for flag, spelling in ((True, "auto"), (False, "none")):
        with pytest.warns(DeprecationWarning, match="norm"):
            d = tx.seqdist(X3, method="EUCLID", step=1, norm=flag)
        np.testing.assert_array_equal(d, tx.seqdist(X3, method="EUCLID", step=1, norm=spelling))


def test_traminex_binding_is_validated_at_import_time():
    code = "import traminex.distances"
    result = subprocess.run(
        [sys.executable, "-c", code],
        env={"PATH": "/usr/bin:/bin", "TRAMINEX_BINDING": "bogus"},
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert "TRAMINEX_BINDING" in result.stderr


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ \ at\ least,\ one\ substitution\ cost\ doesn't\ respect\ the\ triangle\ inequality\.:traminex.distances.python_kernels.CostDiagnosticMessage"
)
def test_om_trate_pairwise_resolves_its_automatic_indel_as_r_does():
    """The package's two TRATE+OM entry points must agree on ``indel=None``.

    ``om_trate_pairwise`` derived its automatic indel as
    ``0.5 * mean(off-diagonal)`` while ``seqdist(subm="TRATE", indel="auto")``
    uses R's ``max(sm)/2`` (``seqdist.R:359-365``). On a three-sequence corpus
    that produced ``[0.0, 0.3125, 0.625]`` against ``[0.0, 1.25, 4.0]`` -- the
    same nominal method giving different numbers depending on which public
    entry point was used (cycle-10 N-17 / 2026-09-20 T-06). The helper has no R
    counterpart, so R's resolution is the one both now share; pass ``indel=``
    explicitly for any other rule.
    """
    import numpy as np

    import traminex as tx
    from traminex.distances import om_trate_pairwise

    sequences = tx.seqdef(
        np.array([["A", "A", "B", "C"], ["A", "B", "B", "C"], ["C", "B", "A", "A"]], dtype=object)
    )
    helper = np.asarray(om_trate_pairwise(sequences.data.astype(np.int32), norm=None), dtype=float)
    dispatcher = np.asarray(
        tx.seqdist(sequences, method="OM", subm="TRATE", indel="auto"), dtype=float
    )
    np.testing.assert_allclose(helper, dispatcher, rtol=0, atol=1e-12)
