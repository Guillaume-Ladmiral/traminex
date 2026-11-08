import warnings

import numpy as np
import pytest
from corestats_fixtures import fixture_b
from r_goldens_utils import r_anchor_values

from traminex import seqdef, seqdist, seqRLCP
from traminex.distances import pairwise


def _toy_sts():
    mat = np.array(
        [
            ["A", "B", "C", "C"],
            ["A", "C", "B", "C"],
            ["C", "C", "B", "A"],
        ],
        dtype=object,
    )
    return seqdef(mat)


def test_rlcp_wrapper_accepts_norm_length():
    x = _toy_sts()

    d_rlcp_raw = seqRLCP(x)
    d_rlcp_len = seqRLCP(x, norm="length")

    assert d_rlcp_raw.shape == d_rlcp_len.shape == (3, 3)
    assert np.mean(d_rlcp_len) <= np.mean(d_rlcp_raw)


def test_rlcp_wrapper_rejects_invalid_norm():
    x = _toy_sts()
    with pytest.raises(ValueError, match="norm"):
        seqRLCP(x, norm="max")


def test_rlcp_wrapper_rejects_legacy_normalize_kwarg():
    x = _toy_sts()
    with pytest.raises(TypeError):
        seqRLCP(x, normalize=True)


def test_lcs_lcp_distance_matrices_still_available_via_seqdist():
    # The LCS/LCP *distance matrix* moved off seqLLCS/seqLLCP (now scalar, R-faithful)
    # onto seqdist(method=...), which is the supported way to get the matrix form.
    x = _toy_sts()
    d_lcs = seqdist(x, method="LCS")
    d_lcp = seqdist(x, method="LCP")
    assert d_lcs.shape == d_lcp.shape == (3, 3)
    assert np.allclose(np.diag(d_lcs), 0.0)
    assert np.allclose(np.diag(d_lcp), 0.0)


#: TraMineR 2.2-14 on corestats fixture B, rows 1 and 2, recorded by
#: scripts/r_goldens/make_goldens_audit20260922_s10.R (tests/r_goldens/inline_r_anchors.txt).
_NORMS = ("auto", "maxlength", "gmean", "maxdist", "YujianBo")
_DENSE_CLASSIC_R_PAIR = {
    (method, norm): float(value)
    for method in ("OM", "LCS", "LCP", "RLCP", "HAM", "DHD")
    for norm, value in zip(_NORMS, r_anchor_values(f"fixture_b_d12_{method}"), strict=True)
}


@pytest.mark.parametrize(("method", "norm"), _DENSE_CLASSIC_R_PAIR)
def test_seqdist_dense_classic_norm_vocabulary_matches_r(method, norm):
    """Dense seqdist keeps the R vocabulary at the R-facing boundary.

    Values are TraMineR 2.2-14's for fixture B, rows 1 and 2 (recorded by
    ``scripts/r_goldens/make_goldens_audit20260922_s10.R``).
    Direct ``pairwise`` retains its older ``length``/``max`` vocabulary.
    """
    kwargs = {"subm": 2.0, "indel": 1.0} if method == "OM" else {}
    actual = seqdist(
        fixture_b(),
        method=method,
        with_missing=True,
        norm=norm,
        **kwargs,
    )
    assert actual[0, 1] == pytest.approx(_DENSE_CLASSIC_R_PAIR[(method, norm)], rel=0.0, abs=1e-10)


_R_NORM_METHODS = {
    "OM": {"subm": 2.0, "indel": 2.0},
    "OMLOC": {"subm": 2.0},
    "OMSLEN": {"subm": 2.0, "indel": 2.0},
    "OMSPELL": {"subm": 2.0, "indel": 2.0},
    "OMSTRAN": {"subm": np.full((3, 3), 2.0) - 2.0 * np.eye(3), "indel": 2.0, "otto": 0.5},
    "TWED": {"subm": 2.0, "nu": 0.5, "h": 0.5},
    "HAM": {},
    "DHD": {},
    "LCS": {},
    "LCP": {},
    "RLCP": {},
}


@pytest.mark.parametrize("method", sorted(_R_NORM_METHODS))
@pytest.mark.parametrize(("alias", "canonical"), [("length", "maxlength"), ("max", "maxdist")])
def test_seqdist_aliases_mean_the_traminer_name_for_every_normalizing_method(
    method, alias, canonical
):
    """``seqdist``'s ``'length'``/``'max'`` are R's ``maxlength``/``maxdist`` everywhere.

    The docs said ``'max'`` was OM-only; ``seqdist`` accepts it for every method
    that takes a TraMineR ``norm`` (audit S1a-13). The indel of 2 separates R's
    ``maxlength`` from the low-level ``'length'``, which omits the indel factor.
    """
    x = seqdef(
        np.array([list("AABBC"), list("ABBCC"), list("CCBAA")], dtype=object),
        alphabet=["A", "B", "C"],
    )
    options = _R_NORM_METHODS[method]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        aliased = seqdist(x, method=method, norm=alias, **options)
        named = seqdist(x, method=method, norm=canonical, **options)
    np.testing.assert_array_equal(aliased, named)


def test_pairwise_dhd_norm_is_applied():
    """The low-level default DHD normalises as seqdist does (TX-AUD-20261010-023).

    It used to ignore ``norm`` and return the raw sums whatever was asked. The
    sequences have no missing cells: seqdist would then apply R's fixed
    missing-state costs, which the code-only low level does not know.
    """
    seqs = seqdef(
        np.array([list("AABBC"), list("ABBCC"), list("CCBAA"), list("BBCCA")], dtype=object),
        alphabet=["A", "B", "C"],
    )
    x = np.asarray(seqs.data, dtype=np.int32)
    np.testing.assert_array_equal(
        pairwise(x, method="DHD", norm="auto"), seqdist(seqs, method="DHD", norm="auto")
    )
    assert not np.array_equal(pairwise(x, method="DHD", norm="auto"), pairwise(x, method="DHD"))


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(r"ignore:\s*Non\-symmetric\ substitution\ cost\ matrix!:UserWarning")
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ \ at\ least,\ one\ substitution\ cost\ doesn't\ respect\ the\ triangle\ inequality\.:traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.parametrize("norm, expected", [("maxdist", 0.1), ("gmean", 0.55)])
def test_seqdist_om_asymmetric_subm_uses_the_whole_off_diagonal_bound(norm, expected):
    """TraMineR scans only ``i < j`` for maxcost, so its bound is not the maximum distance.

    d(A, B) = 1; the bound is min(10, 2 * 10) * 1 = 10 (the whole off-diagonal maximum), so
    maxdist = 0.1 and gmean = 1 - (10 - 1) / (2 * sqrt(10 * 10)) = 0.55; R's upper-triangle
    bound 1 gives 1 and 1 (TX-AUD-20261019-008, an R defect corrected).
    """
    (r_value,) = r_anchor_values(f"om_asymmetric_upper_bound_{norm}")
    x = seqdef(np.array([["A"], ["B"]], dtype=object), alphabet=["A", "B"])
    subm = np.array([[0.0, 1.0], [10.0, 0.0]])

    actual = seqdist(x, method="OM", subm=subm, indel=10.0, norm=norm)

    assert actual[0, 1] == pytest.approx(expected, rel=0.0, abs=1e-12)
    assert r_value == pytest.approx(1.0, abs=1e-12)


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ Empty\ sequence\(s\)\ \ \d+\ \ may\ cause\ inconsistent\ results!:UserWarning"
)
def test_seqdist_twed_empty_border_uses_auto_indel_cost():
    """TWEDdistance.cpp charges ``2*max(sm) + nu + h`` for an empty border."""
    x = seqdef(
        np.array([["%"], ["A"]], dtype=object),
        alphabet=["A"],
        nr="*",
        void="%",
    )

    actual = seqdist(
        x,
        method="TWED",
        subm=2.0,
        nu=0.5,
        h=0.5,
        norm="none",
    )

    assert actual[0, 1] == pytest.approx(5.0, rel=0.0, abs=1e-10)


def test_seqdist_dhd_uses_attached_sequence_weights():
    actual = seqdist(fixture_b(weighted=True), method="DHD", with_missing=True, norm="none")

    np.testing.assert_allclose(
        actual[np.triu_indices(4, k=1)],
        np.asarray([19.8, 10.5, 19.25, 19.8, 24.0, 21.75]),
        rtol=0,
        atol=1e-10,
    )


def test_seqdist_chi2_uses_attached_sequence_weights():
    actual = seqdist(fixture_b(weighted=True), method="CHI2", with_missing=True, norm="none")

    assert actual[0, 1] == pytest.approx(5.789837841, rel=0.0, abs=1e-9)


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ At\ least\ one\ substitution\ cost\ greater\ than\ twice\ the\ indel\ cost\.:traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.parametrize("deduplicate", [False, True])
def test_seqdist_trate_auto_indel_precedes_fixed_missing_cost(deduplicate):
    x = seqdef(
        np.asarray([["A", "*"], ["*", "A"], ["A", "*"]], dtype=object),
        alphabet=["A"],
        nr="*",
    )

    actual = seqdist(
        x,
        method="OM",
        subm="TRATE",
        with_missing=True,
        norm="none",
        deduplicate=deduplicate,
    )

    np.testing.assert_array_equal(actual, np.zeros((3, 3), dtype=float))
