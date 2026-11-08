"""Complete TraMineR 2.2-14 R1 contract matrix for ``seqdist``.

The R generator records every official method on both an interior-missing and
a ragged fixture.  It also freezes normalization, output-shape, reference,
state-dependent-indel, and position-dependent-cost contracts.  Every manifest
case is executed: a case either matches R numerically or raises the same
normalized refusal category/message.

Regenerate with::

    /usr/bin/Rscript scripts/r_goldens/make_goldens_distances_void.R
"""

from __future__ import annotations

import ast
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from corestats_fixtures import ALPH_C, SEQS_C, fixture_b, fixture_c
from r_goldens_utils import R_DIFFERENCE_ATOL, R_GOLDEN_ATOL

from traminex import seqdist
from traminex.format import seqdef

GOLDEN_DIR = Path(__file__).resolve().parent / "r_goldens"

OFFICIAL_METHODS = {
    "OM",
    "OMloc",
    "OMslen",
    "OMspell",
    "OMstran",
    "HAM",
    "DHD",
    "CHI2",
    "EUCLID",
    "LCS",
    "LCP",
    "RLCP",
    "NMS",
    "NMSMST",
    "SVRspell",
    "TWED",
}
ALL_NORMS = {"none", "auto", "maxlength", "gmean", "maxdist", "YujianBo"}


def fixture_f():
    """Ragged fixture whose second sequence is uniquely most frequent."""
    rows = [list(row) for row in SEQS_C]
    rows.append(list(SEQS_C[1]))
    return seqdef(np.asarray(rows, dtype=object), alphabet=ALPH_C, nr="*")


def external_ref():
    """Separately constructed shorter reference with retained missing + void."""
    return seqdef(
        np.asarray([["A", "*", "C", "C", "%", "%"]], dtype=object),
        alphabet=ALPH_C,
        nr="*",
    )


def fixture_g():
    """Complete source whose external reference alone materializes missing."""
    return seqdef(
        np.asarray([["A", "A", "A"], ["A", "B", "A"]], dtype=object),
        alphabet=["A", "B"],
        nr="*",
    )


def asymmetric_external_ref():
    """External reference with one retained interior missing state."""
    return seqdef(
        np.asarray([["A", "*", "A"]], dtype=object),
        alphabet=["A", "B"],
        nr="*",
    )


def fixture_h():
    """Frequency tie whose R numeric-SPS order differs from input order."""
    return seqdef(
        np.asarray([["B", "B"], ["A", "A"], ["C", "C"]], dtype=object),
        alphabet=["A", "B", "C"],
    )


def fixture_j():
    """Two-state fixture for asymmetric substitution-cost normalization."""
    return seqdef(np.asarray([["A"], ["B"]], dtype=object), alphabet=["A", "B"])


def fixture_k():
    """One structurally empty row plus two non-empty transition sequences."""
    return seqdef(
        np.asarray(
            [["%", "%", "%"], ["A", "B", "A"], ["A", "B", "B"]],
            dtype=object,
        ),
        alphabet=["A", "B"],
        nr="*",
    )


def fixture_l():
    """TRATE fixture whose automatic indel is zero before missing-cost fixup."""
    return seqdef(
        np.asarray([["A", "*"], ["*", "A"]], dtype=object),
        alphabet=["A"],
        nr="*",
    )


def fixture_m():
    """Two structural-empty rows that form the most frequent sequence."""
    return seqdef(
        np.asarray(
            [
                ["%", "%", "%"],
                ["%", "%", "%"],
                ["A", "*", "A"],
                ["A", "B", "%"],
            ],
            dtype=object,
        ),
        alphabet=["A", "B"],
        nr="*",
    )


def fixture_n():
    """Two-state fixture for state-dependent indels."""
    return seqdef(
        np.asarray([["A", "A", "B"], ["A", "B", "B"]], dtype=object),
        alphabet=["A", "B"],
    )


def fixture_o():
    """DHD fixture whose declared C state is unobserved."""
    return seqdef(
        np.asarray([["A", "A", "A"], ["A", "B", "A"]], dtype=object),
        alphabet=["A", "B", "C"],
    )


def fixture_p():
    """Equal logical-length DHD rows with interior missing and trailing void."""
    return seqdef(
        np.asarray([["A", "*", "A", "%"], ["B", "A", "B", "%"]], dtype=object),
        alphabet=["A", "B"],
        nr="*",
    )


def short_ref():
    """External reference shorter than the complete source rows."""
    return seqdef(
        np.asarray([["A", "B"]], dtype=object),
        alphabet=["A", "B"],
        nr="*",
    )


def fixture_bw():
    return fixture_b(weighted=True)


def fixture_cw():
    return fixture_c(weighted=True)


def asymmetric_subm() -> np.ndarray:
    return np.asarray([[0.0, 1.0], [10.0, 0.0]])


def soft_prox() -> np.ndarray:
    return np.asarray([[1.0, 0.5], [0.5, 1.0]])


def position_subm() -> np.ndarray:
    """DHD costs: off-diagonal cost ``t + 1`` at zero-based position ``t``."""
    out = np.zeros((6, 4, 4), dtype=float)
    for time in range(6):
        out[time, :, :] = float(time + 1)
        np.fill_diagonal(out[time], 0.0)
    return out


def static_subm() -> np.ndarray:
    """A static matrix, which R-facing DHD explicitly refuses."""
    return np.asarray([[0.0, 2.0], [2.0, 0.0]])


def wrong_time_subm() -> np.ndarray:
    """Four positions for a source whose original stored width is three."""
    out = np.full((4, 2, 2), 2.0, dtype=float)
    for position in range(4):
        np.fill_diagonal(out[position], 0.0)
    return out


def ragged_position_subm() -> np.ndarray:
    """Stored-width DHD costs for a three-state, four-column ragged source."""
    out = np.zeros((4, 3, 3), dtype=float)
    for position in range(4):
        out[position, :, :] = float(position + 1)
        np.fill_diagonal(out[position], 0.0)
    return out


def ragged_canonical_subm() -> np.ndarray:
    """Invalid logical-width costs that omit the stored structural-void column."""
    return ragged_position_subm()[:3]


def negative_dhd_subm() -> np.ndarray:
    out = np.broadcast_to(static_subm(), (3, 2, 2)).copy()
    out[0, 0, 1] = out[0, 1, 0] = -1.0
    return out


def diagonal_dhd_subm() -> np.ndarray:
    out = np.broadcast_to(static_subm(), (3, 2, 2)).copy()
    out[0, 0, 0] = 1.0
    return out


def nan_dhd_subm() -> np.ndarray:
    out = np.broadcast_to(static_subm(), (3, 2, 2)).copy()
    out[0, 0, 1] = out[0, 1, 0] = np.nan
    return out


def infinite_dhd_subm() -> np.ndarray:
    out = np.broadcast_to(static_subm(), (3, 2, 2)).copy()
    out[0, 0, 1] = out[0, 1, 0] = np.inf
    return out


def nan_value() -> float:
    return float("nan")


def infinite_value() -> float:
    return float("inf")


_EXECUTION_NAMESPACE = {
    "seqdist": seqdist,
    "fixture_b": fixture_b,
    "fixture_c": fixture_c,
    "fixture_f": fixture_f,
    "fixture_g": fixture_g,
    "fixture_h": fixture_h,
    "fixture_j": fixture_j,
    "fixture_k": fixture_k,
    "fixture_l": fixture_l,
    "fixture_m": fixture_m,
    "fixture_n": fixture_n,
    "fixture_o": fixture_o,
    "fixture_p": fixture_p,
    "fixture_bw": fixture_bw,
    "fixture_cw": fixture_cw,
    "external_ref": external_ref,
    "asymmetric_external_ref": asymmetric_external_ref,
    "short_ref": short_ref,
    "position_subm": position_subm,
    "static_subm": static_subm,
    "wrong_time_subm": wrong_time_subm,
    "ragged_position_subm": ragged_position_subm,
    "ragged_canonical_subm": ragged_canonical_subm,
    "negative_dhd_subm": negative_dhd_subm,
    "diagonal_dhd_subm": diagonal_dhd_subm,
    "nan_dhd_subm": nan_dhd_subm,
    "infinite_dhd_subm": infinite_dhd_subm,
    "nan_value": nan_value,
    "infinite_value": infinite_value,
    "asymmetric_subm": asymmetric_subm,
    "soft_prox": soft_prox,
}
_FIXTURE_CALLS = {
    "B": "fixture_b",
    "C": "fixture_c",
    "F": "fixture_f",
    "G": "fixture_g",
    "H": "fixture_h",
    "J": "fixture_j",
    "K": "fixture_k",
    "L": "fixture_l",
    "M": "fixture_m",
    "N": "fixture_n",
    "O": "fixture_o",
    "P": "fixture_p",
    "W": "fixture_bw",
    "CW": "fixture_cw",
}
_ALLOWED_AST_NODES = (
    ast.Expression,
    ast.Call,
    ast.Name,
    ast.Load,
    ast.keyword,
    ast.Constant,
    ast.List,
    ast.Tuple,
    ast.UnaryOp,
    ast.USub,
)


def _validated_manifest_call(row: pd.Series) -> ast.Expression:
    """Parse one checked-in call and reject anything beyond the tiny call DSL."""
    tree = ast.parse(row["python_call"], mode="eval")
    assert isinstance(tree, ast.Expression)
    for node in ast.walk(tree):
        assert isinstance(node, _ALLOWED_AST_NODES), ast.dump(node)
        if isinstance(node, ast.Call):
            assert isinstance(node.func, ast.Name)
            assert node.func.id in _EXECUTION_NAMESPACE
            assert all(keyword.arg is not None for keyword in node.keywords)

    call = tree.body
    assert isinstance(call, ast.Call)
    assert isinstance(call.func, ast.Name) and call.func.id == "seqdist"
    assert len(call.args) == 1
    fixture_call = call.args[0]
    assert isinstance(fixture_call, ast.Call)
    assert isinstance(fixture_call.func, ast.Name)
    assert fixture_call.func.id == _FIXTURE_CALLS[row["fixture"]]
    assert not fixture_call.args and not fixture_call.keywords

    keywords = {keyword.arg: keyword.value for keyword in call.keywords}
    assert len(keywords) == len(call.keywords)
    method = keywords.get("method")
    assert isinstance(method, ast.Constant) and method.value == row["method"]
    return tree


def _execute_manifest_call(row: pd.Series):
    tree = _validated_manifest_call(row)
    return eval(  # noqa: S307 - AST and globals are deliberately restricted above.
        compile(tree, "<distance-r1-manifest>", "eval"),
        {"__builtins__": {}, **_EXECUTION_NAMESPACE},
        {},
    )


NUMERIC = pd.read_csv(
    GOLDEN_DIR / "distances_r1_cases.csv",
    dtype=str,
    keep_default_na=False,
)
ERRORS = pd.read_csv(
    GOLDEN_DIR / "distances_r1_errors.txt",
    sep="|",
    dtype=str,
    keep_default_na=False,
)
_MANIFEST_COLUMNS = ["case", "method", "fixture", "options", "python_call"]
MANIFEST = pd.concat(
    [NUMERIC[_MANIFEST_COLUMNS], ERRORS[_MANIFEST_COLUMNS]],
    ignore_index=True,
).drop_duplicates()
CASE_IDS = tuple(MANIFEST["case"])


def _numeric(token: str) -> float:
    lowered = token.strip().lower()
    if lowered in {"na", "nan"}:
        return float("nan")
    if lowered in {"inf", "+inf"}:
        return float("inf")
    if lowered == "-inf":
        return float("-inf")
    return float(token)


def _single(group: pd.DataFrame, column: str) -> str:
    values = group[column].unique().tolist()
    assert len(values) == 1, f"inconsistent {column}: {values}"
    return values[0]


def _expected_numeric(case_id: str) -> tuple[str, np.ndarray]:
    group = NUMERIC.loc[NUMERIC["case"] == case_id].copy()
    assert not group.empty
    kind = _single(group, "output_kind")
    shape_0 = int(_single(group, "shape_0"))
    shape_1_token = _single(group, "shape_1")

    group["index_i"] = group["index"].astype(int)
    group = group.sort_values("index_i")
    assert group["index_i"].tolist() == list(range(len(group)))
    values = np.asarray([_numeric(value) for value in group["value"]], dtype=float)

    if kind == "matrix":
        shape_1 = int(shape_1_token)
        assert len(group) == shape_0 * shape_1
        assert group["row"].astype(int).tolist() == [
            index // shape_1 for index in range(len(group))
        ]
        assert group["col"].astype(int).tolist() == [index % shape_1 for index in range(len(group))]
        return kind, values.reshape(shape_0, shape_1)

    assert kind in {"vector", "condensed"}
    assert shape_1_token == ""
    assert len(group) == shape_0
    assert set(group["row"]) == {""}
    assert set(group["col"]) == {""}
    return kind, values


def _normalize_exception(error: BaseException, method: str) -> tuple[str, str]:
    message = " ".join(str(error).split())
    lowered = message.lower()

    if "indel" in lowered and "finite" in lowered and "positive" in lowered:
        return "invalid_indel", "indel cost must be finite and strictly positive"
    if "indel" in lowered and "positive" in lowered:
        return "invalid_indel", "indel cost must be strictly positive"
    if method == "DHD" and (
        "substitution costs must be finite" in lowered
        or "substitution costs cannot be negative" in lowered
        or "substitution-cost diagonals must be zero" in lowered
    ):
        return (
            "invalid_dhd_costs",
            "DHD substitution costs must be finite, nonnegative, and zero-diagonal",
        )
    if "substitution costs must be finite" in lowered:
        return "invalid_subm", "substitution costs must be finite"
    if method == "DHD" and (
        "dhd subm" in lowered
        or "dhd substitution cost" in lowered
        or "time-varying array" in lowered
    ):
        return (
            "invalid_dhd_subm",
            "DHD subm must match the active states and original stored width",
        )

    if (
        ("with_missing" in lowered or "with.missing" in lowered)
        and "missing" in lowered
        and ("true" in lowered or "opt" in lowered)
    ):
        return (
            "missing_requires_opt_in",
            "with_missing must be true when sequences contain missing values",
        )
    if (
        "different length" in lowered
        or "different lengths" in lowered
        or ("equal length" in lowered and method in {"HAM", "DHD"})
    ):
        return (
            "unequal_length",
            f"{method} is not defined for sequences of different length",
        )
    if "norm" in lowered and any(
        marker in lowered
        for marker in ("not implemented", "isn't implemented", "not supported", "only")
    ):
        return "unsupported_norm", f"norm is not supported for method {method}"
    return "unclassified", message


def test_manifest_schema_and_case_map_are_exact() -> None:
    assert NUMERIC.columns.tolist() == [
        "case",
        "method",
        "fixture",
        "options",
        "python_call",
        "output_kind",
        "shape_0",
        "shape_1",
        "index",
        "row",
        "col",
        "value",
    ]
    assert ERRORS.columns.tolist() == [
        "case",
        "method",
        "fixture",
        "options",
        "python_call",
        "category",
        "message",
    ]

    numeric_cases = set(NUMERIC["case"])
    error_cases = set(ERRORS["case"])
    assert numeric_cases.isdisjoint(error_cases)
    assert numeric_cases | error_cases == set(CASE_IDS)
    assert len(CASE_IDS) == 172
    assert MANIFEST["case"].is_unique
    assert set(MANIFEST["method"]) == OFFICIAL_METHODS
    assert set(MANIFEST["fixture"]) == {
        "B",
        "C",
        "CW",
        "F",
        "G",
        "H",
        "J",
        "K",
        "L",
        "M",
        "N",
        "O",
        "P",
        "W",
    }
    assert MANIFEST["python_call"].str.startswith("seqdist(").all()
    for _, row in MANIFEST.iterrows():
        _validated_manifest_call(row)

    for method in OFFICIAL_METHODS:
        slug = method.lower()
        assert f"interior_{slug}" in CASE_IDS
        assert f"ragged_{slug}" in CASE_IDS

    observed_norms = {
        option.split(";norm=", 1)[1].split(";", 1)[0]
        for option in MANIFEST["options"]
        if ";norm=" in option
    }
    assert observed_norms == ALL_NORMS
    assert {
        "ragged_lcs_full_matrix_false",
        "ragged_lcs_condensed_alias",
        "ragged_lcs_refseq_sets",
        "ragged_lcs_ref_sets",
        "ragged_om_refseq_single",
        "ragged_lcp_refseq_most_frequent",
        "ragged_lcs_refseq_external",
        "ragged_om_vector_indel",
        "interior_dhd_position_subm",
        "interior_lcs_missing_refused",
        "asymmetric_om_maxdist",
        "empty_twed_auto_indel_border",
        "ragged_om_trate_refseq_external",
        "interior_om_trate_zero_pre_missing_indel",
        "interior_dhd_attached_weighted",
        "interior_chi2_attached_weighted",
        "ragged_omstran_prob_previous_true_add_column_true",
        "ragged_omstran_prob_refseq_external",
        "refseq_most_frequent_tie_numeric_sps",
        "refseq_most_frequent_structural_empty",
        "dhd_declared_unobserved_state_auto",
        "dhd_ragged_missing_none",
        "dhd_ragged_missing_auto",
        "chi2_ragged_empty_interval_auto",
        "euclid_ragged_empty_interval_auto",
        "complete_ham_refseq_shorter",
        "complete_dhd_refseq_shorter",
        "vector_indel_omslen",
        "vector_indel_omspell",
        "vector_indel_twed",
        "explicit_sm_auto_indel_om",
        "explicit_sm_auto_indel_omslen",
        "explicit_sm_auto_indel_omspell",
        "explicit_sm_auto_indel_omstran",
        "interior_dhd_trate_method",
        "interior_dhd_indels_method",
        "interior_dhd_indelslog_method",
        "interior_dhd_indels_attached_unweighted",
        "complete_dhd_indels_refseq_shorter",
        "refuse_dhd_static_subm",
        "refuse_dhd_wrong_time_subm",
        "refuse_dhd_wrong_time_subm_refseq_shorter",
        "refuse_nonpositive_indel_omslen",
        "refuse_nonpositive_indel_omspell",
        "refuse_nonpositive_indel_omstran",
        "refuse_nonpositive_indel_twed",
        "refuse_zero_auto_indel_om",
        "dhd_ragged_position_subm_stored_width",
        "refuse_dhd_canonical_width_subm_on_ragged",
        "refuse_dhd_negative_subm",
        "refuse_dhd_diagonal_subm",
        "refuse_dhd_nan_subm",
        "refuse_dhd_infinite_subm",
        "refuse_nan_indel_om",
        "refuse_nan_indel_twed",
        "accept_infinite_indel_om",
        "refuse_om_nan_subm",
        "refuse_om_infinite_subm",
        "infinite_indel_scalar_omslen",
        "infinite_indel_scalar_omspell",
        "refuse_infinite_indel_scalar_omstran",
        "infinite_indel_scalar_twed",
        "infinite_indel_vector_om",
        "infinite_indel_vector_omslen",
        "infinite_indel_vector_omspell",
        "refuse_infinite_indel_vector_omstran",
        "infinite_indel_vector_twed",
    } <= set(CASE_IDS)


#: Degenerate inputs that TraMineR 2.2-14 accepts and traminex deliberately
#: refuses. The goldens record R's actual behaviour -- the oracle artifact must
#: stay truthful -- and the divergence is asserted here instead.
#:
#: 2.2-13 refused all three. 2.2-14's rewritten ``checkcost`` demoted a zero
#: indel from an error to a warning and added division-by-zero guards to
#: ``OMstran``, so the calls now return values. For the two ``OMstran`` cases
#: those values are **NaN** off the diagonal, which is exactly the "silently
#: meaningless matrix" that ``tests/test_seqdist_infinite_cost_contract.py``
#: records traminex as refusing to return. ``refuse_zero_auto_indel_om`` yields
#: an all-zero matrix, which is arithmetically right but carries no information
#: at all; traminex keeps naming the cause rather than returning it.
#:
#: Recorded in ``docs/API_PARITY.md`` under the 2.2-14 divergences.
DELIBERATE_REFUSALS_OF_R_ACCEPTED_INPUT = {
    "refuse_infinite_indel_scalar_omstran": (
        "invalid_indel",
        "indel cost must be finite and strictly positive",
    ),
    "refuse_infinite_indel_vector_omstran": (
        "invalid_indel",
        "indel cost must be finite and strictly positive",
    ),
    "refuse_zero_auto_indel_om": (
        "invalid_indel",
        "indel cost must be strictly positive",
    ),
}


def test_every_deliberate_refusal_is_a_case_r_actually_accepts() -> None:
    """Guard the exception list: each entry must be a *numeric* golden case.

    If a future oracle goes back to refusing one of these, its row moves to the
    errors file and this test fails, forcing the divergence to be re-examined
    rather than silently kept.
    """
    numeric_cases = set(NUMERIC["case"])
    stale = sorted(set(DELIBERATE_REFUSALS_OF_R_ACCEPTED_INPUT) - numeric_cases)
    assert not stale, f"no longer accepted by the oracle, so the exception is stale: {stale}"


#: Cases where R's value is an R defect that traminex corrects, with the value derived from
#: the definition. fixture_j: d(A, B) = 1 under costs [[0, 1], [10, 0]] and indel 10; the
#: maximum possible distance is min(10, 2 * 10) * 1 = 10, not R's upper-triangle 1, so maxdist
#: is 0.1 and gmean 1 - (10 - 1) / (2 * sqrt(10 * 10)) = 0.55 (R: 1 and 1).
R_DEFECT_CORRECTED = {
    "asymmetric_om_maxdist": np.array([[0.0, 0.1], [0.1, 0.0]]),
    "asymmetric_om_gmean": np.array([[0.0, 0.55], [0.55, 0.0]]),
    # OMloc YujianBo (cost 2, expcost 0.25, context 0.5; lengths 6, 3, 4, 2): one local indel
    # costs up to 0.25 * 2 + 0.5 * 2 = 1.5, so the bound is |n - m| * 1.5 + 3 * min(m, n), not
    # R's |n - m| + 2 * min(m, n); R's raw distances (5.5 for lengths 6 and 3) over it
    # (TX-REP-20261021-003). R's values reach 1.077.
    "auto_omloc_ragged": np.array(
        [
            [0.0, 11 / 19, 4 / 7, 10 / 17],
            [11 / 19, 0.0, 4 / 9, 11 / 13],
            [4 / 7, 4 / 9, 0.0, 7 / 8],
            [10 / 17, 11 / 13, 7 / 8, 0.0],
        ]
    ),
}


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(r"ignore:\s*Non\-symmetric\ substitution\ cost\ matrix!:UserWarning")
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ indel\ vector\ not\ supported\ by\ the\ chosen\ method,\ max\(indel\)\ used\ instead!:UserWarning"
)
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ At\ least\ one\ substitution\ cost\ greater\ than\ twice\ the\ indel\ cost\.:traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ Empty\ sequence\(s\)\ \ \d+,\ \d+\ \ may\ cause\ inconsistent\ results!:UserWarning"
)
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ Empty\ sequence\(s\)\ \ \d+\ \ may\ cause\ inconsistent\ results!:UserWarning"
)
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ \ at\ least,\ one\ substitution\ cost\ doesn't\ respect\ the\ triangle\ inequality\.:traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ seqdist.\ 'with\.missing'\ set\ as\ FALSE\ as\ 'seqdata'\ has\ no\ non\-void\ missing\ values:UserWarning"
)
@pytest.mark.parametrize("case_id", CASE_IDS)
def test_seqdist_r1_case_matches_r(case_id: str) -> None:
    row = MANIFEST.loc[MANIFEST["case"] == case_id].iloc[0]
    method = row["method"]
    expected_error = ERRORS.loc[ERRORS["case"] == case_id]

    deliberate = DELIBERATE_REFUSALS_OF_R_ACCEPTED_INPUT.get(case_id)
    if deliberate is not None:
        with pytest.raises(Exception) as captured:
            _execute_manifest_call(row)
        assert _normalize_exception(captured.value, method) == deliberate
        return

    if not expected_error.empty:
        assert len(expected_error) == 1
        with pytest.raises(Exception) as captured:
            _execute_manifest_call(row)
        category, message = _normalize_exception(captured.value, method)
        assert category == expected_error.iloc[0]["category"]
        assert message == expected_error.iloc[0]["message"]
        return

    kind, expected = _expected_numeric(case_id)
    actual = np.asarray(_execute_manifest_call(row), dtype=float)
    corrected = R_DEFECT_CORRECTED.get(case_id)
    if corrected is not None:
        # R's value stays recorded; traminex asserts the definition's (TX-AUD-20261019-008,
        # TX-REP-20261021-003).
        np.testing.assert_allclose(actual, corrected, rtol=0, atol=1e-12)
        assert not np.allclose(actual, expected, rtol=0, atol=R_DIFFERENCE_ATOL)
        return
    if kind == "matrix":
        assert actual.ndim == 2
    else:
        assert actual.ndim == 1
    assert actual.shape == expected.shape
    np.testing.assert_allclose(actual, expected, rtol=0, atol=R_GOLDEN_ATOL, equal_nan=True)
