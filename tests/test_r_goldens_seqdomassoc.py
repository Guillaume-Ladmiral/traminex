"""R-golden tests for seqdomassoc (RP4, rep.method="overall").

Anchors traminex's cross-domain LRT / Cramer's V to live TraMineR 2.2-14
``seqdomassoc`` (see scripts/r_goldens/make_goldens_seqdomassoc.R), on the shared
3-channel multidomain fixture.

The UNWEIGHTED golden is direct ``seqdomassoc``. The WEIGHTED golden is R's own
``loglin`` / Cramer's V on the *correctly* (per-case) weighted cross-table:
traminex intentionally diverges from R's ``seqdomassoc(weighted=TRUE)``, whose
column-major ``unlist`` vs case-major weight-rep assigns each cell the weight of
the wrong case (a latent R bug). See docs/API_PARITY.md and the RP3/RP4 audit.
"""

from __future__ import annotations

import csv
import inspect
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_DIFFERENCE_ATOL, R_GOLDEN_ATOL, load_golden, load_golden_lines

from traminex import seqdef, seqdomassoc

_ALPH = ["1", "2", "3"]
R2_CASES = Path(__file__).with_name("r_goldens") / "seqdomassoc_r2_cases.csv"


def _domains(weighted=False):
    w = load_golden("md_input_weights").ravel() if weighted else None
    doms = []
    for k in (1, 2, 3):
        codes = load_golden(f"md_input_s{k}").astype(int)
        chars = np.array([[_ALPH[c] for c in row] for row in codes], dtype=object)
        doms.append(seqdef(chars, alphabet=_ALPH, weights=w))
    return doms


@pytest.mark.parametrize("tag,weighted", [("uw", False), ("w", True)])
def test_seqdomassoc_overall_matches_live_r(tag, weighted):
    doms = _domains(weighted=weighted)
    out = seqdomassoc(doms, assoc=("LRT", "V"), rep_method="overall", weighted=weighted)
    assert list(out.columns) == ["df", "LRT", "p(LRT)", "v", "p(v)"]
    assert out.shape == (3, 5)
    np.testing.assert_allclose(
        out.to_numpy(dtype=float),
        load_golden(f"seqdomassoc_{tag}"),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


def _distances():
    return [load_golden(f"md_input_d{k}") for k in (1, 2, 3)]


def _missing_domains():
    positions = ((0, 0), (1, 2), (2, 4))
    domains = _domains(weighted=False)
    result = []
    for domain, (row, col) in zip(domains, positions):
        chars = np.asarray(
            [[_ALPH[code] for code in values] for values in domain.data], dtype=object
        )
        chars[row, col] = None
        result.append(seqdef(chars, alphabet=_ALPH))
    return result


def _r2_rows():
    with R2_CASES.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _case_rows(case, section="result"):
    return [row for row in _r2_rows() if row["case"] == case and row["section"] == section]


def _expected_frame(case, section="result"):
    rows = _case_rows(case, section)
    row_numbers = sorted({int(row["row"]) for row in rows})
    col_numbers = sorted({int(row["col"]) for row in rows})
    values = np.full((len(row_numbers), len(col_numbers)), np.nan)
    names = {}
    index = {}
    for row in rows:
        i = row_numbers.index(int(row["row"]))
        j = col_numbers.index(int(row["col"]))
        values[i, j] = float(row["value"])
        names[int(row["col"])] = row["text"]
        index[int(row["row"])] = row["pair"]
    return values, [names[col] for col in col_numbers], [index[row] for row in row_numbers]


_CASES = {
    "overall_uw_default_cross": dict(
        domains="uw", weighted=False, cross_table=True, dnames=["Home", "Work", "Care"]
    ),
    "overall_lrt_no_p_structfalse": dict(
        domains="uw", weighted=False, assoc=("LRT",), p_value=False, struct_zero=False
    ),
    "overall_v_no_p_custom_overlap": dict(
        domains="uw",
        weighted=False,
        assoc=("V",),
        p_value=False,
        dnames=["Life x Work", "Life", "Health"],
    ),
    "overall_wrange_cross": dict(domains="uw", weighted=True, wrange=(-1, 1), cross_table=True),
    "overall_missing_drop": dict(domains="missing", weighted=False, with_missing=False),
    "overall_missing_keep": dict(domains="missing", weighted=False, with_missing=True),
    "rep_uw_density": dict(
        domains="uw",
        weighted=False,
        rep_method="rep",
        diss_dom=True,
        seqrep_args={"coverage": 0.6, "pradius": 0.2},
        cross_table=True,
    ),
    "rep_uw_dist": dict(
        domains="uw",
        weighted=False,
        rep_method="rep",
        diss_dom=True,
        assoc=("V",),
        seqrep_args={"criterion": "dist", "nrep": 3, "pradius": 0.1},
    ),
    "eq_uw_prop": dict(
        domains="uw",
        weighted=False,
        rep_method="eq.group",
        diss_dom=True,
        seqrf_args={"k": 3, "sortv": np.arange(12)},
        cross_table=True,
    ),
    "eq_uw_first": dict(
        domains="uw",
        weighted=False,
        rep_method="eq.group",
        diss_dom=True,
        seqrf_args={"k": 4, "sortv": np.arange(12), "grp_meth": "first"},
    ),
    "eq_weighted_prop": dict(
        domains="w",
        weighted=True,
        rep_method="eq.group",
        diss_dom=True,
        seqrf_args={"k": 3, "sortv": np.arange(12)},
    ),
    "overall_weighted_correct": dict(domains="w", weighted=True),
    "rep_weighted_correct": dict(
        domains="w",
        weighted=True,
        rep_method="rep",
        diss_dom=True,
        seqrep_args={"coverage": 0.6, "pradius": 0.2},
    ),
}


def _run_case(config):
    config = dict(config)
    kind = config.pop("domains")
    domains = _missing_domains() if kind == "missing" else _domains(weighted=kind == "w")
    if config.pop("diss_dom", False):
        config["diss_dom"] = _distances()
    return seqdomassoc(domains, **config)


def test_seqdomassoc_r2_signature_is_complete():
    assert list(inspect.signature(seqdomassoc).parameters) == [
        "seqdata_dom",
        "assoc",
        "rep_method",
        "diss_dom",
        "wrange",
        "cross_table",
        "seqrep_args",
        "seqrf_args",
        "with_missing",
        "weighted",
        "struct_zero",
        "p_value",
        "dnames",
    ]


# R's windowed tables count states of different people (seqdomassoc.R:277-301);
# traminex windows within each person, a recorded correction (TX-AUD-20261011-009).
# This case asserts the definition's values, which R computes position-wise on the
# window's pairs (make_goldens_domassoc_window.R, case md_wrange), and that R's
# recorded values differ.
_WINDOW = pd.read_csv(
    Path(__file__).with_name("r_goldens") / "domassoc_window.csv", dtype=str, keep_default_na=False
)


def _window_definition(section: str, pair: str) -> dict[str, float]:
    rows = _WINDOW[
        (_WINDOW["case"] == "md_wrange")
        & (_WINDOW["kind"] == "definition")
        & (_WINDOW["section"] == section)
        & (_WINDOW["pair"] == pair)
    ]
    return {name: float(value) for name, value in zip(rows["name"], rows["value"])}


@pytest.mark.parametrize("case", list(_CASES))
def test_seqdomassoc_r2_routes_and_controls_match_oracle(case):
    assert {row["case"] for row in _r2_rows()} == set(_CASES)
    out = _run_case(_CASES[case])
    expected, columns, index = _expected_frame(case)
    assert list(out.columns) == columns
    assert list(out.index) == index
    if case == "overall_wrange_cross":
        definition = np.array(
            [[_window_definition("stat", pair)[column] for column in columns] for pair in index]
        )
        np.testing.assert_allclose(out.to_numpy(dtype=float), definition, rtol=0, atol=1e-10)
        assert not np.allclose(expected, definition, rtol=0, atol=R_DIFFERENCE_ATOL)
        return
    np.testing.assert_allclose(
        out.to_numpy(dtype=float), expected, rtol=0, atol=R_GOLDEN_ATOL, equal_nan=True
    )

    upstream = _case_rows(case, "upstream_result")
    if upstream:
        upstream_values, _, _ = _expected_frame(case, "upstream_result")
        assert not np.allclose(expected, upstream_values, rtol=0, atol=R_DIFFERENCE_ATOL)


@pytest.mark.parametrize(
    "case", ["overall_uw_default_cross", "overall_wrange_cross", "rep_uw_density", "eq_uw_prop"]
)
def test_seqdomassoc_r2_cross_tables_match_r(case):
    out = _run_case(_CASES[case])
    expected_rows = _case_rows(case, "cross_table")
    assert "cross_tables" in out.attrs
    assert set(out.attrs["cross_tables"]) == {row["pair"] for row in expected_rows}
    for name, table in out.attrs["cross_tables"].items():
        rows = [row for row in expected_rows if row["pair"] == name]
        expected = np.full(
            (max(int(row["row"]) for row in rows), max(int(row["col"]) for row in rows)),
            np.nan,
        )
        row_names = {}
        col_names = {}
        for row in rows:
            i, j = int(row["row"]) - 1, int(row["col"]) - 1
            expected[i, j] = float(row["value"])
            row_name, col_name = row["text"].split("|", 1)
            row_names[i] = row_name
            col_names[j] = col_name
        if case == "overall_wrange_cross":
            definition = _window_definition("table", name)
            ours = {f"{r}|{c}": float(table.loc[r, c]) for r in table.index for c in table.columns}
            assert ours == definition
            assert not np.allclose(table.to_numpy(), expected, atol=R_DIFFERENCE_ATOL, rtol=0)
        else:
            np.testing.assert_allclose(table.to_numpy(), expected, atol=R_GOLDEN_ATOL, rtol=0)
        assert list(table.index) == [row_names[i] for i in range(len(row_names))]
        assert list(table.columns) == [col_names[j] for j in range(len(col_names))]


def test_seqdomassoc_routes_through_public_representative_primitives(monkeypatch):
    import traminex.cluster as cluster_module
    import traminex.repseq_rf as rf_module

    seqrep_calls = []
    seqrf_calls = []
    original_seqrep = cluster_module.seqrep
    original_seqrf = rf_module.seqrf

    def recording_seqrep(*args, **kwargs):
        seqrep_calls.append(kwargs.copy())
        return original_seqrep(*args, **kwargs)

    def recording_seqrf(*args, **kwargs):
        seqrf_calls.append(kwargs.copy())
        return original_seqrf(*args, **kwargs)

    monkeypatch.setattr(cluster_module, "seqrep", recording_seqrep)
    monkeypatch.setattr(rf_module, "seqrf", recording_seqrf)
    _run_case(_CASES["rep_uw_density"])
    _run_case(_CASES["eq_uw_first"])
    assert len(seqrep_calls) == 3
    assert all(call["coverage"] == 0.6 and call["pradius"] == 0.2 for call in seqrep_calls)
    assert len(seqrf_calls) == 3
    assert all(call["grp_meth"] == "first" and call["k"] == 4 for call in seqrf_calls)


def test_seqdomassoc_r2_refusals_are_frozen_and_python_is_actionable():
    assert load_golden_lines("seqdomassoc_r2_errors") == [
        "not_list=ERROR:invalid 'times' argument",
        "one_domain=ERROR:seqdata.dom should be a list of at least two state sequence objects",
        "rep_missing_diss=ERROR: Except for rep.method='overall', diss.dom must be a list of dissimilarity matrices",
        "eq_diss_wrong_length=ERROR:diss.dom should be a list of same length as seqdata.dom",
        "bad_assoc=ERROR:bad assoc values, allowed are LRT,V",
        "bad_rep_method=ERROR: Unknown rep.method: bogus",
        "reversed_wrange=ERROR:wrange[1] cannot be larger than wrange[2]",
        "incompatible_ncases=ERROR:variable lengths differ (found for 'repseq[[d2]]')",
        "incompatible_width=ERROR:variable lengths differ (found for 'repseq[[d2]]')",
        "seqrep_bad_args=ERROR:Unknown criterion / no score provided",
    ]
    domains = _domains()
    with pytest.raises(ValueError, match="at least two"):
        seqdomassoc(domains[:1])
    with pytest.raises(ValueError, match="diss_dom"):
        seqdomassoc(domains, rep_method="rep")
    with pytest.raises(ValueError, match="same length"):
        seqdomassoc(domains, rep_method="eq.group", diss_dom=_distances()[:2])
    with pytest.raises(ValueError, match="assoc"):
        seqdomassoc(domains, assoc=("tau",))
    with pytest.raises(ValueError, match="rep_method"):
        seqdomassoc(domains, rep_method="bogus", diss_dom=_distances())
    with pytest.raises(ValueError, match="wrange"):
        seqdomassoc(domains, wrange=(1, -1))

    short = _domains()
    short[1] = seqdef(
        np.asarray([[_ALPH[c] for c in row] for row in short[1].data[:-1]], dtype=object),
        alphabet=_ALPH,
    )
    with pytest.raises(ValueError, match="same number of cases"):
        seqdomassoc(short, weighted=False)
