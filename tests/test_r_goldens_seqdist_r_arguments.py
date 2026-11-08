"""seqdist accepts R's arguments for every R method, as R does (TX-AUD-20261010-071).

R's seqdist takes all of its arguments for every method (seqdist.R:4-11): it
checks some for every method (norm, indel, sm, prox, link, step), some only for
the method that uses them, and ignores an argument the method does not use, even
at a value it would refuse for a method that uses it. traminex scoped R's
arguments per method and refused the others: 251 of 378 single-argument calls
with an R default given explicitly, and every call giving all of them. Fixture
``scripts/r_goldens/make_goldens_seqdist_r_arguments.R``; this module rebuilds
each of its calls with traminex's spelling of R's argument names.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "seqdist_r_arguments.csv",
    dtype=str,
    keep_default_na=False,
).set_index("case")
SX = tx.seqdef(
    np.array(
        [["A", "B", "A", "C"], ["B", "B", "C", "A"], ["A", "A", "B", "B"], ["C", "C", "C", "A"]],
        dtype=object,
    ),
    alphabet=["A", "B", "C"],
)
METHODS = ["OM", "OMloc", "OMslen", "OMspell", "OMstran", "HAM", "DHD", "CHI2", "EUCLID", "LCS",
           "LCP", "RLCP", "NMS", "NMSMST", "SVRspell", "TWED"]  # fmt: skip
REQUIRED = {
    "OM": {"subm": "CONSTANT"},
    "OMloc": {"subm": "CONSTANT"},
    "OMslen": {"subm": "CONSTANT"},
    "OMspell": {"subm": "CONSTANT"},
    "OMstran": {"subm": "CONSTANT", "otto": 0.5},
    "TWED": {"subm": "CONSTANT", "nu": 0.5},
}
#: R's argument names in traminex's spelling (R's ``sm`` is ``subm``).
NAME = {"sm": "subm", "with.missing": "with_missing", "full.matrix": "full_matrix",
        "add.column": "add_column", "global.pdotj": "global_pdotj", "check.max.size": "check_max_size",
        "check.costs": "check_costs", "opt.args": "opt_args"}  # fmt: skip
USES = {
    "indel": {"OM", "OMslen", "OMspell", "OMstran", "TWED"},
    "sm": {"OM", "OMloc", "OMslen", "OMspell", "OMstran", "HAM", "DHD", "TWED", "NMS", "SVRspell"},
    "kweights": {"NMS", "NMSMST", "SVRspell"},
    "tpow": {"OMspell", "NMSMST", "SVRspell"},
    "expcost": {"OMloc", "OMspell"},
    "context": {"OMloc"},
    "link": {"OMslen"},
    "h": {"OMslen", "TWED"},
    "nu": {"TWED"},
    "transindel": {"OMstran"},
    "otto": {"OMstran"},
    "previous": {"OMstran"},
    "add.column": {"OMstran"},
    "breaks": {"CHI2", "EUCLID"},
    "step": {"CHI2", "EUCLID"},
    "overlap": {"CHI2", "EUCLID"},
    "global.pdotj": {"CHI2", "EUCLID"},
    "prox": {"NMS", "SVRspell"},
    "opt.args": {"OMspell"},
}
DEFAULTS = {"refseq": None, "norm": "none", "indel": "auto", "sm": None, "with.missing": False,
            "full.matrix": True, "kweights": np.ones(4), "tpow": 1.0, "expcost": 0.5, "link": "mean",
            "h": 0.5, "transindel": "constant", "previous": False, "add.column": True, "breaks": None,
            "step": 1, "overlap": False, "weighted": True, "global.pdotj": None, "prox": None,
            "check.max.size": True, "check.costs": True, "opt.args": {}}  # fmt: skip
IGNORED = {"indel": 3, "sm": "INDELS", "kweights": np.full(4, 2.0), "tpow": 2, "expcost": 0.25,
           "context": 0.3, "link": "gmean", "h": 1, "nu": 0.7, "transindel": "subcost", "otto": 0.3,
           "previous": True, "add.column": False, "breaks": [(1, 2)], "step": 2, "overlap": True,
           "global.pdotj": np.array([0.5, 0.25, 0.25]), "prox": np.eye(3),
           "opt.args": {"tokdep.coeff": [1, 1, 1]}, "check.max.size": False}  # fmt: skip
INVALID_IGNORED = {"expcost": -1, "h": -1, "nu": 0, "otto": 2, "kweights": -1, "tpow": -1}
GLOBAL_STOPS = {
    "global_link_bogus": ("OM", {"subm": "CONSTANT", "link": "bogus"}),
    "global_step_fraction": ("LCS", {"step": 2.5}),
    "global_prox_not_matrix": ("OM", {"subm": "CONSTANT", "prox": 1}),
    "global_sm_bogus": ("CHI2", {"subm": "bogus"}),
    "global_indel_bogus": ("HAM", {"indel": "bogus"}),
    "global_indel_vector_length": ("HAM", {"indel": [1, 2]}),
    "nms_sm_given": ("NMS", {"subm": "CONSTANT"}),
    "svrspell_sm_given": ("SVRspell", {"subm": "CONSTANT"}),
    "nms_norm_auto": ("NMS", {"norm": "auto"}),
    "nmsmst_norm_auto": ("NMSMST", {"norm": "auto"}),
}
#: R's sm is NULL, a matrix, an array or a method name for every method; opt.args
#: counts only as a list with a tokdep.coeff entry (batch-6 review M-2, M-3).
ARGUMENT_FORMS = {
    "global_sm_number": ("CHI2", {"subm": 2}),
    "global_sm_number_lcs": ("LCS", {"subm": 2}),
    "global_sm_vector": ("LCP", {"subm": [1, 2]}),
    "omspell_opt_args_number": ("OMspell", {"subm": "CONSTANT", "opt_args": 5}),
    "omspell_opt_args_unnamed_list": (
        "OMspell",
        {"subm": "CONSTANT", "opt_args": ["tokdep.coeff"]},
    ),
    "omspell_opt_args_null_coeff": (
        "OMspell",
        {"subm": "CONSTANT", "opt_args": {"tokdep.coeff": None}},
    ),
}


def _cases() -> dict[str, tuple[str, dict]]:
    cases: dict[str, tuple[str, dict]] = {}
    for method in METHODS:
        base = dict(REQUIRED.get(method, {}))
        cases[f"{method}_base"] = (method, dict(base))
        every = dict(base)
        for name, value in DEFAULTS.items():
            if not (name == "sm" and "subm" in base):
                every[NAME.get(name, name)] = value
        cases[f"{method}_all_defaults"] = (method, every)
        for prefix, table in (("ignored", IGNORED), ("ignored_invalid", INVALID_IGNORED)):
            for name, value in table.items():
                if method in USES.get(name, set()) or NAME.get(name, name) in base:
                    continue
                cases[f"{method}_{prefix}_{name}"] = (method, {**base, NAME.get(name, name): value})
    cases.update(GLOBAL_STOPS)
    cases.update(ARGUMENT_FORMS)
    return cases


CASES = _cases()


def test_every_r_case_is_rebuilt() -> None:
    assert set(CASES) == set(GOLDEN.index)


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ \ \ \ \ \ for\ example,\ cost\ is\ zero\ for\ the\ pair\ \(\ A\ ,\ C\ \):traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ \[!\]\ At\ pos\ \d+,\ \d+\ substitution\ cost\ value\(s\)\ is\ \d+!:traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.parametrize("case", sorted(CASES))
def test_every_case_matches_r(case: str) -> None:
    method, kwargs = CASES[case]
    row = GOLDEN.loc[case]
    if row["ok"] == "FALSE":
        # Any ValueError used to pass (batch-6 review M-7): the refusal must name the
        # argument R names, the first one R's message quotes.
        named = re.search(r"'([A-Za-z_.]+)'", row["error"]).group(1)
        with pytest.raises(ValueError, match=re.escape(named)):
            tx.seqdist(SX, method=method, **kwargs)
        return
    expected = np.array([float(v) for v in row["value"].split()]).reshape(4, 4, order="F")
    np.testing.assert_allclose(
        tx.seqdist(SX, method=method, **kwargs), expected, atol=1e-12, rtol=0
    )


def test_an_unknown_keyword_still_raises() -> None:
    with pytest.raises(TypeError, match="tpwo"):
        tx.seqdist(SX, method="OM", subm="CONSTANT", tpwo=2)


@pytest.mark.parametrize("key", ["tokdep.coeff", "tokdep_coeff"])
def test_omspell_token_dependent_costs_are_refused(key: str) -> None:
    # R reads opt.args only for OMspell, where tokdep.coeff selects the undocumented
    # development method OMtspell (seqdist.R:690-697; seqdist.Rd:253-255).
    with pytest.raises(ValueError, match="OMtspell"):
        tx.seqdist(SX, method="OMspell", subm="CONSTANT", opt_args={key: [1, 2, 3]})


@pytest.mark.parametrize("case", ["global_sm_number", "global_sm_number_lcs", "global_sm_vector"])
def test_a_number_or_vector_sm_is_refused_where_no_cost_uses_it(case: str) -> None:
    assert GOLDEN.loc[case, "error"].strip() == "'sm' invalid type or unknown value"
    method, kwargs = CASES[case]
    with pytest.raises(ValueError, match="'sm' invalid type or unknown value"):
        tx.seqdist(SX, method=method, **kwargs)
