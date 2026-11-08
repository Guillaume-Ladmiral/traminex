"""Every ``seqdist`` option R documents, checked name by name and value by value (W-05).

``seqdist`` takes its per-method options through ``**kwargs``, so
``inspect.signature`` cannot see them and the cycle-6 argument-name sweep
skipped them entirely -- roughly two dozen options across nineteen methods had
never been checked against ``SOURCES/TraMineR/man/seqdist.Rd``.

Two things are checked here:

1. **Vocabulary.** For every (method, option) pair the Rd documents, traminex
   accepts the option, spelled per ``docs/API_CONSISTENCY.md`` section 1b
   (dots become underscores, so ``add.column`` -> ``add_column`` and
   ``global.pdotj`` -> ``global_pdotj``; ``sm`` -> ``subm`` is one of the
   listed semantic renames).
2. **Value.** Each pair is exercised at a NON-DEFAULT value and compared to R,
   so an option that were accepted and then ignored would still be caught.

Generator: ``scripts/r_goldens/make_goldens_seqdist_rd_options.R``.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDEN_ATOL, R_GOLDENS_DIR, load_golden_lines

import traminex as tx
from traminex.distances._validation import BASE_ALLOWED_PAIRWISE_KWARGS

# --------------------------------------------------------------------------- #
# The Rd's per-method applicability, transcribed from the \arguments section of
# SOURCES/TraMineR/man/seqdist.Rd with the line range each claim rests on.
# R method names on the left, exactly as the Rd spells them.
# --------------------------------------------------------------------------- #
RD_OPTIONS: dict[str, tuple[str, ...]] = {
    # Rd:48-63
    "norm": (
        "OM",
        "OMloc",
        "OMslen",
        "OMspell",
        "OMstran",
        "TWED",
        "HAM",
        "DHD",
        "LCS",
        "LCP",
        "RLCP",
        "CHI2",
        "EUCLID",
    ),
    "indel": ("OM", "OMslen", "OMspell", "OMstran"),  # Rd:65-81
    "sm": ("OM", "OMloc", "OMslen", "OMspell", "OMstran", "HAM", "TWED", "DHD"),  # Rd:83-115
    "kweights": ("NMS", "NMSMST", "SVRspell"),  # Rd:131-141
    "tpow": ("OMspell", "NMSMST", "SVRspell"),  # Rd:143-148
    "expcost": ("OMloc", "OMspell"),  # Rd:150-156
    "context": ("OMloc",),  # Rd:158-162
    "link": ("OMslen",),  # Rd:164-170
    "h": ("OMslen", "TWED"),  # Rd:172-181
    "nu": ("TWED",),  # Rd:183-188
    "transindel": ("OMstran",),  # Rd:190-196
    "otto": ("OMstran",),  # Rd:198-201
    "previous": ("OMstran",),  # Rd:203-208
    "add.column": ("OMstran",),  # Rd:210-216
    "breaks": ("CHI2", "EUCLID"),  # Rd:218-222
    "step": ("CHI2", "EUCLID"),  # Rd:224-230
    "overlap": ("CHI2", "EUCLID"),  # Rd:232-238
    "global.pdotj": ("CHI2",),  # Rd:248-253
    "prox": ("NMS", "SVRspell"),  # Rd:255-260
}

#: R option name -> traminex kwarg. Everything not listed keeps its R spelling.
#: The first two follow the mechanical rule in API_CONSISTENCY.md section 1b
#: (dots become underscores); ``sm`` -> ``subm`` is a semantic rename.
RD_TO_PYTHON = {"add.column": "add_column", "global.pdotj": "global_pdotj", "sm": "subm"}

#: R method name -> the key traminex's validation table uses.
R_METHOD_TO_PYTHON = {
    "OMloc": "OMLOC",
    "OMslen": "OMSLEN",
    "OMspell": "OMSPELL",
    "OMstran": "OMSTRAN",
    "SVRspell": "SVRSPELL",
}


def _python_method(name: str) -> str:
    return R_METHOD_TO_PYTHON.get(name, name.upper())


# --------------------------------------------------------------------------- #
# Corpus and cost objects, rebuilt from the fixtures.
# --------------------------------------------------------------------------- #
ALPHABET = load_golden_lines("seqdist_rd_alphabet")
WEIGHTS = np.array([float(v) for v in load_golden_lines("seqdist_rd_weights")])
KWEIGHTS = np.array([float(v) for v in load_golden_lines("seqdist_rd_kweights")])


def _csv(name: str) -> np.ndarray:
    return pd.read_csv(R_GOLDENS_DIR / f"{name}.csv").to_numpy(dtype=float)


CODES = pd.read_csv(R_GOLDENS_DIR / "seqdist_rd_codes.csv").to_numpy(dtype=int)
SUBM = _csv("seqdist_rd_subm")
PROX = _csv("seqdist_rd_prox")
K = len(ALPHABET)
L = CODES.shape[1]
# R writes the DHD array as (K*K) x L, column t holding vec(sm[,,t]).
SUBM_ARRAY_R = _csv("seqdist_rd_subm_array").reshape((K, K, L), order="F")
# traminex puts time first: R's (K, K, L) becomes (L, K, K). Recorded divergence.
SUBM_ARRAY = np.moveaxis(SUBM_ARRAY_R, 2, 0)


def _sequences():
    states = np.asarray(ALPHABET, dtype=object)[CODES]
    return tx.seqdef(pd.DataFrame(states), alphabet=list(ALPHABET), weights=WEIGHTS)


def _expected() -> dict[str, tuple[str, np.ndarray]]:
    out: dict[str, tuple[str, np.ndarray]] = {}
    for line in (R_GOLDENS_DIR / "seqdist_rd_option_values.csv").read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        label, method, blob = line.split(",", 2)
        out[label] = (method, np.array([float(v) for v in blob.split(";")]))
    return out


EXPECTED = _expected()

#: label -> traminex call. Mirrors the generator case for case.
CASES: dict[str, dict] = {
    "norm|OM": dict(method="OM", subm=SUBM, indel=1.0, norm="maxlength"),
    "norm|OMloc": dict(method="OMLOC", subm=SUBM, expcost=0.5, norm="YujianBo"),
    "norm|OMslen": dict(method="OMSLEN", subm=SUBM, indel=1.0, norm="YujianBo"),
    "norm|OMspell": dict(method="OMSPELL", subm=SUBM, indel=1.0, norm="YujianBo"),
    "norm|OMstran": dict(method="OMSTRAN", subm=SUBM, indel=1.0, otto=0.5, norm="YujianBo"),
    "norm|TWED": dict(method="TWED", subm=SUBM, nu=0.5, h=0.5, norm="YujianBo"),
    "norm|HAM": dict(method="HAM", subm=SUBM, norm="maxlength"),
    "norm|DHD": dict(method="DHD", norm="maxlength"),
    "norm|LCS": dict(method="LCS", norm="gmean"),
    "norm|LCP": dict(method="LCP", norm="gmean"),
    "norm|RLCP": dict(method="RLCP", norm="gmean"),
    "norm|CHI2": dict(method="CHI2", step=1, norm="auto"),
    "norm|EUCLID": dict(method="EUCLID", step=1, norm="auto"),
    "indel|OM": dict(method="OM", subm=SUBM, indel=1.75),
    "indel|OMslen": dict(method="OMSLEN", subm=SUBM, indel=1.75),
    "indel|OMspell": dict(method="OMSPELL", subm=SUBM, indel=1.75),
    "indel|OMstran": dict(method="OMSTRAN", subm=SUBM, indel=1.75, otto=0.5),
    "indel_vec|OM": dict(method="OM", subm=SUBM, indel=[0.5, 1.0, 1.5, 2.0]),
    "sm|OM": dict(method="OM", subm=SUBM, indel=1.0),
    "sm|OMloc": dict(method="OMLOC", subm=SUBM, expcost=0.5),
    "sm|OMslen": dict(method="OMSLEN", subm=SUBM, indel=1.0),
    "sm|OMspell": dict(method="OMSPELL", subm=SUBM, indel=1.0),
    "sm|OMstran": dict(method="OMSTRAN", subm=SUBM, indel=1.0, otto=0.5),
    "sm|HAM": dict(method="HAM", subm=SUBM),
    "sm|TWED": dict(method="TWED", subm=SUBM, nu=0.5, h=0.5),
    "sm_array|DHD": dict(method="DHD", subm=SUBM_ARRAY),
    "sm_str_TRATE|OM": dict(method="OM", subm="TRATE", indel=1.0),
    "sm_str_CONSTANT|OM": dict(method="OM", subm="CONSTANT", indel=1.0),
    "sm_str_INDELS|OM": dict(method="OM", subm="INDELS", indel=1.0),
    "sm_str_INDELSLOG|OM": dict(method="OM", subm="INDELSLOG", indel=1.0),
    "kweights|NMS": dict(method="NMS", kweights=KWEIGHTS, prox=PROX),
    "kweights|NMSMST": dict(method="NMSMST", kweights=KWEIGHTS),
    "kweights|SVRspell": dict(method="SVRSPELL", kweights=KWEIGHTS, prox=PROX),
    "tpow|OMspell": dict(method="OMSPELL", subm=SUBM, indel=1.0, tpow=1.5),
    "tpow|NMSMST": dict(method="NMSMST", tpow=1.5),
    "tpow|SVRspell": dict(method="SVRSPELL", prox=PROX, tpow=1.5),
    "expcost|OMloc": dict(method="OMLOC", subm=SUBM, expcost=0.3),
    "expcost|OMspell": dict(method="OMSPELL", subm=SUBM, indel=1.0, expcost=0.3),
    "context|OMloc": dict(method="OMLOC", subm=SUBM, expcost=0.3, context=0.2),
    "link|OMslen": dict(method="OMSLEN", subm=SUBM, indel=1.0, link="gmean"),
    "h|OMslen": dict(method="OMSLEN", subm=SUBM, indel=1.0, h=0.75),
    "h|TWED": dict(method="TWED", subm=SUBM, nu=0.5, h=0.75),
    "nu|TWED": dict(method="TWED", subm=SUBM, nu=0.25, h=0.5),
    "transindel|OMstran": dict(
        method="OMSTRAN", subm=SUBM, indel=1.0, otto=0.5, transindel="subcost"
    ),
    "transindel_prob|OMstran": dict(
        method="OMSTRAN", subm=SUBM, indel=1.0, otto=0.5, transindel="prob"
    ),
    "otto|OMstran": dict(method="OMSTRAN", subm=SUBM, indel=1.0, otto=0.25),
    "previous|OMstran": dict(method="OMSTRAN", subm=SUBM, indel=1.0, otto=0.5, previous=True),
    "add_column|OMstran": dict(method="OMSTRAN", subm=SUBM, indel=1.0, otto=0.5, add_column=False),
    "step|CHI2": dict(method="CHI2", step=2),
    "step|EUCLID": dict(method="EUCLID", step=2),
    "overlap|CHI2": dict(method="CHI2", step=2, overlap=True),
    "overlap|EUCLID": dict(method="EUCLID", step=2, overlap=True),
    "breaks|CHI2": dict(method="CHI2", breaks=[(1, 3), (4, 8)]),
    "breaks|EUCLID": dict(method="EUCLID", breaks=[(1, 3), (4, 8)]),
    "breaks_overlapping|CHI2": dict(method="CHI2", breaks=[(1, 4), (3, 6), (5, 8)]),
    "global_pdotj_obs|CHI2": dict(method="CHI2", step=2, global_pdotj="obs"),
    "global_pdotj_vec|CHI2": dict(method="CHI2", step=2, global_pdotj=[0.4, 0.3, 0.2, 0.1]),
    "global_pdotj_unnormalised|CHI2": dict(
        method="CHI2", step=2, global_pdotj=[4.0, 3.0, 2.0, 1.0]
    ),
    "weighted_false|CHI2": dict(method="CHI2", step=2, weighted=False),
    "weighted_false|OM_TRATE": dict(method="OM", subm="TRATE", indel=1.0, weighted=False),
    "prox|NMS": dict(method="NMS", prox=PROX),
    "prox|SVRspell": dict(method="SVRSPELL", prox=PROX),
    "breaks_beats_step|CHI2": dict(method="CHI2", breaks=[(1, 3), (4, 8)], step=5, overlap=True),
}


# --------------------------------------------------------------------------- #
# 1. Vocabulary: every documented (method, option) pair is accepted.
# --------------------------------------------------------------------------- #
def test_every_documented_option_is_covered_by_a_value_case():
    """The denominator, stated rather than assumed."""
    pairs = {(method, option) for option, methods in RD_OPTIONS.items() for method in methods}
    assert len(pairs) == 51
    assert len(EXPECTED) == 63
    assert set(CASES) == set(EXPECTED)


def test_r_raised_no_error_for_any_documented_option():
    """The generator records any case R refused in its own fixture.

    That fixture was written but never read, so an R refusal would have
    surfaced only as a case missing from the values file (2026-09-22 audit
    S10-03). R refused none: the file holds its header alone.
    """
    lines = (R_GOLDENS_DIR / "seqdist_rd_option_errors.csv").read_text().splitlines()
    assert lines[0] == "# label,error"
    assert [line for line in lines[1:] if line.strip()] == []


@pytest.mark.parametrize(
    "option,r_method",
    [(option, method) for option, methods in RD_OPTIONS.items() for method in methods],
)
def test_every_documented_option_is_in_the_methods_accepted_vocabulary(option, r_method):
    """The option R documents for this method is spelled per section 1b and accepted."""
    python_option = RD_TO_PYTHON.get(option, option)
    allowed = BASE_ALLOWED_PAIRWISE_KWARGS[_python_method(r_method)]

    assert python_option in allowed, (
        f"R documents {option!r} for method {r_method}; traminex's "
        f"{_python_method(r_method)} accepts {sorted(allowed)}"
    )


def test_the_only_renames_are_the_documented_ones():
    """Guards the mapping above against absorbing a silent divergence."""
    assert RD_TO_PYTHON == {
        "add.column": "add_column",
        "global.pdotj": "global_pdotj",
        "sm": "subm",
    }


# --------------------------------------------------------------------------- #
# 2. Value: each option at a non-default value reproduces R.
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("label", sorted(EXPECTED))
def test_each_documented_option_reproduces_r_at_a_non_default_value(label):
    kwargs = dict(CASES[label])
    method = kwargs.pop("method")
    _, expected = EXPECTED[label]

    with warnings.catch_warnings():
        # A deliberately large substitution matrix trips TraMineR's own
        # "greater than twice the indel cost" notice, and breaks= trips the
        # step/overlap notice; neither is what this test is about.
        warnings.simplefilter("ignore", UserWarning)
        got = np.asarray(tx.seqdist(_sequences(), method=method, **kwargs), dtype=float)

    np.testing.assert_allclose(
        got.ravel(order="F"), expected, rtol=0, atol=R_GOLDEN_ATOL, equal_nan=True
    )


def test_the_option_grid_is_not_degenerate():
    """Distinct options must give distinct answers, or the grid proves nothing.

    Exactly three ties are expected, all properties of R's own output:

    * ``global.pdotj`` is renormalised by its sum, so 4:3:2:1 and 0.4:0.3:0.2:0.1
      are the same vector (``seqdist-CHI2.R:64``);
    * ``breaks`` overrides ``step``/``overlap`` (``seqdist.R:256-258``);
    * with ``indel=1.0``, every ``INDELS`` substitution cost on this corpus
      exceeds twice the indel (the smallest is 6.31), so no substitution is ever
      used and only indels remain; ``CONSTANT`` costs exactly 2 = 2 x indel, so
      a substitution ties with an indel pair and the alignment totals coincide.
      TraMineR emits its own "greater than twice the indel cost" notice here.
    """
    values = {label: tuple(np.round(vector, 12)) for label, (_, vector) in EXPECTED.items()}
    expected_ties = {
        frozenset({"global_pdotj_vec|CHI2", "global_pdotj_unnormalised|CHI2"}),
        frozenset({"breaks|CHI2", "breaks_beats_step|CHI2"}),
        frozenset({"sm_str_CONSTANT|OM", "sm_str_INDELS|OM"}),
    }
    seen: dict[tuple, list[str]] = {}
    for label, vector in values.items():
        seen.setdefault(vector, []).append(label)
    collisions = {frozenset(labels) for labels in seen.values() if len(labels) > 1}

    assert collisions == expected_ties
    assert len(seen) == 60, "60 distinct answers across 63 cases"


# --------------------------------------------------------------------------- #
# 3. breaks -- previously refused outright.
# --------------------------------------------------------------------------- #
def test_breaks_produces_a_different_answer_from_the_step_derived_intervals():
    """Otherwise the breaks cases could pass while breaks was ignored."""
    _, from_breaks = EXPECTED["breaks|CHI2"]
    _, from_step = EXPECTED["step|CHI2"]
    assert not np.allclose(from_breaks, from_step, rtol=0, atol=1e-9)


def test_breaks_overrides_step_and_overlap_with_a_warning_like_r():
    """seqdist.R:256-258 warns and ignores both."""
    sequences = _sequences()

    with pytest.warns(UserWarning, match="'step' has been ignored because 'breaks' is specified"):
        with_step = tx.seqdist(sequences, method="CHI2", breaks=[(1, 3), (4, 8)], step=5)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        without_step = tx.seqdist(sequences, method="CHI2", breaks=[(1, 3), (4, 8)])

    np.testing.assert_allclose(with_step, without_step, rtol=0, atol=0)


def test_overlapping_breaks_are_accepted():
    """The Rd calls the intervals "possibly overlapping" (Rd:218-222)."""
    _, expected = EXPECTED["breaks_overlapping|CHI2"]

    got = tx.seqdist(_sequences(), method="CHI2", breaks=[(1, 4), (3, 6), (5, 8)])

    np.testing.assert_allclose(
        np.asarray(got, dtype=float).ravel(order="F"), expected, rtol=0, atol=R_GOLDEN_ATOL
    )


@pytest.mark.parametrize(
    "bad,message",
    [
        ("1,3", "ordered pairs"),
        ([(3, 1)], "not ordered"),
        ([(0, 3)], "outside the sequence positions"),
        ([(1, 99)], "outside the sequence positions"),
        ([(1, 2, 3)], "exactly two positions"),
        ([], "at least one interval"),
    ],
)
def test_malformed_breaks_are_refused_with_a_message_naming_the_problem(bad, message):
    with pytest.raises(ValueError, match=message):
        tx.seqdist(_sequences(), method="CHI2", breaks=bad)


# --------------------------------------------------------------------------- #
# 4. global_pdotj -- "obs" used to crash, and a vector was not renormalised.
# --------------------------------------------------------------------------- #
def test_global_pdotj_obs_uses_the_overall_state_distribution_like_r():
    """seqdist-CHI2.R:59 resolves "obs" as seqmeant(..., prop=TRUE)."""
    sequences = _sequences()
    _, expected = EXPECTED["global_pdotj_obs|CHI2"]

    got = tx.seqdist(sequences, method="CHI2", step=2, global_pdotj="obs")

    np.testing.assert_allclose(
        np.asarray(got, dtype=float).ravel(order="F"), expected, rtol=0, atol=R_GOLDEN_ATOL
    )
    # and it is genuinely the seqmeant vector, not an accident of the default
    explicit = tx.seqdist(
        sequences,
        method="CHI2",
        step=2,
        global_pdotj=np.asarray(tx.seqmeant(sequences, prop=True), dtype=float).ravel(),
    )
    np.testing.assert_allclose(got, explicit, rtol=0, atol=0)


def test_global_pdotj_obs_differs_from_the_per_interval_default():
    _, obs = EXPECTED["global_pdotj_obs|CHI2"]
    _, default = EXPECTED["step|CHI2"]
    assert not np.allclose(obs, default, rtol=0, atol=1e-9)


def test_a_supplied_global_pdotj_vector_is_renormalised_like_r():
    """seqdist-CHI2.R:64 divides by the sum, so 4:3:2:1 equals 0.4:0.3:0.2:0.1."""
    sequences = _sequences()

    scaled = tx.seqdist(sequences, method="CHI2", step=2, global_pdotj=[4.0, 3.0, 2.0, 1.0])
    unit = tx.seqdist(sequences, method="CHI2", step=2, global_pdotj=[0.4, 0.3, 0.2, 0.1])

    # Not bit-identical: 0.4 + 0.3 + 0.2 + 0.1 sums to 1.0000000000000002, so the
    # already-normalised vector is divided by something fractionally above one
    # while 4 + 3 + 2 + 1 is exactly 10. R does the same division and its two
    # recorded matrices agree to every digit it prints.
    np.testing.assert_allclose(scaled, unit, rtol=0, atol=R_GOLDEN_ATOL)
    _, expected = EXPECTED["global_pdotj_unnormalised|CHI2"]
    np.testing.assert_allclose(
        np.asarray(scaled, dtype=float).ravel(order="F"), expected, rtol=0, atol=R_GOLDEN_ATOL
    )


@pytest.mark.parametrize(
    "bad,message",
    [
        ("observed", "either 'obs' or a vector"),
        ([0.5, 0.5], "length must equal"),
        ([-1.0, 1.0, 1.0, 1.0], "non negative"),
        ([0.0, 0.0, 0.0, 0.0], "non negative"),
    ],
)
def test_malformed_global_pdotj_is_refused_with_a_message_naming_the_option(bad, message):
    with pytest.raises(ValueError, match=message):
        tx.seqdist(_sequences(), method="CHI2", step=2, global_pdotj=bad)


# --------------------------------------------------------------------------- #
# 5. Options traminex accepts that R does not document for that method.
#    Recorded so the superset stays deliberate.
# --------------------------------------------------------------------------- #
def test_traminex_only_kwargs_are_the_expected_superset():
    """Every accepted kwarg not documented by the Rd for that method."""
    extra: dict[str, set[str]] = {}
    for python_method, allowed in BASE_ALLOWED_PAIRWISE_KWARGS.items():
        documented = {
            RD_TO_PYTHON.get(option, option)
            for option, methods in RD_OPTIONS.items()
            if python_method in {_python_method(name) for name in methods}
        }
        difference = set(allowed) - documented
        if difference:
            extra[python_method] = difference

    assert extra == {
        # traminex-only methods with no Rd entry at all
        "WHAM": {"weights"},
        # ``eps`` left this set when it stopped being an accepted kwarg: it is
        # always refused, and advertising it in the "Allowed kwargs" line offered
        # users a parameter the next call rejects (2026-09-20). The refusal and
        # its reason now live in ALWAYS_REFUSED_PAIRWISE_KWARGS.
        "DUR": {"metric"},
        "JSD_TRANS": {"eps"},
        # traminex extensions on R methods
        "OM": {"indel_by_pos", "subm_by_pos"},
        "DHD": {
            "variant",
            "smooth_window",
            "dhd_profile",
            "calibrated_params_path",
            "calibrated_params_paths",
        },
        # R documents global.pdotj for CHI2 only; traminex also accepts it for
        # EUCLID, where the metric ignores it. `weights` is the vector form of
        # R's seqdist-level `weighted`.
        "CHI2": {"weights"},
        "EUCLID": {"weights", "global_pdotj"},
        "OMSTRAN": {"weights"},
        # R's Rd comments indel out for TWED ("%, or \\code{\"TWED\"}"), but the
        # implementation still has an indel; traminex accepts it.
        "TWED": {"indel"},
    }
