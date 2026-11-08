"""Rejection messages must enumerate the set the *called layer* actually accepts.

F-011 -- ``seqdist`` accepts all four TraMineR cost-method strings
(``SOURCES/TraMineR/R/seqdist.R:151`` defines
``sm.methods <- c("TRATE","CONSTANT","INDELS","INDELSLOG")``), yet a user porting
``sm="CONSTANT"`` was told to "use subm=<float|matrix|'TRATE'>" and could
reasonably conclude ``CONSTANT`` was unsupported.  The same text is *accurate*
for the low-level ``pairwise``, whose contract really is ``'TRATE'``-only, so the
message is made layer-aware rather than deleted.

F-011 (second half) / F-020 -- ``seqdist(method="OM")`` accepts nine ``norm``
spellings, all R-exact, but an invalid one was rejected with a message naming
three.  Accepted sets measured on this fixture:

    OM                  None, 'length', 'max', 'none', 'maxlength', 'gmean',
                        'maxdist', 'YujianBo', 'auto'   (case-insensitive)
    LCS / LCP / RLCP    the same minus 'max'

F-010 -- the ``subm`` type error inside ``seqdist``'s classic-normalization
branch is unreachable (``seqdist`` resolves any cost string to an ``ndarray``
before it, and rejects an unknown string earlier).  It is kept as
defence-in-depth, so its message must name all four accepted strings rather than
``'TRATE'`` alone.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

import traminex as tx
from traminex.distances import pairwise

_RAW = np.array([["A", "B", "A", "B"], ["B", "A", "B", "A"]], dtype=object)
COST_METHODS = ("TRATE", "CONSTANT", "INDELS", "INDELSLOG")


def _sequences():
    return tx.seqdef(_RAW.copy(), alphabet=["A", "B"])


def _matrix():
    return np.asarray(_sequences().data, dtype=np.int32)


def test_seqdist_sm_alias_message_names_every_cost_method_it_accepts():
    with pytest.raises(TypeError) as excinfo:
        tx.seqdist(_sequences(), method="OM", sm="CONSTANT")
    message = str(excinfo.value)
    assert "subm" in message
    for name in COST_METHODS:
        assert name in message, (name, message)


def test_pairwise_sm_alias_message_stays_trate_only_because_that_layer_is():
    """``pairwise`` really does accept only 'TRATE'; its narrower text is correct."""
    with pytest.raises(TypeError) as excinfo:
        pairwise(_matrix(), method="OM", sm="CONSTANT")
    message = str(excinfo.value)
    assert "TRATE" in message
    assert "CONSTANT" not in message


@pytest.mark.parametrize("cost_method", COST_METHODS)
def test_seqdist_really_accepts_every_cost_method_the_message_names(cost_method):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = np.asarray(tx.seqdist(_sequences(), method="OM", subm=cost_method))
    assert result.shape == (2, 2)
    assert np.all(np.isfinite(result))


def test_seqdist_norm_message_names_every_spelling_it_accepts():
    with pytest.raises(ValueError) as excinfo:
        tx.seqdist(_sequences(), method="OM", subm="TRATE", norm="bogus")
    message = str(excinfo.value)
    for spelling in ("length", "max", "none", "maxlength", "gmean", "maxdist", "YujianBo", "auto"):
        assert spelling in message, (spelling, message)


@pytest.mark.parametrize(
    "norm",
    [None, "length", "max", "none", "maxlength", "gmean", "maxdist", "YujianBo", "auto"],
)
def test_seqdist_om_really_accepts_every_norm_the_message_names(norm):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = np.asarray(tx.seqdist(_sequences(), method="OM", subm="TRATE", norm=norm))
    assert result.shape == (2, 2)


# --- cycle-6 F-006: the ragged-sequence refusal named a set it does not mean ---
# The guard admits 16 methods -- the 14 order-based ones plus HAM and DHD -- but
# its message enumerated only the 14 and said ragged input is "only supported
# for" those. HAM and DHD are admitted deliberately: on void-bearing input where
# deletion leaves every sequence the same length they succeed, and on genuinely
# ragged input they reach the more specific "not defined for sequences of
# different length" error rather than this generic refusal.
_UNIFORM_VOID_RAW = np.array(
    [["A", "B", "A", "%"], ["B", "A", "B", "%"], ["A", "A", "B", "%"]], dtype=object
)
_RAGGED_RAW = np.array(
    [["A", "B", "A", "B"], ["B", "A", "%", "%"], ["A", "A", "B", "%"]], dtype=object
)


def _uniform_void_sequences():
    return tx.seqdef(_UNIFORM_VOID_RAW.copy(), alphabet=["A", "B"])


def _ragged_sequences():
    return tx.seqdef(_RAGGED_RAW.copy(), alphabet=["A", "B"])


def test_ragged_refusal_message_names_exactly_the_set_the_guard_admits():
    from traminex.distances import _VOID_ADMITTED_METHODS

    with pytest.raises(NotImplementedError) as excinfo:
        tx.seqdist(_ragged_sequences(), method="DUR")
    message = str(excinfo.value)
    named = {name for name in _VOID_ADMITTED_METHODS if name in message}
    assert named == set(_VOID_ADMITTED_METHODS), (
        "the refusal must enumerate every method it admits; missing "
        f"{sorted(set(_VOID_ADMITTED_METHODS) - named)}"
    )
    assert "only supported for" not in message


@pytest.mark.parametrize("method", ["DUR", "WHAM", "JSD_TRANS"])
def test_the_ragged_refusal_is_reachable_from_public_seqdist(method):
    with pytest.raises(NotImplementedError):
        tx.seqdist(_ragged_sequences(), method=method)


@pytest.mark.parametrize("method", ["HAM", "DHD"])
def test_ham_and_dhd_succeed_on_uniform_void_input(method):
    result = np.asarray(tx.seqdist(_uniform_void_sequences(), method=method))
    assert result.shape == (3, 3)
    assert np.all(np.isfinite(result))
    np.testing.assert_allclose(result, result.T)
    assert np.all(np.diag(result) == 0)


@pytest.mark.parametrize("method", ["HAM", "DHD"])
def test_ham_and_dhd_still_raise_the_specific_length_error_on_ragged_input(method):
    with pytest.raises(ValueError, match="different length"):
        tx.seqdist(_ragged_sequences(), method=method)


# --- cycle-6 F-009: the indel diagnostic named neither the argument nor its set ---
# The `indel` vocabulary genuinely differs by layer, exactly as `subm`'s does:
# `seqdist` resolves the string "auto" (and None) against the substitution matrix
# the way seqdist.R:359-365 does, while the low-level `pairwise` has no matrix to
# resolve against and takes a number. Neither layer said so: a mistyped
# indel="atuo" produced a bare "could not convert string to float: 'atuo'".


def test_seqdist_indel_message_names_the_parameter_and_its_accepted_set():
    with pytest.raises(ValueError) as excinfo:
        tx.seqdist(_sequences(), method="OM", subm=2.0, indel="atuo")
    message = str(excinfo.value)
    assert "indel" in message
    assert "auto" in message
    assert "atuo" in message


def test_pairwise_indel_message_states_its_own_narrower_set():
    with pytest.raises(ValueError) as excinfo:
        pairwise(_matrix(), method="OM", subm=2.0, indel="auto")
    message = str(excinfo.value)
    assert "indel" in message
    assert "seqdist" in message, "the narrower layer must say where 'auto' is resolved"


def test_seqdist_still_accepts_the_auto_indel_it_advertises():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        automatic = np.asarray(tx.seqdist(_sequences(), method="OM", subm=4.0, indel="auto"))
        explicit = np.asarray(tx.seqdist(_sequences(), method="OM", subm=4.0, indel=2.0))
    # seqdist.R:359-365 resolves an automatic indel as max(subm)/2.
    np.testing.assert_allclose(automatic, explicit)


def test_pairwise_still_accepts_a_numeric_indel():
    result = np.asarray(pairwise(_matrix(), method="OM", subm=2.0, indel=1.0))
    assert result.shape == (2, 2)


def test_the_indel_message_is_not_mistaken_for_the_positivity_refusal():
    """The void manifest classifies refusals by keyword; keep the two distinct.

    ``tests/test_r_goldens_distances_void.py::_normalize_exception`` maps any
    message containing "indel" plus "positive" to the ``invalid_indel`` category.
    An unknown-spelling diagnostic is a different failure and must not borrow
    that vocabulary.
    """
    with pytest.raises(ValueError) as excinfo:
        tx.seqdist(_sequences(), method="OM", subm=2.0, indel="atuo")
    lowered = str(excinfo.value).lower()
    assert "positive" not in lowered
    assert "finite" not in lowered


def test_the_defence_in_depth_subm_type_error_names_all_four_strings():
    """F-010: unreachable from ``seqdist``, so assert the message text directly."""
    import inspect

    from traminex import distances

    source = inspect.getsource(distances)
    marker = "OM subm accepts float, ndarray, or string 'TRATE'"
    # The low-level ``pairwise`` guard keeps the accurate 'TRATE'-only wording;
    # the ``seqdist`` classic-branch guard must name all four.
    assert source.count(marker) == 1, "only pairwise() may use the TRATE-only wording"
    assert "TRATE, CONSTANT, INDELS, INDELSLOG" in source


# --- cycle-6 §8.2: the same F-009 gap, for DHD's subm at the pairwise layer ---
# The cycle-6 audit refuted "the subm string-cost route is inconsistent between
# layers" -- correctly, for OM, where `pairwise` states its own accepted set
# ("OM subm accepts float, ndarray, or string 'TRATE'"). DHD had no such guard:
# a cost string fell through to np.asarray(subm, dtype=float).


def test_pairwise_dhd_subm_string_names_the_parameter_and_its_layer():
    with pytest.raises((TypeError, ValueError)) as excinfo:
        pairwise(_matrix(), method="DHD", subm="CONSTANT")
    message = str(excinfo.value)
    assert "subm" in message
    assert "seqdist" in message, "the narrower layer must say where cost strings are resolved"
    assert "could not convert string to float" not in message


def test_seqdist_dhd_subm_string_keeps_its_wider_accepted_set():
    """The `seqdist` layer resolves three cost strings and must keep saying so."""
    with pytest.raises(ValueError) as excinfo:
        tx.seqdist(_sequences(), method="DHD", subm="CONSTANT")
    message = str(excinfo.value)
    for name in ("TRATE", "INDELS", "INDELSLOG"):
        assert name in message, (name, message)


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ max\ substitution\ cost\ less\ or\ equal\ to\ \d+!:traminex.distances.python_kernels.CostDiagnosticMessage"
)
def test_dhd_still_accepts_the_arrays_and_strings_each_layer_advertises():
    sequences = _sequences()
    resolved = np.asarray(tx.seqdist(sequences, method="DHD", subm="TRATE"))
    assert resolved.shape == (2, 2)
    matrix = np.array([[0.0, 1.0], [1.0, 0.0]])
    assert np.asarray(pairwise(_matrix(), method="DHD", subm=matrix)).shape == (2, 2)


def test_an_explicit_subm_too_small_for_the_reference_is_refused_by_name():
    """A caller-supplied ``subm`` must be checked against the *combined* codes.

    ``_validate_fixed_substitution_costs`` bounds the matrix by
    ``seqmat.max()``, which covers only the source. An external ``refseq``
    contributing a higher code -- its missing state -- slipped past and reached
    the kernel, which failed with a bare
    ``IndexError: index 3 is out of bounds for axis 1 with size 3`` naming no
    parameter and no remedy (2026-09-20 T-02). The generated-cost route no
    longer produces an undersized matrix (see
    ``tests/test_seqcost_refseq_missing_contract.py``), but an explicit one
    still can.
    """
    import numpy as np

    import traminex as tx

    source = tx.seqdef(
        np.array([["A", "B", "B", "C"], ["A", "A", "C", "C"], ["B", "B", "A", "A"]], dtype=object),
        alphabet=["A", "B", "C"],
    )
    reference = tx.seqdef(np.array([["A", "B", None, "C"]], dtype=object), alphabet=["A", "B", "C"])
    too_small = np.full((3, 3), 2.0)
    np.fill_diagonal(too_small, 0.0)

    # Since 2026-09-22 (audit S1b-01) the declared alphabet fixes the width
    # before the kernel's coverage check can run: 3 states plus the missing one.
    with pytest.raises(ValueError, match=r"subm must be 4x4.*plus the missing state"):
        tx.seqdist(
            source,
            method="OM",
            subm=too_small,
            indel=1.0,
            refseq=reference,
            with_missing=True,
        )


def test_the_allowed_kwarg_list_never_advertises_a_kwarg_that_is_refused():
    """``DUR``'s "Allowed kwargs" line offered ``eps``, which is always refused.

    ``eps`` was removed from ``DUR`` as unused (``DURConfig`` rejects it too),
    and a targeted ``ValueError`` explains that. But ``eps`` was left in
    ``BASE_ALLOWED_PAIRWISE_KWARGS["DUR"]`` so the targeted message would fire
    ahead of the generic one -- with the side effect that an *unrelated* bad
    kwarg reported ``Allowed kwargs: ['eps', 'metric']``, advertising the very
    parameter the next call refuses (2026-09-20). The two messages are now
    driven from one table and cannot disagree.
    """
    import numpy as np

    import traminex as tx

    sequences = tx.seqdef(np.array([["A", "B"], ["B", "A"]], dtype=object))

    with pytest.raises(TypeError) as unknown:
        tx.seqdist(sequences, method="DUR", bogus=1)
    assert "eps" not in str(unknown.value), str(unknown.value)
    assert "metric" in str(unknown.value)

    # The targeted explanation must survive.
    with pytest.raises(ValueError, match="eps is not supported for method DUR"):
        tx.seqdist(sequences, method="DUR", eps=1e-9)
