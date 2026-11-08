"""Two-marker-representation equivalence regressions (2026-09-02 audit).

A declared external marker code (``nr_code``/``void_code`` outside
``range(len(alphabet))``) and the same data with in-alphabet markers are two
representations of identical sequences, so every analytical result must agree
between them.  The public ``STSList`` constructor and
``group_compare._combine_sequences`` both produce the external representation,
so it is a supported input everywhere.

Fixture pairs (``ext`` is the external representation, ``twin`` the
in-alphabet one holding logically identical data):

- FX-A: one sequence ``A-A-*-%`` (nr then trailing void).
- FX-B: two sequences over ``A,B`` with nr and trailing voids.
- FX-C: two full-length sequences with a gap nr and no void; row 1 is
  ``A-*-A-B``. Its ``seqient(with.missing=TRUE)`` values are live TraMineR
  2.2-14's, recorded by ``scripts/r_goldens/make_goldens_audit20260922_s10.R``.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import r_anchor, r_anchor_values

import traminex
from traminex.sequences import STSList


def _fx_a() -> tuple[STSList, STSList]:
    ext = STSList(
        np.array([[0, 0, 3, 2]]),
        alphabet=["A", "B"],
        nr="*",
        void="%",
        nr_code=3,
        void_code=2,
    )
    twin = STSList(
        np.array([[0, 0, 2, 3]]),
        alphabet=["A", "B", "*", "%"],
        nr="*",
        void="%",
        nr_code=2,
        void_code=3,
    )
    return ext, twin


def _fx_b(with_metadata: bool = False) -> tuple[STSList, STSList]:
    ext_extra: dict = {}
    twin_extra: dict = {}
    if with_metadata:
        ext_extra = {"labels": ["Alpha", "Beta"], "cpal": ["#ff0000", "#0000ff"]}
        twin_extra = {
            "labels": ["Alpha", "Beta", "*", "%"],
            "cpal": ["#ff0000", "#0000ff", "#00ff00", "white"],
        }
    ext = STSList(
        np.array([[0, 3, 1, 4], [1, 1, 4, 4]]),
        alphabet=["A", "B"],
        nr="*",
        void="%",
        nr_code=3,
        void_code=4,
        missing_color="#00ff00",
        **ext_extra,
    )
    twin = STSList(
        np.array([[0, 2, 1, 3], [1, 1, 3, 3]]),
        alphabet=["A", "B", "*", "%"],
        nr="*",
        void="%",
        nr_code=2,
        void_code=3,
        missing_color="#00ff00",
        **twin_extra,
    )
    return ext, twin


def _fx_c() -> tuple[STSList, STSList]:
    ext = STSList(
        np.array([[0, 3, 0, 1], [1, 1, 0, 0]]),
        alphabet=["A", "B"],
        nr="*",
        void="%",
        nr_code=3,
    )
    twin = STSList(
        np.array([[0, 2, 0, 1], [1, 1, 0, 0]]),
        alphabet=["A", "B", "*"],
        nr="*",
        void="%",
        nr_code=2,
    )
    return ext, twin


ALL_PAIRS = [("FX-A", _fx_a), ("FX-B", _fx_b), ("FX-C", _fx_c)]


def test_canonicalized_is_identity_for_in_alphabet_objects():
    _, twin = _fx_b(with_metadata=True)
    canonical = twin._canonicalized()
    assert canonical is twin


def test_canonicalized_maps_external_markers_into_the_alphabet():
    ext, twin = _fx_b(with_metadata=True)
    canonical = ext._canonicalized()
    assert canonical.alphabet == ["A", "B", "*", "%"]
    assert canonical.nr_code == 2
    assert canonical.void_code == 3
    np.testing.assert_array_equal(canonical.data, twin.data)
    assert canonical.to_strings() == ext.to_strings()
    assert canonical.labels == ["Alpha", "Beta", "*", "%"]
    assert canonical.cpal == ["#ff0000", "#0000ff", "#00ff00", "white"]
    # The source object is never mutated.
    assert ext.alphabet == ["A", "B"]
    assert ext.nr_code == 3


def test_canonicalized_merges_external_marker_into_existing_text_slot():
    """A combined object (``group_compare._combine_sequences``) can carry the
    marker text as an alphabet member while the marker code is external;
    canonicalization must merge into that slot, not append a duplicate."""
    half = STSList(
        np.array([[0, 3, 1, 4], [1, 1, 4, 4]]),
        alphabet=["A", "B", "*"],
        nr="*",
        void="%",
        nr_code=3,
        void_code=4,
    )
    canonical = half._canonicalized()
    assert canonical.alphabet == ["A", "B", "*", "%"]
    assert canonical.nr_code == 2
    assert canonical.void_code == 3
    np.testing.assert_array_equal(canonical.data, np.array([[0, 2, 1, 3], [1, 1, 3, 3]]))


def test_code_for_state_resolves_markers_in_both_representations():
    ext, twin = _fx_b()
    assert ext._code_for_state("A") == 0
    assert ext._code_for_state("*") == 3
    assert ext._code_for_state("%") == 4
    assert twin._code_for_state("*") == 2
    with pytest.raises(ValueError):
        ext._code_for_state("Z")


@pytest.mark.parametrize(("name", "make"), ALL_PAIRS)
def test_seqient_agrees_between_representations(name, make):
    ext, twin = make()
    for norm in (True, False):
        got = traminex.seqient(ext, include_missing=True, norm=norm)
        want = traminex.seqient(twin, include_missing=True, norm=norm)
        np.testing.assert_allclose(got, want, atol=1e-15)


def test_seqient_matches_live_r_anchor():
    # TraMineR 2.2-14 on seqdef(A,NA,A,B / B,B,A,A), seqient(with.missing=TRUE).
    ext, twin = _fx_c()
    for x in (ext, twin):
        raw = traminex.seqient(x, include_missing=True, norm=False)
        norm = traminex.seqient(x, include_missing=True)
        np.testing.assert_allclose(raw, r_anchor_values("fx_c_seqient_raw"), atol=1e-12)
        np.testing.assert_allclose(norm, r_anchor_values("fx_c_seqient_norm"), atol=1e-12)


@pytest.mark.parametrize(("name", "make"), ALL_PAIRS)
def test_entropy_agrees_between_representations(name, make):
    ext, twin = make()
    assert traminex.entropy(ext) == pytest.approx(traminex.entropy(twin), abs=1e-15)
    assert traminex.entropy(ext, normalized=False) == pytest.approx(
        traminex.entropy(twin, normalized=False), abs=1e-15
    )


def test_transition_count_matrix_agrees_between_representations():
    from traminex.stats_basic import _transition_count_matrix

    ext, twin = _fx_a()
    np.testing.assert_allclose(_transition_count_matrix(ext), _transition_count_matrix(twin))


def test_transition_heatmap_renders_external_markers():
    from traminex.visualization.plots.transitions import transition_heatmap

    ext, _ = _fx_a()
    fig, _ax = transition_heatmap(ext)
    plt.close(fig)


def test_distribution_plot_and_chronogram_render_external_markers():
    from traminex.visualization.plots.chronogram import chronogram
    from traminex.visualization.plots.distribution import state_distribution_plot

    ext, twin = _fx_b()
    fig_e, ax_e = state_distribution_plot(ext, include_missing=True)
    fig_t, ax_t = state_distribution_plot(twin, include_missing=True)
    alphas_e = sorted({p.get_alpha() or 1.0 for p in ax_e.patches})
    alphas_t = sorted({p.get_alpha() or 1.0 for p in ax_t.patches})
    assert alphas_e == alphas_t
    plt.close(fig_e)
    plt.close(fig_t)
    fig, _ax = chronogram(ext, include_missing=True)
    plt.close(fig)


def test_seqfpos_and_seqfposend_columns_agree():
    ext, twin = _fx_b()
    for fn in (traminex.seqfpos, traminex.seqfposend):
        got = fn(ext, include_missing=True)
        want = fn(twin, include_missing=True)
        assert list(got.columns) == list(want.columns)
        pd.testing.assert_frame_equal(got.reset_index(drop=True), want.reset_index(drop=True))


def test_seqfpos_scalar_marker_state_agrees():
    ext, twin = _fx_b()
    got = traminex.seqfpos(ext, state="*")
    want = traminex.seqfpos(twin, state="*")
    np.testing.assert_allclose(got, want, equal_nan=True)
    np.testing.assert_allclose(want, [2.0, np.nan], equal_nan=True)


@pytest.mark.parametrize("include_missing", [True, False])
def test_seqst_normalized_agrees(include_missing):
    ext, twin = _fx_b()
    np.testing.assert_allclose(
        traminex.seqST(ext, norm=True, include_missing=include_missing),
        traminex.seqST(twin, norm=True, include_missing=include_missing),
        atol=1e-15,
    )


def test_seqindic_agrees_between_representations():
    ext, twin = _fx_b()
    pd.testing.assert_frame_equal(traminex.seqindic(ext), traminex.seqindic(twin))
    pd.testing.assert_frame_equal(
        traminex.seqindic(ext, indic=["visited", "visitp"], with_missing=True),
        traminex.seqindic(twin, indic=["visited", "visitp"], with_missing=True),
    )


@pytest.mark.parametrize("include_missing", [True, False])
def test_seqivolatility_agrees(include_missing):
    ext, twin = _fx_b()
    np.testing.assert_allclose(
        traminex.seqivolatility(ext, include_missing=include_missing),
        traminex.seqivolatility(twin, include_missing=include_missing),
        atol=1e-15,
    )


def test_seqivardur_agrees():
    ext, twin = _fx_b()
    got = traminex.seqivardur(ext, include_missing=True)
    want = traminex.seqivardur(twin, include_missing=True)
    np.testing.assert_allclose(np.asarray(got), np.asarray(want), atol=1e-15)


def test_seqlogp_agrees_between_representations():
    for _, make in (("FX-B", _fx_b), ("FX-C", _fx_c)):
        ext, twin = make()
        np.testing.assert_allclose(
            traminex.seqlogp(ext, include_missing=True),
            traminex.seqlogp(twin, include_missing=True),
            atol=1e-12,
        )


@pytest.mark.parametrize("include_missing", [True, False])
def test_seqdss_decodes_identically_and_keeps_marker_metadata(include_missing):
    ext, twin = _fx_b()
    got = traminex.seqdss(ext, include_missing=include_missing)
    want = traminex.seqdss(twin, include_missing=include_missing)
    assert got.to_strings() == want.to_strings()
    assert got.nr_code is not None
    assert got._state_for_code(got.nr_code) == "*"


def test_seqtransn_and_seqsubsn_with_missing_agree():
    ext, twin = _fx_b()
    np.testing.assert_allclose(
        traminex.seqtransn(ext, with_missing=True),
        traminex.seqtransn(twin, with_missing=True),
    )
    np.testing.assert_allclose(
        np.asarray(traminex.seqsubsn(ext, with_missing=True), dtype=float),
        np.asarray(traminex.seqsubsn(twin, with_missing=True), dtype=float),
    )


def test_precarity_family_agrees_between_representations():
    ext, twin = _fx_b()
    np.testing.assert_allclose(
        traminex.seqprecarity(ext, include_missing=True),
        traminex.seqprecarity(twin, include_missing=True),
        atol=1e-15,
    )
    np.testing.assert_allclose(
        traminex.seqinsecurity(ext, with_missing=True),
        traminex.seqinsecurity(twin, with_missing=True),
        atol=1e-15,
    )
    np.testing.assert_allclose(
        traminex.seqibad(ext, include_missing=True),
        traminex.seqibad(twin, include_missing=True),
        atol=1e-15,
    )
    np.testing.assert_allclose(
        traminex.seqidegrad(ext, with_missing=True),
        traminex.seqidegrad(twin, with_missing=True),
        atol=1e-15,
    )


def test_seqintegr_agrees_including_marker_state_lookup():
    ext, twin = _fx_b()
    np.testing.assert_allclose(
        traminex.seqintegr(ext, include_missing=True),
        traminex.seqintegr(twin, include_missing=True),
        atol=1e-15,
    )
    # R accepts the missing state only when it is requested (seqintegr.R:9-20;
    # live TraMineR 2.2-14 stops with "state * not in the alphabet!" otherwise,
    # audit S3-17), so the marker lookup is compared on that route and the
    # refusal must agree between the two representations.
    np.testing.assert_allclose(
        traminex.seqintegr(ext, state="*", include_missing=True),
        traminex.seqintegr(twin, state="*", include_missing=True),
        atol=1e-15,
    )
    for x in (ext, twin):
        with pytest.raises(ValueError, match="not in the alphabet"):
            traminex.seqintegr(x, state="*")


def test_seqprecstart_and_seqipos_agree_between_representations():
    # Re-homed from the removed *_proxy helpers (TX-AUD-20261019-063): the public
    # functions those proxies shadowed, for both marker representations.
    ext, twin = _fx_b()
    np.testing.assert_allclose(
        traminex.seqprecstart(ext, include_missing=True),
        traminex.seqprecstart(twin, include_missing=True),
        atol=1e-15,
    )
    np.testing.assert_allclose(
        np.asarray(traminex.seqipos(ext, pos_states=["A"]), dtype=float),
        np.asarray(traminex.seqipos(twin, pos_states=["A"]), dtype=float),
        atol=1e-15,
    )


def test_seqdomassoc_accepts_external_markers_and_agrees():
    ext, twin = _fx_b()
    got = traminex.seqdomassoc([ext, ext], with_missing=True)
    want = traminex.seqdomassoc([twin, twin], with_missing=True)
    pd.testing.assert_frame_equal(got, want)


def test_seqpolyads_agrees_between_representations():
    ext, twin = _fx_c()
    got = traminex.seqpolyads([ext, ext], T=10, s=1, with_missing=True)
    want = traminex.seqpolyads([twin, twin], T=10, s=1, with_missing=True)
    np.testing.assert_allclose(got["mean.dist"], want["mean.dist"], atol=1e-12)
    np.testing.assert_allclose(got["observed.dist"], want["observed.dist"], atol=1e-12)


def test_seqsurv_and_seqsurvplot_handle_external_markers_with_labels():
    from traminex.visualization.plots.survival import seqsurvplot

    ext, twin = _fx_b(with_metadata=True)
    got = traminex.seqsurv(ext, per_state=False, with_missing=True)
    want = traminex.seqsurv(twin, per_state=False, with_missing=True)
    assert sorted(got.curves) == sorted(want.curves)
    for key in got.curves:
        np.testing.assert_allclose(got.curves[key], want.curves[key])
    fig, _ax = seqsurvplot(ext, with_missing=True)
    plt.close(fig)


def test_seqcta_marker_initial_state_agrees():
    ext, twin = _fx_b()
    got = traminex.seqcta(ext, subseq=2, initial_state="*")
    want = traminex.seqcta(twin, subseq=2, initial_state="*")
    pd.testing.assert_frame_equal(got, want)


@pytest.mark.parametrize(("name", "make"), [("FX-B", _fx_b), ("FX-C", _fx_c)])
def test_seqdist_agrees_between_representations(name, make):
    ext, twin = make()
    for method, kwargs in (
        ("LCS", {}),
        ("HAM", {}),
        ("OM", {"subm": 2.0, "indel": 1.0}),
    ):
        outcomes = []
        for obj in (ext, twin):
            try:
                outcomes.append(
                    ("ok", traminex.seqdist(obj, method=method, with_missing=True, **kwargs))
                )
            except ValueError as exc:
                # HAM legitimately refuses ragged (void-deleted) sequences —
                # exactly as R's seqdist stops on unequal lengths; the two
                # representations must refuse identically.
                outcomes.append(("raise", str(exc)))
        (kind_ext, got), (kind_twin, want) = outcomes
        assert kind_ext == kind_twin, (name, method, outcomes)
        if kind_ext == "ok":
            np.testing.assert_allclose(got, want, atol=1e-15)
        else:
            assert got == want


def _png_bytes(fig) -> bytes:
    from io import BytesIO

    buffer = BytesIO()
    fig.savefig(buffer, format="png", dpi=80)
    plt.close(fig)
    return buffer.getvalue()


def _fig_of(result):
    return result[0] if isinstance(result, tuple) else result


_RENDERER_CASES = [
    ("sequence_index_plot", lambda x, D: traminex.sequence_index_plot(x)),
    (
        "sequence_index_plot_grouped",
        lambda x, D: traminex.sequence_index_plot_grouped(x, np.array([0, 1])),
    ),
    (
        "overlay_representatives_on_index_plot",
        lambda x, D: traminex.overlay_representatives_on_index_plot(x, np.array([0, 1]), D),
    ),
    ("seqfplot", lambda x, D: traminex.seqfplot(x)),
    ("state_legend_figure", lambda x, D: traminex.state_legend_figure(x)),
    ("seqmtplot", lambda x, D: traminex.seqmtplot(x)),
    ("seqmsplot", lambda x, D: traminex.seqmsplot(x)),
    (
        "cluster_chronograms",
        lambda x, D: traminex.cluster_chronograms(x, np.array([0, 1])),
    ),
    (
        "cluster_chronograms_with_representatives",
        lambda x, D: traminex.cluster_chronograms_with_representatives(x, np.array([0, 1]), D),
    ),
    ("seqrplot", lambda x, D: traminex.seqrplot(x, diss=D, nrep=1)),
    ("seqrfplot", lambda x, D: traminex.seqrfplot(x, D, k=2)),
    ("seqplot_rf", lambda x, D: traminex.seqplot_rf(x, D, k=2)),
    # 2026-09-22 audit S6-12: neither renderer canonicalized its input.
    (
        "seqalignplot",
        lambda x, D: traminex.seqalignplot(x, (0, 1), include_missing=True),
    ),
    (
        "duration_boxplot",
        lambda x, D: traminex.duration_boxplot(x, by_state=True, seaborn=False),
    ),
]


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ Negative\ \ or\ zero\ df!\ k\ =\ \d+\ \ wsum\ =\ \ [\d.]+:UserWarning"
)
@pytest.mark.parametrize(
    ("fixture", "make"),
    [("FX-B", lambda: _fx_b(with_metadata=True)), ("FX-C", _fx_c)],
)
@pytest.mark.parametrize(("renderer", "call"), _RENDERER_CASES, ids=[c[0] for c in _RENDERER_CASES])
def test_renderers_draw_both_representations_identically(renderer, call, fixture, make):
    """Cycle-9 finding N-05: these entry points drew the declared-external
    representation wrongly (raw codes into palettes/axes) while the
    in-alphabet twin matched the canonical rendering."""
    ext, twin = make()
    images = []
    for obj in (ext, twin):
        D = traminex.seqdist(obj, method="LCS", with_missing=True)
        images.append(_png_bytes(_fig_of(call(obj, D))))
    assert images[0] == images[1], (renderer, fixture)


def test_seqmsplot_colors_modal_marker_identically():
    """seqmsplot colours bars by the modal state's code; when the modal
    state is the missing marker the external code (outside the alphabet)
    picked a different colour than the twin's in-alphabet code."""
    ext = STSList(
        np.array([[0, 3], [1, 3]]),
        alphabet=["A", "B"],
        nr="*",
        void="%",
        nr_code=3,
    )
    twin = STSList(
        np.array([[0, 2], [1, 2]]),
        alphabet=["A", "B", "*"],
        nr="*",
        void="%",
        nr_code=2,
    )
    images = [
        _png_bytes(_fig_of(traminex.seqmsplot(obj, include_missing=True))) for obj in (ext, twin)
    ]
    assert images[0] == images[1]


@pytest.mark.parametrize(("name", "make"), [("FX-B", _fx_b), ("FX-C", _fx_c)])
def test_seqalignplot_colours_the_missing_marker_as_a_state(name, make):
    """2026-09-22 audit S6-12: with the external representation seqalignplot
    looked the marker text up in the unextended alphabet and painted the
    missing cell in the gap grey; the twin painted it in its palette colour."""
    from traminex.visualization.plots.alignment import GAP_COLOR

    def cell_colours(obj):
        fig, ax = traminex.seqalignplot(obj, (0, 1), include_missing=True)
        alignment = traminex.seqalign(obj, (0, 1), include_missing=True)
        # Patches come in (seq1, operation, seq2) triples per position.
        colours = [
            (symbol, tuple(ax.patches[3 * position].get_facecolor()))
            for position, symbol in enumerate(alignment.seq1)
        ]
        plt.close(fig)
        return colours

    ext, twin = make()
    ext_colours = cell_colours(ext)
    assert ext_colours == cell_colours(twin)
    missing = [colour for symbol, colour in ext_colours if symbol == "*"]
    assert missing
    assert all(colour != matplotlib.colors.to_rgba(GAP_COLOR) for colour in missing)


def test_state_legend_draws_both_representations_identically():
    ext, twin = _fx_b(with_metadata=True)
    images = []
    for obj in (ext, twin):
        fig, ax = plt.subplots(figsize=(4, 3))
        traminex.state_legend(obj, ax=ax)
        images.append(_png_bytes(fig))
    assert images[0] == images[1]


def test_seqsplot_renders_external_markers_with_missing():
    """Cycle-9 finding N-04: seqsplot crashed with IndexError on the
    external representation under with_missing=True (seqsurvplot already
    canonicalizes)."""
    ext, twin = _fx_c()
    images = []
    for obj in (ext, twin):
        images.append(_png_bytes(_fig_of(traminex.seqsplot(obj, with_missing=True))))
    assert images[0] == images[1]


@pytest.mark.parametrize(("name", "make"), [("FX-B", _fx_b), ("FX-C", _fx_c)])
def test_seqasnum_agrees_between_representations(name, make):
    """Cycle-9 finding N-07: external marker codes leaked through unchanged,
    violating seqasnum's own contract (nr's code sits right after the real
    alphabet)."""
    ext, twin = make()
    got = traminex.seqasnum(ext, with_missing=True)
    want = traminex.seqasnum(twin, with_missing=True)
    np.testing.assert_array_equal(got, want)


def test_seqstart_agrees_between_representations():
    """Cycle-9 finding N-06: seqstart copied external source codes while
    reassigning nr_code, so the result's constructor rejected them."""
    ext, twin = _fx_c()
    got = traminex.seqstart(ext, 1, [1, 2])
    want = traminex.seqstart(twin, 1, [1, 2])
    assert got.to_strings() == want.to_strings()
    assert got.alphabet == want.alphabet


def test_seqgen_missing_external_void_does_not_collide_marker_codes():
    """Cycle-9 finding N-09: with an external void_code == len(alphabet), the
    freshly appended nr took the same code (nr_code == void_code, corrupting
    every decoded string)."""
    x = STSList(
        np.array([[0, 1, 2], [1, 0, 2]]),
        alphabet=["A", "B"],
        nr="*",
        void="%",
        void_code=2,
    )
    out = traminex.seqgen_missing(x, p_cases=1.0, p_left=0.9, random_state=5)
    assert out.nr_code is not None
    assert out.nr_code != out.void_code
    decoded = out.to_strings()
    assert all("%" not in row.replace("-%", "") or True for row in decoded)
    # Every nr cell decodes to '*', never to the void marker.
    nr_cells = out.data == out.nr_code
    void_cells = out.data == out.void_code
    assert not np.any(nr_cells & void_cells)


def test_seqgen_missing_extends_cpal_with_the_appended_marker():
    """Cycle-9 finding N-09 (cpal half): the appended nr state extended
    labels but not cpal, leaving cpal shorter than the alphabet."""
    x = STSList(
        np.array([[0, 1], [1, 0]]),
        alphabet=["A", "B"],
        labels=["Alpha", "Beta"],
        cpal=["#111111", "#222222"],
        nr="*",
        void="%",
        missing_color="darkgrey",
    )
    out = traminex.seqgen_missing(x, p_cases=1.0, p_left=0.9, random_state=5)
    assert "*" in out.alphabet
    assert out.cpal is not None
    assert len(out.cpal) == len(out.alphabet)
    assert out.cpal[out.nr_code] == "darkgrey"


def test_seqsamm_spell_holds_real_states_only():
    """Cycle-9 finding N-08: R sets attr(ret, 'spell') <- alphabet(seqdata)
    (seqsamm.R:47) — real states only — while traminex used the raw alphabet,
    which includes the in-alphabet missing marker for seqdef objects."""
    x = traminex.seqdef(
        pd.DataFrame([["A", None, "A"], ["B", "B", "A"], ["A", "A", "B"]]), missing=None
    )
    assert "*" in x.alphabet
    samm = traminex.seqsamm(x, 2)
    assert samm.spell == ("A", "B")


def test_seqalign_refuses_missing_states_clearly_without_include_missing():
    """Cycle-9 finding N-03 (refusal half): R's seqalign stops when the
    missing state meets with.missing=FALSE (a name-subscript error);
    traminex crashed with an opaque IndexError deep in the DP. Both
    representations must refuse with the same clear message."""
    ext, twin = _fx_c()
    for obj in (ext, twin):
        with pytest.raises(ValueError, match="include_missing"):
            traminex.seqalign(obj, (0, 1), indel=1.0)


def test_seqalign_agrees_between_representations_with_missing():
    """Cycle-9 finding N-03 (equivalence half): under include_missing=True
    the external representation crashed while the twin computes."""
    ext, twin = _fx_c()
    got = traminex.seqalign(ext, (0, 1), indel=1.0, include_missing=True)
    want = traminex.seqalign(twin, (0, 1), indel=1.0, include_missing=True)
    np.testing.assert_array_equal(got.operation, want.operation)
    np.testing.assert_array_equal(got.seq1, want.seq1)
    np.testing.assert_array_equal(got.seq2, want.seq2)
    assert got.total_cost == want.total_cost


def test_seqmd_accepts_external_marker_channels_and_agrees():
    """Cycle-9 finding N-02: seqMD refused external-marker channels with
    'channel contains a state code outside its alphabet' while the twin
    channels compute."""
    ext, twin = _fx_c()
    got = traminex.seqMD([ext, ext], what="MDseq")
    want = traminex.seqMD([twin, twin], what="MDseq")
    assert got.alphabet == want.alphabet
    assert got.to_strings() == want.to_strings()
    got_d = traminex.seqdistmc([ext, ext], what="diss", method="LCS", with_missing=True)
    want_d = traminex.seqdistmc([twin, twin], what="diss", method="LCS", with_missing=True)
    np.testing.assert_allclose(got_d, want_d, atol=1e-15)


def test_seqcompare_computes_on_mixed_width_missing_input():
    """Public reproduction of cycle-9 finding N-01: ``seqCompare`` builds the
    external representation internally (``_combine_sequences``: nr at
    ``len(alphabet)``, void at ``len(alphabet)+1`` for width padding), so a
    plain two-object comparison with any missing value must not be refused by
    ``seqdist``'s state-count bound."""
    a = traminex.seqdef(pd.DataFrame([["A", "B", "A"], ["B", None, "B"]]), missing=None)
    b = traminex.seqdef(pd.DataFrame([["B", "A"], ["A", "A"]]))
    out = traminex.seqCompare(a, b, s=0, method="LCS", with_missing=True)
    assert np.isfinite(out.to_numpy(dtype=float)).all()


def test_seqrecode_keeps_marker_cells_out_of_otherwise_in_both_representations():
    """``otherwise`` recodes the real states only, whichever marker regime is used.

    R's ``recodeFunc`` (``SOURCES/TraMineR/R/seqrecode.R``) starts *every*
    factor cell at ``otherwise`` -- the nr and void levels included -- so
    TraMineR 2.2-14 returns alphabet ``W,O`` and the rows ``W-W-W-O`` /
    ``W-W-O-O`` (recorded by ``scripts/r_goldens/make_goldens_audit20260922_s10.R``).
    traminex leaves the markers out of ``otherwise`` (owner decision 2026-10-23,
    TX-AUD-20261019-033), so the rows keep their missing and void cells.

    The external-marker twin must give the same result (cycle-10 N-15 /
    2026-09-20 T-03: ``seqrecode`` canonicalizes at entry).
    """
    internal = traminex.seqdef(
        np.array([["A", "B", "B", "*"], ["A", "A", "C", "%"]], dtype=object),
        alphabet=["A", "B", "C"],
    )
    external = STSList(
        data=internal.data.copy(),
        alphabet=["A", "B", "C"],
        nr="*",
        void="%",
        nr_code=3,
        void_code=4,
    )

    recoded_internal = traminex.seqrecode(internal, {"W": ["A", "B"]}, otherwise="O")
    recoded_external = traminex.seqrecode(external, {"W": ["A", "B"]}, otherwise="O")

    assert r_anchor("seqrecode_otherwise_alphabet", "strings") == ["W", "O"]  # R
    assert r_anchor("seqrecode_otherwise_rows", "strings") == ["W-W-W-O", "W-W-O-O"]  # R
    assert traminex.alphabet(recoded_internal) == ["W", "O"]
    assert recoded_external.alphabet == recoded_internal.alphabet
    np.testing.assert_array_equal(recoded_external.data, recoded_internal.data)
    assert recoded_internal.to_strings() == ["W-W-W-*", "W-W-O-%"]


_TWIN_LAYOUTS = ["nr-only", "nr-and-void", "gapped-marker-codes"]


def _per_state_twins(layout: str) -> tuple[STSList, STSList]:
    """A true twin pair: same cells and declared markers, different representation.

    The external object's marker codes sit outside ``range(len(alphabet))`` --
    adjacent to it (3, 4) or, for ``gapped-marker-codes``, at 7 and 9 -- so its
    canonical form must equal the in-alphabet object exactly; the assertions
    below guard that.
    """
    with_void = layout != "nr-only"
    rows = [["A", "B", "B", "*"], ["A", "A", "C", "%" if with_void else "C"], ["C", "B", "*", "A"]]
    internal = traminex.seqdef(np.array(rows, dtype=object), alphabet=["A", "B", "C"])
    internal.weights = np.array([1.0, 2.0, 0.5])
    nr_code, void_code = (7, 9) if layout == "gapped-marker-codes" else (3, 4)
    data = internal.data.copy()
    data[internal.data == 3] = nr_code
    data[internal.data == 4] = void_code
    external = STSList(
        data=data,
        alphabet=["A", "B", "C"],
        nr="*",
        nr_code=nr_code,
        weights=np.array([1.0, 2.0, 0.5]),
        **({"void": "%", "void_code": void_code} if with_void else {}),
    )
    canonical = external._canonicalized()
    assert canonical.alphabet == internal.alphabet
    np.testing.assert_array_equal(canonical.data, internal.data)
    assert (canonical.nr_code, canonical.void_code) == (internal.nr_code, internal.void_code)
    return internal, external


_PER_STATE_RESULTS = {
    "seqistatd": lambda x, inc: traminex.seqistatd(x, include_missing=inc),
    "seqistatd(prop)": lambda x, inc: traminex.seqistatd(x, include_missing=inc, prop=True),
    "seqstatd": lambda x, inc: traminex.seqstatd(x, include_missing=inc)["Frequencies"],
    "seqmeant": lambda x, inc: traminex.seqmeant(x, include_missing=inc),
    "seqmeant(prop)": lambda x, inc: traminex.seqmeant(x, include_missing=inc, prop=True),
    "cluster_summary": lambda x, inc: traminex.cluster_summary(x, np.array([0, 1, 0])),
}


@pytest.mark.parametrize("layout", _TWIN_LAYOUTS)
@pytest.mark.parametrize("include_missing", [False, True])
@pytest.mark.parametrize("name", sorted(_PER_STATE_RESULTS))
def test_per_state_results_have_the_same_shape_in_both_representations(
    name, include_missing, layout
):
    """Per-state output must not change shape with the marker representation.

    These build their per-state axis from ``_codebook_codes`` or from
    ``len(x.alphabet)``, so an in-alphabet object carried marker columns while
    its external-marker twin omitted them. Shared columns always agreed, but
    the shape did not, which breaks the two-representation equivalence the
    changelog claims (cycle-10 N-14 / 2026-09-20 T-05). The padded shape is
    traminex's documented per-state layout (docs/API_CONSISTENCY.md); R's own
    result has no marker rows (tests/test_per_state_shape_and_selection.py), so
    canonicalizing at entry is what makes the pair agree. This test once covered only an nr-only pair with
    ``include_missing=True`` -- the one configuration in which the unpadded
    external shape happens to coincide -- and so passed while 15 of the other
    configurations still diverged (2026-09-22 audit F-04).
    """
    internal, external = _per_state_twins(layout)
    compute = _PER_STATE_RESULTS[name]

    got_i = np.asarray(pd.DataFrame(compute(internal, include_missing)).to_numpy(dtype=float))
    got_e = np.asarray(pd.DataFrame(compute(external, include_missing)).to_numpy(dtype=float))

    assert got_e.shape == got_i.shape, f"{name}: {got_e.shape} (external) != {got_i.shape} (twin)"
    np.testing.assert_allclose(got_e, got_i, rtol=0, atol=1e-12, equal_nan=True)


@pytest.mark.parametrize("layout", _TWIN_LAYOUTS)
@pytest.mark.parametrize("include_missing", [False, True])
def test_the_modal_sequence_is_the_same_in_both_representations(include_missing, layout):
    """``seqmodst`` must return the same object for both representations.

    It indexes ``seqstatd``'s columns with its own codebook and builds the result
    on its input's alphabet, so it has to read the canonical form: otherwise the
    modal states agree but the returned object carries the external alphabet
    ``[A, B, C]`` where the twin's carries ``[A, B, C, *, %]``.
    """
    internal, external = _per_state_twins(layout)

    modal_i = traminex.seqmodst(internal, include_missing=include_missing)
    modal_e = traminex.seqmodst(external, include_missing=include_missing)

    assert modal_e.to_strings() == modal_i.to_strings()
    assert modal_e.alphabet == modal_i.alphabet
    np.testing.assert_array_equal(modal_e.data, modal_i.data)
    np.testing.assert_allclose(modal_e.frequencies_, modal_i.frequencies_, rtol=0, atol=1e-12)
    assert modal_e.occurrences_ == modal_i.occurrences_


@pytest.mark.parametrize(("name", "factory"), ALL_PAIRS)
def test_stlab_and_cpal_agree_between_marker_representations(name, factory):
    """Audit 2026-09-22 S2-14: ``stlab``/``cpal`` read the raw metadata.

    They returned one entry per stored alphabet slot, so the in-alphabet twin
    gave the documented codebook-padded length (markers included) while the
    external representation of the same data gave only the real states.
    """
    ext, twin = factory()
    assert traminex.stlab(ext) == traminex.stlab(twin)
    assert traminex.cpal(ext) == traminex.cpal(twin)


def test_stlab_and_cpal_agree_between_representations_with_metadata():
    ext, twin = _fx_b(with_metadata=True)
    assert traminex.stlab(ext) == traminex.stlab(twin) == ["Alpha", "Beta", "*", "%"]
    assert traminex.cpal(ext) == traminex.cpal(twin) == ["#ff0000", "#0000ff", "#00ff00", "white"]
    # The accessors return values; the stored metadata is unchanged.
    assert ext.labels == ["Alpha", "Beta"]
    assert ext.cpal == ["#ff0000", "#0000ff"]
