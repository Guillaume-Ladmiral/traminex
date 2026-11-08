from __future__ import annotations

import inspect
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import traminex.disstree_viz as tree_viz
from traminex import (
    disstree2dot,
    disstree2dotp,
    disstreedisplay,
    seqdef,
    seqplot,
    seqrplot,
    seqsplot,
    seqtree,
    seqtree2dot,
    seqtreedisplay,
)


def _tree_fixture():
    seqdata = seqdef(
        np.array(
            [
                ["A", "A", "B"],
                ["A", "B", "B"],
                ["B", "B", "A"],
                ["B", "A", "A"],
            ]
        )
    )
    diss = np.array(
        [
            [0.0, 0.2, 0.8, 1.0],
            [0.2, 0.0, 0.7, 0.9],
            [0.8, 0.7, 0.0, 0.3],
            [1.0, 0.9, 0.3, 0.0],
        ]
    )
    predictors = pd.DataFrame({"group": ["a", "a", "b", "b"]})
    tree = seqtree(seqdata, predictors, diss=diss, R=1)
    return tree, seqdata


def _normalized_dot(path: Path, stem: Path) -> str:
    return path.read_text(encoding="utf-8").replace(str(stem), "STEM").replace(stem.name, "STEM")


def test_disstree2dot_supported_aliases_match_canonical_and_preserve_data(
    tmp_path,
) -> None:
    tree, _ = _tree_fixture()
    image_data = np.arange(4)
    original = image_data.copy()
    legend = tmp_path / "legend.png"
    quality = tmp_path / "quality.png"
    legend.write_bytes(b"legend")
    quality.write_bytes(b"quality")

    def image_fun(subset, output_path):
        Path(output_path).write_bytes(np.asarray(subset).tobytes())

    canonical_stem = tmp_path / "canonical"
    disstree2dot(
        tree,
        str(canonical_stem),
        image_fun=image_fun,
        image_data=image_data,
        only_leaf=True,
        image_format="gif",
        legend_text="Legend",
        image_legend=legend,
        image_quality=quality,
        show_depth=True,
    )

    alias_stem = tmp_path / "alias"
    with pytest.warns(DeprecationWarning) as caught:
        disstree2dot(
            tree,
            str(alias_stem),
            imagefunc=image_fun,
            imagedata=image_data,
            imgLeafOnly=True,
            imageext="gif",
            legendtext="Legend",
            legendimage=legend,
            qualityimage=quality,
            showdepth=True,
        )

    assert len(caught) == 8
    for alias in (
        "imagefunc",
        "imagedata",
        "imgLeafOnly",
        "imageext",
        "legendtext",
        "legendimage",
        "qualityimage",
        "showdepth",
    ):
        assert any(alias in str(item.message) for item in caught)
    assert _normalized_dot(alias_stem.with_suffix(".dot"), alias_stem) == _normalized_dot(
        canonical_stem.with_suffix(".dot"), canonical_stem
    )
    np.testing.assert_array_equal(image_data, original)


def test_disstree2dot_alias_default_precedence_conflict_and_stacklevel(tmp_path) -> None:
    tree, _ = _tree_fixture()
    stem = tmp_path / "precedence"
    with pytest.warns(DeprecationWarning, match=r"showdepth.*show_depth") as caught:
        disstree2dot(tree, str(stem), show_depth=False, showdepth=True)
    # Resolved: through a symbolic link the code's file name and __file__ can name
    # the same file by two paths (TX-AUD-20261010-090 (a)).
    assert Path(caught[0].filename).resolve() == Path(__file__).resolve()
    assert "rank=same" in stem.with_suffix(".dot").read_text(encoding="utf-8")

    with pytest.raises(ValueError, match=r"show_depth.*showdepth"):
        disstree2dot(
            tree,
            str(tmp_path / "conflict"),
            show_depth=True,
            showdepth=False,
        )


def test_tree_dot_wrapper_aliases_match_canonical(tmp_path) -> None:
    tree, seqdata = _tree_fixture()
    canonical_diss = tmp_path / "diss-canonical"
    disstree2dotp(tree, str(canonical_diss), with_quality=False)
    alias_diss = tmp_path / "diss-alias"
    with pytest.warns(DeprecationWarning, match=r"withquality.*with_quality"):
        disstree2dotp(tree, str(alias_diss), withquality=False)
    assert _normalized_dot(alias_diss.with_suffix(".dot"), alias_diss) == _normalized_dot(
        canonical_diss.with_suffix(".dot"), canonical_diss
    )

    canonical_seq = tmp_path / "seq-canonical"
    seqtree2dot(
        tree,
        str(canonical_seq),
        seqdata=seqdata,
        with_quality=True,
        with_legend=False,
    )
    alias_seq = tmp_path / "seq-alias"
    with pytest.warns(DeprecationWarning) as caught:
        seqtree2dot(
            tree,
            str(alias_seq),
            seqdata=seqdata,
            withquality=True,
            withlegend=False,
        )
    assert len(caught) == 2
    assert _normalized_dot(alias_seq.with_suffix(".dot"), alias_seq) == _normalized_dot(
        canonical_seq.with_suffix(".dot"), canonical_seq
    )


def test_tree_dot_wrapper_alias_default_precedence_and_conflicts(tmp_path) -> None:
    tree, seqdata = _tree_fixture()
    with pytest.warns(DeprecationWarning):
        disstree2dotp(
            tree,
            str(tmp_path / "diss-default"),
            with_quality=True,
            withquality=False,
        )
    with pytest.raises(ValueError, match=r"with_quality.*withquality"):
        disstree2dotp(
            tree,
            str(tmp_path / "diss-conflict"),
            with_quality=False,
            withquality=True,
        )

    with pytest.warns(DeprecationWarning):
        seqtree2dot(
            tree,
            str(tmp_path / "seq-default"),
            seqdata=seqdata,
            with_legend=True,
            withlegend=False,
        )
    with pytest.raises(ValueError, match=r"with_legend.*withlegend"):
        seqtree2dot(
            tree,
            str(tmp_path / "seq-conflict"),
            seqdata=seqdata,
            with_legend=False,
            withlegend=True,
        )


def test_tree_display_supported_aliases_match_canonical(tmp_path, monkeypatch) -> None:
    tree, seqdata = _tree_fixture()
    monkeypatch.setattr(tree_viz, "_render_dot", lambda *args: None)

    canonical_diss = tmp_path / "display-diss-canonical"
    expected = disstreedisplay(
        tree,
        str(canonical_diss),
        image_format="svg",
        with_quality=False,
        show_depth=True,
    )
    alias_diss = tmp_path / "display-diss-alias"
    with pytest.warns(DeprecationWarning) as caught:
        actual = disstreedisplay(
            tree,
            str(alias_diss),
            imageformat="svg",
            withquality=False,
            showdepth=True,
        )
    assert len(caught) == 3
    assert actual.endswith(".svg") and expected.endswith(".svg")
    assert _normalized_dot(alias_diss.with_suffix(".dot"), alias_diss) == _normalized_dot(
        canonical_diss.with_suffix(".dot"), canonical_diss
    )

    canonical_seq = tmp_path / "display-seq-canonical"
    expected = seqtreedisplay(
        tree,
        str(canonical_seq),
        seqdata=seqdata,
        image_format="svg",
        with_quality=False,
        with_legend=False,
        show_depth=True,
    )
    alias_seq = tmp_path / "display-seq-alias"
    with pytest.warns(DeprecationWarning) as caught:
        actual = seqtreedisplay(
            tree,
            str(alias_seq),
            seqdata=seqdata,
            imageformat="svg",
            withquality=False,
            withlegend=False,
            showdepth=True,
        )
    assert len(caught) == 4
    assert actual.endswith(".svg") and expected.endswith(".svg")
    assert _normalized_dot(alias_seq.with_suffix(".dot"), alias_seq) == _normalized_dot(
        canonical_seq.with_suffix(".dot"), canonical_seq
    )


def test_tree_display_alias_default_precedence_and_conflict(tmp_path, monkeypatch) -> None:
    tree, _ = _tree_fixture()
    monkeypatch.setattr(tree_viz, "_render_dot", lambda *args: None)
    with pytest.warns(DeprecationWarning):
        actual = disstreedisplay(
            tree,
            str(tmp_path / "precedence"),
            image_format="png",
            imageformat="svg",
        )
    assert actual.endswith(".svg")
    with pytest.raises(ValueError, match=r"image_format.*imageformat"):
        disstreedisplay(
            tree,
            str(tmp_path / "conflict"),
            image_format="jpg",
            imageformat="svg",
        )

    _, seqdata = _tree_fixture()
    with pytest.warns(DeprecationWarning):
        seqtreedisplay(
            tree,
            str(tmp_path / "seq-precedence"),
            seqdata=seqdata,
            with_legend=True,
            withlegend=False,
        )
    with pytest.raises(ValueError, match=r"with_legend.*withlegend"):
        seqtreedisplay(
            tree,
            str(tmp_path / "seq-conflict"),
            seqdata=seqdata,
            with_legend=False,
            withlegend=True,
        )


def test_seqsplot_supported_aliases_match_canonical_and_preserve_data() -> None:
    seqdata = _tree_fixture()[1]
    original = seqdata.data.copy()
    expected_fig, expected_ax = seqsplot(
        seqdata,
        main="Survival",
        cex_axis=1.25,
        with_legend=False,
        xaxis=False,
    )
    with pytest.warns(DeprecationWarning) as caught:
        actual_fig, actual_ax = seqsplot(
            seqdata,
            title="Survival",
            cex_plot=1.25,
            withlegend=False,
            axes=False,
        )

    assert len(caught) == 4
    assert actual_ax.get_title() == expected_ax.get_title()
    assert len(actual_ax.get_xticks()) == len(expected_ax.get_xticks()) == 0
    assert len(actual_fig.legends) == len(expected_fig.legends) == 0
    assert actual_ax.get_yticklabels()[0].get_fontsize() == pytest.approx(
        expected_ax.get_yticklabels()[0].get_fontsize()
    )
    np.testing.assert_array_equal(seqdata.data, original)
    plt.close(actual_fig)
    plt.close(expected_fig)


def test_seqsplot_alias_default_precedence_conflict_and_unknown() -> None:
    seqdata = _tree_fixture()[1]
    expected_fig, expected_ax = seqsplot(seqdata, main="Alias title")
    with pytest.warns(DeprecationWarning):
        fig, ax = seqsplot(seqdata, main="auto", title="Alias title")
    assert ax.get_title() == expected_ax.get_title()
    plt.close(fig)
    plt.close(expected_fig)

    with pytest.raises(ValueError, match=r"main.*title"):
        seqsplot(seqdata, main="Canonical", title="Alias")
    with pytest.raises((AttributeError, TypeError)):
        seqsplot(seqdata, not_an_option=True)


@pytest.mark.parametrize(
    "alias",
    ["title", "cex_plot", "withlegend", "axes", "tlim", "ylas"],
)
def test_seqplot_rejects_unsupported_source_aliases(alias) -> None:
    seqdata = _tree_fixture()[1]
    with pytest.raises(TypeError, match=alias):
        seqplot(seqdata, **{alias: False})


def test_seqrplot_rejects_unsupported_cex_plot_even_with_supplied_distance() -> None:
    _, seqdata = _tree_fixture()
    diss = np.zeros((seqdata.nseq, seqdata.nseq), dtype=float)
    with pytest.raises(TypeError, match="cex_plot"):
        seqrplot(seqdata, diss=diss, cex_plot=1)


@pytest.mark.parametrize(
    ("function", "names"),
    [
        (
            disstree2dot,
            {
                "imagefunc",
                "imagedata",
                "imgLeafOnly",
                "imageext",
                "legendtext",
                "legendimage",
                "qualityimage",
                "showdepth",
            },
        ),
        (disstree2dotp, {"withquality"}),
        (seqtree2dot, {"withlegend", "withquality"}),
        (disstreedisplay, {"imageformat", "withquality", "showdepth"}),
        (
            seqtreedisplay,
            {"imageformat", "withquality", "withlegend", "showdepth"},
        ),
        (seqsplot, {"title", "cex_plot", "withlegend", "axes"}),
    ],
)
def test_phase3_supported_plot_alias_runtime_signatures(function, names) -> None:
    assert names <= set(inspect.signature(function).parameters)
