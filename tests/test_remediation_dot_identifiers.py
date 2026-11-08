"""TX-AUD-20261005-009: Graphviz sees one node per tree/depth identifier."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from traminex.disstree import disstree, disstree2dot


@pytest.fixture(scope="module")
def graphviz_dot():
    configured = os.environ.get("GVPATH")
    executable = Path(configured) if configured else Path(shutil.which("dot") or "")
    if executable.is_dir():
        executable = executable / "bin" / "dot"
    if not executable.is_file() or not os.access(executable, os.X_OK):
        pytest.skip("Graphviz dot is unavailable")
    return str(executable.resolve())


@pytest.fixture(scope="module")
def three_node_tree():
    distances = np.array([[0, 1, 4, 4], [1, 0, 4, 4], [4, 4, 0, 1], [4, 4, 1, 0]], float)
    return disstree(
        distances,
        pd.DataFrame({"group": ["a", "a", "b", "b"]}),
        R=0,
        min_size=1,
        max_depth=2,
    )


@pytest.mark.parametrize(
    "basename",
    ["normal", "with space", 'with-quote"', "with-backslash\\", "trajectory_日本語"],
    ids=["ordinary", "space", "quote", "backslash", "unicode"],
)
@pytest.mark.parametrize("show_depth", [False, True], ids=["no-depth", "depth"])
def test_rendered_dot_has_exact_tree_node_identities(
    tmp_path, graphviz_dot, three_node_tree, basename, show_depth
):
    stem = tmp_path / basename
    disstree2dot(three_node_tree, str(stem), show_depth=show_depth, use_title=False)
    dot_path = Path(f"{stem}.dot")

    # Graphviz's parsed graph is the authority, rather than a text match on the
    # serialized DOT: rank statements can silently invent additional nodes.
    parsed = subprocess.run(
        [graphviz_dot, "-Tjson", str(dot_path)], capture_output=True, text=True, check=False
    )
    assert parsed.returncode == 0, parsed.stderr
    graph = json.loads(parsed.stdout)
    nodes = [item for item in graph["objects"] if "height" in item]
    actual_names = [node["name"] for node in nodes]
    root = f"{basename}_none"
    expected_names = {root, f"{root}_left", f"{root}_right"}
    if show_depth:
        expected_names.update({"1", "2"})
    assert len(actual_names) == len(expected_names)
    assert set(actual_names) == expected_names

    names_by_id = {node["_gvid"]: node["name"] for node in nodes}
    actual_edges = {
        (names_by_id[edge["tail"]], names_by_id[edge["head"]]) for edge in graph["edges"]
    }
    expected_edges = {(root, f"{root}_left"), (root, f"{root}_right")}
    if show_depth:
        expected_edges.add(("1", "2"))
    assert actual_edges == expected_edges

    # Also verify a real exported image, including its node identity titles.
    image_path = tmp_path / "rendered.svg"
    rendered = subprocess.run(
        [graphviz_dot, "-Tsvg", "-o", str(image_path), str(dot_path)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert rendered.returncode == 0, rendered.stderr
    svg = ET.parse(image_path)
    ns = {"svg": "http://www.w3.org/2000/svg"}
    image_nodes = [
        group.find("svg:title", ns).text
        for group in svg.findall(".//svg:g", ns)
        if group.get("class") == "node"
    ]
    assert len(image_nodes) == len(expected_names)
    assert set(image_nodes) == expected_names


def _parse_graph(graphviz_dot, dot_path):
    parsed = subprocess.run(
        [graphviz_dot, "-Tjson", str(dot_path)], capture_output=True, text=True, check=False
    )
    assert parsed.returncode == 0, parsed.stderr
    graph = json.loads(parsed.stdout)
    nodes = [item for item in graph["objects"] if "height" in item]
    names_by_id = {node["_gvid"]: node["name"] for node in nodes}
    edges = {(names_by_id[edge["tail"]], names_by_id[edge["head"]]) for edge in graph["edges"]}
    return nodes, edges


@pytest.mark.parametrize(
    "basename",
    ["x" + "\\" * count + '"q' for count in range(6)] + ["x\nq", "x\rq", '軌跡_日本語\\"q'],
    ids=[
        "quote",
        "slash1-quote",
        "slash2-quote",
        "slash3-quote",
        "slash4-quote",
        "slash5-quote",
        "newline",
        "return",
        "unicode-slash-quote",
    ],
)
@pytest.mark.parametrize("show_depth", [False, True], ids=["no-depth", "depth"])
def test_lexical_filename_boundaries_preserve_tree_structure(
    tmp_path, graphviz_dot, three_node_tree, basename, show_depth
):
    stem = tmp_path / basename
    disstree2dot(three_node_tree, str(stem), show_depth=show_depth, use_title=False)
    dot_path = Path(f"{stem}.dot")
    nodes, edges = _parse_graph(graphviz_dot, dot_path)
    names = [node["name"] for node in nodes]
    depth_names = {"1", "2"} if show_depth else set()
    tree_names = set(names) - depth_names
    assert len(names) == (5 if show_depth else 3)
    assert len(tree_names) == 3
    assert depth_names.issubset(names)
    tree_edges = edges - ({("1", "2")} if show_depth else set())
    assert len(tree_edges) == 2
    roots = {source for source, _ in tree_edges}
    assert len(roots) == 1
    root = roots.pop()
    assert {target for _, target in tree_edges} == tree_names - {root}
    assert len(edges) == (3 if show_depth else 2)

    # Unsafe filename identities may use a documented encoded fallback. Its
    # spelling is internal, but repeated serialization must map consistently.
    repeated_stem = tmp_path / "repeat" / basename
    disstree2dot(three_node_tree, str(repeated_stem), show_depth=show_depth, use_title=False)
    repeated_nodes, repeated_edges = _parse_graph(graphviz_dot, Path(f"{repeated_stem}.dot"))
    assert {node["name"] for node in repeated_nodes} == set(names)
    assert repeated_edges == edges

    image_path = tmp_path / "lexical.svg"
    rendered = subprocess.run(
        [graphviz_dot, "-Tsvg", "-o", str(image_path), str(dot_path)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert rendered.returncode == 0, rendered.stderr
    svg = ET.parse(image_path)
    groups = svg.findall(".//{http://www.w3.org/2000/svg}g")
    assert sum(group.get("class") == "node" for group in groups) == len(names)


@pytest.mark.parametrize(
    "basenames",
    [("x\nq", "x\\nq"), ("x\rq", "x\\rq"), ('x\\"q', 'x\\\\"q')],
    ids=["newline-versus-slash-n", "return-versus-slash-r", "odd-versus-even-slash-quote"],
)
@pytest.mark.parametrize("show_depth", [False, True], ids=["no-depth", "depth"])
def test_distinct_lexical_filenames_have_distinct_tree_identifiers(
    tmp_path, graphviz_dot, three_node_tree, basenames, show_depth
):
    # Escaping a newline as literal backslash-n must not alias an actual
    # backslash-n filename. The same injective mapping must serve every edge
    # and rank reference, regardless of whether fallback encoding is needed.
    names_per_graph = []
    for basename in basenames:
        stem = tmp_path / basename
        disstree2dot(three_node_tree, str(stem), show_depth=show_depth, use_title=False)
        nodes, _ = _parse_graph(graphviz_dot, Path(f"{stem}.dot"))
        tree_names = {node["name"] for node in nodes} - ({"1", "2"} if show_depth else set())
        assert len(tree_names) == 3
        names_per_graph.append(tree_names)
    assert names_per_graph[0].isdisjoint(names_per_graph[1])


def test_reserved_identifier_prefix_cannot_alias_encoded_filename(graphviz_dot):
    from traminex.disstree_viz import _dot_identifier

    # A literal filename can match the transport spelling of another filename;
    # the serializer must keep those inputs distinct in the actual parsed graph.
    identifiers = ["x\nq", "_traminex_id_780a71", "ordinary_control"]
    declarations = [
        f'{_dot_identifier(name)} [label="marker{position}"];'
        for position, name in enumerate(identifiers)
    ]
    parsed = subprocess.run(
        [graphviz_dot, "-Tjson"],
        input="digraph tree {" + "\n".join(declarations) + "}",
        capture_output=True,
        text=True,
        check=False,
    )
    assert parsed.returncode == 0, parsed.stderr
    nodes = json.loads(parsed.stdout)["objects"]
    assert len(nodes) == 3
    assert {node["label"] for node in nodes} == {"marker0", "marker1", "marker2"}
    assert next(node["name"] for node in nodes if node["label"] == "marker2") == "ordinary_control"


@pytest.mark.parametrize(
    "basename", ['image\\"name', "image\nname"], ids=["slash-quote", "newline"]
)
@pytest.mark.parametrize("show_depth", [False, True], ids=["no-depth", "depth"])
def test_callback_png_paths_render_with_lexical_filename_boundaries(
    tmp_path, graphviz_dot, three_node_tree, basename, show_depth
):
    from matplotlib.image import imsave

    stem = tmp_path / basename
    emitted = []

    def emit_png(subset, output_path):
        emitted.append(output_path)
        imsave(output_path, np.ones((8, 8, 3)), format="png")

    disstree2dot(
        three_node_tree,
        str(stem),
        image_fun=emit_png,
        use_title=True,
        show_depth=show_depth,
    )
    assert len(emitted) == 3
    assert all(path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n") for path in emitted)
    dot_path = Path(f"{stem}.dot")
    nodes, _ = _parse_graph(graphviz_dot, dot_path)
    assert len(nodes) == (5 if show_depth else 3)
    image_path = tmp_path / "with-images.svg"
    rendered = subprocess.run(
        [graphviz_dot, "-Tsvg", "-o", str(image_path), str(dot_path)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    assert rendered.returncode == 0, rendered.stderr
    svg = ET.parse(image_path)
    images = svg.findall(".//{http://www.w3.org/2000/svg}image")
    assert len(images) == 3, rendered.stderr


@pytest.mark.parametrize(
    "basename", ['display\\"name', "display\nname"], ids=["slash-quote", "newline"]
)
def test_display_png_output_accepts_lexical_filename_boundaries(
    tmp_path, graphviz_dot, three_node_tree, basename, monkeypatch
):
    from traminex.disstree import disstreedisplay

    monkeypatch.setenv("GVPATH", graphviz_dot)
    filename = tmp_path / f"{basename}.png"
    output = disstreedisplay(
        three_node_tree,
        str(filename),
        with_quality=False,
        show_depth=True,
    )
    assert Path(output) == filename
    assert filename.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")


@pytest.mark.parametrize(
    "column", ['predictor\\"name', "predictor\nname"], ids=["slash-quote", "newline"]
)
def test_predictor_labels_render_with_lexical_boundaries(tmp_path, graphviz_dot, column):
    distances = np.array([[0, 1, 4, 4], [1, 0, 4, 4], [4, 4, 0, 1], [4, 4, 1, 0]], float)
    tree = disstree(
        distances,
        pd.DataFrame({column: ["a", "a", "b", "b"]}),
        R=0,
        min_size=1,
        max_depth=2,
    )
    stem = tmp_path / "label-tree"
    disstree2dot(tree, str(stem), use_title=False)
    dot_path = Path(f"{stem}.dot")
    nodes, _ = _parse_graph(graphviz_dot, dot_path)
    assert len(nodes) == 3
    root = next(node for node in nodes if node["name"] == "label-tree_none")
    assert "Split:" in root["label"]
    image_path = tmp_path / "labels.svg"
    rendered = subprocess.run(
        [graphviz_dot, "-Tsvg", "-o", str(image_path), str(dot_path)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert rendered.returncode == 0, rendered.stderr
    text = "".join(ET.parse(image_path).getroot().itertext())
    assert "predictor" in text and "name" in text
    if "\\" in column:
        assert column in text


@pytest.mark.parametrize("image_option", ["image_legend", "image_quality"])
@pytest.mark.parametrize(
    "relative_image",
    ['panel\\"name.png', 'unsafe\\"directory/panel.png', "panel.png\\", "panel.png"],
    ids=["slash-quote-basename", "slash-quote-parent", "trailing-backslash", "ordinary"],
)
def test_external_png_images_render_from_lexical_paths(
    tmp_path, graphviz_dot, three_node_tree, image_option, relative_image
):
    from matplotlib.image import imsave

    source = tmp_path / "external" / relative_image
    source.parent.mkdir(parents=True, exist_ok=True)
    imsave(source, np.full((8, 8, 3), 0.5), format="png")
    original_image = source.read_bytes()
    output_directory = tmp_path / "output"
    output_directory.mkdir()
    stem = output_directory / "external-panels"
    disstree2dot(three_node_tree, str(stem), **{image_option: source})
    dot_path = Path(f"{stem}.dot")
    nodes, _ = _parse_graph(graphviz_dot, dot_path)
    expected_panel = "node_image.legend" if image_option == "image_legend" else "node_qualityimage"
    assert len(nodes) == 4
    assert expected_panel in {node["name"] for node in nodes}

    image_path = output_directory / "external-panels.svg"
    rendered = subprocess.run(
        [graphviz_dot, "-Tsvg", "-o", str(image_path), str(dot_path)],
        cwd=output_directory,
        capture_output=True,
        text=True,
        check=False,
    )
    assert rendered.returncode == 0, rendered.stderr
    images = ET.parse(image_path).findall(".//{http://www.w3.org/2000/svg}image")
    assert len(images) == 1, rendered.stderr
    loaded_reference = images[0].get("{http://www.w3.org/1999/xlink}href")
    assert loaded_reference is not None
    assert (output_directory / loaded_reference).read_bytes() == original_image
    assert source.read_bytes() == original_image


@pytest.mark.parametrize("existing_parent", [False, True], ids=["new-parent", "existing-parent"])
@pytest.mark.parametrize("mode", ["plain", "callback", "image_legend", "image_quality"])
def test_nested_output_parent_supports_png_generation_and_external_panels(
    tmp_path, graphviz_dot, three_node_tree, existing_parent, mode
):
    from matplotlib.image import imsave

    output_directory = tmp_path / "nested" / "output"
    if existing_parent:
        output_directory.mkdir(parents=True)
    stem = output_directory / "tree"
    arguments = {}
    callback_rows = []
    callback_paths = []
    if mode == "callback":

        def emit_png(subset, output_path):
            callback_rows.append(subset.tolist())
            callback_paths.append(output_path)
            imsave(output_path, np.full((8, 8, 3), 0.5), format="png")

        arguments = {"image_fun": emit_png, "image_data": np.array([10, 11, 12, 13])}
    elif mode != "plain":
        source = tmp_path / 'panel\\"image.png'
        imsave(source, np.full((8, 8, 3), 0.5), format="png")
        arguments = {mode: source}

    disstree2dot(three_node_tree, str(stem), **arguments)
    dot_path = Path(f"{stem}.dot")
    assert dot_path.is_file()
    if mode == "callback":
        assert callback_rows == [[10, 11, 12, 13], [10, 11], [12, 13]]
        assert all(path.parent == output_directory for path in callback_paths)
        assert all(path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n") for path in callback_paths)
    rendered_path = output_directory / "tree.svg"
    rendered = subprocess.run(
        [graphviz_dot, "-Tsvg", "-o", str(rendered_path), str(dot_path)],
        cwd=output_directory,
        capture_output=True,
        text=True,
        check=False,
    )
    assert rendered.returncode == 0, rendered.stderr
    svg = ET.parse(rendered_path)
    image_nodes = svg.findall(".//{http://www.w3.org/2000/svg}image")
    expected_images = 0 if mode == "plain" else 3 if mode == "callback" else 1
    assert len(image_nodes) == expected_images, rendered.stderr
