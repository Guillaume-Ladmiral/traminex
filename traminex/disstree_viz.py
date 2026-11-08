"""Deterministic GraphViz serialization for dissimilarity regression trees."""

from __future__ import annotations

import hashlib
import os
import re
import shutil
import subprocess
import tempfile
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from itertools import count
from pathlib import Path
from typing import Any, List, Optional

import numpy as np
import pandas as pd

from ._deprecated import MISSING, resolve_deprecated_alias
from ._rformat import r_as_character

_DISPLAY_TEMP_DIRECTORY: tempfile.TemporaryDirectory | None = None
_DISPLAY_STEM_COUNTER = count()


# R's KP_MAX (format.c): the largest power of ten whose fixed-notation rounding
# ``scientific()`` still checks.
_R_KP_MAX = 27
# disstree2dotp's quality lines: format(x, digits = getOption("digits") - 2)
# with R's default getOption("digits") of 7.
_R_QUALITY_DIGITS = 5


def _fmt(x: float, digits: int) -> str:
    """R ``format(x, digits = digits)`` for one number (``formatReal``).

    R writes fixed notation unless it would be wider than scientific notation
    (``scipen = 0``); Python's ``"%.{d}g"`` switches to scientific notation as
    soon as a value has more than ``d`` integer digits (2026-09-22 audit S6-03).
    """
    value = float(x)
    if np.isnan(value):
        return "NaN"
    if np.isinf(value):
        return "Inf" if value > 0 else "-Inf"
    negative = int(value < 0)
    if value == 0:
        kpower, nsig, rounding_widens = 0, 1, False
    else:
        magnitude = abs(value)
        mantissa, exponent = format(magnitude, f".{digits - 1}e").split("e")
        nsig = max(len(mantissa.replace(".", "").rstrip("0")), 1)
        kpower = int(exponent)
        # 9996 at 3 digits is 1e+04 in scientific but 9996 in fixed notation.
        right_digits = min(max(digits - kpower, 0), _R_KP_MAX)
        rounding_widens = (
            0 < kpower <= _R_KP_MAX and magnitude < 10.0**kpower - 0.5 / 10.0**right_digits
        )
    left = kpower + 1 - int(rounding_widens)
    right = max(nsig - left, 0)
    fixed_width = negative + (left if left > 0 else 1) + right + int(right > 0)
    # R's third exponent digit; it can never tip the choice, because fixed
    # notation is at least 100 characters wide at such exponents.
    exponent_digits = 2 if kpower >= 100 or kpower <= -99 else 1
    scientific_width = negative + int(nsig > 1) + (nsig - 1) + 4 + exponent_digits
    if fixed_width <= scientific_width:
        return format(value, f".{right}f")
    return format(value, f".{nsig - 1}e")


def _dot_text(text: object) -> str:
    """Quote literal label text, including backslashes before quotes."""
    return str(text).replace("\\", "\\\\").replace('"', '\\"')


def _dot_path_text(path: str) -> str:
    """Retain Graphviz's literal path backslashes rather than label escapes."""
    return path.replace('"', '\\"')


def _dot_identifier(text: object, *, quoted: bool = True) -> str:
    """Encode filename-derived identities consistently, retaining safe ranks."""
    value = str(text)
    # DOT's scanner preserves most backslashes but cannot quote an odd run
    # immediately before a quote. Escaping control characters would also alias
    # literal backslash-n/r filenames. Use an injective reserved transport for
    # these names, including names that already start with the reserved prefix.
    prefix = "_traminex_id_"
    odd_quote = any(len(match.group(1)) % 2 for match in re.finditer(r'(\\+)"', value))
    if value.startswith(prefix) or odd_quote or any(ord(char) < 32 for char in value):
        value = prefix + value.encode("utf-8").hex()
    if not quoted and re.fullmatch(r"[A-Za-z_][A-Za-z_0-9]*", value):
        return value
    # DOT preserves ordinary backslashes in identifiers rather than decoding
    # them as Python/JSON escapes. Only embedded quotes need escaping here.
    escaped = value.replace('"', '\\"').replace("\n", "\\n").replace("\r", "\\r")
    return f'"{escaped}"'


def _stars(p: float) -> str:
    # R: as.character(cut(val, c(0, 0.001, 0.01, 0.05, 1), labels = c("***",
    # "**", "*", ""))) -- values outside (0, 1] fall out of every bin -> "NA".
    value = float(p)
    if not np.isfinite(value) or value <= 0 or value > 1:
        return "NA"
    if value <= 0.001:
        return "***"
    if value <= 0.01:
        return "**"
    if value <= 0.05:
        return "*"
    return ""


def _quality_lines(tree: Any) -> List[str]:
    """R ``disstree2dotp``'s quality text; its digits do not follow ``digits``."""
    stat = tree.adjustment.stat
    lines = []
    for name in ("Pseudo F", "Pseudo R2", "Levene"):
        value = _fmt(float(stat.loc[name, "t0"]), _R_QUALITY_DIGITS)
        stars = _stars(float(stat.loc[name, "p_value"]))
        # R: paste(n, ": ", value, " ", star, sep="") -- the separator space is
        # kept even when the star is empty.
        lines.append(f"{name}: {value} {stars}")
    return lines


@dataclass(frozen=True)
class _DotNode:
    name: str
    schedule: int
    content: str
    is_leaf: bool
    parent: str | None
    identifiers: tuple[Any, ...]
    title: str = " "


def _branch_labels(split: Any, digits: int) -> tuple[str, str]:
    """R's ``split_print`` (disstree2dot.R:400-414): the text of the two branches."""
    if split.breaks is not None:
        breaks = _fmt(split.breaks, digits)
        labels = [f"<= {breaks}", f"> {breaks}"]
    else:
        index = np.asarray(split.index)
        levels = [str(r_as_character(level)) for level in split.labels]
        labels = [
            "[ " + ",".join(level for level, side in zip(levels, index) if side == branch) + " ]"
            for branch in (1, 2)
        ]
    if split.na_group is not None:
        labels[int(split.na_group) - 1] += " with NA"
    return labels[0], labels[1]


def _collect_dot_nodes(tree: Any, filename: str, digits: int) -> list[_DotNode]:
    """Flatten a tree to immutable preorder node events.

    ``content`` is R's label for ``use.title = FALSE`` (node and split text run
    together, disstree2dot.R:392-394); ``title`` is the text R writes into the
    node image's title (branch label and node line in ``main``, split line in
    ``sub``; :365-383), which a node without an image shows as its label.
    """
    stem = Path(filename).name
    columns = list(tree.data.columns)
    records: list[_DotNode] = []
    stack: list[tuple[Any, str, str | None, str]] = [(tree.root, f"{stem}_none", None, "Root")]
    while stack:
        node, name, parent, branch = stack.pop()
        node_text = f"n: {_fmt(node.info['n'], digits)} s2: {_fmt(node.info['vardis'], digits)}"
        split_text = ""
        if node.split is not None:
            variable = _dot_text(columns[int(node.split.varindex)])
            split_text = f"Split: {variable} R2:{_fmt(node.split.info['R2'], digits)}"
        title = f"{_dot_text(branch)}\\n{node_text}" + (f"\\n{split_text}" if split_text else "")
        records.append(
            _DotNode(
                name=name,
                schedule=int(node.info["splitschedule"]),
                content=node_text + split_text,
                is_leaf=node.kids is None,
                parent=parent,
                identifiers=tuple(np.asarray(node.info["ind"]).tolist()),
                title=title,
            )
        )
        if node.kids is not None:
            left, right = _branch_labels(node.split, digits)
            stack.append((node.kids[1], f"{name}_right", name, right))
            stack.append((node.kids[0], f"{name}_left", name, left))
    return records


def _dot_image_node(name: str, image: str) -> str:
    return f'{_dot_identifier(name)}[shape=box, image="{_dot_path_text(image)}", imagescale=true, label=" "];'


def _build_dot_lines(
    tree: Any,
    filename: str,
    digits: int,
    use_title: bool,
    show_depth: bool,
    *,
    legend_label: Optional[str] = None,
    quality_lines: Optional[List[str]] = None,
    node_images: Optional[dict[str, str]] = None,
    legend_text: Optional[str] = None,
    image_legend: Optional[str] = None,
    image_quality: Optional[str] = None,
    node_records: Optional[list[_DotNode]] = None,
) -> List[str]:
    """Serialize already-collected events in TraMineR's stable DOT order."""
    records = _collect_dot_nodes(tree, filename, digits) if node_records is None else node_records
    images = node_images or {}
    lines = ["digraph distree{"]
    if image_legend is not None:
        lines.append(_dot_image_node("node_image.legend", image_legend))
    elif legend_label is not None:
        lines.append(f'"node_image.legend"[shape=box, label="{_dot_text(legend_label)}"];')

    for record in records:
        image = images.get(record.name)
        image = None if image is None else _dot_path_text(image)
        if not use_title:
            image_fragment = f' image="{image}", imagescale=true, ' if image is not None else ""
            lines.append(
                f'{_dot_identifier(record.name)}[shape=box, {image_fragment} label="{record.content}"];'
            )
        elif image is None:
            # R always draws node images and writes this text into their titles.
            # When no node has an image, the text is the box's own label
            # (TX-AUD-20261010-014); with images at the leaves only, the inner
            # boxes stay blank, as in R.
            label = record.title if not images else " "
            lines.append(f'{_dot_identifier(record.name)}[shape=box,  label="{label}"];')
        else:
            lines.append(
                f'{_dot_identifier(record.name)}[shape=box,  image="{image}", imagescale=true, '
                f' label=" ", xlabel="{record.content}"];'
            )

    if image_quality is not None:
        lines.append(_dot_image_node("node_qualityimage", image_quality))
    elif quality_lines is not None:
        label = "Global quality\\n" + "\\n".join(quality_lines)
        lines.append(f'"node_qualityimage"[shape=box, label="{label}"];')
    if legend_text is not None:
        lines.append(f'"node_legend.text"[shape=box, label=<{legend_text}>];')

    if show_depth:
        schedules = list(dict.fromkeys(record.schedule for record in records))
        for position, schedule in enumerate(schedules):
            spacer = " " if position == len(schedules) - 1 else ""
            lines.append(f"{schedule}{spacer}[shape=box, fontsize=20];")
        for schedule in schedules:
            names = "; ".join(
                _dot_identifier(record.name, quoted=False)
                for record in records
                if record.schedule == schedule
            )
            lines.append(f"{{ rank=same ; {schedule} ; {names} ;}}")
        lines.append(" -> ".join(map(str, sorted(schedules))) + " ;")

    for record in records:
        if record.parent is not None:
            lines.append(f"{_dot_identifier(record.parent)}->{_dot_identifier(record.name)};")
    lines.append("}")
    return lines


def _write_dot(lines: List[str], filename: str) -> str:
    dot_path = Path(f"{filename}.dot")
    dot_path.parent.mkdir(parents=True, exist_ok=True)
    dot_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(dot_path)


def _relative_dot_reference(path: str | Path, dot_path: Path) -> str:
    candidate = Path(path)
    if not candidate.is_absolute():
        return candidate.as_posix()
    return Path(os.path.relpath(candidate, dot_path.parent)).as_posix()


def _dot_image_reference(path: str | Path, dot_path: Path) -> str:
    """Copy paths DOT cannot quote to a deterministic adjacent image name."""
    reference = _relative_dot_reference(path, dot_path)
    unsafe = any(len(m.group(1)) % 2 for m in re.finditer(r'(\\+)"', reference))
    trailing = re.search(r"(\\+)$", reference)
    unsafe = unsafe or bool(trailing and len(trailing.group(1)) % 2)
    unsafe = unsafe or any(ord(char) < 32 for char in reference)
    source = Path(reference)
    if not source.is_absolute():
        source = dot_path.parent / source
    # Graphviz writes image filenames into SVG href attributes without XML
    # escaping. Existing images with XML punctuation need the same safe alias.
    unsafe = unsafe or (source.is_file() and any(char in reference for char in '"&<>'))
    if not unsafe:
        return reference
    # The suffix selects the image decoder. Hashing avoids filename-length
    # expansion and distinguishes callback, legend and quality image paths.
    suffix = source.suffix
    if not re.fullmatch(r"\.[A-Za-z0-9]+", suffix):
        suffix = ""
    alias = "_traminex_image_" + hashlib.sha256(reference.encode("utf-8")).hexdigest() + suffix
    shutil.copyfile(source, dot_path.parent / alias)
    return alias


def _image_positions(records: list[_DotNode]) -> dict[str, list[int]]:
    """Resolve stored case identifiers to positional indexes once."""
    if not records:
        return {}
    root_ids = np.asarray(records[0].identifiers, dtype=int)
    n_cases = root_ids.size
    one_based = bool(
        n_cases and np.array_equal(np.sort(root_ids), np.arange(1, n_cases + 1, dtype=int))
    )
    offset = 1 if one_based else 0
    return {
        record.name: [int(identifier) - offset for identifier in record.identifiers]
        for record in records
    }


def _subset_image_data(image_data: object, positions: list[int]) -> object:
    if isinstance(image_data, (pd.DataFrame, pd.Series)):
        return image_data.iloc[positions].copy()
    if isinstance(image_data, np.ndarray):
        return image_data[np.asarray(positions, dtype=int)].copy()
    if isinstance(image_data, Sequence) and not isinstance(image_data, (str, bytes)):
        return [image_data[position] for position in positions]
    raise TypeError("image_data must support positional row selection")


def _row_count(image_data: object) -> int:
    try:
        return len(image_data)  # type: ignore[arg-type]
    except TypeError as exc:
        raise TypeError("image_data must have one row per root observation") from exc


def _emit_dot(
    tree: Any,
    filename: str,
    *,
    digits: int,
    image_fun: Callable[[object, Path], None] | None = None,
    image_data: object | None = None,
    only_leaf: bool = False,
    image_format: str = "png",
    use_title: bool = True,
    show_depth: bool = False,
    legend_text: str | None = None,
    image_legend: str | Path | None = None,
    image_quality: str | Path | None = None,
    legend_label: str | None = None,
    quality_lines: list[str] | None = None,
) -> None:
    records = _collect_dot_nodes(tree, filename, digits)
    dot_path = Path(f"{filename}.dot")
    dot_path.parent.mkdir(parents=True, exist_ok=True)
    positions = _image_positions(records)
    node_images: dict[str, str] = {}
    if image_fun is not None:
        if image_data is None:
            image_data = pd.DataFrame({0: records[0].identifiers})
        expected_rows = len(records[0].identifiers)
        if _row_count(image_data) != expected_rows:
            raise ValueError(
                f"image_data must have {expected_rows} rows, one for each root observation"
            )
        for record in records:
            if only_leaf and not record.is_leaf:
                continue
            output_path = dot_path.parent / f"{record.name}.{image_format}"
            subset = _subset_image_data(image_data, positions[record.name])
            image_fun(subset, output_path)
            if not output_path.is_file():
                raise RuntimeError("image callback must create a regular file")
            node_images[record.name] = _dot_image_reference(output_path, dot_path)

    legend_reference = (
        None if image_legend is None else _dot_image_reference(image_legend, dot_path)
    )
    quality_reference = (
        None if image_quality is None else _dot_image_reference(image_quality, dot_path)
    )
    lines = _build_dot_lines(
        tree,
        filename,
        digits,
        use_title,
        show_depth,
        legend_label=legend_label,
        quality_lines=quality_lines,
        node_images=node_images,
        legend_text=legend_text,
        image_legend=legend_reference,
        image_quality=quality_reference,
        node_records=records,
    )
    _write_dot(lines, filename)


def disstree2dot(
    tree: Any,
    filename: str,
    *,
    digits: int = 3,
    image_fun: Callable[[object, Path], None] | None = None,
    image_data: object | None = None,
    only_leaf: bool = False,
    image_format: str = "png",
    use_title: bool = True,
    show_depth: bool = False,
    legend_text: str | None = None,
    image_legend: str | Path | None = None,
    image_quality: str | Path | None = None,
    imagefunc=MISSING,
    imagedata=MISSING,
    imgLeafOnly=MISSING,
    imageext=MISSING,
    legendtext=MISSING,
    legendimage=MISSING,
    qualityimage=MISSING,
    showdepth=MISSING,
) -> None:
    """Write a deterministic TraMineR-compatible ``.dot`` tree description."""
    image_fun = resolve_deprecated_alias(
        image_fun,
        imagefunc,
        canonical_name="image_fun",
        alias_name="imagefunc",
        function_name="disstree2dot",
        default=None,
    )
    image_data = resolve_deprecated_alias(
        image_data,
        imagedata,
        canonical_name="image_data",
        alias_name="imagedata",
        function_name="disstree2dot",
        default=None,
    )
    only_leaf = resolve_deprecated_alias(
        only_leaf,
        imgLeafOnly,
        canonical_name="only_leaf",
        alias_name="imgLeafOnly",
        function_name="disstree2dot",
        default=False,
    )
    image_format = resolve_deprecated_alias(
        image_format,
        imageext,
        canonical_name="image_format",
        alias_name="imageext",
        function_name="disstree2dot",
        default="png",
    )
    legend_text = resolve_deprecated_alias(
        legend_text,
        legendtext,
        canonical_name="legend_text",
        alias_name="legendtext",
        function_name="disstree2dot",
        default=None,
    )
    image_legend = resolve_deprecated_alias(
        image_legend,
        legendimage,
        canonical_name="image_legend",
        alias_name="legendimage",
        function_name="disstree2dot",
        default=None,
    )
    image_quality = resolve_deprecated_alias(
        image_quality,
        qualityimage,
        canonical_name="image_quality",
        alias_name="qualityimage",
        function_name="disstree2dot",
        default=None,
    )
    show_depth = resolve_deprecated_alias(
        show_depth,
        showdepth,
        canonical_name="show_depth",
        alias_name="showdepth",
        function_name="disstree2dot",
        default=False,
    )
    _emit_dot(
        tree,
        filename,
        digits=digits,
        image_fun=image_fun,
        image_data=image_data,
        only_leaf=bool(only_leaf),
        image_format=str(image_format),
        use_title=use_title,
        show_depth=bool(show_depth),
        legend_text=legend_text,
        image_legend=image_legend,
        image_quality=image_quality,
    )


def disstree2dotp(
    tree: Any,
    filename: str,
    *,
    digits: int = 3,
    with_quality: bool = True,
    use_title: bool = True,
    show_depth: bool = False,
    withquality=MISSING,
) -> None:
    """Write a discrepancy tree to DOT with optional tree-quality lines."""
    from .disstree import DissTree

    if not isinstance(tree, DissTree):
        raise TypeError("tree must be a DissTree (disstree or seqtree)")
    with_quality = resolve_deprecated_alias(
        with_quality,
        withquality,
        canonical_name="with_quality",
        alias_name="withquality",
        function_name="disstree2dotp",
        default=True,
    )
    _emit_dot(
        tree,
        filename,
        digits=digits,
        use_title=use_title,
        show_depth=show_depth,
        quality_lines=_quality_lines(tree) if with_quality else None,
    )


def _legend_label(seqdata: Any) -> str:
    from .sequences import alphabet

    return "Legend: " + ", ".join(map(str, alphabet(seqdata)))


def seqtree2dot(
    tree: Any,
    filename: str,
    *,
    seqdata: Any = None,
    digits: int = 3,
    with_quality: bool = False,
    with_legend: bool = True,
    use_title: bool = True,
    show_depth: bool = False,
    withlegend=MISSING,
    withquality=MISSING,
) -> None:
    """Write a sequence regression tree to DOT with its legend and quality."""
    from .disstree import DissTree

    if not isinstance(tree, DissTree):
        raise TypeError("tree must be a DissTree (disstree or seqtree)")
    with_legend = resolve_deprecated_alias(
        with_legend,
        withlegend,
        canonical_name="with_legend",
        alias_name="withlegend",
        function_name="seqtree2dot",
        default=True,
    )
    with_quality = resolve_deprecated_alias(
        with_quality,
        withquality,
        canonical_name="with_quality",
        alias_name="withquality",
        function_name="seqtree2dot",
        default=False,
    )
    if seqdata is None:
        seqdata = tree.info.get("object")
    if seqdata is None:
        raise ValueError("seqdata is required for a sequence tree")
    _emit_dot(
        tree,
        filename,
        digits=digits,
        use_title=use_title,
        show_depth=show_depth,
        legend_label=_legend_label(seqdata) if with_legend else None,
        quality_lines=_quality_lines(tree) if with_quality else None,
    )


def _dot_candidates() -> List[Path]:
    """The paths GVPATH can name for ``dot``: GraphViz's installation directory
    (``<dir>/bin/dot``), its ``bin`` directory (``<dir>/dot``) or the executable.

    Deliberate divergence: R honours GVPATH only through Windows-specific
    ``\\bin\\dot.exe`` path construction; traminex resolves it on every
    platform (development is Linux-only, AGENTS.md; TX-AUD-20261010-082).
    """
    configured = os.environ.get("GVPATH")
    if not configured:
        return []
    base = Path(configured)
    if base.is_dir():
        return [base / "bin" / "dot", base / "dot"]
    return [base]


def _locate_dot() -> Optional[str]:
    for candidate in _dot_candidates():
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
    return shutil.which("dot")


def _render_dot(dot_path: str, out_path: str, image_format: str) -> None:
    executable = _locate_dot()
    if executable is None:
        configured = os.environ.get("GVPATH")
        setting = f"GVPATH={configured!r}" if configured else "GVPATH is not set"
        tried = [str(path) for path in _dot_candidates()] + ["'dot' on PATH"]
        raise RuntimeError(
            f"GraphViz 'dot' was not found ({setting}; tried {', '.join(tried)}). Install "
            "GraphViz or set GVPATH to its installation directory, its bin directory or the "
            f"dot executable. The DOT file is at {os.path.abspath(dot_path)}"
        )
    # Absolute paths: a relative name starting with "-" was read by dot as an
    # option (2026-09-22 audit S6-17).
    completed = subprocess.run(
        [
            executable,
            f"-T{image_format}",
            "-o",
            os.path.abspath(out_path),
            os.path.abspath(dot_path),
        ],
        # DOT image references are relative to the emitted file. Graphviz
        # otherwise resolves them from the caller's unrelated directory and
        # silently omits callback images even with a successful exit status.
        cwd=Path(dot_path).absolute().parent,
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        detail = completed.stderr.strip() or completed.stdout.strip()
        raise RuntimeError(f"GraphViz rendering failed: {detail}")
    if not Path(out_path).is_file():
        raise RuntimeError("GraphViz did not create the requested output file")


def _display_stem(filename: Optional[str], prefix: str) -> str:
    if filename is not None:
        return filename
    global _DISPLAY_TEMP_DIRECTORY
    if _DISPLAY_TEMP_DIRECTORY is None:
        _DISPLAY_TEMP_DIRECTORY = tempfile.TemporaryDirectory(prefix="traminex-display-")
    directory = Path(_DISPLAY_TEMP_DIRECTORY.name)
    return str(directory / f"{prefix}-{next(_DISPLAY_STEM_COUNTER)}")


def _display_paths(filename: Optional[str], prefix: str, image_format: str) -> tuple[str, str]:
    """``(stem, image path)`` for the display functions.

    As in R (``file.copy(..., filename)``), a ``filename`` that already ends
    with the image format's extension is the image path itself; any other name
    is a stem that gains the extension (2026-09-22 audit S6-18). The ``.dot``
    file is written at ``stem + ".dot"``.
    """
    if filename is not None:
        path = Path(filename)
        if path.suffix.lower() == f".{image_format}".lower():
            return str(path.with_suffix("")), str(path)
    stem = _display_stem(filename, prefix)
    return stem, f"{stem}.{image_format}"


def disstreedisplay(
    tree: Any,
    filename: Optional[str] = None,
    *,
    image_format: str = "png",
    with_quality: bool = True,
    show_depth: bool = False,
    digits: int = 3,
    imageformat=MISSING,
    withquality=MISSING,
    showdepth=MISSING,
) -> str:
    """Render a discrepancy-tree DOT graph and return the output image path.

    ``filename`` names the image, as in R, when it ends with the extension of
    ``image_format`` (``"tree.png"``); otherwise it is a stem that gains the
    extension (``"tree"`` -> ``tree.png``). The ``.dot`` file is kept beside the
    image. Without ``filename`` both go to a session temporary directory.
    ``dot`` renders ``image_format`` directly (R renders PNG or JPEG and
    converts other formats with ImageMagick). ``dot`` is found through the
    ``GVPATH`` environment variable, which may name GraphViz's installation
    directory, its ``bin`` directory or the executable, or else on ``PATH``;
    when it is missing the error names what was tried and the ``.dot`` file
    stays for rendering by hand.
    """
    from .disstree import DissTree

    if not isinstance(tree, DissTree):
        raise TypeError("tree must be a DissTree (disstree or seqtree)")
    image_format = resolve_deprecated_alias(
        image_format,
        imageformat,
        canonical_name="image_format",
        alias_name="imageformat",
        function_name="disstreedisplay",
        default="png",
    )
    with_quality = resolve_deprecated_alias(
        with_quality,
        withquality,
        canonical_name="with_quality",
        alias_name="withquality",
        function_name="disstreedisplay",
        default=True,
    )
    show_depth = resolve_deprecated_alias(
        show_depth,
        showdepth,
        canonical_name="show_depth",
        alias_name="showdepth",
        function_name="disstreedisplay",
        default=False,
    )
    stem, output = _display_paths(filename, "disstree", str(image_format))
    disstree2dotp(
        tree,
        stem,
        digits=digits,
        with_quality=bool(with_quality),
        show_depth=bool(show_depth),
    )
    _render_dot(f"{stem}.dot", output, str(image_format))
    return output


def seqtreedisplay(
    tree: Any,
    filename: Optional[str] = None,
    *,
    seqdata: Any = None,
    image_format: str = "png",
    with_quality: bool = True,
    with_legend: bool = True,
    show_depth: bool = False,
    digits: int = 3,
    imageformat=MISSING,
    withquality=MISSING,
    withlegend=MISSING,
    showdepth=MISSING,
) -> str:
    """Render a sequence-tree DOT graph and return the output image path.

    ``filename`` follows the :func:`disstreedisplay` rule: an image path when
    it ends with the extension of ``image_format``, otherwise a stem; ``dot``
    is found as :func:`disstreedisplay` finds it (``GVPATH`` or ``PATH``).
    """
    from .disstree import DissTree

    if not isinstance(tree, DissTree):
        raise TypeError("tree must be a DissTree (disstree or seqtree)")
    image_format = resolve_deprecated_alias(
        image_format,
        imageformat,
        canonical_name="image_format",
        alias_name="imageformat",
        function_name="seqtreedisplay",
        default="png",
    )
    with_quality = resolve_deprecated_alias(
        with_quality,
        withquality,
        canonical_name="with_quality",
        alias_name="withquality",
        function_name="seqtreedisplay",
        default=True,
    )
    with_legend = resolve_deprecated_alias(
        with_legend,
        withlegend,
        canonical_name="with_legend",
        alias_name="withlegend",
        function_name="seqtreedisplay",
        default=True,
    )
    show_depth = resolve_deprecated_alias(
        show_depth,
        showdepth,
        canonical_name="show_depth",
        alias_name="showdepth",
        function_name="seqtreedisplay",
        default=False,
    )
    stem, output = _display_paths(filename, "seqtree", str(image_format))
    seqtree2dot(
        tree,
        stem,
        seqdata=seqdata,
        digits=digits,
        with_quality=bool(with_quality),
        with_legend=bool(with_legend),
        show_depth=bool(show_depth),
    )
    _render_dot(f"{stem}.dot", output, str(image_format))
    return output
