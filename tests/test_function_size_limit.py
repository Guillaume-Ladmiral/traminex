"""No function of the package grows past its recorded size (TX-AUD-20261019-072).

The 2026-10-19 audit found 22 functions longer than 200 lines or with an approximate
cyclomatic complexity above 50; ``seqdist`` grew from 1,179 to 1,244 lines during the repair
that followed. The owner chose a growth limit now and the restructuring of ``seqdist`` and
``pairwise`` later (2026-10-24). Every other function must stay within 200 lines and
complexity 50. A listed function must measure exactly as recorded: growing fails, and a
shrink fails until the table records the new, lower size. The test cannot stop one commit
that raises a recorded size together with the code: that shows in the diff as an edit to
RECORDED, and review must catch it.
Length and complexity are counted as the audit counted them (lines from ``def`` to the
last line; 1 plus the decision points: if, for, while, try, except, conditional
expressions, comprehensions, with, assert, each extra operand of and/or, match cases).
"""

from __future__ import annotations

import ast
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
MAX_LENGTH = 200
MAX_COMPLEXITY = 50

#: (file, function name) -> (length, complexity), measured at 3579bf12.
RECORDED = {
    ("traminex/align_p2.py", "seqalign"): (223, 39),
    ("traminex/diss.py", "dissassoc"): (328, 39),
    ("traminex/diss_p1.py", "dissrep"): (365, 42),
    ("traminex/disstree.py", "disstree"): (242, 32),
    ("traminex/distances/__init__.py", "pairwise"): (277, 79),
    ("traminex/distances/__init__.py", "seqdist"): (1244, 336),
    ("traminex/distances/_omstran.py", "omstran_transform_rows"): (147, 55),
    ("traminex/distances/omloc.py", "omloc"): (256, 42),
    ("traminex/distances/python_kernels.py", "_pairwise_rows"): (270, 76),
    ("traminex/distances/python_kernels.py", "pairwise"): (412, 112),
    ("traminex/events/_object.py", "seqecreate"): (231, 66),
    ("traminex/format.py", "TSE_to_STS"): (
        231,
        61,
    ),  # 2026-11-07: rows grouped by _rows_by_person (F-02)
    ("traminex/format.py", "_hspell_horizontal_sts"): (180, 61),
    ("traminex/format.py", "seqdef"): (
        433,
        152,
    ),  # 2026-11-07: states renamed by _rename_states (F-01)
    ("traminex/multidomain_cat.py", "seqMD"): (221, 63),
    ("traminex/seqformat.py", "_normalize_spell_input"): (305, 100),
    ("traminex/seqformat.py", "seqformat"): (211, 37),
    ("traminex/stats_complexity.py", "_seqidegrad_correction"): (180, 54),
    ("traminex/stats_p2.py", "seqpolyads"): (229, 49),
    ("traminex/stats_p2.py", "seqrecode"): (243, 63),
    ("traminex/stats_utils.py", "seqindic"): (355, 68),
    ("traminex/visualization/plots/survival.py", "seqsplot"): (297, 76),
}

_DECISIONS = (
    ast.If,
    ast.For,
    ast.AsyncFor,
    ast.While,
    ast.Try,
    ast.ExceptHandler,
    ast.IfExp,
    ast.comprehension,
    ast.With,
    ast.Assert,
)


def _complexity(function: ast.AST) -> int:
    nodes = list(ast.walk(function))
    return (
        1
        + sum(isinstance(node, _DECISIONS) for node in nodes)
        + sum(len(node.values) - 1 for node in nodes if isinstance(node, ast.BoolOp))
        + sum(len(node.cases) for node in nodes if isinstance(node, ast.Match))
    )


def _measured() -> list[tuple[str, str, int, int]]:
    rows = []
    for path in sorted(REPO.glob("traminex/**/*.py")):
        relative = path.relative_to(REPO).as_posix()
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                length = node.end_lineno - node.lineno + 1
                rows.append((relative, node.name, length, _complexity(node)))
    return rows


def test_no_function_passes_its_size_limit() -> None:
    rows = _measured()
    found: dict[tuple[str, str], list[tuple[int, int]]] = {}
    for relative, name, length, complexity in rows:
        found.setdefault((relative, name), []).append((length, complexity))

    problems = []
    for key, sizes in found.items():
        if key in RECORDED:
            if sizes != [RECORDED[key]]:
                problems.append(f"{key}: measured {sizes}, recorded {RECORDED[key]}")
        else:
            problems += [
                f"{key}: {length} lines, complexity {complexity}"
                for length, complexity in sizes
                if length > MAX_LENGTH or complexity > MAX_COMPLEXITY
            ]
    problems += [f"{key}: recorded but not found" for key in RECORDED if key not in found]

    assert problems == []


def test_the_documents_say_review_must_catch_a_raised_recorded_size() -> None:
    # One commit can raise a RECORDED size together with the code and this test passes, so the
    # limit does not only tighten (batch-18 review M-6).
    contributing = (REPO / "CONTRIBUTING.md").read_text(encoding="utf-8")
    for name, text in (("CONTRIBUTING.md", contributing), ("module docstring", __doc__ or "")):
        flat = " ".join(text.split())
        assert "only tightens" not in flat, name
        assert "review must catch it" in flat, name
