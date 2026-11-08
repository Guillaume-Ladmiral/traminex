"""Every comparison in a golden test states its tolerance (TX-AUD-20261019-075, TX-AUD-20261011-065).

``tests/r_goldens_utils.py`` sets the contract for comparing with R: ``R_GOLDEN_ATOL`` (1e-10)
and ``assert_matches_golden``. The 2026-10-19 audit counted comparisons in the golden tests
that applied a library default instead -- ``assert_allclose``'s rtol 1e-7, ``np.allclose``'s
rtol 1e-5 and atol 1e-8, ``pytest.approx``'s rel 1e-6, pandas' rtol 1e-5 -- none of which a
reader sees at the call. Each such call now states its tolerance, the contract's or another
value with its reason beside it.
"""

from __future__ import annotations

import ast
from pathlib import Path

TESTS = Path(__file__).resolve().parent


def _unstated(name: str | None, keywords: set[str], exact: bool = False) -> bool:
    """Whether a comparison call leaves a tolerance at its library default.

    ``exact`` says the call passes ``check_exact=True``: pandas then compares exactly; with
    ``check_exact=False`` its rtol 1e-5 applies unless rtol and atol are given (batch-16
    review M-6)."""
    if name == "assert_allclose":
        return "rtol" not in keywords
    if name in ("allclose", "isclose"):
        return not {"rtol", "atol"} <= keywords
    if name == "approx":
        return not keywords & {"rel", "abs"}
    if name in ("assert_frame_equal", "assert_series_equal"):
        return not exact and not {"rtol", "atol"} <= keywords
    if name == "assert_array_almost_equal":
        return "decimal" not in keywords
    if name == "assert_almost_equal":  # NumPy's decimal 7, or pandas' rtol and atol
        return "decimal" not in keywords and not {"rtol", "atol"} <= keywords
    if name == "assert_approx_equal":  # NumPy's 7 significant digits
        return "significant" not in keywords
    return False


def _default_tolerance_calls(path: Path) -> list[str]:
    found = []
    for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
        if not isinstance(node, ast.Call):
            continue
        function = node.func
        name = (
            function.attr
            if isinstance(function, ast.Attribute)
            else function.id if isinstance(function, ast.Name) else None
        )
        keywords = {keyword.arg for keyword in node.keywords}
        if None in keywords:  # **options: the tolerance is not visible here
            continue
        exact = any(
            keyword.arg == "check_exact"
            and isinstance(keyword.value, ast.Constant)
            and keyword.value.value is True
            for keyword in node.keywords
        )
        if _unstated(name, keywords, exact):
            found.append((node.lineno, f"{path.name}:{node.lineno} {name}"))
    return [entry for _, entry in sorted(found)]


def test_golden_comparisons_state_their_tolerance() -> None:
    found = [
        call
        for path in sorted(TESTS.glob("test_r_goldens*.py"))
        for call in _default_tolerance_calls(path)
    ]

    assert found == [], f"{len(found)} comparisons use a library default tolerance:\n" + "\n".join(
        found
    )


def test_the_census_flags_a_default_tolerance(tmp_path) -> None:
    path = tmp_path / "test_r_goldens_example.py"
    path.write_text(
        "np.testing.assert_allclose(a, b, atol=1e-10)\n"
        "assert x == pytest.approx(1.0)\n"
        "pd.testing.assert_frame_equal(left, right)\n"
        "np.testing.assert_allclose(a, b, rtol=0, atol=R_GOLDEN_ATOL)\n"
        "assert x == pytest.approx(1.0, abs=R_GOLDEN_ATOL)\n"
        "pd.testing.assert_frame_equal(left, right, check_exact=True)\n",
        encoding="utf-8",
    )

    assert _default_tolerance_calls(path) == [
        "test_r_goldens_example.py:1 assert_allclose",
        "test_r_goldens_example.py:2 approx",
        "test_r_goldens_example.py:3 assert_frame_equal",
    ]


def _negated_equality_tolerance_calls(path: Path) -> list[str]:
    """'Values differ' checks (``not np.allclose``/``np.isclose``, ``!= pytest.approx``) that
    use the equality tolerance R_GOLDEN_ATOL: such a check passes on a difference of 1e-10,
    so it is weaker the tighter the tolerance (batch-16 review M-5)."""

    def uses_contract(call: ast.Call) -> bool:
        return any(
            isinstance(keyword.value, ast.Name) and keyword.value.id == "R_GOLDEN_ATOL"
            for keyword in call.keywords
        )

    def name_of(call: ast.Call) -> str | None:
        function = call.func
        return (
            function.attr if isinstance(function, ast.Attribute) else getattr(function, "id", None)
        )

    found = []
    for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
        calls = []
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.Not):
            if isinstance(node.operand, ast.Call) and name_of(node.operand) in (
                "allclose",
                "isclose",
            ):
                calls.append(node.operand)
        if isinstance(node, ast.Compare) and any(isinstance(op, ast.NotEq) for op in node.ops):
            calls += [
                side
                for side in [node.left, *node.comparators]
                if isinstance(side, ast.Call) and name_of(side) == "approx"
            ]
        found += [
            (call.lineno, f"{path.name}:{call.lineno}") for call in calls if uses_contract(call)
        ]
    return [entry for _, entry in sorted(found)]


def test_values_differ_checks_state_a_difference_threshold() -> None:
    found = [
        call
        for path in sorted(TESTS.glob("test_r_goldens*.py"))
        for call in _negated_equality_tolerance_calls(path)
    ]

    assert found == [], "'values differ' checks with the equality tolerance:\n" + "\n".join(found)


def test_the_census_flags_a_negated_equality_tolerance(tmp_path) -> None:
    path = tmp_path / "test_r_goldens_example.py"
    path.write_text(
        "assert not np.allclose(a, b, rtol=0, atol=R_GOLDEN_ATOL)\n"
        "assert x != pytest.approx(y, rel=0, abs=R_GOLDEN_ATOL)\n"
        "assert not np.allclose(a, b, rtol=0, atol=R_DIFFERENCE_ATOL)\n"
        "np.testing.assert_allclose(a, b, rtol=0, atol=R_GOLDEN_ATOL)\n",
        encoding="utf-8",
    )

    assert _negated_equality_tolerance_calls(path) == [
        "test_r_goldens_example.py:1",
        "test_r_goldens_example.py:2",
    ]


def test_the_census_flags_inexact_frames_and_the_other_numpy_helpers(tmp_path) -> None:
    """``check_exact=False`` alone keeps pandas' rtol 1e-5; ``assert_almost_equal`` (decimal 7)
    and ``assert_approx_equal`` (7 significant digits) have defaults too (batch-16 review M-6)."""
    path = tmp_path / "test_r_goldens_example.py"
    path.write_text(
        "pd.testing.assert_frame_equal(left, right, check_exact=False)\n"
        "np.testing.assert_almost_equal(a, b)\n"
        "np.testing.assert_approx_equal(a, b)\n"
        "pd.testing.assert_frame_equal(left, right, check_exact=False, rtol=0, atol=R_GOLDEN_ATOL)\n"
        "pd.testing.assert_series_equal(left, right, check_exact=True)\n"
        "np.testing.assert_almost_equal(a, b, decimal=10)\n"
        "np.testing.assert_approx_equal(a, b, significant=12)\n",
        encoding="utf-8",
    )

    assert _default_tolerance_calls(path) == [
        "test_r_goldens_example.py:1 assert_frame_equal",
        "test_r_goldens_example.py:2 assert_almost_equal",
        "test_r_goldens_example.py:3 assert_approx_equal",
    ]
