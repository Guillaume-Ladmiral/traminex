from __future__ import annotations

import importlib.metadata as metadata
import json
import os
import shutil
import subprocess
import sys
import textwrap
import tomllib
from pathlib import Path

import numpy as np
import pytest
import yaml
from packaging.requirements import Requirement
from packaging.specifiers import SpecifierSet
from packaging.utils import canonicalize_name
from packaging.version import Version

REPO = Path(__file__).resolve().parents[1]

# Tools whose declarations must never sit below the maintained environment,
# keyed by the pre-commit hook repository URL. The environment is the reference;
# the codebase tracks it (see test_no_formatter_declaration_sits_below_the_installed_toolchain).
PINNED_HOOK_TOOLS = {
    "https://github.com/psf/black": "black",
    "https://github.com/astral-sh/ruff-pre-commit": "ruff",
}


def _maintained_python() -> Path | None:
    """The interpreter the Makefile's gates actually run.

    ``AGENTS.md`` makes one environment the reference for development, linting
    and packaging, and the Makefile names it in ``PYTHON ?=``. The
    pre-commit-revision check below is a claim about ``make format-check``, so
    it has to compare against *that* interpreter and not whichever one happens
    to be running the tests -- a throwaway 3.11 matrix environment resolves the
    newest ruff, which is neither the reference nor what the gate runs.

    ``make release-check`` passes its own ``PYTHON`` in ``TRAMINEX_MAKE_PYTHON``,
    which is read first: CI calls it with ``PYTHON=python``, where the default path
    below does not exist and the check skipped (TX-AUD-20261011-076). A command name
    is resolved on ``PATH``, as make's shell resolves it.
    """
    given = os.environ.get("TRAMINEX_MAKE_PYTHON")
    if given:
        located = shutil.which(given) if os.sep not in given else given
        return Path(located) if located and Path(located).is_file() else None
    for line in (REPO / "Makefile").read_text(encoding="utf-8").splitlines():
        if line.startswith("PYTHON ?="):
            # The default is written from make's `$(HOME)`, which make takes from the
            # environment; expand it the same way, or the check skips (2026-11-08).
            value = line.split("=", 1)[1].strip().replace("$(HOME)", os.environ.get("HOME", ""))
            candidate = Path(value)
            return candidate if candidate.is_file() else None
    return None


def _maintained_tool_version(tool: str) -> Version | None:
    interpreter = _maintained_python()
    if interpreter is None:
        return None
    completed = subprocess.run(
        [
            str(interpreter),
            "-c",
            f"import importlib.metadata as m; print(m.version({tool!r}))",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        return None
    return Version(completed.stdout.strip())


EXPECTED_REQUIREMENTS = {
    "core": {
        "numpy": ">=1.26",
        "pandas": ">=2.1",
        "scipy": ">=1.11.1",
    },
    "plot": {"matplotlib": ">=3.6", "seaborn": ">=0.12"},
    "perf": {"numba": ">=0.59"},
    "extras": {"scikit-learn": ">=1.3"},
    # `tree` is deliberately empty: tree rendering shells out to the system
    # GraphViz `dot`, which pip cannot install (cycle-6 F-004).
    "tree": {},
    "dev": {
        # Floors at the versions the maintained `datasci` environment carries.
        # They must not pin below it -- that would force the reference
        # environment to downgrade instead of adapting the codebase to it.
        "black": ">=26.5.1",
        "build": ">=1.2",
        "mypy": ">=1.10",
        "pre-commit": ">=3.7.0",
        "pyyaml": ">=6.0",
        "pytest": ">=7.0",
        "pytest-cov": ">=4.0",
        "ruff": ">=0.16.9",
        "scikit-build-core": ">=0.11.0",
    },
}


def _requirements_by_extra() -> dict[str, dict[str, Requirement]]:
    dist = metadata.distribution("traminex")
    requirements = [Requirement(raw) for raw in dist.requires or ()]
    grouped: dict[str, dict[str, Requirement]] = {key: {} for key in EXPECTED_REQUIREMENTS}
    for requirement in requirements:
        matching_extras = [
            extra
            for extra in EXPECTED_REQUIREMENTS
            if extra != "core"
            and requirement.marker is not None
            and requirement.marker.evaluate({"extra": extra})
        ]
        group = matching_extras[0] if matching_extras else "core"
        grouped[group][canonicalize_name(requirement.name)] = requirement
    return grouped


def _run_isolated(code: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-I", "-c", code],
        check=False,
        capture_output=True,
        text=True,
    )


def test_installed_metadata_declares_supported_python_and_dependency_floors():
    dist = metadata.distribution("traminex")
    requires_python = dist.metadata["Requires-Python"]

    assert requires_python == ">=3.11"
    assert SpecifierSet(requires_python).contains(
        f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
    )

    grouped = _requirements_by_extra()
    assert {group: set(requirements) for group, requirements in grouped.items()} == {
        group: set(requirements) for group, requirements in EXPECTED_REQUIREMENTS.items()
    }
    for group, expected in EXPECTED_REQUIREMENTS.items():
        for name, declared in expected.items():
            requirement = grouped[group][name]
            assert str(requirement.specifier) == declared
            # Every declaration is a floor -- never a pin below the environment,
            # and never a bare upper bound.
            assert not any(spec.operator in {"<", "<="} for spec in requirement.specifier)


def _dev_requirements_from_pyproject() -> dict[str, Requirement]:
    pyproject = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))
    return {
        canonicalize_name(Requirement(raw).name): Requirement(raw)
        for raw in pyproject["project"]["optional-dependencies"]["dev"]
    }


def test_no_formatter_declaration_sits_below_the_installed_toolchain():
    """The maintained environment is the reference; the codebase tracks it.

    ``AGENTS.md`` makes the ``datasci`` environment the standard for development,
    linting and packaging. A declaration that pins *below* what that environment
    carries inverts the relationship: it forces a downgrade of the reference
    instead of adapting the codebase to it. An exact ``==`` pin did exactly that
    for black and ruff, and the earlier fix for cycle-6 F-002 chose that
    direction; this test pins the correct one.

    Two separate things are asserted:

    * the ``[dev]`` specifier must admit the installed version, so
      ``pip install -e '.[dev]'`` never downgrades the reference environment; and
    * the ``.pre-commit-config.yaml`` revision must not be *older* than the
      installed version, so ``pre-commit install`` does not reformat the tree
      differently from ``make format-check`` -- which runs the environment's own
      black, exactly as ``make lint`` runs the environment's own ruff.

    A newer revision than the installed one is allowed: it means the environment
    is behind, which upgrading it fixes without touching the repository.
    """
    config = yaml.safe_load((REPO / ".pre-commit-config.yaml").read_text(encoding="utf-8"))
    revisions = {repo["repo"]: str(repo["rev"]) for repo in config["repos"]}
    dev_requirements = _dev_requirements_from_pyproject()

    for url, tool in PINNED_HOOK_TOOLS.items():
        assert url in revisions, f".pre-commit-config.yaml no longer declares {url}"
        installed = Version(metadata.version(tool))
        requirement = dev_requirements[tool]

        assert requirement.specifier.contains(installed, prereleases=True), (
            f"[dev] declares {tool}{requirement.specifier!s} but the maintained "
            f"environment has {tool} {installed}. The environment is the "
            f"reference: raise the declaration, do not downgrade the environment."
        )
        hook_version = Version(revisions[url].lstrip("v"))
        # `make format-check` runs the Makefile's PYTHON, so the revision has to
        # be compared against THAT environment's tool -- not the interpreter
        # running this test, which under the 3.11 matrix cell is a throwaway
        # venv carrying whatever version uv last resolved.
        reference = _maintained_tool_version(tool)
        if reference is None:
            pytest.skip(f"the maintained environment is not available to read {tool} from")
        assert hook_version >= reference, (
            f".pre-commit-config.yaml runs {tool} {hook_version} but the "
            f"maintained environment has {reference}; the commit hook would "
            f"format differently from `make format-check`. Bump the revision."
        )


def test_installed_runtime_dependencies_satisfy_declared_requirements():
    grouped = _requirements_by_extra()
    for group in ("core", "plot", "perf", "extras", "tree"):
        for name, requirement in grouped[group].items():
            installed = metadata.version(name)
            assert requirement.specifier.contains(
                installed
            ), f"installed {name}=={installed} does not satisfy {requirement.specifier}"


def _package_import_roots() -> dict[str, list[str]]:
    """Top-level module roots imported anywhere in the tracked package sources.

    Covers ``import x``, ``from x import ...`` and a literal
    ``importlib.import_module("x")`` / ``__import__("x")``.
    """
    import ast

    listed = subprocess.run(
        ["git", "-C", str(REPO), "ls-files", "traminex"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.split()
    sources = [name for name in listed if name.endswith(".py")]
    assert sources, "no tracked traminex sources found"

    roots: dict[str, list[str]] = {}
    for relative in sources:
        tree = ast.parse((REPO / relative).read_text(encoding="utf-8"), relative)
        for node in ast.walk(tree):
            found: list[str] = []
            if isinstance(node, ast.Import):
                found = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                found = [node.module]
            elif isinstance(node, ast.Call):
                function = node.func
                is_import_call = (
                    isinstance(function, ast.Attribute) and function.attr == "import_module"
                ) or (isinstance(function, ast.Name) and function.id in {"__import__"})
                if (
                    is_import_call
                    and node.args
                    and isinstance(node.args[0], ast.Constant)
                    and isinstance(node.args[0].value, str)
                ):
                    found = [node.args[0].value]
            for name in found:
                roots.setdefault(name.split(".")[0], []).append(f"{relative}:{node.lineno}")
    return roots


def test_every_optional_dependency_extra_has_a_real_import_site():
    """A mapped distribution the package never imports advertises nothing.

    ``__getattr__`` reads ``_OPTIONAL_DEPENDENCY_EXTRAS`` only from the module
    name inside a caught ``ImportError``. If traminex never imports that module,
    the handler can never see the name and the advertised "install the '<extra>'
    extra" message is unreachable. ``lifelines`` and ``graphviz`` were mapped in
    exactly that state (cycle-6 F-004); the install-profile verifier and the
    declared-floor test both certified them because they check distribution
    *presence*, never *use*.
    """
    import traminex

    roots = _package_import_roots()
    unreachable = {
        module: extra
        for module, extra in traminex._OPTIONAL_DEPENDENCY_EXTRAS.items()
        if not roots.get(module)
    }
    assert not unreachable, (
        "these _OPTIONAL_DEPENDENCY_EXTRAS keys have no import site in traminex/, "
        f"so their extras message can never fire: {sorted(unreachable)}"
    )


def test_optional_dependency_extras_name_declared_extras():
    """Every mapped extra must exist in [project.optional-dependencies]."""
    import traminex

    pyproject = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))
    declared = set(pyproject["project"]["optional-dependencies"])
    mapped = set(traminex._OPTIONAL_DEPENDENCY_EXTRAS.values())
    assert (
        mapped <= declared
    ), f"undeclared extras named in the mapping: {sorted(mapped - declared)}"


def test_core_import_is_lazy_with_respect_to_optional_dependencies():
    optional_roots = {
        "graphviz",
        "lifelines",
        "matplotlib",
        "numba",
        "seaborn",
        "sklearn",
        "traminex._nb_kernels",
        "traminex._pb_kernels",
    }
    result = _run_isolated(
        "import json, sys; import traminex; "
        f"names = {json.dumps(sorted(optional_roots))}; "
        "print(json.dumps(sorted(name for name in names if name in sys.modules)))"
    )

    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout) == []


def test_missing_optional_dependency_error_is_feature_local():
    result = _run_isolated(textwrap.dedent("""
            import importlib.abc
            import sys

            class BlockSklearn(importlib.abc.MetaPathFinder):
                def find_spec(self, fullname, path=None, target=None):
                    if fullname == "sklearn" or fullname.startswith("sklearn."):
                        raise ModuleNotFoundError("blocked optional dependency")
                    return None

            sys.meta_path.insert(0, BlockSklearn())
            import numpy as np
            import traminex

            assert "sklearn" not in sys.modules
            try:
                traminex.cluster_agglomerative(np.zeros((2, 2)), 1)
            except RuntimeError as exc:
                assert str(exc) == "scikit-learn is required for clustering"
            else:
                raise AssertionError("missing optional dependency did not fail locally")

            try:
                traminex.disstree_fit(np.zeros((2, 2)), [0, 1])
            except RuntimeError as exc:
                assert str(exc) == "scikit-learn is required for disstree_fit"
            else:
                raise AssertionError("missing disstree dependency did not fail locally")

            sequences = traminex.seqdef(np.array([["A", "B"], ["B", "A"]]))
            try:
                traminex.seqtree_fit(sequences, [0, 1])
            except RuntimeError as exc:
                assert str(exc) == "scikit-learn is required for seqtree_fit"
            else:
                raise AssertionError("missing seqtree dependency did not fail locally")
            """))

    assert result.returncode == 0, result.stderr


def test_exploratory_tree_failures_remain_local_to_requested_feature(monkeypatch):
    import traminex as tx
    from traminex import disstree_explore, seqtree_explore

    # scikit-learn is imported where it is used (TX-AUD-20261011-034): block the import.
    monkeypatch.setitem(sys.modules, "sklearn", None)
    with pytest.raises(RuntimeError, match="scikit-learn is required for disstree_fit"):
        disstree_explore.disstree_fit(np.zeros((2, 2)), [0, 1])

    # A real sequence object: the type of x is checked first (TX-AUD-20261019-036).
    sequences = tx.seqdef(np.array([["A", "B"], ["B", "A"]], dtype=object))
    with pytest.raises(RuntimeError, match="scikit-learn is required for seqtree_fit"):
        seqtree_explore.seqtree_fit(sequences, [0, 1])


def test_the_maintained_interpreter_is_the_one_release_check_runs(monkeypatch, tmp_path) -> None:
    """``make release-check`` passes its PYTHON on (TX-AUD-20261011-076).

    CI runs ``make release-check PYTHON=python``; the Makefile's default path does not
    exist there, so the formatter check skipped. A command name is resolved on PATH, as
    make's shell resolves it.
    """
    fake = tmp_path / "python"
    fake.write_text("#!/bin/sh\n", encoding="utf-8")
    fake.chmod(0o755)
    monkeypatch.setenv("PATH", str(tmp_path))
    monkeypatch.setenv("TRAMINEX_MAKE_PYTHON", "python")

    assert _maintained_python() == fake

    monkeypatch.setenv("TRAMINEX_MAKE_PYTHON", str(tmp_path / "absent-python"))
    assert _maintained_python() is None


def test_release_check_passes_its_python_to_the_formatter_check() -> None:
    makefile = (REPO / "Makefile").read_text(encoding="utf-8")

    assert (
        'TRAMINEX_MAKE_PYTHON="$(PYTHON)" $(PYTHON) -m pytest tests/test_supported_versions.py'
        in (makefile)
    )
