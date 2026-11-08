from __future__ import annotations

import os
import subprocess
import sys
from importlib import import_module
from importlib.abc import MetaPathFinder
from importlib.machinery import PathFinder
from importlib.util import spec_from_file_location
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
ACTIVE_UBUNTU_DOCS = (
    ROOT / "README.md",
    ROOT / "CONTRIBUTING.md",
    ROOT / "docs" / "DEV_BUILD.md",
    ROOT / "docs" / "INDEX.md",
    ROOT / "stubs" / "README.md",
)


def test_no_powershell_development_helpers_remain():
    assert list((ROOT / "scripts").glob("*.ps1")) == []


def test_pybind11_test_gate_is_directly_executable():
    launcher = ROOT / "scripts" / "run_pybind11_tests.sh"
    assert os.access(launcher, os.X_OK)


def test_checkout_import_preference_defeats_an_editable_redirect(tmp_path, monkeypatch):
    from scripts._checkout_imports import prefer_checkout_imports

    package_name = "checkout_import_probe"
    current_package = tmp_path / "current" / package_name
    redirected_package = tmp_path / "redirected" / package_name
    current_package.mkdir(parents=True)
    redirected_package.mkdir(parents=True)
    (current_package / "__init__.py").write_text("origin = 'current'\n", encoding="utf-8")
    (redirected_package / "__init__.py").write_text("origin = 'redirected'\n", encoding="utf-8")

    class EditableRedirect(MetaPathFinder):
        def find_spec(self, fullname, path=None, target=None):
            if fullname == package_name:
                return spec_from_file_location(fullname, redirected_package / "__init__.py")
            return None

    original_meta_path = sys.meta_path.copy()
    monkeypatch.syspath_prepend(str(current_package.parent))
    monkeypatch.delitem(sys.modules, package_name, raising=False)
    sys.meta_path.insert(0, EditableRedirect())
    try:
        prefer_checkout_imports()
        imported = import_module(package_name)
    finally:
        sys.meta_path[:] = original_meta_path
        sys.modules.pop(package_name, None)

    assert imported.origin == "current"
    assert sys.meta_path.index(PathFinder) == original_meta_path.index(PathFinder)


def test_checkout_import_preference_preserves_builtin_module_priority(tmp_path, monkeypatch):
    from scripts._checkout_imports import prefer_checkout_imports

    module_name = "_symtable"
    assert module_name in sys.builtin_module_names
    (tmp_path / f"{module_name}.py").write_text(
        "origin = 'path-shadow'\n",
        encoding="utf-8",
    )

    original_meta_path = sys.meta_path.copy()
    monkeypatch.syspath_prepend(str(tmp_path))
    monkeypatch.delitem(sys.modules, module_name, raising=False)
    try:
        prefer_checkout_imports()
        imported = import_module(module_name)
    finally:
        sys.meta_path[:] = original_meta_path
        sys.modules.pop(module_name, None)

    assert imported.__spec__.origin == "built-in"
    assert not hasattr(imported, "origin")


def test_make_typecheck_uses_the_canonical_bare_mypy_command():
    makefile = (ROOT / "Makefile").read_text(encoding="utf-8")
    typecheck = makefile.split("typecheck:\n", maxsplit=1)[1].split("\n\n", maxsplit=1)[0]

    assert typecheck.strip() == "$(PYTHON) -m mypy"
    assert "follow-imports" not in makefile


def test_format_check_runs_the_environment_black_in_read_only_mode():
    """``format-check`` must check with the same black the developer runs.

    It used to invoke pre-commit's own isolated, pinned black. That made the gate
    structurally unable to see a disagreement with the maintained environment --
    and there was one: the gate passed on black 24.8.0 while the environment's
    26.5.1 reformatted eight files, so the workflow ``CONTRIBUTING.md`` documents
    broke the gate. Running the environment's black, exactly as ``lint`` runs the
    environment's ruff, makes that divergence unrepresentable rather than merely
    detected.

    The pre-commit hooks stay for ``pre-commit install`` users;
    ``tests/test_supported_versions.py`` keeps their revisions from falling
    behind the environment.
    """
    makefile = (ROOT / "Makefile").read_text(encoding="utf-8")
    release_target = makefile.split("release-check:", maxsplit=1)[1].splitlines()[0]
    format_target = makefile.split("format-check:\n", maxsplit=1)[1].split("\n\n", maxsplit=1)[0]
    lint_target = makefile.split("lint:\n", maxsplit=1)[1].split("\n\n", maxsplit=1)[0]

    assert "format-check" in release_target.split()

    normalized = " ".join(format_target.replace("\\\n", " ").split())
    # black's notebook support needs IPython and tokenize_rt; the prerequisite check
    # names them (and the install command) before black runs (TX-AUD-20261010-046).
    assert normalized.startswith("$(PYTHON) scripts/check_dev_environment.py"), normalized
    sys.path.insert(0, str(ROOT / "scripts"))
    import check_dev_environment

    assert {"IPython", "tokenize_rt", "black"} <= set(check_dev_environment.REQUIRED)
    assert "$(PYTHON) -m black --check" in normalized, normalized
    assert "pre_commit" not in normalized, normalized
    # Read-only: the gate must never rewrite the tree.
    assert " -w" not in normalized and "--write" not in normalized

    # Both formatter gates run the environment's own tool, from the same Python.
    assert lint_target.strip().startswith("$(PYTHON) -m ruff check")

    # Every tracked Python file is covered by the paths the gate passes.
    listed = subprocess.run(
        ["git", "-C", str(ROOT), "ls-files", "*.py"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.split()
    scope = tuple(
        token for token in normalized.split() if not token.startswith(("$", "-", "black", "python"))
    )
    uncovered = [
        name
        for name in listed
        if not any(name == s or name.startswith(s.rstrip("/") + "/") for s in scope)
    ]
    assert not uncovered, f"format-check does not cover: {uncovered[:10]}"

    config = yaml.safe_load((ROOT / ".pre-commit-config.yaml").read_text(encoding="utf-8"))
    black_repo = next(repo for repo in config["repos"] if repo["repo"].endswith("/psf/black"))
    black_check = next(hook for hook in black_repo["hooks"] if hook.get("alias") == "black-check")

    assert black_check["id"] == "black"
    assert "--check" in black_check["args"]
    assert black_check["stages"] == ["manual"]


def test_active_build_docs_use_only_supported_ubuntu_commands():
    text = "\n".join(path.read_text(encoding="utf-8") for path in ACTIVE_UBUNTU_DOCS)

    for forbidden in (".ps1", "```powershell", "follow-imports", "TRAMINEX_BUILD_NANOBIND"):
        assert forbidden not in text
    assert "Python 3.11" in text
    assert "Python 3.12" in text
    assert "make typecheck" in text
    assert "make test-pybind11" in text


def test_lint_covers_every_tracked_python_file():
    """``make lint`` must cover what ``make format-check`` covers (S7-22).

    ``lint`` passed ``traminex tests scripts benchmarks`` while ``format-check``
    also passed ``examples``, ``docs`` and ``verify_audit_improvements.py``, so
    ruff never saw the example scripts or the release verifier (audit S7-22,
    2026-09-22).
    """
    makefile = (ROOT / "Makefile").read_text(encoding="utf-8")
    lint_target = makefile.split("lint:\n", maxsplit=1)[1].split("\n\n", maxsplit=1)[0]
    normalized = " ".join(lint_target.replace("\\\n", " ").split())
    assert normalized.startswith("$(PYTHON) -m ruff check"), normalized
    scope = tuple(
        token
        for token in normalized.split()
        if not token.startswith(("$", "-")) and token not in {"ruff", "check"}
    )

    listed = subprocess.run(
        ["git", "-C", str(ROOT), "ls-files", "*.py"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.split()
    uncovered = [
        name
        for name in listed
        if not any(name == s or name.startswith(s.rstrip("/") + "/") for s in scope)
    ]
    assert not uncovered, f"lint does not cover: {uncovered[:10]}"


def _run_native_gate(tmp_path: Path, body: str) -> subprocess.CompletedProcess:
    """Run pytest with the native gate's skip plugin on one temporary test module."""
    (tmp_path / "test_probe.py").write_text(
        "import pytest\n\n\ndef test_passes():\n    assert True\n\n\n" + body, encoding="utf-8"
    )
    program = (
        "import sys, pytest\n"
        "from scripts._native_skip_gate import NativeSkipGate\n"
        "sys.exit(pytest.main(sys.argv[1:], plugins=[NativeSkipGate()]))\n"
    )
    return subprocess.run(
        [sys.executable, "-c", program, str(tmp_path), "-q", "-p", "no:cacheprovider",
         "--rootdir", str(tmp_path)],  # fmt: skip
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )


def test_native_runner_fails_on_a_compiled_skip(tmp_path) -> None:
    """docs/DEV_BUILD.md: a skipped compiled differential module is a failed pybind11 gate
    even when pytest exits zero (TX-AUD-20261011-075, TX-AUD-20261019-077)."""
    result = _run_native_gate(
        tmp_path,
        "def test_kernel():\n    pytest.skip('compiled _pb_kernels not built')\n",
    )

    assert result.returncode == 1, result.stdout + result.stderr
    assert "compiled _pb_kernels not built" in result.stdout


def test_native_runner_accepts_an_environment_skip(tmp_path) -> None:
    result = _run_native_gate(
        tmp_path, "def test_dot():\n    pytest.skip('Graphviz dot is unavailable')\n"
    )

    assert result.returncode == 0, result.stdout + result.stderr


def test_the_native_runner_loads_the_skip_gate() -> None:
    launcher = (ROOT / "scripts" / "run_pybind11_tests.sh").read_text(encoding="utf-8")

    assert "plugins=[NativeSkipGate()]" in launcher


def test_every_compiled_skip_reason_is_known_to_the_native_gate() -> None:
    """A skip reason in tests/ naming the compiled kernels is one the gate refuses, and the
    gate lists no reason the tests no longer use."""
    import ast

    from scripts._native_skip_gate import COMPILED_SKIP_REASONS

    reasons = set()
    for path in sorted((ROOT / "tests").glob("test_*.py")):
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if not isinstance(node, ast.Call):
                continue
            function = node.func
            name = function.attr if isinstance(function, ast.Attribute) else None
            if name not in ("skip", "skipif", "importorskip"):
                continue
            given = [keyword.value for keyword in node.keywords if keyword.arg == "reason"]
            if name == "skip" and node.args:
                given.append(node.args[0])
            reasons.update(
                value.value
                for value in given
                if isinstance(value, ast.Constant) and isinstance(value.value, str)
            )
    skip_reasons = {
        text
        for text in reasons
        if any(word in text for word in ("_pb_kernels", "pybind11", "compiled kernel"))
    }

    assert skip_reasons == set(COMPILED_SKIP_REASONS)


# TX-AUD-20261010-088: a run killed before its EXIT trap leaves its link in traminex/, which
# every later run refused. The link is live while the run's build directory exists (until /tmp
# is cleared) and dead after. The runner records its process ID in its build directory and,
# before building, removes a link whose target is gone, or whose run's process is gone and
# whose module no process has loaded (with that build directory). Anything else is left alone
# and refused (batch-17 review I-1).
_NATIVE_LINK = ROOT / "scripts" / "_native_link.sh"
_MODULE = "_pb_kernels.cpython-312-x86_64-linux-gnu.so"


def _run_link_function(package: Path, call: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["bash", "-c", f'set -euo pipefail; source "$1"; {call}', "bash", str(_NATIVE_LINK)],
        cwd=package.parent,
        env={**os.environ, "PACKAGE": str(package)},
        capture_output=True,
        text=True,
        check=False,
    )


def _no_such_process() -> int:
    """A process ID no process can have: PIDs are below the kernel's pid_max."""
    return int(Path("/proc/sys/kernel/pid_max").read_text(encoding="utf-8"))


def _package_with(tmp_path: Path, kind: str, runner_build: Path | None = None) -> Path:
    package = tmp_path / "traminex"
    package.mkdir()
    module = package / _MODULE
    if kind in {"killed run link", "running run link", "runner link without a process ID"}:
        assert runner_build is not None
        (runner_build / _MODULE).write_bytes(b"built")
        module.symlink_to(runner_build / _MODULE)
        if kind != "runner link without a process ID":
            pid = _no_such_process() if kind == "killed run link" else os.getpid()
            (runner_build / "runner.pid").write_text(f"{pid}\n", encoding="utf-8")
    elif kind == "dead runner link":
        module.symlink_to(f"/tmp/traminex-r1-pybind11.no-such-build/{_MODULE}")
    elif kind == "dead foreign link":
        module.symlink_to(tmp_path / "elsewhere" / _MODULE)
    elif kind == "live link":
        (tmp_path / "build").mkdir()
        (tmp_path / "build" / _MODULE).write_bytes(b"built")
        module.symlink_to(tmp_path / "build" / _MODULE)
    else:
        module.write_bytes(b"copied")
    return package


def test_the_native_runner_removes_its_own_dead_link(tmp_path: Path) -> None:
    package = _package_with(tmp_path, "dead runner link")

    result = _run_link_function(package, 'clear_dead_native_links "$PACKAGE"')

    assert result.returncode == 0, result.stderr
    assert not (package / _MODULE).is_symlink()
    assert "dead link" in result.stderr
    refused = _run_link_function(package, f'refuse_existing_module "$PACKAGE/{_MODULE}"')
    assert refused.returncode == 0, refused.stderr


@pytest.fixture
def runner_build():
    import shutil
    import tempfile

    directory = Path(tempfile.mkdtemp(prefix="traminex-r1-pybind11.", dir="/tmp"))
    yield directory
    shutil.rmtree(directory, ignore_errors=True)


def test_the_native_runner_removes_a_killed_runs_link_and_build_directory(
    tmp_path: Path, runner_build: Path
) -> None:
    package = _package_with(tmp_path, "killed run link", runner_build)

    result = _run_link_function(package, 'clear_dead_native_links "$PACKAGE"')

    assert result.returncode == 0, result.stderr
    assert not (package / _MODULE).is_symlink()
    assert not runner_build.exists()
    assert "killed run" in result.stderr


def test_the_native_runner_keeps_a_link_whose_module_a_process_has_loaded(
    tmp_path: Path, runner_build: Path
) -> None:
    import mmap

    package = _package_with(tmp_path, "killed run link", runner_build)
    with (
        open(runner_build / _MODULE, "rb") as handle,
        mmap.mmap(handle.fileno(), 0, access=mmap.ACCESS_READ),
    ):
        cleared = _run_link_function(package, 'clear_dead_native_links "$PACKAGE"')

    assert cleared.returncode == 0, cleared.stderr
    assert (package / _MODULE).is_symlink()
    assert (runner_build / "runner.pid").is_file()


@pytest.mark.parametrize(
    "kind",
    [
        "running run link",
        "runner link without a process ID",
        "dead foreign link",
        "live link",
        "regular file",
    ],
)
def test_the_native_runner_keeps_and_refuses_any_other_module(
    tmp_path: Path, runner_build: Path, kind: str
) -> None:
    package = _package_with(tmp_path, kind, runner_build)
    module = package / _MODULE
    before = (module.is_symlink(), os.readlink(module) if module.is_symlink() else None)

    cleared = _run_link_function(package, 'clear_dead_native_links "$PACKAGE"')
    refused = _run_link_function(package, f'refuse_existing_module "$PACKAGE/{_MODULE}"')

    assert cleared.returncode == 0, cleared.stderr
    assert (module.is_symlink(), os.readlink(module) if module.is_symlink() else None) == before
    assert module.is_symlink() or module.read_bytes() == b"copied"
    assert refused.returncode == 1
    assert "Refusing to replace existing module" in refused.stderr


def test_the_native_runner_uses_the_link_functions() -> None:
    launcher = (ROOT / "scripts" / "run_pybind11_tests.sh").read_text(encoding="utf-8")

    assert 'source "$ROOT/scripts/_native_link.sh"' in launcher
    assert launcher.index('clear_dead_native_links "$ROOT/traminex"') < launcher.index("cmake -S")
    assert 'refuse_existing_module "$PB_CANDIDATE"' in launcher
    record = 'echo "$$" > "$BUILD_DIR/$NATIVE_PID_FILE"'
    assert (
        launcher.index("BUILD_DIR=$(mktemp") < launcher.index(record) < launcher.index("cmake -S")
    )
