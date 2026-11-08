from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CHECKER = ROOT / "scripts" / "check_stub_signatures.py"
SYNC = ROOT / "scripts" / "sync_stub_public_api.py"


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _fixture_package(tmp_path: Path) -> tuple[Path, Path, Path]:
    runtime = tmp_path / "fixturepkg"
    stubs = tmp_path / "stubs" / "fixturepkg"
    manifest = runtime / "api_manifest.py"
    _write(
        runtime / "__init__.py",
        'from .public import Widget, public\n\n__all__ = ["public", "Widget"]\n',
    )
    _write(manifest, 'PUBLIC_API = ["public", "Widget"]\n')
    _write(
        runtime / "public.py",
        "def public(a, /, b=1, *, flag, option=None, **kwargs):\n"
        "    return a, b, flag, option, kwargs\n\n"
        "class Widget:\n"
        "    pass\n",
    )
    _write(
        stubs / "__init__.pyi",
        "from .public import Widget as Widget\n"
        "from .public import public as public\n\n"
        '__all__ = ["public", "Widget"]\n',
    )
    _write(stubs / "api_manifest.pyi", "PUBLIC_API: list[str]\n")
    _write(
        stubs / "public.pyi",
        "def public(a, /, b=..., *, flag, option=..., **kwargs): ...\n\n" "class Widget: ...\n",
    )
    return runtime, stubs, manifest


def _run_checker(runtime: Path, stubs: Path) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join([str(runtime.parent), env.get("PYTHONPATH", "")]).rstrip(
        os.pathsep
    )
    return subprocess.run(
        [
            sys.executable,
            str(CHECKER),
            "--runtime-root",
            str(runtime),
            "--stub-root",
            str(stubs),
            "--package",
            "fixturepkg",
        ],
        cwd=runtime.parent,
        env=env,
        check=False,
        capture_output=True,
        text=True,
    )


def test_stub_checker_accepts_clean_bidirectional_package(tmp_path: Path):
    runtime, stubs, _manifest = _fixture_package(tmp_path)

    result = _run_checker(runtime, stubs)

    assert result.returncode == 0, result.stdout + result.stderr
    assert "Runtime modules: 3" in result.stdout
    assert "Stub modules: 3" in result.stdout
    assert "Public functions: 2" in result.stdout
    assert "Public classes: 2" in result.stdout


def test_stub_checker_detects_missing_module(tmp_path: Path):
    runtime, stubs, _manifest = _fixture_package(tmp_path)
    (stubs / "public.pyi").unlink()

    result = _run_checker(runtime, stubs)

    assert result.returncode == 1
    assert "missing stub module fixturepkg.public" in result.stdout


def test_stub_checker_detects_orphan_stub(tmp_path: Path):
    runtime, stubs, _manifest = _fixture_package(tmp_path)
    _write(stubs / "orphan.pyi", "def ghost(): ...\n")

    result = _run_checker(runtime, stubs)

    assert result.returncode == 1
    assert "orphan stub module fixturepkg.orphan" in result.stdout


def test_stub_checker_detects_missing_and_extra_symbols(tmp_path: Path):
    runtime, stubs, _manifest = _fixture_package(tmp_path)
    _write(stubs / "public.pyi", "class Widget: ...\n\ndef ghost(): ...\n")

    result = _run_checker(runtime, stubs)

    assert result.returncode == 1
    assert "stub missing public function 'public'" in result.stdout
    assert "runtime module fixturepkg.public missing stub function 'ghost'" in result.stdout


def test_stub_checker_detects_a_reexport_the_runtime_module_lacks(tmp_path: Path):
    # A stub without __all__ re-exports ("from x import y as y") a name its runtime module does not
    # have (stubs/traminex/seqtree_explore.pyi kept seqtransn_count after the module stopped
    # importing it; TX-REP-20261023-006): the name comparison reads definitions only, so it passed.
    runtime, stubs, _manifest = _fixture_package(tmp_path)
    _write(runtime / "other.py", "def other():\n    return 1\n")
    _write(
        stubs / "other.pyi", "from fixturepkg.public import Widget as Widget\n\ndef other(): ...\n"
    )

    result = _run_checker(runtime, stubs)

    assert result.returncode == 1
    assert (
        "fixturepkg.other: 'Widget' re-exports fixturepkg.public.Widget, but the runtime module "
        "has no 'Widget'" in result.stdout
    )


def test_stub_checker_detects_parameter_kind_and_default_drift(tmp_path: Path):
    runtime, stubs, _manifest = _fixture_package(tmp_path)
    _write(
        stubs / "public.pyi",
        "def public(a, b, *, flag=..., option=..., **kwargs): ...\n\n" "class Widget: ...\n",
    )

    result = _run_checker(runtime, stubs)

    assert result.returncode == 1
    assert "signature mismatch" in result.stdout
    assert "POSITIONAL_ONLY" in result.stdout
    assert "required" in result.stdout


_CLASS_RUNTIME = (
    "from dataclasses import dataclass, field\n\n"
    "@dataclass\n"
    "class Config:\n"
    "    metric: str = 'l1'\n\n"
    "@dataclass\n"
    "class Result:\n"
    "    values: list\n"
    "    labels: list | None = None\n"
    "    max_time: int = 0\n\n"
    "class Holder:\n"
    "    def __init__(self, values, *, meta):\n"
    "        self.values = values\n"
    "        self.meta = meta\n\n"
    "class Plain:\n"
    "    pass\n\n"
    "@dataclass(kw_only=True)\n"
    "class Options:\n"
    "    width: int\n"
    "    height: int = 1\n"
)

_CLASS_STUB_CLEAN = (
    "from dataclasses import dataclass\n\n"
    "@dataclass\n"
    "class Config:\n"
    "    metric: str = 'l1'\n\n"
    "@dataclass\n"
    "class Result:\n"
    "    values: list\n"
    "    labels: list | None = ...\n"
    "    max_time: int = 0\n\n"
    "class Holder:\n"
    "    values: list\n"
    "    def __init__(self, values, *, meta) -> None: ...\n\n"
    "class Plain: ...\n\n"
    "@dataclass(kw_only=True)\n"
    "class Options:\n"
    "    width: int\n"
    "    height: int = ...\n"
)


def _class_fixture(tmp_path: Path, stub_text: str) -> tuple[Path, Path]:
    """Fixture package whose ``classes`` module exercises constructor contracts.

    Audit S7-02 (2026-09-22): the checker compared only top-level function
    signatures and class presence, so ``DURConfig.eps`` (a stub field the runtime
    dataclass had dropped) and nine result classes whose stubs declared no
    constructor passed the gate.
    """
    runtime, stubs, _manifest = _fixture_package(tmp_path)
    _write(runtime / "classes.py", _CLASS_RUNTIME)
    _write(stubs / "classes.pyi", stub_text)
    return runtime, stubs


def test_stub_checker_accepts_matching_class_constructors(tmp_path: Path):
    runtime, stubs = _class_fixture(tmp_path, _CLASS_STUB_CLEAN)

    result = _run_checker(runtime, stubs)

    assert result.returncode == 0, result.stdout + result.stderr
    # Config, Result, Holder and Options; Plain declares no constructor on
    # either side.
    assert "Constructors checked: 4" in result.stdout


def test_stub_checker_detects_dataclass_keyword_only_drift(tmp_path: Path):
    stub = _CLASS_STUB_CLEAN.replace("@dataclass(kw_only=True)\n", "@dataclass\n")
    runtime, stubs = _class_fixture(tmp_path, stub)

    result = _run_checker(runtime, stubs)

    assert result.returncode == 1
    assert "fixturepkg.classes.Options constructor mismatch" in result.stdout


def test_stub_checker_detects_extra_dataclass_field(tmp_path: Path):
    stub = _CLASS_STUB_CLEAN.replace(
        "    metric: str = 'l1'\n", "    metric: str = 'l1'\n    eps: float = ...\n"
    )
    runtime, stubs = _class_fixture(tmp_path, stub)

    result = _run_checker(runtime, stubs)

    assert result.returncode == 1
    assert "fixturepkg.classes.Config constructor mismatch" in result.stdout
    assert "'eps'" in result.stdout


def test_stub_checker_detects_undeclared_constructor(tmp_path: Path):
    stub = _CLASS_STUB_CLEAN.replace(
        "    def __init__(self, values, *, meta) -> None: ...\n", "    meta: dict\n"
    )
    runtime, stubs = _class_fixture(tmp_path, stub)

    result = _run_checker(runtime, stubs)

    assert result.returncode == 1
    assert (
        "fixturepkg.classes.Holder: stub declares no constructor, but the runtime "
        "class defines __init__(values, *, meta)" in result.stdout
    )


def test_stub_checker_detects_dataclass_stub_without_decorator(tmp_path: Path):
    stub = _CLASS_STUB_CLEAN.replace("@dataclass\nclass Result:\n", "class Result:\n")
    runtime, stubs = _class_fixture(tmp_path, stub)

    result = _run_checker(runtime, stubs)

    assert result.returncode == 1
    assert "fixturepkg.classes.Result: stub declares no constructor" in result.stdout


def test_stub_checker_detects_missing_and_phantom_dataclass_fields(tmp_path: Path):
    # An explicit __init__ that matches the runtime does not excuse a missing
    # attribute declaration (the SeqSurv.max_time case), and a declared
    # attribute the runtime dataclass does not have is reported (SeqImplic).
    stub = _CLASS_STUB_CLEAN.replace(
        "@dataclass\nclass Result:\n"
        "    values: list\n"
        "    labels: list | None = ...\n"
        "    max_time: int = 0\n",
        "class Result:\n"
        "    values: list\n"
        "    labels: list | None\n"
        "    duration: int\n"
        "    def __init__(self, values, labels=..., max_time=...) -> None: ...\n",
    )
    runtime, stubs = _class_fixture(tmp_path, stub)

    result = _run_checker(runtime, stubs)

    assert result.returncode == 1
    assert "fixturepkg.classes.Result: stub does not declare field 'max_time'" in result.stdout
    assert (
        "fixturepkg.classes.Result: stub declares attribute 'duration' that the "
        "runtime dataclass does not have" in result.stdout
    )
    assert "constructor mismatch" not in result.stdout


def test_stub_checker_does_not_require_private_fields_outside_the_constructor(tmp_path: Path):
    # A private ``field(init=False)`` is implementation state, not public API:
    # Eseq._insertion_order (audit S5-01) made the real-tree check fail once the
    # S7-02 field rule met it. A missing public field is still reported (above).
    runtime, stubs = _class_fixture(tmp_path, _CLASS_STUB_CLEAN)
    _write(
        runtime / "classes.py",
        _CLASS_RUNTIME.replace(
            "    max_time: int = 0\n\n",
            "    max_time: int = 0\n"
            "    _order: list | None = field(default=None, init=False, repr=False)\n\n",
        ),
    )

    result = _run_checker(runtime, stubs)

    assert result.returncode == 0, result.stdout + result.stderr
    assert "_order" not in result.stdout


def test_stub_checker_detects_dataclass_default_value_drift(tmp_path: Path):
    stub = _CLASS_STUB_CLEAN.replace("    metric: str = 'l1'\n", "    metric: str = 'l2'\n")
    runtime, stubs = _class_fixture(tmp_path, stub)

    result = _run_checker(runtime, stubs)

    assert result.returncode == 1
    assert (
        "fixturepkg.classes.Config.__init__ default mismatch for 'metric': "
        "stub='l2' runtime='l1'" in result.stdout
    )


def test_real_stub_tree_passes_the_contract_checker():
    """The maintained ``stubs/traminex`` tree satisfies every contract above."""
    result = subprocess.run(
        [sys.executable, str(CHECKER)],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert "Stub contract check PASSED." in result.stdout


_REDIRECTED_RUN = """
import runpy
import sys
from importlib.abc import MetaPathFinder
from importlib.machinery import PathFinder


class {finder}(MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == "fixturepkg" or fullname.startswith("fixturepkg."):
            search = [{shadow_parent!r}] if fullname == "fixturepkg" else path
            return PathFinder.find_spec(fullname, search)
        return None


sys.meta_path.insert(0, {finder}())
sys.path.insert(0, {scripts!r})
sys.argv = [{checker!r}, "--runtime-root", {runtime!r}, "--stub-root", {stubs!r},
            "--package", "fixturepkg"]
runpy.run_path({checker!r}, run_name="__main__")
"""


def _run_checker_behind_a_redirect(tmp_path: Path, finder: str) -> subprocess.CompletedProcess[str]:
    """Check a drifted copy while a meta-path finder serves a clean shadow copy.

    Audit S7-12 (2026-09-22): the checker parsed ``--runtime-root`` from disk but
    imported the package through ``importlib``, where an editable-install finder
    outranks ``sys.path``; it then compared the copy's source with another tree's
    objects. Here the checked copy drifts from its stubs while the shadow copy
    the finder serves matches them, so a mixed run passes silently.
    """
    runtime, stubs, _manifest = _fixture_package(tmp_path)
    shadow = tmp_path / "shadow" / "fixturepkg"
    shutil.copytree(runtime, shadow)
    _write(
        runtime / "public.py",
        "def public(a, /, b=1, *, flag, option=None, extra=None, **kwargs):\n"
        "    return a, b, flag, option, extra, kwargs\n\n"
        "class Widget:\n"
        "    pass\n",
    )
    code = _REDIRECTED_RUN.format(
        finder=finder,
        shadow_parent=str(shadow.parent),
        scripts=str(CHECKER.parent),
        checker=str(CHECKER),
        runtime=str(runtime),
        stubs=str(stubs),
    )
    return subprocess.run(
        [sys.executable, "-c", code],
        cwd=tmp_path,
        check=False,
        capture_output=True,
        text=True,
    )


def test_stub_checker_imports_the_runtime_root_despite_an_editable_redirect(tmp_path: Path):
    result = _run_checker_behind_a_redirect(tmp_path, "_EditableRedirectFinder")

    assert result.returncode == 1, result.stdout + result.stderr
    assert "fixturepkg.public.public signature mismatch" in result.stdout
    assert "'extra'" in result.stdout


def test_stub_checker_refuses_a_runtime_imported_from_another_tree(tmp_path: Path):
    result = _run_checker_behind_a_redirect(tmp_path, "_OpaqueFinder")

    assert result.returncode == 1, result.stdout + result.stderr
    assert "outside --runtime-root" in result.stdout


def test_stub_sync_check_is_non_mutating_and_update_is_deterministic(tmp_path: Path):
    _runtime, stubs, manifest = _fixture_package(tmp_path)
    stub_init = stubs / "__init__.pyi"
    original = stub_init.read_text(encoding="utf-8").replace(
        '__all__ = ["public", "Widget"]', '__all__ = ["Widget", "public"]'
    )
    stub_init.write_text(original, encoding="utf-8")
    base_args = [
        sys.executable,
        str(SYNC),
        "--manifest",
        str(manifest),
        "--stub-init",
        str(stub_init),
    ]

    checked = subprocess.run([*base_args, "--check"], check=False, capture_output=True, text=True)
    assert checked.returncode == 1
    assert stub_init.read_text(encoding="utf-8") == original

    updated = subprocess.run(base_args, check=False, capture_output=True, text=True)
    assert updated.returncode == 0, updated.stdout + updated.stderr
    first = stub_init.read_bytes()
    repeated = subprocess.run(base_args, check=False, capture_output=True, text=True)
    assert repeated.returncode == 0, repeated.stdout + repeated.stderr
    assert stub_init.read_bytes() == first
    final_check = subprocess.run(
        [*base_args, "--check"], check=False, capture_output=True, text=True
    )
    assert final_check.returncode == 0, final_check.stdout + final_check.stderr
