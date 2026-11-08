from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tarfile
import tomllib
import zipfile
from pathlib import Path

from scripts import sync_sdist_manifest

ROOT = Path(__file__).resolve().parents[1]


def test_sdist_policy_is_an_exact_existing_source_manifest():
    settings = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["tool"][
        "scikit-build"
    ]
    declared = settings["sdist"]["include"]

    # No inclusion-mode: scikit-build-core 0.11.0, the declared floor, does not
    # know it (TX-AUD-20261010-047); an entry anchored with "/" is a root path.
    assert "inclusion-mode" not in settings["sdist"]
    assert settings["sdist"]["exclude"] == ["**"]
    assert declared == sorted(set(declared), key=str.encode)
    assert not [path for path in declared if any(char in path for char in "*?[")]
    assert all((ROOT / path.lstrip("/")).is_file() for path in declared)
    assert "traminex/_pb_kernels.cp312-win_amd64.pyd" not in declared
    assert "AUDIT_ONLY_PROMPT.md" not in declared
    assert "FULL_CODEBASE_REMEDIATION_PROMPT.md" not in declared

    # The list is generated from `git ls-files` and the policy in
    # scripts/sync_sdist_manifest.py (`make sdist-manifest`; TX-AUD-20261019-078).
    assert declared == sync_sdist_manifest.sdist_entries(ROOT)


def _build_sdist(
    source: Path, output_dir: Path
) -> tuple[subprocess.CompletedProcess[str], set[str]]:
    output_dir.mkdir()
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "from scikit_build_core.build import build_sdist; "
            "print(build_sdist(__import__('sys').argv[1]))",
            str(output_dir),
        ],
        cwd=source,
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        return result, set()
    sdist = next(output_dir.glob("traminex-0.2.2.tar.gz"))
    with tarfile.open(sdist, "r:gz") as archive:
        members = {
            member.name.removeprefix("traminex-0.2.2/")
            for member in archive.getmembers()
            if member.isfile()
        }
    return result, members


def _copy_tracked_source(destination: Path) -> None:
    destination.mkdir()
    tracked = subprocess.run(
        ["git", "-C", str(ROOT), "ls-files", "-z"],
        check=True,
        capture_output=True,
    ).stdout.split(b"\0")
    source_paths = {encoded.decode() for encoded in tracked if encoded}
    for source_path in sorted(source_paths, key=str.encode):
        relative = Path(source_path)
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / relative, target)


def test_sdist_excludes_internal_control_and_reference_artifacts(tmp_path: Path):
    source = tmp_path / "source"
    source.mkdir()
    for name in ("CMakeLists.txt", "pyproject.toml", "README.md"):
        shutil.copy2(ROOT / name, source / name)
    package = source / "traminex"
    package.mkdir()
    shutil.copy2(ROOT / "traminex" / "__init__.py", package / "__init__.py")
    shutil.copy2(ROOT / "traminex" / "py.typed", package / "py.typed")
    internal_paths = [
        ".coverage",
        ".coverage.worker",
        "docs/superpowers/prompts/control.md",
        "docs/superpowers/plans/internal-plan.md",
        "docs/superpowers/research/internal-notes.md",
        "graphify-out/graph.json",
        ".codebase-memory/index.json",
        ".planning/state.json",
        ".external/Sequenzo/reference.py",
        "AUDIT_AND_FIX_PROMPT.md",
        "reports/audit/full_codebase_audit_20990101.md",
        "reports/audit/full_codebase_audit_20990101_inventory.csv",
        "vendor/Sequenzo/reference.py",
    ]
    for relative_path in internal_paths:
        artifact = source / relative_path
        artifact.parent.mkdir(parents=True, exist_ok=True)
        artifact.write_text("internal project control or reference material\n", encoding="utf-8")
    output_dir = tmp_path / "sdist"
    output_dir.mkdir()
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "from scikit_build_core.build import build_sdist; "
            "print(build_sdist(__import__('sys').argv[1]))",
            str(output_dir),
        ],
        cwd=source,
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    sdist = next(output_dir.glob("traminex-0.2.2.tar.gz"))
    with tarfile.open(sdist, "r:gz") as archive:
        names = archive.getnames()
    members = {name.removeprefix("traminex-0.2.2/") for name in names}
    assert members.isdisjoint(internal_paths)
    assert any(name.endswith("/pyproject.toml") for name in names)
    assert any(name.endswith("/traminex/py.typed") for name in names)


def test_sdist_membership_is_unchanged_by_arbitrary_untracked_files(tmp_path: Path):
    clean_source = tmp_path / "clean-source"
    dirty_source = tmp_path / "dirty-source"
    _copy_tracked_source(clean_source)
    shutil.copytree(clean_source, dirty_source)
    sentinels = {
        "AUDIT_ONLY_PROMPT.md",
        "FULL_CODEBASE_REMEDIATION_PROMPT.md",
        "surprise-root.txt",
        "scratch/random/nested_payload.py",
        "docs/generated/random-note.txt",
        "traminex/arbitrary_untracked_module.py",
        "traminex/arbitrary_untracked_binary.pyd",
        "tests/r_goldens/arbitrary_untracked_fixture.csv",
    }
    for relative_name in sentinels:
        path = dirty_source / relative_name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("must not be distributed\n", encoding="utf-8")

    clean_result, clean_members = _build_sdist(clean_source, tmp_path / "clean-dist")
    dirty_result, dirty_members = _build_sdist(dirty_source, tmp_path / "dirty-dist")

    assert clean_result.returncode == 0, clean_result.stdout + clean_result.stderr
    assert dirty_result.returncode == 0, dirty_result.stdout + dirty_result.stderr
    assert dirty_members == clean_members
    assert dirty_members.isdisjoint(sentinels)
    required = {
        "pyproject.toml",
        "CMakeLists.txt",
        "README.md",
        "cpp/bindings/pybind11_module.cpp",
        "traminex/py.typed",
        "tests/r_goldens/OMspell_exp05.csv",
        "SOURCES/TraMineR/DESCRIPTION",
        "docs/API_PARITY_DECISIONS.csv",
        "scripts/_checkout_imports.py",
        "scripts/r_goldens/README.md",
        "tests/test_sequence_weight_validation.py",
        "tests/test_remediation_codebook.py",
        "tests/test_remediation_discrepancy.py",
        "tests/test_remediation_distances.py",
        "tests/test_remediation_dot_identifiers.py",
        "tests/test_remediation_gower_weights.py",
        "tests/test_remediation_native_codes.py",
        "tests/test_remediation_reference_route.py",
        "tests/test_remediation_representatives.py",
        "tests/test_remediation_vertical_spells.py",
    }
    assert required <= clean_members


def test_default_wheel_build_is_compiler_free(tmp_path: Path):
    source = tmp_path / "source"
    source.mkdir()
    for name in ("CMakeLists.txt", "pyproject.toml", "README.md"):
        shutil.copy2(ROOT / name, source / name)
    shutil.copytree(
        ROOT / "traminex",
        source / "traminex",
        ignore=shutil.ignore_patterns("*.dylib", "*.pyd", "*.so", "__pycache__"),
    )
    output_dir = tmp_path / "wheelhouse"
    output_dir.mkdir()
    env = os.environ.copy()
    env.update(
        {
            "CC": "/definitely/missing-cc",
            "CXX": "/definitely/missing-cxx",
            "SKBUILD_BUILD_DIR": str(tmp_path / "build"),
        }
    )
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "from scikit_build_core.build import build_wheel; "
            "print(build_wheel(__import__('sys').argv[1]))",
            str(output_dir),
        ],
        cwd=source,
        env=env,
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    wheel = next(output_dir.glob("traminex-0.2.2-*.whl"))
    assert wheel.name.endswith("-py3-none-any.whl")
    with zipfile.ZipFile(wheel) as archive:
        names = archive.namelist()
        metadata_name = next(name for name in names if name.endswith(".dist-info/METADATA"))
        wheel_metadata = archive.read(metadata_name).decode("utf-8")
    assert not any(name.endswith((".so", ".pyd", ".dylib")) for name in names)
    assert "traminex/py.typed" in names
    assert "Requires-Python: >=3.11" in wheel_metadata

    venv = tmp_path / "typing-venv"
    subprocess.run([sys.executable, "-m", "venv", str(venv)], check=True)
    venv_python = venv / "bin" / "python"
    subprocess.run(
        [str(venv_python), "-m", "pip", "install", "--no-deps", str(wheel)],
        check=True,
        capture_output=True,
        text=True,
    )
    consumer = tmp_path / "consumer"
    consumer.mkdir()
    (consumer / "typing_smoke.py").write_text(
        "from traminex.sequences import STSList\n\nvalue: type[STSList] = STSList\n",
        encoding="utf-8",
    )
    (consumer / "mypy.ini").write_text(
        "[mypy]\nignore_missing_imports = True\n",
        encoding="utf-8",
    )
    typing_env = os.environ.copy()
    typing_env.pop("MYPYPATH", None)
    typing_env["PYTHONPATH"] = ""
    typing_result = subprocess.run(
        [
            sys.executable,
            "-m",
            "mypy",
            "--config-file",
            str(consumer / "mypy.ini"),
            "--python-executable",
            str(venv_python),
            str(consumer / "typing_smoke.py"),
        ],
        cwd=consumer,
        env=typing_env,
        check=False,
        capture_output=True,
        text=True,
    )
    assert typing_result.returncode == 0, typing_result.stdout + typing_result.stderr
    assert "Skipping analyzing" not in typing_result.stdout


def test_cmake_exposes_only_an_opt_in_pybind11_profile():
    cmake = (ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    launcher = (ROOT / "scripts" / "run_pybind11_tests.sh").read_text(encoding="utf-8")
    runtime = (ROOT / "traminex" / "distances" / "__init__.py").read_text(encoding="utf-8")
    binding_sources = {
        path.relative_to(ROOT).as_posix() for path in (ROOT / "cpp" / "bindings").glob("*.cpp")
    }

    assert "project(traminex_kernels LANGUAGES NONE)" in cmake
    assert 'option(TRAMINEX_BUILD_PYBIND11 "Build pybind11 bindings" OFF)' in cmake
    assert "enable_language(CXX)" in cmake
    assert binding_sources == {"cpp/bindings/pybind11_module.cpp"}
    assert "nanobind" not in cmake.lower()
    assert "nanobind" not in pyproject.lower()
    assert "nanobind" not in launcher.lower()
    assert "nanobind" not in runtime.lower()
