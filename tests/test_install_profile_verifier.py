from __future__ import annotations

import io
import json
import os
import sys
import tarfile
import zipfile
from pathlib import Path

import pytest

from scripts import verify_install_profiles as verifier

METADATA = (
    "Metadata-Version: 2.3\n" "Name: traminex\n" "Version: 0.2.2\n" "Requires-Python: >=3.11\n"
)


def _write_wheel(
    path: Path,
    *,
    metadata: tuple[str, ...] = (METADATA,),
    py_typed: bool = True,
    binaries: tuple[str, ...] = (),
) -> None:
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("traminex/__init__.py", "")
        if py_typed:
            archive.writestr("traminex/py.typed", "")
        for index, content in enumerate(metadata):
            suffix = "" if index == 0 else f"-duplicate-{index}"
            archive.writestr(
                f"traminex{suffix}-0.2.2.dist-info/METADATA",
                content,
            )
        for name in binaries:
            archive.writestr(name, b"ELF")


def _write_sdist(path: Path, *, py_typed: bool = True) -> None:
    with tarfile.open(path, "w:gz") as archive:
        for name, content in (
            ("traminex-0.2.2/PKG-INFO", METADATA.encode()),
            ("traminex-0.2.2/traminex/__init__.py", b""),
        ):
            info = tarfile.TarInfo(name)
            info.size = len(content)
            archive.addfile(info, io.BytesIO(content))
        if py_typed:
            info = tarfile.TarInfo("traminex-0.2.2/traminex/py.typed")
            info.size = 0
            archive.addfile(info, io.BytesIO())


def test_owned_root_cleanup_requires_matching_marker(tmp_path):
    owned = verifier.create_owned_root(tmp_path)
    child = owned.path / "profile"
    child.mkdir()
    (child / "proof.txt").write_text("owned")

    verifier.remove_owned_root(owned)
    assert not owned.path.exists()

    foreign = tmp_path / "foreign"
    foreign.mkdir()
    with pytest.raises(verifier.UnsafeCleanupError, match="ownership marker"):
        verifier.remove_owned_root(verifier.OwnedRoot(foreign, "not-the-marker"))
    assert foreign.is_dir()

    mismatched = verifier.create_owned_root(tmp_path)
    with pytest.raises(verifier.UnsafeCleanupError, match="does not match"):
        verifier.remove_owned_root(verifier.OwnedRoot(mismatched.path, "wrong-token"))
    assert mismatched.path.is_dir()
    verifier.remove_owned_root(mismatched)


def test_owned_root_refuses_symlinks_and_dangerous_paths(tmp_path):
    target = tmp_path / "target"
    target.mkdir()
    link = tmp_path / "link"
    link.symlink_to(target, target_is_directory=True)

    with pytest.raises(verifier.UnsafeCleanupError):
        verifier.remove_owned_root(verifier.OwnedRoot(Path("/"), "token"))
    with pytest.raises(verifier.UnsafeCleanupError):
        verifier.remove_owned_root(verifier.OwnedRoot(link, "token"))


def test_artifact_inspection_hashes_and_classifies_wheels(tmp_path):
    pure = tmp_path / "traminex-0.2.2-py3-none-any.whl"
    compiled = tmp_path / "traminex-0.2.2-cp312-cp312-linux_x86_64.whl"
    _write_wheel(pure)
    _write_wheel(
        compiled,
        binaries=("traminex/_pb_kernels.cpython-312-x86_64-linux-gnu.so",),
    )

    pure_info = verifier.inspect_wheel(pure, expect_compiled=False)
    compiled_info = verifier.inspect_wheel(compiled, expect_compiled=True)

    assert len(pure_info["sha256"]) == 64
    assert pure_info["compiled_extension"] is False
    assert compiled_info["compiled_extension"] is True
    assert pure_info["version"] == compiled_info["version"] == "0.2.2"
    assert pure_info["requires_python"] == ">=3.11"


def test_artifact_inspection_rejects_missing_corrupt_and_wrong_extension(tmp_path):
    missing = tmp_path / "missing.whl"
    corrupt = tmp_path / "corrupt.whl"
    corrupt.write_bytes(b"not a zip archive")
    wrong_extension = tmp_path / "traminex-0.2.2.zip"
    _write_wheel(wrong_extension)

    with pytest.raises(FileNotFoundError):
        verifier.inspect_wheel(missing, expect_compiled=False)
    with pytest.raises(zipfile.BadZipFile):
        verifier.inspect_wheel(corrupt, expect_compiled=False)
    with pytest.raises(ValueError, match="wheel artifact is invalid"):
        verifier.inspect_wheel(wrong_extension, expect_compiled=False)


@pytest.mark.parametrize("metadata_count", [0, 2])
def test_wheel_inspection_requires_exactly_one_metadata_member(tmp_path, metadata_count):
    wheel = tmp_path / "traminex-0.2.2-py3-none-any.whl"
    _write_wheel(wheel, metadata=(METADATA,) * metadata_count)

    with pytest.raises(ValueError, match="exactly one METADATA"):
        verifier.inspect_wheel(wheel, expect_compiled=False)


@pytest.mark.parametrize(
    "metadata",
    [
        METADATA.replace("Version: 0.2.2", "Version: 0.0.9"),
        METADATA.replace("Requires-Python: >=3.11", "Requires-Python: >=3.12"),
    ],
)
def test_wheel_inspection_rejects_metadata_contract_drift(tmp_path, metadata):
    wheel = tmp_path / "traminex-0.2.2-py3-none-any.whl"
    _write_wheel(wheel, metadata=(metadata,))

    with pytest.raises(ValueError, match="metadata drift"):
        verifier.inspect_wheel(wheel, expect_compiled=False)


def test_wheel_inspection_rejects_malformed_binary_profiles(tmp_path):
    pure_with_binary = tmp_path / "traminex-0.2.2-py3-none-any.whl"
    compiled_without_kernel = tmp_path / "traminex-0.2.2-cp312-cp312-linux_x86_64.whl"
    compiled_with_extra_binary = tmp_path / "traminex-extra-cp312-cp312-linux_x86_64.whl"
    _write_wheel(pure_with_binary, binaries=("traminex/unexpected.so",))
    _write_wheel(compiled_without_kernel)
    _write_wheel(
        compiled_with_extra_binary,
        binaries=(
            "traminex/_pb_kernels.cpython-312-x86_64-linux-gnu.so",
            "traminex/unexpected.so",
        ),
    )

    with pytest.raises(ValueError, match="pure wheel"):
        verifier.inspect_wheel(pure_with_binary, expect_compiled=False)
    with pytest.raises(ValueError, match="compiled wheel"):
        verifier.inspect_wheel(compiled_without_kernel, expect_compiled=True)
    with pytest.raises(ValueError, match="compiled wheel"):
        verifier.inspect_wheel(compiled_with_extra_binary, expect_compiled=True)


def test_artifact_inspection_requires_packaged_typing_marker(tmp_path):
    wheel = tmp_path / "traminex-0.2.2-py3-none-any.whl"
    sdist = tmp_path / "traminex-0.2.2.tar.gz"
    _write_wheel(wheel, py_typed=False)
    _write_sdist(sdist, py_typed=False)

    with pytest.raises(ValueError, match="py.typed"):
        verifier.inspect_wheel(wheel, expect_compiled=False)
    with pytest.raises(ValueError, match="py.typed"):
        verifier.inspect_sdist(sdist)


def test_paths_with_spaces_and_shell_characters_remain_literal_arguments(tmp_path):
    hostile = tmp_path / "artifact space;$(touch SHOULD_NOT_EXIST).whl"
    evidence = tmp_path / "evidence space;$(touch ALSO_NOT_CREATED)" / "report.json"
    requirement = verifier._artifact_requirement(hostile, "plot")
    record = verifier._run_command(
        [sys.executable, "-c", "import sys; print(sys.argv[1])", requirement],
        cwd=tmp_path,
        env=os.environ.copy(),
    )

    verifier.write_report(evidence, {"argument": requirement})

    assert record["exit_code"] == 0
    assert record["stdout"].strip() == requirement
    assert record["command"][-1] == requirement
    assert evidence.is_file()
    assert not (tmp_path / "SHOULD_NOT_EXIST").exists()
    assert not (tmp_path / "ALSO_NOT_CREATED").exists()


def test_profile_matrix_is_exact_and_reuses_all_three_artifacts():
    profiles = verifier.profile_definitions()

    assert [profile.name for profile in profiles] == [
        "core",
        "performance",
        "plotting",
        "tree",
        "compiled",
        "fallback",
    ]
    assert {profile.artifact_kind for profile in profiles} == {
        "sdist",
        "pure_wheel",
        "compiled_wheel",
    }
    assert {profile.python_slot for profile in profiles} == {"python311", "python312"}
    assert profiles[1].extra == "perf"
    assert profiles[2].extra == "plot"
    assert profiles[3].extra == "tree"


def test_report_writer_is_deterministic_and_creates_only_parent_dirs(tmp_path):
    destination = tmp_path / "evidence" / "profiles.json"
    report = {"profiles": [{"name": "core", "status": "passed"}], "schema_version": 1}

    verifier.write_report(destination, report)
    first = destination.read_bytes()
    verifier.write_report(destination, report)

    assert destination.read_bytes() == first
    assert json.loads(first) == report


def test_make_install_profile_target_uses_the_maintained_verifier():
    makefile = (Path(__file__).resolve().parents[1] / "Makefile").read_text(encoding="utf-8")

    assert "install-profiles:" in makefile
    assert "scripts/verify_install_profiles.py" in makefile
    for variable in ("SDIST", "PURE_WHEEL", "COMPILED_WHEEL", "PYTHON311", "DOT"):
        assert f"$({variable})" in makefile


def _required_parser_args(tmp_path: Path) -> list[str]:
    executable = str(Path(sys.executable).resolve())
    return [
        "--sdist",
        executable,
        "--pure-wheel",
        executable,
        "--compiled-wheel",
        executable,
        "--python311",
        executable,
        "--python312",
        executable,
        "--uv",
        executable,
        "--dot",
        executable,
        "--evidence",
        str(tmp_path / "evidence.json"),
    ]


def test_main_reports_missing_file_as_concise_argparse_error(tmp_path, capsys):
    missing = tmp_path / "missing artifact.whl"
    arguments = _required_parser_args(tmp_path)
    arguments[1] = str(missing)

    with pytest.raises(SystemExit) as caught:
        verifier.main(arguments)

    stderr = capsys.readouterr().err
    assert caught.value.code == 2
    assert str(missing) in stderr
    assert "expected an existing file" in stderr
    assert "Traceback" not in stderr


def test_main_reports_missing_directory_as_concise_argparse_error(tmp_path, capsys):
    missing = tmp_path / "missing temp parent"
    arguments = [*_required_parser_args(tmp_path), "--temp-parent", str(missing)]

    with pytest.raises(SystemExit) as caught:
        verifier.main(arguments)

    stderr = capsys.readouterr().err
    assert caught.value.code == 2
    assert str(missing) in stderr
    assert "expected an existing directory" in stderr
    assert "Traceback" not in stderr


def _smoke_distance_block() -> str:
    """The smoke code from the seqdef call to the end of its distance checks."""
    code = verifier._common_smoke(verifier.profile_definitions()[0])
    start = code.index("x = seqdef(")
    end = code.index("result = {")
    return code[start:end]


def test_the_install_smoke_checks_the_distance_values(monkeypatch) -> None:
    """A wrong distance must fail the smoke test, not only a wrong shape (TX-AUD-20261019-058).

    A-B-B and A-A-B differ at one position: one substitution costs 2, a deletion plus an
    insertion 1 + 1, so OM with indel 1 and subm 2 gives [[0, 2], [2, 0]].
    """
    import numpy as np

    import traminex

    block = _smoke_distance_block()
    exec(block, {"np": np, "seqdef": traminex.seqdef, "seqdist": traminex.seqdist})

    def wrong(*args, **kwargs):
        return np.array([[0.0, 1.0], [1.0, 0.0]])

    with pytest.raises(AssertionError):
        exec(block, {"np": np, "seqdef": traminex.seqdef, "seqdist": wrong})


def test_make_test_pybind11_passes_its_python_to_the_runner() -> None:
    """``make test-pybind11 PYTHON=...`` ran the runner's own default (TX-AUD-20261011-078)."""
    import subprocess

    root = Path(__file__).resolve().parents[1]
    dry_run = subprocess.run(
        ["make", "-n", "test-pybind11", "PYTHON=/opt/example/bin/python"],
        cwd=root,
        capture_output=True,
        text=True,
        check=True,
    ).stdout

    assert 'TRAMINEX_PYTHON="/opt/example/bin/python" scripts/run_pybind11_tests.sh' in dry_run


def test_the_native_runner_falls_back_to_python3_on_path() -> None:
    """The fallback was one user's interpreter path (TX-AUD-20261019-058)."""
    launcher = (
        Path(__file__).resolve().parents[1] / "scripts" / "run_pybind11_tests.sh"
    ).read_text(encoding="utf-8")

    assert "PY=${TRAMINEX_PYTHON:-python3}" in launcher
    assert "/home/" not in launcher


def test_the_default_install_evidence_is_git_ignored() -> None:
    """The default report path was the repository root, which git does not ignore
    (TX-AUD-20261010-088)."""
    import re
    import subprocess

    root = Path(__file__).resolve().parents[1]
    makefile = (root / "Makefile").read_text(encoding="utf-8")
    default = re.search(r"^INSTALL_PROFILE_EVIDENCE \?= (\S+)$", makefile, re.MULTILINE)

    assert default is not None
    ignored = subprocess.run(
        ["git", "-C", str(root), "check-ignore", "-q", "--no-index", default.group(1)],
        check=False,
    )
    assert ignored.returncode == 0, default.group(1)
