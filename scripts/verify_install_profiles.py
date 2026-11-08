#!/usr/bin/env python3
"""Verify the six maintained traminex install profiles in fresh environments."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import secrets
import shutil
import subprocess
import tarfile
import tempfile
import textwrap
import zipfile
from dataclasses import dataclass
from datetime import datetime, timezone
from email.parser import Parser
from pathlib import Path
from typing import Any, Sequence

ROOT = Path(__file__).resolve().parents[1]
MARKER_NAME = ".traminex-install-profile-owner"
FLOOR_CONSTRAINTS = (
    "numpy==1.26.0",
    "pandas==2.1.0",
    "scipy==1.11.1",
    "matplotlib==3.6.0",
    "seaborn==0.12.0",
    "numba==0.59.0",
    "scikit-learn==1.3.0",
)


class UnsafeCleanupError(RuntimeError):
    """Raised when a path cannot be proven to be owned by this invocation."""


@dataclass(frozen=True)
class OwnedRoot:
    path: Path
    token: str


@dataclass(frozen=True)
class ProfileDefinition:
    name: str
    artifact_kind: str
    python_slot: str
    extra: str | None = None
    floor_constraints: bool = False
    binding: str = "python"


def profile_definitions() -> tuple[ProfileDefinition, ...]:
    """Return the frozen six-profile R5 install matrix."""
    return (
        ProfileDefinition("core", "pure_wheel", "python311", floor_constraints=True),
        ProfileDefinition(
            "performance",
            "sdist",
            "python311",
            extra="perf",
            floor_constraints=True,
        ),
        ProfileDefinition("plotting", "pure_wheel", "python312", extra="plot"),
        ProfileDefinition("tree", "sdist", "python312", extra="tree"),
        ProfileDefinition(
            "compiled",
            "compiled_wheel",
            "python312",
            binding="pybind11",
        ),
        ProfileDefinition("fallback", "pure_wheel", "python312", binding="auto"),
    )


def create_owned_root(parent: Path | None = None) -> OwnedRoot:
    """Create a new root and an unpredictable ownership marker."""
    if parent is not None:
        parent = parent.expanduser().resolve(strict=True)
        if not parent.is_dir() or parent.is_symlink():
            raise ValueError("temporary parent must be a real directory")
    path = Path(tempfile.mkdtemp(prefix="traminex-r5-install-", dir=parent))
    token = secrets.token_hex(32)
    (path / MARKER_NAME).write_text(token + "\n", encoding="utf-8")
    return OwnedRoot(path=path, token=token)


def remove_owned_root(owned: OwnedRoot) -> None:
    """Delete only a non-symlink directory carrying the matching marker."""
    path = owned.path.expanduser()
    forbidden = {Path("/"), Path.home(), ROOT, ROOT.parent}
    if path in forbidden or len(path.parts) < 3:
        raise UnsafeCleanupError(f"refusing dangerous cleanup path: {path}")
    if path.is_symlink():
        raise UnsafeCleanupError(f"refusing symlink cleanup path: {path}")
    if not path.exists():
        return
    if not path.is_dir():
        raise UnsafeCleanupError(f"owned root is not a directory: {path}")
    marker = path / MARKER_NAME
    try:
        marker_token = marker.read_text(encoding="utf-8").strip()
    except OSError as error:
        raise UnsafeCleanupError(f"ownership marker is missing: {marker}") from error
    if not secrets.compare_digest(marker_token, owned.token):
        raise UnsafeCleanupError(f"ownership marker does not match: {marker}")
    shutil.rmtree(path)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _metadata_fields(content: str) -> tuple[str | None, str | None]:
    metadata = Parser().parsestr(content)
    return metadata.get("Version"), metadata.get("Requires-Python")


def inspect_wheel(path: Path, *, expect_compiled: bool) -> dict[str, Any]:
    path = path.expanduser().resolve(strict=True)
    if not path.is_file() or path.suffix != ".whl":
        raise ValueError(f"wheel artifact is invalid: {path}")
    with zipfile.ZipFile(path) as archive:
        names = archive.namelist()
        metadata_names = [name for name in names if name.endswith(".dist-info/METADATA")]
        if len(metadata_names) != 1:
            raise ValueError(f"wheel must contain exactly one METADATA file: {path}")
        metadata_text = archive.read(metadata_names[0]).decode("utf-8")
    version, requires_python = _metadata_fields(metadata_text)
    binaries = [name for name in names if name.endswith((".so", ".pyd", ".dylib"))]
    kernel_binaries = [name for name in binaries if name.startswith("traminex/_pb_kernels")]
    if expect_compiled:
        if len(kernel_binaries) != 1 or len(binaries) != 1:
            raise ValueError(
                "compiled wheel must contain exactly one _pb_kernels binary, "
                f"found binaries={binaries}"
            )
    elif binaries:
        raise ValueError(f"pure wheel must not contain binaries: {binaries}")
    compiled = expect_compiled
    if "traminex/py.typed" not in names:
        raise ValueError("wheel is missing traminex/py.typed")
    if version != "0.2.3" or requires_python != ">=3.11":
        raise ValueError(
            f"wheel metadata drift: version={version!r}, requires_python={requires_python!r}"
        )
    return {
        "path": str(path),
        "filename": path.name,
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
        "version": version,
        "requires_python": requires_python,
        "compiled_extension": compiled,
        "binary_members": binaries,
        "member_count": len(names),
        "py_typed": True,
    }


def inspect_sdist(path: Path) -> dict[str, Any]:
    path = path.expanduser().resolve(strict=True)
    if not path.is_file() or not path.name.endswith(".tar.gz"):
        raise ValueError(f"sdist artifact is invalid: {path}")
    with tarfile.open(path, "r:gz") as archive:
        members = archive.getmembers()
        metadata_members = [member for member in members if member.name.endswith("/PKG-INFO")]
        if len(metadata_members) != 1:
            raise ValueError(f"sdist must contain exactly one PKG-INFO: {path}")
        extracted = archive.extractfile(metadata_members[0])
        if extracted is None:
            raise ValueError(f"cannot read sdist metadata: {path}")
        metadata_text = extracted.read().decode("utf-8")
    version, requires_python = _metadata_fields(metadata_text)
    names = [member.name for member in members]
    if not any(name.endswith("/traminex/py.typed") for name in names):
        raise ValueError("sdist is missing traminex/py.typed")
    if version != "0.2.3" or requires_python != ">=3.11":
        raise ValueError(
            f"sdist metadata drift: version={version!r}, requires_python={requires_python!r}"
        )
    if any(name.endswith((".so", ".pyd", ".dylib")) for name in names):
        raise ValueError("default sdist unexpectedly contains a compiled extension")
    return {
        "path": str(path),
        "filename": path.name,
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
        "version": version,
        "requires_python": requires_python,
        "compiled_extension": False,
        "member_count": len(names),
        "py_typed": True,
    }


def write_report(destination: Path, report: dict[str, Any]) -> None:
    destination = destination.expanduser()
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _run_command(
    command: Sequence[str],
    *,
    cwd: Path,
    env: dict[str, str],
) -> dict[str, Any]:
    started = time_monotonic()
    completed = subprocess.run(
        list(command),
        cwd=cwd,
        env=env,
        check=False,
        capture_output=True,
        text=True,
    )
    return {
        "command": list(command),
        "cwd": str(cwd),
        "exit_code": completed.returncode,
        "duration_seconds": time_monotonic() - started,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
    }


def time_monotonic() -> float:
    # Isolated for deterministic unit substitution without patching subprocess.
    import time

    return time.monotonic()


def _checked(
    command: Sequence[str],
    *,
    cwd: Path,
    env: dict[str, str],
    records: list[dict[str, Any]],
) -> dict[str, Any]:
    record = _run_command(command, cwd=cwd, env=env)
    records.append(record)
    if record["exit_code"] != 0:
        rendered = " ".join(command)
        raise RuntimeError(
            f"command failed ({record['exit_code']}): {rendered}\n"
            f"{record['stdout']}{record['stderr']}"
        )
    return record


def _common_smoke(profile: ProfileDefinition) -> str:
    # `graphviz` (the Python distribution) is expected absent everywhere,
    # including the `tree` profile: that profile's real requirement is the
    # *system* `dot` executable, which pip cannot install, and traminex never
    # imports the Python package (cycle-6 F-004). Asserting its absence in the
    # tree profile is what keeps the declaration honest.
    expected = {
        "core": {"matplotlib": False, "numba": False, "graphviz": False},
        "performance": {"matplotlib": False, "numba": True, "graphviz": False},
        "plotting": {"matplotlib": True, "numba": False, "graphviz": False},
        "tree": {"matplotlib": False, "numba": False, "graphviz": False},
        "compiled": {"matplotlib": False, "numba": False, "graphviz": False},
        "fallback": {"matplotlib": False, "numba": False, "graphviz": False},
    }[profile.name]
    return textwrap.dedent(f"""
        import importlib.metadata
        import importlib.util
        import json
        import os
        import subprocess
        import sys
        from pathlib import Path

        import numpy as np
        import traminex
        from traminex import distances, seqdef, seqdist
        from traminex.api_manifest import PUBLIC_API

        checkout = Path({str(ROOT)!r}).resolve()
        origin = Path(traminex.__file__).resolve()
        assert not origin.is_relative_to(checkout), (origin, checkout)
        assert importlib.metadata.version("traminex") == "0.2.3"
        assert traminex.__all__ == PUBLIC_API
        assert len(PUBLIC_API) == len(set(PUBLIC_API)) == 226

        presence = {{
            name: importlib.util.find_spec(name) is not None
            for name in ("matplotlib", "numba", "graphviz", "traminex._pb_kernels")
        }}
        expected = {expected!r}
        for name, present in expected.items():
            assert presence[name] is present, (name, presence, expected)

        x = seqdef(np.array([["A", "B", "B"], ["A", "A", "B"]], dtype=object))
        distance = np.asarray(seqdist(x, method="OM", indel=1.0, subm=2.0))
        # A-B-B and A-A-B differ at one position: a substitution (2) or a deletion and an
        # insertion (1 + 1) (TX-AUD-20261019-058).
        assert distance.tolist() == [[0.0, 2.0], [2.0, 0.0]], distance.tolist()

        result = {{
            "profile": {profile.name!r},
            "python": sys.version.split()[0],
            "executable": sys.executable,
            "version": importlib.metadata.version("traminex"),
            "origin": str(origin),
            "backend": distances._impl_name(),
            "api_count": len(PUBLIC_API),
            "optional_modules": presence,
            "distance": distance.tolist(),
        }}
        """)


def _profile_smoke(profile: ProfileDefinition, dot: Path) -> str:
    code = _common_smoke(profile)
    if profile.name == "performance":
        code += textwrap.dedent("""
            from traminex.distances import python_kernels
            assert python_kernels._HAVE_NUMBA is True
            matrix = np.array([[0, 1, 2], [0, 2, 2], [2, 1, 0]], dtype=np.int32)
            python_kernels._HAVE_NUMBA = False
            reference = python_kernels.pairwise(matrix, method="LCS")
            python_kernels._HAVE_NUMBA = True
            accelerated = python_kernels.pairwise(matrix, method="LCS")
            assert np.allclose(reference, accelerated, rtol=1e-10, atol=1e-12)
            result["numba_agreement"] = True
            """)
    elif profile.name == "plotting":
        code += textwrap.dedent("""
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from traminex import sequence_index_plot
            fig, axis = sequence_index_plot(x)
            assert axis is not None and plt.fignum_exists(fig.number)
            plt.close(fig)
            assert not plt.fignum_exists(fig.number)
            result["figure_created_and_closed"] = True
            """)
    elif profile.name == "tree":
        # Exercise the path traminex actually uses: its own `dot` locator and
        # subprocess renderer, not the Python `graphviz` distribution, which the
        # profile deliberately does not install (cycle-6 F-004).
        code += textwrap.dedent(f"""
            import tempfile

            from traminex.disstree_viz import _locate_dot, _render_dot

            assert Path({str(dot)!r}).is_file()
            located = _locate_dot()
            assert located is not None, "traminex could not locate the system dot"
            assert Path(located).resolve() == Path({str(dot)!r}).resolve(), located
            with tempfile.TemporaryDirectory() as workdir:
                source = Path(workdir) / "tree.dot"
                target = Path(workdir) / "tree.svg"
                source.write_text("digraph traminex {{ start -> done }}\\n")
                _render_dot(str(source), str(target), "svg")
                rendered = target.read_bytes()
            assert b"<svg" in rendered and b"start" in rendered
            result["graphviz_svg_bytes"] = len(rendered)
            result["dot_executable"] = located
            """)
    elif profile.name == "compiled":
        code += textwrap.dedent("""
            from traminex import _pb_kernels
            from traminex.distances import python_kernels
            assert distances._impl_name() == "pybind11"
            assert _pb_kernels.accel_openmp_enabled() is True
            rng = np.random.default_rng(20260715)
            for n, length, states in ((8, 6, 3), (5, 0, 2), (4, 9, 5)):
                matrix = rng.integers(0, states, size=(n, length), dtype=np.int32)
                assert np.array_equal(
                    _pb_kernels.pairwise_ham_unit(matrix),
                    python_kernels.pairwise(matrix, method="HAM"),
                )
                assert np.array_equal(
                    _pb_kernels.pairwise_lcs(matrix),
                    python_kernels.pairwise(matrix, method="LCS"),
                )
                assert np.array_equal(
                    _pb_kernels.pairwise_om(matrix, 1.0, 2.0),
                    python_kernels.pairwise(matrix, method="OM", indel=1.0, subm=2.0),
                )
            result["openmp_enabled"] = True
            result["byte_identity_corpora"] = 3
            """)
    elif profile.name == "fallback":
        code += textwrap.dedent("""
            assert presence["traminex._pb_kernels"] is False
            assert distances._impl_name() == "python"
            explicit_env = os.environ.copy()
            explicit_env["TRAMINEX_BINDING"] = "pybind11"
            explicit = subprocess.run(
                [
                    sys.executable,
                    "-I",
                    "-c",
                    "from traminex import distances; distances._impl_name()",
                ],
                cwd=Path(sys.executable).parent,
                env=explicit_env,
                check=False,
                capture_output=True,
                text=True,
            )
            assert explicit.returncode != 0
            assert "pybind11 backend was requested but is unavailable" in explicit.stderr
            result["explicit_failure"] = explicit.stderr.strip().splitlines()[-1]
            """)
    code += "\nprint(json.dumps(result, sort_keys=True))\n"
    return code


def _artifact_requirement(path: Path, extra: str | None) -> str:
    if extra is None:
        return str(path)
    return f"traminex[{extra}] @ {path.as_uri()}"


def _base_environment(binding: str, dot: Path) -> dict[str, str]:
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    env.pop("MYPYPATH", None)
    env["PYTHONPATH"] = ""
    env["MPLBACKEND"] = "Agg"
    env["TRAMINEX_BINDING"] = binding
    env["PATH"] = os.pathsep.join([str(dot.parent), env.get("PATH", "")])
    return env


def verify_profiles(
    *,
    artifacts: dict[str, Path],
    artifact_info: dict[str, dict[str, Any]],
    python311: Path,
    python312: Path,
    uv: Path,
    dot: Path,
    evidence: Path,
    temp_parent: Path | None = None,
    keep: bool = False,
) -> dict[str, Any]:
    owned = create_owned_root(temp_parent)
    report: dict[str, Any] = {
        "schema_version": 1,
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "owned_root": str(owned.path),
        "kept_owned_root": keep,
        "artifacts": artifact_info,
        "profiles": [],
    }
    constraints = owned.path / "floor-constraints.txt"
    constraints.write_text("\n".join(FLOOR_CONSTRAINTS) + "\n", encoding="utf-8")
    interpreters = {"python311": python311, "python312": python312}
    try:
        for definition in profile_definitions():
            profile_dir = owned.path / definition.name
            venv = profile_dir / "venv"
            profile_dir.mkdir()
            commands: list[dict[str, Any]] = []
            entry: dict[str, Any] = {
                "name": definition.name,
                "definition": {
                    "artifact_kind": definition.artifact_kind,
                    "python_slot": definition.python_slot,
                    "extra": definition.extra,
                    "floor_constraints": definition.floor_constraints,
                    "binding": definition.binding,
                },
                "artifact_sha256": artifact_info[definition.artifact_kind]["sha256"],
                "commands": commands,
                "status": "running",
            }
            report["profiles"].append(entry)
            env = _base_environment(definition.binding, dot)
            try:
                _checked(
                    [
                        str(uv),
                        "venv",
                        "--python",
                        str(interpreters[definition.python_slot]),
                        "--seed",
                        str(venv),
                    ],
                    cwd=profile_dir,
                    env=env,
                    records=commands,
                )
                venv_python = venv / "bin" / "python"
                version_record = _checked(
                    [str(venv_python), "--version"],
                    cwd=profile_dir,
                    env=env,
                    records=commands,
                )
                install_command = [str(venv_python), "-m", "pip", "install"]
                if definition.floor_constraints:
                    install_command.extend(["--constraint", str(constraints)])
                install_command.append(
                    _artifact_requirement(artifacts[definition.artifact_kind], definition.extra)
                )
                _checked(
                    install_command,
                    cwd=profile_dir,
                    env=env,
                    records=commands,
                )
                check_record = _checked(
                    [str(venv_python), "-m", "pip", "check"],
                    cwd=profile_dir,
                    env=env,
                    records=commands,
                )
                freeze_record = _checked(
                    [str(venv_python), "-m", "pip", "freeze", "--all"],
                    cwd=profile_dir,
                    env=env,
                    records=commands,
                )
                smoke_path = profile_dir / "smoke.py"
                smoke_path.write_text(_profile_smoke(definition, dot), encoding="utf-8")
                smoke_record = _checked(
                    [str(venv_python), "-I", str(smoke_path)],
                    cwd=profile_dir,
                    env=env,
                    records=commands,
                )
                smoke_output = json.loads(smoke_record["stdout"].strip().splitlines()[-1])
                entry.update(
                    {
                        "status": "passed",
                        "interpreter": version_record["stdout"].strip()
                        or version_record["stderr"].strip(),
                        "pip_check": check_record["stdout"].strip(),
                        "pip_freeze": freeze_record["stdout"].splitlines(),
                        "smoke": smoke_output,
                    }
                )
            except Exception as error:
                entry["status"] = "failed"
                entry["error"] = str(error)
                raise
            finally:
                write_report(evidence, report)
        report["completed_utc"] = datetime.now(timezone.utc).isoformat()
        report["status"] = "passed"
        write_report(evidence, report)
        return report
    finally:
        if not keep:
            remove_owned_root(owned)


def _existing_file(value: str) -> Path:
    candidate = Path(value).expanduser()
    try:
        path = candidate.resolve(strict=True)
    except OSError:
        raise argparse.ArgumentTypeError(f"expected an existing file: {candidate}") from None
    if not path.is_file():
        raise argparse.ArgumentTypeError(f"expected an existing file: {candidate}")
    return path


def _existing_directory(value: str) -> Path:
    candidate = Path(value).expanduser()
    try:
        path = candidate.resolve(strict=True)
    except OSError:
        raise argparse.ArgumentTypeError(f"expected an existing directory: {candidate}") from None
    if not path.is_dir():
        raise argparse.ArgumentTypeError(f"expected an existing directory: {candidate}")
    return path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sdist", required=True, type=_existing_file)
    parser.add_argument("--pure-wheel", required=True, type=_existing_file)
    parser.add_argument("--compiled-wheel", required=True, type=_existing_file)
    parser.add_argument("--python311", required=True, type=_existing_file)
    parser.add_argument("--python312", required=True, type=_existing_file)
    parser.add_argument("--uv", required=True, type=_existing_file)
    parser.add_argument("--dot", required=True, type=_existing_file)
    parser.add_argument("--evidence", required=True, type=Path)
    parser.add_argument("--temp-parent", type=_existing_directory)
    parser.add_argument("--keep", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    artifacts = {
        "sdist": args.sdist,
        "pure_wheel": args.pure_wheel,
        "compiled_wheel": args.compiled_wheel,
    }
    artifact_info = {
        "sdist": inspect_sdist(args.sdist),
        "pure_wheel": inspect_wheel(args.pure_wheel, expect_compiled=False),
        "compiled_wheel": inspect_wheel(args.compiled_wheel, expect_compiled=True),
    }
    verify_profiles(
        artifacts=artifacts,
        artifact_info=artifact_info,
        python311=args.python311,
        python312=args.python312,
        uv=args.uv,
        dot=args.dot,
        evidence=args.evidence,
        temp_parent=args.temp_parent,
        keep=args.keep,
    )
    print(f"verified six isolated install profiles; evidence: {args.evidence}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
